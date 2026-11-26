"""AP17/AP18: reproduzierbarer REAL-HEOS-Test mit Auswertung aller vier Fluide.

python 03_Tests_Implementierung/validate_ap17_ap18.py --output 03_Tests_Implementierung/results/ap17_ap18.json
Exitcode 2 = Pflichtabhaengigkeit nicht verfuegbar; 1 = Test nicht bestanden.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import platform
import sys
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path
from time import perf_counter
import numpy as np

PROGRAM_DIR = Path(__file__).resolve().parents[1] / "01_Programm"
sys.path.insert(0, str(PROGRAM_DIR))

from homogeneous_pvt import (solve_spinodal_branch, build_homogeneous_surfaces,
                             validate_homogeneous_boundaries, count_renderable_triangles)


def env_report():
    report={'python':platform.python_version(),'platform':platform.platform()}
    for name in ('numpy','scipy','CoolProp','pyvista','pyvistaqt','PySide6','PyYAML'):
        try: report[name]=version(name)
        except PackageNotFoundError: report[name]=None
    return report


def run_fluid(stoff, *, nT=24, nscan=160, npoints=20):
    from CoolProp import AbstractState
    import CoolProp.CoolProp as cp
    state=AbstractState('HEOS',stoff)
    Tc=float(state.T_critical());pc=float(state.p_critical());rhoc=float(state.rhomolar_critical())
    tmin=max(.9*Tc,1.001*state.Ttriple(),1.001*state.Tmin())
    tmax=min(.995*Tc,state.Tmax())
    if not tmin<tmax:raise ValueError('Keine subkritischen Temperaturstuetzstellen')
    temps=np.linspace(tmin,tmax,nT)
    liquid=np.empty((nT,3));vapor=np.empty((nT,3))
    for i,T in enumerate(temps):
        state.update(cp.QT_INPUTS,0.,float(T));pl=state.p();vl=1/state.rhomolar()
        state.update(cp.QT_INPUTS,1.,float(T));pg=state.p();vg=1/state.rhomolar()
        if not np.isclose(pl,pg,rtol=1e-8,atol=1e-3):
            raise AssertionError('Sättigungsdruckdifferenz')
        liquid[i]=(T,vl,pl); vapor[i]=(T,vg,pl)
    branches={};times={};diagnostics={}
    for branch in ('liquid','vapor'):
        start=perf_counter()
        branches[branch],diagnostics[branch]=solve_spinodal_branch(
           stoff,temps,branch,n_scan=nscan)
        times[branch]=perf_counter()-start
    start=perf_counter()
    meshes, mesh_report=build_homogeneous_surfaces(
        stoff,liquid,vapor,branches['liquid'],branches['vapor'],
        n_points=npoints,include_unstable=True)
    times['surfaces']=perf_counter()-start
    boundary=validate_homogeneous_boundaries(meshes,liquid,vapor,
                  branches['liquid'],branches['vapor'],(Tc,pc,1/rhoc))
    # Exakte Grenztransformation; zur Diagnose max Dimensionless Difference.
    coordinate_errors={}
    for key,a,b in (('metastable_liquid_surface',liquid,branches['liquid']),
                    ('metastable_vapor_surface',branches['vapor'],vapor),
                    ('unstable_surface',branches['liquid'],branches['vapor'])):
        x=meshes[key]
        rows=np.isfinite(x).all(axis=(1,2))
        if np.any(rows):
            def xyz(z):
                return np.stack((np.log10(z[...,1]*rhoc),z[...,0]/Tc,
                                 np.log10(z[...,2]/pc)),axis=-1)
            sample=np.concatenate((x[rows,0,:],x[rows,-1,:]),axis=0)
            reference=np.concatenate((a[rows],b[rows]),axis=0)
            mask=(sample[:,2]>0)&(reference[:,2]>0)
            diff=float(np.max(np.abs(xyz(sample[mask])-xyz(reference[mask])))) if np.any(mask) else None
        else:diff=None
        coordinate_errors[key]=diff
    pvreport={}
    if importlib.util.find_spec('pyvista'):
        from visualization import homogeneous_surface_mesh
        for key,mesh in meshes.items():
            surface=homogeneous_surface_mesh(mesh,(Tc,pc,1/rhoc))
            if surface.n_cells:
                cell_area=surface.compute_cell_sizes(length=False,volume=False).cell_data['Area']
                valid=bool(np.isfinite(cell_area).all() and np.all(cell_area>0))
            else: valid=False
            pvreport[key]={'points':int(surface.n_points), 'cells':int(surface.n_cells),
                           'positive_finite_area':valid}
    mesh_counts={name:count_renderable_triangles(mesh, 'log') for name,mesh in meshes.items()}
    available={name:n>0 for name,n in mesh_counts.items()}
    return {'fluid':stoff,'Tc_K':Tc,'pc_Pa':pc,'n_T':nT,'n_scan':nscan,
            'n_points_per_row':npoints,'spinodal':diagnostics,'surfaces':mesh_report,
            'edge_validation':boundary,'transformed_max_error':coordinate_errors,
            'pyvista':pvreport,'mesh_bytes':int(sum(v.nbytes for v in meshes.values())),
            'elapsed_s':times,'available_surfaces':available,
            'passes_core':bool(available['metastable_liquid_surface'] and
                               available['metastable_vapor_surface'] and
                               all(v is None or v <= 1e-10 for v in coordinate_errors.values()))}


def run_all():
    report={'environment':env_report(), 'fluids':[], 'benchmarks':[], 'status':'INCOMPLETE'}
    if not importlib.util.find_spec('CoolProp'):
        report['status']='BLOCKED_NO_COOLPROP'
        report['reason']='CoolProp HEOS für reale Stoffprüfungen nicht installiert'
        return report
    for stoff in ('Water','R134a','CO2','Argon'):
        try:
            report['fluids'].append(run_fluid(stoff))
        except Exception as exc:
            report['fluids'].append({'fluid':stoff,'status':'FAILED','error':repr(exc)})
    # Die Benchmark-Zeit ist numerisch gemessen, nicht eine Schaetzung.
    for points in ((16,64,12),(32,128,20),(48,256,28)):
        try:
            result=run_fluid('Water',nT=points[0],nscan=points[1],npoints=points[2])
            report['benchmarks'].append({'settings':points,'status':'MEASURED',
               'spinodal_and_surface_s':sum(result['elapsed_s'].values()),
               'array_bytes':result['mesh_bytes'],
               'valid_rows':{k:v['valid_rows'] for k,v in result['surfaces'].items()}})
        except Exception as exc:
            report['benchmarks'].append({'settings':points,'status':'FAILED','error':repr(exc)})
    report['status']='HEOS_CORE_PASS' if all(item.get('passes_core',False)
                      for item in report['fluids']) and all(item['status']=='MEASURED'
                      for item in report['benchmarks']) else 'HEOS_CORE_INCOMPLETE'
    report['note']='Qt-Klicktests, komplette p-v-T-Regression und fehlende externe Baselines separat erforderlich.'
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default=str(Path(__file__).resolve().parent / 'results' /
                                                  'ap17_ap18_validation.json'))
    options=parser.parse_args()
    output=Path(options.output);output.parent.mkdir(parents=True,exist_ok=True)
    report=run_all()
    output.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')
    print(json.dumps({'status':report['status'],'report':str(output)},ensure_ascii=False))
    return 0 if report['status']=='HEOS_CORE_PASS' else (2 if report['status'].startswith('BLOCKED') else 1)

if __name__=='__main__':sys.exit(main())
