"""Reproduzierbare HEOS-Aufloesungsdiagnose, ohne GUI/PyVista.

Aufruf: python 02_Tests_Archiv/diagnose_resolution.py --fluid Water --nT 24 --nscan 160
Erfordert CoolProp und scipy. Nur Diagnose, keine Modellparameter veraendern.
"""
from __future__ import annotations
import argparse
import json
import platform
import sys
from pathlib import Path
import numpy as np

PROGRAM_DIR = Path(__file__).resolve().parents[1] / "01_Programm"
sys.path.insert(0, str(PROGRAM_DIR))

from homogeneous_pvt import (solve_spinodal_branch, build_homogeneous_surfaces,
                             validate_homogeneous_boundaries, count_renderable_triangles)


def inspect(fluid='Water', nT=24, nscan=160, tmin=None, tmax=None, samples=(5,20,40)):
    from CoolProp import AbstractState
    import CoolProp.CoolProp as cp
    state=AbstractState('HEOS',fluid)
    Tc=float(state.T_critical()); pc=float(state.p_critical())
    rhoc=float(state.rhomolar_critical())
    low=max(.9*Tc, 1.001*state.Ttriple(), 1.001*state.Tmin()) if tmin is None else tmin
    high=min(.995*Tc,state.Tmax()) if tmax is None else tmax
    if not max(state.Ttriple(),state.Tmin())<low<high<Tc:
        raise ValueError('Verlange einen gueltigen unterkritischen Temperaturbereich')
    temps=np.linspace(low,high,nT)
    sat_l=np.empty((nT,3));sat_g=np.empty((nT,3))
    for i,T in enumerate(temps):
        state.update(cp.QT_INPUTS,0.,float(T));psat=float(state.p());vl=1./float(state.rhomolar())
        state.update(cp.QT_INPUTS,1.,float(T));pv=float(state.p());vg=1./float(state.rhomolar())
        if not np.isclose(psat,pv,rtol=1e-8,atol=1e-3):
            raise AssertionError('Sattigungsdruck inkonsistent')
        sat_l[i]=(T,vl,psat);sat_g[i]=(T,vg,psat)
    sl,rl=solve_spinodal_branch(fluid,temps,'liquid',n_scan=nscan)
    sg,rg=solve_spinodal_branch(fluid,temps,'vapor',n_scan=nscan)
    output={'fluid':fluid,'python':platform.python_version(),
            'critical':{'Tc':Tc,'pc':pc,'vc':1/rhoc},
            'T_range':[float(low),float(high)],
            'last_subcritical_gap_K':float(Tc-temps[-1]),
            'nT':nT,'nscan':nscan,'spinodal':{'liquid':rl,'vapor':rg},
            'resolutions':[]}
    for n in samples:
        surfaces,reports=build_homogeneous_surfaces(fluid,sat_l,sat_g,sl,sg,
                                   n_points=n,include_unstable=True)
        checks=validate_homogeneous_boundaries(surfaces,sat_l,sat_g,sl,sg,(Tc,pc,1/rhoc))
        counts={name:{axis:count_renderable_triangles(mesh,axis)
                      for axis in ('log','asinh')} for name,mesh in surfaces.items()}
        output['resolutions'].append({'n_points':n,'surface_report':reports,
                 'boundary_check':checks,'potential_triangles':counts})
    return output


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fluid',default='Water')
    p.add_argument('--nT',type=int,default=24)
    p.add_argument('--nscan',type=int,default=160)
    p.add_argument('--tmin',type=float,default=None)
    p.add_argument('--tmax',type=float,default=None)
    p.add_argument('--output', default=str(Path(__file__).resolve().parent / 'results' /
                                             'resolution_diagnosis.json'))
    args=p.parse_args()
    result=inspect(args.fluid,args.nT,args.nscan,args.tmin,args.tmax)
    path=Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    for block in result['resolutions']:
        print(block['n_points'],'Volumenpunkte:',block['potential_triangles'])
    print('Diagnose gespeichert:',path)

if __name__=='__main__': main()
