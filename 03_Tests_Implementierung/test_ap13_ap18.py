"""AP13--AP17: unabhaengige analytische EOS- und Geometrietests ohne CoolProp/Qt.

Start: python -m unittest test_ap13_ap18 -v
"""
from __future__ import annotations
import sys
import types
import unittest
from scipy.optimize import brentq
import numpy as np
from homogeneous_pvt import (solve_spinodal_branch, build_homogeneous_surfaces,
                             validate_homogeneous_boundaries, homogeneous_point)


class FakeCP:
    iphase_gas=1; iphase_liquid=2; DmolarT_INPUTS=3; QT_INPUTS=4
    iP=5; iDmolar=6; iT=7


def primitive(rho):
    return rho**3/3-5*rho*rho+21*rho


RHO_G=2.
RHO_L=brentq(lambda r:primitive(r)-primitive(RHO_G),7.1,10.)


class AnalyticEOS:
    """Eine S-Isotherme mit lokalen Extrema bei rho=3 und rho=7.

    Siededichten 2 und RHO_L wurden so gewaehlt, dass die p-Werte identisch
    sind (keine HEOS-Referenz oder Behauptung eines echten Phasengleichgewichts).
    """
    def __init__(self, *, bad_temperature=None, base=2e6, dynamic=False):
        self.base=base; self.rho=1.; self.temp=400.; self.dynamic=dynamic
        self.bad_temperature=bad_temperature; self.phase=None
    def T_critical(self):return 500.
    def p_critical(self):return 2e6
    def rhomolar_critical(self):return 5.
    def Ttriple(self):return 300.
    def Tmin(self):return 250.
    def Tmax(self):return 600.
    def specify_phase(self,code):self.phase=code
    def unspecify_phase(self):self.phase=None
    def update(self,mode,x,T):
        self.temp=T
        if self.bad_temperature is not None and abs(T-self.bad_temperature)<1e-6:
            raise ValueError('Gesperrte Temperatur')
        self.rho=(RHO_L if x==0 else RHO_G) if mode==FakeCP.QT_INPUTS else x
    def rhomolar(self):return self.rho
    def roots(self):
        d=.005*(self.temp-400.) if self.dynamic else 0.
        return 3.+d,7.-d
    def p(self):
        a,b=self.roots(); r=self.rho
        return self.base+1000*(r**3/3-(a+b)*r*r/2+a*b*r)
    def first_partial_deriv(self,*_):
        a,b=self.roots()
        return 1000*(self.rho-a)*(self.rho-b)


class HomogeneousTests(unittest.TestCase):
    def setUp(self):
        self.T=np.linspace(400.,440.,8)
        self.factory=AnalyticEOS
        self.l=np.stack((self.T,np.full(len(self.T),1/RHO_L),
                         np.full(len(self.T),2e6+1000*primitive(RHO_L))),axis=1)
        self.g=np.stack((self.T,np.full(len(self.T),1/RHO_G),
                         np.full(len(self.T),2e6+1000*primitive(RHO_G))),axis=1)
    def solve(self,branch,**kw):
        return solve_spinodal_branch('Synthetic',self.T,branch,n_scan=64,
                  cp_module=FakeCP,state_factory=lambda:AnalyticEOS(**kw))
    def test_spinodal_both_branches(self):
        for branch,root in (('vapor',3.),('liquid',7.)):
            with self.subTest(branch=branch):
                points,report=self.solve(branch)
                self.assertEqual(report['valid_points'],len(self.T),report)
                self.assertLess(report['max_residual'],1e-9)
                np.testing.assert_allclose(points[:,1],1/root,rtol=1e-10)
                np.testing.assert_array_equal(points[:,0],self.T)
                self.assertEqual(report['nonpositive_pressures'],0)
    def test_temperature_dependent_branches_move_toward_critical(self):
        for branch in ('liquid','vapor'):
            pts,rep=self.solve(branch,dynamic=True)
            self.assertEqual(rep['valid_points'],len(self.T),rep)
            d=np.diff(pts[:,1])
            self.assertTrue(np.all(d>0) if branch=='liquid' else np.all(d<0))

    def test_gap_not_bridged(self):
        bad=self.T[3]
        points,report=self.solve('liquid',bad_temperature=bad)
        self.assertTrue(np.isnan(points[3,1:]).all())
        self.assertEqual(report['valid_points'],len(self.T)-1)
        self.assertIn(bad,report['missing_temperatures'])
    def test_all_surfaces_physics_and_exact_boundaries(self):
        sl,_=self.solve('liquid'); sg,_=self.solve('vapor')
        meshes,reports=build_homogeneous_surfaces('Synthetic',self.l,self.g,sl,sg,
                       n_points=15,include_unstable=True,state_factory=AnalyticEOS,cp_module=FakeCP)
        for key in meshes:
            self.assertEqual(reports[key]['valid_rows'],len(self.T),(key,reports[key]))
        result=validate_homogeneous_boundaries(meshes,self.l,self.g,sl,sg,(500.,2e6,.2))
        self.assertEqual(result['unstable_surface']['valid_rows'],len(self.T))
        for key,x in meshes.items():
            for row in x:
                for T,v,p in row[1:-1]:
                    self.assertAlmostEqual(p,2e6+1000*primitive(1/v),places=5)
    def test_nonpositive_pressure_physics_not_deleted(self):
        base=-20000.
        state_factory=lambda:AnalyticEOS(base=base)
        sl, rep=solve_spinodal_branch('Synthetic',self.T,'liquid',n_scan=64,
                                     cp_module=FakeCP,state_factory=state_factory)
        self.assertEqual(rep['valid_points'],len(self.T))
        self.assertLess(sl[0,2],0)
        self.assertEqual(rep['nonpositive_pressures'],len(self.T))
    def test_inconsistent_boundary_rejected(self):
        sl,_=self.solve('liquid');sg,_=self.solve('vapor')
        bad=self.g.copy();bad[:,2]*=1.10
        meshes,reports=build_homogeneous_surfaces('Synthetic',self.l,bad,sl,sg,
                 n_points=10,state_factory=AnalyticEOS,cp_module=FakeCP)
        self.assertEqual(reports['metastable_vapor_surface']['valid_rows'],0)
        self.assertTrue(np.isnan(meshes['metastable_vapor_surface']).all())
    def test_strict_monotone_temperature_validation(self):
        with self.assertRaises(ValueError):
            solve_spinodal_branch('Synthetic',[401,400],'liquid',state_factory=AnalyticEOS,
                                  cp_module=FakeCP)
    def test_render_mesh_skips_invalid_rows_and_negative_p(self):
        # Minimaler PV-Double; VTK/Qt wird fuer diesen Test nicht benoetigt.
        original=sys.modules.get('pyvista')
        original_visualization=sys.modules.get('visualization')
        fake=types.ModuleType('pyvista')
        class Mesh:
            def __init__(self,points=None,faces=None):
                self.points=np.asarray(points if points is not None else np.empty((0,3)))
                self.faces=np.asarray(faces if faces is not None else np.empty(0),dtype=int)
                self.n_cells=len(self.faces)//4
        fake.PolyData=Mesh
        sys.modules['pyvista']=fake
        try:
            import importlib
            if 'visualization' in sys.modules:
                del sys.modules['visualization']
            import visualization
            x=np.empty((5,5,3),float)
            for i,T in enumerate(np.linspace(400,440,5)):
                for j,v in enumerate(np.linspace(.1,.5,5)):
                    x[i,j]=T,v,1e6+i*1e4+j*1e3
            x[2]=np.nan
            mesh=visualization.homogeneous_surface_mesh(x,(500,2e6,.2))
            self.assertEqual(mesh.n_cells,2*(2*4))
            self.assertTrue(np.isfinite(mesh.points).all())
            x[1,1,2]=-1.0
            mesh=visualization.homogeneous_surface_mesh(x,(500,2e6,.2))
            self.assertLess(mesh.n_cells,16)
        finally:
            if original is None:
                sys.modules.pop('pyvista',None)
            else:
                sys.modules['pyvista']=original
            if original_visualization is None:
                sys.modules.pop('visualization',None)
            else:
                sys.modules['visualization']=original_visualization

if __name__=='__main__':unittest.main()
