"""Integration ohne CoolProp/Qt: main-Verzweigung, Keys und AP15-Verkabelung.

Nicht als Ersatz fuer echte HEOS-/GUI-Tests gedacht.
"""
import importlib.util
from pathlib import Path
from types import ModuleType
from unittest.mock import patch
import sys
import unittest
import numpy as np
from test_ap13_ap18 import AnalyticEOS, FakeCP, RHO_G, RHO_L, primitive
from homogeneous_pvt import solve_spinodal_branch, build_homogeneous_surfaces


class MainWiring(unittest.TestCase):
    def test_result_dictionary_contains_old_and_new_keys(self):
        pvt=ModuleType('coolprop_heos_pvt')
        pvt.validate_rasters=lambda *_:None
        pvt.validate_heos_temperature_ranges=lambda *_:None
        T=np.array([400.,420.,500.,550.])
        satT=T[:2]
        p_sat=2e6+1000*primitive(RHO_G)
        l=np.column_stack((satT,np.full(2,1/RHO_L),np.full(2,p_sat)))
        g=np.column_stack((satT,np.full(2,1/RHO_G),np.full(2,p_sat)))
        received_surface_cfg=[]
        def fake_build_temperature_grid(cfg, *_):
            received_surface_cfg.append(cfg)
            return T,satT
        pvt.build_pvt_temperature_grid=fake_build_temperature_grid
        requested_isotherm_temperatures=[]
        def fake_multi_isotherm(**kwargs):
            if kwargs.get('T_values') is not None:
                requested_isotherm_temperatures.append(kwargs['T_values'])
            return np.stack((np.repeat(T[:,None],5,axis=1),
                             np.broadcast_to(np.linspace(.1,.6,5),(4,5)),
                             np.full((4,5),p_sat)),axis=-1)
        pvt.multi_isotherm=fake_multi_isotherm
        pvt.multi_binodale=lambda **_:(l,g)
        pvt.berechne_zweiphasenmaske=lambda *_:np.zeros((4,5),bool)
        pvt.multi_isochore=lambda **_:np.empty((0,2,3))
        pvt.multi_isobare=lambda **_:[]
        pvt.calc_two_phase_surface=lambda a,b,n_quality:(np.stack((
                  np.broadcast_to(a[:,None,0],(len(a),n_quality)),
                  a[:,None,1]*(1-np.linspace(0,1,n_quality))+b[:,None,1]*np.linspace(0,1,n_quality),
                  np.broadcast_to(a[:,None,2],(len(a),n_quality))),axis=-1),
                  np.linspace(0,1,n_quality))
        pvt.calc_spinodale=lambda *_:(np.empty((0,3)),np.empty((0,3)))
        pvt.prepare_spinodal_branch=lambda raw,*_:(raw,{'geometry_complete':False})
        pvt.validate_pvt_phase_boundaries=lambda *_:None
        common=ModuleType('coolprop_heos_common')
        common.critical_point=lambda *_:(500.,2e6,.2)
        common.fluid_limits=lambda *_:(300.,250.,600.)
        gui=ModuleType('gui');gui.ThermoSurfaceWindow=object
        pyside=ModuleType('PySide6');widgets=ModuleType('PySide6.QtWidgets')
        widgets.QApplication=object;pyside.QtWidgets=widgets
        doubles={'coolprop_heos_common':common,'coolprop_heos_pvt':pvt,
                 'gui':gui,'PySide6':pyside,'PySide6.QtWidgets':widgets}
        for name in ('pht','pst','hsp','uvt'):
            doubles['coolprop_heos_'+name]=ModuleType('coolprop_heos_'+name)
        with patch.dict(sys.modules,doubles):
            spec=importlib.util.spec_from_file_location(
                '_isolated_main_ap15', Path(__file__).resolve().parents[1] / '01_Programm' / 'main.py')
            app=importlib.util.module_from_spec(spec)
            spec.loader.exec_module(app)
            app.DEBUG_MODE=False
            app.test_spinodale_ableitung=lambda *_:{'liquid':{'passed':False},'vapor':{'passed':False}}
            app.solve_spinodal_branch=lambda stoff,temps,branch,**kw:solve_spinodal_branch(
                stoff,temps,branch,n_scan=48,state_factory=AnalyticEOS,cp_module=FakeCP)
            app.build_homogeneous_surfaces=lambda stoff,*args,**kw:build_homogeneous_surfaces(
                stoff,*args,**kw,state_factory=AnalyticEOS,cp_module=FakeCP)
            from zustandsflaeche import build_raw_field
            app.build_raw_field=lambda stoff,*args,**kw:build_raw_field(
                stoff,*args,**kw,state_factory=AnalyticEOS,cp_module=FakeCP)
            # 'Synthetic' ist kein CoolProp-Stoff: Kritische Zusatzzeilen entfallen im Mock.
            app.critical_refinement_rows=lambda *a,**k:(np.empty(0),np.empty((0,3)),np.empty((0,3)))
            raster={'isotherm_grid':{'T_start_isotherm':500.,'T_stop_isotherm':900.,
                                     'T_step_isotherm':10.},
                    'spinodale_grid':{'n_scan':48,'residual_tolerance':1e-8}}
            visu={'binodale_visu':{},'two_phase_visu':{'Q_points':5},
                  'spinodale_visu':{'T_start_spinodal':400.,'T_stop_spinodal':490.},
                  'isotherm_visu':{'T_start_isotherm':420.,'T_stop_isotherm':480.,
                                   'T_step_isotherm':20.,'n_isotherms':4},
                  'isochore_visu':{},'isobare_visu':{}}
            data=app.berechne_daten('Synthetic',raster,visu)
        self.assertEqual(data['two_phase_surface'].shape,(2,5,3))
        self.assertEqual(received_surface_cfg[-1]['T_start_isotherm'],500.)
        self.assertEqual(received_surface_cfg[-1]['T_stop_isotherm'],900.)
        np.testing.assert_array_equal(
            requested_isotherm_temperatures[-1], np.array([420., 440., 460., 480.]))
        self.assertTrue(data['homogeneous_validation']['available'])
        for key in ('surface','phase_mask','two_phase_surface','Q_values',
                   'binodal_liquid','binodal_vapor','spinodale_liquid',
                   'spinodale_vapor','critical_point','isothermen','isobaren','isochoren',
                   'metastable_liquid_surface','metastable_vapor_surface',
                   'unstable_surface','homogeneous_spinodal_liquid','homogeneous_spinodal_vapor'):
            self.assertIn(key,data)
        np.testing.assert_array_equal(data['metastable_liquid_surface'][:,0],l)
        np.testing.assert_array_equal(data['metastable_vapor_surface'][:,-1],g)

if __name__=='__main__': unittest.main()
