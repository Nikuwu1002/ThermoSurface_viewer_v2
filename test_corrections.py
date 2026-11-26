"""Gezielte Regressionen fuer Punktverlust, Dreiecksbildung und Druckachsen.

Ohne CoolProp/VTK ausfuehrbar: python -m unittest test_corrections -v
"""
import importlib.util
from pathlib import Path
from unittest.mock import patch
import types
import unittest
import sys
import numpy as np

from test_ap13_ap18 import (FakeCP, AnalyticEOS, RHO_G, RHO_L, primitive)
from homogeneous_pvt import (solve_spinodal_branch, build_homogeneous_surfaces,
                             validate_homogeneous_boundaries, count_renderable_triangles)


class Corrections(unittest.TestCase):
    def setUp(self):
        self.T = np.linspace(400., 440., 8)
        psat = 2e6 + 1000*primitive(RHO_G)
        self.l = np.column_stack((self.T, np.full(len(self.T), 1/RHO_L), np.full(len(self.T), psat)))
        self.g = np.column_stack((self.T, np.full(len(self.T), 1/RHO_G), np.full(len(self.T), psat)))
        self.sl, _ = solve_spinodal_branch('Synthetic', self.T, 'liquid', n_scan=64,
                          state_factory=AnalyticEOS, cp_module=FakeCP)
        self.sg, _ = solve_spinodal_branch('Synthetic', self.T, 'vapor', n_scan=64,
                          state_factory=AnalyticEOS, cp_module=FakeCP)

    def test_point_error_does_not_remove_entire_row(self):
        target_v = np.geomspace(self.sl[0, 1], self.sg[0, 1], 20)[9]
        rho_bad = 1./target_v
        class InjectedFault(AnalyticEOS):
            def update(self, mode, x, T):
                if mode == FakeCP.DmolarT_INPUTS and abs(x-rho_bad) < 1e-7:
                    raise ValueError('Gezielt gesperrter einzelner Stuetzpunkt')
                super().update(mode,x,T)
        result, report = build_homogeneous_surfaces(
            'Synthetic', self.l, self.g, self.sl, self.sg, n_points=20,
            include_unstable=True, cp_module=FakeCP, state_factory=InjectedFault)
        info = report['unstable_surface']
        self.assertEqual(info['partial_rows'], len(self.T))
        self.assertEqual(info['valid_rows'], 0)
        self.assertEqual(info['boundary_rows'], len(self.T))
        self.assertEqual(len(info['point_reasons']), len(self.T))
        self.assertGreater(count_renderable_triangles(result['unstable_surface']), 0)
        edges = validate_homogeneous_boundaries(result, self.l, self.g,
                        self.sl, self.sg, (500.,2e6,.2))
        self.assertTrue(edges['unstable_surface']['exact_boundary'])
        self.assertEqual(edges['unstable_surface']['partial_rows'], len(self.T))

    def test_no_false_success_for_empty_surface(self):
        empty = {key: np.full((len(self.T), 5, 3), np.nan) for key in (
            'metastable_liquid_surface','metastable_vapor_surface','unstable_surface')}
        checked = validate_homogeneous_boundaries(empty, self.l, self.g, self.sl, self.sg,
                                                  (500.,2e6,.2))
        for v in checked.values():
            self.assertEqual(v['boundary_rows'], 0)
            self.assertFalse(v['exact_boundary'])
            self.assertEqual(v['valid_nodes'], 0)

    def test_no_triangle_between_isolated_rows(self):
        x = np.full((5,5,3),np.nan)
        for i in (0, 4):
            x[i,:,0] = 400 + i*10
            x[i,:,1] = np.linspace(.1,.5,5)
            x[i,:,2] = 2e6
        self.assertEqual(count_renderable_triangles(x),0)

    def test_3_corner_mesh_and_signed_pressure_axis(self):
        stub = types.ModuleType('pyvista')
        class PolyData:
            def __init__(self, points=None, faces=None):
                self.points = np.asarray(points if points is not None else np.empty((0,3)))
                self.faces = np.asarray(faces if faces is not None else [], dtype=np.int64)
                self.n_cells = len(self.faces)//4
        stub.PolyData = PolyData
        with patch.dict(sys.modules, {'pyvista':stub}):
            spec = importlib.util.spec_from_file_location(
                '_viz_correction', Path(__file__).resolve().parents[1] / '01_Programm' / 'visualization.py')
            viz = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(viz)
            x = np.array([[[400.,.1,2e6],[400.,.2,2e6]],
                          [[420.,.1,2e6],[420.,.2,np.nan]]])
            self.assertEqual(viz.homogeneous_surface_mesh(x, (500.,2e6,.2)).n_cells, 1)
            x[1,1,2] = -2e6
            self.assertEqual(viz.homogeneous_surface_mesh(x, (500.,2e6,.2),'log').n_cells,1)
            self.assertEqual(viz.homogeneous_surface_mesh(x, (500.,2e6,.2),'asinh').n_cells,2)
            self.assertEqual(count_renderable_triangles(x, 'asinh'),2)
            self.assertTrue(np.isfinite(viz.to_xyz(x, (500.,2e6,.2),'asinh')).all())


if __name__ == '__main__': unittest.main()
