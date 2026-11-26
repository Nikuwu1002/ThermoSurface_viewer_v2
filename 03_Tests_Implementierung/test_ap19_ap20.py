"""AP19/AP20: Tests fuer instabile Fortsetzung, Netzkontinuitaet und Abdeckung.

Start: python -m unittest test_ap19_ap20 -v

Synthetische Tests laufen ohne CoolProp/VTK. Die HEOS-Tests (Wasser, CO2, Argon) werden
uebersprungen, wenn CoolProp nicht installiert ist.
"""
from __future__ import annotations
import importlib.util
import sys
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.optimize import brentq

import homogeneous_pvt as hp
from homogeneous_pvt import (solve_spinodal_branch, build_homogeneous_surfaces, grid_triangles,
                             node_validity, mesh_topology_report, surface_coverage,
                             combined_surface_report, validate_unstable_nodes, flag_triangles,
                             critical_approach_report, critical_refinement_rows,
                             unstable_row_layout)
from test_ap13_ap18 import FakeCP, AnalyticEOS, RHO_G, RHO_L, primitive

ROOT = Path(__file__).resolve().parent
KEYS = ('metastable_liquid_surface', 'unstable_surface', 'metastable_vapor_surface')

# dp/drho = K (r-1)(r-3)(r-6)(r-8): Dampfspinodale bei 1 (Maximum), Fluessigkeitsspinodale bei 8
# (Minimum) und ein zusaetzliches Extremapaar bei 3 (Minimum) und 6 (Maximum). Die Extrema
# liegen ausserhalb von [p_spin,l; p_spin,g] (Mehrparameter-EOS-Artefakt in Reinform).
POLY = np.poly1d(np.poly([1., 3., 6., 8.]))
ANTI = POLY.integ()
SCALE = 1e4
BASE = 1e7
ART_RHO_G = 0.8
ART_RHO_L = brentq(lambda r: ANTI(r) - ANTI(ART_RHO_G), 8.0001, 40.)


def art_p(rho):
    return BASE + SCALE * (ANTI(rho) - ANTI(5.))


class ArtifactEOS:
    def __init__(self):
        self.rho, self.temp, self.phase = 1., 400., None
    def T_critical(self): return 500.
    def p_critical(self): return 2e6
    def rhomolar_critical(self): return 2.5
    def Ttriple(self): return 300.
    def Tmin(self): return 250.
    def Tmax(self): return 600.
    def specify_phase(self, code): self.phase = code
    def unspecify_phase(self): self.phase = None
    def update(self, mode, x, T):
        self.temp = T
        self.rho = (ART_RHO_L if x == 0 else ART_RHO_G) if mode == FakeCP.QT_INPUTS else x
    def rhomolar(self): return self.rho
    def p(self): return float(art_p(self.rho))
    def first_partial_deriv(self, *_): return float(SCALE * POLY(self.rho))


def load_legacy():
    path = ROOT / 'regression' / 'ap19_vorher' / 'homogeneous_pvt_legacy.py'
    spec = importlib.util.spec_from_file_location('homogeneous_pvt_legacy_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Fixture:
    """Gemeinsame Randkurven fuer die kuenstliche (Artefakt-)EOS und die regulaere Analytik-EOS."""
    @staticmethod
    def artifact(n_T=8, n_points=21):
        T = np.linspace(400., 440., n_T)
        psat = float(art_p(ART_RHO_G))
        l = np.column_stack((T, np.full(n_T, 1 / ART_RHO_L), np.full(n_T, psat)))
        g = np.column_stack((T, np.full(n_T, 1 / ART_RHO_G), np.full(n_T, psat)))
        kw = dict(n_scan=64, state_factory=ArtifactEOS, cp_module=FakeCP)
        sl, rl = solve_spinodal_branch('Synthetic', T, 'liquid', **kw)
        sg, rg = solve_spinodal_branch('Synthetic', T, 'vapor', **kw)
        return T, l, g, sl, sg, (rl, rg)

    @staticmethod
    def regular(n_T=8):
        T = np.linspace(400., 440., n_T)
        psat = 2e6 + 1000 * primitive(RHO_G)
        l = np.column_stack((T, np.full(n_T, 1 / RHO_L), np.full(n_T, psat)))
        g = np.column_stack((T, np.full(n_T, 1 / RHO_G), np.full(n_T, psat)))
        kw = dict(n_scan=64, state_factory=AnalyticEOS, cp_module=FakeCP)
        sl, _ = solve_spinodal_branch('Synthetic', T, 'liquid', **kw)
        sg, _ = solve_spinodal_branch('Synthetic', T, 'vapor', **kw)
        return T, l, g, sl, sg


def build(module, fixture, factory, n_points=21, **kw):
    T, l, g, sl, sg = fixture[:5]
    return module.build_homogeneous_surfaces('Synthetic', l, g, sl, sg, n_points=n_points,
                                             include_unstable=True, state_factory=factory,
                                             cp_module=FakeCP, **kw)


class SyntheticArtifactEOS(unittest.TestCase):
    """AP19: Zusaetzliche EOS-Extrema duerfen nie als instabile Flaeche erscheinen."""

    def setUp(self):
        self.fix = Fixture.artifact()
        self.T, self.l, self.g, self.sl, self.sg, _ = self.fix
        self.result, self.reports = build(hp, self.fix, ArtifactEOS)

    def test_boundary_order_and_spinodal_values(self):
        self.assertTrue(np.all(self.l[:, 1] < self.sl[:, 1]))
        self.assertTrue(np.all(self.sl[:, 1] < self.sg[:, 1]))
        self.assertTrue(np.all(self.sg[:, 1] < self.g[:, 1]))
        np.testing.assert_allclose(self.sl[:, 1], 1 / 8., rtol=1e-9)
        np.testing.assert_allclose(self.sg[:, 1], 1 / 1., rtol=1e-9)
        self.assertTrue(np.all(self.sl[:, 2] < self.sg[:, 2]))

    def test_spinodal_uses_root_nearest_to_saturation(self):
        # Zwei Minima (rho=3 und 8) im Liquid-Suchintervall: gueltig ist die saettigungsnaechste.
        rep = Fixture.artifact()[5][0]
        self.assertEqual(rep['valid_points'], len(self.T))
        self.assertEqual(len(rep['ambiguous_temperatures']), len(self.T))

    def test_no_node_outside_spinodal_pressure_window(self):
        u = self.result['unstable_surface']
        finite = np.isfinite(u).all(axis=2)
        for i in range(len(self.T)):
            p = u[i, finite[i], 2]
            tol = 1e-9 * max(abs(self.sl[i, 2]), abs(self.sg[i, 2]), 2e6)
            self.assertTrue(np.all(p >= self.sl[i, 2] - tol), (i, p.min(), self.sl[i, 2]))
            self.assertTrue(np.all(p <= self.sg[i, 2] + tol), (i, p.max(), self.sg[i, 2]))

    def test_legacy_reproduces_outliers_and_fix_removes_them(self):
        legacy = load_legacy()
        old, _ = build(legacy, self.fix, ArtifactEOS)
        windows = np.column_stack((self.sl[:, 2], self.sg[:, 2]))
        valid = node_validity(old['unstable_surface'], 'asinh')
        tri, _cells = grid_triangles(valid)
        flagged = flag_triangles(old['unstable_surface'], tri, windows)
        self.assertGreater(len(flagged), 0, 'Legacy-Ausreisser muessen reproduzierbar sein')
        first = flagged[0]
        self.assertEqual(len(first['nodes']), 3)
        self.assertTrue(all(len(n) == 5 for n in first['nodes_outside']))  # (i, j, T, v, p)
        new_valid = node_validity(self.result['unstable_surface'], 'asinh')
        new_tri, _ = grid_triangles(new_valid)
        self.assertEqual(flag_triangles(self.result['unstable_surface'], new_tri, windows), [])

    def test_all_rows_have_gap_and_are_not_bridged(self):
        info = self.reports['unstable_surface']
        self.assertEqual(info['regular_rows'], 0)
        self.assertEqual(len(info['gap_rows']), len(self.T))
        n_left, _ = unstable_row_layout(21)
        u = self.result['unstable_surface']
        self.assertTrue(np.isnan(u[:, n_left, 1:]).all(), 'Trennspalte muss NaN bleiben')
        valid = node_validity(u, 'asinh')
        tri, _cells = grid_triangles(valid)
        columns = tri % u.shape[1]
        crossing = [(c.min() < n_left) and (c.max() > n_left) for c in columns]
        self.assertFalse(any(crossing), 'Dreiecke duerfen die EOS-Artefaktluecke nicht ueberbruecken')
        top = mesh_topology_report(tri, _cells, valid)
        self.assertEqual(top['components'], 2)
        self.assertEqual(top['open_edges_unexplained'], 0)

    def test_branch_ends_exactly_at_window_edge(self):
        u = self.result['unstable_surface']
        n_left, _ = unstable_row_layout(21)
        for i in range(len(self.T)):
            pG, pL = self.sg[i, 2], self.sl[i, 2]
            self.assertAlmostEqual(u[i, n_left - 1, 2], pG, delta=1e-6 * pG)
            self.assertAlmostEqual(u[i, n_left + 1, 2], pL, delta=1e-6 * pG)

    def test_independent_node_validation_passes_and_detects_tampering(self):
        state = ArtifactEOS()
        u = self.result['unstable_surface']
        ok = validate_unstable_nodes(state, FakeCP, u, self.sl, self.sg, l=self.l, g=self.g)
        self.assertTrue(ok['passed'], ok['failures'][:2])
        self.assertGreater(ok['checked'], 50)
        # Knoten im spurioeus-stabilen Bereich (rho=4.5, dp/dv<0) einschleusen
        bad = u.copy()
        T = self.T[0]
        bad[0, 5] = (T, 1 / 4.5, float(art_p(4.5)))               # dp/dv < 0
        bad[0, 6] = (T, 1 / 1.5, float(art_p(1.5)) + 5e6)         # p ausserhalb des Fensters
        bad[0, 7] = (T, 1 / 1.4, float(art_p(1.4)) * 1.01)        # Druck passt nicht zur EOS
        res = validate_unstable_nodes(state, FakeCP, bad, self.sl, self.sg)
        self.assertFalse(res['passed'])
        reasons = ' '.join(r for f in res['failures'] for r in f['reasons'])
        self.assertIn('dp/dv', reasons)
        self.assertIn('p ausserhalb', reasons)
        self.assertIn('weicht von der EOS', reasons)

    def test_report_documents_gap_causes(self):
        gap = self.reports['unstable_surface']['gap_rows']['0']
        self.assertIn('Dampfspinodaldruck', gap['left_stop'])
        self.assertIn('Fluessigkeitsspinodaldruck', gap['right_stop'])
        self.assertGreaterEqual(gap['slope_sign_changes'], 2)
        self.assertLess(gap['scan_p_min'], gap['window'][0])
        self.assertGreater(gap['scan_p_max'], gap['window'][1])

    def test_metastable_surfaces_stay_monotone_and_complete(self):
        for key in ('metastable_liquid_surface', 'metastable_vapor_surface'):
            x = self.result[key]
            self.assertTrue(np.isfinite(x).all(), key)
            self.assertTrue(np.all(np.diff(x[..., 2], axis=1) < 0), key)

    def test_boundary_order_violation_blocks_row(self):
        T, l, g, sl, sg, _ = self.fix
        sl_bad = sl.copy()
        sl_bad[3, 1], sg_bad = sg[3, 1] * 1.1, sg
        _, rep = build_homogeneous_surfaces(
            'Synthetic', l, g, sl_bad, sg_bad, n_points=11, include_unstable=True,
            state_factory=ArtifactEOS, cp_module=FakeCP)
        for key in KEYS:
            self.assertIn('3', rep[key]['order_violations'])
            self.assertIn(T[3], rep[key]['missing_rows'])

    def test_unstable_independent_check_runs_inside_build(self):
        check = self.reports['unstable_surface']['independent_check']
        self.assertTrue(check['passed'], check)
        self.assertGreater(check['checked'], 0)


class RegularEOS(unittest.TestCase):
    """AP20: regulaere Isothermen bleiben unveraendert, Raender sind exakt gemeinsam."""

    def setUp(self):
        self.fix = Fixture.regular()
        self.T, self.l, self.g, self.sl, self.sg = self.fix
        self.result, self.reports = build(hp, self.fix, AnalyticEOS)

    def test_regular_rows_use_unchanged_geometric_grid(self):
        u = self.result['unstable_surface']
        for i in range(len(self.T)):
            np.testing.assert_allclose(u[i, :, 1], np.geomspace(self.sl[i, 1], self.sg[i, 1], 21),
                                       rtol=1e-14)
        self.assertEqual(self.reports['unstable_surface']['regular_rows'], len(self.T))
        self.assertEqual(self.reports['unstable_surface']['gap_rows'], {})
        self.assertTrue(np.isfinite(u).all())

    def test_shared_boundary_nodes_are_bit_identical(self):
        liquid, unstable, vapor = (self.result[k] for k in KEYS)
        np.testing.assert_array_equal(liquid[:, -1], unstable[:, 0])
        np.testing.assert_array_equal(unstable[:, -1], vapor[:, 0])
        np.testing.assert_array_equal(liquid[:, 0], self.l)
        np.testing.assert_array_equal(vapor[:, -1], self.g)
        np.testing.assert_array_equal(unstable[:, 0], self.sl)
        np.testing.assert_array_equal(unstable[:, -1], self.sg)

    def test_combined_surface_is_one_closed_component(self):
        comb = combined_surface_report(self.result, 'asinh')
        self.assertEqual(comb['components'], 1)
        self.assertEqual(comb['components_expected'], 1)
        self.assertEqual(comb['open_edges_between_shared_nodes'], 0)
        self.assertEqual(comb['nonmanifold_edges'], 0)
        for name, seam in comb['seams'].items():
            self.assertTrue(seam['identical_nodes'], name)
            self.assertEqual(seam['missing_left'], [], name)
            self.assertEqual(seam['missing_right'], [], name)
            self.assertEqual(seam['expected_edges'], len(self.T) - 1)

    def test_each_surface_topology_is_clean(self):
        for key in KEYS:
            x = self.result[key]
            valid = node_validity(x, 'asinh')
            tri, cells = grid_triangles(valid)
            xy = np.stack((np.log10(x[..., 1]), x[..., 0]), axis=-1).reshape(-1, 2)
            top = mesh_topology_report(tri, cells, valid, xy=xy)
            self.assertEqual(top['triangles'], 2 * (len(self.T) - 1) * 20, key)
            self.assertEqual(top['components'], 1, key)
            self.assertEqual(top['open_edges_unexplained'], 0, key)
            self.assertEqual(top['open_edges_gap'], 0, key)
            self.assertEqual(top['orientation_conflicts'], 0, key)
            self.assertEqual(top['projected_orientation_flips'], 0, key)
            self.assertEqual(top['nonmanifold_edges'], 0, key)

    def test_coverage_counts_missing_node(self):
        x = self.result['metastable_vapor_surface'].copy()
        full = surface_coverage(x)
        self.assertEqual(full['coverage'], 1.0)
        self.assertEqual(full['cells_evaluable'], (len(self.T) - 1) * 20)
        x[3, 7] = np.nan  # Ein Knoten betrifft vier Zellen
        part = surface_coverage(x)
        self.assertEqual(part['cells_covered'], full['cells_covered'] - 4)
        self.assertEqual(part['cells_partial'], 4)
        self.assertAlmostEqual(part['coverage'], 1 - 4 / full['cells_evaluable'])

    def test_missing_row_is_not_interpolated_or_bridged(self):
        x = self.result['unstable_surface'].copy()
        x[3] = np.nan
        cov = surface_coverage(x)
        self.assertEqual(cov['rows_missing'], 1)
        self.assertEqual(cov['cells_evaluable'], (len(self.T) - 3) * 20)
        tri, cells = grid_triangles(node_validity(x, 'asinh'))
        rows = tri // x.shape[1]
        self.assertFalse(np.any((rows.min(axis=1) < 3) & (rows.max(axis=1) > 3)))

    def test_critical_cap_helpers_and_approach_report(self):
        T, l, g, sl, sg = self.fix
        # Kuenstliche, streng schrumpfende Kurven um (vc=0.2): Schliessung plausibel.
        n = 6
        s = np.linspace(.5, .1, n)
        vl, vg = .2 - s * .05, .2 + s * .05
        Tr = np.linspace(480., 499., n)
        mk = lambda v, p: np.column_stack((Tr, v, p))
        ll, gg = mk(.2 - s * .08, np.full(n, 1.99e6)), mk(.2 + s * .08, np.full(n, 1.99e6))
        sll, sgg = mk(vl, np.full(n, 1.95e6)), mk(vg, np.full(n, 2.03e6))
        rep = critical_approach_report(ll, sll, sgg, gg, (500., 2e6, .2))
        self.assertTrue(rep['closure_supported'], rep['reasons'])
        # Nicht monoton schrumpfend: keine Schliessung
        sgg_bad = sgg.copy()
        sgg_bad[-2, 1] = sgg[-3, 1] * 0.99 + 0.2 * 0.01
        sgg_bad[-1, 1] = sgg[-2, 1] * 1.5
        self.assertFalse(critical_approach_report(ll, sll, sgg_bad, gg, (500., 2e6, .2))['closure_supported'])

    def test_refinement_rows_between_last_row_and_tc(self):
        T_extra, liq, vap = critical_refinement_rows(
            'Synthetic', 440., 4, state_factory=AnalyticEOS, cp_module=FakeCP)
        self.assertEqual(len(T_extra), 4)
        np.testing.assert_allclose(T_extra, 500. - 60. / 2. ** np.arange(1, 5))
        self.assertTrue(np.all(np.diff(T_extra) > 0) and T_extra[0] > 440. and T_extra[-1] < 500.)
        self.assertTrue(np.all(liq[:, 1] < vap[:, 1]))
        self.assertEqual(critical_refinement_rows('Synthetic', 440., 0, state_factory=AnalyticEOS,
                                                  cp_module=FakeCP)[0].size, 0)


class TopologyTools(unittest.TestCase):
    def test_unexplained_open_edge_and_orientation_and_fold_detection(self):
        valid = np.ones((4, 4), dtype=bool)
        tri, cells = grid_triangles(valid)
        base = mesh_topology_report(tri, cells, valid)
        self.assertEqual((base['open_edges_unexplained'], base['orientation_conflicts']), (0, 0))
        self.assertEqual(base['open_edges_frame'], 12)
        # Ein Dreieck einer voll gueltigen Zelle fehlt -> Nahtfehler (unerklaerte offene Kante)
        broken = mesh_topology_report(tri[1:], cells[1:], valid)
        self.assertGreater(broken['open_edges_unexplained'], 0)
        # Umgekehrte Orientierung eines Dreiecks
        flipped = tri.copy()
        flipped[0] = flipped[0][::-1]
        self.assertGreater(mesh_topology_report(flipped, cells, valid)['orientation_conflicts'], 0)
        # Gefaltetes Gitter (Spalten vertauscht) -> Projektionsfaltung
        grid = np.stack(np.meshgrid(np.arange(4.), np.arange(4.), indexing='ij'), axis=-1)
        grid[:, 2], grid[:, 3] = grid[:, 3].copy(), grid[:, 2].copy()
        xy = np.stack((grid[..., 1], grid[..., 0]), axis=-1).reshape(-1, 2)
        self.assertGreater(mesh_topology_report(tri, cells, valid, xy=xy)['projected_orientation_flips'], 0)

    def test_hole_is_classified_as_gap_not_seam_error(self):
        valid = np.ones((5, 5), dtype=bool)
        valid[2, 2] = False
        tri, cells = grid_triangles(valid)
        top = mesh_topology_report(tri, cells, valid)
        self.assertEqual(top['open_edges_unexplained'], 0)
        self.assertGreater(top['open_edges_gap'], 0)
        self.assertEqual(top['components'], 1)
        self.assertEqual(top['triangles'], 2 * 12 + 4)

    def test_isolated_rows_form_separate_components(self):
        valid = np.ones((6, 4), dtype=bool)
        valid[3] = False
        tri, cells = grid_triangles(valid)
        top = mesh_topology_report(tri, cells, valid)
        self.assertEqual(top['components'], 2)


class MeshWithStubbedPyvista(unittest.TestCase):
    def setUp(self):
        stub = types.ModuleType('pyvista')

        class PolyData:
            def __init__(self, points=None, faces=None):
                self.points = np.asarray(points if points is not None else np.empty((0, 3)))
                self.faces = np.asarray(faces if faces is not None else [], dtype=np.int64)
                self.n_cells = len(self.faces) // 4
        stub.PolyData = PolyData
        with patch.dict(sys.modules, {'pyvista': stub}):
            spec = importlib.util.spec_from_file_location(
                '_viz_ap19', ROOT.parent / '01_Programm' / 'visualization.py')
            self.viz = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(self.viz)
        self.critical = (500., 2e6, .2)

    def grid(self):
        x = np.empty((4, 5, 3))
        for i, T in enumerate(np.linspace(450, 480, 4)):
            for j, v in enumerate(np.linspace(.1, .3, 5)):
                x[i, j] = T, v, 1.5e6 + 1e4 * j + 1e3 * i
        return x

    def test_mesh_matches_topology_tool_and_skips_gap_column(self):
        x = self.grid()
        x[:, 2] = np.nan
        mesh = self.viz.homogeneous_surface_mesh(x, self.critical, 'asinh')
        self.assertEqual(mesh.n_cells, 2 * 3 * 2)
        diag = self.viz.homogeneous_mesh_diagnostics(x, self.critical, 'asinh')
        self.assertEqual(diag['components'], 2)
        self.assertEqual(diag['triangles'], mesh.n_cells)
        self.assertEqual(diag['open_edges_unexplained'], 0)

    def test_degenerate_triangles_are_not_rendered(self):
        x = self.grid()
        x[1, 1] = x[0, 1]  # zwei identische Knoten -> Dreiecke mit Flaeche 0
        diag = self.viz.homogeneous_mesh_diagnostics(x, self.critical, 'asinh')
        self.assertGreater(diag['degenerate_removed'], 0)

    def test_critical_cap_adds_fan_to_single_apex(self):
        x = self.grid()
        plain = self.viz.homogeneous_surface_mesh(x, self.critical, 'asinh')
        capped = self.viz.homogeneous_surface_mesh(x, self.critical, 'asinh', critical_cap=True)
        self.assertEqual(capped.n_cells, plain.n_cells + 4)
        apex = self.viz.to_xyz(np.array([[500., .2, 2e6]]), self.critical, 'asinh')[0]
        hits = np.all(np.isclose(capped.points, apex), axis=1)
        self.assertEqual(int(hits.sum()), 1)


class HeosValidation(unittest.TestCase):
    """AP19.6/AP20.8: echte HEOS-Stoffe, drei Aufloesungen."""

    @classmethod
    def setUpClass(cls):
        try:
            import CoolProp  # noqa: F401
        except ImportError:
            raise unittest.SkipTest('CoolProp nicht installiert')
        import validate_ap19_ap20 as v
        cls.v = v

    def check(self, fluid, resolution, require_outliers=True):
        row = self.v.run(fluid, resolution, load_legacy())
        after, before = row['nachher'], row['vorher']
        u = after['surfaces']['unstable_surface']
        self.assertEqual(u['nodes_outside_spinodal_pressure_window'], 0, (fluid, resolution))
        self.assertEqual(after['unstable_independent_check']['failures'], 0,
                         after['unstable_independent_check'])
        if require_outliers:
            self.assertGreater(before['surfaces']['unstable_surface']
                               ['nodes_outside_spinodal_pressure_window'], 0, (fluid, resolution))
        for key in KEYS:
            s = after['surfaces'][key]
            self.assertEqual(s['topology']['open_edges_unexplained'], 0, (fluid, resolution, key))
            self.assertEqual(s['topology']['orientation_conflicts'], 0)
            self.assertEqual(s['topology']['projected_orientation_flips'], 0)
            if s['coverage']['coverage'] is not None:
                self.assertGreaterEqual(s['coverage']['coverage'], 0.98, (fluid, resolution, key))
        comb = after['combined']
        self.assertEqual(comb['open_edges_between_shared_nodes'], 0)
        for name, seam in comb['seams'].items():
            self.assertTrue(seam['identical_nodes'], name)
            self.assertEqual(seam['missing_left'], [], (fluid, name))
            self.assertEqual(seam['missing_right'], [], (fluid, name))
        if row['rows_all_four_boundaries'] == row['n_T']:
            self.assertEqual(comb['components'], comb['components_expected'], (fluid, resolution))
        return row

    def test_water_three_resolutions(self):
        for resolution in ('niedrig', 'mittel', 'hoch'):
            with self.subTest(resolution=resolution):
                self.check('Water', resolution)

    def test_co2_three_resolutions(self):
        for resolution in ('niedrig', 'mittel', 'hoch'):
            with self.subTest(resolution=resolution):
                self.check('CO2', resolution)

    def test_argon_and_methane_medium(self):
        for fluid in ('Argon', 'Methane'):
            with self.subTest(fluid=fluid):
                self.check(fluid, 'mittel')

    def test_water_recovers_all_spinodal_rows(self):
        row = self.check('Water', 'mittel')
        self.assertEqual(row['rows_all_four_boundaries'], row['n_T'])
        self.assertEqual(row['missing_boundary_temperatures'], [])


if __name__ == '__main__':
    unittest.main()
