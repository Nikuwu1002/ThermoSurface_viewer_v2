"""Tests fuer zustandsflaeche.py (Datenmodell, gemeinsame Randpunkte, Ersatzflaeche, Netz, Druckachse).

Start: python -m unittest test_zustandsflaeche -v      (auch mit pytest lauffaehig)
Synthetische Tests ohne CoolProp/VTK; HEOS-Tests (Wasser, Argon, CO2) werden ohne CoolProp uebersprungen.
"""
from __future__ import annotations
import unittest
import numpy as np

import zustandsflaeche as zf
from test_ap19_ap20 import ArtifactEOS, Fixture, FakeCP, art_p


def artifact_field(n_s=31, n_T=8):
    fix = Fixture.artifact(n_T=n_T)
    T, l, g, sl, sg = fix[:5]
    fld = zf.build_raw_field('Synthetic', T, sl, sg, n_s, state_factory=ArtifactEOS, cp_module=FakeCP)
    return fix, fld


class P1Datenmodell(unittest.TestCase):
    def setUp(self):
        self.fix, self.fld = artifact_field()
        self.T, self.l, self.g, self.sl, self.sg = self.fix[:5]
        self.pc = 2e6

    def test_raw_field_has_no_nan_without_eval_failure(self):
        self.assertTrue(np.all(self.fld.row_status == zf.ROW_OK))
        self.assertTrue(np.isfinite(self.fld.p_raw).all())
        self.assertTrue(np.all(self.fld.eval_status == zf.EVAL_OK))

    def test_filter_sets_only_display_stability_and_code(self):
        before = {k: getattr(self.fld, k).copy() for k in ('v', 'p_raw', 'dpdv', 'eval_status', 'source')}
        digest = zf.raw_digest(self.fld)
        zf.apply_display_filter(self.fld, self.sl, self.sg, self.pc)
        self.assertEqual(zf.raw_digest(self.fld), digest, 'p_raw muss bitgleich bleiben')
        for key, ref in before.items():
            np.testing.assert_array_equal(getattr(self.fld, key), ref)
        c = self.fld.counts()
        self.assertGreater(c['hidden'], 0)
        self.assertGreater(c['shown'], 0)
        # Ausgeblendete Knoten behalten ihren endlichen Rohwert (Filter setzen nie NaN in (a))
        self.assertTrue(np.isfinite(self.fld.p_raw[~self.fld.display]).all())

    def test_raw_values_bitwise_equal_after_triangulation_and_display_view(self):
        digest = zf.raw_digest(self.fld)
        zf.apply_display_filter(self.fld, self.sl, self.sg, self.pc)
        view = self.fld.display_points()
        self.assertTrue(np.isnan(view[~self.fld.display]).all())
        shared = zf.assemble_shared_vertices({'u': self.fld.raw_points()}, [])
        zf.surface_triangles(self.fld.display, shared.ids['u'])
        self.assertEqual(zf.raw_digest(self.fld), digest)

    def test_every_hidden_node_has_a_filter_code_and_shown_nodes_none(self):
        zf.apply_display_filter(self.fld, self.sl, self.sg, self.pc)
        shown = self.fld.display
        self.assertTrue(np.all(self.fld.filter_code[shown] == zf.FILTER_NONE))
        self.assertTrue(np.all(self.fld.filter_code[~shown] != zf.FILTER_NONE))
        codes = set(np.unique(self.fld.filter_code[~shown]))
        self.assertIn(zf.FILTER_STABLE_SIGN, codes)
        self.assertIn(zf.FILTER_PRESSURE_WINDOW, codes)

    def test_shown_nodes_satisfy_stability_and_window(self):
        zf.apply_display_filter(self.fld, self.sl, self.sg, self.pc)
        for i in range(len(self.T)):
            idx = np.flatnonzero(self.fld.display[i])
            inner = idx[(idx > 0) & (idx < self.fld.shape[1] - 1)]
            self.assertTrue(np.all(self.fld.dpdv[i, inner] > 0))
            p = self.fld.p_raw[i, idx]
            tol = 1e-9 * max(abs(self.sl[i, 2]), abs(self.sg[i, 2]), self.pc)
            self.assertTrue(np.all((p >= self.sl[i, 2] - tol) & (p <= self.sg[i, 2] + tol)))
            self.assertEqual(self.fld.stability[i, inner].tolist(), [zf.STAB_UNSTABLE] * len(inner))

    def test_stability_classes_cover_metastable_islands(self):
        zf.apply_display_filter(self.fld, self.sl, self.sg, self.pc)
        self.assertGreater(np.count_nonzero(self.fld.stability == zf.STAB_METASTABLE), 0)
        self.assertGreater(np.count_nonzero(self.fld.stability == zf.STAB_UNSTABLE), 0)
        self.assertEqual(np.count_nonzero(self.fld.stability == zf.STAB_STABLE), 0)

    def test_invalid_evaluation_is_distinct_from_deliberate_hiding(self):
        fix = Fixture.artifact(n_T=4)
        T, l, g, sl, sg = fix[:5]
        target_v = np.geomspace(sl[0, 1], sg[0, 1], 31)[12]

        class Fault(ArtifactEOS):
            def update(self, mode, x, Tk):
                if mode == FakeCP.DmolarT_INPUTS and abs(x - 1. / target_v) < 1e-9:
                    raise ValueError('Homogene EOS nicht auswertbar (gezielt)')
                super().update(mode, x, Tk)
        fld = zf.build_raw_field('Synthetic', T, sl, sg, 31, state_factory=Fault, cp_module=FakeCP)
        self.assertEqual(fld.eval_status[0, 12], zf.EVAL_NONFINITE)
        self.assertTrue(np.isnan(fld.p_raw[0, 12]))
        self.assertIn((0, 12), fld.eval_reason)
        self.assertEqual(int(np.count_nonzero(fld.eval_status != zf.EVAL_OK)), 4)  # eine Stelle je Zeile
        zf.apply_display_filter(fld, sl, sg, 2e6)
        self.assertEqual(fld.filter_code[0, 12], zf.FILTER_NOT_EVALUATED)
        deliberately_hidden = (~fld.display) & (fld.eval_status == zf.EVAL_OK)
        self.assertGreater(int(deliberately_hidden.sum()), 4)
        self.assertTrue(np.all(np.isfinite(fld.p_raw[deliberately_hidden])))

    def test_missing_rows_are_marked_not_filled(self):
        fix = Fixture.artifact(n_T=5)
        T, l, g, sl, sg = fix[:5]
        sl = sl.copy()
        sl[2, 1:] = np.nan
        fld = zf.build_raw_field('Synthetic', T, sl, sg, 21, state_factory=ArtifactEOS, cp_module=FakeCP)
        self.assertEqual(fld.row_status[2], zf.ROW_NO_SPINODAL)
        self.assertTrue(np.isnan(fld.p_raw[2]).all())
        zf.apply_display_filter(fld, sl, sg, 2e6)
        self.assertFalse(fld.display[2].any())
        self.assertTrue(np.all(fld.filter_code[2] == zf.FILTER_NO_ROW))

    def test_source_sum_equals_node_count(self):
        zf.apply_display_filter(self.fld, self.sl, self.sg, self.pc)
        c = self.fld.counts()
        self.assertEqual(c['source_sum'], c['nodes'])
        self.assertEqual(c['by_source']['HEOS'], c['nodes'])
        self.assertEqual(sum(c['by_filter'].values()), c['nodes'])
        self.assertEqual(sum(c['by_stability'].values()), c['nodes'])


class P2RandpunkteUndSpinodalzuordnung(unittest.TestCase):
    def setUp(self):
        self.fix = Fixture.artifact()
        self.T, self.l, self.g, self.sl, self.sg = self.fix[:5]
        self.fld = zf.build_raw_field('Synthetic', self.T, self.sl, self.sg, 25,
                                      state_factory=ArtifactEOS, cp_module=FakeCP)
        self.sub = zf.build_substitute_field(self.T, self.sl, self.sg, (500., 2e6, .4), 25)

    def grids(self):
        from homogeneous_pvt import build_homogeneous_surfaces
        res, _ = build_homogeneous_surfaces('Synthetic', self.l, self.g, self.sl, self.sg, n_points=25,
                                            include_unstable=True, state_factory=ArtifactEOS, cp_module=FakeCP)
        return {'metastable_liquid_surface': res['metastable_liquid_surface'],
                'unstable_heos': self.fld.raw_points(), 'replacement': self.sub.raw_points(),
                'metastable_vapor_surface': res['metastable_vapor_surface']}

    LINKS = [('metastable_liquid_surface', 'unstable_heos'), ('metastable_liquid_surface', 'replacement'),
             ('unstable_heos', 'metastable_vapor_surface'), ('replacement', 'metastable_vapor_surface')]

    def test_boundary_nodes_share_identical_vertex_indices(self):
        grids = self.grids()
        sv = zf.assemble_shared_vertices(grids, self.LINKS)
        liquid, vapor = sv.ids['metastable_liquid_surface'], sv.ids['metastable_vapor_surface']
        np.testing.assert_array_equal(liquid[:, -1], sv.ids['unstable_heos'][:, 0])
        np.testing.assert_array_equal(liquid[:, -1], sv.ids['replacement'][:, 0])
        np.testing.assert_array_equal(vapor[:, 0], sv.ids['unstable_heos'][:, -1])
        np.testing.assert_array_equal(vapor[:, 0], sv.ids['replacement'][:, -1])
        self.assertTrue(all(not m for m in sv.mismatches.values()), sv.mismatches)

    def test_mesh_boundary_equals_spinodal_exactly(self):
        grids = self.grids()
        sv = zf.assemble_shared_vertices(grids, self.LINKS)
        for key in ('replacement', 'unstable_heos'):
            first = sv.vertices[sv.ids[key][:, 0]]
            last = sv.vertices[sv.ids[key][:, -1]]
            np.testing.assert_array_equal(first, self.sl)
            np.testing.assert_array_equal(last, self.sg)
            rel = np.abs(first - self.sl) / np.maximum(np.abs(self.sl), 1e-300)
            self.assertLessEqual(float(rel.max()), 1e-12)

    def test_mismatching_values_keep_separate_vertices_and_are_reported(self):
        grids = self.grids()
        grids['replacement'] = grids['replacement'].copy()
        grids['replacement'][3, 0, 2] *= 1.0000001
        sv = zf.assemble_shared_vertices(grids, self.LINKS)
        self.assertEqual(sv.mismatches['metastable_liquid_surface|replacement'], [3])
        self.assertNotEqual(sv.ids['replacement'][3, 0], sv.ids['metastable_liquid_surface'][3, -1])
        self.assertEqual(sv.ids['replacement'][2, 0], sv.ids['metastable_liquid_surface'][2, -1])

    def test_spinodal_assignment_first_zero_from_saturation_synthetic(self):
        state = ArtifactEOS()
        from test_ap19_ap20 import ART_RHO_G, ART_RHO_L
        rho_liq = zf.first_zero_from_saturation(state, FakeCP, 400., ART_RHO_L, 2.5)
        rho_vap = zf.first_zero_from_saturation(state, FakeCP, 400., ART_RHO_G, 2.5)
        self.assertAlmostEqual(rho_liq, 8.0, places=8)     # nicht die zusaetzliche Nullstelle bei 3
        self.assertAlmostEqual(rho_vap, 1.0, places=8)
        np.testing.assert_allclose(1. / self.sl[:, 1], rho_liq, rtol=1e-9)
        np.testing.assert_allclose(1. / self.sg[:, 1], rho_vap, rtol=1e-9)


class P2HeosSpinodalzuordnung(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import CoolProp.CoolProp as cp
            from coolprop_heos_common import _new_state
        except ImportError:
            raise unittest.SkipTest('CoolProp nicht installiert')
        cls.cp, cls._new_state = cp, staticmethod(_new_state)

    def test_code_spinodals_match_first_zero_rule_water_argon(self):
        import homogeneous_pvt as hp
        for fluid in ('Water', 'Argon'):
            st, sat = self._new_state(fluid), self._new_state(fluid)
            Tc = st.T_critical()
            for rel in (0.8, 0.95, 0.99):
                T = rel * Tc
                sat.update(self.cp.QT_INPUTS, 0., T)
                rl = sat.rhomolar()
                sat.update(self.cp.QT_INPUTS, 1., T)
                rg = sat.rhomolar()
                sl, _ = hp.solve_spinodal_branch(fluid, np.array([T, T * 1.0000001]), 'liquid')
                sg, _ = hp.solve_spinodal_branch(fluid, np.array([T, T * 1.0000001]), 'vapor')
                z_l = zf.first_zero_from_saturation(st, self.cp, T, rl, st.rhomolar_critical())
                z_g = zf.first_zero_from_saturation(st, self.cp, T, rg, st.rhomolar_critical())
                with self.subTest(fluid=fluid, rel=rel):
                    self.assertAlmostEqual(1. / sl[0, 1] / z_l, 1., places=6)
                    self.assertAlmostEqual(1. / sg[0, 1] / z_g, 1., places=6)


class P1HeosPhasenvorgabe(unittest.TestCase):
    def test_explicit_phase_returns_raw_eos_not_flash(self):
        try:
            import CoolProp.CoolProp as cp
            from coolprop_heos_common import _new_state
        except ImportError:
            raise unittest.SkipTest('CoolProp nicht installiert')
        import homogeneous_pvt as hp
        st = _new_state('Water')
        p_raw, _ = hp.homogeneous_point(st, cp, 600., 1. / 21372.7, 'vapor')
        flash = _new_state('Water')
        flash.update(cp.DmolarT_INPUTS, 21372.7, 600.)
        self.assertGreater(abs(p_raw - flash.p()), 1e9)       # roh 1,16 GPa vs. Flash ~12,3 MPa
        self.assertLess(flash.p(), 2e7)


class P3Ersatzflaeche(unittest.TestCase):
    N = 41

    def test_monotone_instabil_for_pV_greater_pL_with_exact_boundaries(self):
        st, reason, v, p, dpdv = zf.substitute_isotherm(2e-5, 3e6, 1.5e-4, 1.4e7, self.N)
        self.assertEqual(st, zf.ROW_OK, reason)
        self.assertEqual((v[0], v[-1], p[0], p[-1]), (2e-5, 1.5e-4, 3e6, 1.4e7))   # bitgleiche Randwerte
        self.assertTrue(np.all(np.diff(v) > 0))
        self.assertTrue(np.all(np.diff(p) >= 0))
        self.assertTrue(np.all(dpdv >= 0))        # (dp/dv)_T >= 0: mechanisch instabil
        self.assertTrue(np.all(dpdv[1:-1] > 0))

    def test_boundary_derivative_is_zero_and_matches_finite_difference(self):
        st, _, v, p, dpdv = zf.substitute_isotherm(1e-5, -2e6, 2e-4, 5e6, self.N)
        self.assertEqual(st, zf.ROW_OK)
        self.assertEqual((dpdv[0], dpdv[-1]), (0., 0.))
        fd = np.gradient(p, v)
        np.testing.assert_allclose(fd[1:-1], dpdv[1:-1], rtol=2e-2)

    def test_curvature_signs_fit_the_spinodals(self):
        # d2p/dv2 > 0 an der Fluessigkeitsspinodale, < 0 an der Dampfspinodale (Hermite: 6 dp (1-2s)/dv^2)
        s = np.linspace(0, 1, 2001)
        d2 = np.gradient(zf.hermite_dh(s), s)
        self.assertGreater(d2[0], 0)
        self.assertLess(d2[-1], 0)

    def test_negative_pressures_are_data_not_errors(self):
        st, _, v, p, _ = zf.substitute_isotherm(2e-5, -9e6, 1.5e-4, 1.1e7, self.N)
        self.assertEqual(st, zf.ROW_OK)
        self.assertLess(p.min(), 0)

    def test_pV_not_greater_than_pL_gives_no_points(self):
        for pL, pV in ((5e6, 5e6), (5e6, 4e6)):
            st, reason, v, p, dpdv = zf.substitute_isotherm(2e-5, pL, 1.5e-4, pV, self.N)
            self.assertEqual(st, zf.ROW_ORDER)
            self.assertIn('nicht konstruierbar', reason)
            self.assertIsNone(v)
            self.assertIsNone(p)

    def test_missing_wrong_or_nonfinite_spinodal_gives_no_points(self):
        self.assertEqual(zf.substitute_isotherm(np.nan, 1e6, 1e-4, 2e6, 9)[0], zf.ROW_NO_SPINODAL)
        self.assertEqual(zf.substitute_isotherm(2e-4, 1e6, 1e-4, 2e6, 9)[0], zf.ROW_ORDER)
        self.assertEqual(zf.substitute_isotherm(1e-5, 1e6, 1e-4, np.inf, 9)[0], zf.ROW_NO_SPINODAL)

    def test_collapse_near_critical_point_by_width_and_pressure(self):
        vc, pc = 5e-5, 2e6
        st, reason, *_ = zf.substitute_isotherm(vc * 0.999, 1.9e6, vc * 1.001, 2.0e6, 9, v_c=vc, p_c=pc)
        self.assertEqual(st, zf.ROW_COLLAPSED)
        self.assertIn('Breite', reason)
        st, reason, *_ = zf.substitute_isotherm(vc * 0.5, 1.99e6, vc * 1.5, 1.99e6 * (1 + 1e-10), 9, v_c=vc, p_c=pc)
        self.assertEqual(st, zf.ROW_COLLAPSED)
        self.assertIn('Druckdifferenz', reason)
        self.assertEqual(zf.substitute_isotherm(vc * 0.7, 1.5e6, vc * 1.3, 2.0e6, 9, v_c=vc, p_c=pc)[0], zf.ROW_OK)

    def test_field_marks_source_substitute_and_hides_nonconstructible_rows(self):
        fix = Fixture.artifact(n_T=6)
        T, l, g, sl, sg = fix[:5]
        sg = sg.copy()
        sg[2, 2] = sl[2, 2] - 1.          # p_V <= p_L in Zeile 2
        sl = sl.copy()
        sl[4, 1:] = np.nan                # fehlende Spinodale in Zeile 4
        fld = zf.build_substitute_field(T, sl, sg, (500., 2e6, .4), 15)
        ok_rows = fld.row_status == zf.ROW_OK
        self.assertTrue(np.all(fld.source[ok_rows][:, 1:-1] == zf.SRC_SUBSTITUTE))
        self.assertTrue(np.all(fld.source[ok_rows][:, [0, -1]] == zf.SRC_HEOS))   # Randknoten = HEOS-Spinodalpunkte
        self.assertEqual(fld.row_status.tolist(), [0, 0, zf.ROW_ORDER, 0, zf.ROW_NO_SPINODAL, 0])
        for i in (2, 4):
            self.assertFalse(fld.display[i].any())
            self.assertTrue(np.isnan(fld.p_raw[i]).all())
            self.assertIn(i, fld.row_reason)
        for i in (0, 1, 3, 5):
            np.testing.assert_array_equal(fld.raw_points()[i, 0], sl[i])
            np.testing.assert_array_equal(fld.raw_points()[i, -1], sg[i])
            self.assertTrue(fld.display[i].all())
        self.assertEqual(fld.counts()['by_source']['HEOS'], 2 * 4)      # 4 OK-Zeilen x 2 Randknoten
        self.assertEqual(fld.counts()['by_source']['Ersatz'], 6 * 15 - 8)

    def test_critical_and_supercritical_rows_get_no_points(self):
        T = np.array([480., 499., 500., 520.])
        sl = np.array([[480., .1, 1.5e6], [499., .18, 1.99e6], [500, np.nan, np.nan], [520, np.nan, np.nan]])
        sg = np.array([[480., .5, 2.5e6], [499., .22, 2.01e6], [500, np.nan, np.nan], [520, np.nan, np.nan]])
        fld = zf.build_substitute_field(T, sl, sg, (500., 2e6, .2), 9, eps_width=1e-3)
        self.assertEqual(fld.row_status.tolist(), [zf.ROW_OK, zf.ROW_OK, zf.ROW_CRITICAL, zf.ROW_SUPERCRITICAL])
        self.assertFalse(fld.display[2:].any())

    def test_substitute_has_no_derived_thermodynamic_fields(self):
        fld = zf.build_substitute_field(np.array([400.]), np.array([[400., .1, 1e6]]),
                                        np.array([[400., .4, 3e6]]), (500., 2e6, .2), 9)
        forbidden = {'s', 'h', 'u', 'g', 'cp', 'entropy', 'enthalpy'}
        self.assertFalse(forbidden & set(vars(fld)))


class P4Vernetzung(unittest.TestCase):
    CRIT = (500., 2e6, .4)

    def substitute(self, n_T=8, n_s=21):
        fix = Fixture.artifact(n_T=n_T)
        T, l, g, sl, sg = fix[:5]
        return zf.build_substitute_field(T, sl, sg, self.CRIT, n_s)

    def report(self, fld, shown=None):
        shown = fld.display if shown is None else shown
        sv = zf.assemble_shared_vertices({'s': fld.raw_points()}, [])
        tri = zf.surface_triangles(shown, sv.ids['s'])
        xyz = zf.vertices_to_xyz(sv.vertices, self.CRIT, 'asinh')
        v = sv.vertices
        xy = np.column_stack((np.log10(np.where(np.isfinite(v[:, 1]), v[:, 1], 1.)), v[:, 0]))
        return tri, sv, zf.mesh_report(tri, xyz, xy=xy)

    def test_full_grid_has_only_the_outer_boundary_loop(self):
        fld = self.substitute()
        tri, sv, rep = self.report(fld)
        self.assertEqual(rep['triangles'], 2 * 7 * 20)
        self.assertEqual((rep['components'], rep['boundary_loops'], rep['inner_holes']), (1, 1, 0))
        self.assertEqual(rep['open_edges'], 2 * (7 + 20))
        self.assertEqual(rep['degenerate'], 0)
        self.assertEqual((rep['orientation_conflicts'], rep['nonmanifold_edges'],
                          rep['projected_orientation_flips']), (0, 0, 0))

    def test_cells_use_the_fixed_diagonal(self):
        fld = self.substitute(6, 9)
        n_s = fld.shape[1]
        local = zf.grid_triangles_strict(fld.display)
        self.assertEqual(len(local), 2 * 5 * 8)
        for t in local:
            r, c = t // n_s, t % n_s
            diag = [(a, b) for a in range(3) for b in range(3) if r[b] - r[a] == 1 and c[b] - c[a] == 1]
            self.assertEqual(len(diag), 1, t)

    def test_triangle_only_if_all_three_corners_are_shown(self):
        fld = self.substitute(6, 9)
        shown = fld.display.copy()
        shown[3, 4] = False
        sv = zf.assemble_shared_vertices({'s': fld.raw_points()}, [])
        local = zf.grid_triangles_strict(shown)
        self.assertNotIn(3 * 9 + 4, set(local.ravel().tolist()))
        for t in local:
            self.assertTrue(shown.ravel()[t].all())
        tri = zf.surface_triangles(shown, sv.ids['s'])
        rep = zf.mesh_report(tri, zf.vertices_to_xyz(sv.vertices, self.CRIT, 'asinh'))
        self.assertEqual(rep['inner_holes'], 1)       # definierter innerer Rand um den fehlenden Punkt
        self.assertEqual(rep['components'], 1)

    def test_missing_row_gives_two_components_not_a_bridge(self):
        fld = self.substitute(8, 11)
        shown = fld.display.copy()
        shown[4] = False
        tri, sv, rep = self.report(fld, shown)
        rows = (tri // fld.shape[1])
        self.assertFalse(np.any((rows.min(axis=1) < 4) & (rows.max(axis=1) > 4)))
        self.assertEqual((rep['components'], rep['inner_holes']), (2, 0))

    def test_open_edges_match_pyvista_boundary_edges(self):
        try:
            import pyvista as pv
        except ImportError:
            self.skipTest('pyvista fehlt')
        fld = self.substitute(6, 9)
        shown = fld.display.copy()
        shown[2, 3] = False
        tri, sv, rep = self.report(fld, shown)
        xyz = zf.vertices_to_xyz(sv.vertices, self.CRIT, 'asinh')
        faces = np.column_stack((np.full(len(tri), 3), tri)).ravel()
        mesh = pv.PolyData(xyz, faces)
        edges = mesh.extract_feature_edges(boundary_edges=True, feature_edges=False,
                                           manifold_edges=False, non_manifold_edges=False)
        self.assertEqual(edges.n_cells, rep['open_edges'])

    def test_critical_fan_closes_to_one_apex_without_degenerate_triangles(self):
        fld = self.substitute(6, 9)
        sv = zf.assemble_shared_vertices({'s': fld.raw_points()}, [], self.CRIT)
        status = np.concatenate((fld.row_status, [zf.ROW_COLLAPSED, zf.ROW_CRITICAL, zf.ROW_SUPERCRITICAL]))
        ok, last = zf.fan_applicable(status)
        self.assertTrue(ok)
        self.assertEqual(last, 5)
        fan = zf.critical_fan(fld.display[last], sv.ids['s'][last], sv.apex_id)
        self.assertEqual(len(fan), 8)
        tri = np.vstack((zf.surface_triangles(fld.display, sv.ids['s']), fan))
        rep = zf.mesh_report(tri, zf.vertices_to_xyz(sv.vertices, self.CRIT, 'asinh'))
        self.assertEqual((rep['components'], rep['inner_holes'], rep['degenerate']), (1, 0, 0))
        self.assertEqual(rep['orientation_conflicts'], 0)
        self.assertFalse(zf.fan_applicable(np.array([0, 0, zf.ROW_NO_SPINODAL]))[0])
        self.assertFalse(zf.fan_applicable(np.array([0, 0, 0]))[0])


class P5KritischerPunktHeos(unittest.TestCase):
    """T/T_c in {0,6; 0,8; 0,95; 0,99; 1,0; 1,05}: Fächer zum kritischen Punkt ohne degenerierte Dreiecke."""
    REL = (0.6, 0.8, 0.95, 0.99, 1.0, 1.05)

    @classmethod
    def setUpClass(cls):
        try:
            import CoolProp.CoolProp  # noqa: F401
        except ImportError:
            raise unittest.SkipTest('CoolProp nicht installiert')

    def build(self, fluid):
        import homogeneous_pvt as hp
        from coolprop_heos_common import _new_state
        st = _new_state(fluid)
        Tc, pc, vc = float(st.T_critical()), float(st.p_critical()), 1. / float(st.rhomolar_critical())
        T = np.array([r * Tc for r in self.REL])
        sl, _ = hp.solve_spinodal_branch(fluid, T, 'liquid')
        sg, _ = hp.solve_spinodal_branch(fluid, T, 'vapor')
        crit = (Tc, pc, vc)
        return T, sl, sg, crit, zf.build_substitute_field(T, sl, sg, crit, 41)

    def test_fan_has_no_degenerate_triangles_for_water_argon_co2(self):
        for fluid in ('Water', 'Argon', 'CO2'):
            with self.subTest(fluid=fluid):
                T, sl, sg, crit, fld = self.build(fluid)
                self.assertEqual(fld.row_status[4], zf.ROW_CRITICAL)
                self.assertEqual(fld.row_status[5], zf.ROW_SUPERCRITICAL)
                self.assertFalse(fld.display[4:].any())
                ok, last = zf.fan_applicable(fld.row_status)
                self.assertTrue(ok)
                self.assertEqual(last, 3)                    # T/T_c = 0,99
                sv = zf.assemble_shared_vertices({'s': fld.raw_points()}, [], crit)
                tri = np.vstack((zf.surface_triangles(fld.display, sv.ids['s']),
                                 zf.critical_fan(fld.display[last], sv.ids['s'][last], sv.apex_id)))
                for axis in ('asinh', 'symlog', 'linear'):
                    rep = zf.mesh_report(tri, zf.vertices_to_xyz(sv.vertices, crit, axis))
                    self.assertEqual(rep['degenerate'], 0, (fluid, axis, rep))
                    self.assertEqual((rep['components'], rep['inner_holes']), (1, 0))

    def test_rows_closer_than_eps_collapse_instead_of_degenerate_cells(self):
        import homogeneous_pvt as hp
        from coolprop_heos_common import _new_state
        st = _new_state('Water')
        Tc, pc, vc = float(st.T_critical()), float(st.p_critical()), 1. / float(st.rhomolar_critical())
        T = Tc * (1 - np.array([1e-2, 1e-3, 1e-5, 1e-6, 1e-7]))
        sl, _ = hp.solve_spinodal_branch('Water', T, 'liquid')
        sg, _ = hp.solve_spinodal_branch('Water', T, 'vapor')
        fld = zf.build_substitute_field(T, sl, sg, (Tc, pc, vc), 21)
        self.assertEqual(fld.row_status.tolist(), [zf.ROW_OK] * 4 + [zf.ROW_COLLAPSED])
        # Breiten (Water): 1e-6 -> 0,015 v_c (>= eps=1e-2, bleibt Zeile), 1e-7 -> 0,0049 v_c (< eps, kollabiert)
        self.assertTrue(np.isnan(fld.p_raw[4:]).all())
        sv = zf.assemble_shared_vertices({'s': fld.raw_points()}, [], (Tc, pc, vc))
        ok, last = zf.fan_applicable(fld.row_status)
        self.assertTrue(ok)
        tri = np.vstack((zf.surface_triangles(fld.display, sv.ids['s']),
                         zf.critical_fan(fld.display[last], sv.ids['s'][last], sv.apex_id)))
        rep = zf.mesh_report(tri, zf.vertices_to_xyz(sv.vertices, (Tc, pc, vc), 'asinh'))
        self.assertEqual((rep['degenerate'], rep['components'], rep['inner_holes']), (0, 1, 0))


class P6Druckachse(unittest.TestCase):
    PC = 2e6

    def test_negative_zero_and_positive_pressures_stay_finite_and_ordered(self):
        p = np.array([-9e6, -2e5, 0., 2e5, 3e6, 8e7])
        for axis in ('asinh', 'symlog', 'linear'):
            z = zf.pressure_to_z(p, self.PC, axis)
            self.assertTrue(np.isfinite(z).all(), axis)
            self.assertTrue(np.all(np.diff(z) > 0), axis)
            self.assertEqual(float(zf.pressure_to_z(0., self.PC, axis)), 0.)
        # symlog ist punktsymmetrisch und jenseits der Linearschwelle logarithmisch
        s = zf.pressure_to_z(np.array([-1e6, 1e6]), self.PC, 'symlog')
        self.assertAlmostEqual(s[0], -s[1])
        hi = zf.pressure_to_z(np.array([1e7, 1e8]), self.PC, 'symlog')
        self.assertAlmostEqual(hi[1] - hi[0], 1.0, places=2)

    def test_log_axis_marks_nonpositive_as_nan_without_clipping(self):
        z = zf.pressure_to_z(np.array([-1e6, 0., 1e6]), self.PC, 'log')
        self.assertTrue(np.isnan(z[0]) and np.isnan(z[1]) or np.isneginf(z[1]))
        self.assertTrue(np.isfinite(z[2]))
        with self.assertRaises(ValueError):
            zf.pressure_to_z(1., self.PC, 'foo')

    def test_raw_pressure_is_never_modified_by_the_axis_transform(self):
        v = np.array([[400., .1, -3e6], [400., .2, 0.], [400., .3, 4e6]])
        before = v.copy()
        xyz = zf.vertices_to_xyz(v, (500., self.PC, .2), 'symlog')
        np.testing.assert_array_equal(v, before)
        self.assertTrue(np.isfinite(xyz).all())

    def test_ticks_are_labelled_with_real_pressure_values(self):
        for axis in ('asinh', 'symlog'):
            z, labels = zf.pressure_ticks(-5e6, 2e7, self.PC, axis)
            self.assertIn('0 Pa', labels)
            self.assertTrue(any('MPa' in s for s in labels))
            self.assertTrue(any(s.startswith('-') for s in labels))       # negative Werte beschriftet
            for zi, label in zip(z, labels):
                val, unit = label.split()
                value = float(val) * {'Pa': 1., 'kPa': 1e3, 'MPa': 1e6, 'GPa': 1e9}[unit]
                self.assertAlmostEqual(zi, float(zf.pressure_to_z(value, self.PC, axis)), places=9)
        z, labels = zf.pressure_ticks(-5e6, 2e7, self.PC, 'log')
        self.assertTrue(all(not s.startswith('-') and s != '0 Pa' for s in labels))

    def test_mesh_keeps_all_signed_pressure_states_on_symlog_and_drops_them_on_log(self):
        import visualization as viz
        n_T, n_s = 4, 6
        T = np.linspace(400, 460, n_T)
        v = np.linspace(.1, .5, n_s)
        p = np.tile(np.linspace(-4e6, 6e6, n_s), (n_T, 1)) + np.arange(n_T)[:, None] * 1e4
        p[1, 2] = 0.                                   # p = 0 exakt
        grid = np.stack((np.repeat(T[:, None], n_s, 1), np.tile(v, (n_T, 1)), p), axis=-1)
        daten = {'critical_point': (500., self.PC, .2), 'replacement_surface': grid,
                 'pressure_axis': 'symlog'}
        meshes, info = viz.homogeneous_shared_meshes(daten)
        self.assertEqual(meshes['replacement_surface'].n_cells, 2 * (n_T - 1) * (n_s - 1))
        daten['pressure_axis'] = 'log'
        meshes_log, _ = viz.homogeneous_shared_meshes(daten)
        self.assertLess(meshes_log['replacement_surface'].n_cells, 2 * (n_T - 1) * (n_s - 1))
        # Rohwerte bleiben in der Mesh-Punktdatei unveraendert (auch p<0 und p=0)
        np.testing.assert_array_equal(
            np.sort(meshes['replacement_surface'].point_data['p_raw'][:n_T * n_s]), np.sort(p.ravel()))


class P7Kennzeichnung(unittest.TestCase):
    def setUp(self):
        from homogeneous_pvt import build_homogeneous_surfaces
        self.fix = Fixture.artifact(n_T=6)
        T, l, g, sl, sg = self.fix[:5]
        res, _ = build_homogeneous_surfaces('Synthetic', l, g, sl, sg, n_points=15, include_unstable=True,
                                            state_factory=ArtifactEOS, cp_module=FakeCP)
        raw = zf.build_raw_field('Synthetic', T, sl, sg, 15, state_factory=ArtifactEOS, cp_module=FakeCP)
        zf.apply_display_filter(raw, sl, sg, 2e6)
        crit = (500., 2e6, .4)
        sub = zf.build_substitute_field(T, sl, sg, crit, 17)
        self.daten = dict(res, homogeneous_raw_field=raw, replacement_field=sub,
                          replacement_surface=sub.raw_points(), critical_point=crit,
                          pressure_axis='symlog')

    def test_export_has_source_column_and_origin_sum_equals_total(self):
        import csv, tempfile, os
        path = os.path.join(tempfile.mkdtemp(), 'ex.csv')
        stats = zf.export_homogeneous_csv(self.daten, path)
        with open(path, encoding='utf-8') as handle:
            rows = list(csv.DictReader(handle))
        self.assertIn('source', rows[0])
        self.assertEqual(len(rows), stats['rows'])
        self.assertEqual(stats['source_sum'], stats['rows'])
        self.assertEqual({r['source'] for r in rows}, {'HEOS', 'Ersatz'})
        n_sub = sum(1 for r in rows if r['surface'] == 'replacement_surface' and r['source'] == 'Ersatz')
        sub = self.daten['replacement_field']
        self.assertEqual(n_sub, int(np.count_nonzero(sub.source == zf.SRC_SUBSTITUTE)))
        # Ersatzpunkte nie als HEOS exportiert; Randknoten der Ersatzflaeche sind unveraenderte HEOS-Werte
        sub_rows = [r for r in rows if r['surface'] == 'replacement_surface']
        n_s = self.daten['replacement_field'].shape[1]
        for k in range(0, len(sub_rows), n_s):
            block = sub_rows[k:k + n_s]
            if block[0]['display'] == 'True':
                self.assertEqual((block[0]['source'], block[-1]['source']), ('HEOS', 'HEOS'))
                self.assertTrue(all(r['source'] == 'Ersatz' for r in block[1:-1]))
        for r in rows:
            if r['surface'] != 'replacement_surface':
                self.assertEqual(r['source'], 'HEOS')

    def test_all_substitute_vertices_carry_label_in_hover_and_export(self):
        import visualization as viz
        meshes, _ = viz.homogeneous_shared_meshes(self.daten)
        mesh = meshes['replacement_surface']
        pd = mesh.point_data
        sub = self.daten['replacement_field']
        n_labelled = 0
        used = np.unique(mesh.faces.reshape(-1, 4)[:, 1:])
        for pid in used:
            text = viz.hover_text_for_point(mesh, int(pid))
            is_sub = int(pd['source'][pid]) == zf.SRC_SUBSTITUTE
            self.assertEqual(zf.SUBSTITUTE_HOVER_LABEL in text, is_sub)
            n_labelled += is_sub
        interior_shown = int(np.count_nonzero((sub.source == zf.SRC_SUBSTITUTE) & sub.display))
        self.assertEqual(n_labelled, interior_shown)               # 100 % der Ersatzpunkte gekennzeichnet
        for key in ('metastable_liquid_surface', 'metastable_vapor_surface'):
            other = meshes[key]
            for pid in np.unique(other.faces.reshape(-1, 4)[:, 1:]):
                self.assertNotIn(zf.SUBSTITUTE_HOVER_LABEL, viz.hover_text_for_point(other, int(pid)))

    def test_heos_and_substitute_are_never_shown_together(self):
        states = {'unstable_surface': True, 'replacement_surface': False}
        out = zf.resolve_homogeneous_choice({'unstable_surface': True, 'replacement_surface': True}, 'replacement_surface')
        self.assertEqual((out['unstable_surface'], out['replacement_surface']), (False, True))
        out = zf.resolve_homogeneous_choice({'unstable_surface': True, 'replacement_surface': True}, 'unstable_surface')
        self.assertEqual((out['unstable_surface'], out['replacement_surface']), (True, False))
        self.assertEqual(zf.resolve_homogeneous_choice(states, 'metastable_liquid_surface'), states)

    def test_pick_hover_text_offscreen(self):
        import pyvista as pv
        import visualization as viz
        meshes, _ = viz.homogeneous_shared_meshes(self.daten)
        plotter = pv.Plotter(off_screen=True, window_size=(400, 400))
        plotter.add_mesh(meshes['replacement_surface'], color='gold')
        plotter.reset_camera()
        plotter.render()
        w, h = plotter.window_size
        found = None
        for x in range(w // 4, 3 * w // 4, 20):
            for y in range(h // 4, 3 * h // 4, 20):
                text = viz.pick_hover_text(plotter, x, y)
                if text:
                    found = text
                    break
            if found:
                break
        plotter.close()
        self.assertIsNotNone(found)
        self.assertTrue('T=' in found and 'p=' in found)


def _load_script(name):
    import importlib.util
    from pathlib import Path
    path = Path(__file__).resolve().parent / 'scripts' / f'{name}.py'
    spec = importlib.util.spec_from_file_location(f'_script_{name}', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class P8Validierungsmatrix(unittest.TestCase):
    """Matrix aus scripts/validate_zustandsflaeche.py (niedrig/mittel; 'hoch' laeuft im Skript)."""

    @classmethod
    def setUpClass(cls):
        try:
            import CoolProp.CoolProp  # noqa: F401
        except ImportError:
            raise unittest.SkipTest('CoolProp nicht installiert')
        cls.vz = _load_script('validate_zustandsflaeche')

    def check(self, r, label):
        s, c = r['replacement_surface'], r['combined']
        self.assertEqual(s['inner_holes'], 0, label)
        self.assertEqual(c['inner_holes'], 0, label)
        self.assertEqual((s['degenerate'], c['degenerate']), (0, 0), label)
        self.assertEqual((s['orientation_conflicts'], s['nonmanifold_edges'], s['projected_orientation_flips']),
                         (0, 0, 0), label)
        self.assertLessEqual(r['boundary_rel_deviation'], 1e-12, label)
        self.assertTrue(r['shared_vertex_indices_identical'], label)
        self.assertTrue(r['boundary_polyline_monotone_in_T'], label)
        self.assertTrue(r['raw_bitwise_unchanged'], label)
        for key in ('origin_sum_substitute', 'origin_sum_raw'):
            self.assertEqual(r[key][0], r[key][1], label)
        self.assertEqual(sum(r['shared_index_mismatches'].values()), 0, label)

    def test_matrix_low_and_medium_resolution_linear_and_log(self):
        for fluid in ('Water', 'Argon', 'CO2'):
            for name in ('niedrig', 'mittel'):
                n_T, n_s = self.vz.RESOLUTIONS[name]
                ctx = self.vz.prepare(fluid, n_T, n_s)
                for axis in ('linear', 'log'):
                    with self.subTest(fluid=fluid, resolution=name, axis=axis):
                        r = self.vz.evaluate(ctx, axis)
                        self.check(r, (fluid, name, axis))
                        if axis == 'log':
                            self.assertEqual(r['axis_effective'], 'symlog')   # p<=0 vorhanden -> Umschaltung
                        self.assertTrue(r['fan'])

    def test_true_log_axis_when_all_pressures_positive(self):
        for fluid in ('Water', 'Argon', 'CO2'):
            with self.subTest(fluid=fluid):
                ctx = self.vz.prepare(fluid, 60, 120, t_lo=0.97)
                r = self.vz.evaluate(ctx, 'log')
                self.assertEqual(r['axis_effective'], 'log')
                self.assertFalse(r['axis_notice'])
                self.check(r, (fluid, 'log-p>0'))

    def test_single_isotherms_status(self):
        rows = {f: self.vz.single_isotherms(f) for f in ('Water', 'Argon')}
        for fluid, res in rows.items():
            status = {r['T_over_Tc']: r['status'] for r in res}
            self.assertEqual([status[x] for x in (0.6, 0.8, 0.95, 0.99)], ['ok'] * 4, fluid)
            self.assertEqual(status[1.0], zf.ROW_LABELS[zf.ROW_CRITICAL])
            self.assertEqual(status[1.05], zf.ROW_LABELS[zf.ROW_SUPERCRITICAL])


class P8EndeZuEnde(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import CoolProp.CoolProp  # noqa: F401
            import main as app
        except ImportError:
            raise unittest.SkipTest('CoolProp/PySide6 nicht installiert')
        cls.app = app

    def compute(self, **overrides):
        app = self.app
        stoff, raster, visu = app.lade_konfiguration()
        original = app.lade_diagramm_optionen
        options = original('pvt')
        options.update(overrides)
        app.lade_diagramm_optionen = lambda t: dict(options) if t == 'pvt' else original(t)
        try:
            return app.berechne_auswahl(stoff, raster, visu, 'pvt')
        finally:
            app.lade_diagramm_optionen = original

    def test_log_request_switches_to_symlog_with_notice_and_keeps_all_nodes(self):
        import visualization as viz
        daten = self.compute(pressure_axis='log', auto_pressure_axis=True, include_replacement=True)
        self.assertEqual(daten['pressure_axis'], 'symlog')
        self.assertIn('symlog', daten['pressure_axis_notice'])
        sub = daten['replacement_field']
        self.assertGreater(int(np.count_nonzero(sub.p_raw[sub.display] <= 0)), 0)
        meshes, info = viz.homogeneous_shared_meshes(daten)
        used = np.unique(meshes['replacement_surface'].faces.reshape(-1, 4)[:, 1:])
        # jeder gezeigte Ersatzknoten (p<0, p=0, p>0) ist Teil des Netzes (+1 Apex-Vertex bei Faecher)
        expected = int(np.count_nonzero(sub.display)) + (1 if info['replacement_surface']['fan'] else 0)
        self.assertEqual(len(used), expected)

    def test_auto_axis_can_be_disabled_and_log_then_drops_nonpositive_nodes(self):
        import visualization as viz
        daten = self.compute(pressure_axis='log', auto_pressure_axis=False)
        self.assertEqual(daten['pressure_axis'], 'log')
        self.assertIsNone(daten['pressure_axis_notice'])
        meshes, _ = viz.homogeneous_shared_meshes(daten)
        full = self.compute(pressure_axis='asinh')
        meshes_full, _ = viz.homogeneous_shared_meshes(full)
        self.assertLess(meshes['replacement_surface'].n_cells, meshes_full['replacement_surface'].n_cells)

    def test_replacement_default_state_and_origin_sums_in_result(self):
        daten = self.compute(pressure_axis='asinh')
        st = daten['replacement_status']
        self.assertEqual(st['state'], 'ok')
        c = daten['replacement_field'].counts()
        self.assertEqual(c['source_sum'], c['nodes'])
        self.assertEqual(daten['homogeneous_raw_field'].counts()['source_sum'], daten['homogeneous_raw_field'].counts()['nodes'])

    def test_regression_against_p0_baseline_is_bitwise_equal(self):
        from pathlib import Path
        import json
        baseline = Path(__file__).resolve().parent / 'regression' / 'p0_baseline.json'
        if not baseline.exists():
            self.skipTest('P0-Baseline fehlt')
        p0 = _load_script('p0_baseline')
        reference = json.loads(baseline.read_text(encoding='utf-8'))
        current = p0.collect(reference['stoff'])
        for name in ('equilibrium_surface', 'single_phase_surface'):
            self.assertEqual(current[name]['n_points'], reference[name]['n_points'], name)
            self.assertEqual(current[name]['n_cells'], reference[name]['n_cells'], name)
            np.testing.assert_array_equal(current[name]['bounds'], reference[name]['bounds'])
        for diagram, entries in reference['diagrams'].items():
            for key, ref in entries.items():
                self.assertEqual(current['diagrams'][diagram].get(key), ref, f'{diagram}.{key}')


if __name__ == '__main__':
    unittest.main()
