"""g-p-T und f-v-T: thermodynamische Validierung, Integration (Registry, Adapter, Resolver), GUI.

Alle Referenzwerte stammen aus unabhaengigen CoolProp-Aufrufen (PropsSI bzw. eigene AbstractState-Objekte), nicht aus
den getesteten Modulen. Fluide: Water, CO2, Nitrogen, Methane. Die GUI-Tests brauchen ein OpenGL-Fenster
(nicht QT_QPA_PLATFORM=offscreen) und werden sonst uebersprungen.
"""
import os
import unittest

import numpy as np

import test_state_comparison as tsc
from test_gui_multiview import GUI_SKIP, app_main, spin

try:
    import CoolProp
    import CoolProp.CoolProp as CP
    from CoolProp import AbstractState
    import coolprop_heos_fvt as FVT
    import coolprop_heos_gpt as GPT
    import diagram_coords as dc
    from comparison_view import build_geometry
    from diagram_registry import DIAGRAMS
    from selection_manager import SelectionManager
    from state_store import NodeRef
    import thermo_state as ts
    CP_SKIP = None
except Exception as exc:  # pragma: no cover
    CP_SKIP = f"CoolProp/Module nicht verfuegbar: {exc}"

FLUIDS = ("Water", "CO2", "Nitrogen", "Methane")
_DATA = {}


def data(kind, fluid, **kwargs):
    key = (kind, fluid, tuple(sorted(kwargs.items())))
    if key not in _DATA:
        _DATA[key] = (GPT.berechne_gpt_daten if kind == "gpt" else FVT.berechne_fvt_daten)(fluid, **kwargs)
    return _DATA[key]


def molar_mass(fluid):
    return float(AbstractState("HEOS", fluid).molar_mass())


def ref_state(fluid):
    return AbstractState("HEOS", fluid)


def ref_gf_dmass(fluid, T, rho_mass):
    """(p, g, f, h, s, u) mit unabhaengigem PropsSI bei (T, rho_mass) in massenspezifischen SI-Einheiten."""
    props = {k: CP.PropsSI(k, "T", T, "D", rho_mass, fluid) for k in ("P", "G", "H", "S", "U")}
    return props["P"], props["G"], props["U"] - T * props["S"], props["H"], props["S"], props["U"]


def sample_nodes(surface, rows=(0, 7, 20, 44), cols=(0, 4, 19, 38, 40)):
    out = []
    for i in rows:
        for j in cols:
            if i < surface.shape[0] and j < surface.shape[1]:
                out.append(surface[i, j])
    return out


# Stuetzknoten ohne die Sattigungsrandknoten und ohne den kritischen Zustand (T = Tc, p = pc): dort ist ein
# unabhaengiger (T,p)-Aufruf von PropsSI nicht definiert; diese Knoten werden ueber Q bzw. (T, rho) geprueft.
GPT_SAMPLES = {"surface_liquid": dict(rows=(0, 7, 20, 44), cols=(1, 4, 19, 38, 40)),
               "surface_vapor": dict(rows=(0, 7, 20, 44), cols=(0, 4, 19, 38, 39)),
               "surface_high_temperature": dict(rows=(2, 7, 13, 24), cols=(0, 4, 19, 40, 60, 79))}


def gas_scale(fluid, T):
    """Natuerliche Skala R*T [J/kg] fuer absolute Toleranzen bei Werten nahe null (Referenzzustand)."""
    return 8.31446261815324 / molar_mass(fluid) * T


@unittest.skipIf(CP_SKIP, CP_SKIP or "")
class GptThermodynamics(unittest.TestCase):
    """g-p-T: g = h - T s, Sattigungslinie, Knick, Schwalbenschwanz, Spinodalen."""

    def test_surface_nodes_match_independent_coolprop(self):
        for fluid in FLUIDS:
            d = data("gpt", fluid)
            for name, sel in GPT_SAMPLES.items():
                with self.subTest(fluid=fluid, surface=name):
                    for T, p, g in sample_nodes(d[name], **sel):
                        g_ref = CP.PropsSI("G", "T", T, "P", p, fluid)
                        h = CP.PropsSI("H", "T", T, "P", p, fluid)
                        s = CP.PropsSI("S", "T", T, "P", p, fluid)
                        tol = 1e-8 * max(abs(g_ref), gas_scale(fluid, T))
                        self.assertAlmostEqual(g, g_ref, delta=tol, msg=f"T={T} p={p}")
                        self.assertAlmostEqual(h - T * s, g_ref, delta=tol)  # g = h - T s

    def test_g_equals_f_plus_p_v_on_same_state(self):
        for fluid in FLUIDS:
            d = data("gpt", fluid)
            for name, sel in GPT_SAMPLES.items():
                for T, p, g in sample_nodes(d[name], **sel):
                    with self.subTest(fluid=fluid, surface=name, T=T, p=p):
                        rho = CP.PropsSI("D", "T", T, "P", p, fluid)
                        _p, _g, f, *_ = ref_gf_dmass(fluid, T, rho)
                        self.assertAlmostEqual(g, f + p / rho, delta=1e-8 * max(abs(g), gas_scale(fluid, T)))

    def test_dg_dp_is_specific_volume_and_dg_dT_is_minus_entropy(self):
        for fluid in FLUIDS:
            d = data("gpt", fluid)
            T, p, _ = d["surface_high_temperature"][8, 30]
            dp, dT = 1e-4 * p, 1e-5 * T
            dg_dp = (CP.PropsSI("G", "T", T, "P", p + dp, fluid) - CP.PropsSI("G", "T", T, "P", p - dp, fluid)) / (2 * dp)
            dg_dT = (CP.PropsSI("G", "T", T + dT, "P", p, fluid) - CP.PropsSI("G", "T", T - dT, "P", p, fluid)) / (2 * dT)
            self.assertAlmostEqual(dg_dp * CP.PropsSI("D", "T", T, "P", p, fluid), 1.0, delta=1e-5)
            self.assertAlmostEqual(-dg_dT / CP.PropsSI("S", "T", T, "P", p, fluid), 1.0, delta=1e-5)

    def test_saturation_line_has_equal_g_and_is_a_curve(self):
        for fluid in FLUIDS:
            d = data("gpt", fluid)
            l, g = d["binodal_liquid"], d["binodal_vapor"]
            with self.subTest(fluid=fluid):
                self.assertNotIn("two_phase_surface", d)
                np.testing.assert_array_equal(l[:, :2], g[:, :2])
                np.testing.assert_allclose(l[:, 2], g[:, 2], rtol=1e-7)
                for T, p, gl in l[[0, 10, 30]]:
                    self.assertAlmostEqual(gl / CP.PropsSI("G", "T", T, "Q", 0, fluid), 1.0, delta=1e-8)
                    self.assertAlmostEqual(gl / CP.PropsSI("G", "T", T, "Q", 1, fluid), 1.0, delta=1e-7)
                    self.assertAlmostEqual(p / CP.PropsSI("P", "T", T, "Q", 0, fluid), 1.0, delta=1e-8)

    def test_equilibrium_surface_has_kink_at_saturation(self):
        for fluid in FLUIDS:
            d = data("gpt", fluid)
            M = molar_mass(fluid)
            for i in (2, 20):
                T, psat, gs = d["binodal_liquid"][i]
                v_l = 1.0 / CP.PropsSI("D", "T", T, "Q", 0, fluid)
                v_g = 1.0 / CP.PropsSI("D", "T", T, "Q", 1, fluid)
                p_v, g_v = d["surface_vapor"][i, -2, 1], d["surface_vapor"][i, -2, 2]
                p_l, g_l = d["surface_liquid"][i, 1, 1], d["surface_liquid"][i, 1, 2]
                slope_vapor = (gs - g_v) / (psat - p_v)       # -> v_g
                slope_liquid = (g_l - gs) / (p_l - psat)      # -> v_l
                with self.subTest(fluid=fluid, row=i):
                    self.assertAlmostEqual(slope_vapor / v_g, 1.0, delta=0.15)
                    self.assertAlmostEqual(slope_liquid / v_l, 1.0, delta=0.15)
                    self.assertGreater(slope_vapor, slope_liquid)   # konkaver Knick: g(T,p) = min(g_l, g_g)

    def test_equilibrium_surface_is_the_lower_envelope(self):
        """Die Gleichgewichtsflaeche ist die untere Huelle: g_stabil <= g jedes homogenen Zustands bei gleichem (T,p)."""
        for fluid in ("Water", "CO2"):
            d = data("gpt", fluid, pressure_axis="symlog")
            checked = 0
            for key, sl in (("metastable_liquid_surface", np.s_[:, 1:]), ("metastable_vapor_surface", np.s_[:, :-1])):
                for T, p, g in d[key][sl].reshape(-1, 3):     # Randspalte = Saettigungsknoten (dort g_stabil = g)
                    if not (np.isfinite([T, p, g]).all() and p > 0 and d["raster"]["p_min"] < p < d["raster"]["p_max"]):
                        continue
                    g_eq = CP.PropsSI("G", "T", T, "P", p, fluid)
                    self.assertGreaterEqual(g, g_eq - 1e-6 * abs(g_eq), (fluid, key, T, p))
                    checked += 1
            self.assertGreater(checked, 200)

    def test_metastable_surfaces_are_homogeneous_heos_values(self):
        """Metastabile Knoten: identisch zur direkten HEOS-Auswertung (T, v) mit vorgegebener Phase."""
        for fluid in ("Water", "CO2", "Nitrogen"):
            d = data("gpt", fluid)
            state = ref_state(fluid)
            for key, phase in (("metastable_liquid_surface", CP.iphase_liquid), ("metastable_vapor_surface", CP.iphase_gas),
                               ("unstable_surface", CP.iphase_gas)):
                arr, volumes = d[key], d["homogeneous_v"][key]
                finite = np.argwhere(np.isfinite(arr).all(axis=-1))
                for i, j in finite[::max(1, len(finite) // 25)]:
                    T, p, g = arr[i, j]
                    state.specify_phase(phase)
                    state.update(CP.DmolarT_INPUTS, 1.0 / volumes[i, j], T)
                    ref_p, ref_g = state.p(), state.gibbsmass()
                    state.unspecify_phase()
                    with self.subTest(fluid=fluid, key=key, node=(i, j)):
                        self.assertAlmostEqual(p, ref_p, delta=1e-7 * max(abs(ref_p), 1e3))
                        self.assertAlmostEqual(g, ref_g, delta=1e-7 * max(abs(ref_g), 1.0))

    def test_unstable_surface_pressure_rises_with_volume_and_metastable_falls(self):
        for fluid in ("Water", "CO2"):
            d = data("gpt", fluid)
            pc = d["critical_point"][1]
            for key, expected in (("metastable_liquid_surface", -1), ("metastable_vapor_surface", -1), ("unstable_surface", +1)):
                p, v = d[key][..., 1], d["homogeneous_v"][key]
                for i in range(p.shape[0]):
                    ok = np.isfinite(p[i]) & np.isfinite(v[i])
                    for j in np.flatnonzero(ok[:-1] & ok[1:]):
                        dp, dv = p[i, j + 1] - p[i, j], v[i, j + 1] - v[i, j]
                        with self.subTest(fluid=fluid, key=key, row=i, col=j):
                            self.assertTrue(expected * np.sign(dp) * np.sign(dv) >= 0 or abs(dp) <= 1e-9 * pc,
                                            (dp, dv))

    def test_spinodal_curves_vanish_dpdrho_and_agree_with_f_v_t(self):
        """Spinodalen: (dp/drho)_T = 0 (unabhaengig ausgewertet); g-p-T und f-v-T beschreiben denselben Zustand."""
        for fluid in ("Water", "CO2", "Nitrogen"):
            dg, df = data("gpt", fluid), data("fvt", fluid)
            state = ref_state(fluid)
            rhoc, pc = state.rhomolar_critical(), state.p_critical()
            for key, phase in (("spinodal_liquid", CP.iphase_liquid), ("spinodal_vapor", CP.iphase_gas)):
                cg, cf = dg[key][:-1], df[key][:-1]          # letzte Zeile = kritischer Punkt
                good = np.flatnonzero(np.isfinite(cg).all(axis=1) & np.isfinite(cf).all(axis=1))
                self.assertGreater(len(good), 20)
                for k in good[::5]:
                    T, p, g = cg[k]
                    _T, v, f = cf[k]
                    state.specify_phase(phase)
                    state.update(CP.DmolarT_INPUTS, 1.0 / v, T)
                    slope = state.first_partial_deriv(CP.iP, CP.iDmolar, CP.iT)
                    ref = (state.p(), state.gibbsmass(), state.helmholtzmass())
                    state.unspecify_phase()
                    with self.subTest(fluid=fluid, key=key, T=T):
                        self.assertLess(abs(slope) * rhoc / pc, 1e-6)
                        self.assertAlmostEqual(p / ref[0], 1.0, delta=1e-8)
                        self.assertAlmostEqual(g / ref[1], 1.0, delta=1e-8)
                        self.assertAlmostEqual(f / ref[2], 1.0, delta=1e-8)

    def test_critical_point_is_shared_by_all_curves(self):
        for fluid in FLUIDS:
            d = data("gpt", fluid)
            st = ref_state(fluid)
            st.update(CP.DmolarT_INPUTS, st.rhomolar_critical(), st.T_critical())
            Tc, pc, gc = d["critical_point"]
            self.assertAlmostEqual(gc / st.gibbsmass(), 1.0, delta=1e-9)
            self.assertAlmostEqual(pc / st.p_critical(), 1.0, delta=1e-9)
            for key in ("binodal_liquid", "binodal_vapor", "spinodal_liquid", "spinodal_vapor"):
                np.testing.assert_array_equal(d[key][-1], [Tc, pc, gc])

    def test_corrupted_data_is_rejected(self):
        d = data("gpt", "Water")
        bad = {**d, "surface_liquid": d["surface_liquid"].copy()}
        bad["surface_liquid"][5, 3, 2] = np.nan
        with self.assertRaises(ValueError):
            GPT.validate_gpt_data(bad)
        bad = {**d, "surface_vapor": d["surface_vapor"].copy()}
        bad["surface_vapor"][5, -1, 2] += 1.0          # Rand weicht von der Binodale ab
        with self.assertRaises(ValueError):
            GPT.validate_gpt_data(bad)
        bad = {**d, "surface_liquid": d["surface_liquid"].copy()}
        bad["surface_liquid"][5, 4, 2] = bad["surface_liquid"][5, 3, 2] - 1.0   # dg/dp = v > 0 verletzt
        with self.assertRaises(ValueError):
            GPT.validate_gpt_data(bad)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            GPT.berechne_gpt_daten("Water", pressure_axis="kubisch")
        with self.assertRaises(ValueError):
            GPT.berechne_gpt_daten("Water", n_p=1)
        with self.assertRaises(ValueError):
            GPT.berechne_gpt_daten("Water", p_min=5e8)
        with self.assertRaises(ValueError):
            GPT.berechne_gpt_daten("kein_stoff")
        self.assertTrue(np.isnan(GPT.g_eos(-1.0, 1e5, "Water")))

    def test_resolution_options(self):
        coarse = GPT.berechne_gpt_daten("Water", n_T_sub=12, n_T_super=6, n_p=9, include_homogeneous=False,
                                        include_isolines=False)
        fine = data("gpt", "Water")
        self.assertEqual(coarse["surface_liquid"].shape, (13, 9, 3))
        self.assertGreater(fine["surface_liquid"].size, coarse["surface_liquid"].size)
        self.assertNotIn("metastable_liquid_surface", coarse)
        self.assertEqual(coarse["isothermen"], [])

    def test_homogeneous_info_and_pressure_axis_notice(self):
        d = data("gpt", "Water")
        info = d["homogeneous_info"]
        self.assertEqual(info["mode"], "gpt")
        self.assertGreater(info["outside_pressure_window_nodes"], 0)
        self.assertIn("Druckfenster", d["pressure_axis_notice"])
        auto = data("gpt", "Water", auto_pressure_axis=True)
        self.assertEqual(auto["pressure_axis"], "symlog")


@unittest.skipIf(CP_SKIP, CP_SKIP or "")
class FvtThermodynamics(unittest.TestCase):
    """f-v-T: f = u - T s, Regelflaeche (gemeinsame Tangente), homogene Fortsetzung oberhalb der Tangente."""

    def test_surface_nodes_match_independent_coolprop(self):
        for fluid in FLUIDS:
            d = data("fvt", fluid)
            M = molar_mass(fluid)
            for name in ("surface_liquid", "surface_vapor", "surface_high_temperature", "two_phase_surface"):
                with self.subTest(fluid=fluid, surface=name):
                    for T, v, f in sample_nodes(d[name], cols=(0, 4, 19, 30, 39)):
                        rho = M / v                                         # kg/m3
                        _p, _g, f_ref, *_ = ref_gf_dmass(fluid, T, rho)
                        self.assertAlmostEqual(f / f_ref, 1.0, delta=1e-7, msg=f"T={T} v={v}")

    def test_g_equals_f_plus_p_v_with_mass_specific_quantities(self):
        for fluid in FLUIDS:
            d = data("fvt", fluid)
            M = molar_mass(fluid)
            for name in ("surface_liquid", "surface_vapor", "surface_high_temperature"):
                for T, v, f in sample_nodes(d[name], cols=(0, 4, 19, 30, 39)):
                    with self.subTest(fluid=fluid, surface=name, T=T):
                        rho = M / v
                        p, g, f_ref, *_ = ref_gf_dmass(fluid, T, rho)
                        self.assertAlmostEqual(g, f + p / rho, delta=1e-7 * max(abs(g), 1.0))   # g = f + p v

    def test_p_is_minus_df_dv_and_s_is_minus_df_dT(self):
        for fluid in FLUIDS:
            d = data("fvt", fluid)
            M = molar_mass(fluid)
            T, v, f = d["surface_high_temperature"][8, 30]
            rho = M / v

            def f_at(T_, rho_):
                return ref_gf_dmass(fluid, T_, rho_)[2]

            dv_mass = 1e-5 / rho
            df_dv = (f_at(T, 1.0 / (1.0 / rho + dv_mass)) - f_at(T, 1.0 / (1.0 / rho - dv_mass))) / (2 * dv_mass)
            dT = 1e-4 * T
            df_dT = (f_at(T + dT, rho) - f_at(T - dT, rho)) / (2 * dT)
            p, _g, _f, _h, s, _u = ref_gf_dmass(fluid, T, rho)
            self.assertAlmostEqual(-df_dv / p, 1.0, delta=1e-4)
            self.assertAlmostEqual(-df_dT / s, 1.0, delta=1e-4)

    def test_two_phase_surface_is_the_common_tangent(self):
        for fluid in FLUIDS:
            d = data("fvt", fluid)
            M = molar_mass(fluid)
            mix, Q = d["two_phase_surface"], d["Q_values"]
            for i in (3, 20, 40):
                T = mix[i, 0, 0]
                psat = CP.PropsSI("P", "T", T, "Q", 0, fluid)
                slopes = np.diff(mix[i, :, 2]) / np.diff(mix[i, :, 1])   # df/dv_molar = -p_sat / M
                with self.subTest(fluid=fluid, row=i):
                    np.testing.assert_allclose(slopes, -psat / M, rtol=1e-6)
                    for j in (1, len(Q) // 2, len(Q) - 2):
                        rho = M / mix[i, j, 1]
                        f_ref = ref_gf_dmass(fluid, T, rho)[2]          # unabhaengiger Gleichgewichtszustand (T, rho)
                        self.assertAlmostEqual(mix[i, j, 2] / f_ref, 1.0, delta=1e-7)
                    # g = f + p v ist ueber das ganze Zweiphasengebiet konstant (= g_sat)
                    g_sat = CP.PropsSI("G", "T", T, "Q", 0, fluid)
                    g_mix = mix[i, :, 2] + psat * mix[i, :, 1] / M
                    np.testing.assert_allclose(g_mix, g_sat, rtol=1e-7, atol=1e-5 * abs(g_sat))

    def test_homogeneous_surfaces_lie_above_the_common_tangent(self):
        for fluid in ("Water", "CO2", "Nitrogen"):
            d = data("fvt", fluid)
            l, g = d["binodal_liquid"], d["binodal_vapor"]
            checked = 0
            for key in ("metastable_liquid_surface", "metastable_vapor_surface", "unstable_surface"):
                arr = d[key]
                for i in range(arr.shape[0]):
                    for T, v, f in arr[i][np.isfinite(arr[i]).all(axis=-1)]:
                        if not (l[i, 1] < v < g[i, 1]):
                            continue
                        tangent = l[i, 2] + (v - l[i, 1]) / (g[i, 1] - l[i, 1]) * (g[i, 2] - l[i, 2])
                        self.assertGreaterEqual(f, tangent - 1e-9 * abs(tangent), (fluid, key, i, v))
                        checked += 1
            self.assertGreater(checked, 500)

    def test_metastable_and_unstable_curvature_classes(self):
        """Metastabil: (d2f/dv2)_T = -(dp/dv) > 0; instabil: < 0 (aus dp/drho der homogenen EOS, unabhaengig)."""
        for fluid in ("Water", "CO2", "Nitrogen"):
            d = data("fvt", fluid)
            state = ref_state(fluid)
            for key, phase, positive in (("metastable_liquid_surface", CP.iphase_liquid, True),
                                         ("metastable_vapor_surface", CP.iphase_gas, True),
                                         ("unstable_surface", CP.iphase_gas, False)):
                arr = d[key]
                finite = np.argwhere(np.isfinite(arr).all(axis=-1))
                for i, j in finite[1::max(1, len(finite) // 20)]:
                    T, v, f = arr[i, j]
                    if j in (0, arr.shape[1] - 1):      # Randknoten: Binodale (dp/drho > 0) bzw. Spinodale (dp/drho = 0)
                        continue
                    state.specify_phase(phase)
                    state.update(CP.DmolarT_INPUTS, 1.0 / v, T)
                    slope = state.first_partial_deriv(CP.iP, CP.iDmolar, CP.iT)
                    state.unspecify_phase()
                    with self.subTest(fluid=fluid, key=key, node=(i, j)):
                        self.assertEqual(slope > 0, positive)

    def test_homogeneous_f_is_direct_heos_value(self):
        for fluid in ("Water", "CO2"):
            d = data("fvt", fluid)
            state = ref_state(fluid)
            for key, phase in (("metastable_liquid_surface", CP.iphase_liquid), ("metastable_vapor_surface", CP.iphase_gas),
                               ("unstable_surface", CP.iphase_gas)):
                arr = d[key]
                finite = np.argwhere(np.isfinite(arr).all(axis=-1))
                for i, j in finite[::max(1, len(finite) // 25)]:
                    T, v, f = arr[i, j]
                    state.specify_phase(phase)
                    state.update(CP.DmolarT_INPUTS, 1.0 / v, T)
                    f_ref = state.umass() - T * state.smass()
                    state.unspecify_phase()
                    self.assertAlmostEqual(f / f_ref, 1.0, delta=1e-8)

    def test_critical_point_and_shared_curves(self):
        for fluid in FLUIDS:
            d = data("fvt", fluid)
            st = ref_state(fluid)
            st.update(CP.DmolarT_INPUTS, st.rhomolar_critical(), st.T_critical())
            Tc, vc, fc = d["critical_point"]
            self.assertAlmostEqual(fc / (st.umass() - Tc * st.smass()), 1.0, delta=1e-9)
            self.assertAlmostEqual(vc * st.rhomolar_critical(), 1.0, delta=1e-12)
            for key in ("binodal_liquid", "binodal_vapor", "spinodal_liquid", "spinodal_vapor"):
                np.testing.assert_array_equal(d[key][-1], [Tc, vc, fc])

    def test_corrupted_data_is_rejected(self):
        d = data("fvt", "Water")
        for key, edit in (("surface_liquid", lambda a: a.__setitem__((5, 3, 2), np.nan)),
                          ("two_phase_surface", lambda a: a.__setitem__((5, 0, 2), a[5, 0, 2] + 1.0)),
                          ("binodal_vapor", lambda a: a.__setitem__((5, 1), a[5, 1] * 0.5))):
            bad = {**d, key: d[key].copy()}
            edit(bad[key])
            with self.subTest(key=key), self.assertRaises(ValueError):
                FVT.validate_fvt_data(bad)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            FVT.berechne_fvt_daten("Water", n_v=1)
        with self.assertRaises(ValueError):
            FVT.berechne_fvt_daten("Water", v_min=1.0, v_max=2.0)
        with self.assertRaises(ValueError):
            FVT.berechne_fvt_daten("kein_stoff")
        self.assertTrue(np.isnan(FVT.f_eos(-5.0, 1e-3, "Water")))

    def test_resolution_options(self):
        coarse = FVT.berechne_fvt_daten("Water", n_T_sub=12, n_T_super=6, n_v=9, n_quality=5,
                                        include_homogeneous=False, include_isolines=False)
        self.assertEqual(coarse["surface_liquid"].shape, (13, 9, 3))
        self.assertEqual(coarse["two_phase_surface"].shape, (13, 5, 3))
        self.assertNotIn("unstable_surface", coarse)


@unittest.skipIf(CP_SKIP, CP_SKIP or "")
class RegistryAndAdapters(unittest.TestCase):
    def test_registry_and_dispatch(self):
        self.assertEqual(DIAGRAMS.ids()[-2:], ["gpt", "fvt"])
        self.assertEqual(DIAGRAMS.get("gpt").label, "g-p-T (3D)")
        self.assertEqual(DIAGRAMS.get("fvt").label, "f-v-T (3D)")
        self.assertNotIn("zweiphasengebiet", DIAGRAMS.get("gpt").layers)     # Zweiphasengebiet = Sattigungslinie
        self.assertIn("zweiphasengebiet", DIAGRAMS.get("fvt").layers)
        for spec in (DIAGRAMS.get("gpt"), DIAGRAMS.get("fvt")):
            for layer in ("spinodale", "metastable_liquid_surface", "metastable_vapor_surface", "unstable_surface"):
                self.assertIn(layer, spec.layers)
                self.assertFalse(spec.initial_layer_state()[layer])            # Gleichgewichtsflaeche ohne Vermischung
            self.assertTrue(spec.initial_layer_state()["surface"])
            self.assertNotIn("replacement_surface", spec.layers)                # pvt-spezifisch
        self.assertEqual(set(app_main.BERECHNUNGEN), set(DIAGRAMS.ids()))

    def test_berechne_auswahl_returns_registered_diagram_data(self):
        _, raster, visu = app_main.lade_konfiguration()
        for diagram in ("gpt", "fvt"):
            daten = app_main.berechne_auswahl("Water", raster, visu, diagram_type=diagram)
            for key in ("surface_liquid", "surface_vapor", "surface_high_temperature", "binodal_liquid", "binodal_vapor",
                        "critical_point", "raster", "meta", "isothermen", "isobaren", "isochoren"):
                self.assertIn(key, daten, (diagram, key))
        with self.assertRaises(ValueError):
            app_main.berechne_auswahl("Water", raster, visu, diagram_type="xyz")

    def test_cache_parameters_follow_settings_and_q_points(self):
        _, raster, visu = app_main.lade_konfiguration()
        for diagram in ("gpt", "fvt"):
            params = DIAGRAMS.get(diagram).cache_params(raster, visu)
            self.assertIn("Q_points", params)
            self.assertEqual(params["diagram_settings"], app_main.lade_diagramm_optionen(diagram))
        self.assertEqual(app_main.lade_diagramm_optionen("gpt")["pressure_axis"], "log")

    def test_coordinate_roundtrip(self):
        for fluid in ("Water", "CO2"):
            d_g, d_f = data("gpt", fluid), data("fvt", fluid)
            for axis in ("log", "asinh", "symlog", "linear"):
                pts = d_g["surface_liquid"].reshape(-1, 3)
                xyz = dc.to_xyz_gpt(pts, d_g["critical_point"], axis)
                np.testing.assert_allclose(dc.from_xyz_gpt(xyz, d_g["critical_point"], axis), pts, rtol=1e-9)
            pts = d_f["surface_vapor"].reshape(-1, 3)
            np.testing.assert_allclose(dc.from_xyz_fvt(dc.to_xyz_fvt(pts, d_f["critical_point"]), d_f["critical_point"]),
                                       pts, rtol=1e-9)
            # Achsenzuordnung: x = log10(p/pc), y = T/Tc, z = (g-gc)/1e6
            Tc, pc, gc = d_g["critical_point"]
            np.testing.assert_allclose(dc.to_xyz_gpt([[Tc, pc, gc]], d_g["critical_point"]), [[0.0, 1.0, 0.0]], atol=1e-12)
            Tc, vc, fc = d_f["critical_point"]
            np.testing.assert_allclose(dc.to_xyz_fvt([[Tc, vc, fc]], d_f["critical_point"]), [[0.0, 1.0, 0.0]], atol=1e-12)

    def test_pressure_window_hides_but_does_not_change_data(self):
        d = data("gpt", "Water")
        arr = d["metastable_liquid_surface"]
        window = (d["raster"]["p_min"], d["raster"]["p_max"])
        shown = dc.to_xyz_gpt(arr, d["critical_point"], "log", window)
        everything = dc.to_xyz_gpt(arr, d["critical_point"], "symlog")
        self.assertGreater(np.isnan(shown).any(axis=-1).sum(), 0)
        self.assertTrue(np.isfinite(everything[np.isfinite(arr).all(axis=-1)]).all())
        stable = dc.to_xyz_gpt(d["surface_liquid"], d["critical_point"], "log", window)
        self.assertTrue(np.isfinite(stable).all())                                # Gleichgewichtsflaeche nie abgeschnitten

    def test_every_diagram_has_an_adapter_and_unique_node_references(self):
        from state_store import ADAPTERS
        self.assertEqual(set(ADAPTERS), set(DIAGRAMS.ids()))
        for fluid in ("Water", "CO2"):
            for diagram in ("gpt", "fvt"):
                store, index, _ = tsc.make(fluid, diagram)
                refs = list(index.iter_refs())
                self.assertEqual(len(refs), len(set(refs)))
                self.assertGreater(len(refs), 5000)
                # Rasterarrays werden referenziert, nicht kopiert
                self.assertTrue(np.shares_memory(index.surfaces["surface_liquid"].array, index.daten["surface_liquid"]))

    def test_fluid_matrix_for_ten_gui_fluids(self):
        _, raster, visu = app_main.lade_konfiguration()
        fluids = ["Water", "CO2", "Methane", "n-Propane", "Argon", "Nitrogen", "Ammonia", "Ethanol", "R134a", "n-Hexane"]
        for fluid in fluids:
            r, v, _meta = app_main.generate_stoff_grid(fluid, raster, visu)
            for diagram in ("gpt", "fvt"):
                with self.subTest(fluid=fluid, diagram=diagram):
                    daten = app_main.berechne_auswahl(fluid, r, v, diagram_type=diagram)
                    self.assertTrue(np.isfinite(daten["surface_liquid"]).all())
                    self.assertTrue(np.isfinite(daten["surface_high_temperature"]).all())


@unittest.skipIf(CP_SKIP, CP_SKIP or "")
class PointPickerStates(unittest.TestCase):
    """Zustandsrekonstruktion fuer Punkte auf g-p-T- und f-v-T-Flaechen (gegen direkte HEOS-Aufrufe)."""

    def test_exact_nodes_reproduce_source_data_in_correct_units(self):
        cases = (("gpt", "surface_liquid", (12, 8)), ("gpt", "surface_vapor", (12, 8)),
                 ("gpt", "surface_high_temperature", (5, 20)), ("fvt", "surface_liquid", (12, 8)),
                 ("fvt", "surface_vapor", (12, 8)), ("fvt", "surface_high_temperature", (5, 20)),
                 ("fvt", "two_phase_surface", (12, 20)))
        for fluid in ("Water", "CO2"):
            M = molar_mass(fluid)
            for diagram, surface, node in cases:
                with self.subTest(fluid=fluid, diagram=diagram, surface=surface):
                    store, index, resolver = tsc.make(fluid, diagram)
                    state = resolver.resolve_node(index, NodeRef(surface, node))
                    T, a, b = index.native(NodeRef(surface, node))
                    self.assertEqual((state.fluid, state.source_diagram, state.eos_source), (fluid, diagram, ts.SOURCE_HEOS))
                    self.assertIn(state.validity, (ts.VALID, ts.PARTIALLY_VALID), state.notes)
                    self.assertAlmostEqual(state.T / T, 1.0, delta=1e-6)
                    if diagram == "gpt":
                        self.assertAlmostEqual(state.p / a, 1.0, delta=1e-6)
                        self.assertAlmostEqual(state.g_molar / (b * M), 1.0, delta=1e-6)       # J/kg -> J/mol
                    else:
                        self.assertAlmostEqual(state.v_molar / a, 1.0, delta=1e-6)
                        self.assertAlmostEqual(state.f_molar / (b * M), 1.0, delta=1e-6)
                    # g = f + p v aus den molaren Zustandsgroessen desselben Zustands
                    self.assertAlmostEqual(state.g_molar, state.f_molar + state.p * state.v_molar,
                                           delta=1e-8 * max(abs(state.g_molar), 1.0))
                    # g = h - T s, f = u - T s
                    self.assertAlmostEqual(state.g_molar, state.h_molar - state.T * state.s_molar,
                                           delta=1e-8 * max(abs(state.g_molar), 1.0))
                    self.assertAlmostEqual(state.f_molar, state.u_molar - state.T * state.s_molar,
                                           delta=1e-8 * max(abs(state.f_molar), 1.0))

    def test_interior_points_match_direct_heos(self):
        for fluid in ("Water", "CO2", "Nitrogen"):
            M = molar_mass(fluid)
            for surface in ("surface_liquid", "surface_vapor", "surface_high_temperature"):
                store, index, resolver = tsc.make(fluid, "gpt")
                pick = tsc.cell(index, surface, [(6, 7), (6, 8), (7, 7)], (0.3, 0.3, 0.4))
                state = resolver.resolve(pick, index)
                T, p, _g = pick.physical_hint["native_interpolated"]
                with self.subTest(fluid=fluid, diagram="gpt", surface=surface):
                    self.assertIn(state.validity, (ts.VALID, ts.PARTIALLY_VALID), state.notes)
                    self.assertAlmostEqual(state.g_molar / (M * CP.PropsSI("G", "T", T, "P", p, fluid)), 1.0, delta=1e-8)
                    self.assertAlmostEqual(state.T / T, 1.0, delta=1e-9)
                    self.assertAlmostEqual(state.p / p, 1.0, delta=1e-9)
                store, index, resolver = tsc.make(fluid, "fvt")
                pick = tsc.cell(index, surface, [(6, 7), (6, 8), (7, 7)], (0.3, 0.3, 0.4))
                state = resolver.resolve(pick, index)
                T, v, _f = pick.physical_hint["native_interpolated"]
                with self.subTest(fluid=fluid, diagram="fvt", surface=surface):
                    f_ref = ref_gf_dmass(fluid, T, M / v)[2]
                    self.assertIn(state.validity, (ts.VALID, ts.PARTIALLY_VALID), state.notes)
                    self.assertAlmostEqual(state.f_molar / (M * f_ref), 1.0, delta=1e-8)

    def test_gpt_saturation_edge_nodes_are_one_sided_states_without_q(self):
        for fluid in ("Water", "CO2"):
            store, index, resolver = tsc.make(fluid, "gpt")
            M = molar_mass(fluid)
            for surface, node, phase in (("surface_liquid", (15, 0), "liquid"), ("surface_vapor", (15, -1), "gas")):
                ref = NodeRef(surface, (node[0], node[1] % index.surfaces[surface].array.shape[1]))
                state = resolver.resolve_node(index, ref)
                T, p, g = index.native(ref)
                with self.subTest(fluid=fluid, surface=surface):
                    self.assertEqual(state.phase, phase)
                    self.assertIsNone(state.Q)
                    self.assertEqual(state.stability, "saturated")
                    self.assertTrue(any("Saettigungslinie" in n for n in state.notes))
                    self.assertAlmostEqual(state.p / CP.PropsSI("P", "T", T, "Q", 0, fluid), 1.0, delta=1e-6)
                    self.assertAlmostEqual(state.g_molar / (g * M), 1.0, delta=1e-7)

    def test_gpt_equilibrium_point_across_saturation_gets_stable_phase_and_flag(self):
        """Zwischen zwei Sattigungsknoten liegt der interpolierte Punkt (konkav in log p) unter p_sat: Gleichgewicht = Dampf."""
        store, index, resolver = tsc.make("Water", "gpt")
        pick = tsc.cell(index, "surface_liquid", [(10, 0), (11, 0), (10, 1)], (0.5, 0.5, 0.0))
        state = resolver.resolve(pick, index)
        T, p, _ = pick.physical_hint["native_interpolated"]
        self.assertLess(p, CP.PropsSI("P", "T", T, "Q", 0, "Water"))
        self.assertEqual(state.validity, ts.PARTIALLY_VALID)
        self.assertIn(state.phase, ("gas", "twophase", "supercritical_gas"))
        self.assertTrue(any("Phasenrandnaehe" in n for n in state.notes), state.notes)
        self.assertAlmostEqual(state.g_molar / (molar_mass("Water") * CP.PropsSI("G", "T", T, "P", p, "Water")), 1.0,
                               delta=1e-7)

    def test_gpt_homogeneous_points_use_stored_volume_and_are_marked(self):
        for fluid in ("Water", "CO2"):
            M = molar_mass(fluid)
            store, index, resolver = tsc.make(fluid, "gpt")
            direct = ref_state(fluid)
            for surface, stability, phase in (("metastable_vapor_surface", "metastable", CP.iphase_gas),
                                              ("metastable_liquid_surface", "metastable", CP.iphase_liquid),
                                              ("unstable_surface", "unstable", CP.iphase_gas)):
                s = index.surface(surface)
                nodes = None
                for i in range(10, 30):
                    for j in range(5, 12):
                        cand = [(i, j), (i + 1, j), (i + 1, j + 1)]
                        xyz = index.adapter.to_display(np.array([s.node(n) for n in cand]), index.daten)
                        if np.isfinite(xyz).all():
                            nodes = cand
                            break
                    if nodes:
                        break
                self.assertIsNotNone(nodes, surface)
                w = (0.3, 0.3, 0.4)
                pick = tsc.cell(index, surface, nodes, w)
                state = resolver.resolve(pick, index)
                v = float(np.dot(w, [s.node_aux("v", n) for n in nodes]))
                T = pick.physical_hint["native_interpolated"][0]
                direct.specify_phase(phase)
                direct.update(CP.DmolarT_INPUTS, 1.0 / v, T)
                with self.subTest(fluid=fluid, surface=surface):
                    self.assertEqual(state.eos_source, ts.SOURCE_HEOS_HOMOGENEOUS)
                    self.assertEqual(state.stability, stability)
                    self.assertIn(state.validity, (ts.VALID, ts.PARTIALLY_VALID))
                    self.assertAlmostEqual(state.p / direct.p(), 1.0, delta=1e-8)
                    self.assertAlmostEqual(state.g_molar / (M * direct.gibbsmass()), 1.0, delta=1e-8)
                    self.assertAlmostEqual(state.f_molar / (M * direct.helmholtzmass()), 1.0, delta=1e-8)
                direct.unspecify_phase()

    def test_fvt_homogeneous_and_two_phase_points(self):
        for fluid in ("Water", "CO2"):
            M = molar_mass(fluid)
            store, index, resolver = tsc.make(fluid, "fvt")
            two = index.surface("two_phase_surface")
            for j in (0, 20, 39):
                state = resolver.resolve_node(index, NodeRef("two_phase_surface", (12, j)))
                T, v, f = index.native(NodeRef("two_phase_surface", (12, j)))
                with self.subTest(fluid=fluid, Q_index=j):
                    self.assertAlmostEqual(state.Q, two.q_values[j], delta=1e-7)
                    self.assertAlmostEqual(state.f_molar / (f * M), 1.0, delta=1e-7)
                    self.assertAlmostEqual(state.v_molar / v, 1.0, delta=1e-7)
            for surface in ("metastable_vapor_surface", "unstable_surface"):
                pick = tsc.cell(index, surface, [(15, 6), (16, 6), (16, 7)], (0.3, 0.3, 0.4))
                state = resolver.resolve(pick, index)
                T, v, _f = pick.physical_hint["native_interpolated"]
                phase = CP.iphase_gas
                direct = ref_state(fluid)
                direct.specify_phase(phase)
                direct.update(CP.DmolarT_INPUTS, 1.0 / v, T)
                with self.subTest(fluid=fluid, surface=surface):
                    self.assertEqual(state.eos_source, ts.SOURCE_HEOS_HOMOGENEOUS)
                    self.assertAlmostEqual(state.f_molar / (M * direct.helmholtzmass()), 1.0, delta=1e-8)

    def test_critical_node_is_valid_for_both_diagrams(self):
        for diagram in ("gpt", "fvt"):
            store, index, resolver = tsc.make("Water", diagram)
            state = resolver.resolve_node(index, NodeRef("critical_point", (0,)))
            self.assertEqual(state.validity, ts.VALID)
            self.assertEqual(state.phase, "critical_point")

    def test_comparison_view_accepts_g_and_f_axes_without_recomputation(self):
        manager = SelectionManager()
        states = []
        for diagram, surface in (("gpt", "surface_vapor"), ("fvt", "surface_liquid"), ("gpt", "surface_high_temperature")):
            store, index, resolver = tsc.make("Water", diagram)
            pick = tsc.cell(index, surface, [(6, 7), (6, 8), (7, 7)], (0.3, 0.3, 0.4))
            states.append(resolver.resolve(pick, index))
            manager.add(states[-1], pick.display_xyz, {"T": 1.0})
        evaluations = sum(1 for _ in states)
        geometry = build_geometry(manager.points, ("g", "f", "T"))
        self.assertEqual(len(geometry.valid_rows), 3)
        self.assertEqual({r.source_diagram for r in geometry.valid_rows}, {"gpt", "fvt"})
        geometry2 = build_geometry(manager.points, ("T", "g", "f"))
        self.assertEqual(len(geometry2.valid_rows), evaluations)


# ---------------------------------------------------------------- GUI (nativer OpenGL-Kontext)

if not GUI_SKIP:
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QApplication
    from state_picker import PickError, pick_from_viewport


@unittest.skipIf(GUI_SKIP or CP_SKIP, GUI_SKIP or CP_SKIP or "")
class GptFvtGui(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def make_window(self):
        stoff, raster, visu = app_main.lade_konfiguration()
        self.calls = []

        def berechnen(stoff_, raster_, visu_, diagram_type="pvt", progress_callback=None):
            self.calls.append(diagram_type)
            return app_main.berechne_auswahl(stoff_, raster_, visu_, diagram_type=diagram_type,
                                             progress_callback=progress_callback)

        window = app_main.ThermoSurfaceWindow(
            start_stoff=stoff, surface_cfg=raster["isotherm_grid"], raster=raster, visu_raster=visu,
            berechnen=berechnen, grid_generator=app_main.generate_stoff_grid)
        window.show()
        self.addCleanup(self.close_window, window)
        return window

    def close_window(self, window):
        window.manager.shutdown()
        window.close()
        self.app.processEvents()

    def calculate(self, window):
        window.on_calculate_clicked()
        spin(self.app, lambda: not window.manager.is_running)
        self.app.processEvents()

    def pump(self, seconds):
        import time
        end = time.time() + seconds
        while time.time() < end:
            self.app.processEvents()
            time.sleep(0.01)

    def find_hit(self, window, viewport, preferred=("surface",)):
        key = viewport.result_key
        index = window.state_store.ensure(key.fluid, viewport.diagram_id, viewport.daten, key)
        width, height = viewport.plotter.width(), viewport.plotter.height()
        ratio = viewport.plotter.devicePixelRatioF()
        for fy in np.linspace(0.1, 0.9, 41):
            for fx in np.linspace(0.05, 0.95, 46):
                qx, qy = int(width * fx), int(height * fy)
                try:
                    pick = pick_from_viewport(viewport, index, int(round(qx * ratio)), int(round((height - 1 - qy) * ratio)))
                except PickError:
                    continue
                if pick.actor_key in preferred:
                    return qx, qy
        self.fail("Kein pickbarer Pixel gefunden")

    def test_diagram_selection_contains_all_seven_types_in_order(self):
        window = self.make_window()
        ids = [window.diagram_combo.itemData(i) for i in range(window.diagram_combo.count())]
        self.assertEqual(ids, ["pvt", "pht", "pst", "hsp", "uvt", "gpt", "fvt"])
        texts = [window.diagram_combo.itemText(i) for i in range(window.diagram_combo.count())]
        self.assertIn("g-p-T (3D)", texts)
        self.assertIn("f-v-T (3D)", texts)
        self.assertEqual(sorted(window.diagram_checks), sorted(ids))

    def test_single_mode_layers_toggle_and_diagram_switch_without_restart(self):
        window = self.make_window()
        window.on_confirm_stoff()
        for diagram in ("gpt", "fvt", "pht", "gpt"):
            window.diagram_combo.setCurrentIndex(window.diagram_combo.findData(diagram))
            window.on_confirm_stoff()
            self.calculate(window)
            viewport = window.workspace.viewports[diagram]
            self.assertTrue(viewport.has_result)
            layers = viewport.available_layers()
            self.assertTrue(layers["surface"] and layers["critical_point"] and layers["binodale"])
            if diagram in ("gpt", "fvt"):
                for key in ("spinodale", "metastable_liquid_surface", "metastable_vapor_surface", "unstable_surface",
                            "isothermen", "isobaren", "isochoren"):
                    self.assertTrue(layers[key], (diagram, key))
                    self.assertFalse(window.layer_checks[key].isChecked() and key not in ("critical_point",))
                self.assertEqual(layers["zweiphasengebiet"], diagram == "fvt")
                # Ebene ein-/ausblenden wirkt nur auf diesen Viewport
                box = window.layer_checks["metastable_vapor_surface"]
                box.setChecked(True)
                actor = viewport.actors["metastable_vapor_surface"]
                self.assertTrue(actor.visibility)
                box.setChecked(False)
                self.assertFalse(actor.visibility)
                window.layer_checks["spinodale"].setChecked(True)
                self.assertTrue(all(a.visibility for a in viewport.actors["spinodale"]))
                window.layer_checks["spinodale"].setChecked(False)
                window.layer_checks["surface"].setChecked(False)
                self.assertTrue(all(not a.visibility for a in viewport.actors["surface"]))
                window.layer_checks["surface"].setChecked(True)
        self.assertEqual(self.calls, ["gpt", "fvt", "pht"])        # Wiederholung: Cache-Treffer

    def test_multi_view_with_all_seven_diagrams_and_independent_views(self):
        window = self.make_window()
        window.on_confirm_stoff()
        window.mode_combo.setCurrentIndex(window.mode_combo.findData("multi"))
        for diagram_id, box in window.diagram_checks.items():
            box.setChecked(diagram_id in ("pvt", "gpt", "fvt", "uvt"))
        self.calculate(window)
        self.assertEqual(sorted(self.calls), sorted(["pvt", "gpt", "fvt", "uvt"]))
        views = window.workspace.viewports
        self.assertEqual(set(views), {"pvt", "gpt", "fvt", "uvt"})
        self.assertTrue(all(v.has_result for v in views.values()))
        views["gpt"].set_layer_state({"metastable_liquid_surface": True})
        views["gpt"].apply_layers()
        self.assertTrue(views["gpt"].actors["metastable_liquid_surface"].visibility)
        self.assertFalse(views["fvt"].actors["metastable_liquid_surface"].visibility)
        window.workspace.maximize(views["gpt"])
        window.workspace.restore()
        self.assertEqual(len(self.calls), 4)          # reiner Ansichtswechsel: keine Neuberechnung
        for box in window.diagram_checks.values():
            box.setChecked(True)
        self.calculate(window)
        self.assertEqual(len(window.workspace.viewports), 7)
        self.assertEqual(len(self.calls), 7)          # nur die drei neuen Diagramme wurden berechnet

    def test_point_picker_on_gpt_and_fvt_and_comparison(self):
        window = self.make_window()
        window.on_confirm_stoff()
        window.mode_combo.setCurrentIndex(window.mode_combo.findData("multi"))
        for diagram_id, box in window.diagram_checks.items():
            box.setChecked(diagram_id in ("gpt", "fvt"))
        self.calculate(window)
        window.compare_check.setChecked(True)
        picks = {}
        for diagram in ("gpt", "fvt"):
            viewport = window.workspace.viewports[diagram]
            qx, qy = self.find_hit(window, viewport)
            QTest.mouseClick(viewport.plotter, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, QPoint(qx, qy))
            self.pump(0.8)
            picks[diagram] = viewport
        self.assertEqual({p.source_diagram for p in window.selection.points}, {"gpt", "fvt"})
        self.assertEqual(picks["gpt"].marker_count, 1)
        self.assertEqual(picks["fvt"].marker_count, 1)
        for p in window.selection.points:
            self.assertIn(p.state.validity, (ts.VALID, ts.PARTIALLY_VALID))
            self.assertIsNotNone(p.state.g_molar)
            self.assertIsNotNone(p.state.f_molar)
            self.assertAlmostEqual(p.state.g_molar, p.state.f_molar + p.state.p * p.state.v_molar,
                                   delta=1e-8 * abs(p.state.g_molar))
        window.on_open_comparison()
        window.comparison_view.set_axes("g", "f", "T")
        self.assertEqual(len([r for r in window.comparison_view.geometry.rows if r.xyz is not None]), 2)
        # Doppelklick maximiert, waehlt aber keinen Punkt
        viewport = picks["gpt"]
        QTest.mouseDClick(viewport.plotter, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, QPoint(qx, qy))
        self.pump(0.8)
        self.assertEqual(len(window.selection), 2)

    def test_picker_inactive_without_compare_mode_and_unpickable_layers(self):
        window = self.make_window()
        window.on_confirm_stoff()
        window.diagram_combo.setCurrentIndex(window.diagram_combo.findData("gpt"))
        window.on_confirm_stoff()
        self.calculate(window)
        viewport = window.workspace.viewports["gpt"]
        qx, qy = self.find_hit(window, viewport)
        before = window.state_store.evaluations
        QTest.mouseClick(viewport.plotter, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, QPoint(qx, qy))
        self.pump(0.8)
        self.assertEqual((len(window.selection), window.state_store.evaluations), (0, before))

    def test_homogeneous_layer_points_can_be_picked(self):
        window = self.make_window()
        window.on_confirm_stoff()
        for diagram in ("fvt", "gpt"):
            window.diagram_combo.setCurrentIndex(window.diagram_combo.findData(diagram))
            window.on_confirm_stoff()
            self.calculate(window)
            viewport = window.workspace.viewports[diagram]
            for key in ("surface", "metastable_vapor_surface", "unstable_surface"):
                window.layer_checks[key].setChecked(key != "surface")
            self.pump(0.3)
            window.compare_check.setChecked(True)
            qx, qy = self.find_hit(window, viewport, preferred=("metastable_vapor_surface", "unstable_surface"))
            before = len(window.selection)
            QTest.mouseClick(viewport.plotter, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, QPoint(qx, qy))
            self.pump(0.8)
            self.assertEqual(len(window.selection), before + 1, diagram)
            point = window.selection.points[-1]
            self.assertEqual(point.state.eos_source, ts.SOURCE_HEOS_HOMOGENEOUS)
            self.assertIn(point.state.stability, ("metastable", "unstable"))
            window.on_compare_reset()
            window.compare_check.setChecked(False)
            window.layer_checks["surface"].setChecked(True)


    def test_renderers_for_all_pressure_axes_and_minimal_options(self):
        import pyvista as pv
        import visualization as V
        plotter = pv.Plotter(off_screen=True)
        self.addCleanup(plotter.close)
        for axis in ("log", "asinh", "symlog", "linear"):
            actors = V.render_gpt(plotter, GPT.berechne_gpt_daten("Water", pressure_axis=axis))
            self.assertEqual(len(actors["surface"]), 3, axis)
            self.assertIsNone(actors["zweiphasengebiet"])                  # g-p-T: keine Zweiphasenflaeche
            self.assertEqual(len(actors["binodale"]), 1)                   # Siede- und Taulinie fallen zusammen
            self.assertIsNotNone(actors["unstable_surface"], axis)
        minimal = V.render_gpt(plotter, GPT.berechne_gpt_daten("CO2", include_homogeneous=False, include_isolines=False))
        self.assertNotIn("metastable_liquid_surface", minimal)
        self.assertEqual(minimal["isothermen"], [])
        minimal = V.render_fvt(plotter, FVT.berechne_fvt_daten("CO2", include_homogeneous=False, include_isolines=False))
        self.assertIsNotNone(minimal["zweiphasengebiet"])
        self.assertEqual(len(minimal["binodale"]), 2)


if __name__ == "__main__":
    unittest.main()
