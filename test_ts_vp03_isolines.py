"""TS-VP03/TS-VP07: geometrisch verteilte Isolinien in allen sieben Berechnungspfaden (Testmatrix L04, L05, C02).

Echte HEOS-Laeufe (CoolProp); ohne CoolProp/PyVista werden die Tests uebersprungen.
"""
import os
import unittest

import numpy as np

os.environ.setdefault("QT_API", "pyside6")
SKIP = None
try:
    import CoolProp.CoolProp as CP
    from CoolProp import AbstractState
    import main as M
    from calculation_manager import DiagramTask, cache_key_for
    from diagram_registry import DIAGRAMS
    from isoline_levels import generate_levels
except Exception as exc:  # pragma: no cover
    SKIP = f"CoolProp/main nicht verfuegbar: {exc}"

FLUIDS = ("Water", "CO2", "R134a", "Argon")
FAMILIES = ("isothermen", "isochoren", "isobaren")
_CACHE = {}


def _grids(fluid):
    _, raster, visu = M.lade_konfiguration()
    return M.generate_stoff_grid(fluid, raster, visu)[:2]


def _critical(fluid):
    state = AbstractState("HEOS", fluid)
    return state.T_critical(), state.p_critical(), 1.0 / state.rhomolar_critical()


def _config(fluid, factor_T=1.06, factor_p=1.5, factor_v=1.6):
    Tc, pc, vc = _critical(fluid)
    return {"isothermen": {"start": 0.8 * Tc, "stop": 1.1 * Tc, "factor": factor_T},
            "isobaren": {"start": 0.2 * pc, "stop": 1.5 * pc, "factor": factor_p},
            "isochoren": {"start": 0.7 * vc, "stop": 5.0 * vc, "factor": factor_v}}


def compute(fluid, diagram_id, isolines=None):
    key = (fluid, diagram_id, repr(isolines))
    if key not in _CACHE:
        raster, visu = _grids(fluid)
        _CACHE[key] = M.berechne_auswahl(fluid, raster, visu, diagram_type=diagram_id, isolines=isolines)
    return _CACHE[key]


def finite_rows(curve):
    curve = np.asarray(curve, dtype=float)
    return curve[np.isfinite(curve).all(axis=1)]


def state_of(diagram_id, fluid, row):
    """HEOS-Zustand zu einem nativen Datenpunkt dieses Diagramms (oder None, wenn nicht rekonstruierbar)."""
    state = AbstractState("HEOS", fluid)
    try:
        if diagram_id == "pvt" or diagram_id in ("uvt", "fvt"):
            state.update(CP.DmolarT_INPUTS, 1.0 / row[1], row[0])
        elif diagram_id == "pht":
            state.update(CP.HmassP_INPUTS, row[1], row[2])
        elif diagram_id == "pst":
            state.update(CP.PSmass_INPUTS, row[2], row[1])
        elif diagram_id == "hsp":
            state.update(CP.HmassP_INPUTS, row[0], row[2])
        elif diagram_id == "gpt":
            state.update(CP.PT_INPUTS, row[1], row[0])
    except (ValueError, RuntimeError):
        return None
    return state


def constant_values(diagram_id, fluid, family, curve):
    """Werte der konstant gehaltenen Groesse entlang der Kurve (T, v oder p), aus Daten bzw. HEOS rekonstruiert."""
    rows = finite_rows(curve)
    columns = {"pvt": {"isothermen": 0, "isochoren": 1, "isobaren": 2}, "pht": {"isothermen": 0, "isobaren": 2},
               "pst": {"isothermen": 0, "isobaren": 2}, "hsp": {"isobaren": 2},
               "uvt": {"isothermen": 0, "isochoren": 1}, "gpt": {"isothermen": 0, "isobaren": 1},
               "fvt": {"isothermen": 0, "isochoren": 1}}[diagram_id]
    if family in columns:
        return rows[:, columns[family]], len(rows)
    values = []
    for row in rows:
        state = state_of(diagram_id, fluid, row)
        if state is None:
            continue
        values.append({"isothermen": lambda: state.T(), "isochoren": lambda: 1.0 / state.rhomolar(),
                       "isobaren": lambda: state.p()}[family]())
    return np.array(values), len(rows)


@unittest.skipIf(SKIP, SKIP or "")
class LevelsAcrossDiagrams(unittest.TestCase):
    def test_same_geometric_levels_in_all_seven_diagrams(self):
        for fluid in ("Water", "CO2"):
            Tc, pc, vc = _critical(fluid)
            cfg = _config(fluid)
            expected = {fam: generate_levels(cfg[fam]["start"], cfg[fam]["stop"], cfg[fam]["factor"]) for fam in FAMILIES}
            for diagram_id in DIAGRAMS.ids():
                with self.subTest(fluid=fluid, diagram=diagram_id):
                    daten = compute(fluid, diagram_id, cfg)
                    report = daten["isoline_report"]
                    for family in FAMILIES:
                        reported = np.array(report[family]["levels"])
                        # innerhalb des Diagrammbereichs exakt die Formel x_n = start * f**n (zusammenhaengende n)
                        self.assertTrue(len(reported) >= 1, (family, report[family]))
                        idx = [int(np.argmin(np.abs(expected[family] / r - 1.0))) for r in reported]
                        np.testing.assert_allclose(reported, expected[family][idx], rtol=1e-12)
                        self.assertEqual(idx, list(range(idx[0], idx[0] + len(idx))))
                        self.assertEqual(report[family]["n_levels"], len(reported))
                        # Levelanzahl und Abtastpunktzahl werden getrennt protokolliert
                        self.assertEqual(len(report[family]["samples_per_curve"]), report[family]["n_curves"])
                        self.assertTrue(all(s >= 2 for s in report[family]["samples_per_curve"]))

    def test_domain_filtering_reports_dropped_levels(self):
        daten = compute("Water", "gpt", {"isobaren": {"start": 1.0, "stop": 1e12, "factor": 10.0}})
        entry = daten["isoline_report"]["isobaren"]
        self.assertTrue(entry["dropped_outside_domain"])
        lo, hi = entry["domain"]
        self.assertTrue(all(lo * (1 - 1e-9) <= x <= hi * (1 + 1e-9) for x in entry["levels"]))


@unittest.skipIf(SKIP, SKIP or "")
class ConstantAlongCurve(unittest.TestCase):
    """L04: Residuen der konstant gehaltenen Zustandsgroesse entlang jeder Kurve."""

    def test_T_v_p_are_constant_on_every_curve(self):
        for fluid in ("Water", "CO2"):
            cfg = _config(fluid)
            for diagram_id in DIAGRAMS.ids():
                daten = compute(fluid, diagram_id, cfg)
                report = daten["isoline_report"]
                for family in FAMILIES:
                    curves = daten[family]
                    levels = report[family].get("computed_levels") or report[family]["levels"]
                    self.assertEqual(len(curves), len(levels), (fluid, diagram_id, family))
                    for level, curve in zip(levels, curves):
                        with self.subTest(fluid=fluid, diagram=diagram_id, family=family, level=level):
                            values, total = constant_values(diagram_id, fluid, family, curve)
                            # g-p-T: Isochorenpunkte im Nassdampfgebiet liegen auf der Saettigungslinie, dort ist (T, p)
                            # kein eindeutiger Zustand (g_l = g_g); nur die einphasigen Punkte sind rekonstruierbar.
                            fraction = 0 if (diagram_id, family) == ("gpt", "isochoren") else total // 3
                            self.assertGreaterEqual(len(values), max(2, fraction))
                            tol = 2e-3 if family == "isothermen" else 1e-5                 # K absolut bzw. relativ
                            if family == "isothermen":
                                np.testing.assert_allclose(values, level, rtol=0, atol=tol)
                            else:
                                np.testing.assert_allclose(values, level, rtol=tol, atol=0)

    def test_no_curve_bridges_over_nan(self):
        for diagram_id in DIAGRAMS.ids():
            daten = compute("Water", diagram_id, _config("Water"))
            for family in FAMILIES:
                for curve in daten[family]:
                    curve = np.asarray(curve, dtype=float)
                    self.assertEqual(curve.shape[1], 3)
                    nan_rows = ~np.isfinite(curve).all(axis=1)
                    # NaN steht immer zeilenweise (Linien werden dort getrennt), nie teilweise innerhalb einer Zeile
                    partial = np.isnan(curve).any(axis=1) & np.isfinite(curve).any(axis=1)
                    self.assertTrue(np.all(~partial | (np.isnan(curve[:, 1:]).all(axis=1))))
                    self.assertTrue(np.isfinite(curve[~nan_rows]).all())


@unittest.skipIf(SKIP, SKIP or "")
class SeparationOfLevelsAndGrids(unittest.TestCase):
    def test_surface_arrays_are_independent_of_isoline_levels(self):
        fine = compute("Water", "pvt", {"isothermen": {"factor": 1.05}, "isochoren": {"factor": 1.5},
                                         "isobaren": {"factor": 1.3}})
        coarse = compute("Water", "pvt", {"isothermen": {"factor": 3.14}})
        for key in ("surface", "binodal_liquid", "binodal_vapor", "two_phase_surface", "phase_mask"):
            np.testing.assert_array_equal(fine[key], coarse[key])                       # Oberflaechenrasterpfad unveraendert
        self.assertGreater(fine["isoline_report"]["isothermen"]["n_levels"],
                           coarse["isoline_report"]["isothermen"]["n_levels"])
        for diagram_id in ("uvt", "gpt", "fvt", "pst"):
            a = compute("Water", diagram_id, {"isothermen": {"factor": 1.05}})
            b = compute("Water", diagram_id, {"isothermen": {"factor": 3.14}})
            for key in ("surface_liquid", "surface_vapor", "surface_high_temperature", "binodal_liquid", "binodal_vapor"):
                np.testing.assert_array_equal(a[key], b[key])

    def test_include_isolines_false_keeps_empty_lists(self):
        import coolprop_heos_gpt as GPT
        daten = GPT.berechne_gpt_daten("Water", include_homogeneous=False, include_isolines=False)
        self.assertEqual(daten["isothermen"], [])
        import coolprop_heos_pht as PHT
        self.assertNotIn("isothermen", PHT.berechne_pht_daten("Water"))


@unittest.skipIf(SKIP, SKIP or "")
class CacheIdentity(unittest.TestCase):
    """C02: Aenderungen der Isolinien-Berechnungsparameter veraendern den Cache-Schluessel der betroffenen Diagramme."""

    def key(self, diagram_id, isolines=None):
        _, raster, visu = M.lade_konfiguration()
        return cache_key_for(DiagramTask(diagram_id, "Water", raster, visu, isolines=isolines))

    def test_isoline_parameters_are_part_of_every_key(self):
        for diagram_id in DIAGRAMS.ids():
            base = self.key(diagram_id)
            changed = self.key(diagram_id, {"isothermen": {"factor": 1.2}})
            self.assertNotEqual(base, changed, diagram_id)
            self.assertEqual(changed, self.key(diagram_id, {"isothermen": {"factor": 1.2}}))
            self.assertNotEqual(changed, self.key(diagram_id, {"isothermen": {"factor": 1.2, "start": 500.0}}))
            # Vorgabe == config.yaml-Vorgabe: identischer Schluessel (kein unnoetiges Neuberechnen)
            self.assertEqual(base, self.key(diagram_id, M.lade_isolinien_vorgaben()[diagram_id]))


@unittest.skipIf(SKIP, SKIP or "")
class FluidMatrix(unittest.TestCase):
    """TS-VP07-1: sieben Diagrammtypen fuer Wasser, CO2, R134a und Argon mit Standardisolinien (Faktor 3.14)."""

    def test_default_isolines_for_four_fluids_and_seven_diagrams(self):
        for fluid in FLUIDS:
            for diagram_id in DIAGRAMS.ids():
                with self.subTest(fluid=fluid, diagram=diagram_id):
                    daten = compute(fluid, diagram_id, None)
                    for family in FAMILIES:
                        self.assertIn(family, daten)
                    report = daten["isoline_report"]
                    self.assertTrue(any(report[f]["n_levels"] for f in FAMILIES))
                    for family in FAMILIES:
                        self.assertLessEqual(report[family]["n_levels"], 200)


@unittest.skipIf(SKIP, SKIP or "")
class RendererActors(unittest.TestCase):
    """L05: keine fehlenden Isolinien-Actors (insbesondere p-h-T), Layer einzeln schaltbar."""

    def test_isoline_actors_exist_for_all_seven_diagrams(self):
        import pyvista as pv
        import visualization as V
        for diagram_id in DIAGRAMS.ids():
            with self.subTest(diagram=diagram_id):
                daten = compute("Water", diagram_id, _config("Water", 1.1, 1.5, 1.6))
                plotter = pv.Plotter(off_screen=True)
                actors = getattr(V, DIAGRAMS.get(diagram_id).render_name)(plotter, daten)
                for family in FAMILIES:
                    flat = actors[family] if isinstance(actors[family], list) else [actors[family]]
                    self.assertTrue(any(a is not None for a in flat), (diagram_id, family))
                    self.assertIn(family, DIAGRAMS.get(diagram_id).layers)
                plotter.close()


if __name__ == "__main__":
    unittest.main()
