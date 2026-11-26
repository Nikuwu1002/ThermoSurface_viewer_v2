"""Numeriktests des Dampfspinodal-Solvers; ohne CoolProp zunaechst Fake-EOS.

Echte HEOS-Faelle werden bei fehlendem CoolProp explizit uebersprungen.
Zum Ausfuehren: python -m unittest test_spinodal_extension -v
"""
from __future__ import annotations

import importlib.util
import unittest
from types import SimpleNamespace

import numpy as np

from spinodal_vapor_solver import (
    _brackets, _check_options, calc_spinodale_vapor_temperature,
    compare_spinodal_overlap,
)


CP = SimpleNamespace(iphase_gas=1, DmolarT_INPUTS=2, QT_INPUTS=3,
                     iP=4, iDmolar=5, iT=6)


class AnalyticState:
    """Analytische Druckisothermen mit genau zwei stationaeren Stellen.

    Dies ist bewusst keine Stoffdatenreferenz: Es testet die algorithmische
    Unterscheidung des Maximums vom Minimum, unabhaengig von CoolProp.
    """
    def __init__(self, *, nan_at=None, bad_T=None):
        self.phase = None
        self.T = 400.0
        self.rho = 1.0
        self.nan_at = nan_at
        self.bad_T = bad_T
        self.unset_calls = 0

    @staticmethod
    def roots(T):
        return 3.0 + 0.02*(T-400.0), 7.0 + 0.015*(T-400.0)

    def primitive(self, rho, T):
        a, b = self.roots(T)
        return rho**3/3 - (a+b)*rho**2/2 + a*b*rho

    def base(self, T):
        return 1e5 + (T-400)*8e3

    def p(self):
        return self.base(self.T) + 1000*self.primitive(self.rho, self.T)

    def update(self, code, x, T):
        if code == CP.QT_INPUTS:
            self.T = T
            self.rho = 1.
        else:
            self.rho = x
            self.T = T
        if self.bad_T is not None and self.bad_T[0] <= T <= self.bad_T[1]:
            raise ValueError("Simulierter Zustand nicht auswertbar")

    def rhomolar(self): return self.rho
    def first_partial_deriv(self, *_):
        a, b = self.roots(self.T)
        if self.nan_at is not None and self.rho < self.nan_at:
            return float('nan')
        return 1000*(self.rho-a)*(self.rho-b)
    def specify_phase(self, phase): self.phase = phase
    def unspecify_phase(self):
        self.phase = None
        self.unset_calls += 1
    def T_critical(self): return 500.
    def p_critical(self): return 5e6
    def rhomolar_critical(self): return 10.
    def Ttriple(self): return 200.
    def Tmin(self): return 180.
    def Tmax(self): return 600.


class VaporSolverSyntheticTests(unittest.TestCase):
    def solve(self, **kwargs):
        states = []
        def factory():
            obj = AnalyticState(**kwargs)
            states.append(obj)
            return obj
        results = calc_spinodale_vapor_temperature(
            "Synthetic", 400, 450, n_T=30, n_scan=64,
            with_diagnostics=True, state_factory=factory, cp_module=CP)
        self.assertIsNone(states[0].phase)
        return (*results, states)

    def test_chooses_maximum_not_minimum(self):
        points, report, _ = self.solve()
        self.assertEqual(report["valid_points"], 30)
        self.assertEqual(report["success_fraction"], 1.)
        self.assertEqual(len(report["missing_temperatures"]), 0)
        self.assertGreater(report["found_minima"], 0)
        self.assertLess(report["max_residual"], 1e-9)
        rho = 1/points[:, 1]
        roots = np.array([AnalyticState.roots(t)[0] for t in points[:, 0]])
        np.testing.assert_allclose(rho, roots, rtol=1e-10)
        self.assertTrue(np.all(np.diff(points[:, 2]) > 0))

    def test_exact_node_and_nan_separation(self):
        roots = _brackets(np.array([1., 2., 3., 4., 5.]),
                          np.array([4., 1., 0., -1., 2.]))
        self.assertIn((3., 3., "maximum"), roots)
        self.assertIn((4., 5., "minimum"), roots)
        self.assertEqual(_brackets(np.array([1., 2., 3.]),
                                   np.array([2., np.nan, -2.])), [])

    def test_failed_saturation_rows_stay_nan(self):
        points, report, _ = self.solve(bad_T=(418.0, 424.0))
        failed = (points[:, 0] >= 418.) & (points[:, 0] <= 424.)
        self.assertTrue(np.isnan(points[failed, 1:]).all())
        self.assertFalse(report["is_complete"])
        self.assertLess(report["success_fraction"], 1.)
        self.assertTrue(report["large_gaps"])
        # Nicht einfach den naechsten vorhandenen Punkt duplizieren.
        self.assertEqual(len(np.unique(points[:, 0])), len(points))

    def test_multiple_credible_maxima_are_rejected(self):
        class MultiPeakState(AnalyticState):
            deriv_poly = -np.poly1d(np.poly([2., 3., 4., 5., 6.]))
            primitive_poly = np.polyint(deriv_poly)
            def first_partial_deriv(self, *_):
                return 1000. * float(self.deriv_poly(self.rho))
            def p(self):
                return self.base(self.T) + 1000.*float(self.primitive_poly(self.rho))
        points, info = calc_spinodale_vapor_temperature(
            "SyntheticMulti", 400., 405., n_T=8, n_scan=200,
            state_factory=MultiPeakState, cp_module=CP,
            with_diagnostics=True)
        self.assertEqual(info["valid_points"], 0)
        self.assertEqual(len(info["ambiguous_temperatures"]), 8)
        self.assertTrue(np.isnan(points[:, 1:]).all())

    def test_invalid_settings(self):
        for settings in [
            (410, 400, 20, 40, 1e-8),
            (400, 450, 1, 40, 1e-8),
            (400, 450, 20, 10, 1e-8),
            (400, 450, 20, 40, 0.),
            (400, 450, 10000, 10000, 1e-8),
        ]:
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                _check_options(*settings)

    def test_comparison_never_adds_points(self):
        points, report, _ = self.solve()
        old = points[::4].copy()
        original = points.copy()
        comparison = compare_spinodal_overlap(points, old, critical_T=500.)
        self.assertGreater(comparison["compared"], 0)
        self.assertLess(comparison["max_rel_volume"], 1e-2)
        np.testing.assert_array_equal(points, original)

    def test_critical_temperature_cap_and_mismatch(self):
        points, report = calc_spinodale_vapor_temperature(
            "Synthetic", 450, 520, n_T=12, n_scan=40,
            with_diagnostics=True, state_factory=AnalyticState, cp_module=CP)
        self.assertLess(points[-1, 0], 500.)
        self.assertLessEqual(report["evaluated_range"][1], 500.)
        with self.assertRaises(ValueError):
            calc_spinodale_vapor_temperature(
                "Synthetic", 501, 520, n_T=12, n_scan=40,
                state_factory=AnalyticState, cp_module=CP)


@unittest.skipUnless(importlib.util.find_spec("CoolProp") is not None,
                     "HEOS-Bibliothek auf diesem Rechner nicht vorhanden")
class VaporSolverHEOSTests(unittest.TestCase):
    WATER_REF = {  # Nutzer-Messungen vom 19.10.2026; dimensionslos, gerundet.
        580.: (.274950, 3.637028, .515090, 1.202898),
        600.: (.355207, 2.815258, .628857, 1.123961),
        620.: (.464580, 2.152482, .762474, 1.058026),
        630.: (.541137, 1.847960, .839532, 1.030857),
        640.: (.658638, 1.518284, .926810, 1.009076),
    }

    def test_water_reference_points(self):
        import CoolProp.CoolProp as CP_real
        from coolprop_heos_common import _new_state
        state = _new_state("Water")
        rho_c, pc = state.rhomolar_critical(), state.p_critical()
        for T, ref in self.WATER_REF.items():
            with self.subTest(T=T):
                data, report = calc_spinodale_vapor_temperature(
                    "Water", T, T+0.0001, n_T=2, n_scan=300,
                    with_diagnostics=True)
                self.assertEqual(report["valid_points"], 2)
                state.update(CP_real.QT_INPUTS, 1., T)
                psat = state.p()
                rho = 1/data[0, 1]
                actual = (rho/rho_c, data[0, 1]*rho_c, data[0, 2]/pc,
                          data[0, 2]/psat)
                np.testing.assert_allclose(actual, ref, rtol=1e-4, atol=2e-5)
                self.assertLess(report["max_residual"], 1e-8)

    def test_water_full_and_r134a(self):
        from CoolProp import AbstractState
        for fluid in ("Water", "R134a"):
            with self.subTest(fluid=fluid):
                state = AbstractState("HEOS", fluid)
                Tc = state.T_critical()
                T_start = 580. if fluid == "Water" else 0.9*Tc
                T_end = 647. if fluid == "Water" else 0.997*Tc
                points, report = calc_spinodale_vapor_temperature(
                    fluid, T_start, T_end, n_T=250, n_scan=200,
                    with_diagnostics=True)
                self.assertGreater(report["valid_points"], 100)
                self.assertTrue(report["covers_start"])
                self.assertTrue(np.all(points[np.isfinite(points).all(axis=1), 1] > 0))
                old = AbstractState("HEOS", fluid)
                old.build_spinodal()




class VaporSolverIntegrationTests(unittest.TestCase):
    """Testet die echte main.py-Verzweigung gegen kontrollierte Moduldoubles."""

    @classmethod
    def setUpClass(cls):
        try:
            from Test import AP3SyntheticTests
        except ModuleNotFoundError as exc:
            if exc.name == "Test":
                raise unittest.SkipTest("Nicht mitgelieferte Testabhaengigkeit: Test.py") from exc
            raise
        AP3SyntheticTests.setUpClass()
        cls.pvt = AP3SyntheticTests.pvt_ap3

    def _run(self, *, fraction=1., min_points=20):
        import importlib.util
        import sys
        from types import ModuleType
        from pathlib import Path
        from unittest.mock import patch
        pvt = self.pvt
        critical = (500., 1e6, 1e-4)
        gui_stub = ModuleType("gui")
        gui_stub.ThermoSurfaceWindow = type("ThermoSurfaceWindow", (), {})
        pyside = ModuleType("PySide6")
        widgets = ModuleType("PySide6.QtWidgets")
        widgets.QApplication = type("QApplication", (), {})
        pyside.QtWidgets = widgets
        common = ModuleType("coolprop_heos_common")
        common.critical_point = lambda _: critical
        common.fluid_limits = lambda _: (200., 200., 700.)
        deps = {"gui": gui_stub, "PySide6": pyside,
                "PySide6.QtWidgets": widgets,
                "coolprop_heos_common": common,
                "coolprop_heos_pvt": pvt}
        for n in ("coolprop_heos_pht", "coolprop_heos_pst", "coolprop_heos_hsp", "coolprop_heos_uvt"):
            deps[n] = ModuleType(n)
        spec = importlib.util.spec_from_file_location(
            "_spinodal_integration_mock", Path(__file__).resolve().parents[1] / "01_Programm" / "main.py")
        app = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, deps):
            spec.loader.exec_module(app)
        app.DEBUG_MODE = False
        raster = {"isotherm_grid": {}, "spinodale_grid": {
            "T_start_spinodal": 400., "T_stop_spinodal": 490.,
            "n_T": 20, "n_scan": 40, "residual_tolerance": 1e-8}}
        visu = {"binodale_visu": {}, "two_phase_visu": {"Q_points": 3},
                "spinodale_visu": {"T_start_spinodal": 400., "T_stop_spinodal": 490.},
                "isotherm_visu": {}, "isochore_visu": {}, "isobare_visu": {}}
        liquid = np.array([[400., 5e-5, 1e5], [420., 6e-5, 2e5], [490., 8e-5, 8e5]])
        old = np.array([[400., 2e-4, 1e5], [420., 1.8e-4, 2e5], [490., 1.4e-4, 8e5]])
        t = np.linspace(400., 490., 20)
        new = np.column_stack((t, np.linspace(2e-4, 1.4e-4, 20),
                               np.linspace(1e5, 8e5, 20)))
        report = {"valid_points": min_points, "requested_points": 20,
                  "covers_start": fraction >= 0.9,
                  "success_fraction": fraction, "ambiguous_temperatures": []}
        arr = np.zeros((4, 3, 3))
        checks = {key: {"passed": True} for key in ("liquid", "vapor")}
        with patch.multiple(pvt,
                            validate_rasters=lambda *_: None,
                            validate_heos_temperature_ranges=lambda *_: None,
                            build_pvt_temperature_grid=lambda *_: (np.array([400., 420., 500., 550.]), np.array([400.,420.])),
                            multi_isotherm=lambda **_: arr,
                            multi_binodale=lambda **_: (liquid[:2], old[:2]),
                            berechne_zweiphasenmaske=lambda *_: np.zeros((4,3), dtype=bool),
                            multi_isochore=lambda **_: arr,
                            multi_isobare=lambda **_: [],
                            calc_spinodale=lambda *_: (liquid, old),
                            calc_spinodale_vapor_temperature=lambda *_a, **_k: (new, report),
                            compare_spinodal_overlap=lambda *_a, **_k: {
                                "compared": 4, "max_rel_volume": 0.005, "max_rel_pressure": 0.005},
                            validate_pvt_phase_boundaries=lambda *_: None), \
             patch.object(app, "test_spinodale_ableitung", return_value=checks):
            return app.berechne_daten("Synthetic", raster, visu)

    def test_new_solver_selected_if_qualified(self):
        data = self._run()
        self.assertEqual(data["spinodal_solver"]["selected"], "temperature_heos")
        self.assertEqual(len(data["spinodale_vapor"]), 20)
        self.assertTrue(data["spinodal_validation"]["vapor"]["derivative"]["passed"])

    def test_weak_solver_keeps_legacy_branch(self):
        data = self._run(fraction=0.10, min_points=2)
        self.assertEqual(data["spinodal_solver"]["selected"], "legacy")
        self.assertEqual(len(data["spinodale_vapor"]), 3)
        self.assertIsNotNone(data["spinodal_solver"]["fallback_reason"])


if __name__ == "__main__":
    unittest.main()
