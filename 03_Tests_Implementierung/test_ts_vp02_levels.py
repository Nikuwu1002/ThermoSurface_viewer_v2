"""TS-VP02: zentraler geometrischer Isolinien-Levelgenerator (Testmatrix L01-L03, Konfiguration, Cache-Signatur)."""
import math
import sys
import unittest

import numpy as np

from isoline_levels import (FACTOR_DEFAULT, IsolineLevelError, config_signature, generate_levels,
                            normalise_isoline_config, resolve_levels)


class GeneratorFormula(unittest.TestCase):
    def test_L01_example_10_1000_314(self):
        levels = generate_levels(10, 1000, 3.14)
        np.testing.assert_allclose(levels, [10, 31.4, 98.596, 309.59144, 972.1171216], rtol=1e-14)
        self.assertLess(levels[-1], 1000.0)               # stop wird nicht zusaetzlich erzwungen
        self.assertEqual(len(levels), 5)

    def test_exact_formula_and_constant_ratio(self):
        for start, stop, factor in ((1e-5, 1e-2, 3.14), (1e3, 1e8, 3.14), (273.15, 1300.0, 1.1), (1.0, 1e12, 7.5)):
            levels = generate_levels(start, stop, factor)
            expected = start * factor ** np.arange(len(levels))
            np.testing.assert_array_equal(levels, expected)            # x_n = x_0 * f**n, keine kumulative Multiplikation
            np.testing.assert_allclose(levels[1:] / levels[:-1], factor, rtol=1e-13)
            self.assertTrue(np.all(np.diff(levels) > 0))
            self.assertLessEqual(levels[-1], stop * (1 + 1e-12))
            self.assertGreater(levels[-1] * factor, stop)            # der naechste Wert liegt ausserhalb

    def test_deterministic_for_identical_si_input(self):
        np.testing.assert_array_equal(generate_levels(8e6, 2.8e7, 3.14), generate_levels(8e6, 2.8e7, 3.14))

    def test_boundary_equality_included_numerically(self):
        factor = 1000.0 ** (1.0 / 3.0)
        levels = generate_levels(1.0, 1000.0, factor)
        self.assertEqual(len(levels), 4)                   # 1000 ist numerisch erreichbar und wird nicht verworfen
        self.assertAlmostEqual(levels[-1], 1000.0, places=9)

    def test_start_equals_stop_gives_single_level(self):
        np.testing.assert_array_equal(generate_levels(5.0, 5.0, 3.14), [5.0])

    def test_no_linear_interpolation_or_forced_endpoint(self):
        levels = generate_levels(10, 1000, 3.14)
        self.assertNotIn(1000.0, levels)
        self.assertFalse(np.allclose(np.diff(levels), np.diff(levels)[0]))       # nicht aequidistant

    def test_fluid_orders_of_magnitude(self):
        # T [K], v [m3/mol], p [Pa] verschiedener Stoffe: Anzahl = floor(log_f(stop/start)) + 1
        for start, stop in ((100.0, 1500.0), (3e-5, 2e-2), (1e3, 1e9), (1e-9, 1e9)):
            count = len(generate_levels(start, stop, 3.14, max_levels=200))
            self.assertEqual(count, int(math.floor(math.log(stop / start) / math.log(3.14) + 1e-9)) + 1)


class GeneratorValidation(unittest.TestCase):
    def test_L02_invalid_inputs_raise_controlled_errors(self):
        bad = ((0.0, 10, 3), (-1, 10, 3), (1, 10, 1.0), (1, 10, 0.5), (1, 10, 1 + 1e-12), (10, 1, 3),
               (float("nan"), 10, 3), (1, float("inf"), 3), (1, 10, float("nan")), (1, 10, float("-inf")),
               ("abc", 10, 3), (None, 10, 3))
        for args in bad:
            with self.subTest(args=args):
                with self.assertRaises(IsolineLevelError):
                    generate_levels(*args)

    def test_error_is_a_valueerror(self):
        self.assertTrue(issubclass(IsolineLevelError, ValueError))

    def test_L03_factor_near_one_is_guarded_without_hanging(self):
        with self.assertRaises(IsolineLevelError):
            generate_levels(1.0, 1e6, 1.0000001)
        with self.assertRaises(IsolineLevelError):
            generate_levels(1.0, 1e300, 1.001, max_levels=200)
        # knapp unter der Grenze bleibt gueltig: floor(ln 2 / ln 1.0035) + 1 = 199 Werte
        self.assertEqual(len(generate_levels(1.0, 2.0, 1.0035, max_levels=200)), 199)

    def test_max_levels_is_respected_exactly(self):
        self.assertEqual(len(generate_levels(1.0, 2.0 ** 199, 2.0, max_levels=200)), 200)
        with self.assertRaises(IsolineLevelError):
            generate_levels(1.0, 2.0 ** 200, 2.0, max_levels=200)
        for bad in (0, -1, 1.5, True, 10 ** 7):
            with self.assertRaises(IsolineLevelError):
                generate_levels(1.0, 10.0, 2.0, max_levels=bad)

    def test_overflow_is_reported(self):
        with self.assertRaises(IsolineLevelError):
            generate_levels(1.0, sys.float_info.max, 1e100, max_levels=200)   # Ueberlauf an der Gleitkommagrenze


class ConfigAndResolution(unittest.TestCase):
    def test_defaults_and_auto_markers(self):
        cfg = normalise_isoline_config(None)
        for family in ("isothermen", "isochoren", "isobaren"):
            self.assertEqual(cfg[family]["factor"], FACTOR_DEFAULT)
            self.assertIsNone(cfg[family]["start"])
            self.assertIsNone(cfg[family]["stop"])
        cfg = normalise_isoline_config({"isobaren": {"start": "auto", "stop": 0, "factor": 2}})
        self.assertIsNone(cfg["isobaren"]["start"])
        self.assertIsNone(cfg["isobaren"]["stop"])

    def test_unknown_keys_and_families_rejected(self):
        with self.assertRaises(IsolineLevelError):
            normalise_isoline_config({"isohypsen": {}})
        with self.assertRaises(IsolineLevelError):
            normalise_isoline_config({"isobaren": {"step": 2}})
        with self.assertRaises(IsolineLevelError):
            normalise_isoline_config({"isobaren": {"start": -1}})

    def test_resolve_uses_domain_for_auto_and_drops_outside_levels(self):
        domains = {"isothermen": (400.0, 1200.0), "isobaren": (1e6, 1e8), "isochoren": (1e-5, 1e-3)}
        levels, report = resolve_levels({"isobaren": {"start": 1e5, "factor": 10.0}}, domains)
        np.testing.assert_allclose(levels["isobaren"], [1e6, 1e7, 1e8], rtol=1e-12)      # 1e5 liegt ausserhalb
        self.assertEqual(report["isobaren"]["dropped_outside_domain"], [1e5])
        self.assertEqual(report["isobaren"]["n_levels"], 3)
        np.testing.assert_allclose(levels["isothermen"], [400.0])                       # 1256 K > 1200 K
        self.assertEqual(report["isochoren"]["start"], 1e-5)

    def test_resolution_error_names_the_family(self):
        with self.assertRaisesRegex(IsolineLevelError, "isochoren"):
            resolve_levels({"isochoren": {"start": 5e-4, "stop": 1e-4}}, {"isothermen": (1, 2), "isochoren": (1e-5, 1e-3),
                                                                           "isobaren": (1, 2)})

    def test_signature_changes_with_levels_only(self):
        base = config_signature({"isobaren": {"factor": 3.14}})
        self.assertEqual(base, config_signature({"isobaren": {"factor": 3.14, "start": None}}))
        self.assertNotEqual(base, config_signature({"isobaren": {"factor": 3.0}}))
        self.assertNotEqual(base, config_signature({"isobaren": {"factor": 3.14, "start": 1e6}}))


if __name__ == "__main__":
    unittest.main()
