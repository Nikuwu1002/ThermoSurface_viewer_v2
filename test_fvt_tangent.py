import unittest
from unittest.mock import patch

import numpy as np

import CoolProp.CoolProp as CP
from CoolProp import AbstractState

from fvt_tangent import calculate_fvt_tangent, homogeneous_fvt_isotherm
from fvt_tangent import tangent_display_points, tangent_plot_data
import diagram_coords


class FvtTangentTests(unittest.TestCase):
    def test_tangent_segments_follow_the_physical_curve_on_logarithmic_volume_axis(self):
        tangent = calculate_fvt_tangent("Water", 300.0)
        volumes = tangent.line_v
        physical_midpoints = np.sqrt(volumes[:-1] * volumes[1:])
        chord_midpoints = (volumes[:-1] + volumes[1:]) / 2.0
        max_error = np.max(np.abs(chord_midpoints - physical_midpoints) / np.diff(volumes))
        self.assertLess(max_error, 0.05)

    def test_tangent_display_points_use_mass_specific_f_and_logarithmic_volume_transform(self):
        tangent = calculate_fvt_tangent("Water", 300.0, n_points=7)
        data = {"meta": {"molar_mass": 0.01801528},
                "critical_point": (647.096, 0.0031, 200000.0)}
        actual = tangent_display_points(tangent, data)
        native = np.column_stack((np.full(7, 300.0), tangent.line_v,
                                  tangent.line_f_molar / data["meta"]["molar_mass"]))
        np.testing.assert_allclose(actual, diagram_coords.to_xyz_fvt(native, data["critical_point"]))
        self.assertNotAlmostEqual(actual[3, 2], (actual[0, 2] + actual[-1, 2]) / 2.0)

    def test_two_dimensional_plot_payload_exposes_the_zero_volume_intercept(self):
        tangent = calculate_fvt_tangent("Water", 300.0, n_points=7)
        payload = tangent_plot_data(tangent)
        np.testing.assert_array_equal(payload["volume"], tangent.line_v)
        np.testing.assert_array_equal(payload["helmholtz"], tangent.line_f_molar)
        self.assertAlmostEqual(payload["intercept"], tangent.mu)

    def test_homogeneous_isotherm_samples_only_outside_the_two_phase_interval(self):
        tangent = calculate_fvt_tangent("Water", 300.0, n_points=7)
        data = {"meta": {"molar_mass": 0.01801528, "v_min": tangent.v_l * 0.8,
                         "v_max": tangent.v_v * 1.2}}
        branches = homogeneous_fvt_isotherm(tangent, data, n_points=24)
        self.assertEqual(len(branches), 2)
        for volumes, energies in branches:
            self.assertTrue(np.isfinite(volumes).all())
            self.assertTrue(np.isfinite(energies).all())
            self.assertTrue(np.all((volumes < tangent.v_l) | (volumes > tangent.v_v)))

    def test_homogeneous_isotherm_does_not_connect_across_invalid_samples(self):
        tangent = calculate_fvt_tangent("Water", 300.0)
        data = {"meta": {"molar_mass": 0.01801528, "v_min": tangent.v_l * 0.8,
                         "v_max": tangent.v_v * 1.2}}
        liquid_samples = np.geomspace(data["meta"]["v_min"], tangent.v_l * (1.0 - 1e-8), 6)
        vapor_samples = np.geomspace(tangent.v_v * (1.0 + 1e-8), data["meta"]["v_max"], 6)

        def evaluate(_state, _temperature, volume):
            if np.isclose(volume, liquid_samples[2], rtol=1e-10, atol=0) or np.isclose(
                    volume, vapor_samples[2], rtol=1e-10, atol=0):
                raise ValueError("invalid EOS sample")
            return 0.0, 0.0, volume

        with patch("coolprop_heos_fvt._f_at_tv", side_effect=evaluate):
            segments = homogeneous_fvt_isotherm(tangent, data, n_points=12)
        self.assertEqual(len(segments), 4)
        for (volumes, _energies), expected in zip(
                segments, (liquid_samples[:2], liquid_samples[3:], vapor_samples[:2], vapor_samples[3:])):
            np.testing.assert_allclose(volumes, expected)

    def test_homogeneous_isotherm_reports_when_no_segment_can_be_calculated(self):
        tangent = calculate_fvt_tangent("Water", 300.0)
        data = {"meta": {"molar_mass": 0.01801528, "v_min": tangent.v_l * 0.8,
                         "v_max": tangent.v_v * 1.2}}
        with patch("coolprop_heos_fvt._f_at_tv", side_effect=ValueError("invalid EOS sample")):
            with self.assertRaisesRegex(ValueError, "Keine zusammenhängenden"):
                homogeneous_fvt_isotherm(tangent, data, n_points=12)

    def test_tangent_matches_heos_saturation_for_multiple_pure_fluids(self):
        for fluid, temperature in (("Water", 300.0), ("CO2", 250.0), ("R134a", 300.0)):
            with self.subTest(fluid=fluid):
                result = calculate_fvt_tangent(fluid, temperature, n_points=21)
                state = AbstractState("HEOS", fluid)
                state.update(CP.QT_INPUTS, 0.0, temperature)
                expected_pressure = state.p()
                expected_v_l = 1.0 / state.rhomolar()
                expected_g_l = state.gibbsmolar()
                state.update(CP.QT_INPUTS, 1.0, temperature)
                expected_v_v = 1.0 / state.rhomolar()
                expected_g_v = state.gibbsmolar()

                self.assertAlmostEqual(result.p_sat, expected_pressure, delta=expected_pressure * 1e-9)
                self.assertAlmostEqual(result.v_l, expected_v_l, delta=expected_v_l * 1e-9)
                self.assertAlmostEqual(result.v_v, expected_v_v, delta=expected_v_v * 1e-9)
                self.assertAlmostEqual(result.g_l_molar, expected_g_l, delta=max(1e-5, abs(expected_g_l) * 1e-9))
                self.assertAlmostEqual(result.g_v_molar, expected_g_v, delta=max(1e-5, abs(expected_g_v) * 1e-9))
                self.assertLess(result.pressure_relative_error, 1e-8)
                self.assertLess(abs(result.chemical_potential_delta), max(1e-4, abs(result.mu) * 1e-8))
                np.testing.assert_allclose(result.line_v[[0, -1]], [result.v_l, result.v_v])
                np.testing.assert_allclose(result.line_f_molar[[0, -1]], [result.f_l_molar, result.f_v_molar])
                self.assertEqual(len(result.line_v), 21)

    def test_rejects_non_two_phase_temperatures_and_invalid_point_counts(self):
        state = AbstractState("HEOS", "Water")
        for temperature in (state.Tmin() - 1.0, state.T_critical(), state.T_critical() + 1.0):
            with self.subTest(temperature=temperature):
                with self.assertRaises(ValueError):
                    calculate_fvt_tangent("Water", temperature)
        with self.assertRaises(ValueError):
            calculate_fvt_tangent("Water", 300.0, n_points=1)


if __name__ == "__main__":
    unittest.main()
