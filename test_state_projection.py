import unittest

import numpy as np

import diagram_coords
import thermo_state as ts
from state_projection import project_state


class StateProjectionTests(unittest.TestCase):
    def setUp(self):
        self.state = ts.ThermoState(
            id="state-1",
            fluid="Nitrogen",
            source_diagram="fvt",
            point_id="pick",
            T=320.0,
            p=2e5,
            v_molar=0.03,
            h_molar=10000.0,
            s_molar=200.0,
            u_molar=8000.0,
            f_molar=7000.0,
            g_molar=9000.0,
            provenance={"molar_mass_kg_per_mol": 0.028},
        )
        self.data = {
            "pvt": {"critical_point": (600.0, 10e6, 0.001)},
            "pht": {"critical_point": (600.0, 10e6, 400000.0)},
            "pst": {"critical_point": (600.0, 10e6, 3000.0)},
            "hsp": {"critical_point": (400000.0, 3000.0, 10e6)},
            "uvt": {"critical_point": (600.0, 0.001, 100000.0)},
            "gpt": {"critical_point": (600.0, 10e6, 350000.0)},
            "fvt": {"critical_point": (600.0, 0.001, 250000.0)},
        }

    def test_projects_molar_state_into_each_diagram_native_axes_and_renderer_coordinates(self):
        expected_native = {
            "pvt": (320.0, 0.03, 2e5),
            "pht": (320.0, 10000.0 / 0.028, 2e5),
            "pst": (320.0, 200.0 / 0.028, 2e5),
            "hsp": (10000.0 / 0.028, 200.0 / 0.028, 2e5),
            "uvt": (320.0, 0.03, 8000.0 / 0.028),
            "gpt": (320.0, 2e5, 9000.0 / 0.028),
            "fvt": (320.0, 0.03, 7000.0 / 0.028),
        }
        for diagram_id, native in expected_native.items():
            with self.subTest(diagram=diagram_id):
                result = project_state(self.state, diagram_id, self.data[diagram_id])
                self.assertTrue(result.valid, result.reason)
                np.testing.assert_allclose(result.native_values, native)
                adapter = __import__("state_store").adapter_for(diagram_id)
                expected_display = adapter.to_display(np.asarray(native), self.data[diagram_id])
                np.testing.assert_allclose(result.display_xyz, expected_display)

    def test_refuses_states_without_required_values_or_mass_conversion(self):
        incomplete = self.state.with_changes(v_molar=None)
        result = project_state(incomplete, "pvt", self.data["pvt"])
        self.assertFalse(result.valid)
        self.assertIn("v", result.reason)

        no_mass = self.state.with_changes(provenance={})
        result = project_state(no_mass, "pht", self.data["pht"])
        self.assertFalse(result.valid)
        self.assertIn("Molare Masse", result.reason)

    def test_refuses_synthetic_and_nonpositive_log_pressure_states(self):
        synthetic = self.state.with_changes(eos_source=ts.SOURCE_SYNTHETIC, h_molar=None, s_molar=None,
                                            u_molar=None, f_molar=None, g_molar=None)
        self.assertFalse(project_state(synthetic, "pvt", self.data["pvt"]).valid)

        nonpositive = self.state.with_changes(p=0.0)
        result = project_state(nonpositive, "pvt", self.data["pvt"])
        self.assertFalse(result.valid)
        self.assertIn("Druck", result.reason)

    def test_ambiguous_gpt_saturation_cannot_project_an_assumed_phase(self):
        state = self.state.with_changes(source_diagram="gpt", stability="saturated", Q=None)
        result = project_state(state, "pvt", self.data["pvt"])
        self.assertFalse(result.valid)
        self.assertIn("Q", result.reason)
        self.assertTrue(project_state(state, "gpt", self.data["gpt"]).valid)
        self.assertTrue(project_state(state.with_changes(Q=0.0), "pvt", self.data["pvt"]).valid)

    def test_nonpositive_pressure_is_allowed_on_nonlogarithmic_target_axes(self):
        state = self.state.with_changes(p=-1e5)
        for diagram_id in ("pvt", "gpt"):
            for axis in ("asinh", "symlog", "linear"):
                with self.subTest(diagram=diagram_id, axis=axis):
                    target = {**self.data[diagram_id], "pressure_axis": axis}
                    result = project_state(state, diagram_id, target)
                    self.assertTrue(result.valid, result.reason)
                    self.assertTrue(np.isfinite(result.display_xyz).all())


if __name__ == "__main__":
    unittest.main()
