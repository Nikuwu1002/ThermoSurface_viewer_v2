import unittest

try:
    import hover_data
except ImportError:
    hover_data = None


class HoverDataTests(unittest.TestCase):
    def test_formats_each_diagram_in_its_physical_axis_units(self):
        self.assertIsNotNone(hover_data, "hover_data module should exist")
        format_hover_values = hover_data.format_hover_values
        examples = {
            "pvt": ((300.0, 0.02, 2e6), "pvt: T=300 K, v=0.02 m³/mol, p=2 MPa"),
            "pht": ((300.0, 120000.0, 2e6), "pht: T=300 K, h=120 kJ/kg, p=2 MPa"),
            "pst": ((300.0, 1500.0, 2e6), "pst: T=300 K, s=1.5 kJ/(kg K), p=2 MPa"),
            "hsp": ((120000.0, 1500.0, 2e6), "hsp: h=120 kJ/kg, s=1.5 kJ/(kg K), p=2 MPa"),
            "uvt": ((300.0, 0.02, 90000.0), "uvt: T=300 K, v=0.02 m³/mol, u=90 kJ/kg"),
            "gpt": ((300.0, 2e6, 110000.0), "gpt: T=300 K, p=2 MPa, g=110 kJ/kg"),
            "fvt": ((300.0, 0.02, 70000.0), "fvt: T=300 K, v=0.02 m³/mol, f=70 kJ/kg"),
        }
        for diagram_id, (values, expected) in examples.items():
            with self.subTest(diagram_id=diagram_id):
                self.assertEqual(format_hover_values(diagram_id, values), expected + " (Gitterpunkt; HEOS)")

    def test_rejects_nonfinite_or_unknown_hover_values(self):
        self.assertIsNotNone(hover_data, "hover_data module should exist")
        format_hover_values = hover_data.format_hover_values
        self.assertIsNone(format_hover_values("pvt", (300.0, float("nan"), 2e6)))
        self.assertIsNone(format_hover_values("unknown", (1.0, 2.0, 3.0)))

    def test_marks_homogeneous_and_synthetic_sources(self):
        self.assertIsNotNone(hover_data, "hover_data module should exist")
        format_hover_values = hover_data.format_hover_values
        values = (300.0, 0.02, 2e6)
        self.assertIn("HEOS homogen", format_hover_values("pvt", values, "HEOS_homogeneous"))
        self.assertIn("Modell", format_hover_values("pvt", values, "synthetic_substitute"))

    def test_structured_surface_metadata_matches_vtk_point_order(self):
        import numpy as np
        import visualization

        native = np.array([[[300.0, 0.01, 1e5], [300.0, 0.02, 2e5]],
                           [[400.0, 0.01, 3e5], [400.0, 0.02, 4e5]]])
        mesh = visualization._new_surface_mesh(native, (1.0, 1.0, 1.0), lambda values, _critical: values, "fvt")
        np.testing.assert_array_equal(mesh.point_data["ts_hover_native"], mesh.points)

    def test_pvt_two_phase_hover_retains_physical_values_with_and_without_critical_cap(self):
        import numpy as np
        import visualization

        class CapturePlotter:
            def add_mesh(self, mesh, **_kwargs):
                return mesh

        surface = np.array([[[300.0, 0.001, 1e5], [300.0, 0.01, 1e5]],
                            [[400.0, 0.001, 2e5], [400.0, 0.01, 2e5]]])
        for reaches_critical in (False, True):
            with self.subTest(reaches_critical=reaches_critical):
                data = {"two_phase_surface": surface, "critical_point": (600.0, 5e6, 0.002),
                        "two_phase_reaches_critical": reaches_critical}
                mesh = visualization.add_two_phase_surface(CapturePlotter(), data)
                expected_points = list(surface.reshape(-1, 3))
                if reaches_critical:
                    expected_points.append(np.array([600.0, 0.002, 5e6]))
                for values in expected_points:
                    xyz = visualization.to_xyz(values, data["critical_point"])
                    point_id = np.argmin(np.linalg.norm(mesh.points - xyz, axis=1))
                    text = hover_data.format_mesh_hover_point(mesh, point_id)
                    self.assertEqual(text, hover_data.format_hover_values("pvt", values))

    def test_clipped_hover_retains_native_coordinates_and_marks_interpolated_points(self):
        import numpy as np
        import pyvista as pv
        import visualization
        from viewport_clipping import ClipRegistry

        native = np.array([[[1.0, 1.0, 1.0], [1.0, 2.0, 1.0]],
                           [[2.0, 1.0, 2.0], [2.0, 2.0, 2.0]]])
        mesh = visualization._new_surface_mesh(native, (1.0, 1.0, 1.0), lambda values, _critical: values, "fvt")
        actor = pv.Actor(mapper=pv.DataSetMapper(dataset=mesh))
        registry = ClipRegistry()
        registry.register({"surface": actor})
        registry.apply((1.5, 2.1, 0.0, 3.0, 0.0, 3.0))
        clipped = pv.wrap(actor.mapper.GetInput())
        for i, point in enumerate(clipped.points):
            text = hover_data.format_mesh_hover_point(clipped, i)
            self.assertIsNotNone(text)
            np.testing.assert_allclose(clipped.point_data["ts_hover_native"][i], point)
            self.assertIn("Schnittpunkt" if point[0] == 1.5 else "Gitterpunkt", text)


if __name__ == "__main__":
    unittest.main()
