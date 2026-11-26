"""TS-VP04/05/06: Ansichtsbereich ohne Qt/HEOS - Achsenmatrix, Grenzenpruefung, geometrisches Clipping, Actor-Registry.

Testmatrix V02-V04, V08-V10 sowie die Renderer-Seite von C01/V03/V11. Die GUI-gebundenen Faelle stehen in
test_ts_vp_gui.py (benoetigt ein OpenGL-Fenster).
"""
import unittest

import numpy as np
import pyvista as pv

import axes_view
import diagram_axes as A
import geometry_clip as G
from diagram_registry import DIAGRAMS
from viewport_bounds import ViewportState, parse_field
from viewport_clipping import ORIG_CELL, ClipRegistry
from zustandsflaeche import SRC_CLIPPED

# Synthetische kritische Zustaende (Struktur wie daten["critical_point"] je Diagramm).
CRITICAL = {
    "pvt": (647.1, 22.06e6, 5.6e-5), "pht": (647.1, 22.06e6, 2.1e6), "pst": (647.1, 22.06e6, 4.4e3),
    "hsp": (2.1e6, 4.4e3, 22.06e6), "uvt": (647.1, 5.6e-5, 2.0e6), "gpt": (647.1, 22.06e6, -1.0e6),
    "fvt": (647.1, 5.6e-5, -2.0e5),
}
# Gueltige Grenzen in Anzeigeeinheiten (x, y, z) je Diagramm.
BOUNDS = {
    "pvt": [(2e-5, 2e-3), (300.0, 800.0), (0.1, 60.0)], "pht": [(-1500.0, 2500.0), (300.0, 800.0), (0.1, 60.0)],
    "pst": [(-3.0, 4.0), (300.0, 800.0), (0.1, 60.0)], "hsp": [(-1500.0, 2500.0), (-3.0, 4.0), (0.1, 60.0)],
    "uvt": [(2e-5, 2e-3), (300.0, 800.0), (-1500.0, 2500.0)], "gpt": [(0.1, 60.0), (300.0, 800.0), (-500.0, 500.0)],
    "fvt": [(2e-5, 2e-3), (300.0, 800.0), (-500.0, 500.0)],
}


def daten_for(diagram_id, pressure_axis="log"):
    return {"critical_point": CRITICAL[diagram_id], "pressure_axis": pressure_axis}


def flat_surface(n=21, z=0.3):
    x, y = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n), indexing="ij")
    pts = np.column_stack((x.ravel(), y.ravel(), np.full(n * n, z)))
    idx = np.arange(n * n).reshape(n, n)
    t1 = np.column_stack((idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[:-1, 1:].ravel()))
    t2 = np.column_stack((idx[1:, :-1].ravel(), idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()))
    return pts, np.vstack((t1, t2))


def area_xy(points, triangles):
    a, b, c = (points[triangles[:, k]] for k in range(3))
    return 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])).sum()


class AxisMatrix(unittest.TestCase):
    def test_coordinate_contract_of_all_seven_diagrams(self):
        expected = {"pvt": "vTp", "pht": "hTp", "pst": "sTp", "hsp": "hsp", "uvt": "vTu", "gpt": "pTg", "fvt": "vTf"}
        self.assertEqual(set(A.DIAGRAM_AXES), set(DIAGRAMS.ids()))
        for diagram_id, keys in expected.items():
            self.assertEqual("".join(spec.key for spec in A.axes_for(diagram_id)), keys)

    def test_forward_transform_is_the_one_of_the_meshes_and_invertible(self):
        for diagram_id in DIAGRAMS.ids():
            d = daten_for(diagram_id)
            for index, (lo, hi) in enumerate(BOUNDS[diagram_id]):
                spec = A.axes_for(diagram_id)[index]
                si = np.linspace(lo, hi, 7) / spec.scale
                coords = A.to_display(diagram_id, d, index, si)
                self.assertTrue(np.all(np.diff(coords) > 0), (diagram_id, index))          # monoton steigend
                np.testing.assert_allclose(A.to_physical(diagram_id, d, index, coords), si, rtol=1e-9, atol=1e-12)

    def test_axis_values_match_mesh_transform(self):
        import diagram_coords as dc
        d = daten_for("pvt", "asinh")
        T, v, p = 500.0, 1e-4, 5e6
        xyz = dc.to_xyz(np.array([[T, v, p]]), d["critical_point"], "asinh")[0]
        self.assertAlmostEqual(A.to_display("pvt", d, 0, [v])[0], xyz[0])
        self.assertAlmostEqual(A.to_display("pvt", d, 1, [T])[0], xyz[1])
        self.assertAlmostEqual(A.to_display("pvt", d, 2, [p])[0], xyz[2])

    def test_units_and_titles(self):
        d = daten_for("pvt", "log")
        self.assertEqual([A.axis_title("pvt", i, d) for i in range(3)], ["v [m³/mol] (log)", "T [K]", "p [MPa] (log)"])
        self.assertEqual([A.axis_title("hsp", i, daten_for("hsp")) for i in range(3)],
                         ["h [kJ/kg]", "s [kJ/(kg K)]", "p [MPa] (log)"])
        self.assertEqual(A.axis_title("pvt", 2, daten_for("pvt", "symlog")), "p [MPa] (symlog)")


class BoundsValidation(unittest.TestCase):
    def test_valid_bounds_round_trip_for_every_diagram(self):
        for diagram_id in DIAGRAMS.ids():
            d = daten_for(diagram_id)
            bounds_si, box = A.validate_bounds(diagram_id, d, BOUNDS[diagram_id])
            self.assertEqual(len(box), 6)
            self.assertTrue(np.all(np.array(box[0::2]) < np.array(box[1::2])))
            back = A.display_bounds_to_physical(diagram_id, d, box)
            np.testing.assert_allclose(back, BOUNDS[diagram_id], rtol=1e-9)
            self.assertEqual(box, A.bounds_to_box(diagram_id, d, bounds_si))

    def test_V08_negative_energy_axes_are_accepted(self):
        for diagram_id, index in (("pht", 0), ("pst", 0), ("hsp", 0), ("hsp", 1), ("uvt", 2), ("gpt", 2), ("fvt", 2)):
            bounds = [list(b) for b in BOUNDS[diagram_id]]
            bounds[index] = (-900.0, -100.0)
            _si, box = A.validate_bounds(diagram_id, daten_for(diagram_id), [tuple(b) for b in bounds])
            phys = A.display_bounds_to_physical(diagram_id, daten_for(diagram_id), box)
            np.testing.assert_allclose(phys[index], (-900.0, -100.0), rtol=1e-9)
            ticks = A.physical_ticks("linear", -900.0, -100.0)
            self.assertTrue(all(t < 0 for t in ticks) and len(ticks) >= 3)

    def test_V09_nonpositive_pressure_only_with_valid_transformations(self):
        bounds = [(2e-5, 2e-3), (300.0, 800.0), (-5.0, 60.0)]
        with self.assertRaisesRegex(A.ViewportBoundsError, "Logarithmische Achse"):
            A.validate_bounds("pvt", daten_for("pvt", "log"), bounds)
        for axis in ("asinh", "symlog", "linear"):
            _si, box = A.validate_bounds("pvt", daten_for("pvt", axis), bounds)
            self.assertLess(box[4], 0.0 if axis in ("asinh", "symlog", "linear") else 1.0)
        # Diagramme mit fester Log-Druckachse lehnen p <= 0 immer ab (g-p-T: Druck auf der x-Achse)
        for diagram_id, index in (("pht", 2), ("pst", 2), ("hsp", 2)):
            b = [tuple(x) for x in BOUNDS[diagram_id]]
            b[index] = (-1.0, 10.0)
            with self.assertRaises(A.ViewportBoundsError):
                A.validate_bounds(diagram_id, daten_for(diagram_id), b)
        b = [tuple(x) for x in BOUNDS["gpt"]]
        b[0] = (-1.0, 10.0)
        with self.assertRaises(A.ViewportBoundsError):
            A.validate_bounds("gpt", daten_for("gpt", "log"), b)
        A.validate_bounds("gpt", daten_for("gpt", "symlog"), b)

    def test_temperature_and_volume_must_be_positive_and_finite(self):
        for bad in ([(0.0, 800.0)], [(-5.0, 800.0)], [(float("nan"), 800.0)], [(300.0, float("inf"))],
                    [(800.0, 300.0)], [(300.0, 300.0)]):
            b = [tuple(x) for x in BOUNDS["pvt"]]
            b[1] = bad[0]
            with self.assertRaises(A.ViewportBoundsError):
                A.validate_bounds("pvt", daten_for("pvt"), b)
        b = [tuple(x) for x in BOUNDS["pvt"]]
        b[0] = (-1e-5, 1e-3)
        with self.assertRaises(A.ViewportBoundsError):
            A.validate_bounds("pvt", daten_for("pvt"), b)

    def test_state_applies_atomically_and_keeps_last_valid_bounds(self):
        d = daten_for("pvt")
        state = ViewportState("pvt")
        good = ["2e-5", "2e-3", "300", "800", "0,1", "60"]            # Dezimalkomma wird akzeptiert
        box = state.apply_fields(good, d)
        self.assertTrue(state.enabled and state.user_defined and len(box) == 6)
        before = (dict(state.bounds), list(state.fields))
        for bad in (["2e-5", "2e-3", "300", "800", "0.1", "abc"], ["2e-5", "2e-3", "800", "300", "0.1", "60"],
                    ["2e-5", "2e-3", "300", "800", "0.1"], ["", "2e-3", "300", "800", "0.1", "60"]):
            with self.assertRaises(A.ViewportBoundsError):
                state.apply_fields(bad, d)
            self.assertEqual((dict(state.bounds), list(state.fields)), before)       # unveraendert
        state.reset()
        self.assertFalse(state.enabled)
        self.assertIsNone(state.box(d))
        self.assertEqual(dict(state.bounds), before[0])                               # Eingaben bleiben als Vorschlag
        self.assertEqual(parse_field(" 1,5e-3 "), 1.5e-3)

    def test_copies_are_independent_per_diagram(self):
        a, b = ViewportState("pvt"), ViewportState("hsp")
        a.apply_fields(["2e-5", "2e-3", "300", "800", "0.1", "60"], daten_for("pvt"))
        self.assertFalse(b.enabled)
        c = a.copy()
        c.reset()
        self.assertTrue(a.enabled and not c.enabled)


class Ticks(unittest.TestCase):
    def test_hover_label_marks_clipped_points_as_render_interpolation(self):
        from zustandsflaeche import hover_label
        text = hover_label(SRC_CLIPPED, 500.0, 1e-4, 2e6)
        self.assertIn("Schnittpunkt", text)
        self.assertIn("keine HEOS-Auswertung", text)
        self.assertTrue(hover_label(0, 500.0, 1e-4, 2e6).startswith("HEOS"))
        self.assertIn("Ersatz", hover_label(1, 500.0, 1e-4, 2e6))
    def test_linear_log_symlog_ticks_are_physical_numbers(self):
        self.assertEqual(A.physical_ticks("linear", 4.0, 26.0, target=5), [5.0, 10.0, 15.0, 20.0, 25.0])
        self.assertEqual(A.physical_ticks("log", 1e-5, 1e-2), [1e-5, 1e-4, 1e-3, 1e-2])
        narrow = A.physical_ticks("log", 1.5e-4, 7e-4)
        self.assertEqual(narrow, [2e-4, 5e-4])
        sym = A.physical_ticks("symlog", -5.0, 30.0, min_abs=0.2)
        self.assertIn(0.0, sym)
        self.assertTrue(min(sym) < 0 < max(sym) and all(abs(t) >= 0.2 or t == 0 for t in sym))
        self.assertEqual(A.physical_ticks("log", -1.0, 10.0), [])
        self.assertEqual(A.format_tick(1e-5), "1e-5")
        self.assertEqual(A.format_tick(0.0), "0")
        self.assertEqual(A.format_tick(972.1171216), "972.117")

    def test_tick_positions_use_the_mesh_transformation(self):
        d = daten_for("pvt", "log")
        box = A.validate_bounds("pvt", d, BOUNDS["pvt"])[1]
        edges = axes_view.choose_edges(box)
        ticks = axes_view._tick_points("pvt", d, box, edges)
        points, labels = ticks[2]
        self.assertTrue(len(points) >= 3)
        expected = A.to_display("pvt", d, 2, np.array([float(t) for t in labels]) * 1e6)
        np.testing.assert_allclose(points[:, 2], expected, rtol=1e-12)

    def test_xy_axis_titles_stay_below_the_box_even_when_camera_is_underneath(self):
        box = A.validate_bounds("pvt", daten_for("pvt"), BOUNDS["pvt"])[1]
        lo, hi = np.array(box[0::2]), np.array(box[1::2])
        center = 0.5 * (lo + hi)
        camera = (center + np.array([0.0, -3.0, -3.0]), center, np.array([0.0, 0.0, 1.0]))

        edges = axes_view.choose_edges(box, camera)

        for index in (0, 1):
            fixed, outward = edges[index]
            title_point = fixed.copy()
            title_point[index] = center[index]
            title_point += outward * 0.11 * (hi - lo)
            self.assertLess(title_point[2], box[4])

    def test_default_view_looks_from_negative_y(self):
        plotter = pv.Plotter(off_screen=True)
        self.addCleanup(plotter.close)
        box = (-2.0, 4.0, -3.0, 5.0, 10.0, 30.0)

        axes_view.set_default_view(plotter, box)

        position, focal, up = (np.asarray(point) for point in plotter.camera_position)
        np.testing.assert_allclose(focal, (1.0, 1.0, 20.0))
        self.assertLess(position[1], focal[1])
        self.assertGreater(position[2], focal[2])
        np.testing.assert_allclose(up, axes_view.DEFAULT_VIEW_UP)

    def test_default_xy_axis_titles_project_below_the_plot(self):
        box = (-1.0, 1.0, -2.0, 2.0, 0.0, 1.0)
        plotter = pv.Plotter(off_screen=True, window_size=(600, 400))
        self.addCleanup(plotter.close)
        plotter.add_mesh(pv.Cube(bounds=box), opacity=0.1)
        axes_view.set_default_view(plotter, box)
        edges = axes_view.choose_edges(box, axes_view._camera_vectors(plotter))
        size = np.maximum(np.array(box[1::2]) - np.array(box[0::2]), 1e-12)
        title_points = []
        for index in (0, 1):
            fixed, outward = edges[index]
            point = fixed.copy()
            point[index] = 0.5 * (box[2 * index] + box[2 * index + 1])
            title_points.append(point + outward * 0.11 * size)
        plotter.render()

        def display_y(point):
            plotter.renderer.SetWorldPoint(*point, 1.0)
            plotter.renderer.WorldToDisplay()
            return plotter.renderer.GetDisplayPoint()[1]

        frame_y = [display_y((x, y, z)) for x in box[0::2] for y in box[1::2] for z in box[4::2]]
        self.assertTrue(all(display_y(point) < min(frame_y) for point in title_points))


class TriangleClipping(unittest.TestCase):
    def setUp(self):
        self.points, self.triangles = flat_surface()

    def test_area_is_exact_and_all_points_lie_inside_the_box(self):
        box = (0.2, 0.7, 0.1, 0.55, 0.0, 1.0)
        result = G.clip_triangles(self.points, self.triangles, box)
        self.assertAlmostEqual(area_xy(result["points"], result["triangles"]), 0.5 * 0.45, places=12)
        lo, hi = np.array(box[0::2]), np.array(box[1::2])
        self.assertTrue(np.all(result["points"] >= lo - 1e-12) and np.all(result["points"] <= hi + 1e-12))
        self.assertTrue(np.any(result["new"]))                                # Schnittpunkte sind markiert
        # Neue Dreieckskanten liegen auf den Boxflaechen: alle neuen Punkte haben mindestens eine Boxkoordinate
        new = result["points"][result["new"]]
        on_face = np.isclose(new, lo, atol=1e-12) | np.isclose(new, hi, atol=1e-12)
        self.assertTrue(np.all(on_face.any(axis=1)))

    def test_attributes_are_interpolated_and_original_cells_tracked(self):
        box = (0.25, 0.75, 0.0, 1.0, 0.0, 1.0)
        data = {"T": self.points[:, 0] * 100.0}
        result = G.clip_triangles(self.points, self.triangles, box, data)
        np.testing.assert_allclose(result["point_data"]["T"], result["points"][:, 0] * 100.0, atol=1e-10)
        self.assertEqual(len(result["orig_cell"]), len(result["triangles"]))
        for tri, orig in zip(result["triangles"][:50], result["orig_cell"][:50]):          # Ursprungsdreieck enthaelt Punkt
            orig_xy = self.points[self.triangles[orig]][:, :2]
            centroid = result["points"][tri].mean(axis=0)[:2]
            self.assertTrue(orig_xy[:, 0].min() - 1e-12 <= centroid[0] <= orig_xy[:, 0].max() + 1e-12)

    def test_box_containing_everything_keeps_every_triangle_bitwise(self):
        result = G.clip_triangles(self.points, self.triangles, (-1, 2, -1, 2, -1, 2))
        self.assertEqual(len(result["triangles"]), len(self.triangles))
        self.assertFalse(result["new"].any())
        self.assertAlmostEqual(area_xy(result["points"], result["triangles"]), 1.0, places=12)

    def test_V04_empty_window_gives_empty_result_without_errors(self):
        result = G.clip_triangles(self.points, self.triangles, (5, 6, 5, 6, 5, 6), {"T": self.points[:, 0]})
        self.assertEqual(len(result["triangles"]), 0)
        self.assertEqual(result["points"].shape, (0, 3))
        self.assertEqual(len(result["point_data"]["T"]), 0)

    def test_V10_nan_gaps_are_never_bridged(self):
        points = self.points.copy()
        points[10] = np.nan
        result = G.clip_triangles(points, self.triangles, (-1, 2, -1, 2, -1, 2))
        used = set(self.triangles[np.isfinite(points[self.triangles]).all(axis=(1, 2))].ravel())
        self.assertLess(len(result["triangles"]), len(self.triangles))
        self.assertTrue(np.isfinite(result["points"]).all())
        self.assertEqual(len(result["triangles"]), int(np.isfinite(points[self.triangles]).all(axis=(1, 2)).sum()))
        self.assertNotIn(10, used)

    def test_z_clipping_of_an_inclined_plane(self):
        pts, tris = flat_surface()
        pts[:, 2] = pts[:, 0] + pts[:, 1]                                     # z = x + y
        result = G.clip_triangles(pts, tris, (0, 1, 0, 1, 0.5, 1.5))
        self.assertTrue(np.all(result["points"][:, 2] >= 0.5 - 1e-12) and np.all(result["points"][:, 2] <= 1.5 + 1e-12))
        self.assertAlmostEqual(area_xy(result["points"], result["triangles"]), 1.0 - 0.125 - 0.125, places=12)

    def test_repeated_clipping_always_starts_from_the_original(self):
        original = (self.points.copy(), self.triangles.copy())
        narrow = G.clip_triangles(*original, (0.4, 0.6, 0.4, 0.6, 0, 1))
        wide = G.clip_triangles(*original, (0.0, 1.0, 0.0, 1.0, 0, 1))
        self.assertLess(area_xy(narrow["points"], narrow["triangles"]), 0.05)
        self.assertAlmostEqual(area_xy(wide["points"], wide["triangles"]), 1.0, places=12)       # V03: wieder voll
        np.testing.assert_array_equal(original[0], self.points)
        np.testing.assert_array_equal(original[1], self.triangles)


class LineClipping(unittest.TestCase):
    def test_polyline_is_cut_at_the_box_and_connected_inside(self):
        t = np.linspace(0.0, 2.0 * np.pi, 400)
        line = np.column_stack((np.cos(t), np.sin(t), t / 6.0))
        result = G.clip_polylines(line, [np.arange(400)], (-0.5, 0.5, -2, 2, 0, 1))
        self.assertEqual(len(result["lines"]), 2)
        pts = result["points"]
        self.assertTrue(np.all(np.abs(pts[:, 0]) <= 0.5 + 1e-12))
        for idx in result["lines"]:
            ends = pts[[idx[0], idx[-1]]]
            self.assertTrue(np.all(np.isclose(np.abs(ends[:, 0]), 0.5, atol=1e-12)))    # Enden liegen auf der Boxflaeche

    def test_nan_split_curves_stay_separate_and_outside_curves_vanish(self):
        pts = np.array([[0.1, 0.1, 0.1], [0.2, 0.2, 0.2], [0.9, 0.9, 0.9], [0.95, 0.95, 0.95]])
        result = G.clip_polylines(pts, [np.array([0, 1]), np.array([2, 3])], (0, 0.5, 0, 0.5, 0, 0.5))
        self.assertEqual(len(result["lines"]), 1)
        self.assertEqual(G.clip_polylines(pts, [np.arange(4)], (5, 6, 5, 6, 5, 6))["lines"], [])

    def test_points_inside_mask(self):
        pts = np.array([[0.5, 0.5, 0.5], [2.0, 0.5, 0.5], [1.0, 1.0, 1.0]])
        np.testing.assert_array_equal(G.points_in_box(pts, (0, 1, 0, 1, 0, 1)), [True, False, True])


def make_actor_scene():
    plotter = pv.Plotter(off_screen=True)
    pts, tris = flat_surface(15, 0.3)
    faces = np.column_stack((np.full(len(tris), 3), tris)).ravel()
    mesh = pv.PolyData(pts, faces)
    mesh.point_data["T"] = pts[:, 0] * 100.0
    mesh.point_data["v"] = pts[:, 1]
    mesh.point_data["p_raw"] = pts[:, 0] + 1.0
    mesh.point_data["source"] = np.zeros(len(pts), dtype=np.uint8)
    line = pv.lines_from_points(np.column_stack((np.linspace(0, 1, 30), np.linspace(0, 1, 30), np.full(30, 0.3))))
    marker = pv.PolyData(np.array([[0.1, 0.1, 0.3], [0.9, 0.9, 0.3]]))
    actors = {"surface": plotter.add_mesh(mesh, color="lightblue", opacity=0.6, show_edges=True),
              "zweiphasengebiet": None, "isothermen": [plotter.add_mesh(line, color="red")],
              "critical_point": plotter.add_mesh(marker, color="black", point_size=10)}
    return plotter, actors, mesh


class ActorRegistry(unittest.TestCase):
    def test_clip_restore_and_repeat_leave_raw_data_untouched(self):
        plotter, actors, mesh = make_actor_scene()
        raw_before = (mesh.points.copy(), mesh.faces.copy(), mesh.point_data["T"].copy())
        registry = ClipRegistry()
        registry.register(actors)
        self.assertEqual(len(registry.entries), 3)
        stats = registry.apply((0.2, 0.6, 0.2, 0.6, 0.0, 1.0))
        self.assertEqual((stats["actors"], stats["empty"]), (3, 1))        # nur der Marker liegt komplett ausserhalb
        surface = actors["surface"].GetMapper().GetInput()
        pts = np.array([surface.GetPoint(i) for i in range(surface.GetNumberOfPoints())])
        self.assertTrue(np.all(pts[:, 0] >= 0.2 - 1e-12) and np.all(pts[:, 0] <= 0.6 + 1e-12))
        self.assertIsNotNone(surface.GetCellData().GetArray(ORIG_CELL))
        # Marker (Vorbild: kritischer Punkt) ausserhalb der Box wird ausgeblendet, nicht geloescht
        self.assertEqual(actors["critical_point"].GetMapper().GetInput().GetNumberOfPoints(), 0)
        for box in ((0.4, 0.5, 0.4, 0.5, 0, 1), (0.0, 1.0, 0.0, 1.0, 0, 1), (0.1, 0.9, 0.1, 0.9, 0, 1)):
            registry.apply(box)                                  # mehrfache aufeinanderfolgende Bereichswechsel
        registry.apply(None)
        restored = actors["surface"].GetMapper().GetInput()
        self.assertEqual(restored.GetNumberOfCells(), mesh.n_cells)       # vollstaendig originale Darstellung
        np.testing.assert_array_equal(mesh.points, raw_before[0])
        np.testing.assert_array_equal(mesh.faces, raw_before[1])
        np.testing.assert_array_equal(mesh.point_data["T"], raw_before[2])
        self.assertEqual(actors["critical_point"].GetMapper().GetInput().GetNumberOfPoints(), 2)

    def test_new_points_are_flagged_as_render_interpolation_not_heos_nodes(self):
        plotter, actors, mesh = make_actor_scene()
        registry = ClipRegistry()
        registry.register(actors)
        registry.apply((0.33, 0.61, 0.0, 1.0, 0.0, 1.0),
                       lambda xyz: {"T": xyz[:, 0] * 100.0, "v": xyz[:, 1], "p_raw": xyz[:, 0] + 1.0})
        data = pv.wrap(actors["surface"].GetMapper().GetInput())
        source = np.asarray(data.point_data["source"])
        self.assertTrue((source == SRC_CLIPPED).any() and (source == 0).any())
        new = source == SRC_CLIPPED
        np.testing.assert_allclose(data.points[new, 0], np.array([0.33, 0.61])[np.abs(data.points[new, 0][:, None] -
                                   np.array([0.33, 0.61])).argmin(axis=1)], atol=1e-12)
        np.testing.assert_allclose(data.point_data["T"], data.points[:, 0] * 100.0, atol=1e-9)
        # Alle echten Originalknoten bleiben exakt erhalten (HEOS-Knoten werden nie neu erzeugt)
        originals = {tuple(p) for p in mesh.points}
        for point, flag in zip(data.points, source):
            if flag == 0:
                self.assertIn(tuple(point), originals)

    def test_hook_without_point_data_and_diagrams_without_arrays(self):
        plotter, actors, mesh = make_actor_scene()
        registry = ClipRegistry()
        registry.register(actors)
        stats = registry.apply((0.33, 0.61, 0.0, 1.0, 0.0, 1.0), lambda xyz: None)   # nicht-pvt-Diagramme: kein Hook-Ergebnis
        self.assertEqual(stats["empty"], 1)
        bare = pv.PolyData(mesh.points, mesh.faces)                                     # Netz ohne Punktdaten
        actor = plotter.add_mesh(bare)
        registry.register({"surface": [actor]})
        registry.apply((0.1, 0.5, 0.1, 0.5, 0.0, 1.0))
        self.assertGreater(actor.GetMapper().GetInput().GetNumberOfCells(), 0)

    def test_pick_geometry_refers_to_the_original_triangle(self):
        plotter, actors, mesh = make_actor_scene()
        registry = ClipRegistry()
        registry.register(actors)
        registry.apply((0.33, 0.61, 0.0, 1.0, 0.0, 1.0))
        surface = actors["surface"].GetMapper().GetInput()
        orig = surface.GetCellData().GetArray(ORIG_CELL)
        for cell in range(0, surface.GetNumberOfCells(), 37):
            base = registry.base_triangle(actors["surface"], int(orig.GetValue(cell)))
            clipped = [surface.GetPoint(surface.GetCell(cell).GetPointId(k)) for k in range(3)]
            lo, hi = base.min(axis=0) - 1e-12, base.max(axis=0) + 1e-12
            self.assertTrue(np.all(np.array(clipped) >= lo) and np.all(np.array(clipped) <= hi))
        self.assertIsNone(registry.base_triangle(actors["surface"], 10 ** 9))

    def test_V04_empty_window_and_axes_without_geometry(self):
        plotter, actors, _mesh = make_actor_scene()
        registry = ClipRegistry()
        registry.register(actors)
        stats = registry.apply((5, 6, 5, 6, 5, 6))
        self.assertEqual(stats["empty"], stats["actors"])
        for key in ("surface", "isothermen", "critical_point"):
            for actor in (actors[key] if isinstance(actors[key], list) else [actors[key]]):
                self.assertEqual(actor.GetMapper().GetInput().GetNumberOfPoints(), 0)
        d = daten_for("pvt")
        box = A.validate_bounds("pvt", d, BOUNDS["pvt"])[1]
        used = axes_view.apply_axes(plotter, "pvt", d, actors, box=box)           # Rahmen/Achsen trotz leerem Fenster
        self.assertEqual(used, box)
        plotter.render()
        self.assertIsNone(axes_view.visible_bounds({k: v for k, v in actors.items()}))
        registry.apply(None)
        self.assertGreater(actors["surface"].GetMapper().GetInput().GetNumberOfCells(), 0)

    def test_axes_replace_instead_of_duplicate_on_repeated_calls(self):
        plotter, actors, _mesh = make_actor_scene()
        d = daten_for("pvt")
        box = A.validate_bounds("pvt", d, BOUNDS["pvt"])[1]
        for _ in range(3):
            axes_view.apply_axes(plotter, "pvt", d, actors, box=box)
        names = [name for name in plotter.renderer.actors if str(name).startswith("ts_axes")]
        self.assertEqual(sorted(names), sorted(set(names)))
        self.assertLessEqual(len(names), 3)


if __name__ == "__main__":
    unittest.main()
