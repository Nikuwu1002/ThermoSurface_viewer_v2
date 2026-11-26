"""TS-VP04..VP07: Ansichtsbereich und Isolinien in der GUI (Testmatrix V02, V05-V07, V11, C01, C02).

Wie test_gui_multiview.py: braucht ein echtes (OpenGL-faehiges) Qt-Fenster, CoolProp und PyVistaQt; ohne diese
werden die Tests uebersprungen.
"""
import os
import time
import unittest
from unittest.mock import patch

import numpy as np

GUI_SKIP = None
try:
    os.environ.setdefault("QT_API", "pyside6")
    import main as app_main
except Exception as exc:  # pragma: no cover
    GUI_SKIP = f"main.py nicht importierbar: {exc}"
if GUI_SKIP is None and os.environ.get("QT_QPA_PLATFORM", "").lower() in ("offscreen", "minimal"):
    GUI_SKIP = "GUI-Tests brauchen eine OpenGL-faehige Plattform (nicht offscreen/minimal)"
if GUI_SKIP is None:
    try:
        import CoolProp  # noqa: F401
        import pyvistaqt  # noqa: F401
        from PySide6.QtWidgets import QApplication, QMessageBox
        import diagram_axes as axes
    except Exception as exc:  # pragma: no cover
        GUI_SKIP = f"GUI-Abhaengigkeiten fehlen: {exc}"


def spin(app, condition, timeout=240.0):
    end = time.time() + timeout
    while not condition():
        app.processEvents()
        time.sleep(0.01)
        if time.time() > end:
            raise TimeoutError("Bedingung nicht rechtzeitig erfuellt")


def mapper_points(actor):
    data = actor.GetMapper().GetInput()
    return np.array([data.GetPoint(i) for i in range(data.GetNumberOfPoints())]).reshape(-1, 3)


def mapper_cells(actor):
    return actor.GetMapper().GetInput().GetNumberOfCells()


def all_actors(viewport, keys=("surface", "zweiphasengebiet", "isothermen", "isochoren", "isobaren", "binodale")):
    out = []
    for key in keys:
        objects = viewport.actors.get(key)
        for actor in (objects if isinstance(objects, (list, tuple)) else [objects]):
            if actor is not None:
                out.append(actor)
    return out


# Gueltige Grenzen je Diagramm fuer Wasser (Anzeigeeinheiten x, y, z); jeweils ein echter Teilausschnitt.
BOUNDS_WATER = {
    "uvt": [(2e-5, 3e-4), (500.0, 700.0), (-200.0, 2800.0)],
    "pht": [(-1000.0, 1500.0), (500.0, 700.0), (0.5, 30.0)],
    "hsp": [(-1000.0, 1500.0), (-3.0, 3.0), (0.5, 30.0)],
}


@unittest.skipIf(GUI_SKIP, GUI_SKIP or "")
class ViewportGuiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def make_window(self):
        stoff, raster, visu = app_main.lade_konfiguration()
        self.calls = []

        def berechnen(stoff_, raster_, visu_, diagram_type="pvt", progress_callback=None, isolines=None):
            self.calls.append((diagram_type, None if isolines is None else dict(isolines)))
            return app_main.berechne_auswahl(stoff_, raster_, visu_, diagram_type=diagram_type,
                                             progress_callback=progress_callback, isolines=isolines)

        window = app_main.ThermoSurfaceWindow(
            start_stoff=stoff, surface_cfg=raster["isotherm_grid"], raster=raster, visu_raster=visu,
            berechnen=berechnen, grid_generator=app_main.generate_stoff_grid,
            isoline_defaults=app_main.lade_isolinien_vorgaben())
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

    def single(self, window, diagram_id):
        window.diagram_combo.setCurrentIndex(window.diagram_combo.findData(diagram_id))
        window.on_confirm_stoff()
        self.calculate(window)
        return window.workspace.viewports[diagram_id]

    def multi(self, window, ids):
        window.mode_combo.setCurrentIndex(window.mode_combo.findData("multi"))
        for diagram_id, box in window.diagram_checks.items():
            box.setChecked(diagram_id in ids)
        window.on_confirm_stoff()
        self.calculate(window)
        return window.workspace.viewports

    def enter(self, window, bounds):
        window.view_check.setChecked(True)
        for edit, value in zip(window.view_fields, [x for pair in bounds for x in pair]):
            edit.setText(repr(value))

    # ------------------------------------------------------------------ V01/V02
    def test_controls_follow_the_axis_matrix_and_start_disabled(self):
        window = self.make_window()
        viewport = self.single(window, "uvt")
        position, focal, up = viewport.plotter.camera_position
        self.assertLess(position[1], focal[1])
        self.assertGreater(position[2], focal[2])
        np.testing.assert_allclose(up, (0.0, 0.0, 1.0))
        self.assertFalse(window.view_check.isChecked())                           # zunaechst deaktiviert
        self.assertTrue(all(not e.isEnabled() for e in window.view_fields))
        self.assertFalse(window.view_apply_button.isEnabled())
        self.assertFalse(viewport.view_state.enabled)
        texts = [label.text() for label in window.view_labels]
        self.assertEqual(texts, ["x: v von [m³/mol]", "x: v bis [m³/mol]", "y: T von [K]", "y: T bis [K]",
                                 "z: u von [kJ/kg]", "z: u bis [kJ/kg]"])
        window.workspace.set_active(viewport)
        window.on_diagram_changed()
        self.assertEqual(viewport.manual_box(), None)                              # V01: kein Clipping

    def test_V02_apply_clips_original_geometry_without_any_recomputation(self):
        window = self.make_window()
        viewport = self.single(window, "uvt")
        calls = list(self.calls)
        key = viewport.result_key
        daten = viewport.daten
        raw = {k: np.array(daten[k], copy=True) for k in ("surface_liquid", "surface_vapor", "surface_high_temperature")}
        stats_before = window.cache.stats()
        original_cells = [mapper_cells(a) for a in all_actors(viewport)]

        self.enter(window, BOUNDS_WATER["uvt"])
        self.assertTrue(all(e.isEnabled() for e in window.view_fields))
        self.assertTrue(window.view_apply_button.isEnabled())
        # Eingabe allein aendert nichts (nur "Ansicht anwenden" wirkt)
        self.assertFalse(viewport.view_state.enabled)
        self.assertEqual([mapper_cells(a) for a in all_actors(viewport)], original_cells)

        window.view_apply_button.click()
        self.app.processEvents()
        _si, box = axes.validate_bounds("uvt", daten, BOUNDS_WATER["uvt"])
        self.assertTrue(viewport.view_state.enabled)
        self.assertEqual(viewport.manual_box(), box)
        lo, hi = np.array(box[0::2]), np.array(box[1::2])
        for actor in all_actors(viewport):
            pts = mapper_points(actor)
            if len(pts):
                self.assertTrue(np.all(pts >= lo - 1e-9) and np.all(pts <= hi + 1e-9))
        self.assertLess(sum(mapper_cells(a) for a in all_actors(viewport)), sum(original_cells))
        # Achsen kommen aus den sechs Eingaben (nicht aus den Actor-Bounds)
        used = viewport.plotter.renderer.actors["ts_axes_frame"].GetBounds()
        np.testing.assert_allclose(used, box, atol=1e-9)

        # C01: keine neue HEOS-Berechnung, identisches Ergebnisobjekt, identische Rohdaten und Cache-ID
        self.assertEqual(self.calls, calls)
        self.assertIs(viewport.daten, daten)
        self.assertEqual(viewport.result_key, key)
        for name, before in raw.items():
            np.testing.assert_array_equal(daten[name], before)
        self.assertEqual(window.cache.stats()["entries"], stats_before["entries"])
        self.assertEqual(window.cache.stats()["misses"], stats_before["misses"])

    def test_invalid_input_reports_error_and_keeps_last_valid_view(self):
        window = self.make_window()
        viewport = self.single(window, "uvt")
        self.enter(window, BOUNDS_WATER["uvt"])
        window.view_apply_button.click()
        box = viewport.manual_box()
        cells = [mapper_cells(a) for a in all_actors(viewport)]
        for bad in (["abc"], ["-1"], ["1e-3"]):
            window.view_fields[0].setText(bad[0])                  # Start v: Text / negativ / groesser als Stop
            window.view_apply_button.click()
            self.assertIn("color", window.view_message.styleSheet())
            self.assertTrue(window.view_message.text())
            self.assertEqual(viewport.manual_box(), box)
            self.assertEqual([mapper_cells(a) for a in all_actors(viewport)], cells)
        self.assertEqual(len(self.calls), 1)

    def test_V03_reset_restores_full_original_geometry(self):
        window = self.make_window()
        viewport = self.single(window, "pht")
        original = [mapper_cells(a) for a in all_actors(viewport)]
        points = [mapper_points(a).copy() for a in all_actors(viewport)]
        self.enter(window, BOUNDS_WATER["pht"])
        window.view_apply_button.click()
        narrow = [mapper_cells(a) for a in all_actors(viewport)]
        self.assertLess(sum(narrow), sum(original))
        # Grenzen wieder erweitern und danach zuruecksetzen
        bounds = [(-5000.0, 8000.0), (300.0, 1000.0), (0.001, 300.0)]
        self.enter(window, bounds)
        window.view_apply_button.click()
        self.assertEqual(sum(mapper_cells(a) for a in all_actors(viewport)), sum(original))
        window.view_reset_button.click()
        self.assertFalse(viewport.view_state.enabled)
        self.assertIsNone(viewport.manual_box())
        self.assertEqual([mapper_cells(a) for a in all_actors(viewport)], original)
        for actor, before in zip(all_actors(viewport), points):
            np.testing.assert_array_equal(mapper_points(actor), before)
        self.assertFalse(window.view_check.isChecked())

    def test_V04_empty_window_keeps_frame_and_axes(self):
        window = self.make_window()
        viewport = self.single(window, "pht")
        messages = []
        viewport.status_message.connect(messages.append)
        self.enter(window, [(2.5e3, 2.6e3), (900.0, 950.0), (150.0, 160.0)])
        window.view_apply_button.click()
        self.assertTrue(all(mapper_points(a).shape[0] == 0 for a in all_actors(viewport)))
        self.assertIn("ts_axes_frame", viewport.plotter.renderer.actors)
        self.assertGreaterEqual(viewport.plotter.renderer.GetActors2D().GetNumberOfItems(), 2)    # Ticks und Titel
        self.assertTrue(any("keine Geometrie" in m for m in messages))
        self.assertTrue(viewport.scene_is_ready())
        self.calculate(window)
        self.assertIs(window.workspace.viewports["pht"], viewport)
        window.on_view_rebuild()
        rebuilt = window.workspace.viewports["pht"]
        self.assertTrue(rebuilt.scene_is_ready())
        self.assertTrue(rebuilt.view_state.enabled)
        self.assertTrue(all(mapper_points(a).shape[0] == 0 for a in all_actors(rebuilt)))

    # ------------------------------------------------------------------ V05
    def test_V05_layer_toggles_keep_box_and_camera(self):
        window = self.make_window()
        viewport = self.single(window, "uvt")
        window.workspace.set_active(viewport)
        self.enter(window, BOUNDS_WATER["uvt"])
        window.view_apply_button.click()
        box = viewport.manual_box()
        camera = viewport.camera_state()
        frame = viewport.plotter.renderer.actors["ts_axes_frame"].GetBounds()
        for key in ("binodale", "isothermen", "surface", "zweiphasengebiet", "critical_point"):
            checkbox = window.layer_checks[key]
            checkbox.setChecked(not checkbox.isChecked())
            self.app.processEvents()
            self.assertEqual(viewport.manual_box(), box)
            self.assertEqual(viewport.camera_state(), camera)
            self.assertEqual(viewport.plotter.renderer.actors["ts_axes_frame"].GetBounds(), frame)
        window.grid_checkbox.setChecked(False)
        self.assertEqual(viewport.camera_state(), camera)
        self.assertEqual(viewport.manual_box(), box)

    # ------------------------------------------------------------------ V06/V07
    def test_V06_V07_independent_viewports_survive_maximize_and_mode_switch(self):
        window = self.make_window()
        viewports = self.multi(window, ["pht", "uvt"])
        workspace = window.workspace
        workspace.set_active(viewports["pht"])
        self.enter(window, BOUNDS_WATER["pht"])
        window.view_apply_button.click()
        workspace.set_active(viewports["uvt"])
        self.enter(window, BOUNDS_WATER["uvt"])
        window.view_apply_button.click()
        boxes = {i: v.manual_box() for i, v in viewports.items()}
        self.assertNotEqual(boxes["pht"], boxes["uvt"])
        self.assertTrue(all(b is not None for b in boxes.values()))
        calls = len(self.calls)

        workspace.maximize(viewports["uvt"])                               # Maximieren / Entmaximieren
        self.app.processEvents()
        workspace.restore()
        self.assertEqual({i: v.manual_box() for i, v in viewports.items()}, boxes)
        workspace.set_active(viewports["pht"])                             # GUI zeigt den Zustand des aktiven Diagramms
        self.assertTrue(window.view_check.isChecked())
        self.assertEqual(window.view_fields[0].text(), repr(BOUNDS_WATER["pht"][0][0]))

        # Wechsel Mehrfach -> Einzel -> Mehrfach: neu gebildete Viewports desselben Typs uebernehmen die Einstellung
        window.mode_combo.setCurrentIndex(window.mode_combo.findData("single"))
        window.diagram_combo.setCurrentIndex(window.diagram_combo.findData("uvt"))
        self.calculate(window)
        self.assertEqual(window.workspace.viewports["uvt"].manual_box(), boxes["uvt"])
        window.mode_combo.setCurrentIndex(window.mode_combo.findData("multi"))
        for diagram_id, check in window.diagram_checks.items():
            check.setChecked(diagram_id in ("pht", "uvt"))
        self.calculate(window)
        self.assertEqual({i: v.manual_box() for i, v in window.workspace.viewports.items()}, boxes)
        self.assertEqual(len(self.calls), calls)                           # alles aus dem Cache: keine HEOS-Berechnung

    def test_markers_outside_the_box_are_hidden_but_kept(self):
        window = self.make_window()
        viewport = self.single(window, "uvt")
        _si, box = axes.validate_bounds("uvt", viewport.daten, BOUNDS_WATER["uvt"])
        inside = tuple(0.5 * (box[2 * i] + box[2 * i + 1]) for i in range(3))
        outside = (box[1] + 1.0, inside[1], inside[2])
        viewport.set_selection_markers([(1, inside, "red"), (2, outside, "blue")])
        self.assertEqual((viewport.marker_count, len(viewport._marker_actors)), (2, 4))
        self.enter(window, BOUNDS_WATER["uvt"])
        window.view_apply_button.click()
        self.assertEqual((viewport.marker_count, len(viewport._marker_actors)), (2, 2))     # nur der innere sichtbar
        window.view_reset_button.click()
        self.assertEqual(len(viewport._marker_actors), 4)

    # ------------------------------------------------------------------ V11
    def test_V11_picking_on_a_clipped_surface_uses_the_original_triangle(self):
        from state_picker import pick_from_viewport
        import vtkmodules.vtkRenderingCore as vtkrc
        window = self.make_window()
        viewport = self.single(window, "uvt")
        window.compare_check.setChecked(True)
        for key in ("zweiphasengebiet", "binodale", "critical_point", "isothermen", "isochoren", "isobaren"):
            viewport.layer_state[key] = False
        viewport.apply_layers()
        self.enter(window, BOUNDS_WATER["uvt"])
        window.view_apply_button.click()
        index = window.state_store.ensure("Water", "uvt", viewport.daten, viewport.result_key)
        plotter = viewport.plotter
        picked = 0
        for actor in viewport.actors["surface"]:
            data = actor.GetMapper().GetInput()
            orig = data.GetCellData().GetArray("ts_orig_cell")
            if data.GetNumberOfCells() == 0:
                continue
            self.assertIsNotNone(orig)
            for cell_id in range(0, data.GetNumberOfCells(), max(1, data.GetNumberOfCells() // 6)):
                cell = data.GetCell(cell_id)
                centroid = np.mean([data.GetPoint(cell.GetPointId(k)) for k in range(3)], axis=0)
                coordinate = vtkrc.vtkCoordinate()
                coordinate.SetCoordinateSystemToWorld()
                coordinate.SetValue(*centroid)
                x, y = coordinate.GetComputedDoubleDisplayValue(plotter.renderer)
                try:
                    pick = pick_from_viewport(viewport, index, int(round(x)), int(round(y)))
                except Exception:
                    continue                                         # verdeckt/ausserhalb: kein gueltiger Treffer
                picked += 1
                self.assertTrue(pick.clipped)
                self.assertTrue(pick.fully_located)
                self.assertAlmostEqual(sum(pick.weights), 1.0, places=9)
                state = window.resolver.resolve(pick, index)
                self.assertTrue(state.provenance.get("viewport_clipped"))
        self.assertGreater(picked, 0)

    # ------------------------------------------------------------------ Isolinien in der GUI (C02)
    def test_C02_isoline_changes_recompute_only_the_affected_diagram(self):
        window = self.make_window()
        viewport = self.single(window, "uvt")
        self.assertEqual(self.calls, [("uvt", None)])
        default_levels = viewport.daten["isoline_report"]["isothermen"]["n_levels"]
        # unveraenderter Faktor: Cache-Treffer
        self.calculate(window)
        self.assertEqual(len(self.calls), 1)
        # Faktor aendern -> neue Isolinien-Berechnung mit anderem Cache-Schluessel
        window.isoline_fields["isothermen"]["factor"].setValue(1.1)
        self.calculate(window)
        self.assertEqual(len(self.calls), 2)
        self.assertEqual(self.calls[-1][1]["isothermen"]["factor"], 1.1)
        viewport = window.workspace.viewports["uvt"]
        self.assertGreater(viewport.daten["isoline_report"]["isothermen"]["n_levels"], default_levels)
        self.assertIn("Isothermen", window.isoline_report_label.text())
        # Zurueck zum Faktor 3.14 -> Cache-Treffer; Viewportgrenzen aendern nie den Cache
        window.isoline_fields["isothermen"]["factor"].setValue(3.14)
        self.calculate(window)
        self.assertEqual(len(self.calls), 2)                              # identische Parameter: Cache-Treffer

    def test_invalid_isoline_input_is_rejected_before_calculation(self):
        window = self.make_window()
        window.diagram_combo.setCurrentIndex(window.diagram_combo.findData("uvt"))
        window.on_confirm_stoff()
        window.isoline_fields["isobaren"]["start"].setValue(5.0)       # 5 MPa
        window.isoline_fields["isobaren"]["stop"].setValue(1.0)        # < start
        with patch.object(QMessageBox, "warning") as warning:
            window.on_calculate_clicked()
        warning.assert_called_once()
        self.assertFalse(window.manager.is_running)
        self.assertEqual(self.calls, [])

    def test_isoline_units_are_converted_to_si(self):
        window = self.make_window()
        window.diagram_combo.setCurrentIndex(window.diagram_combo.findData("gpt"))
        window.isoline_fields["isobaren"]["start"].setValue(2.0)       # MPa
        window.isoline_fields["isobaren"]["stop"].setValue(20.0)
        window.isoline_fields["isochoren"]["start"].setValue(5e-5)     # m3/mol
        config = window.isoline_config_si("gpt")
        self.assertEqual((config["isobaren"]["start"], config["isobaren"]["stop"]), (2e6, 2e7))
        self.assertEqual(config["isochoren"]["start"], 5e-5)
        self.assertIsNone(config["isothermen"]["start"])               # 0 = automatisch


if __name__ == "__main__":
    unittest.main()
