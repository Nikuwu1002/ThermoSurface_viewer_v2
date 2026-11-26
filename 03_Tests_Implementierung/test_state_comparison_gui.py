"""TS-B07 (GUI): Punktauswahl per echtem Mausklick, Marker, Vergleichsmodus, Loeschen/Reihenfolge, Stoffwechsel.

Braucht ein echtes OpenGL-Fenster (nicht QT_QPA_PLATFORM=offscreen) und wird sonst uebersprungen.
"""
import time
import unittest
from unittest.mock import patch

from test_gui_multiview import GUI_SKIP, app_main, spin

if not GUI_SKIP:
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QApplication, QMessageBox

    from state_picker import PickError, pick_from_viewport


def pump(app, seconds):
    end = time.time() + seconds
    while time.time() < end:
        app.processEvents()
        time.sleep(0.01)


@unittest.skipIf(GUI_SKIP, GUI_SKIP or "")
class B07GuiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def make_window(self):
        stoff, raster, visu = app_main.lade_konfiguration()
        window = app_main.ThermoSurfaceWindow(
            start_stoff=stoff, surface_cfg=raster["isotherm_grid"], raster=raster, visu_raster=visu,
            berechnen=app_main.berechne_auswahl, grid_generator=app_main.generate_stoff_grid)
        window.show()
        self.addCleanup(self.close_window, window)
        return window

    def close_window(self, window):
        window.manager.shutdown()
        window.close()
        self.app.processEvents()

    def calculate(self, window, ids):
        window.mode_combo.setCurrentIndex(window.mode_combo.findData("multi"))
        for diagram_id, box in window.diagram_checks.items():
            box.setChecked(diagram_id in ids)
        window.on_calculate_clicked()
        spin(self.app, lambda: not window.manager.is_running)
        self.app.processEvents()

    def find_hit(self, window, viewport, preferred=None):
        """Qt-Pixel, an dem der VTK-Picker eine echte Zustandsflaeche trifft (gleiche Umrechnung wie der Viewport)."""
        key = viewport.result_key
        index = window.state_store.ensure(key.fluid, viewport.diagram_id, viewport.daten, key)
        width, height = viewport.plotter.width(), viewport.plotter.height()
        ratio = viewport.plotter.devicePixelRatioF()
        for fy in (0.5, 0.4, 0.6, 0.3, 0.7):
            for fx in (0.5, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8):
                qx, qy = int(width * fx), int(height * fy)
                try:
                    pick = pick_from_viewport(viewport, index, int(round(qx * ratio)),
                                              int(round((height - 1 - qy) * ratio)))
                except PickError:
                    continue
                if preferred is None or pick.actor_key in preferred:
                    return qx, qy
        self.fail("Kein pickbarer Pixel gefunden")

    def click(self, viewport, qx, qy, wait=0.8):
        QTest.mouseClick(viewport.plotter, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, QPoint(qx, qy))
        pump(self.app, wait)   # Einfachklick wird erst nach dem Doppelklick-Intervall gemeldet

    def test_picking_inactive_until_mode_enabled(self):
        # Vergleichsmodus aus: keine Auswertung
        window = self.make_window()
        window.on_confirm_stoff()
        self.calculate(window, ["pht"])
        viewport = window.workspace.viewports["pht"]
        qx, qy = self.find_hit(window, viewport)
        before = window.state_store.evaluations
        self.click(viewport, qx, qy)
        self.assertEqual(len(window.selection), 0)
        self.assertEqual(window.state_store.evaluations, before)

    def test_surface_pick_is_not_blocked_by_a_nearby_isoline(self):
        from vtkmodules.vtkRenderingCore import vtkCellPicker

        from state_picker import find_actor

        window = self.make_window()
        window.on_confirm_stoff()
        self.calculate(window, ["pvt"])
        viewport = window.workspace.viewports["pvt"]
        index = window.state_store.ensure(viewport.result_key.fluid, "pvt", viewport.daten, viewport.result_key)
        exact, tolerant = vtkCellPicker(), vtkCellPicker()
        exact.SetTolerance(0.0)
        tolerant.SetTolerance(0.005)
        width, height = viewport.plotter.render_window.GetSize()
        for y in range(height // 4, 3 * height // 4, 4):
            for x in range(width // 4, 3 * width // 4, 4):
                if not exact.Pick(x, y, 0, viewport.plotter.renderer):
                    continue
                actor_key, _ = find_actor(viewport.actors, exact.GetActor(), exact.GetMapper())
                if actor_key != "surface" or not tolerant.Pick(x, y, 0, viewport.plotter.renderer):
                    continue
                other_key, _ = find_actor(viewport.actors, tolerant.GetActor(), tolerant.GetMapper())
                if other_key not in ("isochoren", "isothermen", "binodale"):
                    continue
                pick = pick_from_viewport(viewport, index, x, y)
                self.assertEqual(pick.actor_key, "surface")
                self.assertEqual(pick.cell_id, exact.GetCellId())
                return
        self.fail("Keinen sichtbaren Oberflaechentreffer nahe einer Isolinie gefunden")

    def test_picker_rejects_holes_background_hidden_and_occluded_surfaces(self):
        from types import SimpleNamespace

        import numpy as np
        import pyvista as pv
        from vtkmodules.vtkRenderingCore import vtkCellPicker, vtkCoordinate

        from diagram_coords import from_xyz_pht
        from state_store import StateStore

        x_values = (-1.0, -0.05, 0.05, 1.0)
        y_values = (0.1, 1.05, 1.15, 2.1)
        points = np.array([[(x, y, 0.0) for x in x_values] for y in y_values])
        faces = []
        for row in range(3):
            for col in range(3):
                if (row, col) == (1, 1):
                    continue
                a = 4 * row + col
                faces.extend((3, a, a + 1, a + 4, 3, a + 1, a + 5, a + 4))
        mesh = pv.PolyData(points.reshape(-1, 3), np.array(faces))
        critical = (647.0, 22e6, 2e6)
        index = StateStore().register("Water", "pht", {
            "critical_point": critical,
            "surface_liquid": from_xyz_pht(points, critical),
        })
        plotter = pv.Plotter(off_screen=True, window_size=(600, 600))
        self.addCleanup(plotter.close)
        actor = plotter.add_mesh(mesh)
        plotter.camera_position = [(0, 1.1, 5), (0, 1.1, 0), (0, 1, 0)]
        plotter.enable_parallel_projection()
        plotter.camera.parallel_scale = 1.2
        plotter.show(auto_close=False)
        viewport = SimpleNamespace(plotter=plotter, actors={"surface": [actor]})
        coordinate = vtkCoordinate()
        coordinate.SetCoordinateSystemToWorld()
        coordinate.SetValue(0.04, 1.1, 0.0)
        x, y = coordinate.GetComputedDisplayValue(plotter.renderer)
        exact = vtkCellPicker()
        exact.SetTolerance(0.0)
        self.assertFalse(exact.Pick(x, y, 0, plotter.renderer))
        with self.assertRaises(PickError):
            pick_from_viewport(viewport, index, x, y)

        coordinate.SetValue(1.01, 1.1, 0.0)
        x, y = coordinate.GetComputedDisplayValue(plotter.renderer)
        self.assertFalse(exact.Pick(x, y, 0, plotter.renderer))
        pick = pick_from_viewport(viewport, index, x, y)
        self.assertEqual(pick.actor_key, "surface")
        self.assertTrue(pick.fully_located)
        self.assertIn("fallback", pick.message)

        blocker = plotter.add_mesh(pv.Plane(center=(0.985, 1.1, 0.2), i_size=0.03, j_size=0.12),
                                   pickable=False)
        plotter.render()
        self.assertEqual(plotter.renderer.GetZ(x, y), 1.0)
        with self.assertRaises(PickError):
            pick_from_viewport(viewport, index, x, y)
        plotter.remove_actor(blocker)
        plotter.render()

        coordinate.SetValue(1.1, 1.1, 0.0)
        x, y = coordinate.GetComputedDisplayValue(plotter.renderer)
        self.assertFalse(exact.Pick(x, y, 0, plotter.renderer))
        with self.assertRaises(PickError):
            pick_from_viewport(viewport, index, x, y)

        coordinate.SetValue(0.5, 1.1, 0.0)
        x, y = coordinate.GetComputedDisplayValue(plotter.renderer)
        actor.SetVisibility(False)
        plotter.render()
        with self.assertRaises(PickError):
            pick_from_viewport(viewport, index, x, y)
        actor.SetVisibility(True)
        foreground = plotter.add_mesh(pv.Plane(center=(0.5, 1.1, 0.2), i_size=0.3, j_size=0.3))
        viewport.actors["other"] = foreground
        plotter.render()
        self.assertTrue(exact.Pick(x, y, 0, plotter.renderer))
        self.assertEqual(exact.GetActor(), foreground)
        with self.assertRaises(PickError):
            pick_from_viewport(viewport, index, x, y)

        foreground.SetPickable(False)
        plotter.render()
        self.assertTrue(exact.Pick(x, y, 0, plotter.renderer))
        self.assertEqual(exact.GetActor(), actor)
        with self.assertRaises(PickError):
            pick_from_viewport(viewport, index, x, y)

        from vtkmodules.vtkCommonDataModel import vtkCellArray
        vertices = vtkCellArray()
        vertices.InsertNextCell(1)
        vertices.InsertCellPoint(0)
        foreground.GetMapper().GetInput().SetVerts(vertices)
        plotter.render()
        self.assertEqual(foreground.GetMapper().GetInput().GetCell(0).GetCellDimension(), 0)
        with self.assertRaises(PickError):
            pick_from_viewport(viewport, index, x, y)

        plotter.camera.SetClippingRange(0.01, 100000.0)
        plotter.render()
        self.assertTrue(exact.Pick(x, y, 0, plotter.renderer))
        plotter.renderer.SetWorldPoint(*exact.GetPickPosition(), 1.0)
        plotter.renderer.WorldToDisplay()
        self.assertLess(abs(plotter.renderer.GetZ(x, y) - plotter.renderer.GetDisplayPoint()[2]), 1e-5)
        with self.assertRaises(PickError):
            pick_from_viewport(viewport, index, x, y)

    def test_click_selects_marks_and_double_click_only_maximizes(self):
        window = self.make_window()
        window.on_confirm_stoff()
        self.calculate(window, ["pvt", "hsp"])
        window.compare_check.setChecked(True)
        pvt, hsp = window.workspace.viewports["pvt"], window.workspace.viewports["hsp"]

        qx, qy = self.find_hit(window, hsp)
        self.click(hsp, qx, qy)
        self.assertEqual(len(window.selection), 1)
        self.assertEqual(window.selection.points[0].source_diagram, "hsp")
        self.assertEqual(hsp.marker_count, 1)
        self.assertEqual(pvt.marker_count, 0)
        self.assertEqual(window.compare_list.count(), 1)

        qx, qy = self.find_hit(window, pvt)
        self.click(pvt, qx, qy)
        self.assertEqual(len(window.selection), 2)
        self.assertEqual({p.source_diagram for p in window.selection.points}, {"hsp", "pvt"})
        self.assertEqual((pvt.marker_count, hsp.marker_count), (1, 1))
        self.assertEqual({p.state.fluid for p in window.selection.points}, {window.confirmed_stoff})

        # Doppelklick maximiert und waehlt keinen Punkt
        QTest.mouseDClick(hsp.plotter, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, QPoint(qx, qy))
        pump(self.app, 0.8)
        self.assertEqual(len(window.selection), 2)
        self.assertIs(window.workspace.maximized, hsp)
        window.workspace.restore()
        pump(self.app, 0.8)
        # Vergleichsdiagramm: Achsentausch ohne neue Thermodynamikauswertung
        window.on_open_comparison()
        view = window.comparison_view
        self.assertEqual(len([r for r in view.geometry.rows if r.xyz is not None]), 2)
        view.close()
        self.app.processEvents()
        self.assertIs(window.comparison_view, view)
        self.assertFalse(view.isVisible())
        self.assertFalse(view._disposed)
        window.on_open_comparison()
        self.assertIs(window.comparison_view, view)
        self.assertTrue(view.isVisible())
        evaluations = window.state_store.evaluations
        view.set_axes("T", "p", "v")
        view.set_axes("v", "T", "p")
        self.assertEqual(window.state_store.evaluations, evaluations)
        self.assertEqual(len(set(view.axes)), 3)

        # Reihenfolge, Loeschen (Nummern bleiben eindeutig), Reset
        first, second = [p.selection_index for p in window.selection.points]
        window.compare_list.setCurrentRow(0)
        window.on_compare_move(+1)
        self.assertEqual([p.selection_index for p in window.selection.points], [second, first])
        window.compare_list.setCurrentRow(0)
        window.on_compare_delete()
        self.assertEqual([p.selection_index for p in window.selection.points], [first])
        self.assertEqual(sum(v.marker_count for v in window.workspace.ordered()), 1)
        qx, qy = self.find_hit(window, hsp)
        self.click(hsp, qx, qy)
        self.assertNotIn(window.selection.points[-1].selection_index, (first, second))
        window.on_compare_reset()
        self.assertEqual(len(window.selection), 0)
        self.assertEqual(sum(v.marker_count for v in window.workspace.ordered()), 0)

    def test_cross_diagram_selection_projects_to_target_surface_without_replacing_it(self):
        window = self.make_window()
        window.on_confirm_stoff()
        self.calculate(window, ["pvt", "hsp"])
        window.compare_check.setChecked(True)
        source = window.workspace.viewports["hsp"]
        target = window.workspace.viewports["pvt"]
        qx, qy = self.find_hit(window, source)
        self.click(source, qx, qy)
        original_target_actors = dict(target.actors)

        window.on_toggle_projection()
        if window.selection.points[0].state.phase == "twophase":
            self.assertEqual(target.projected_marker_count, 0)
            self.assertIn("Zweiphasen", window.projection_info.toolTip())
            window.layer_checks["zweiphasengebiet"].setChecked(True)

        self.assertEqual(target.projected_marker_count, 1, window.projection_info.toolTip())
        self.assertEqual(target.marker_count, 0)
        self.assertIs(target.actors["surface"], original_target_actors["surface"])
        self.assertIn("1 von 1", window.projection_info.text())

        selected = window.selection.points[0]
        camera = target.camera_state()
        prop_count = target.plotter.renderer.GetViewProps().GetNumberOfItems()
        for _ in range(25):
            window.selection.reset()
            self.assertEqual(target.projected_marker_count, 0)
            window.selection.add(selected.state, selected.display_xyz, selected.original)
            self.assertEqual(target.projected_marker_count, 1)
            self.assertEqual(target.plotter.renderer.GetViewProps().GetNumberOfItems(), prop_count)
            self.assertEqual(target.camera_state(), camera)
            self.assertIs(target.actors["surface"], original_target_actors["surface"])

        selected_v = window.selection.points[0].state.v_molar
        fields = window._suggest_view_fields(target)
        fields[0] = str(selected_v * 1.1)
        fields[1] = str(selected_v * 1.5)
        window.view_check.setChecked(True)
        for edit, value in zip(window.view_fields, fields):
            edit.setText(value)
        window.on_view_apply()
        self.assertIn("0 von 1", window.projection_info.text())
        self.assertIn("außerhalb", window.projection_info.toolTip())
        self.assertEqual(len(target._marker_actors), 0)

        window.on_view_reset()
        self.assertIn("1 von 1", window.projection_info.text())
        self.assertEqual(len(target._marker_actors), 2)
        window.on_toggle_projection()
        self.assertEqual(target.projected_marker_count, 0)

    def test_old_fluid_viewports_cannot_select_or_project_after_fluid_change(self):
        window = self.make_window()
        window.on_confirm_stoff()
        self.calculate(window, ["pvt", "hsp"])
        window.compare_check.setChecked(True)
        source = window.workspace.viewports["hsp"]
        target = window.workspace.viewports["pvt"]
        qx, qy = self.find_hit(window, source)
        self.click(source, qx, qy)
        old_state = window.selection.points[0].state

        window.stoff_combo.setCurrentIndex(window.stoff_combo.findData("CO2"))
        window.on_confirm_stoff()
        self.assertEqual(len(window.selection), 0)
        self.click(source, qx, qy)
        self.assertEqual(len(window.selection), 0)
        self.assertIn("CO2", window.statusBar().currentMessage())

        window.selection.add(old_state.with_changes(fluid="CO2"), (0.0, 0.0, 0.0))
        window._projection_enabled = True
        window._sync_markers()
        self.assertEqual(target.projected_marker_count, 0)
        self.assertIn("Fluid", window.projection_info.toolTip())

    def test_projection_reports_out_of_range_and_homogeneous_states(self):
        import thermo_state as ts

        window = self.make_window()
        window.on_confirm_stoff()
        self.calculate(window, ["pvt", "hsp"])
        window.compare_check.setChecked(True)
        source = window.workspace.viewports["hsp"]
        target = window.workspace.viewports["pvt"]
        self.click(source, *self.find_hit(window, source))
        state = window.selection.points[0].state
        window.selection.reset()
        window.selection.add(state.with_changes(T=1e7), (0.0, 0.0, 0.0))
        window.selection.add(state.with_changes(stability="unstable", eos_source=ts.SOURCE_HEOS_HOMOGENEOUS),
                             (0.0, 0.0, 0.0))
        window.on_toggle_projection()
        self.assertEqual(target.projected_marker_count, 0)
        self.assertIn("außerhalb", window.projection_info.toolTip())
        self.assertIn("homogen", window.projection_info.toolTip().lower())

    def test_projection_prepares_pvt_grid_for_the_confirmed_fluid(self):
        window = self.make_window()
        window.on_confirm_stoff()
        window.diagram_combo.setCurrentIndex(window.diagram_combo.findData("fvt"))
        window.stoff_combo.setCurrentIndex(window.stoff_combo.findData("CO2"))
        window.on_confirm_stoff()
        window.on_calculate_clicked()
        spin(self.app, lambda: not window.manager.is_running)
        source = window.workspace.viewports["fvt"]
        window.compare_check.setChecked(True)
        self.click(source, *self.find_hit(window, source))
        with patch.object(QMessageBox, "warning") as warning:
            window.on_toggle_projection()
        self.assertFalse(warning.called)
        spin(self.app, lambda: not window.manager.is_running)
        target = window.workspace.viewports["pvt"]
        self.assertEqual(target.result_key.fluid, "CO2")
        self.assertEqual(target.projected_marker_count, 1)

    def test_projection_requested_during_calculation_is_not_lost(self):
        window = self.make_window()
        window.on_confirm_stoff()
        window.diagram_combo.setCurrentIndex(window.diagram_combo.findData("pht"))
        window.on_confirm_stoff()
        window.on_calculate_clicked()
        spin(self.app, lambda: not window.manager.is_running)
        source = window.workspace.viewports["pht"]
        window.compare_check.setChecked(True)
        self.click(source, *self.find_hit(window, source))
        window.on_calculate_clicked()
        self.assertTrue(window.manager.is_running)
        window.on_toggle_projection()
        spin(self.app, lambda: (not window.manager.is_running and "pvt" in window.workspace.viewports
                                and window.workspace.viewports["pvt"].has_result), timeout=30.0)
        if window.selection.points[0].state.phase == "twophase":
            self.assertIn("Zweiphasen", window.projection_info.toolTip())
            window.layer_checks["zweiphasengebiet"].setChecked(True)
        self.assertEqual(window.workspace.viewports["pvt"].projected_marker_count, 1,
                         window.projection_info.toolTip())

    def test_exit_releases_threads_views_and_persists_no_selection(self):
        # B-T07: Verlassen der Anwendung
        import glob
        import os
        import threading
        window = self.make_window()
        window.on_confirm_stoff()
        self.calculate(window, ["pht", "pst"])
        window.compare_check.setChecked(True)
        viewport = window.workspace.viewports["pht"]
        qx, qy = self.find_hit(window, viewport)
        self.click(viewport, qx, qy)
        self.assertEqual(len(window.selection), 1)
        window.on_open_comparison()
        view = window.comparison_view
        viewports = list(window.workspace.ordered())
        files_before = set(glob.glob("cache/*")) if os.path.isdir("cache") else set()

        window.close()
        self.app.processEvents()
        self.assertFalse(window.manager.is_running)
        self.assertIsNone(window.comparison_view)
        self.assertEqual(window.workspace.viewports, {})
        self.assertTrue(all(v.plotter._closed for v in viewports + [view]))
        self.assertFalse(window.state_store.config.disk_enabled)
        files_after = set(glob.glob("cache/*")) if os.path.isdir("cache") else set()
        self.assertEqual(files_before, files_after)   # keine Auswahl/kein Zustand auf der Platte
        self.assertFalse([t for t in threading.enumerate() if t.name.startswith("QThread")])

    def test_fluid_change_clears_selection(self):
        window = self.make_window()
        window.on_confirm_stoff()
        self.calculate(window, ["pht"])
        window.compare_check.setChecked(True)
        viewport = window.workspace.viewports["pht"]
        qx, qy = self.find_hit(window, viewport)
        self.click(viewport, qx, qy)
        self.assertEqual(len(window.selection), 1)
        other = "CO2" if window.confirmed_stoff != "CO2" else "Water"
        window.stoff_combo.setCurrentIndex(window.stoff_combo.findData(other))
        window.on_confirm_stoff()
        self.assertEqual(len(window.selection), 0)
        self.assertEqual(window.compare_list.count(), 0)


if __name__ == "__main__":
    unittest.main()
