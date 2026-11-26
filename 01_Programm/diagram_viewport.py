"""
Datei: diagram_viewport.py
Projekt: ThermoSurface
Kategorie: GUI / 3D-Viewport

Beschreibung:
Verwaltet eine interaktive PyVista-Ansicht für genau einen Diagrammtyp.

Aufgaben:
- Steuert Rendering, Kamera, Layer, Hover und Punktauswahl.
- Wendet manuelle Viewportgrenzen und geometrisches Clipping an.

Abhängigkeiten:
- Projektmodule (direkt): axes_view.py, diagram_axes.py, diagram_coords.py,
  diagram_registry.py, fvt_tangent.py, geometry_clip.py, visualization.py
- Weitere Imports: viewport_bounds.py, viewport_clipping.py, viewport_recovery.py
- Bibliotheken: PySide6, numpy, pyvista, pyvistaqt

Verwendet von:
- Projektmodule: diagram_workspace.py
"""

import logging
import math
import time

import numpy as np
from PySide6.QtCore import QEvent, QObject, QTimer, Qt, Signal
from PySide6.QtWidgets import QApplication, QFrame, QLabel, QVBoxLayout
from pyvistaqt import QtInteractor

import axes_view
import diagram_coords
import geometry_clip
from diagram_axes import ViewportBoundsError, pressure_axis_of
from diagram_registry import DIAGRAMS, LAYER_LABELS
from viewport_bounds import ViewportState
from viewport_clipping import ClipRegistry
from viewport_recovery import DIRTY, DISPOSED, EMPTY, FAILED, READY, RENDERING
from visualization import pick_hover_text

LOGGER = logging.getLogger("thermosurface.viewport")

GRID_ACTOR_KEYS = ("surface", "surface_single_phase", "metastable_liquid_surface", "metastable_vapor_surface",
                   "unstable_surface", "replacement_surface")


class _Header(QLabel):
    clicked = Signal()
    double_clicked = Signal()

    def mousePressEvent(self, event):
        self.clicked.emit()
        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event):
        self.double_clicked.emit()
        super().mouseDoubleClickEvent(event)


class _MouseFilter(QObject):
    """Meldet Doppelklicks und Einfachklicks (ohne Ziehen) auf den Interactor, ohne Ereignisse zu verbrauchen.

    Rotation (Ziehen) und Doppelklick-Maximierung bleiben unberuehrt; ein Einfachklick wird erst gemeldet,
    wenn kein Doppelklick folgt.
    """

    double_clicked = Signal()
    clicked = Signal(int, int)     # Qt-Widgetkoordinaten

    DRAG_TOLERANCE_PX = 4

    def __init__(self, parent=None):
        super().__init__(parent)
        self.click_enabled = False
        self._press = None
        self._pending = None
        self._after_double = False
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self._emit_click)

    def _emit_click(self):
        if self._pending is not None:
            x, y = self._pending
            self._pending = None
            self.clicked.emit(x, y)

    def eventFilter(self, obj, event):
        kind = event.type()
        if kind == QEvent.Type.MouseButtonDblClick and event.button() == Qt.MouseButton.LeftButton:
            self._timer.stop()
            self._pending = None
            self._after_double = True
            self.double_clicked.emit()
        elif kind == QEvent.Type.MouseButtonPress and event.button() == Qt.MouseButton.LeftButton:
            self._press = (event.position().x(), event.position().y())
            self._after_double = False
        elif kind == QEvent.Type.MouseButtonRelease and event.button() == Qt.MouseButton.LeftButton:
            press, self._press = self._press, None
            if self.click_enabled and press is not None and not self._after_double:
                dx, dy = event.position().x() - press[0], event.position().y() - press[1]
                if dx * dx + dy * dy <= self.DRAG_TOLERANCE_PX ** 2:
                    self._pending = (int(event.position().x()), int(event.position().y()))
                    self._timer.start(QApplication.doubleClickInterval())
        return False


class DiagramViewport(QFrame):
    activated = Signal(object)
    maximize_requested = Signal(object)
    status_message = Signal(str)
    point_clicked = Signal(object, int, int)   # (viewport, x, y) in VTK-Displaykoordinaten

    def __init__(self, diagram_id, parent=None, show_header=True):
        super().__init__(parent)
        self.spec = DIAGRAMS.get(diagram_id)
        self.diagram_id = diagram_id
        self.daten = None
        self.actors = {}
        self.result_key = None
        self.layer_state = self.spec.initial_layer_state()
        self.grid_on = True
        self.render_count = 0
        self.info_text = ""
        self.tooltips = {}
        self._hover_installed = False
        self._hover_observer_id = None
        self._last_hover_text = None
        self._last_hover_at = 0.0
        self._axes_observer_installed = False
        self._axes_observer_id = None
        self._disposed = False
        self._render_state = EMPTY
        self.view_state = ViewportState(diagram_id)
        self.clip_registry = ClipRegistry()
        self.default_camera = None
        self._manual_box = None

        self.setObjectName("diagramViewport")
        self.setFrameShape(QFrame.Shape.NoFrame)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(2)

        self.header = _Header(self.spec.label)
        self.header.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.header.setCursor(Qt.CursorShape.PointingHandCursor)
        self.header.setToolTip("Klick: Diagramm aktivieren (Layersteuerung) – Doppelklick: maximieren")
        self.header.clicked.connect(lambda: self.activated.emit(self))
        self.header.double_clicked.connect(lambda: self.maximize_requested.emit(self))
        layout.addWidget(self.header)
        self.header.setVisible(show_header)

        self.plotter = QtInteractor(self)
        layout.addWidget(self.plotter, 1)
        self.plotter.add_camera_orientation_widget(animate=True, n_frames=20)

        self._filter = _MouseFilter(self)
        self._filter.double_clicked.connect(lambda: self.maximize_requested.emit(self))
        self._filter.clicked.connect(self._on_qt_click)
        self.plotter.installEventFilter(self._filter)
        self._markers = []          # [(selection_index, display_xyz, color)]
        self._projected_markers = []
        self._marker_actors = []
        self._tangent_data = None
        self._tangent_show_slice = False
        self._tangent_actors = []
        self.set_active(False)

    # ------------------------------------------------------------------ Zustand
    @property
    def has_result(self):
        return self.daten is not None

    def set_header_visible(self, visible):
        self.header.setVisible(bool(visible))

    def set_active(self, active):
        color = "#2d7ff9" if active else "#8a8a8a"
        width = 2 if active else 1
        self.setStyleSheet(f"#diagramViewport {{ border: {width}px solid {color}; }}")
        self.header.setStyleSheet("font-weight: bold;" + ("color: #2d7ff9;" if active else ""))

    # ------------------------------------------------------------------ Rendern
    def show_result(self, daten, key=None):
        """Rendert `daten` mit dem registrierten Renderer in genau diesen Plotter."""
        if self._disposed:
            raise RuntimeError("Ein freigegebener Viewport kann nicht neu gerendert werden.")
        self._render_state = RENDERING
        LOGGER.debug("event=render_begin diagram=%s viewport=%x key=%s",
                     self.diagram_id, id(self), key)
        try:
            self.actors = self.spec.render(self.plotter, daten)
            self.daten = daten
            self.result_key = key
            self.render_count += 1
            self._manual_box = None
            if self.spec.hover:
                self._install_hover()
            self._install_axes_observer()
            self.clip_registry.register(self.actors)
            self.apply_layers(render=False)
            self.apply_grid(render=False)
            default_bounds = axes_view.auto_bounds(self.diagram_id, self.actors)
            if default_bounds is not None:
                axes_view.set_default_view(self.plotter, default_bounds)
                self.update_axes(render=False)
            self.default_camera = self.camera_state()
            self._marker_actors = []
            self._tangent_actors = []
            self._apply_view(fit_camera=self.view_state.enabled, render=False)
            recovery_camera = getattr(self, "_recovery_camera", None)
            if recovery_camera is not None:
                self.plotter.camera_position = recovery_camera
                self.update_axes(render=False)
                del self._recovery_camera
            if not self._markers:
                self.plotter.render()
            if not self.scene_is_ready():
                raise RuntimeError("Renderer lieferte keine gueltige Zustandsoberflaeche.")
        except Exception:
            self._render_state = FAILED
            LOGGER.exception("event=render_error diagram=%s viewport=%x key=%s",
                             self.diagram_id, id(self), key)
            raise
        self._render_state = READY
        LOGGER.debug("event=render_end diagram=%s viewport=%x actors=%d",
                     self.diagram_id, id(self), self._scene_actor_count())

    @property
    def render_state(self):
        if self._disposed:
            return DISPOSED
        if self._render_state == READY and not self.scene_is_ready():
            return FAILED
        return self._render_state

    def mark_scene_dirty(self):
        if not self._disposed:
            self._render_state = DIRTY

    def _surface_actors(self):
        surfaces = self.actors.get("surface")
        return list(surfaces) if isinstance(surfaces, (list, tuple)) else [surfaces]

    def _scene_actor_count(self):
        return sum(actor is not None for actor in self._surface_actors())

    def scene_is_ready(self):
        if (self._disposed or self.daten is None or self.plotter._closed
                or self.plotter.iren is None or self.plotter.renderer.GetRenderWindow() is None):
            return False
        valid_surface = False
        for actor in self._surface_actors():
            if actor is None:
                continue
            try:
                if not self.plotter.renderer.HasViewProp(actor):
                    return False
                mapper = actor.GetMapper()
                dataset = mapper.GetInput() if mapper is not None else None
                if dataset is None:
                    return False
                if self.clip_registry.is_clipped(actor):
                    dataset = self.clip_registry.original_dataset(actor)
                bounds = dataset.GetBounds() if dataset is not None else None
                if (dataset is not None and dataset.GetNumberOfCells() > 0 and bounds is not None
                        and all(math.isfinite(value) for value in bounds)):
                    valid_surface = True
            except (AttributeError, RuntimeError, TypeError):
                LOGGER.debug("Oberflaechengeometrie fuer %s nicht verfuegbar", self.diagram_id, exc_info=True)
                return False
        return valid_surface

    def _install_axes_observer(self):
        """Tick-/Titelkanten folgen der Kamera: nach jeder Interaktion neu zeichnen (nur Beschriftung, keine Kamera)."""
        if self._axes_observer_installed:
            return

        def on_end(_obj, _event):
            try:
                if self.daten is not None and not self._disposed:
                    self.update_axes(render=True)
            except Exception:
                LOGGER.debug("Achsenaktualisierung nach Interaktion fehlgeschlagen", exc_info=True)

        self._axes_observer_id = self.plotter.iren.add_observer("EndInteractionEvent", on_end)
        self._axes_observer_installed = True

    def _install_hover(self):
        if self._hover_installed:
            return

        def on_move(_obj, _event):
            now = time.monotonic()
            if now - self._last_hover_at < 1 / 30:
                return
            self._last_hover_at = now
            try:
                x, y = self.plotter.iren.get_event_position()
                text = pick_hover_text(self.plotter, x, y)
            except Exception:
                return
            if text != self._last_hover_text:
                self._last_hover_text = text
                self.status_message.emit(text or "")

        self._hover_observer_id = self.plotter.iren.add_observer("MouseMoveEvent", on_move)
        self._hover_installed = True

    # ------------------------------------------------------------------ Layer
    def layer_available(self, key):
        if key not in self.spec.layers:
            return False
        objects = self.actors.get(key)
        if isinstance(objects, (list, tuple)):
            return any(actor is not None for actor in objects)
        return objects is not None

    def available_layers(self):
        return {key: self.layer_available(key) for key in LAYER_LABELS}

    def effective_layers(self):
        return {key: bool(self.layer_state.get(key, False)) and self.layer_available(key) for key in LAYER_LABELS}

    def set_layer_state(self, state):
        for key, value in state.items():
            if key in self.layer_state:
                self.layer_state[key] = bool(value)

    def _set_visibility(self, key, visible):
        objects = self.actors.get(key)
        if objects is None:
            return
        for actor in (objects if isinstance(objects, (list, tuple)) else [objects]):
            if actor is not None:
                actor.visibility = bool(visible)

    def apply_layers(self, render=True):
        if not self.actors:
            return
        effective = self.effective_layers()
        self._set_visibility("surface", effective["surface"])
        self._set_visibility("zweiphasengebiet", effective["zweiphasengebiet"])
        if self.diagram_id == "pvt":
            # Die randkonforme Einphasenflaeche bleibt unabhaengig sichtbar (nur 'surface' steuert sie).
            self._set_visibility("surface_single_phase", False)
        for key, visible in effective.items():
            if key not in ("surface", "zweiphasengebiet"):
                self._set_visibility(key, visible)
        # Zentrale Achsenaktualisierung: Layerwechsel aendert weder Kamera noch manuelle Grenzen.
        self.update_axes(render=False)
        if render:
            self.plotter.render()

    # ------------------------------------------------------------------ Ansichtsbereich (TS-VP04..06)
    def update_axes(self, render=False):
        """Achsenrahmen, physikalische Ticks und Titel; manuelle Box hat Vorrang vor automatischen Grenzen."""
        if not self.actors or self.daten is None:
            return
        axes_view.apply_axes(self.plotter, self.diagram_id, self.daten, self.actors, box=self._manual_box)
        if render:
            self.plotter.render()

    def _new_point_values(self, xyz):
        """Punktdaten neuer Schnittpunkte aus den Darstellungskoordinaten (Rueckwaertstransformation, keine HEOS)."""
        if self.diagram_id != "pvt":
            return None
        T, v, p = diagram_coords.from_xyz(xyz, self.daten["critical_point"], pressure_axis_of(self.daten)).T
        return {"T": T, "v": v, "p_raw": p}

    def _apply_view(self, fit_camera=False, render=True):
        """Clippt die Originalgeometrie an die manuelle Box (oder stellt sie wieder her) und aktualisiert Achsen.

        Rein darstellungsbezogen: Berechnungsarrays, Ergebnis-Dictionary und Cache bleiben unberuehrt.
        """
        if not self.actors or self.daten is None:
            return None
        try:
            box = self.view_state.box(self.daten)
        except ViewportBoundsError as exc:
            LOGGER.warning("Gespeicherte Grenzen nicht anwendbar: %s", exc)
            self.status_message.emit(f"Ansichtsbereich nicht anwendbar, automatische Ansicht: {exc}")
            self.view_state.reset()
            box = None
        stats = self.clip_registry.apply(box, self._new_point_values)
        self._manual_box = box
        if fit_camera:
            target = box if box is not None else axes_view.auto_bounds(self.diagram_id, self.actors)
            if target is not None:
                axes_view.fit_camera(self.plotter, target)
        self.update_axes(render=False)
        if box is not None and stats["empty"] == stats["actors"]:
            self.status_message.emit("Ansichtsbereich enthält keine Geometrie (Rahmen und Achsen bleiben sichtbar).")
        if self._markers or self._projected_markers:
            self._draw_markers()
        elif render:
            self.plotter.render()
        return stats

    def apply_manual_view(self, texts):
        """Prueft die sechs Grenzen (Anzeigeeinheiten) und wendet sie an; bei Fehler bleibt alles unveraendert.

        Wirft ViewportBoundsError. Die Kamera wird auf die neue Box eingepasst (Blickrichtung bleibt erhalten).
        """
        if self.daten is None:
            raise ViewportBoundsError("Noch kein Ergebnis in diesem Diagramm: zuerst die Berechnung starten.")
        candidate = self.view_state.copy()
        candidate.apply_fields(texts, self.daten)
        previous, self.view_state = self.view_state, candidate
        try:
            stats = self._apply_view(fit_camera=True, render=True)
        except Exception as exc:
            LOGGER.exception("Ansichtsbereich konnte nicht dargestellt werden")
            self.view_state = previous                       # letzter gueltiger Ausschnitt bleibt erhalten
            self._apply_view(fit_camera=False, render=True)
            raise ViewportBoundsError(f"Darstellung des Ausschnitts fehlgeschlagen: {exc}") from exc
        self.view_state.camera_state = self.camera_state()
        return stats

    def reset_manual_view(self):
        """Manuelle Begrenzung aus: vollstaendige Originalgeometrie, automatische Grenzen, Kamera auf alles."""
        self.view_state.reset()
        return self._apply_view(fit_camera=True, render=True)

    def set_view_state(self, state, fit_camera=False):
        """Uebernimmt einen gespeicherten Zustand (z. B. nach Neuerzeugung des Viewports)."""
        self.view_state = state
        if self.daten is not None:
            self._apply_view(fit_camera=fit_camera, render=True)

    def restore_default_camera(self):
        """Stellt die Standardkamera des Renderers wieder her (unabhaengig vom Darstellungsfenster)."""
        if self.default_camera is not None:
            self.restore_camera(self.default_camera)
            self.update_axes(render=True)

    def manual_box(self):
        return self._manual_box

    def point_in_manual_view(self, xyz):
        return self._manual_box is None or bool(geometry_clip.points_in_box([xyz], self._manual_box)[0])

    def point_in_visible_surface_bounds(self, xyz):
        bounds = axes_view.visible_bounds({key: self.actors.get(key) for key in ("surface", "zweiphasengebiet")})
        return bounds is not None and bool(geometry_clip.points_in_box([xyz], bounds)[0])

    def clip_stats(self):
        return {"enabled": self._manual_box is not None, "box": self._manual_box,
                "actors": len(self.clip_registry.entries)}

    def apply_grid(self, render=True):
        if not self.actors:
            return
        for key in GRID_ACTOR_KEYS:
            objects = self.actors.get(key)
            if objects is None:
                continue
            for actor in (objects if isinstance(objects, (list, tuple)) else [objects]):
                if actor is not None:
                    actor.prop.show_edges = bool(self.grid_on)
        if render:
            self.plotter.render()

    def set_grid(self, on):
        self.grid_on = bool(on)
        self.apply_grid()

    # ------------------------------------------------------------------ Punktauswahl (Zustandsvergleich)
    @property
    def picking_enabled(self):
        return self._filter.click_enabled

    @picking_enabled.setter
    def picking_enabled(self, enabled):
        self._filter.click_enabled = bool(enabled)

    def _on_qt_click(self, x, y):
        ratio = self.plotter.devicePixelRatioF()
        self.point_clicked.emit(self, int(round(x * ratio)), int(round((self.plotter.height() - 1 - y) * ratio)))

    def set_selection_markers(self, markers, projected_markers=None):
        """Aktualisiert Quellen- und diagrammübergreifende Auswahlmarker ohne Szenen-Neuaufbau."""
        self._markers = [(int(i), tuple(float(c) for c in xyz), color) for i, xyz, color in markers if xyz is not None]
        if projected_markers is not None:
            self._projected_markers = [(int(i), tuple(float(c) for c in xyz), color)
                                       for i, xyz, color in projected_markers if xyz is not None]
        self._draw_markers()

    def _draw_markers(self):
        import pyvista as pv
        for actor in self._marker_actors:
            try:
                self.plotter.remove_actor(actor, render=False)
            except Exception:
                LOGGER.debug("Markierung konnte nicht entfernt werden", exc_info=True)
        self._marker_actors = []
        markers_by_index = {index: (index, xyz, color) for index, xyz, color in self._projected_markers}
        markers_by_index.update({index: (index, xyz, color) for index, xyz, color in self._markers})
        for index, xyz, color in markers_by_index.values():
            # Marker ausserhalb des manuellen Ansichtsbereichs werden ausgeblendet, aber nicht aus den Daten geloescht.
            if not self.point_in_manual_view(xyz):
                continue
            point = pv.PolyData([xyz])
            self._marker_actors.append(self.plotter.add_mesh(
                point, color=color, point_size=16, render_points_as_spheres=True, pickable=False, reset_camera=False))
            self._marker_actors.append(self.plotter.add_point_labels(
                [xyz], [str(index)], font_size=16, bold=True, text_color="black", show_points=False, shape=None,
                always_visible=True, pickable=False, reset_camera=False, render=False))
        self.plotter.render()

    @property
    def marker_count(self):
        return len(self._markers)

    @property
    def projected_marker_count(self):
        return len(self._projected_markers)

    def set_fvt_tangent(self, tangent, show_slice=False):
        """Aktualisiert die nicht-pickbare gemeinsame Tangente, ohne die Szene oder Kamera zurueckzusetzen."""
        self._tangent_data = tangent
        self._tangent_show_slice = bool(show_slice)
        if self.diagram_id != "fvt" or not self.has_result or self._disposed:
            return
        import pyvista as pv

        for actor in self._tangent_actors:
            try:
                self.plotter.remove_actor(actor, render=False)
            except (RuntimeError, ValueError):
                LOGGER.debug("Tangenten-Overlay konnte nicht entfernt werden", exc_info=True)
        self._tangent_actors = []
        if tangent is not None:
            from fvt_tangent import tangent_display_points

            xyz = tangent_display_points(tangent, self.daten)
            line = pv.lines_from_points(xyz)
            self._tangent_actors.append(self.plotter.add_mesh(
                line, color="lime", line_width=5, pickable=False, reset_camera=False, render=False))
            endpoints = pv.PolyData(xyz[[0, -1]])
            self._tangent_actors.append(self.plotter.add_mesh(
                endpoints, color="lime", point_size=18, render_points_as_spheres=True, pickable=False,
                reset_camera=False, render=False))
            if self._tangent_show_slice:
                self._add_tangent_slice_plane(pv, tangent.temperature)
        self.plotter.render()

    def _add_tangent_slice_plane(self, pv, temperature):
        bounds = self.plotter.bounds
        if bounds is None or len(bounds) != 6 or not all(math.isfinite(value) for value in bounds):
            return
        y = float(temperature) / float(self.daten["critical_point"][0])
        points = np.array(((bounds[0], y, bounds[4]), (bounds[1], y, bounds[4]),
                           (bounds[1], y, bounds[5]), (bounds[0], y, bounds[5])))
        plane = pv.PolyData(points, np.array((4, 0, 1, 2, 3)))
        self._tangent_actors.append(self.plotter.add_mesh(
            plane, color="lime", opacity=0.12, show_edges=False, pickable=False, reset_camera=False, render=False))

    # ------------------------------------------------------------------ Kamera
    def camera_state(self):
        return [tuple(float(x) for x in point) for point in self.plotter.camera_position]

    def restore_camera(self, state):
        self.plotter.camera_position = state
        self.plotter.render()

    # ------------------------------------------------------------------ Lebenszyklus
    def dispose(self):
        """Gibt Plotter, Beobachter und Qt-Widget frei (idempotent)."""
        if self._disposed:
            return
        LOGGER.debug("event=viewport_dispose diagram=%s viewport=%x actors=%d",
                     self.diagram_id, id(self), self._scene_actor_count())
        self._disposed = True
        self._render_state = DISPOSED
        try:
            self.plotter.removeEventFilter(self._filter)
            for observer_id in (self._hover_observer_id, self._axes_observer_id):
                if observer_id is not None:
                    self.plotter.iren.remove_observer(observer_id)
            self.plotter.clear()
            self.plotter.close()
        except Exception:
            LOGGER.exception("Viewport %s konnte nicht sauber geschlossen werden", self.diagram_id)
        self.actors = {}
        self.daten = None
        self.hide()
        self.setParent(None)
        self.deleteLater()
