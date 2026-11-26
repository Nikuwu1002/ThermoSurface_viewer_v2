"""
Datei: diagram_workspace.py
Projekt: ThermoSurface
Kategorie: GUI / Mehrfachansicht

Beschreibung:
Ordnet mehrere DiagramViewport-Instanzen in einer Qt-Kachelansicht an.

Aufgaben:
- Erstellt, aktiviert, maximiert und entfernt Diagramm-Viewports.
- Erhält Ansichtseinstellungen bei gezieltem Austausch eines Viewports.

Abhängigkeiten:
- Projektmodule (direkt): diagram_viewport.py
- Weitere Imports: tile_layout.py
- Bibliotheken: PySide6

Verwendet von:
- Projektmodule: gui.py
"""

import logging
import math

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QGridLayout, QHBoxLayout, QPushButton, QVBoxLayout, QWidget

from diagram_viewport import DiagramViewport
from tile_layout import tile_grid_shape

LOGGER = logging.getLogger("thermosurface.workspace")


class DiagramWorkspace(QWidget):
    active_changed = Signal(object)      # DiagramViewport | None
    maximized_changed = Signal(object)   # DiagramViewport | None
    status_message = Signal(str)
    point_clicked = Signal(object, int, int)   # (viewport, x, y): Einfachklick fuer die Punktauswahl

    def __init__(self, parent=None):
        super().__init__(parent)
        self.viewports = {}
        self._active = None
        self._maximized = None
        self._single_mode = True
        self._grid_extent = (0, 0)
        self._picking = False

        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(2)
        bar = QHBoxLayout()
        bar.setContentsMargins(0, 0, 0, 0)
        self.return_button = QPushButton("◀ Zurück zur Mehrfachansicht")
        self.return_button.setVisible(False)
        self.return_button.clicked.connect(self.restore)
        bar.addWidget(self.return_button)
        bar.addStretch(1)
        root.addLayout(bar)

        self._grid_host = QWidget()
        self._grid = QGridLayout(self._grid_host)
        self._grid.setContentsMargins(0, 0, 0, 0)
        self._grid.setSpacing(4)
        root.addWidget(self._grid_host, 1)

    # ------------------------------------------------------------------ Zugriff
    @property
    def active(self):
        return self._active

    @property
    def maximized(self):
        return self._maximized

    def ordered(self):
        return list(self.viewports.values())

    def set_single_mode(self, single):
        self._single_mode = bool(single)
        for viewport in self.viewports.values():
            viewport.set_header_visible(not self._single_mode)

    @property
    def picking_enabled(self):
        return self._picking

    @picking_enabled.setter
    def picking_enabled(self, enabled):
        self._picking = bool(enabled)
        for viewport in self.viewports.values():
            viewport.picking_enabled = self._picking

    # ------------------------------------------------------------------ Aufbau
    def ensure(self, diagram_ids):
        """Haelt genau die Viewports fuer `diagram_ids` (Reihenfolge = Kachelreihenfolge) vor.

        Vorhandene Viewports gleicher ID bleiben erhalten (Kamera, Layer); nicht mehr
        benoetigte werden freigegeben, fehlende neu erzeugt. Gibt (viewports, neu) zurueck.
        """
        diagram_ids = list(dict.fromkeys(diagram_ids))
        previous = self._active
        if self._maximized is not None:
            self.restore()

        for diagram_id in list(self.viewports):
            if diagram_id not in diagram_ids:
                self._remove(diagram_id)

        created = set()
        ordered = {}
        for diagram_id in diagram_ids:
            viewport = self.viewports.get(diagram_id)
            if viewport is None:
                viewport = self._new_viewport(diagram_id)
                created.add(diagram_id)
            ordered[diagram_id] = viewport
        self.viewports = ordered

        if self._active is None:
            self._active = next(iter(self.viewports.values()), None)
            if self._active is not None:
                self._active.set_active(True)
        if self._active is not previous:
            self.active_changed.emit(self._active)
        self._relayout()
        return ordered, created

    def _new_viewport(self, diagram_id):
        viewport = DiagramViewport(diagram_id, self._grid_host, show_header=not self._single_mode)
        viewport.activated.connect(self.set_active)
        viewport.maximize_requested.connect(self.maximize)
        viewport.status_message.connect(self.status_message)
        viewport.point_clicked.connect(self.point_clicked)
        viewport.picking_enabled = self._picking
        return viewport

    def replace_viewport(self, diagram_id):
        """Replace one renderer while preserving its diagram-local presentation preferences."""
        old = self.viewports.get(diagram_id)
        if old is None:
            return None

        was_active = old is self._active
        was_maximized = old is self._maximized
        replacement = self._new_viewport(diagram_id)
        replacement.layer_state = dict(old.layer_state)
        replacement.grid_on = old.grid_on
        replacement.tooltips = dict(old.tooltips)
        replacement.view_state = old.view_state.copy()
        replacement._markers = list(old._markers)
        try:
            camera = old.camera_state()
            if len(camera) != 3 or not all(len(point) == 3 for point in camera):
                raise ValueError("ungueltige Kamerakoordinaten")
            if not all(math.isfinite(value) for point in camera for value in point):
                raise ValueError("nicht endliche Kamerakoordinaten")
            replacement._recovery_camera = camera
        except (AttributeError, RuntimeError, TypeError, ValueError):
            LOGGER.warning("Kamerazustand fuer %s nicht uebernehmbar; Standardkamera wird verwendet.",
                           diagram_id, exc_info=True)
        self._grid.removeWidget(old)
        self.viewports[diagram_id] = replacement
        if was_active:
            self._active = replacement
        if was_maximized:
            self._maximized = replacement
        old.dispose()
        replacement.set_active(was_active)
        self._relayout()
        if was_active:
            self.active_changed.emit(replacement)
        if was_maximized:
            self.maximized_changed.emit(replacement)
        return replacement

    def _remove(self, diagram_id):
        viewport = self.viewports.pop(diagram_id)
        self._grid.removeWidget(viewport)
        if viewport is self._active:
            self._active = None
        viewport.dispose()

    def clear(self):
        """Ansicht zuruecksetzen (Moduswechsel); Berechnungsdaten liegen im Cache, nicht hier."""
        previous = self._active
        if self._maximized is not None:
            self.restore()
        for diagram_id in list(self.viewports):
            self._remove(diagram_id)
        self._active = None
        if previous is not None:
            self.active_changed.emit(None)
        self._relayout()

    def close_all(self):
        self.clear()

    # ------------------------------------------------------------------ Aktivierung
    def set_active(self, viewport):
        self._set_active(viewport)

    def _set_active(self, viewport):
        if viewport is self._active and (viewport is None or viewport in self.viewports.values()):
            if viewport is not None:
                viewport.set_active(True)
            return
        if self._active is not None and self._active in self.viewports.values():
            self._active.set_active(False)
        self._active = viewport
        if viewport is not None:
            viewport.set_active(True)
        self.active_changed.emit(viewport)

    # ------------------------------------------------------------------ Maximieren
    def maximize(self, viewport):
        """Zeigt nur diesen bereits vorhandenen Viewport; es wird nichts berechnet."""
        if viewport is None or viewport not in self.viewports.values():
            return
        if self._single_mode or len(self.viewports) < 2 or self._maximized is not None:
            return
        self._maximized = viewport
        self._set_active(viewport)
        self.return_button.setVisible(True)
        self._relayout()
        self.maximized_changed.emit(viewport)

    def restore(self):
        if self._maximized is None:
            return
        self._maximized = None
        self.return_button.setVisible(False)
        self._relayout()
        self.maximized_changed.emit(None)

    # ------------------------------------------------------------------ Layout
    def _relayout(self):
        for viewport in self.viewports.values():
            self._grid.removeWidget(viewport)

        old_rows, old_cols = self._grid_extent
        for row in range(old_rows):
            self._grid.setRowStretch(row, 0)
        for col in range(old_cols):
            self._grid.setColumnStretch(col, 0)

        if self._maximized is not None:
            for viewport in self.viewports.values():
                viewport.setVisible(viewport is self._maximized)
            self._grid.addWidget(self._maximized, 0, 0)
            self._grid.setRowStretch(0, 1)
            self._grid.setColumnStretch(0, 1)
            self._grid_extent = (1, 1)
            return

        rows, cols = tile_grid_shape(len(self.viewports))
        for index, viewport in enumerate(self.viewports.values()):
            self._grid.addWidget(viewport, index // cols, index % cols)
            viewport.setVisible(True)
        # Gleichmaessige Kacheln; nicht belegte Kacheln bleiben leer (kein Renderer).
        for row in range(rows):
            self._grid.setRowStretch(row, 1)
        for col in range(cols):
            self._grid.setColumnStretch(col, 1)
        self._grid_extent = (rows, cols)
