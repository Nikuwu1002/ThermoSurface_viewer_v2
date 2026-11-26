"""
Datei: comparison_view.py
Projekt: ThermoSurface
Kategorie: GUI / Zustandsvergleich

Beschreibung:
Zeigt ausgewählte thermodynamische Zustände in einem separaten 3D-Vergleich.

Aufgaben:
- Erzeugt die Punkt- und Verbindungsgeometrie für drei wählbare Achsen.
- Stellt Punkte, Legende und Warnungen mit Qt und PyVista dar.

Abhängigkeiten:
- Projektmodule (direkt): keine
- Weitere Imports: thermo_state.py
- Bibliotheken: PySide6, numpy, pyvista, pyvistaqt

Verwendet von:
- Projektmodule: gui.py
"""

import logging
from dataclasses import dataclass, field

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QComboBox, QHBoxLayout, QLabel, QScrollArea, QSplitter, QVBoxLayout, QWidget)

from thermo_state import QUANTITIES

LOGGER = logging.getLogger("thermosurface.comparison")

AXIS_ORDER = ("T", "p", "v", "h", "s", "u", "f", "g")
DEFAULT_AXES = ("p", "v", "T")          # sinnvolle Startbelegung: x = p, y = v, z = T

# Anzeigeeinheiten (Werte bleiben SI in den Datensaetzen); Faktor: Anzeige = SI * Faktor
_DISPLAY_FACTOR = {"h": (1e-3, "kJ/mol"), "u": (1e-3, "kJ/mol"), "f": (1e-3, "kJ/mol"), "g": (1e-3, "kJ/mol")}


@dataclass
class AxisScale:
    quantity: str
    mode: str          # linear | log10
    title: str
    factor: float = 1.0

    def apply(self, value):
        if value is None or not np.isfinite(value):
            return None
        scaled = float(value) * self.factor
        if self.mode == "log10":
            return float(np.log10(scaled)) if scaled > 0 else None
        return scaled


def axis_scale(quantity, values):
    """Skala je Achse: p und v logarithmisch, solange alle Werte > 0 sind; sonst linear (mit Hinweis im Titel)."""
    label, unit = QUANTITIES[quantity]
    finite = [float(x) for x in values if x is not None and np.isfinite(x)]
    if quantity in ("p", "v"):
        if finite and min(finite) > 0:
            return AxisScale(quantity, "log10", f"log10({label} / {unit})")
        if quantity == "p":
            return AxisScale(quantity, "linear", f"{label} [MPa] (linear, da p ≤ 0 vorkommt)", 1e-6)
        return AxisScale(quantity, "linear", f"{label} [{unit}] (linear)")
    factor, shown_unit = _DISPLAY_FACTOR.get(quantity, (1.0, unit))
    return AxisScale(quantity, "linear", f"{label} [{shown_unit}]", factor)


@dataclass
class ComparisonRow:
    selection_index: int
    selection_order: int
    color: str
    source_diagram: str
    xyz: tuple = None              # None: nicht darstellbar
    reason: str = ""
    description: str = ""
    status: str = ""


@dataclass
class ComparisonGeometry:
    axes: tuple
    scales: tuple
    rows: list = field(default_factory=list)
    segments: list = field(default_factory=list)    # (row_position_a, row_position_b), nur unmittelbar aufeinanderfolgende gueltige Punkte
    warnings: list = field(default_factory=list)

    @property
    def valid_rows(self):
        return [r for r in self.rows if r.xyz is not None]


def build_geometry(points, axes):
    """Reine Geometrie aus ausgewaehlten Punkten (in Auswahlreihenfolge) und drei verschiedenen Achsen."""
    axes = tuple(axes)
    if len(axes) != 3 or len(set(axes)) != 3 or any(a not in QUANTITIES for a in axes):
        raise ValueError("Drei verschiedene Zustandsgroessen aus T, p, v, h, s, u, f, g erforderlich")
    points = list(points)
    scales = tuple(axis_scale(a, [p.state.value(a) for p in points]) for a in axes)
    geometry = ComparisonGeometry(axes, scales)
    for point in points:
        state = point.state
        coords, missing = [], []
        for axis, scale in zip(axes, scales):
            value = state.value(axis)
            shown = scale.apply(value)
            if shown is None:
                missing.append(axis if value is None else f"{axis} (nicht darstellbar auf {scale.mode}-Achse)")
            coords.append(shown)
        row = ComparisonRow(point.selection_index, point.selection_order, point.color, point.source_diagram,
                            description=point.description(), status=state.validity)
        if missing:
            row.reason = ", ".join(missing) + f" nicht berechnet ({state.validity})"
            geometry.warnings.append(f"Punkt {point.selection_index}: {row.reason} – nicht dargestellt, Verbindung unterbrochen")
        else:
            row.xyz = tuple(coords)
        geometry.rows.append(row)
    for k in range(len(geometry.rows) - 1):
        a, b = geometry.rows[k], geometry.rows[k + 1]
        if a.xyz is not None and b.xyz is not None:
            geometry.segments.append((k, k + 1))
    return geometry


class ComparisonView(QWidget):
    """Eigenstaendiges Fenster: drei Achsen-Comboboxen, 3D-Ansicht, Legende, Warnungen."""

    def __init__(self, selection, parent=None):
        super().__init__(parent, Qt.WindowType.Window)
        from pyvistaqt import QtInteractor
        self.selection = selection
        self.setWindowTitle("ThermoSurface – Zustandsvergleich (3D)")
        self.resize(1000, 700)
        self._drawn_axes = None
        self.refresh_count = 0
        self.geometry = None
        self._disposed = False

        root = QVBoxLayout(self)
        bar = QHBoxLayout()
        self.axis_combos = []
        for name, default in zip(("x", "y", "z"), DEFAULT_AXES):
            bar.addWidget(QLabel(f"{name}-Achse:"))
            combo = QComboBox()
            for quantity in AXIS_ORDER:
                label, unit = QUANTITIES[quantity]
                combo.addItem(f"{label} [{unit}]", quantity)
            combo.setCurrentIndex(combo.findData(default))
            combo.currentIndexChanged.connect(lambda _i, c=combo: self._on_axis_changed(c))
            self.axis_combos.append(combo)
            bar.addWidget(combo)
        bar.addStretch(1)
        root.addLayout(bar)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        self.plotter = QtInteractor(self)
        splitter.addWidget(self.plotter)
        self.legend = QLabel()
        self.legend.setTextFormat(Qt.TextFormat.RichText)
        self.legend.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
        self.legend.setWordWrap(True)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.legend)
        splitter.addWidget(scroll)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)
        root.addWidget(splitter, 1)
        self.warning_label = QLabel()
        self.warning_label.setWordWrap(True)
        self.warning_label.setStyleSheet("color: #b36b00;")
        root.addWidget(self.warning_label)

        self._previous_axes = self.axes
        self.selection.subscribe(self._on_selection_event)

    # ------------------------------------------------------------------ Achsen
    @property
    def axes(self):
        return tuple(c.currentData() for c in self.axis_combos)

    def set_axes(self, x, y, z):
        if len({x, y, z}) != 3:
            raise ValueError("Achsen muessen verschieden sein")
        for combo, value in zip(self.axis_combos, (x, y, z)):
            combo.blockSignals(True)
            combo.setCurrentIndex(combo.findData(value))
            combo.blockSignals(False)
        self._previous_axes = self.axes
        self.refresh()

    def _on_axis_changed(self, changed):
        # Dieselbe Variable nie doppelt: bei Kollision tauscht die andere Combobox auf den bisherigen Wert.
        position = self.axis_combos.index(changed)
        new_value = changed.currentData()
        old_value = self._previous_axes[position]
        for other_position, other in enumerate(self.axis_combos):
            if other is not changed and other.currentData() == new_value:
                other.blockSignals(True)
                other.setCurrentIndex(other.findData(old_value))
                other.blockSignals(False)
        self._previous_axes = self.axes
        self.refresh()

    # ------------------------------------------------------------------ Darstellung
    def _on_selection_event(self, _event):
        if self.isVisible():
            self.refresh()

    def showEvent(self, event):
        super().showEvent(event)
        self.refresh()

    def refresh(self):
        """Baut die Geometrie neu auf (nur Anzeige; keine Thermodynamikberechnung)."""
        if self._disposed:
            return
        import pyvista as pv
        self.refresh_count += 1
        geometry = build_geometry(self.selection.points, self.axes)
        self.geometry = geometry
        self.plotter.clear()
        valid = geometry.valid_rows
        if valid:
            xyz = np.array([r.xyz for r in valid], dtype=float)
            position = {id(r): k for k, r in enumerate(valid)}
            for color in dict.fromkeys(r.color for r in valid):
                members = [position[id(r)] for r in valid if r.color == color]
                self.plotter.add_mesh(pv.PolyData(xyz[members]), color=color, point_size=18,
                                      render_points_as_spheres=True, reset_camera=False)
            if geometry.segments:
                lines = []
                for a, b in geometry.segments:
                    lines += [2, position[id(geometry.rows[a])], position[id(geometry.rows[b])]]
                poly = pv.PolyData(xyz)
                poly.lines = np.array(lines, dtype=np.int64)
                self.plotter.add_mesh(poly, color="dimgray", line_width=3, reset_camera=False)
            self.plotter.add_point_labels(xyz, [str(r.selection_index) for r in valid], font_size=16, bold=True,
                                          text_color="black", show_points=False, shape=None, always_visible=True,
                                          reset_camera=False, render=False)
            lo, hi = xyz.min(axis=0), xyz.max(axis=0)
            pad = np.maximum((hi - lo) * 0.1, np.maximum(np.abs(hi + lo) * 0.005, 1e-3))
            bounds = (lo[0] - pad[0], hi[0] + pad[0], lo[1] - pad[1], hi[1] + pad[1], lo[2] - pad[2], hi[2] + pad[2])
            self.plotter.show_bounds(bounds=bounds, xtitle=geometry.scales[0].title, ytitle=geometry.scales[1].title,
                                     ztitle=geometry.scales[2].title)
            if self._drawn_axes != geometry.axes:
                self.plotter.reset_camera()
                self._drawn_axes = geometry.axes
        else:
            self.plotter.remove_bounds_axes()
        self.plotter.render()
        self.legend.setText(self._legend_html(geometry))
        self.warning_label.setText("\n".join(geometry.warnings))

    @staticmethod
    def _legend_html(geometry):
        if not geometry.rows:
            return "<i>Keine Punkte ausgewählt.</i>"
        lines = ["<b>Punkte (Auswahlreihenfolge)</b><br>"]
        for row in sorted(geometry.rows, key=lambda r: r.selection_order):
            marker = f"<span style='color:{row.color}'>&#9679;</span>"
            shown = "" if row.xyz is not None else f" <i>(nicht dargestellt: {row.reason})</i>"
            lines.append(f"{marker} <b>{row.selection_index}</b> – {row.description} [{row.status}]{shown}<br>")
        diagrams = dict.fromkeys((r.source_diagram, r.color) for r in geometry.rows)
        lines.append("<br><b>Herkunft (Farbe)</b><br>")
        for diagram, color in diagrams:
            lines.append(f"<span style='color:{color}'>&#9679;</span> {diagram}<br>")
        return "".join(lines)

    def dispose(self):
        if self._disposed:
            return
        self._disposed = True
        self.selection.unsubscribe(self._on_selection_event)
        try:
            self.plotter.close()
        except Exception:
            LOGGER.debug("Vergleichsplotter konnte nicht geschlossen werden", exc_info=True)
        self.deleteLater()

    def closeEvent(self, event):
        LOGGER.debug("event=compare_close selection_count=%d view=%x", len(self.selection), id(self))
        event.ignore()
        self.hide()
