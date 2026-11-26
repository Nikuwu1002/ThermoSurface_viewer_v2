"""
Datei: fvt_slice_view.py
Projekt: ThermoSurface
Kategorie: GUI / f-v-Isotherme

Beschreibung:
Zeigt die gemeinsame f-v-T-Tangente in einer ergänzenden 2D-Matplotlib-Ansicht.

Aufgaben:
- Zeichnet Tangente, Sättigungspunkte und homogene Isothermenabschnitte.
- Berechnet zusätzliche Kurvendaten asynchron über einen Einzel-Worker.

Abhängigkeiten:
- Projektmodule (direkt): fvt_tangent.py
- Bibliotheken: PySide6, matplotlib, numpy

Verwendet von:
- Projektmodule: gui.py
"""

from concurrent.futures import ThreadPoolExecutor
import logging

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import QObject, Signal, Slot
from PySide6.QtWidgets import QVBoxLayout, QWidget

from fvt_tangent import tangent_plot_data

LOGGER = logging.getLogger("thermosurface.fvt_slice")


class FvtIsothermSampler(QObject):
    completed = Signal(object, object, object)
    _finished = Signal(int, object, object, object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="fvt-isotherm")
        self._generation = 0
        self._running = False
        self._pending = None
        self._closed = False
        self._finished.connect(self._on_finished)

    def invalidate(self):
        self._generation += 1
        self._pending = None

    def request(self, tangent, data):
        self.invalidate()
        self._pending = (self._generation, tangent, {"meta": dict(data["meta"])})
        if not self._running:
            self._start_pending()

    def _start_pending(self):
        generation, tangent, data = self._pending
        self._pending = None
        self._running = True

        def sample():
            from fvt_tangent import homogeneous_fvt_isotherm
            return homogeneous_fvt_isotherm(tangent, data)

        def finished(future):
            try:
                branches = future.result()
            except Exception as exc:
                LOGGER.exception("Homogene f-v-Isotherme konnte nicht berechnet werden")
                self._finished.emit(generation, tangent, None, str(exc) or type(exc).__name__)
            else:
                self._finished.emit(generation, tangent, branches, None)

        self._executor.submit(sample).add_done_callback(finished)

    @Slot(int, object, object, object)
    def _on_finished(self, generation, tangent, branches, error):
        self._running = False
        if not self._closed and generation == self._generation:
            self.completed.emit(tangent, branches, error)
        if not self._closed and self._pending is not None:
            self._start_pending()

    def shutdown(self):
        self._closed = True
        self.invalidate()
        self._executor.shutdown(wait=True, cancel_futures=True)


class FvtSliceView(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure()
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

    def show_error(self, message):
        self.figure.clear()
        axes = self.figure.add_subplot(111)
        axes.text(0.5, 0.5, message, ha="center", va="center", wrap=True, transform=axes.transAxes)
        axes.set_axis_off()
        self.toolbar.update()
        self.canvas.draw_idle()

    def update_tangent(self, tangent, homogeneous=None):
        self.figure.clear()
        axes = self.figure.add_subplot(111)
        data = tangent_plot_data(tangent)
        axes.plot(data["volume"], data["helmholtz"], color="green", linewidth=2.2, label="Gemeinsame Tangente")
        axes.plot((0.0, tangent.v_l),
                  (tangent.mu, tangent.mu - tangent.p_sat * tangent.v_l),
                  color="green", linewidth=1.2, linestyle="--")
        axes.scatter((tangent.v_l, tangent.v_v), (tangent.f_l_molar, tangent.f_v_molar),
                     color=("blue", "red"), zorder=3)
        axes.scatter((data["intercept_volume"],), (data["intercept"],), color="black", marker="x", zorder=4)
        axes.annotate("μ = gₘ", (0.0, data["intercept"]), xytext=(6, 8), textcoords="offset points")
        if homogeneous is not None:
            for index, (volumes, energies) in enumerate(homogeneous):
                axes.plot(np.asarray(volumes), np.asarray(energies), color="gray", alpha=0.8,
                          label="Homogene HEOS-Fortsetzung" if index == 0 else None)
        axes.set_xlabel("Molares Volumen vₘ [m³/mol]")
        axes.set_ylabel("Molare Helmholtz-Energie fₘ [J/mol]")
        axes.set_title(f"f-v-Isotherme bei T = {tangent.temperature:.2f} K")
        axes.grid(True, alpha=0.25)
        axes.legend(loc="best")
        axes.set_xlim(left=0.0)
        self.figure.tight_layout()
        self.toolbar.update()
        self.canvas.draw_idle()
