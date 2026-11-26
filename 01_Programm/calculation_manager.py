"""
Datei: calculation_manager.py
Projekt: ThermoSurface
Kategorie: Berechnungssteuerung

Beschreibung:
Koordiniert die sequenzielle Hintergrundberechnung ausgewählter Diagramme.

Aufgaben:
- Verwaltet DiagramTask, Fortschritt und Fehler über QThread-Signale.
- Nutzt den Ergebnis-Cache und berechnet nur fehlende Datensätze.

Abhängigkeiten:
- Projektmodule (direkt): diagram_cache.py, diagram_registry.py
- Bibliotheken: PySide6

Verwendet von:
- Projektmodule: gui.py
"""
import logging
from dataclasses import dataclass, field
from typing import Any, Optional

from PySide6.QtCore import QObject, QThread, QTimer, Signal, Slot

from diagram_cache import DiagramKey
from diagram_registry import DIAGRAMS

LOGGER = logging.getLogger("thermosurface.calculation")


@dataclass
class DiagramTask:
    diagram_id: str
    stoff: str
    raster: dict
    visu_raster: dict
    label: str = ""
    isolines: Optional[dict] = None   # TS-VP03: start/stop/factor je Isolinienfamilie; None = config.yaml

    def __post_init__(self):
        if not self.label:
            self.label = DIAGRAMS.get(self.diagram_id).label if self.diagram_id in DIAGRAMS else self.diagram_id


@dataclass
class TaskOutcome:
    task: DiagramTask
    daten: Any = None
    error: Optional[str] = None
    from_cache: bool = False
    key: Optional[DiagramKey] = None

    @property
    def ok(self):
        return self.error is None


@dataclass
class BatchResult:
    outcomes: list = field(default_factory=list)
    calculations: int = 0
    cache_hits: int = 0

    @property
    def failed(self):
        return [o for o in self.outcomes if not o.ok]

    @property
    def succeeded(self):
        return [o for o in self.outcomes if o.ok]


def cache_key_for(task):
    spec = DIAGRAMS.get(task.diagram_id)
    if spec.cache_params is not None:
        extra = {"isolines": task.isolines} if task.isolines is not None else {}
        params = spec.cache_params(task.raster, task.visu_raster, **extra)
    else:
        params = {"raster": task.raster, "visu_raster": task.visu_raster}
        if task.isolines is not None:
            params["isolines"] = task.isolines
    return DiagramKey.build(task.stoff, task.diagram_id, params)


class CalculationWorker(QObject):
    """Arbeitet eine Aufgabenliste nacheinander ab; ein Fehler stoppt die Folgeaufgaben nicht."""

    progress = Signal(int, str)
    task_done = Signal(int, object)
    task_failed = Signal(int, str)
    finished_all = Signal()

    def __init__(self, berechnen, tasks):
        super().__init__()
        self.berechnen = berechnen
        self.tasks = list(tasks)
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    @Slot()
    def run(self):
        total = len(self.tasks)
        try:
            for index, task in enumerate(self.tasks):
                if self._cancelled:
                    break

                def report(percent, phase, index=index, task=task):
                    fraction = (index + min(max(float(percent), 0.0), 95.0) / 95.0) / total
                    self.progress.emit(int(fraction * 95), f"[{index + 1}/{total}] {task.label}: {phase}")

                try:
                    extra = {"isolines": task.isolines} if task.isolines is not None else {}
                    daten = self.berechnen(task.stoff, task.raster, task.visu_raster,
                                           diagram_type=task.diagram_id, progress_callback=report, **extra)
                except Exception as exc:
                    LOGGER.exception("Berechnung %s fehlgeschlagen", task.diagram_id)
                    self.task_failed.emit(index, str(exc) or type(exc).__name__)
                else:
                    self.task_done.emit(index, daten)
        finally:
            self.finished_all.emit()


class CalculationManager(QObject):
    progress = Signal(int, str)
    finished = Signal(object)  # BatchResult

    def __init__(self, berechnen, cache, parent=None):
        super().__init__(parent)
        self.berechnen = berechnen
        self.cache = cache
        self.compute_calls = 0  # tatsaechlich an `berechnen` uebergebene Aufgaben (Cache-Misses)
        self._thread = None
        self._worker = None
        self._pending = []
        self._batch = None

    @property
    def is_running(self):
        return self._batch is not None

    def start(self, tasks):
        """Startet die Warteschlange; False, wenn bereits eine laeuft oder keine Aufgaben vorliegen."""
        tasks = list(tasks)
        if self.is_running or not tasks:
            return False

        batch = BatchResult()
        outcomes = []
        misses = []
        for task in tasks:
            outcome = TaskOutcome(task)
            outcomes.append(outcome)
            try:
                outcome.key = cache_key_for(task)
            except Exception as exc:
                outcome.error = str(exc)
                continue
            cached = self.cache.get(outcome.key)
            if cached is not None:
                outcome.daten = cached
                outcome.from_cache = True
                batch.cache_hits += 1
            else:
                misses.append(outcome)
        batch.outcomes = outcomes
        self._batch = batch
        self._pending = misses

        if not misses:
            self.progress.emit(95, "Alle Diagramme aus dem Cache")
            QTimer.singleShot(0, self._complete)
            return True

        self._worker = CalculationWorker(self.berechnen, [o.task for o in misses])
        self._thread = QThread(self)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.progress.connect(self.progress)
        self._worker.task_done.connect(self._on_task_done)
        self._worker.task_failed.connect(self._on_task_failed)
        self._worker.finished_all.connect(self._on_worker_done)
        self._thread.start()
        return True

    @Slot(int, object)
    def _on_task_done(self, index, daten):
        outcome = self._pending[index]
        outcome.daten = daten
        self.compute_calls += 1
        self._batch.calculations += 1
        evicted = self.cache.put(outcome.key, daten)
        if evicted:
            LOGGER.info("%d aeltere Cache-Eintraege verdraengt", len(evicted))

    @Slot(int, str)
    def _on_task_failed(self, index, message):
        self._pending[index].error = message

    @Slot()
    def _on_worker_done(self):
        thread, worker = self._thread, self._worker
        if thread is not None:
            thread.quit()
            thread.wait()
            thread.deleteLater()
        if worker is not None:
            worker.deleteLater()
        self._thread = None
        self._worker = None
        self._complete()

    def _complete(self):
        batch, self._batch, self._pending = self._batch, None, []
        if batch is not None:
            self.finished.emit(batch)

    def shutdown(self):
        """Wartet auf eine laufende Berechnung (bricht zwischen Aufgaben ab) und beendet den QThread."""
        if self._worker is not None:
            self._worker.cancel()
        if self._thread is not None:
            self._thread.quit()
            self._thread.wait()
