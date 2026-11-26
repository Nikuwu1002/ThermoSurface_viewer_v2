
"""
Datei: selection_manager.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Verwaltung ausgewählter thermodynamischer Zustandspunkte
für den diagrammübergreifenden Zustandsvergleich.

Aufgaben:
- Speichern ausgewählter Zustandspunkte innerhalb einer Sitzung
- Vergabe eindeutiger Punktnummern und Diagrammfarben
- Hinzufügen, Entfernen und Umsortieren von Punkten
- Prüfung der Stoffidentität und Zustandsgültigkeit
- Benachrichtigung registrierter Listener bei Änderungen

Abhängigkeiten:
- thermo_state.py: Thermodynamische Zustandsdatensätze

Verwendet von:
- gui.py: Verwaltung und Steuerung der Punktauswahl
- comparison_view.py: Funktionale Nutzung der Auswahlpunkte
"""

from dataclasses import dataclass, field
from typing import Callable, Optional

from thermo_state import ThermoState

# Farbe je Ursprungsdiagramm (Herkunft); unbekannte Diagramme erhalten reihum Zusatzfarben.
DIAGRAM_COLORS = {
    "pvt": "#1f77b4",
    "pht": "#d62728",
    "pst": "#2ca02c",
    "hsp": "#ff7f0e",
    "uvt": "#9467bd",
    "gpt": "#8c564b",
    "fvt": "#17becf",
}
_EXTRA_COLORS = ("#e377c2", "#7f7f7f", "#bcbd22", "#393b79", "#637939")


class FluidMismatch(ValueError):
    """Punkt eines anderen Fluids als die bestehende Auswahl."""


class NotAddable(ValueError):
    """Zustand nicht rekonstruierbar (not_unique, out_of_range, numerical_failure, unsupported_eos)."""


ADDABLE_STATUS = ("valid", "partially_valid", "synthetic_substitute")


def color_for_diagram(diagram_id):
    if diagram_id in DIAGRAM_COLORS:
        return DIAGRAM_COLORS[diagram_id]
    return _EXTRA_COLORS[sum(map(ord, diagram_id)) % len(_EXTRA_COLORS)]


@dataclass
class SelectedPoint:
    selection_index: int                 # laufende Punktnummer (1, 2, 3, ...), nie wiederverwendet
    selection_order: int                 # Rang in der Verbindungsreihenfolge
    color: str
    source_diagram: str
    state_id: str
    state: ThermoState
    display_xyz: Optional[tuple] = None  # Position im Ursprungsdiagramm (Markierung)
    original: dict = field(default_factory=dict)  # ausgewaehlte Originalkoordinaten (Quelleinheiten)

    @property
    def label(self):
        return f"Punkt {self.selection_index}"

    def description(self):
        """'aus p-v-T bei T=..., p=..., v=...' (Originalkoordinaten der Auswahl)."""
        parts = ", ".join(f"{k}={v:.6g}" for k, v in self.original.items() if v is not None)
        return f"aus {self.source_diagram} bei {parts}" if parts else f"aus {self.source_diagram}"


class SelectionManager:
    """Verwaltet die Vergleichspunkte; Aenderungen werden ueber Listener gemeldet."""

    def __init__(self):
        self._points = []
        self._next_index = 1
        self._next_order = 1
        self._listeners = []
        self.fluid = None

    # ---------------------------------------------------------------- Listener
    def subscribe(self, callback: Callable[[str], None]):
        self._listeners.append(callback)

    def unsubscribe(self, callback: Callable[[str], None]):
        if callback in self._listeners:
            self._listeners.remove(callback)
            return True
        return False

    def _notify(self, event):
        for callback in list(self._listeners):
            callback(event)

    # ---------------------------------------------------------------- Zugriff
    @property
    def points(self):
        """Punkte in Verbindungsreihenfolge (selection_order)."""
        return sorted(self._points, key=lambda p: p.selection_order)

    def __len__(self):
        return len(self._points)

    def get(self, selection_index):
        for point in self._points:
            if point.selection_index == selection_index:
                return point
        raise KeyError(selection_index)

    def points_of(self, diagram_id):
        return [p for p in self.points if p.source_diagram == diagram_id]

    # ---------------------------------------------------------------- Aenderungen
    def add(self, state: ThermoState, display_xyz=None, original=None):
        """Fuegt einen Zustand hinzu. Anderes Fluid als die bestehende Auswahl: FluidMismatch."""
        if state.validity not in ADDABLE_STATUS:
            raise NotAddable(f"Zustand nicht verwendbar ({state.validity}): " + "; ".join(state.notes))
        if self._points and state.fluid != self.fluid:
            raise FluidMismatch(f"Vergleich nur fuer ein Fluid: Auswahl enthaelt {self.fluid}, Punkt ist {state.fluid}")
        self.fluid = state.fluid
        point = SelectedPoint(self._next_index, self._next_order, color_for_diagram(state.source_diagram),
                              state.source_diagram, state.id, state, display_xyz, dict(original or {}))
        self._next_index += 1
        self._next_order += 1
        self._points.append(point)
        self._notify("added")
        return point

    def remove(self, selection_index):
        point = self.get(selection_index)
        self._points.remove(point)
        if not self._points:
            self.fluid = None
        self._notify("removed")
        return point

    def reset(self):
        """Alle Punkte entfernen; die Nummerierung laeuft weiter (keine Wiederverwendung von IDs)."""
        if self._points:
            self._points.clear()
            self.fluid = None
            self._notify("reset")

    def move(self, selection_index, delta):
        """Verschiebt einen Punkt in der Verbindungsreihenfolge (delta -1 nach vorn, +1 nach hinten)."""
        ordered = self.points
        position = [p.selection_index for p in ordered].index(selection_index)
        target = position + int(delta)
        if not 0 <= target < len(ordered):
            return False
        ordered[position], ordered[target] = ordered[target], ordered[position]
        for rank, point in enumerate(ordered, start=1):
            point.selection_order = rank
        self._next_order = len(ordered) + 1
        self._notify("moved")
        return True

    def clear_for_fluid_change(self, new_fluid):
        """Beim Fluidwechsel: Auswahl loeschen (Vergleich nie ueber Fluide hinweg)."""
        if self._points and new_fluid != self.fluid:
            self.reset()
            return True
        return False
