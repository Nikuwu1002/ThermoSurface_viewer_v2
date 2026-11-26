
"""
Datei: viewport_recovery.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Hilfsmodul zur Verwaltung und Prüfung
des Renderzustands einzelner 3D-Viewports.

Aufgaben:
- Definition der möglichen Renderzustände
- Prüfung der Wiederverwendbarkeit vorhandener Darstellungen
- Vermeidung unnötiger Renderaufrufe bei unveränderten Daten

Abhängigkeiten:
- Keine projektspezifischen oder externen Module

Verwendet von:
- diagram_viewport.py: Renderzustandsverwaltung
- gui.py: Prüfung der Wiederverwendbarkeit
"""


READY = "ready"
DIRTY = "dirty"
RENDERING = "rendering"
FAILED = "failed"
DISPOSED = "disposed"
EMPTY = "empty"


def can_skip_render(has_result, current_key, requested_key, render_state):
    return bool(has_result and current_key == requested_key and render_state == READY)
