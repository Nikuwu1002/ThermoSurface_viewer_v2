
"""
Datei: tile_layout.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Berechnung einer automatischen Kachelanordnung
für die gleichzeitige Darstellung mehrerer 3D-Diagramme.

Aufgaben:
- Bestimmung der benötigten Zeilen und Spalten
- Anpassung der Kachelanordnung an die Diagrammanzahl
- Bereitstellung der Layoutparameter für den Workspace

Abhängigkeiten:
- Python-Standardbibliothek: math
- Keine projektspezifischen Module

Verwendet von:
- diagram_workspace.py: Anordnung der Diagramm-Viewports
"""

from math import ceil


def tile_grid_shape(count):
    """1 = 1x1, 2 = 1x2, 3-4 = 2x2, 5-6 = 2x3, darueber 3 Spalten (erweiterbar)."""
    if count <= 1:
        return 1, 1
    if count == 2:
        return 1, 2
    if count <= 4:
        return 2, 2
    return ceil(count / 3), 3
