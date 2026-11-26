
"""
Datei: state_projection.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Projektion thermodynamischer Zustandspunkte
in die Koordinatensysteme anderer 3D-Diagramme.

Aufgaben:
- Auswahl der benötigten Zustandsgrößen
- Umrechnung molarer und massenspezifischer Größen
- Transformation in diagrammspezifische 3D-Koordinaten
- Prüfung der Darstellbarkeit im Zieldiagramm
- Bereitstellung von Projektionsdaten und Fehlermeldungen

Abhängigkeiten:
- thermo_state.py: Zustandsdaten und Herkunftskennzeichnung
- state_store.py: Diagrammspezifische Koordinatenadapter
- diagram_axes.py: Achsenskalierung und Darstellungsregeln
- NumPy: Numerische Verarbeitung

Verwendet von:
- gui.py: Projektion ausgewählter Zustände
"""

from dataclasses import dataclass

import numpy as np

import diagram_axes
import thermo_state as ts
from state_store import adapter_for


_FIELD_BY_AXIS = {
    "T": "T",
    "p": "p",
    "v": "v_molar",
    "h": "h_molar",
    "s": "s_molar",
    "u": "u_molar",
    "g": "g_molar",
    "f": "f_molar",
}


@dataclass(frozen=True)
class ProjectionResult:
    target_diagram: str
    native_values: tuple = ()
    display_xyz: tuple = ()
    valid: bool = False
    reason: str = ""


def _failure(target_diagram, reason):
    return ProjectionResult(target_diagram=target_diagram, reason=reason)


def project_state(state, target_diagram, target_data):
    """Map a complete resolved state through the target adapter without another HEOS flash."""
    if state.eos_source == ts.SOURCE_SYNTHETIC:
        return _failure(target_diagram, "Synthetische Ersatzwerte sind nicht als HEOS-Zustand projizierbar.")
    if not state.is_usable:
        return _failure(target_diagram, f"Zustand nicht verwendbar ({state.validity}).")
    if (state.source_diagram == "gpt" and state.stability == "saturated"
            and state.Q is None and target_diagram != "gpt"):
        return _failure(target_diagram, "Dampfgehalt Q ist auf der g-p-T-Sättigungslinie nicht eindeutig.")

    try:
        adapter = adapter_for(target_diagram)
    except ValueError as exc:
        return _failure(target_diagram, str(exc))

    molar_mass = state.provenance.get("molar_mass_kg_per_mol")
    if adapter.mass_axes and (molar_mass is None or not np.isfinite(molar_mass) or molar_mass <= 0):
        return _failure(target_diagram, "Molare Masse fehlt fuer die massenspezifischen Zielachsen.")

    native = []
    for axis in adapter.axes:
        field = _FIELD_BY_AXIS[axis]
        value = getattr(state, field, None)
        if value is None or not np.isfinite(value):
            return _failure(target_diagram, f"Zustandsgröße {axis} fehlt oder ist nicht endlich.")
        if axis in adapter.mass_axes:
            value = float(value) / float(molar_mass)
        native.append(float(value))

    by_axis = dict(zip(adapter.axes, native))
    if "v" in by_axis and by_axis["v"] <= 0:
        return _failure(target_diagram, "Volumen muss positiv sein.")
    if "p" in by_axis and by_axis["p"] <= 0:
        pressure_index = next(i for i, spec in enumerate(diagram_axes.axes_for(target_diagram)) if spec.key == "p")
        if diagram_axes.axis_scale(target_diagram, pressure_index, target_data) == "log":
            return _failure(target_diagram, "Druck muss auf einer logarithmischen Achse positiv sein.")

    try:
        xyz = np.asarray(adapter.to_display(np.asarray(native, dtype=float), target_data), dtype=float)
    except (KeyError, TypeError, ValueError, FloatingPointError) as exc:
        return _failure(target_diagram, f"Zielkoordinaten nicht berechenbar: {exc}")
    if xyz.shape != (3,) or not np.isfinite(xyz).all():
        return _failure(target_diagram, "Zieltransformation lieferte keine endlichen 3D-Koordinaten.")

    return ProjectionResult(target_diagram, tuple(native), tuple(float(value) for value in xyz), True)
