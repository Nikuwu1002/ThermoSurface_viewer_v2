"""
Datei: hover_data.py
Projekt: ThermoSurface
Kategorie: Visualisierung / Hover-Daten

Beschreibung:
Speichert und formatiert physikalische Mesh-Daten für 3D-Hoveranzeigen.

Aufgaben:
- Hinterlegt Zustandsgrößen, Diagrammtypen und Datenherkunft im Mesh.
- Formatiert gepickte Gitter- und interpolierte Schnittpunkte.

Abhängigkeiten:
- Projektmodule (direkt): keine
- Bibliotheken: numpy

Verwendet von:
- Projektmodule: visualization.py
"""

import numpy as np


_AXES = {
    "pvt": (("T", "K", 1.0), ("v", "m³/mol", 1.0), ("p", "MPa", 1e-6)),
    "pht": (("T", "K", 1.0), ("h", "kJ/kg", 1e-3), ("p", "MPa", 1e-6)),
    "pst": (("T", "K", 1.0), ("s", "kJ/(kg K)", 1e-3), ("p", "MPa", 1e-6)),
    "hsp": (("h", "kJ/kg", 1e-3), ("s", "kJ/(kg K)", 1e-3), ("p", "MPa", 1e-6)),
    "uvt": (("T", "K", 1.0), ("v", "m³/mol", 1.0), ("u", "kJ/kg", 1e-3)),
    "gpt": (("T", "K", 1.0), ("p", "MPa", 1e-6), ("g", "kJ/kg", 1e-3)),
    "fvt": (("T", "K", 1.0), ("v", "m³/mol", 1.0), ("f", "kJ/kg", 1e-3)),
}
DIAGRAM_CODES = {diagram_id: code for code, diagram_id in enumerate(_AXES)}
CODE_DIAGRAMS = {code: diagram_id for diagram_id, code in DIAGRAM_CODES.items()}

SOURCE_HEOS = 0
SOURCE_HOMOGENEOUS = 1
SOURCE_SYNTHETIC = 2
SOURCE_CODES = {
    "HEOS": SOURCE_HEOS,
    "HEOS_homogeneous": SOURCE_HOMOGENEOUS,
    "synthetic_substitute": SOURCE_SYNTHETIC,
}
CODE_SOURCES = {value: key for key, value in SOURCE_CODES.items()}


def attach_hover_metadata(mesh, native_points, diagram_id, source_kind="HEOS"):
    """Attach physical axis coordinates before any display-coordinate normalization."""
    if diagram_id not in DIAGRAM_CODES:
        raise ValueError(f"Unbekannter Hover-Diagrammtyp: {diagram_id}")
    if source_kind not in SOURCE_CODES:
        raise ValueError(f"Unbekannte Hover-Datenquelle: {source_kind}")
    values = np.asarray(native_points, dtype=float).reshape((-1, 3))
    if len(values) != mesh.n_points or not np.isfinite(values).all():
        raise ValueError("Hover-Koordinaten müssen endlich sein und jedem Mesh-Punkt entsprechen.")
    mesh.point_data["ts_hover_native"] = values
    mesh.point_data["ts_hover_diagram"] = np.full(mesh.n_points, DIAGRAM_CODES[diagram_id], dtype=np.uint8)
    mesh.point_data["ts_hover_source"] = np.full(mesh.n_points, SOURCE_CODES[source_kind], dtype=np.uint8)
    return mesh


def format_hover_values(diagram_id, values, source="HEOS", interpolated=False):
    """Format the nearest source mesh vertex, explicitly identifying it as a grid point."""
    axes = _AXES.get(diagram_id)
    source_label = {
        "HEOS": "HEOS",
        "HEOS_homogeneous": "HEOS homogen",
        "synthetic_substitute": "Modell (kein HEOS)",
    }.get(source)
    try:
        native = np.asarray(values, dtype=float)
    except (TypeError, ValueError):
        return None
    if axes is None or source_label is None or native.shape != (3,) or not np.isfinite(native).all():
        return None
    items = (f"{label}={value * factor:g} {unit}"
             for value, (label, unit, factor) in zip(native, axes))
    point_label = "interpolierter Schnittpunkt" if interpolated else "Gitterpunkt"
    return f"{diagram_id}: " + ", ".join(items) + f" ({point_label}; {source_label})"


def format_mesh_hover_point(mesh, point_id):
    """Read and format physical metadata from a picked VTK/PyVista mesh point."""
    try:
        point_data = mesh.GetPointData()
        native = point_data.GetArray("ts_hover_native")
        diagrams = point_data.GetArray("ts_hover_diagram")
        sources = point_data.GetArray("ts_hover_source")
        if native is None or diagrams is None or sources is None:
            return None
        diagram_id = CODE_DIAGRAMS.get(int(diagrams.GetValue(point_id)))
        source = CODE_SOURCES.get(int(sources.GetValue(point_id)))
        clipped = point_data.GetArray("ts_hover_interpolated")
        return format_hover_values(diagram_id, native.GetTuple(point_id), source,
                                   interpolated=clipped is not None and bool(clipped.GetValue(point_id)))
    except (AttributeError, IndexError, TypeError, ValueError):
        return None
