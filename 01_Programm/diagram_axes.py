"""
Datei: diagram_axes.py
Projekt: ThermoSurface
Kategorie: Visualisierung / Achsenlogik

Beschreibung:
Definiert die physikalischen Achsen und Eingabeeinheiten der sieben Diagramme.

Aufgaben:
- Ordnet Darstellungsachsen den Zustandsgrößen und Skalierungen zu.
- Prüft Viewportgrenzen und berechnet Tickwerte und Koordinaten.
Intern SI; Achsentitel und Eingaben verwenden ausgewiesene Anzeigeeinheiten.

Abhängigkeiten:
- Projektmodule (direkt): diagram_coords.py, zustandsflaeche.py
- Bibliotheken: numpy

Verwendet von:
- Projektmodule: axes_view.py, diagram_viewport.py, gui.py
"""

import math
from dataclasses import dataclass

import numpy as np

import diagram_coords as dc
from zustandsflaeche import PRESSURE_AXES


class ViewportBoundsError(ValueError):
    """Ungueltige manuelle Viewportgrenzen (Validierungsmeldung fuer die GUI)."""


@dataclass(frozen=True)
class AxisSpec:
    key: str          # T, v, p, h, s, u, g, f
    name: str         # Klartextname
    unit: str         # Anzeigeeinheit
    scale: float      # Anzeigewert = SI-Wert * scale
    column: int       # Spalte der nativen Zustandsdaten dieses Diagramms
    positive: bool    # physikalisch nur > 0 zulaessig (T, v)

    @property
    def title(self):
        return f"{self.key} [{self.unit}]"


def _ax(key, column):
    name, unit, scale, positive = {
        "T": ("Temperatur", "K", 1.0, True),
        "v": ("Molares Volumen", "m³/mol", 1.0, True),
        "p": ("Druck", "MPa", 1e-6, False),
        "h": ("Spezifische Enthalpie", "kJ/kg", 1e-3, False),
        "s": ("Spezifische Entropie", "kJ/(kg K)", 1e-3, False),
        "u": ("Spezifische innere Energie", "kJ/kg", 1e-3, False),
        "g": ("Spezifische freie Enthalpie", "kJ/kg", 1e-3, False),
        "f": ("Spezifische freie Energie", "kJ/kg", 1e-3, False),
    }[key]
    return AxisSpec(key, name, unit, scale, column, positive)


# Koordinatenvertrag (TS-VP01, Abschn. 1.3): Darstellungsachsen (x, y, z) -> Zustandsgroesse und Spalte der
# nativen Daten. Native Spaltenreihenfolge: pvt (T,v,p), pht (T,h,p), pst (T,s,p), hsp (h,s,p), uvt (T,v,u),
# gpt (T,p,g), fvt (T,v,f).
DIAGRAM_AXES = {
    "pvt": (_ax("v", 1), _ax("T", 0), _ax("p", 2)),
    "pht": (_ax("h", 1), _ax("T", 0), _ax("p", 2)),
    "pst": (_ax("s", 1), _ax("T", 0), _ax("p", 2)),
    "hsp": (_ax("h", 0), _ax("s", 1), _ax("p", 2)),
    "uvt": (_ax("v", 1), _ax("T", 0), _ax("u", 2)),
    "gpt": (_ax("p", 1), _ax("T", 0), _ax("g", 2)),
    "fvt": (_ax("v", 1), _ax("T", 0), _ax("f", 2)),
}
AXIS_LETTERS = ("x", "y", "z")


def axes_for(diagram_id):
    try:
        return DIAGRAM_AXES[diagram_id]
    except KeyError:
        raise ValueError(f"Unbekannter Diagrammtyp: {diagram_id}") from None


def _native_reference(diagram_id, daten):
    """Kritischer Zustand in nativer Spaltenreihenfolge (Fuellwerte fuer die nicht variierten Spalten)."""
    c = daten["critical_point"]
    return {"pvt": (c[0], c[2], c[1]), "pht": (c[0], c[2], c[1]), "pst": (c[0], c[2], c[1]), "hsp": tuple(c),
            "uvt": tuple(c), "gpt": tuple(c), "fvt": tuple(c)}[diagram_id]


def pressure_axis_of(daten):
    axis = (daten or {}).get("pressure_axis", "log")
    return axis if axis in PRESSURE_AXES else "log"


def axis_scale(diagram_id, index, daten):
    """Skalenart der Darstellungsachse `index`: 'linear' | 'log' | 'asinh' | 'symlog'."""
    spec = axes_for(diagram_id)[index]
    if spec.key == "v":
        return "log"
    if spec.key == "p":
        return pressure_axis_of(daten) if diagram_id in ("pvt", "gpt") else "log"
    return "linear"


def _forward(diagram_id, daten, points):
    critical = daten["critical_point"]
    if diagram_id == "pvt":
        return dc.to_xyz(points, critical, pressure_axis_of(daten))
    if diagram_id == "gpt":
        return dc.to_xyz_gpt(points, critical, pressure_axis_of(daten))
    return {"pht": dc.to_xyz_pht, "pst": dc.to_xyz_pst, "hsp": dc.to_xyz_hsp, "uvt": dc.to_xyz_uvt,
            "fvt": dc.to_xyz_fvt}[diagram_id](points, critical)


def _inverse(diagram_id, daten, xyz):
    critical = daten["critical_point"]
    if diagram_id == "pvt":
        return dc.from_xyz(xyz, critical, pressure_axis_of(daten))
    if diagram_id == "gpt":
        return dc.from_xyz_gpt(xyz, critical, pressure_axis_of(daten))
    return {"pht": dc.from_xyz_pht, "pst": dc.from_xyz_pst, "hsp": dc.from_xyz_hsp, "uvt": dc.from_xyz_uvt,
            "fvt": dc.from_xyz_fvt}[diagram_id](xyz, critical)


def to_display(diagram_id, daten, index, values_si):
    """SI-Werte der Achse `index` -> Darstellungskoordinate (dieselbe Vorwaertstransformation wie die Meshes)."""
    spec = axes_for(diagram_id)[index]
    values = np.atleast_1d(np.asarray(values_si, dtype=float))
    points = np.tile(np.asarray(_native_reference(diagram_id, daten), dtype=float), (len(values), 1))
    points[:, spec.column] = values
    with np.errstate(all="ignore"):
        return _forward(diagram_id, daten, points)[:, index]


def to_physical(diagram_id, daten, index, coords):
    """Darstellungskoordinaten der Achse `index` -> SI-Werte (Umkehrung von `to_display`)."""
    spec = axes_for(diagram_id)[index]
    coords = np.atleast_1d(np.asarray(coords, dtype=float))
    xyz = np.zeros((len(coords), 3))
    reference = _forward(diagram_id, daten, np.asarray([_native_reference(diagram_id, daten)], dtype=float))[0]
    xyz[:] = reference
    xyz[:, index] = coords
    with np.errstate(all="ignore"):
        return _inverse(diagram_id, daten, xyz)[:, spec.column]


# ---------------------------------------------------------------- Grenzenpruefung

def validate_bounds(diagram_id, daten, bounds_display):
    """Prueft sechs Grenzen in ANZEIGEEINHEITEN: [(lo, hi)_x, (lo, hi)_y, (lo, hi)_z].

    Rueckgabe: (bounds_si, box) mit bounds_si = [(lo, hi), ...] in SI und box = (xmin, xmax, ymin, ymax, zmin, zmax)
    in Darstellungskoordinaten. Wirft ViewportBoundsError mit lesbarer Meldung.
    """
    specs = axes_for(diagram_id)
    bounds_si, box = [], []
    for index, ((lo, hi), spec) in enumerate(zip(bounds_display, specs)):
        letter = AXIS_LETTERS[index]
        label = f"{letter}-Achse ({spec.key})"
        for name, value in (("Start", lo), ("Stop", hi)):
            if value is None or not math.isfinite(value):
                raise ViewportBoundsError(f"{label}: {name} muss eine endliche Zahl sein.")
        if not lo < hi:
            raise ViewportBoundsError(f"{label}: Start ({lo:g}) muss kleiner als Stop ({hi:g}) sein.")
        lo_si, hi_si = lo / spec.scale, hi / spec.scale
        scale_kind = axis_scale(diagram_id, index, daten)
        if spec.positive and lo_si <= 0:
            raise ViewportBoundsError(f"{label}: {spec.name} muss > 0 sein (Start = {lo:g} {spec.unit}).")
        if scale_kind == "log" and lo_si <= 0:
            raise ViewportBoundsError(
                f"{label}: Logarithmische Achse erfordert Start > 0 (Start = {lo:g} {spec.unit}); fuer p <= 0 "
                "pressure_axis 'asinh', 'symlog' oder 'linear' verwenden.")
        coords = to_display(diagram_id, daten, index, [lo_si, hi_si])
        if not np.isfinite(coords).all() or not coords[0] < coords[1]:
            raise ViewportBoundsError(f"{label}: Grenzen {lo:g}..{hi:g} {spec.unit} sind nicht darstellbar.")
        bounds_si.append((lo_si, hi_si))
        box.extend((float(coords[0]), float(coords[1])))
    return bounds_si, tuple(box)


def bounds_to_box(diagram_id, daten, bounds_si):
    """Gespeicherte SI-Grenzen -> VTK-Box (xmin, xmax, ymin, ymax, zmin, zmax) fuer die aktuellen Daten."""
    box = []
    for index, (lo, hi) in enumerate(bounds_si):
        coords = to_display(diagram_id, daten, index, [lo, hi])
        if not np.isfinite(coords).all() or not coords[0] < coords[1]:
            raise ViewportBoundsError(
                f"Gespeicherte Grenzen der {AXIS_LETTERS[index]}-Achse sind fuer die aktuellen Daten nicht darstellbar.")
        box.extend((float(coords[0]), float(coords[1])))
    return tuple(box)


def display_bounds_to_physical(diagram_id, daten, box):
    """VTK-Box -> physikalische Grenzen in Anzeigeeinheiten [(lo, hi)_x, (lo, hi)_y, (lo, hi)_z]."""
    out = []
    for index, spec in enumerate(axes_for(diagram_id)):
        lo, hi = to_physical(diagram_id, daten, index, [box[2 * index], box[2 * index + 1]]) * spec.scale
        out.append((float(lo), float(hi)))
    return out


# ---------------------------------------------------------------- Ticks (physikalische Werte)

def format_tick(value):
    if value == 0:
        return "0"
    magnitude = abs(value)
    if magnitude >= 1e6 or magnitude < 1e-3:
        mantissa, exponent = f"{value:.4e}".split("e")
        mantissa = mantissa.rstrip("0").rstrip(".")
        return f"{mantissa}e{int(exponent)}"
    return f"{value:.6g}"


def _linear_ticks(lo, hi, target):
    span = hi - lo
    if not span > 0:
        return [lo]
    raw = span / max(target, 1)
    exponent = math.floor(math.log10(raw))
    fraction = raw / 10.0 ** exponent
    step = next(m for m in (1.0, 2.0, 2.5, 5.0, 10.0) if fraction <= m + 1e-12) * 10.0 ** exponent
    first = math.ceil(lo / step - 1e-9)
    last = math.floor(hi / step + 1e-9)
    return [k * step for k in range(first, last + 1)]


def _decade_ticks(lo, hi, mantissas):
    values = []
    for exponent in range(math.floor(math.log10(lo)) - 1, math.ceil(math.log10(hi)) + 2):
        for m in mantissas:
            value = m * 10.0 ** exponent
            if lo * (1 - 1e-12) <= value <= hi * (1 + 1e-12):
                values.append(value)
    return values


def _thin(values, limit):
    values = list(values)
    while len(values) > limit:
        values = values[::2]
    return values


def physical_ticks(scale_kind, lo, hi, target=6, min_abs=0.0):
    """Tickwerte (in Anzeigeeinheiten) fuer den physikalischen Bereich [lo, hi] der Achse.

    linear: 1-2-2.5-5-Raster; log: Dekaden (bei weniger als 3 Dekaden zusaetzlich 2 und 5);
    asinh/symlog: 0 sowie symmetrische Dekaden bzw. 1-2-5-Werte der positiven und negativen Seite; `min_abs`
    (Anzeigeeinheiten) blendet Ticks im nahezu linearen Bereich um 0 aus, die dort optisch zusammenfielen.
    """
    if not (math.isfinite(lo) and math.isfinite(hi)) or not lo < hi:
        return []
    if scale_kind == "linear":
        return _linear_ticks(lo, hi, target)
    if scale_kind == "log":
        if lo <= 0:
            return []
        values = _decade_ticks(lo, hi, (1.0,))
        if len(values) < 3:
            values = _decade_ticks(lo, hi, (1.0, 2.0, 5.0))
        return _thin(values, 2 * target)
    if scale_kind in ("asinh", "symlog"):
        positive_hi = hi if hi > 0 else 0.0
        negative_hi = -lo if lo < 0 else 0.0
        floor_value = max(min_abs, 1e-30)

        def candidates(mantissas):
            pos = _decade_ticks(floor_value, positive_hi, mantissas) if positive_hi > floor_value else []
            neg = [-x for x in _decade_ticks(floor_value, negative_hi, mantissas)] if negative_hi > floor_value else []
            return sorted(neg + ([0.0] if lo <= 0 <= hi else []) + pos)

        values = candidates((1.0,))
        if len([x for x in values if x != 0.0]) < 5:
            values = candidates((1.0, 2.0, 5.0))
        zero = [0.0] if 0.0 in values else []
        rest = [x for x in values if x != 0.0]
        while len(rest) + len(zero) > 2 * target:
            rest = rest[::2]
        return sorted(rest + zero)
    raise ValueError(f"Unbekannte Skalenart: {scale_kind}")


def tick_floor(diagram_id, index, daten):
    """Kleinster angezeigter |Tick| (Anzeigeeinheiten) fuer asinh/symlog-Druckachsen, sonst 0."""
    spec = axes_for(diagram_id)[index]
    if spec.key == "p" and axis_scale(diagram_id, index, daten) in ("asinh", "symlog"):
        return 1e-2 * float(daten["critical_point"][1]) * spec.scale
    return 0.0


def axis_title(diagram_id, index, daten):
    spec = axes_for(diagram_id)[index]
    kind = axis_scale(diagram_id, index, daten)
    suffix = {"linear": "", "log": " (log)", "asinh": " (asinh)", "symlog": " (symlog)"}[kind]
    return spec.title + suffix
