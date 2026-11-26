"""
Datei: diagram_coords.py
Projekt: ThermoSurface
Kategorie: Visualisierung / Koordinaten

Beschreibung:
Transformiert Zustandsgrößen der sieben Diagramme in 3D-Anzeigekoordinaten.

Aufgaben:
- Stellt Vorwärts- und Rücktransformationen für die Diagrammtypen bereit.
- Berücksichtigt normierte sowie logarithmische Druck- und Volumenachsen.

Abhängigkeiten:
- Projektmodule (direkt): zustandsflaeche.py
- Bibliotheken: numpy

Verwendet von:
- Projektmodule: diagram_axes.py, diagram_viewport.py, fvt_tangent.py, visualization.py
"""

import numpy as np

from zustandsflaeche import PRESSURE_AXES, pressure_to_z, z_to_pressure


# ---------------------------------------------------------------- Vorwaertstransformationen

def to_xyz(points, critical, pressure_axis="log"):
    """Einheitliche p-v-T-Koordinaten fuer ALLE Actors, inklusive Phasengrenzen.

    p_raw bleibt unveraendert; nur die Darstellungskoordinate z wird transformiert
    (zustandsflaeche.pressure_to_z: log nur fuer p>0, asinh/symlog/linear fuer alle p).
    """
    Tc, pc, vc = critical
    pts = np.asarray(points, dtype=float)
    if pressure_axis not in PRESSURE_AXES:
        raise ValueError("pressure_axis muss eines von " + ", ".join(PRESSURE_AXES) + " sein")
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        z = pressure_to_z(pts[..., 2], pc, pressure_axis)
        return np.stack((np.log10(pts[..., 1] / vc), pts[..., 0] / Tc, z), axis=-1)


def to_xyz_pht(points, critical):
    """(T, h, p) -> PyVista-XYZ; critical: (Tc, pc, hc)."""
    Tc, pc, hc = critical
    points = np.asarray(points, dtype=float)
    T = points[..., 0]
    h = points[..., 1]
    p = points[..., 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.stack(((h - hc) / 1e6, T / Tc, np.log10(p / pc)), axis=-1)


def to_xyz_pst(points, critical):
    """(T [K], s [J/(kg K)], p [Pa]) -> (delta_s/1e3, T/Tc, log10(p/pc))."""
    Tc, pc, sc = critical
    pts = np.asarray(points, dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        return np.stack(((pts[..., 1] - sc) / 1e3,
                         pts[..., 0] / Tc,
                         np.log10(pts[..., 2] / pc)), axis=-1)


def to_xyz_hsp(points, critical):
    """(h [J/kg], s [J/kg/K], p [Pa]) -> (delta h [MJ/kg], delta s [kJ/kg/K], log p_r)."""
    hc, sc, pc = critical
    points = np.asarray(points, dtype=np.float64)

    with np.errstate(divide="ignore", invalid="ignore"):
        return np.stack(((points[..., 0] - hc) / 1e6, (points[..., 1] - sc) / 1e3, np.log10(points[..., 2] / pc)), axis=-1)


def to_xyz_uvt(points, critical):
    """(T [K], v [m3/mol], u [J/kg]) -> (log v_r, T_r, delta u [MJ/kg])."""
    Tc, vc, uc = critical
    points = np.asarray(points, dtype=np.float64)

    with np.errstate(divide="ignore", invalid="ignore"):
        return np.stack((np.log10(points[..., 1] / vc), points[..., 0] / Tc, (points[..., 2] - uc) / 1e6), axis=-1)


def to_xyz_gpt(points, critical, pressure_axis="log", p_window=None):
    """(T [K], p [Pa], g [J/kg]) -> (p-Achse (log10 p_r), T_r, delta g [MJ/kg]); critical: (Tc, pc, gc).

    p_window=(p_min, p_max): nur fuer die Log-Achse; Knoten ausserhalb des berechneten Druckfensters der
    Gleichgewichtsflaeche (z. B. homogene Zustaende bei p -> 0+) werden nicht dargestellt (NaN), die Daten bleiben
    unveraendert.
    """
    Tc, pc, gc = critical
    pts = np.asarray(points, dtype=np.float64)

    if pressure_axis not in PRESSURE_AXES:
        raise ValueError("pressure_axis muss eines von " + ", ".join(PRESSURE_AXES) + " sein")

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        z = pressure_to_z(pts[..., 1], pc, pressure_axis)

        if p_window is not None and pressure_axis == "log":
            lo, hi = p_window
            z = np.where((pts[..., 1] >= lo * (1.0 - 1e-9)) & (pts[..., 1] <= hi * (1.0 + 1e-9)), z, np.nan)

        return np.stack((z, pts[..., 0] / Tc, (pts[..., 2] - gc) / 1e6), axis=-1)


def to_xyz_fvt(points, critical):
    """(T [K], v [m3/mol], f [J/kg]) -> (log v_r, T_r, delta f [MJ/kg]); critical: (Tc, vc, fc)."""
    Tc, vc, fc = critical
    pts = np.asarray(points, dtype=np.float64)

    with np.errstate(divide="ignore", invalid="ignore"):
        return np.stack((np.log10(pts[..., 1] / vc), pts[..., 0] / Tc, (pts[..., 2] - fc) / 1e6), axis=-1)


# ---------------------------------------------------------------- Umkehrung (Darstellung -> native Werte)

def from_xyz(xyz, critical, pressure_axis="log"):
    """Umkehrung von to_xyz: -> (T, v, p)."""
    Tc, pc, vc = critical
    xyz = np.asarray(xyz, dtype=float)
    return np.stack((xyz[..., 1] * Tc, vc * 10.0 ** xyz[..., 0], z_to_pressure(xyz[..., 2], pc, pressure_axis)), axis=-1)


def from_xyz_pht(xyz, critical):
    """Umkehrung von to_xyz_pht: -> (T, h, p)."""
    Tc, pc, hc = critical
    xyz = np.asarray(xyz, dtype=float)
    return np.stack((xyz[..., 1] * Tc, hc + xyz[..., 0] * 1e6, pc * 10.0 ** xyz[..., 2]), axis=-1)


def from_xyz_pst(xyz, critical):
    """Umkehrung von to_xyz_pst: -> (T, s, p)."""
    Tc, pc, sc = critical
    xyz = np.asarray(xyz, dtype=float)
    return np.stack((xyz[..., 1] * Tc, sc + xyz[..., 0] * 1e3, pc * 10.0 ** xyz[..., 2]), axis=-1)


def from_xyz_hsp(xyz, critical):
    """Umkehrung von to_xyz_hsp: -> (h, s, p)."""
    hc, sc, pc = critical
    xyz = np.asarray(xyz, dtype=float)
    return np.stack((hc + xyz[..., 0] * 1e6, sc + xyz[..., 1] * 1e3, pc * 10.0 ** xyz[..., 2]), axis=-1)


def from_xyz_uvt(xyz, critical):
    """Umkehrung von to_xyz_uvt: -> (T, v, u)."""
    Tc, vc, uc = critical
    xyz = np.asarray(xyz, dtype=float)
    return np.stack((xyz[..., 1] * Tc, vc * 10.0 ** xyz[..., 0], uc + xyz[..., 2] * 1e6), axis=-1)


def from_xyz_gpt(xyz, critical, pressure_axis="log"):
    """Umkehrung von to_xyz_gpt: -> (T, p, g)."""
    Tc, pc, gc = critical
    xyz = np.asarray(xyz, dtype=float)
    return np.stack((xyz[..., 1] * Tc, z_to_pressure(xyz[..., 0], pc, pressure_axis), gc + xyz[..., 2] * 1e6), axis=-1)


def from_xyz_fvt(xyz, critical):
    """Umkehrung von to_xyz_fvt: -> (T, v, f)."""
    Tc, vc, fc = critical
    xyz = np.asarray(xyz, dtype=float)
    return np.stack((xyz[..., 1] * Tc, vc * 10.0 ** xyz[..., 0], fc + xyz[..., 2] * 1e6), axis=-1)
