"""
Datei: fvt_tangent.py
Projekt: ThermoSurface
Kategorie: Thermodynamik / Tangentenkonstruktion

Beschreibung:
Bestimmt die gemeinsame Tangente der f-v-Isotherme aus HEOS-Sättigungsdaten.

Aufgaben:
- Berechnet Sättigungspunkte, Tangentensteigung und Gleichgewichtsprüfung.
- Liefert 2D-Kurvendaten und transformierte 3D-Tangentenpunkte.

Abhängigkeiten:
- Projektmodule (direkt): coolprop_heos_common.py, coolprop_heos_fvt.py,
  diagram_coords.py
- Bibliotheken: CoolProp, numpy

Verwendet von:
- Projektmodule: diagram_viewport.py, fvt_slice_view.py, gui.py
"""

from dataclasses import dataclass

import CoolProp.CoolProp as CP
from CoolProp import AbstractState
import numpy as np

import diagram_coords


@dataclass(frozen=True)
class FvtTangent:
    fluid: str
    temperature: float
    p_sat: float
    v_l: float
    v_v: float
    f_l_mass: float
    f_v_mass: float
    f_l_molar: float
    f_v_molar: float
    g_l_mass: float
    g_v_mass: float
    g_l_molar: float
    g_v_molar: float
    mu: float
    mu_heos: float
    line_v: np.ndarray
    line_f_molar: np.ndarray
    tangent_pressure: float
    pressure_relative_error: float
    chemical_potential_delta: float


def tangent_display_points(tangent, fvt_data):
    """Transform molar tangent points using the exact mass-specific f-v-T renderer contract."""
    molar_mass = float(fvt_data["meta"]["molar_mass"])
    if not np.isfinite(molar_mass) or molar_mass <= 0:
        raise ValueError("f-v-T-Daten enthalten keine gueltige molare Masse.")
    native = np.column_stack((
        np.full(len(tangent.line_v), tangent.temperature),
        tangent.line_v,
        tangent.line_f_molar / molar_mass,
    ))
    return diagram_coords.to_xyz_fvt(native, fvt_data["critical_point"])


def tangent_plot_data(tangent):
    """Return linear-axis 2D data and the mathematically extrapolated zero-volume intercept."""
    return {
        "volume": tangent.line_v,
        "helmholtz": tangent.line_f_molar,
        "intercept_volume": 0.0,
        "intercept": tangent.mu,
    }


def homogeneous_fvt_isotherm(tangent, fvt_data, n_points=120):
    """Sample actual HEOS f(v) branches outside the equilibrium two-phase interval."""
    if isinstance(n_points, bool) or int(n_points) != n_points or n_points < 4:
        raise ValueError("n_points muss eine ganze Zahl >= 4 sein.")
    molar_mass = float(fvt_data["meta"]["molar_mass"])
    v_min = float(fvt_data["meta"]["v_min"])
    v_max = float(fvt_data["meta"]["v_max"])
    if not np.isfinite([molar_mass, v_min, v_max]).all() or molar_mass <= 0:
        raise ValueError("f-v-T-Daten enthalten keine gueltige molare Masse oder Volumengrenzen.")
    if not v_min < tangent.v_l < tangent.v_v < v_max:
        raise ValueError("Das f-v-T-Volumenraster schliesst die Sättigungsvolumina nicht ein.")

    from coolprop_heos_common import _new_state
    from coolprop_heos_fvt import _f_at_tv

    state = _new_state(tangent.fluid)
    samples_per_branch = max(2, int(n_points) // 2)
    volumes = (
        np.geomspace(v_min, tangent.v_l * (1.0 - 1e-8), samples_per_branch),
        np.geomspace(tangent.v_v * (1.0 + 1e-8), v_max, samples_per_branch),
    )
    branches = []
    for branch_volumes in volumes:
        energies = []
        valid_volumes = []

        def finish_segment():
            if len(valid_volumes) >= 2:
                branches.append((np.asarray(valid_volumes), np.asarray(energies)))
            valid_volumes.clear()
            energies.clear()

        for volume in branch_volumes:
            try:
                energy = _f_at_tv(state, tangent.temperature, float(volume))[2] * molar_mass
            except ValueError:
                finish_segment()
                continue
            if np.isfinite(energy):
                valid_volumes.append(float(volume))
                energies.append(float(energy))
            else:
                finish_segment()
        finish_segment()
    if not branches:
        raise ValueError("Keine zusammenhängenden homogenen HEOS-Werte verfügbar.")
    return branches


def calculate_fvt_tangent(stoff, T, n_points=100):
    """Calculate saturation endpoints and the physical linear f(v) common tangent."""
    if isinstance(n_points, bool) or int(n_points) != n_points or n_points < 2:
        raise ValueError("n_points muss eine ganze Zahl >= 2 sein.")
    try:
        temperature = float(T)
        state = AbstractState("HEOS", str(stoff))
        critical_temperature = float(state.T_critical())
        minimum_temperature = float(state.Tmin())
    except (ValueError, RuntimeError) as exc:
        raise ValueError(f"HEOS-Stoff oder Temperatur nicht auswertbar: {exc}") from exc

    if not np.isfinite(temperature) or not minimum_temperature < temperature < critical_temperature:
        raise ValueError(f"T={temperature:g} K liegt nicht im unterkritischen Sättigungsbereich "
                         f"({minimum_temperature:g} K < T < {critical_temperature:g} K).")

    phases = []
    try:
        for quality in (0.0, 1.0):
            state.update(CP.QT_INPUTS, quality, temperature)
            phases.append({
                "p": float(state.p()),
                "v": 1.0 / float(state.rhomolar()),
                "f_mass": float(state.helmholtzmass()),
                "g_mass": float(state.gibbsmass()),
                "molar_mass": float(state.molar_mass()),
            })
    except (ValueError, RuntimeError, ZeroDivisionError) as exc:
        raise ValueError(f"HEOS-Sättigungszustand bei T={temperature:g} K fehlgeschlagen: {exc}") from exc

    liquid, vapor = phases
    molar_mass = liquid["molar_mass"]
    values = np.array([liquid["p"], liquid["v"], vapor["v"], liquid["f_mass"], vapor["f_mass"],
                       liquid["g_mass"], vapor["g_mass"], molar_mass], dtype=float)
    if not np.isfinite(values).all() or molar_mass <= 0 or liquid["p"] <= 0 or not liquid["v"] < vapor["v"]:
        raise ValueError("HEOS lieferte unplausible oder nicht endliche Sättigungswerte.")
    if not np.isclose(liquid["p"], vapor["p"], rtol=1e-8, atol=1e-6):
        raise ValueError("HEOS-Flüssigkeits- und Dampfdrücke stimmen nicht überein.")

    f_l = liquid["f_mass"] * molar_mass
    f_v = vapor["f_mass"] * molar_mass
    g_l = liquid["g_mass"] * molar_mass
    g_v = vapor["g_mass"] * molar_mass
    mu = f_l + liquid["p"] * liquid["v"]
    mu_heos = g_l
    potential_tolerance = max(1e-4, abs(mu_heos) * 1e-8, 8.31446261815324 * temperature * 1e-7)
    chemical_potential_delta = max(abs(mu - g_l), abs(mu - g_v), abs(g_l - g_v))
    if chemical_potential_delta > potential_tolerance:
        raise ValueError(f"HEOS-Gleichgewichtspotentiale weichen um {chemical_potential_delta:.6g} J/mol ab "
                         f"(Toleranz {potential_tolerance:.6g} J/mol).")

    tangent_pressure = -(f_v - f_l) / (vapor["v"] - liquid["v"])
    pressure_relative_error = abs(tangent_pressure - liquid["p"]) / liquid["p"]
    if not np.isfinite(pressure_relative_error) or pressure_relative_error > 1e-8:
        raise ValueError(f"Tangentensteigung widerspricht dem HEOS-Sättigungsdruck "
                         f"(relativer Fehler {pressure_relative_error:.6g}).")

    line_v = np.geomspace(liquid["v"], vapor["v"], int(n_points))
    quality = (line_v - liquid["v"]) / (vapor["v"] - liquid["v"])
    return FvtTangent(
        fluid=str(stoff), temperature=temperature, p_sat=liquid["p"], v_l=liquid["v"], v_v=vapor["v"],
        f_l_mass=liquid["f_mass"], f_v_mass=vapor["f_mass"], f_l_molar=f_l, f_v_molar=f_v,
        g_l_mass=liquid["g_mass"], g_v_mass=vapor["g_mass"], g_l_molar=g_l, g_v_molar=g_v,
        mu=mu, mu_heos=mu_heos, line_v=line_v,
        line_f_molar=(1.0 - quality) * f_l + quality * f_v, tangent_pressure=tangent_pressure,
        pressure_relative_error=pressure_relative_error, chemical_potential_delta=chemical_potential_delta)
