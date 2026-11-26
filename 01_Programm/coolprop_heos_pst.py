
"""
Datei: coolprop_heos_pst.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Berechnung dreidimensionaler p-s-T-Zustandsflächen
auf Grundlage des CoolProp-HEOS-Backends.

Aufgaben:
- Berechnung spezifischer Entropien und kritischer Zustände
- Erstellung einphasiger Flüssigkeits- und Dampfflächen
- Berechnung der Siede- und Taulinie
- Erzeugung der Gleichgewichts-Zweiphasenfläche
- Berechnung von Isothermen, Isobaren und Isochoren
- Nutzung gemeinsamer Temperatur- und Druckraster

Abhängigkeiten:
- coolprop_heos_pht.py: Raster und Zweiphasen-Hilfsfunktionen
- coolprop_heos_common.py: Gemeinsame HEOS-Funktionen
- isoline_levels.py: Isoliniengenerierung
- CoolProp: Entropie- und Phasengleichgewichtsberechnung

Verwendet von:
- main.py: Berechnung des p-s-T-Diagramms
- coolprop_heos_hsp.py: Bereitstellung der Entropiedaten
"""

from __future__ import annotations
from coolprop_heos_pht import (_check_count, _geom, calc_two_phase_surface, generate_pht_grid, validate_pht_data)
from coolprop_heos_common import _new_state, pt_update_robust

import CoolProp.CoolProp as CP
import numpy as np
import logging

LOGGER = logging.getLogger("thermosurface.coolprop_heos_pst")



def _critical_from_state(state):
    #(Tc [K], pc [Pa], sc [J/(kg K)])
    Tc = float(state.T_critical())
    pc = float(state.p_critical())
    rho_c = float(state.rhomolar_critical())

    if not np.isfinite([Tc, pc, rho_c]).all() or min(Tc, pc, rho_c) <= 0:
        raise ValueError("Ungueltige kritische HEOS-Stoffdaten.")
    
    state.update(CP.DmolarT_INPUTS, rho_c, Tc)
    sc = float(state.smass())
    
    if not np.isfinite(sc):
        raise ValueError("Kritische spezifische Entropie ist nicht endlich.")
    
    return Tc, pc, sc


def critical_point_pst(stoff: str):
    return _critical_from_state(_new_state(stoff))


def s_eos(T: float, p: float, stoff: str) -> float:
    """Spezifische Entropie fuer eindeutige (T,p)-Zustaende [J/(kg K)]."""
    if not np.isfinite([T, p]).all() or T <= 0 or p <= 0:
        return float("nan")
    try:
        state = _new_state(stoff)
        pt_update_robust(state, T, p)
        return float(state.smass())
    except (ValueError, RuntimeError):
        return float("nan")


def _pt_point(state, T: float, p: float, phase=None):
    """Aus (T,p) ein Punkt (T,s,p); auf dem Phasenrand NICHT verwenden."""
    try:
        pt_update_robust(state, T, p, phase)
        s = float(state.smass())

    except (ValueError, RuntimeError) as exc:
        raise ValueError(f"HEOS-PT-Flash fehlgeschlagen bei T={T:.9g} K, " f"p={p:.9g} Pa: {exc}") from exc
    
    if not np.isfinite(s):  # Entropie darf negativ sein.
        raise ValueError(f"HEOS lieferte ungueltiges s: T={T}, p={p}.")
    
    return (float(T), s, float(p))


def _saturation_from_state(state, T: float):
    """((T,s_l,p_sat),(T,s_g,p_sat)) oder None ausserhalb der VLE."""
    if not np.isfinite(T) or not state.Ttriple() < T < state.T_critical():
        return None
    try:
        state.update(CP.QT_INPUTS, 0.0, float(T))
        p_l, s_l = float(state.p()), float(state.smass())

        state.update(CP.QT_INPUTS, 1.0, float(T))
        p_g, s_g = float(state.p()), float(state.smass())

    except (ValueError, RuntimeError) as exc:
        LOGGER.debug("HEOS-Saettigung T=%s nicht berechenbar: %s", T, exc)
        return None
    
    if (not np.isfinite([p_l, p_g, s_l, s_g]).all() or p_l <= 0 or s_l >= s_g or not np.isclose(p_l, p_g, rtol=1e-8, atol=1e-3)):
        return None
    
    return (float(T), s_l, p_l), (float(T), s_g, p_l)


def calc_saturation_point(T: float, stoff: str):
    return _saturation_from_state(_new_state(stoff), float(T))


def multi_binodale(stoff: str, T_values, include_critical: bool = True):
    """Zwei Arrays (n_T,3) mit identischen T- und p-Werten."""
    state = _new_state(stoff)
    Tc, pc, sc = _critical_from_state(state)
    temperatures = np.asarray(T_values, dtype=float)

    if (temperatures.ndim != 1 or len(temperatures) < 2
            or not np.isfinite(temperatures).all()
            or np.any(np.diff(temperatures) <= 0)
            or np.any(temperatures <= state.Ttriple())
            or np.any(temperatures >= Tc)):
        raise ValueError("T_values muessen streng aufsteigend in (Ttriple,Tc) liegen.")
    
    liquid, vapor = [], []

    for T in temperatures:
        sat = _saturation_from_state(state, float(T))
        if sat is None:
            raise ValueError(f"Siedepunkt und Taupunkt bei T={T:.9g} K fehlen.")
        l, g = sat
        liquid.append(l)
        vapor.append(g)

    if include_critical:
        liquid.append((Tc, sc, pc))
        vapor.append((Tc, sc, pc))

    return np.asarray(liquid, dtype=np.float64), np.asarray(vapor, dtype=np.float64)


def _single_phase_surfaces(stoff: str, raster: dict, liquid: np.ndarray, vapor: np.ndarray):
    """Drei randkonforme Teilnetze: fluessig, dampffoermig, T >= Tc."""
    T_sub, T_high = raster["T_sub"], raster["T_high"]
    p_min, p_max, n_p = raster["p_min"], raster["p_max"], raster["n_p"]
    Tc, sc, pc = liquid[-1]
    n_sub = len(T_sub)
    liquid_surface = np.empty((n_sub + 1, n_p, 3), dtype=np.float64)
    vapor_surface = np.empty_like(liquid_surface)
    state = _new_state(stoff)

    for i, T in enumerate(T_sub):
        p_sat = float(liquid[i, 2])
        p_v = _geom(p_min, p_sat, n_p)
        p_l = _geom(p_sat, p_max, n_p)

        for j, p in enumerate(p_v[:-1]):
            vapor_surface[i, j] = _pt_point(state, T, float(p), CP.iphase_gas)
        
        for j, p in enumerate(p_l[1:], start=1):
            liquid_surface[i, j] = _pt_point(state, T, float(p), CP.iphase_liquid)

        vapor_surface[i, -1] = vapor[i]       # exakt derselbe Randpunkt
        liquid_surface[i, 0] = liquid[i]     # exakt derselbe Randpunkt

    p_v_crit = _geom(p_min, pc, n_p)
    p_l_crit = _geom(pc, p_max, n_p)

    for j, p in enumerate(p_v_crit[:-1]):
        vapor_surface[-1, j] = _pt_point(state, Tc, float(p))
    for j, p in enumerate(p_l_crit[1:], start=1):
        liquid_surface[-1, j] = _pt_point(state, Tc, float(p))

    vapor_surface[-1, -1] = (Tc, sc, pc)
    liquid_surface[-1, 0] = (Tc, sc, pc)

    high = np.empty((len(T_high), 2 * n_p - 1, 3), dtype=np.float64)
    high[0] = np.concatenate((vapor_surface[-1], liquid_surface[-1, 1:]), axis=0)
    pressures = high[0, :, 2].copy()

    for i, T in enumerate(T_high[1:], start=1):
        for j, p in enumerate(pressures):
            high[i, j] = _pt_point(state, float(T), float(p))

    return liquid_surface, vapor_surface, high


def calc_isotherme(T: float, stoff: str, p_min: float, p_max: float, n_p: int = 41, n_quality: int = 41):
    """Isotherme inklusive Q-Abschnitt bei subkritischem T."""
    _check_count(n_p, "n_p")
    _check_count(n_quality, "n_quality")
    state = _new_state(stoff)

    if not (np.isfinite([T, p_min, p_max]).all() and state.Tmin() <= T <= state.Tmax() and 0 < p_min < p_max):
        raise ValueError("Ungueltiger Isothermenbereich.")
    
    sat = _saturation_from_state(state, float(T)) if T < state.T_critical() else None

    if sat is None or not p_min < sat[0][2] < p_max:
        return np.array([_pt_point(state, T, p) for p in _geom(p_min, p_max, n_p)])
    
    l, g = sat
    q = np.linspace(0.0, 1.0, n_quality)
    vapor = [_pt_point(state, T, p) for p in _geom(p_min, g[2], n_p)[:-1]] + [g]
    mix = np.column_stack((np.full(n_quality, T), (1-q) * l[1] + q * g[1], np.full(n_quality, l[2])))
    liquid = [l] + [_pt_point(state, T, p) for p in _geom(l[2], p_max, n_p)[1:]]

    return np.concatenate((vapor, mix[-2:0:-1], liquid), axis=0)


def calc_isobare(p: float, stoff: str, T_min: float, T_max: float, n_T: int = 81, n_quality: int = 41):
    """Isobare mit Dampfgehalt-Pfad beim Siedetemperaturpunkt."""
    _check_count(n_T, "n_T")
    _check_count(n_quality, "n_quality")
    state = _new_state(stoff)

    if not (np.isfinite([p, T_min, T_max]).all() and p > 0 and state.Tmin() <= T_min < T_max <= state.Tmax()):
        raise ValueError("Ungueltiger Isobarenbereich.")
    
    sat = None
    if p < state.p_critical():
        try:
            state.update(CP.PQ_INPUTS, float(p), 0.0)
            sat = _saturation_from_state(state, float(state.T()))
        except (ValueError, RuntimeError):
            pass

    rows = []
    for T in np.linspace(T_min, T_max, n_T):
        if sat is not None and np.isclose(T, sat[0][0], rtol=0, atol=1e-10):
            continue

        rows.append(_pt_point(state, float(T), float(p)))
    if sat is not None and T_min <= sat[0][0] <= T_max:
        l, g = sat
        for q in np.linspace(0, 1, n_quality):
            rows.append((l[0], (1-q)*l[1] + q*g[1], p))

    rows.sort(key=lambda row: (row[0], row[1]))

    return np.asarray(rows, dtype=np.float64)


def calc_isochore(v_molar: float, stoff: str, T_min: float, T_max: float, n_T: int = 81, p_min: float = 0.0, p_max: float = float("inf")):
    """Isochore mit konstantem molarem Volumen [m3/mol], T/s/p-Ausgabe."""
    _check_count(n_T, "n_T")
    state = _new_state(stoff)

    if not (v_molar > 0 and state.Tmin() <= T_min < T_max <= state.Tmax()):
        raise ValueError("Ungueltige Isochorenparameter.")
    
    result = np.empty((n_T, 3), dtype=np.float64)

    for i, T in enumerate(np.linspace(T_min, T_max, n_T)):
        try:
            state.update(CP.DmolarT_INPUTS, 1.0 / v_molar, float(T))
            s, p = float(state.smass()), float(state.p())
            if not np.isfinite([s, p]).all() or p <= 0:
                raise ValueError("Ungueltiger Zustand")
            result[i] = (T, s, p) if p_min <= p <= p_max else (T, np.nan, np.nan)
        except (ValueError, RuntimeError):
            result[i] = (T, np.nan, np.nan)

    return result


def _isolinien(stoff: str, raster: dict, critical, isolines=None):
    """Isothermen/Isobaren/Isochoren entlang der geometrisch verteilten Levels (TS-VP03).

    `isolines`: start/stop/factor je Familie (None = Standard: Faktor 3.14, Bereich = Diagrammbereich).
    Rueckgabe: (isothermen, isobaren, isochoren, bericht). Das Flaechenraster bleibt unveraendert.
    """
    from isoline_levels import resolve_levels, build_isolines, V_DOMAIN_FACTORS
    Tc, pc, sc = critical
    T_sub, T_high = raster["T_sub"], raster["T_high"]
    T_min, T_max = float(T_sub[0]), float(T_high[-1])
    p_min, p_max = float(raster["p_min"]), float(raster["p_max"])
    n_Q = int(raster["n_quality"])
    state = _new_state(stoff)
    vc = 1.0 / float(state.rhomolar_critical())
    levels, report = resolve_levels(isolines, {
        "isothermen": (T_min, T_max), "isobaren": (p_min, p_max),
        "isochoren": (V_DOMAIN_FACTORS[0] * vc, V_DOMAIN_FACTORS[1] * vc)})
    curves = build_isolines(levels, report, {
        "isothermen": lambda T: calc_isotherme(T, stoff, p_min, p_max, n_p=41, n_quality=n_Q),
        "isobaren": lambda p: calc_isobare(p, stoff, T_min, T_max, n_T=81, n_quality=n_Q),
        "isochoren": lambda v: calc_isochore(v, stoff, T_min, T_max, n_T=121, p_min=p_min, p_max=p_max)},
        log=LOGGER.warning)

    return curves["isothermen"], curves["isobaren"], curves["isochoren"], report


def berechne_pst_daten(stoff: str = "Water", *, T_min=None, T_max=None, p_min=None, p_max=None, n_T_sub: int = 45, n_T_super: int = 25, n_p: int = 41, n_quality: int = 41, critical_refinement: float = 1.5, include_isolines: bool = True, isolines=None) -> dict:
    """Vier Flaechenteile und Phasengrenzen; Zustandsdaten stets (T,s,p)."""
    # Gemeinsames, bewaehrtes Temperatur-/Druckraster. Nur fuer Geometrie.
    raster, meta_h = generate_pht_grid(stoff, T_min=T_min, T_max=T_max, p_min=p_min, p_max=p_max, n_T_sub=n_T_sub, n_T_super=n_T_super, n_p=n_p, n_quality=n_quality, critical_refinement=critical_refinement)
    Tc, pc, sc = critical_point_pst(stoff)

    if not np.isclose([Tc, pc], [meta_h["Tc"], meta_h["pc"]], rtol=1e-12).all():
        raise ValueError("Uneinheitlicher kritischer HEOS-Zustand.")
    
    meta = dict(meta_h)
    meta.pop("hc", None)
    meta["sc"] = sc

    liquid, vapor = multi_binodale(stoff, raster["T_sub"], include_critical=True)
    two_phase, Q = calc_two_phase_surface(liquid, vapor, n_quality=n_quality)
    surf_l, surf_g, surf_high = _single_phase_surfaces(stoff, raster, liquid, vapor)

    daten = {
        "surface_liquid": surf_l,
        "surface_vapor": surf_g,
        "surface_high_temperature": surf_high,
        "binodal_liquid": liquid,
        "binodal_vapor": vapor,
        "two_phase_surface": two_phase,
        "Q_values": Q,
        "critical_point": (Tc, pc, sc),
        "raster": raster,
        "meta": meta,
    }
    # Die p-h-T-Validierung testet Geometrie/Kontinuitaet, nicht die h-Physik.
    validate_pht_data(daten)
    if include_isolines:
        (daten["isothermen"], daten["isobaren"], daten["isochoren"],
         daten["isoline_report"]) = _isolinien(stoff, raster, (Tc, pc, sc), isolines)
    else:
        daten.update(isothermen=[], isobaren=[], isochoren=[])

    return daten
