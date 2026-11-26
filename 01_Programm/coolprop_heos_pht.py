
"""
Datei: coolprop_heos_pht.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Berechnung dreidimensionaler p-h-T-Zustandsflächen
auf Grundlage des CoolProp-HEOS-Backends.

Aufgaben:
- Berechnung spezifischer Enthalpien und kritischer Zustände
- Erstellung stoffabhängiger Temperatur- und Druckraster
- Berechnung von Flüssigkeits-, Dampf- und Hochtemperaturflächen
- Berechnung von Binodalen und Zweiphasenflächen
- Erzeugung von Isothermen, Isobaren und Isochoren
- Validierung der thermodynamischen Zustandsdaten

Abhängigkeiten:
- coolprop_heos_common.py: Gemeinsame HEOS-Funktionen
- isoline_levels.py: Isoliniengenerierung
- limits.py: Maximale Rasterpunktzahl
- CoolProp: Enthalpie- und Phasengleichgewichtsberechnung

Verwendet von:
- main.py: Berechnung des p-h-T-Diagramms
- coolprop_heos_pst.py / coolprop_heos_hsp.py
- coolprop_heos_uvt.py / coolprop_heos_fvt.py
- coolprop_heos_gpt.py
"""


from __future__ import annotations
from coolprop_heos_common import _new_state, fluid_limits, pt_update_robust
from limits import MAX_POINTS

import CoolProp.CoolProp as CP
import numpy as np
import logging


LOGGER = logging.getLogger("thermosurface.coolprop_heos_pht")


# 1. Kritischer Punkt und einzelne HEOS-Zustaende

def critical_point_pht(stoff: str) -> tuple[float, float, float]:
    """(Tc [K], pc [Pa], hc [J/kg]) aus dem kritischen HEOS-Zustand."""
    state = _new_state(stoff)
    return _critical_from_state(state)


def _critical_from_state(state) -> tuple[float, float, float]:
    Tc = float(state.T_critical())
    pc = float(state.p_critical())
    rho_c = float(state.rhomolar_critical())

    if not (np.isfinite([Tc, pc, rho_c]).all() and min(Tc, pc, rho_c) > 0):
        raise ValueError("HEOS liefert ungueltige kritische Stoffdaten.")
    
    # Dichte und Temperatur definieren den kritischen Zustand eindeutig.
    state.update(CP.DmolarT_INPUTS, rho_c, Tc)
    hc = float(state.hmass())

    if not np.isfinite(hc):
        raise ValueError("Kritische spezifische Enthalpie nicht berechenbar.")
    
    return Tc, pc, hc


def h_eos(T: float, p: float, stoff: str) -> float:
    """Spezifische Enthalpie [J/kg] aus (T,p); saettigungsnahe PT undefiniert."""
    if not (np.isfinite([T, p]).all() and T > 0 and p > 0):
        return float("nan")
    
    state = _new_state(stoff)

    try:
        pt_update_robust(state, T, p)
        h = float(state.hmass())
        return h if np.isfinite(h) else float("nan")
    
    except (ValueError, RuntimeError):
        return float("nan")


def _pt_point(state, T: float, p: float, phase=None) -> tuple[float, float, float]:
    """Einphasiger PT-Flash; optionale stabile Phasenvorgabe."""
    try:
        pt_update_robust(state, T, p, phase)
        h = float(state.hmass())

    except (ValueError, RuntimeError) as exc:
        raise ValueError(f"HEOS-PT-Flash fehlgeschlagen: T={T:.9g} K, " f"p={p:.9g} Pa: {exc}") from exc

    if not np.isfinite(h):
        raise ValueError(f"HEOS lieferte ungueltiges h bei T={T}, p={p}.")
    
    return float(T), h, float(p)


# 2. Saettigungspunkte und Binodalen

def _saturation_from_state(state, T: float):
    """Gibt ((T,h_l,p_sat),(T,h_g,p_sat)) zurueck; nur Ttr < T < Tc."""
    if not np.isfinite(T) or not (state.Ttriple() < T < state.T_critical()):
        return None
    
    try:
        state.update(CP.QT_INPUTS, 0.0, float(T))
        p_l, h_l = float(state.p()), float(state.hmass())
        state.update(CP.QT_INPUTS, 1.0, float(T))
        p_g, h_g = float(state.p()), float(state.hmass())

    except (ValueError, RuntimeError) as exc:
        LOGGER.debug("Saettigungsberechnung bei T=%s K: %s", T, exc)

        return None
    
    if not (np.isfinite([p_l, p_g, h_l, h_g]).all() and p_l > 0 and h_l < h_g and np.isclose(p_l, p_g, rtol=1e-8, atol=1e-3)):
        return None
    # Die beiden Zweiphasenrandpunkte benutzen EXAKT denselben Druckwert.
    return (float(T), h_l, p_l), (float(T), h_g, p_l)


def calc_saturation_point(T: float, stoff: str):
    """(Siedepunkt, Taupunkt), je (T,h,p), oder None ausserhalb VLE."""
    return _saturation_from_state(_new_state(stoff), float(T))


def multi_binodale(T_start_sat: float | None, T_stop_sat: float | None, T_step_sat: float | None, stoff: str, T_values=None, include_critical: bool = True):
    """Zwei Arrays (n,3); bei T_values werden Start/Stop/Step ignoriert."""
    state = _new_state(stoff)
    Tc, pc, hc = _critical_from_state(state)

    if T_values is None:
        if not (np.isfinite([T_start_sat, T_stop_sat, T_step_sat]).all() and T_step_sat > 0 and T_start_sat < T_stop_sat):
            raise ValueError("Ungueltiges Saettigungstemperaturraster.")
        
        n_est = int(np.ceil((T_stop_sat - T_start_sat) / T_step_sat))

        if n_est > MAX_POINTS:
            raise ValueError("Saettigungstemperaturraster zu gross.")
        
        temperatures = np.arange(T_start_sat, T_stop_sat, T_step_sat)
    else:
        temperatures = np.asarray(T_values, dtype=float)

    if (temperatures.ndim != 1 or len(temperatures) < 1
            or not np.isfinite(temperatures).all()
            or np.any(np.diff(temperatures) <= 0)
            or not np.all((temperatures > state.Ttriple()) & (temperatures < Tc))):
        raise ValueError("T_values muessen streng steigend zwischen Ttriple und Tc liegen.")
    
    if len(temperatures) + int(include_critical) > MAX_POINTS:
        raise ValueError("Binodalen ueberschreiten das Punktlimit.")
    
    liquid, vapor = [], []

    for T in temperatures:
        sat = _saturation_from_state(state, float(T))
        if sat is None:
            raise ValueError(f"HEOS-Saettigung bei T={T:.9g} K nicht auswertbar.")
        
        l, g = sat
        liquid.append(l)
        vapor.append(g)

    if include_critical:
        liquid.append((Tc, hc, pc))
        vapor.append((Tc, hc, pc))

    return (np.asarray(liquid, dtype=np.float64).reshape(-1, 3),
            np.asarray(vapor, dtype=np.float64).reshape(-1, 3))


def calc_two_phase_surface(binodal_liquid, binodal_vapor, n_quality: int = 40):
    """(n_T,n_Q,3), Q(n_Q,); h=(1-Q)*h_l+Q*h_g, p=p_sat(T)."""
    l = np.asarray(binodal_liquid, dtype=np.float64)
    g = np.asarray(binodal_vapor, dtype=np.float64)

    if (l.ndim != 2 or l.shape[1:] != (3,) or l.shape != g.shape or len(l) < 2 or not np.isfinite(l).all() or not np.isfinite(g).all()):
        raise ValueError("Binodalen: gleiche endliche Arrays (n_T>=2, 3) erforderlich.")
    
    _check_count(n_quality, "n_quality")

    if len(l) * n_quality > MAX_POINTS:
        raise ValueError("Zweiphasenoberflaeche ueberschreitet MAX_POINTS.")
    
    if (np.any(np.diff(l[:, 0]) <= 0) or np.any(l[:, 0] <= 0)
            or np.any(l[:, 2] <= 0)
            or not np.array_equal(l[:, 0], g[:, 0])
            or not np.array_equal(l[:, 2], g[:, 2])
            or np.any(l[:, 1] > g[:, 1])):
        raise ValueError("Binodalen sind thermodynamisch/strukturell inkonsistent.")
    
    # Gleiche Enthalpien nur im letzten gemeinsamen kritischen Endpunkt zulaessig.
    if np.any(l[:-1, 1] >= g[:-1, 1]):
        raise ValueError("Unterkritisch muss h_l < h_g gelten.")
    
    Q = np.linspace(0.0, 1.0, int(n_quality))
    T = np.broadcast_to(l[:, None, 0], (len(l), len(Q)))
    p = np.broadcast_to(l[:, None, 2], (len(l), len(Q)))
    h = (1.0 - Q[None, :]) * l[:, None, 1] + Q[None, :] * g[:, None, 1]
    # Die Grenzspalten werden bitgenau auf die Originalbinodalen gesetzt.
    h[:, 0] = l[:, 1]
    h[:, -1] = g[:, 1]

    if l[-1, 1] == g[-1, 1]:
        h[-1, :] = l[-1, 1]  # Kritische Zeile als ein einziger Punkt

    return np.stack((T, h, p), axis=-1), Q


# 3. Stoffabhaengige Berechnungsraster

def _check_count(value: int, name: str) -> None:
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or not 2 <= value <= MAX_POINTS):
        raise ValueError(f"{name} muss ganzzahlig zwischen 2 und {MAX_POINTS} sein.")


def _geom(a: float, b: float, n: int) -> np.ndarray:
    if not (0 < a < b):
        raise ValueError("Druckintervall muss 0 < p_start < p_ende erfuellen.")
    
    values = np.geomspace(a, b, n)
    values[0], values[-1] = a, b

    return values


def generate_pht_grid(stoff: str = "Water", *, T_min=None, T_max=None, p_min=None, p_max=None, n_T_sub: int = 45, n_T_super: int = 25, n_p: int = 41, n_quality: int = 41, critical_refinement: float = 1.5):
    """Liefert (raster, metadaten); bestehende p-v-T-config bleibt unberuehrt."""
    for name, value in (("n_T_sub", n_T_sub), ("n_T_super", n_T_super), ("n_p", n_p), ("n_quality", n_quality)):
        _check_count(value, name)

    if not np.isfinite(critical_refinement) or critical_refinement < 1:
        raise ValueError("critical_refinement muss endlich und >= 1 sein.")
    
    state = _new_state(stoff)
    Tc, pc, hc = _critical_from_state(state)
    Ttr, Tmin_backend, Tmax_backend = fluid_limits(stoff)
    Tlo = (max(0.75 * Tc, 1.005 * Ttr, 1.001 * Tmin_backend) if T_min is None else float(T_min))
    Thi = (min(1.20 * Tc, 0.999 * Tmax_backend) if T_max is None else float(T_max))

    if not (np.isfinite([Ttr, Tmin_backend, Tmax_backend, Tlo, Thi]).all()
            and max(Ttr, Tmin_backend) < Tlo < Tc < Thi <= Tmax_backend):
        raise ValueError("T-Bereich muss oberhalb des Tripelpunkts und beidseits Tc liegen.")
    
    sat0 = _saturation_from_state(state, Tlo)

    if sat0 is None:
        raise ValueError("Keine Saettigungswerte am unteren Temperaturrand verfuegbar.")
    
    psat_min = sat0[0][2]
    plo = min(0.02 * pc, 0.25 * psat_min) if p_min is None else float(p_min)
    phi = 2.0 * pc if p_max is None else float(p_max)

    if (not np.isfinite([plo, phi]).all() or not (0 < plo < psat_min < pc < phi)):
        raise ValueError("Druckbereich muss 0 < p_min < psat(T_min) < pc < p_max erfuellen.")
    
    # Gegebenenfalls Datenbankgrenze fuer den hoechsten Druck beachten.
    if hasattr(state, "pmax") and phi > float(state.pmax()):
        raise ValueError("p_max liegt oberhalb der HEOS-Modellgrenze.")
    
    n_total = (2 * (n_T_sub + 1) * n_p + (n_T_super + 1) * (2 * n_p - 1) + (n_T_sub + 1) * (n_quality + 2))

    if n_total > MAX_POINTS:
        raise ValueError(f"Gesamt: {n_total:,} Punkte > MAX_POINTS={MAX_POINTS:,}.")
    
    # T-Stuetzstellen zur kritischen Temperatur verdichten, Tc selbst separat.
    f = np.arange(n_T_sub, dtype=float) / n_T_sub
    T_sub = Tc - (Tc - Tlo) * np.power(1.0 - f, critical_refinement)
    T_sub[0] = Tlo
    T_high = np.linspace(Tc, Thi, n_T_super + 1)

    if (np.any(np.diff(T_sub) <= 0) or not T_sub[-1] < Tc or np.any(np.diff(T_high) <= 0)):
        raise ValueError("Temperaturraster wegen Zahlenaufloesung entartet.")
    
    raster = {"T_sub": T_sub, 
              "T_high": T_high, 
              "p_min": plo, 
              "p_max": phi, 
              "n_p": int(n_p), 
              "n_quality": int(n_quality)}
    
    meta = {"stoff": stoff, 
            "backend": "HEOS", 
            "Tc": Tc, 
            "pc": pc,
            "hc": hc, 
            "Ttr": Ttr, 
            "Tmin_HEOS": Tmin_backend,
            "Tmax_HEOS": Tmax_backend, 
            "n_points_max": MAX_POINTS}
    
    return raster, meta


# 4. Randkonforme einphasige Flaechen (T,p)->h

def _single_phase_surfaces(stoff: str, raster: dict, liquid: np.ndarray, vapor: np.ndarray):
    T_sub, T_high = raster["T_sub"], raster["T_high"]
    p_min, p_max, n_p = raster["p_min"], raster["p_max"], raster["n_p"]
    n_sub = len(T_sub)
    Tc, hc, pc = liquid[-1]
    l_surf = np.empty((n_sub + 1, n_p, 3), dtype=np.float64)
    g_surf = np.empty_like(l_surf)
    state = _new_state(stoff)

    for i, T in enumerate(T_sub):
        psat = float(liquid[i, 2])
        p_g = _geom(p_min, psat, n_p)
        p_l = _geom(psat, p_max, n_p)

        for j, p in enumerate(p_g[:-1]):
            g_surf[i, j] = _pt_point(state, T, float(p), CP.iphase_gas)

        for j, p in enumerate(p_l[1:], start=1):
            l_surf[i, j] = _pt_point(state, T, float(p), CP.iphase_liquid)

        # Exakt identische Randknoten, nicht PT-Flash bei p=psat.
        g_surf[i, -1] = vapor[i]
        l_surf[i, 0] = liquid[i]
    # Kritische T-Zeile: zwei an pc/hc zusammenlaufende Druckzweige.
    p_g_crit = _geom(p_min, pc, n_p)
    p_l_crit = _geom(pc, p_max, n_p)

    for j, p in enumerate(p_g_crit[:-1]):
        g_surf[-1, j] = _pt_point(state, Tc, float(p))

    for j, p in enumerate(p_l_crit[1:], start=1):
        l_surf[-1, j] = _pt_point(state, Tc, float(p))

    g_surf[-1, -1] = (Tc, hc, pc)
    l_surf[-1, 0] = (Tc, hc, pc)
    # Das Hochtemperaturnetz besitzt dieselbe kritische T-Zeile.
    high = np.empty((len(T_high), 2 * n_p - 1, 3), dtype=np.float64)
    high[0] = np.concatenate((g_surf[-1], l_surf[-1, 1:]), axis=0)
    pressures_high = high[0, :, 2].copy()

    for i, T in enumerate(T_high[1:], start=1):
        for j, p in enumerate(pressures_high):
            high[i, j] = _pt_point(state, float(T), float(p))

    return l_surf, g_surf, high


# 5. Gesamte p-h-T-Berechnung

def berechne_pht_daten(stoff: str = "Water", *, T_min=None, T_max=None, p_min=None, p_max=None, n_T_sub: int = 45, n_T_super: int = 25, n_p: int = 41, n_quality: int = 41, critical_refinement: float = 1.5, include_isolines: bool = False, isolines=None) -> dict:
    """Berechnet Teilflaechen und Phasengrenzen im Format (T,h,p).

    Keys: surface_liquid, surface_vapor, surface_high_temperature,
    binodal_liquid, binodal_vapor, two_phase_surface, Q_values,
    critical_point, raster, meta. Keine PyVista-Abhaengigkeit.
    Mit include_isolines zusaetzlich isothermen/isobaren/isochoren (geometrisch verteilte Levels, TS-VP03) und
    isoline_report; `isolines`: start/stop/factor je Familie (None = Standard).
    """
    raster, meta = generate_pht_grid(stoff, T_min=T_min, T_max=T_max, p_min=p_min, p_max=p_max, n_T_sub=n_T_sub, n_T_super=n_T_super, n_p=n_p, n_quality=n_quality, critical_refinement=critical_refinement)
    T_sub = raster["T_sub"]
    bin_l, bin_g = multi_binodale(None, None, None, stoff, T_values=T_sub, include_critical=True)
    two_phase, Q = calc_two_phase_surface(bin_l, bin_g, n_quality=n_quality)
    surf_l, surf_g, surf_high = _single_phase_surfaces(stoff, raster, bin_l, bin_g)
    daten = {
        "surface_liquid": surf_l,
        "surface_vapor": surf_g,
        "surface_high_temperature": surf_high,
        "binodal_liquid": bin_l,
        "binodal_vapor": bin_g,
        "two_phase_surface": two_phase,
        "Q_values": Q,
        "critical_point": (meta["Tc"], meta["pc"], meta["hc"]),
        "raster": raster,
        "meta": meta,}
    
    validate_pht_data(daten)

    if include_isolines:
        (daten["isothermen"], daten["isobaren"], daten["isochoren"],
         daten["isoline_report"]) = _isolinien(stoff, raster, isolines)

    return daten


# 6. Isolinien (geometrisch verteilte Levels, TS-VP03)

def _isolinien(stoff: str, raster: dict, isolines=None):
    """Kurven (T,h,p) entlang der Levels; Rueckgabe (isothermen, isobaren, isochoren, bericht)."""
    from isoline_levels import resolve_levels, build_isolines, V_DOMAIN_FACTORS
    T_min, T_max = float(raster["T_sub"][0]), float(raster["T_high"][-1])
    p_min, p_max = float(raster["p_min"]), float(raster["p_max"])
    n_Q = int(raster["n_quality"])
    vc = 1.0 / float(_new_state(stoff).rhomolar_critical())
    levels, report = resolve_levels(isolines, {
        "isothermen": (T_min, T_max), "isobaren": (p_min, p_max),
        "isochoren": (V_DOMAIN_FACTORS[0] * vc, V_DOMAIN_FACTORS[1] * vc)})

    def isochore(v):
        curve = calc_isochore(v, stoff, T_min, T_max, n_T=121)
        outside = ~((curve[:, 2] >= p_min) & (curve[:, 2] <= p_max))   # wie p-s-T: nur im Druckfenster der Flaeche
        curve[outside] = np.column_stack((curve[outside, 0], np.full((int(outside.sum()), 2), np.nan)))
        return curve

    curves = build_isolines(levels, report, {
        "isothermen": lambda T: calc_isotherme(T, stoff, p_min, p_max, n_p=41, n_quality=n_Q),
        "isobaren": lambda p: calc_isobare(p, stoff, T_min, T_max, n_T=81, n_quality=n_Q),
        "isochoren": isochore}, log=LOGGER.warning)

    return curves["isothermen"], curves["isobaren"], curves["isochoren"], report


# Weitere Isolinien-Bausteine (einzelne Kurven)

def calc_isotherme(T: float, stoff: str, p_min: float, p_max: float, n_p: int = 41, n_quality: int = 41) -> np.ndarray:
    """Isotherme mit ggf. Q-Abschnitt: (n,3) in T,h,p-Reihenfolge."""
    _check_count(n_p, "n_p")
    _check_count(n_quality, "n_quality")

    state = _new_state(stoff)
    Tc = state.T_critical()

    if not (np.isfinite([T, p_min, p_max]).all() and state.Tmin() <= T <= state.Tmax() and 0 < p_min < p_max):
        raise ValueError("Ungueltiger Isothermenbereich.")
    
    sat = _saturation_from_state(state, float(T)) if T < Tc else None

    if sat is None or not p_min < sat[0][2] < p_max:
        return np.asarray([_pt_point(state, T, p) for p in _geom(p_min, p_max, n_p)])
    
    l, g = sat
    pg = _geom(p_min, g[2], n_p)
    pl = _geom(l[2], p_max, n_p)
    q = np.linspace(0.0, 1.0, n_quality)
    vapor = np.asarray([_pt_point(state, T, p) for p in pg[:-1]] + [g])
    mixture = np.stack((np.full(n_quality, T), (1 - q) * l[1] + q * g[1], np.full(n_quality, l[2])), axis=-1)
    liquid = np.asarray([l] + [_pt_point(state, T, p) for p in pl[1:]])
    # Druckpfad: Dampf -> Nassdampf (Q absteigend) -> Fluessigkeit.
    return np.concatenate((vapor, mixture[-2:0:-1], liquid), axis=0)


def calc_isobare(p: float, stoff: str, T_min: float, T_max: float, n_T: int = 81, n_quality: int = 41) -> np.ndarray:
    """Isobare mit exaktem Q-Abschnitt bei Tsat, sortiert nach T und h."""
    _check_count(n_T, "n_T")
    _check_count(n_quality, "n_quality")
    state = _new_state(stoff)

    if not (np.isfinite([p, T_min, T_max]).all() and p > 0 and state.Tmin() <= T_min < T_max <= state.Tmax()):
        raise ValueError("Ungueltiger Isobarenbereich.")
    
    sat = None
    
    if p < state.p_critical():
        try:
            state.update(CP.PQ_INPUTS, float(p), 0.0)
            Tsat = float(state.T())
            sat = _saturation_from_state(state, Tsat)
    
        except (ValueError, RuntimeError):
            sat = None

    T_values = np.linspace(T_min, T_max, n_T)
    
    points = []
    
    for T in T_values:
        if sat is not None and np.isclose(T, sat[0][0], rtol=0, atol=1e-10):
            continue
        points.append(_pt_point(state, float(T), float(p)))

    if sat is not None and T_min <= sat[0][0] <= T_max:
        l, g = sat
        for q in np.linspace(0.0, 1.0, n_quality):
            points.append((l[0], (1 - q) * l[1] + q * g[1], p))

    points.sort(key=lambda row: (row[0], row[1]))

    return np.asarray(points, dtype=np.float64).reshape(-1, 3)


def calc_isochore(v_molar: float, stoff: str, T_min: float, T_max: float, n_T: int = 81) -> np.ndarray:
    """Isochore (molares v [m3/mol]) via DmolarT; Ausgabe (T,h,p)."""
    _check_count(n_T, "n_T")
    state = _new_state(stoff)
    
    if not (np.isfinite([v_molar, T_min, T_max]).all() and v_molar > 0 and state.Tmin() <= T_min < T_max <= state.Tmax()):
        raise ValueError("Ungueltiger Isochorenbereich.")
    
    result = np.empty((n_T, 3), dtype=np.float64)
    
    for i, T in enumerate(np.linspace(T_min, T_max, n_T)):
        try:
            state.update(CP.DmolarT_INPUTS, 1.0 / float(v_molar), float(T))
            h, p = float(state.hmass()), float(state.p())
    
            if not np.isfinite([h, p]).all() or p <= 0:
    
                raise ValueError("ungueltiger Zustand")
    
            result[i] = (T, h, p)
    
        except (ValueError, RuntimeError) as exc:
            LOGGER.debug("Isochore T=%s nicht auswertbar: %s", T, exc)
    
            result[i] = (T, np.nan, np.nan)
    
    return result


# 7. Reine numerische Konsistenzpruefung (ohne HEOS-Neuberechnung)

def validate_pht_data(daten: dict) -> None:
    """Dimensionen und exakte gemeinsame Phasengrenzen validieren."""
    l, g = daten["binodal_liquid"], daten["binodal_vapor"]
    mix, Q = daten["two_phase_surface"], daten["Q_values"]
    sl, sg, sh = (daten["surface_liquid"], daten["surface_vapor"], daten["surface_high_temperature"])
    for name, values in (("Binodale liquid", l), ("Binodale vapor", g), ("Zweiphasenflaeche", mix), ("Liquid", sl), ("Vapor", sg), ("High-T", sh)):
        
        if (values.ndim < 2 or values.shape[-1] != 3 or not np.isfinite(values).all() or np.any(values[..., 0] <= 0) or np.any(values[..., 2] <= 0)):
            raise ValueError(f"{name}: ungueltige Koordinaten/Arrayform.")
        
    if (l.shape != g.shape or sl.shape[0] != len(l) or sg.shape[0] != len(g) or mix.shape != (len(l), len(Q), 3) or sh.shape[1] != sl.shape[1] + sg.shape[1] - 1 or not np.array_equal(Q[[0, -1]], [0.0, 1.0])):
        raise ValueError("Ungueltige Dimensionen der p-h-T-Teilflaechen.")
    
    for name, actual, expected in (
        ("Fluessigkeitsrand", sl[:, 0, :], l),
        ("Dampfrand", sg[:, -1, :], g),
        ("2-Phasen-Siedelinie", mix[:, 0, :], l),
        ("2-Phasen-Taulinie", mix[:, -1, :], g),
        ("Kritische Hochtemperaturzeile", sh[0],
         np.concatenate((sg[-1], sl[-1, 1:]), axis=0)),):
        
        if not np.array_equal(actual, expected):
            raise ValueError(f"Randkontinuitaet verletzt: {name}.")
    
    Tc, pc, hc = daten["critical_point"]

    if not np.array_equal(l[-1], [Tc, hc, pc]) or not np.array_equal(g[-1], l[-1]):
        raise ValueError("Binodalen enden nicht am kritischen Punkt.")
    
    if not (np.all(np.diff(l[:, 0]) > 0) and np.all(l[:-1, 1] < g[:-1, 1])):
        raise ValueError("Phasengrenze unplausibel.")
