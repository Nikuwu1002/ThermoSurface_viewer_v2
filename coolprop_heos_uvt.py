
"""
Datei: coolprop_heos_uvt.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Berechnung dreidimensionaler u-v-T-Zustandsflächen
mit spezifischer innerer Energie und molarem Volumen.

Aufgaben:
- Berechnung der spezifischen inneren Energie mit HEOS
- Erzeugung von Flüssigkeits-, Dampf- und Hochtemperaturflächen
- Berechnung der Binodalen und der Zweiphasenfläche
- Erzeugung von Isothermen, Isobaren und Isochoren
- Validierung gemeinsamer Flächen- und Phasengrenzen

Abhängigkeiten:
- coolprop_heos_common.py: HEOS-Zustandserzeugung
- coolprop_heos_pvt.py: Sättigungsvolumina und Binodalen
- coolprop_heos_pht.py: Temperaturraster und Isobaren
- isoline_levels.py: Isoliniengenerierung
- limits.py: Maximale Rasterpunktzahl

Verwendet von:
- main.py: Berechnung des u-v-T-Diagramms
"""


from __future__ import annotations
from coolprop_heos_common import _new_state
from coolprop_heos_pvt import multi_binodale
from coolprop_heos_pht import generate_pht_grid, _check_count, calc_isobare
from limits import MAX_POINTS

import CoolProp.CoolProp as CP
import numpy as np
import logging


LOGGER = logging.getLogger("thermosurface.coolprop_heos_uvt")


def _u_at_tv(state, T, v):
    """Ein HEOS-Zustand (T,v_molar) -> (T,v_molar,u_mass)."""
    T, v = float(T), float(v)
    
    if not np.isfinite([T, v]).all() or T <= 0 or v <= 0:
        raise ValueError("u-v-T: T und molares Volumen muessen endlich und positiv sein.")
    
    try:
        state.update(CP.DmolarT_INPUTS, 1.0 / v, T)
        u = float(state.umass())
    
    except (ValueError, RuntimeError) as exc:
        raise ValueError(f"u-v-T: HEOS bei T={T:.8g} K, v={v:.8g} m3/mol: {exc}") from exc
    
    if not np.isfinite(u):
        raise ValueError(f"u-v-T: ungueltige innere Energie bei T={T}, v={v}.")
    
    return (T, v, u)


def _sat_uvt(stoff, T_sub, state):
    """Alte p-v-T-Saettigungsvolumina mit HEOS-spezifischen u-Werten."""
    sat_l, sat_g = multi_binodale(None, None, None, stoff, T_values=T_sub)
    liquid = np.empty((len(T_sub) + 1, 3), dtype=np.float64)
    vapor = np.empty_like(liquid)
    
    for i, T in enumerate(T_sub):
        state.update(CP.QT_INPUTS, 0.0, float(T))
        u_l = float(state.umass())
        v_l = 1.0 / float(state.rhomolar())
        state.update(CP.QT_INPUTS, 1.0, float(T))
        u_g = float(state.umass())
        v_g = 1.0 / float(state.rhomolar())
        # p-v-T nutzt ebenfalls CoolProp-Saettigungsdichten: Konsistenzpruefung.
        np.testing.assert_allclose([v_l, v_g], [sat_l[i, 1], sat_g[i, 1]], rtol=1e-8, atol=1e-12)
        liquid[i] = (T, sat_l[i, 1], u_l)
        vapor[i] = (T, sat_g[i, 1], u_g)

    Tc = float(state.T_critical())
    vc = 1.0 / float(state.rhomolar_critical())
    state.update(CP.DmolarT_INPUTS, 1.0 / vc, Tc)
    uc = float(state.umass())

    liquid[-1] = (Tc, vc, uc)
    vapor[-1] = (Tc, vc, uc)
  
    return liquid, vapor


def validate_uvt_data(daten):
    """Topologie- und Konsistenzpruefung auf den u-v-T-Datenarrays."""
    l, g = daten["binodal_liquid"], daten["binodal_vapor"]
    sl, sg = daten["surface_liquid"], daten["surface_vapor"]
    high, mix, q = (daten["surface_high_temperature"], daten["two_phase_surface"], daten["Q_values"])

    for name in ("surface_liquid", "surface_vapor", "surface_high_temperature","binodal_liquid", "binodal_vapor", "two_phase_surface"):
        arr = np.asarray(daten[name])

        if arr.ndim not in (2, 3) or arr.shape[-1] != 3 or not np.isfinite(arr).all():
            raise ValueError(f"u-v-T: ungueltiges Array {name}: {arr.shape}")
        
        if np.any(arr[..., 0] <= 0) or np.any(arr[..., 1] <= 0):
            raise ValueError(f"u-v-T: nichtpositive T/v-Koordinaten in {name}.")
        
    if (l.shape != g.shape or sl.shape[0] != len(l) or sg.shape[0] != len(g) or mix.shape != (len(l), len(q), 3) or high.shape[1] != sl.shape[1] + sg.shape[1] - 1):
        raise ValueError("u-v-T: Teilflaechen haben inkompatible Dimensionen.")
    
    for name, actual, expected in (
        ("Fluessigkeitsrand", sl[:, -1], l),
        ("Dampfrand", sg[:, 0], g),
        ("Zweiphasenrand Q=0", mix[:, 0], l),
        ("Zweiphasenrand Q=1", mix[:, -1], g),
        ("Hochtemperatur-Anschluss", high[0],
         np.concatenate((sl[-1], sg[-1, 1:]), axis=0)),):
        if not np.array_equal(actual, expected):
            raise ValueError(f"u-v-T: Randkontinuitaet verletzt: {name}.")
        
    if not np.array_equal(l[-1], g[-1]) or not np.array_equal(l[-1], np.asarray(daten["critical_point"])):
        raise ValueError("u-v-T: kein gemeinsamer kritischer Endpunkt.")
    
    if (np.any(l[:-1, 1] >= g[:-1, 1]) or not np.all(np.diff(l[:, 0]) > 0)):
        raise ValueError("u-v-T: Binodale thermodynamisch unplausibel.")
    
    if len(q) < 2 or q[0] != 0.0 or q[-1] != 1.0:
        raise ValueError("u-v-T: Dampfgehaltsraster ungueltig.")


def _saturation_uvt_at(state, T):
    """((v_l, u_l), (v_g, u_g)) bei beliebigem T < Tc oder None."""
    try:
        state.update(CP.QT_INPUTS, 0.0, float(T))
        liquid = (1.0 / float(state.rhomolar()), float(state.umass()))
        state.update(CP.QT_INPUTS, 1.0, float(T))
        vapor = (1.0 / float(state.rhomolar()), float(state.umass()))
    except (ValueError, RuntimeError):
        return None
    return liquid, vapor


def _isobare_uvt(state, stoff, p, T_lo, T_hi, n_quality):
    """u-v-T-Isobare aus der p-h-T-Isobare per h-p-Flash (Zustandspunkte, die nicht rueckgerechnet werden: NaN)."""
    hp = calc_isobare(float(p), stoff, T_lo, T_hi, n_T=61, n_quality=n_quality)
    curve = np.full(hp.shape, np.nan, dtype=np.float64)

    for i, (T, h, pressure) in enumerate(hp):
        try:
            state.update(CP.HmassP_INPUTS, float(h), float(pressure))

            if not np.isclose(float(state.T()), T, rtol=0, atol=1e-3):
                raise ValueError("HmassP-Flash ergibt abweichende Temperatur.")

            curve[i] = (float(state.T()), 1.0 / float(state.rhomolar()), float(state.umass()))

        except (ValueError, RuntimeError) as exc:
            LOGGER.debug("u-v-T-Isobarenpunkt %d: %s", i, exc)

    return curve


def _calculate_isolines(daten, stoff, state, n_quality, n_v, isolines=None):
    """Isolinien entlang geometrisch verteilter Levels (TS-VP03); das Flaechenraster bleibt unveraendert."""
    from isoline_levels import resolve_levels, build_isolines, isotherm_over_volume, P_DOMAIN_FACTORS
    Tc, vc, _ = daten["critical_point"]
    T_sub, T_high = daten["raster"]["T_sub"], daten["raster"]["T_high"]
    T_lo, T_hi = float(T_sub[0]), float(T_high[-1])
    temperatures = np.r_[T_sub, T_high]
    v_min, v_max = daten["meta"]["v_min"], daten["meta"]["v_max"]
    pc = float(state.p_critical())
    levels, report = resolve_levels(isolines, {
        "isothermen": (T_lo, T_hi), "isochoren": (float(v_min), float(v_max)),
        "isobaren": (P_DOMAIN_FACTORS[0] * pc, P_DOMAIN_FACTORS[1] * pc)})

    def isochore(v):
        curve = np.full((len(temperatures), 3), np.nan, dtype=float)
        for i, T in enumerate(temperatures):
            try:
                curve[i] = _u_at_tv(state, T, v)
            except ValueError as exc:
                LOGGER.debug("u-v-T-Isochore: %s", exc)
        return curve

    def isobare(p):
        curve = _isobare_uvt(state, stoff, p, T_lo, T_hi, n_quality)
        outside = ~((curve[:, 1] >= v_min) & (curve[:, 1] <= v_max))
        curve[outside, 1:] = np.nan
        return curve

    curves = build_isolines(levels, report, {
        "isothermen": lambda T: isotherm_over_volume(
            T, Tc, float(v_min), float(v_max), n_v, n_quality, lambda t, v, phase: _u_at_tv(state, t, v),
            lambda t: _saturation_uvt_at(state, t)),
        "isochoren": isochore, "isobaren": isobare}, log=LOGGER.warning)
    daten["isothermen"], daten["isobaren"], daten["isochoren"] = (
        curves["isothermen"], curves["isobaren"], curves["isochoren"])
    daten["isoline_report"] = report


def berechne_uvt_daten(stoff: str = "Water", *, T_min=None, T_max=None,v_min=None, v_max=None, n_T_sub: int = 45, n_T_super: int = 25, n_v: int = 41, n_quality: int = 41, critical_refinement: float = 1.5, include_isolines: bool = True, isolines=None) -> dict:
    """Vollstaendige stabile 3D-u-v-T-Daten, mit separater Zweiphasenflaeche.

    Oberflaechenkoordinaten: (T [K], v_molar [m3/mol], u_mass [J/kg]).
    Metastabile/instabile Zustandszweige und Spinodalen sind nicht enthalten.
    """
    _check_count(n_v, "n_v")
    # p-h-T stellt das bereits erprobte Tc-verdichtete Temperaturraster bereit.
    # n_p=2: es wird KEINE p-h-T-Flache und auch kein Druckraster ausgewertet.
    raster, _ = generate_pht_grid(stoff, T_min=T_min, T_max=T_max, n_T_sub=n_T_sub, n_T_super=n_T_super, n_p=2, n_quality=n_quality, critical_refinement=critical_refinement)
    T_sub, T_high = raster["T_sub"], raster["T_high"]
    state = _new_state(stoff)
    Tc, pc = float(state.T_critical()), float(state.p_critical())
    liquid, vapor = _sat_uvt(stoff, T_sub, state)
    _, vc, uc = liquid[-1]

    # Am fluessigen Rand nur moderat verdichten. Zu kleine v-Werte
    # koennen ausserhalb der HEOS-Modellgueltigkeit liegen.
    default_v_min = 0.985 * float(np.min(liquid[:-1, 1]))
    default_v_max = max(30.0 * vc, 1.2 * float(np.max(vapor[:-1, 1])))
    v_min = default_v_min if v_min is None else float(v_min)
    v_max = default_v_max if v_max is None else float(v_max)

    if (not np.isfinite([v_min, v_max]).all() or not (0 < v_min < np.min(liquid[:, 1]) and v_max > np.max(vapor[:, 1]))):
        raise ValueError("u-v-T: v_min/v_max muessen alle Binodalen einschliessen.")

    n_total = (2 * len(liquid) * n_v + len(T_high) * (2 * n_v - 1) + len(liquid) * n_quality)

    if n_total > MAX_POINTS:
        raise ValueError(f"u-v-T: {n_total:,} Punkte > MAX_POINTS={MAX_POINTS:,}.")

    sl = np.empty((len(liquid), n_v, 3), dtype=np.float64)
    sg = np.empty_like(sl)

    for i, T in enumerate(np.r_[T_sub, Tc]):
        v_l = np.geomspace(v_min, liquid[i, 1], n_v)
        v_g = np.geomspace(vapor[i, 1], v_max, n_v)

        for j, v in enumerate(v_l[:-1]):
            sl[i, j] = _u_at_tv(state, T, v)

        for j, v in enumerate(v_g[1:], start=1):
            sg[i, j] = _u_at_tv(state, T, v)
        # Exakt gleiche Binodal-/Flaechenknoten.
        sl[i, -1] = liquid[i]
        sg[i, 0] = vapor[i]

    high = np.empty((len(T_high), 2 * n_v - 1, 3), dtype=np.float64)
    high[0] = np.concatenate((sl[-1], sg[-1, 1:]), axis=0)
    # dieselben Volumenstuetzstellen bei allen T >= Tc
    for i, T in enumerate(T_high[1:], start=1):
        for j, v in enumerate(high[0, :, 1]):
            high[i, j] = _u_at_tv(state, T, v)

    Q = np.linspace(0.0, 1.0, n_quality)
    T_mix = np.broadcast_to(liquid[:, None, 0], (len(liquid), len(Q)))
    v_mix = ((1.0 - Q[None, :]) * liquid[:, None, 1] + Q[None, :] * vapor[:, None, 1])
    u_mix = ((1.0 - Q[None, :]) * liquid[:, None, 2] + Q[None, :] * vapor[:, None, 2])

    mix = np.stack((T_mix, v_mix, u_mix), axis=-1)
    
    mix[:, 0] = liquid
    mix[:, -1] = vapor
    mix[-1, :] = liquid[-1]  # kritische Zeile ist geometrisch ein Punkt

    result = {
        "surface_liquid": sl,
        "surface_vapor": sg,
        "surface_high_temperature": high,
        "binodal_liquid": liquid,
        "binodal_vapor": vapor,
        "two_phase_surface": mix,
        "Q_values": Q,
        "critical_point": (Tc, vc, uc),
        "raster": raster,
        "meta": {"stoff": stoff, 
                 "backend": "HEOS", 
                 "Tc": Tc, 
                 "pc": pc,
                 "vc": vc, 
                 "uc": uc, 
                 "v_min": v_min, 
                 "v_max": v_max},
        "isothermen": [], 
        "isobaren": [], 
        "isochoren": [],}
    
    validate_uvt_data(result)

    if include_isolines:
        _calculate_isolines(result, stoff, state, n_quality, n_v, isolines)
        
    return result
