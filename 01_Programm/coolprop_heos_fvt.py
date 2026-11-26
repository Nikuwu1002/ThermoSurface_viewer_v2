"""HEOS-Zustandsdaten fuer 3D-f-v-T: (T [K], v_molar [m3/mol], f [J/kg]).

f = u - T*s ist die massenspezifische freie Energie (Helmholtz-Energie), direkt aus CoolProp
(`helmholtzmass`, identisch zu umass - T*smass). Das Volumen steht - wie in p-v-T und u-v-T - als molares
Volumen v = 1/rho_molar [m3/mol] in den Daten; das massenspezifische Volumen ist v/M (M in kg/mol).

Thermodynamische Besonderheiten gegenueber u-v-T:

* f(T,v) ist die natuerliche Helmholtz-Potentialflaeche der Zustandsgleichung: p = -(df/dv)_T,
  s = -(df/dT)_v. Sie ist auch im metastabilen und instabilen Bereich einwertig (HEOS wertet
  alpha_r(delta, tau) aus).
* Im Gleichgewicht ersetzt die gemeinsame Tangente (Maxwell-Konstruktion) die nichtkonvexe Fortsetzung: Die
  Zweiphasenflaeche ist eine Regelflaeche mit f = (1-Q)*f_l + Q*f_g und v = (1-Q)*v_l + Q*v_g und
  df/dv = -p_sat. Anders als in g-p-T gilt hier f_l != f_g.
* Die homogene HEOS-Fortsetzung (metastabil, instabil) liegt ueber dieser Tangentenflaeche und ist eine
  eigene optionale Ebene, nie Teil der Gleichgewichtsflaeche.

Oberflaechenpunkte: (T, v, f); Arrays (..., 3), float64. Keine GUI-/Visualisierungsfunktionen.
"""

from __future__ import annotations

from coolprop_heos_common import (_new_state, critical_gf, gibbs_helmholtz_tv, homogeneous_gf_layers,
                                  pt_update_robust)
from coolprop_heos_pht import _check_count, generate_pht_grid
from limits import MAX_POINTS

import CoolProp.CoolProp as CP
import numpy as np
import logging


LOGGER = logging.getLogger("thermosurface.coolprop_heos_fvt")


def critical_point_fvt(stoff: str) -> tuple[float, float, float]:
    """(Tc [K], vc [m3/mol], fc [J/kg]) aus dem kritischen HEOS-Zustand."""
    Tc, _pc, vc, _gc, fc = critical_gf(stoff)
    return Tc, vc, fc


def _f_at_tv(state, T, v, phase=None):
    """Ein HEOS-Zustand (T, v_molar) -> (T, v_molar, f_mass)."""
    try:
        _p, _g, f = gibbs_helmholtz_tv(state, T, v, phase)
    except ValueError as exc:
        raise ValueError(f"f-v-T: {exc}") from exc

    return float(T), float(v), f


def f_eos(T: float, v: float, stoff: str) -> float:
    """Spezifische freie Energie [J/kg] aus (T, v_molar); NaN, wenn nicht berechenbar."""
    try:
        return _f_at_tv(_new_state(stoff), T, v)[2]
    except ValueError:
        return float("nan")


def _sat_fvt(stoff, T_sub, state):
    """Siede-/Taulinie (T, v_sat, f_sat) inkl. kritischem Endpunkt; Spalten fuer f aus HEOS-Saettigungszustaenden."""
    liquid = np.empty((len(T_sub) + 1, 3), dtype=np.float64)
    vapor = np.empty_like(liquid)

    for i, T in enumerate(T_sub):
        try:
            for row, q in ((liquid, 0.0), (vapor, 1.0)):
                state.update(CP.QT_INPUTS, q, float(T))
                row[i] = (float(T), 1.0 / float(state.rhomolar()), float(state.helmholtzmass()))
        except (ValueError, RuntimeError) as exc:
            raise ValueError(f"f-v-T: HEOS-Saettigung bei T={T:.9g} K nicht auswertbar: {exc}") from exc

    Tc, vc, fc = critical_point_fvt(stoff)
    liquid[-1] = (Tc, vc, fc)
    vapor[-1] = (Tc, vc, fc)

    return liquid, vapor


def validate_fvt_data(daten):
    """Topologie- und Konsistenzpruefung auf den f-v-T-Datenarrays."""
    l, g = daten["binodal_liquid"], daten["binodal_vapor"]
    sl, sg = daten["surface_liquid"], daten["surface_vapor"]
    high, mix, q = daten["surface_high_temperature"], daten["two_phase_surface"], daten["Q_values"]

    for name in ("surface_liquid", "surface_vapor", "surface_high_temperature", "binodal_liquid", "binodal_vapor",
                 "two_phase_surface"):
        arr = np.asarray(daten[name])

        if arr.ndim not in (2, 3) or arr.shape[-1] != 3 or not np.isfinite(arr).all():
            raise ValueError(f"f-v-T: ungueltiges Array {name}: {arr.shape}")

        if np.any(arr[..., 0] <= 0) or np.any(arr[..., 1] <= 0):
            raise ValueError(f"f-v-T: nichtpositive T/v-Koordinaten in {name}.")

    if (l.shape != g.shape or sl.shape[0] != len(l) or sg.shape[0] != len(g) or mix.shape != (len(l), len(q), 3)
            or high.shape[1] != sl.shape[1] + sg.shape[1] - 1):
        raise ValueError("f-v-T: Teilflaechen haben inkompatible Dimensionen.")

    for name, actual, expected in (("Fluessigkeitsrand", sl[:, -1], l), ("Dampfrand", sg[:, 0], g),
                                   ("Zweiphasenrand Q=0", mix[:, 0], l), ("Zweiphasenrand Q=1", mix[:, -1], g),
                                   ("Hochtemperatur-Anschluss", high[0], np.concatenate((sl[-1], sg[-1, 1:]), axis=0))):
        if not np.array_equal(actual, expected):
            raise ValueError(f"f-v-T: Randkontinuitaet verletzt: {name}.")

    if not np.array_equal(l[-1], g[-1]) or not np.array_equal(l[-1], np.asarray(daten["critical_point"])):
        raise ValueError("f-v-T: kein gemeinsamer kritischer Endpunkt.")

    if np.any(l[:-1, 1] >= g[:-1, 1]) or not np.all(np.diff(l[:, 0]) > 0):
        raise ValueError("f-v-T: Binodale thermodynamisch unplausibel.")

    if len(q) < 2 or q[0] != 0.0 or q[-1] != 1.0:
        raise ValueError("f-v-T: Dampfgehaltsraster ungueltig.")


def _isobare(state, p, T_min, T_max, n_T, n_quality):
    """f(v) entlang einer Isobare; Fluessig-/Dampfast per PT-Flash mit Phasenvorgabe, Nassdampfstrecke linear in Q."""
    Tsat, sat = None, None

    if p < float(state.p_critical()):
        try:
            state.update(CP.PQ_INPUTS, float(p), 0.0)
            Tsat = float(state.T())
            v_l, f_l = 1.0 / float(state.rhomolar()), float(state.helmholtzmass())
            state.update(CP.PQ_INPUTS, float(p), 1.0)
            v_g, f_g = 1.0 / float(state.rhomolar()), float(state.helmholtzmass())
            sat = (v_l, f_l, v_g, f_g)
        except (ValueError, RuntimeError):
            Tsat = None

    points = []

    for T in np.linspace(T_min, T_max, n_T):
        if Tsat is not None and np.isclose(T, Tsat, rtol=1e-9, atol=0.0):
            continue

        phase = None if Tsat is None else (CP.iphase_liquid if T < Tsat else CP.iphase_gas)
        pt_update_robust(state, float(T), float(p), phase)
        points.append((float(T), 1.0 / float(state.rhomolar()), float(state.helmholtzmass())))

    if sat is not None and T_min <= Tsat <= T_max:
        v_l, f_l, v_g, f_g = sat

        for q in np.linspace(0.0, 1.0, n_quality):
            points.append((Tsat, (1.0 - q) * v_l + q * v_g, (1.0 - q) * f_l + q * f_g))

    # nach T, im Nassdampfstueck nach v sortiert (v steigt mit Q)
    return np.asarray(sorted(points), dtype=np.float64)


def _saturation_fvt_at(state, T):
    """((v_l, f_l), (v_g, f_g)) bei beliebigem T < Tc oder None."""
    try:
        state.update(CP.QT_INPUTS, 0.0, float(T))
        liquid = (1.0 / float(state.rhomolar()), float(state.helmholtzmass()))
        state.update(CP.QT_INPUTS, 1.0, float(T))
        vapor = (1.0 / float(state.rhomolar()), float(state.helmholtzmass()))
    except (ValueError, RuntimeError):
        return None
    return liquid, vapor


def _calculate_isolines(daten, stoff, state, n_quality, n_v, isolines=None):
    """Isolinien entlang geometrisch verteilter Levels (TS-VP03); das Flaechenraster bleibt unveraendert."""
    from isoline_levels import resolve_levels, build_isolines, isotherm_over_volume, P_DOMAIN_FACTORS
    Tc, vc, _fc = daten["critical_point"]
    T_sub, T_high = daten["raster"]["T_sub"], daten["raster"]["T_high"]
    T_lo, T_hi = float(T_sub[0]), float(T_high[-1])
    temperatures = np.r_[T_sub, T_high]
    v_min, v_max = daten["meta"]["v_min"], daten["meta"]["v_max"]
    pc = float(state.p_critical())
    levels, report = resolve_levels(isolines, {
        "isothermen": (T_lo, T_hi), "isochoren": (float(v_min), float(v_max)),
        "isobaren": (P_DOMAIN_FACTORS[0] * pc, P_DOMAIN_FACTORS[1] * pc)})

    def isochore(v):
        curve = np.full((len(temperatures), 3), np.nan, dtype=np.float64)

        for i, T in enumerate(temperatures):
            try:
                curve[i] = _f_at_tv(state, T, v)
            except ValueError as exc:
                LOGGER.debug("f-v-T-Isochore: %s", exc)

        return curve

    def isobare(p):
        curve = _isobare(state, p, T_lo, T_hi, 61, n_quality)
        keep = (curve[:, 1] >= v_min) & (curve[:, 1] <= v_max)
        return np.where(keep[:, None], curve, np.nan)

    curves = build_isolines(levels, report, {
        "isothermen": lambda T: isotherm_over_volume(
            T, Tc, float(v_min), float(v_max), n_v, n_quality,
            lambda t, v, phase: _f_at_tv(state, t, v, {"liquid": CP.iphase_liquid, "gas": CP.iphase_gas}.get(phase)),
            lambda t: _saturation_fvt_at(state, t)),
        "isochoren": isochore, "isobaren": isobare}, log=LOGGER.warning)
    daten["isothermen"], daten["isobaren"], daten["isochoren"] = (
        curves["isothermen"], curves["isobaren"], curves["isochoren"])
    daten["isoline_report"] = report


def berechne_fvt_daten(stoff: str = "Water", *, T_min=None, T_max=None, v_min=None, v_max=None, n_T_sub: int = 45,
                       n_T_super: int = 25, n_v: int = 41, n_quality: int = 41, critical_refinement: float = 1.5,
                       include_isolines: bool = True, include_homogeneous: bool = True, homogeneous_points: int = 20,
                       include_unstable: bool = True, isolines=None) -> dict:
    """Vollstaendige 3D-f-v-T-Daten: stabile Flaechen, Zweiphasen-Regelflaeche, optionale homogene Ebenen.

    Oberflaechenkoordinaten: (T [K], v_molar [m3/mol], f_mass [J/kg]).
    """
    _check_count(n_v, "n_v")
    raster, _ = generate_pht_grid(stoff, T_min=T_min, T_max=T_max, n_T_sub=n_T_sub, n_T_super=n_T_super, n_p=2,
                                  n_quality=n_quality, critical_refinement=critical_refinement)
    T_sub, T_high = raster["T_sub"], raster["T_high"]
    state = _new_state(stoff)
    Tc, pc = float(state.T_critical()), float(state.p_critical())
    liquid, vapor = _sat_fvt(stoff, T_sub, state)
    _, vc, fc = liquid[-1]

    # Wie in u-v-T: am fluessigen Rand nur moderat verdichten (zu kleine v liegen ausserhalb der HEOS-Gueltigkeit).
    default_v_min = 0.985 * float(np.min(liquid[:-1, 1]))
    default_v_max = max(30.0 * vc, 1.2 * float(np.max(vapor[:-1, 1])))
    v_min = default_v_min if v_min is None else float(v_min)
    v_max = default_v_max if v_max is None else float(v_max)

    if not np.isfinite([v_min, v_max]).all() or not (0 < v_min < np.min(liquid[:, 1]) and v_max > np.max(vapor[:, 1])):
        raise ValueError("f-v-T: v_min/v_max muessen alle Binodalen einschliessen.")

    n_total = 2 * len(liquid) * n_v + len(T_high) * (2 * n_v - 1) + len(liquid) * n_quality

    if n_total > MAX_POINTS:
        raise ValueError(f"f-v-T: {n_total:,} Punkte > MAX_POINTS={MAX_POINTS:,}.")

    sl = np.empty((len(liquid), n_v, 3), dtype=np.float64)
    sg = np.empty_like(sl)

    for i, T in enumerate(np.r_[T_sub, Tc]):
        v_l = np.geomspace(v_min, liquid[i, 1], n_v)
        v_g = np.geomspace(vapor[i, 1], v_max, n_v)

        for j, v in enumerate(v_l[:-1]):
            sl[i, j] = _f_at_tv(state, T, v, CP.iphase_liquid if T < Tc else None)

        for j, v in enumerate(v_g[1:], start=1):
            sg[i, j] = _f_at_tv(state, T, v, CP.iphase_gas if T < Tc else None)

        sl[i, -1] = liquid[i]
        sg[i, 0] = vapor[i]

    high = np.empty((len(T_high), 2 * n_v - 1, 3), dtype=np.float64)
    high[0] = np.concatenate((sl[-1], sg[-1, 1:]), axis=0)

    for i, T in enumerate(T_high[1:], start=1):
        for j, v in enumerate(high[0, :, 1]):
            high[i, j] = _f_at_tv(state, T, v)

    Q = np.linspace(0.0, 1.0, n_quality)
    T_mix = np.broadcast_to(liquid[:, None, 0], (len(liquid), len(Q)))
    v_mix = (1.0 - Q[None, :]) * liquid[:, None, 1] + Q[None, :] * vapor[:, None, 1]
    f_mix = (1.0 - Q[None, :]) * liquid[:, None, 2] + Q[None, :] * vapor[:, None, 2]
    mix = np.stack((T_mix, v_mix, f_mix), axis=-1)
    mix[:, 0] = liquid
    mix[:, -1] = vapor
    mix[-1, :] = liquid[-1]

    meta = {"stoff": stoff, "backend": "HEOS", "Tc": Tc, "pc": pc, "vc": vc, "fc": fc, "v_min": v_min, "v_max": v_max,
            "molar_mass": float(state.molar_mass()),
            "f_definition": "f = u - T s (CoolProp helmholtzmass), massenspezifisch; v molar [m3/mol]"}
    daten = {
        "surface_liquid": sl, "surface_vapor": sg, "surface_high_temperature": high,
        "binodal_liquid": liquid, "binodal_vapor": vapor, "two_phase_surface": mix, "Q_values": Q,
        "critical_point": (Tc, vc, fc), "raster": raster, "meta": meta,
        "isothermen": [], "isobaren": [], "isochoren": [],
    }
    validate_fvt_data(daten)

    if include_homogeneous:
        daten.update(homogeneous_gf_layers(stoff, T_sub, liquid[:-1], vapor[:-1], "fvt", n_points=homogeneous_points,
                                           include_unstable=include_unstable))

    if include_isolines:
        _calculate_isolines(daten, stoff, state, n_quality, n_v, isolines)

    return daten
