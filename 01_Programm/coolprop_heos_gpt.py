
"""
Datei: coolprop_heos_gpt.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Berechnung dreidimensionaler g-p-T-Zustandsflächen mit CoolProp HEOS.

Aufgaben:
- Berechnung der Gibbs-Energie und des kritischen Zustands
- Erzeugung von Flüssigkeits-, Dampf- und Hochtemperaturflächen
- Berechnung der Sättigungslinie und der Isolinien
- Bereitstellung homogener metastabiler und instabiler Zustände
- Validierung der berechneten Zustandsflächen

Abhängigkeiten:
- coolprop_heos_common.py: HEOS-Zustandsberechnungen
- coolprop_heos_pht.py: Rastererzeugung und Hilfsfunktionen
- isoline_levels.py: Isoliniengenerierung
- zustandsflaeche.py: Druckachsenskalierung
- limits.py: Begrenzung der Rasterpunktzahl

Verwendet von:
- main.py: Dynamischer Aufruf der g-p-T-Berechnung
"""


from __future__ import annotations

from coolprop_heos_common import (_new_state, critical_gf, gibbs_helmholtz_tv, homogeneous_gf_layers,
                                  pt_update_robust)
from coolprop_heos_pht import _check_count, _geom, generate_pht_grid
from limits import MAX_POINTS
from zustandsflaeche import PRESSURE_AXES

import CoolProp.CoolProp as CP
import numpy as np
import logging


LOGGER = logging.getLogger("thermosurface.coolprop_heos_gpt")

# Relative Toleranz fuer g_l = g_g auf der Sattigungslinie (Rundung der HEOS-Iteration)
SATURATION_G_RTOL = 1e-7


def critical_point_gpt(stoff: str) -> tuple[float, float, float]:
    """(Tc [K], pc [Pa], gc [J/kg]) aus dem kritischen HEOS-Zustand."""
    Tc, pc, _vc, gc, _fc = critical_gf(stoff)
    return Tc, pc, gc


def _g_pt(state, T: float, p: float, phase=None) -> tuple[float, float, float]:
    """Einphasiger PT-Flash -> (T, p, g_mass); optionale stabile Phasenvorgabe."""
    try:
        pt_update_robust(state, T, p, phase)
        g = float(state.gibbsmass())
    except (ValueError, RuntimeError) as exc:
        raise ValueError(f"g-p-T: HEOS-PT-Flash fehlgeschlagen bei T={T:.9g} K, p={p:.9g} Pa: {exc}") from exc

    if not np.isfinite(g):
        raise ValueError(f"g-p-T: ungueltiges g bei T={T}, p={p}.")

    return float(T), float(p), g


def g_eos(T: float, p: float, stoff: str) -> float:
    """Spezifische freie Enthalpie [J/kg] aus (T,p); NaN, wenn nicht berechenbar."""
    if not (np.isfinite([T, p]).all() and T > 0 and p > 0):
        return float("nan")

    try:
        return _g_pt(_new_state(stoff), T, p)[2]
    except ValueError:
        return float("nan")


# ---------------------------------------------------------------- Sattigungslinie

def multi_binodale_gpt(stoff: str, T_values, include_critical: bool = True):
    """Sattigungslinie als (liquid, vapor), je (n,3) mit (T, p_sat, g_sat).

    Beide Aeste benutzen EXAKT denselben Druck. g wird fuer Siede- und Taulinie getrennt aus HEOS gelesen;
    die Gleichheit g_l = g_g (Gleichgewicht) wird geprueft, nicht erzwungen.
    """
    state = _new_state(stoff)
    Tc, pc, gc = critical_point_gpt(stoff)
    temperatures = np.asarray(T_values, dtype=float)

    if (temperatures.ndim != 1 or len(temperatures) < 1 or not np.isfinite(temperatures).all()
            or np.any(np.diff(temperatures) <= 0)
            or not np.all((temperatures > state.Ttriple()) & (temperatures < Tc))):
        raise ValueError("g-p-T: T_values muessen streng steigend zwischen Ttriple und Tc liegen.")

    gas_constant = 8.31446261815324 / float(state.molar_mass())
    liquid, vapor = [], []

    for T in temperatures:
        try:
            state.update(CP.QT_INPUTS, 0.0, float(T))
            p_l, g_l = float(state.p()), float(state.gibbsmass())
            state.update(CP.QT_INPUTS, 1.0, float(T))
            p_g, g_g = float(state.p()), float(state.gibbsmass())
        except (ValueError, RuntimeError) as exc:
            raise ValueError(f"g-p-T: HEOS-Saettigung bei T={T:.9g} K nicht auswertbar: {exc}") from exc

        if not (np.isfinite([p_l, p_g, g_l, g_g]).all() and p_l > 0 and np.isclose(p_l, p_g, rtol=1e-8, atol=1e-3)):
            raise ValueError(f"g-p-T: ungueltige Saettigungswerte bei T={T:.9g} K.")

        if abs(g_l - g_g) > SATURATION_G_RTOL * max(abs(g_l), gas_constant * float(T)):
            raise ValueError(f"g-p-T: g_l != g_g auf der Saettigungslinie bei T={T:.9g} K "
                             f"({g_l} vs {g_g}); Gleichgewichtsbedingung verletzt.")

        liquid.append((float(T), p_l, g_l))
        vapor.append((float(T), p_l, g_g))

    if include_critical:
        liquid.append((Tc, pc, gc))
        vapor.append((Tc, pc, gc))

    return (np.asarray(liquid, dtype=np.float64).reshape(-1, 3),
            np.asarray(vapor, dtype=np.float64).reshape(-1, 3))


# ---------------------------------------------------------------- randkonforme einphasige Flaechen

def _single_phase_surfaces(stoff: str, raster: dict, liquid: np.ndarray, vapor: np.ndarray):
    """Fluessigkeits-, Dampf- und Hochtemperaturnetz (T,p)->g; Randknoten identisch zu den Binodalen."""
    T_sub, T_high = raster["T_sub"], raster["T_high"]
    p_min, p_max, n_p = raster["p_min"], raster["p_max"], raster["n_p"]
    Tc, pc, gc = liquid[-1]
    l_surf = np.empty((len(T_sub) + 1, n_p, 3), dtype=np.float64)
    g_surf = np.empty_like(l_surf)
    state = _new_state(stoff)

    for i, T in enumerate(T_sub):
        psat = float(liquid[i, 1])
        p_g = _geom(p_min, psat, n_p)
        p_l = _geom(psat, p_max, n_p)

        for j, p in enumerate(p_g[:-1]):
            g_surf[i, j] = _g_pt(state, float(T), float(p), CP.iphase_gas)

        for j, p in enumerate(p_l[1:], start=1):
            l_surf[i, j] = _g_pt(state, float(T), float(p), CP.iphase_liquid)

        g_surf[i, -1] = vapor[i]
        l_surf[i, 0] = liquid[i]

    p_g_crit = _geom(p_min, pc, n_p)
    p_l_crit = _geom(pc, p_max, n_p)

    for j, p in enumerate(p_g_crit[:-1]):
        g_surf[-1, j] = _g_pt(state, Tc, float(p))

    for j, p in enumerate(p_l_crit[1:], start=1):
        l_surf[-1, j] = _g_pt(state, Tc, float(p))

    g_surf[-1, -1] = (Tc, pc, gc)
    l_surf[-1, 0] = (Tc, pc, gc)

    high = np.empty((len(T_high), 2 * n_p - 1, 3), dtype=np.float64)
    high[0] = np.concatenate((g_surf[-1], l_surf[-1, 1:]), axis=0)
    pressures_high = high[0, :, 1].copy()

    for i, T in enumerate(T_high[1:], start=1):
        for j, p in enumerate(pressures_high):
            high[i, j] = _g_pt(state, float(T), float(p))

    return l_surf, g_surf, high


# ---------------------------------------------------------------- Isolinien

def _isobare(state, p, T_min, T_max, n_T):
    """g(T) bei konstantem Druck; Fluessig-/Dampfast mit vorgegebener Phase, exakter Sattigungsknoten."""
    Tsat, sat_row = None, None

    if p < float(state.p_critical()):
        try:
            state.update(CP.PQ_INPUTS, float(p), 0.0)
            Tsat = float(state.T())
            sat_row = (Tsat, float(p), float(state.gibbsmass()))
        except (ValueError, RuntimeError):
            Tsat = None

    points = []

    for T in np.linspace(T_min, T_max, n_T):
        if Tsat is not None and np.isclose(T, Tsat, rtol=1e-9, atol=0.0):
            continue
        phase = None if Tsat is None else (CP.iphase_liquid if T < Tsat else CP.iphase_gas)
        points.append(_g_pt(state, float(T), float(p), phase))

    if Tsat is not None and T_min <= Tsat <= T_max:
        points.append(sat_row)

    return np.asarray(sorted(points), dtype=np.float64)


def _isotherme(state, T, p_min, p_max, n_p):
    """g(p) bei konstanter Temperatur; Dampfast -> Sattigungsknoten (Knick) -> Fluessigkeitsast (T < Tc), sonst ein Ast."""
    sat = None

    if T < float(state.T_critical()):
        try:
            state.update(CP.QT_INPUTS, 0.0, float(T))
            p_sat = float(state.p())
            if p_min < p_sat < p_max:
                sat = p_sat
        except (ValueError, RuntimeError):
            sat = None

    def branch(pressures, phase):
        rows = []
        for p in pressures:
            try:
                rows.append(_g_pt(state, float(T), float(p), phase))
            except ValueError:
                rows.append((float(T), float(p), np.nan))
        return rows

    if sat is None:
        return np.asarray(branch(_geom(p_min, p_max, 2 * n_p - 1), None), dtype=np.float64)

    # Randknoten exakt wie in den Gleichgewichtsflaechen: ein gemeinsamer Sattigungsknoten (g_l = g_g).
    vapor = branch(_geom(p_min, sat, n_p)[:-1], CP.iphase_gas)
    liquid = branch(_geom(sat, p_max, n_p)[1:], CP.iphase_liquid)
    state.update(CP.QT_INPUTS, 0.0, float(T))
    node = (float(T), sat, float(state.gibbsmass()))
    return np.asarray(vapor + [node] + liquid, dtype=np.float64)


def _calculate_isolines(daten, stoff, state, isolines=None):
    """Isolinien entlang geometrisch verteilter Levels (TS-VP03): Isothermen ueber p, Isobaren per PT-Flash mit
    Phasenvorgabe, Isochoren per DmolarT. Das Flaechenraster bleibt unveraendert."""
    from isoline_levels import resolve_levels, build_isolines, V_DOMAIN_FACTORS
    Tc, pc, _gc = daten["critical_point"]
    T_sub, T_high = daten["raster"]["T_sub"], daten["raster"]["T_high"]
    p_min, p_max, n_p = daten["raster"]["p_min"], daten["raster"]["p_max"], daten["raster"]["n_p"]
    T_lo, T_hi = float(T_sub[0]), float(T_high[-1])
    vc = 1.0 / float(state.rhomolar_critical())
    temperatures = np.r_[T_sub, T_high]
    levels, report = resolve_levels(isolines, {
        "isothermen": (T_lo, T_hi), "isobaren": (float(p_min), float(p_max)),
        "isochoren": (V_DOMAIN_FACTORS[0] * vc, V_DOMAIN_FACTORS[1] * vc)})

    def isochore(v):
        curve = np.full((len(temperatures), 3), np.nan, dtype=np.float64)

        for i, T in enumerate(temperatures):
            try:
                p, g, _f = gibbs_helmholtz_tv(state, T, v)
            except ValueError as exc:
                LOGGER.debug("g-p-T-Isochore: %s", exc)
                continue

            if p_min <= p <= p_max:
                curve[i] = (T, p, g)

        return curve

    curves = build_isolines(levels, report, {
        "isothermen": lambda T: _isotherme(state, T, float(p_min), float(p_max), int(n_p)),
        "isobaren": lambda p: _isobare(state, p, T_lo, T_hi, 81),
        "isochoren": isochore}, log=LOGGER.warning)
    daten["isothermen"], daten["isobaren"], daten["isochoren"] = (
        curves["isothermen"], curves["isobaren"], curves["isochoren"])
    daten["isoline_report"] = report


# ---------------------------------------------------------------- Validierung

def validate_gpt_data(daten):
    """Topologie- und Konsistenzpruefung auf den g-p-T-Datenarrays."""
    l, g = daten["binodal_liquid"], daten["binodal_vapor"]
    sl, sg, high = daten["surface_liquid"], daten["surface_vapor"], daten["surface_high_temperature"]

    for name in ("surface_liquid", "surface_vapor", "surface_high_temperature", "binodal_liquid", "binodal_vapor"):
        arr = np.asarray(daten[name])

        if arr.ndim not in (2, 3) or arr.shape[-1] != 3 or not np.isfinite(arr).all():
            raise ValueError(f"g-p-T: ungueltiges Array {name}: {arr.shape}")

        if np.any(arr[..., 0] <= 0) or np.any(arr[..., 1] <= 0):
            raise ValueError(f"g-p-T: nichtpositive T/p-Koordinaten in {name}.")

    if (l.shape != g.shape or sl.shape[0] != len(l) or sg.shape[0] != len(g)
            or high.shape[1] != sl.shape[1] + sg.shape[1] - 1):
        raise ValueError("g-p-T: Teilflaechen haben inkompatible Dimensionen.")

    for name, actual, expected in (("Fluessigkeitsrand", sl[:, 0], l), ("Dampfrand", sg[:, -1], g),
                                   ("Hochtemperatur-Anschluss", high[0], np.concatenate((sg[-1], sl[-1, 1:]), axis=0))):
        if not np.array_equal(actual, expected):
            raise ValueError(f"g-p-T: Randkontinuitaet verletzt: {name}.")

    if not np.array_equal(l[-1], g[-1]) or not np.array_equal(l[-1], np.asarray(daten["critical_point"])):
        raise ValueError("g-p-T: kein gemeinsamer kritischer Endpunkt.")

    if not np.array_equal(l[:, :2], g[:, :2]):
        raise ValueError("g-p-T: Siede- und Taulinie muessen dieselben (T,p) besitzen.")

    if not np.all(np.diff(l[:, 0]) > 0) or not np.all(np.diff(l[:, 1]) > 0):
        raise ValueError("g-p-T: Sattigungslinie nicht streng monoton in T und p.")

    # (dg/dp)_T = v > 0: g waechst auf jedem Isothermenast mit dem Druck
    for name, surface in (("surface_liquid", sl), ("surface_vapor", sg), ("surface_high_temperature", high)):
        if np.any(np.diff(surface[..., 2], axis=1) <= 0) or np.any(np.diff(surface[..., 1], axis=1) <= 0):
            raise ValueError(f"g-p-T: (dg/dp)_T = v > 0 verletzt in {name}.")


# ---------------------------------------------------------------- Gesamtberechnung

def berechne_gpt_daten(stoff: str = "Water", *, T_min=None, T_max=None, p_min=None, p_max=None, n_T_sub: int = 45,
                       n_T_super: int = 25, n_p: int = 41, n_quality: int = 41, critical_refinement: float = 1.5,
                       include_isolines: bool = True, include_homogeneous: bool = True, homogeneous_points: int = 20,
                       include_unstable: bool = True, pressure_axis: str = "log",
                       auto_pressure_axis: bool = False, isolines=None) -> dict:
    """Vollstaendige 3D-g-p-T-Daten: stabile Gleichgewichtsflaechen plus optionale homogene Ebenen.

    Oberflaechenkoordinaten: (T [K], p [Pa], g_mass [J/kg]).
    n_quality wird nur fuer das gemeinsame Temperaturraster angenommen (kein Zweiphasengitter in g-p-T).
    """
    if pressure_axis not in PRESSURE_AXES:
        raise ValueError("pressure_axis muss eines von " + ", ".join(PRESSURE_AXES) + " sein")

    _check_count(n_quality, "n_quality")
    raster, meta = generate_pht_grid(stoff, T_min=T_min, T_max=T_max, p_min=p_min, p_max=p_max, n_T_sub=n_T_sub,
                                     n_T_super=n_T_super, n_p=n_p, n_quality=n_quality,
                                     critical_refinement=critical_refinement)
    T_sub = raster["T_sub"]
    Tc, pc, gc = critical_point_gpt(stoff)
    bin_l, bin_g = multi_binodale_gpt(stoff, T_sub, include_critical=True)
    surf_l, surf_g, surf_high = _single_phase_surfaces(stoff, raster, bin_l, bin_g)
    state = _new_state(stoff)
    meta = {**meta, "gc": gc, "vc": 1.0 / float(state.rhomolar_critical()), "molar_mass": float(state.molar_mass()),
            "g_definition": "g = h - T s (CoolProp gibbsmass), massenspezifisch"}

    daten = {
        "surface_liquid": surf_l, "surface_vapor": surf_g, "surface_high_temperature": surf_high,
        "binodal_liquid": bin_l, "binodal_vapor": bin_g,
        "critical_point": (Tc, pc, gc), "raster": raster, "meta": meta,
        "isothermen": [], "isobaren": [], "isochoren": [],
        "pressure_axis": pressure_axis, "pressure_axis_notice": None,
        "two_phase_note": ("In g-p-T entartet das Zweiphasengebiet zur Sattigungslinie (g_l = g_g bei p_sat(T)); "
                           "es existiert keine Zweiphasenflaeche."),
    }
    validate_gpt_data(daten)

    if include_homogeneous:
        daten.update(homogeneous_gf_layers(stoff, T_sub, bin_l[:-1], bin_g[:-1], "gpt", n_points=homogeneous_points,
                                           include_unstable=include_unstable))
        info = daten["homogeneous_info"]
        outside = 0

        for key in ("metastable_liquid_surface", "metastable_vapor_surface", "unstable_surface"):
            if key in daten:
                arr = daten[key]
                finite = np.isfinite(arr).all(axis=-1)
                p = arr[..., 1]
                outside += int(np.count_nonzero(finite & ((p < raster["p_min"]) | (p > raster["p_max"]))))

        nonpositive = sum(s["nonpositive_pressure_nodes"] for s in info["surfaces"].values())
        info["nonpositive_pressure_nodes"] = nonpositive
        info["outside_pressure_window_nodes"] = outside

        if outside and pressure_axis == "log":
            if auto_pressure_axis:
                daten["pressure_axis"] = "symlog"
                daten["pressure_axis_notice"] = (
                    f"Druckachse automatisch von 'log' auf 'symlog' umgestellt: {outside} homogene Knoten liegen "
                    f"ausserhalb des Druckfensters [{raster['p_min']:.4g}, {raster['p_max']:.4g}] Pa "
                    f"({nonpositive} davon p <= 0).")
            else:
                daten["pressure_axis_notice"] = (
                    f"{outside} homogene Knoten ({nonpositive} mit p <= 0) liegen ausserhalb des Druckfensters "
                    f"[{raster['p_min']:.4g}, {raster['p_max']:.4g}] Pa der Log-Achse und werden nicht dargestellt "
                    "(Daten unveraendert; diagram_settings.gpt.pressure_axis: symlog zeigt alle).")

    if include_isolines:
        _calculate_isolines(daten, stoff, state, isolines)

    return daten
