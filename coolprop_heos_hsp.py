
"""
Datei: coolprop_heos_hsp.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Berechnung dreidimensionaler h-s-p-Zustandsflächen
durch Zusammenführung von p-h-T- und p-s-T-Daten.

Aufgaben:
- Berechnung von Enthalpie, Entropie und Druck
- Zusammenführung identischer thermodynamischer Zustände
- Erzeugung von Phasengrenzen und Zweiphasenflächen
- Umrechnung der Isolinien in h-s-p-Koordinaten
- Validierung der Flächengeometrie und Randkontinuität

Abhängigkeiten:
- coolprop_heos_pht.py: Enthalpiedaten und Berechnungsraster
- coolprop_heos_pst.py: Entropiedaten und Isolinien
- coolprop_heos_common.py: HEOS-Zustandserzeugung
- CoolProp: Thermodynamische Zustandsberechnungen

Verwendet von:
- main.py: Berechnung des h-s-p-Diagramms
"""


from __future__ import annotations
from coolprop_heos_common import _new_state
import coolprop_heos_pht as pht
import coolprop_heos_pst as pst

import CoolProp.CoolProp as CP
import numpy as np
import logging


LOGGER = logging.getLogger("thermosurface.coolprop_heos_hsp")


def _merge(h_points, s_points):
    """(T,h,p) und (T,s,p) desselben Zustands -> (h,s,p)."""
    h = np.asarray(h_points, dtype=np.float64)
    s = np.asarray(s_points, dtype=np.float64)

    if h.shape != s.shape or h.ndim < 2 or h.shape[-1] != 3:
        raise ValueError("h-s-p: p-h-T und p-s-T haben verschiedene Rasterformen.")

    np.testing.assert_allclose(h[..., 0], s[..., 0], rtol=0, atol=1e-8, err_msg="h-s-p: Temperaturen sind nicht identisch.")
    np.testing.assert_allclose(h[..., 2], s[..., 2], rtol=1e-8, atol=1e-3, err_msg="h-s-p: Druckstuetzstellen sind nicht identisch.")

    return np.stack((h[..., 1], s[..., 1], h[..., 2]), axis=-1)


def _convert_pst_curve(curve, state):
    """Eine p-s-T-Linie (T,s,p) mit PS-Flash nach (h,s,p) uebertragen."""
    curve = np.asarray(curve, dtype=np.float64)

    if curve.ndim != 2 or curve.shape[1] != 3:
        raise ValueError("h-s-p-Isolinie muss die Form (n,3) besitzen.")
    
    out = np.full(curve.shape, np.nan, dtype=np.float64)
    
    for i, (T, s, p) in enumerate(curve):
        if not np.isfinite([T, s, p]).all() or T <= 0 or p <= 0:
    
            continue
        try:
            state.update(CP.PSmass_INPUTS, float(p), float(s))
      
            if not np.isclose(float(state.T()), T, rtol=0.0, atol=1e-3):
                raise ValueError("PS-Flash ergibt abweichende Temperatur.")
      
            h = float(state.hmass())
      
            if np.isfinite(h):
                out[i] = (h, s, p)
       
        except (ValueError, RuntimeError) as exc:
            LOGGER.debug("h-s-p: Isolinienpunkt %d nicht auswertbar: %s", i, exc)
    
    return out


def validate_hsp_data(daten):
    """Form, endliche Werte und gemeinsame Netzraender ueberpruefen."""
    keys = ("surface_liquid", "surface_vapor", "surface_high_temperature","binodal_liquid", "binodal_vapor", "two_phase_surface")
    for key in keys:
        arr = np.asarray(daten[key])
        
        if arr.ndim not in (2, 3) or arr.shape[-1] != 3 or not np.isfinite(arr).all():
            raise ValueError(f"h-s-p: ungueltiges Array {key}: {arr.shape}")
        
        if np.any(arr[..., 2] <= 0):
            raise ValueError(f"h-s-p: nichtpositiver Druck in {key}.")
    
    l, g = daten["binodal_liquid"], daten["binodal_vapor"]
    sl, sg = daten["surface_liquid"], daten["surface_vapor"]
    mix = daten["two_phase_surface"]
    high = daten["surface_high_temperature"]
    q = daten["Q_values"]
    
    if (l.shape != g.shape or sl.shape[0] != len(l) or sg.shape[0] != len(g) or mix.shape != (len(l), len(q), 3) or high.shape[1] != sl.shape[1] + sg.shape[1] - 1):
        raise ValueError("h-s-p: Teilflaechen besitzen inkompatible Dimensionen.")
    
    if len(q) < 2 or q[0] != 0.0 or q[-1] != 1.0:
        raise ValueError("h-s-p: ungueltiges Dampfgehaltsraster.")
    
    for name, actual, expected in (
        ("Fluessigkeitsrand", sl[:, 0], l),
        ("Dampfrand", sg[:, -1], g),
        ("Zweiphasenrand Q=0", mix[:, 0], l),
        ("Zweiphasenrand Q=1", mix[:, -1], g),
        ("Hochtemperatur-Anschluss", high[0],
         np.concatenate((sg[-1], sl[-1, 1:]), axis=0)),):
        
        if not np.array_equal(actual, expected):
            raise ValueError(f"h-s-p: Randkontinuitaet verletzt: {name}.")
    
    if not np.array_equal(l[-1], np.asarray(daten["critical_point"])):
        raise ValueError("h-s-p: Binodale endet nicht im kritischen Punkt.")
    
    if not np.array_equal(g[-1], l[-1]):
        raise ValueError("h-s-p: Beide Binodalzweige muessen bei Tc zusammentreffen.")
    
    if np.any(l[:-1, 0] >= g[:-1, 0]) or np.any(l[:-1, 1] >= g[:-1, 1]):
        raise ValueError("h-s-p: Binodalen zeigen keine positive Verdampfungsdifferenz.")


def berechne_hsp_daten(stoff: str = "Water", *, T_min=None, T_max=None, p_min=None, p_max=None, n_T_sub: int = 45, n_T_super: int = 25, n_p: int = 41, n_quality: int = 41, critical_refinement: float = 1.5, include_isolines: bool = True, isolines=None) -> dict:
    """Vollstaendiges h-s-p-Datenpaket fuer den gemeinsamen PyVista-Renderer."""
    common = dict(T_min=T_min, T_max=T_max, p_min=p_min, p_max=p_max, n_T_sub=n_T_sub, n_T_super=n_T_super, n_p=n_p, n_quality=n_quality, critical_refinement=critical_refinement)
    h_data = pht.berechne_pht_daten(stoff=stoff, **common)
    s_data = pst.berechne_pst_daten(stoff=stoff, include_isolines=include_isolines, isolines=isolines, **common)

    np.testing.assert_array_equal(h_data["Q_values"], s_data["Q_values"])
    np.testing.assert_allclose([h_data["meta"]["Tc"], h_data["meta"]["pc"]], [s_data["meta"]["Tc"], s_data["meta"]["pc"]], rtol=1e-12)

    hc = float(h_data["meta"]["hc"])
    sc = float(s_data["meta"]["sc"])
    pc = float(h_data["meta"]["pc"])
    Tc = float(h_data["meta"]["Tc"])

    keys = ("surface_liquid", "surface_vapor", "surface_high_temperature", "binodal_liquid", "binodal_vapor", "two_phase_surface")
    result = {key: _merge(h_data[key], s_data[key]) for key in keys}

    result.update(
        critical_point=(hc, sc, pc),
        Q_values=h_data["Q_values"],
        raster=h_data["raster"],
        meta={**h_data["meta"], "sc": sc},
        isothermen=[],
        isobaren=[],
        isochoren=[],)

    if include_isolines:
        state = _new_state(stoff)
        for key in ("isothermen", "isobaren", "isochoren"):
            result[key] = [_convert_pst_curve(curve, state) for curve in s_data[key]]
        result["isoline_report"] = s_data.get("isoline_report")

    validate_hsp_data(result)
    LOGGER.debug("h-s-p: %s, Tc=%g K, %d Saettigungszeilen", stoff, Tc, len(result["binodal_liquid"]))
    
    return result
