"""
Datei: coolprop_heos_common.py
Projekt: ThermoSurface
Kategorie: Thermodynamik / HEOS-Basis

Beschreibung:
Stellt gemeinsam genutzte CoolProp-HEOS-Zustandsfunktionen bereit.

Aufgaben:
- Erzeugt AbstractState-Objekte und ermittelt Stoff- und kritische Daten.
- Berechnet g und f sowie homogene Fortsetzungen für Energieflächen.
Einheiten: T [K], v [m³/mol], p [Pa], g und f [J/kg].

Abhängigkeiten:
- Projektmodule (direkt): homogeneous_pvt.py
- Bibliotheken: CoolProp, numpy

Verwendet von:
- Projektmodule: coolprop_heos_fvt.py, coolprop_heos_pvt.py, fvt_tangent.py,
  homogeneous_pvt.py, main.py, spinodal_vapor_solver.py, zustandsflaeche.py
"""

from CoolProp import AbstractState
import CoolProp.CoolProp as CP
import numpy as np
import logging


LOGGER = logging.getLogger("thermosurface.coolprop_heos")

# ============================================================
# Hilfsfunktion: neuen HEOS-Zustand für einen Stoff anlegen
# ============================================================

def _new_state(stoff):
    """Erzeugt ein neues AbstractState-Objekt mit dem HEOS-Backend.

    Ein AbstractState kapselt den thermodynamischen Zustand eines Stoffes
    und wird im Folgenden über .update(...) immer wieder auf einen neuen
    Zustandspunkt (z.B. gegebene Dichte+Temperatur) gesetzt. Das ist
    deutlich effizienter, als für jeden Punkt CP.PropsSI(...) neu
    aufzurufen, da das AbstractState-Objekt den Stoff nur einmal lädt.
    """
    return AbstractState("HEOS", stoff)


# Relative Naehe zu Tc, ab der (T,p) als ueberkritisch/kritisch behandelt wird. Gerasterte Tc-Werte
# (Rundung, float-Darstellung) liegen minimal unter oder ueber dem HEOS-Wert.
CRITICAL_T_TOLERANCE = 1e-6


def pt_update_robust(state, T, p, phase=None):
    """PT-Flash mit Rueckfall fuer den kritischen Rand; wirft ValueError/RuntimeError nur bei Totalausfall.

    Mit explizit vorgegebener Phase wird genau diese verwendet. Ohne Phase scheitert HEOS bei
    T ~ Tc und p > pc stoffabhaengig (z.B. CO2, Argon, Methan: 'p is not a valid number');
    dort wird auf die ueberkritische Phasenvorgabe zurueckgefallen. Unterkritische Zustaende
    werden nie auf eine andere Phase umgelenkt.
    """
    candidates = [phase]
    if phase is None and T >= state.T_critical() * (1.0 - CRITICAL_T_TOLERANCE):
        candidates += [CP.iphase_supercritical, CP.iphase_supercritical_gas,
                       CP.iphase_supercritical_liquid]
    last_error = None
    for candidate in candidates:
        try:
            if candidate is not None:
                state.specify_phase(candidate)
            state.update(CP.PT_INPUTS, float(p), float(T))
            return
        except (ValueError, RuntimeError) as exc:
            last_error = exc
        finally:
            if candidate is not None:
                state.unspecify_phase()
    raise last_error

# ============================================================
# 1. Stoffdaten
# ============================================================

def stoffdaten(stoff):
    """Stoffdaten mit kompatibler Rückgabe (T_crit, p_crit, M, omega).

    HEOS benötigt omega nicht für die p-v-T-Auswertung.
    """
    state = _new_state(stoff)

    T_crit = state.T_critical()
    p_crit = state.p_critical()
    M = state.molar_mass()

    # omega bleibt aus Gründen der bestehenden Rückgabeschnittstelle erhalten.
    omega = CP.PropsSI("acentric", stoff)

    return T_crit, p_crit, M, omega


def critical_point(stoff):
    """Liefert den kritischen Punkt direkt aus HEOS: T_crit, p_crit, v_crit.

    v_crit ergibt sich aus der kritischen molaren Dichte rhomolar_critical,
    die HEOS aus der Reinstoff-Referenzgleichung kennt:
        v_crit = 1 / rho_molar_crit
    Die kritischen Größen stammen aus dem gewählten HEOS-Stoffmodell.
    """
    state = _new_state(stoff)
    T_crit = state.T_critical()
    p_crit = state.p_critical()
    v_crit = 1.0 / state.rhomolar_critical()

    return T_crit, p_crit, v_crit

def fluid_limits(stoff):
    """(T_triple, T_min_HEOS, T_max_HEOS) [K]; HEOS ist alleinige Quelle."""
    state = _new_state(stoff)
    
    return state.Ttriple(), state.Tmin(), state.Tmax()


# ============================================================
# Freie Energien g und f (massenspezifisch) und homogene Ebenen
# fuer die Diagramme g-p-T und f-v-T
# ============================================================

def gibbs_helmholtz_tv(state, T, v_molar, phase=None):
    """Ein HEOS-Zustand (T, v_molar) -> (p [Pa], g [J/kg], f [J/kg]).

    g und f kommen direkt aus CoolProp (gibbsmass, helmholtzmass) und stammen damit aus demselben
    Zustand wie p. Mit `phase` wird die homogene Zustandsgleichung ohne Gleichgewichts-Flash
    ausgewertet (metastabile/instabile Zustaende); ohne `phase` liefert HEOS im Nassdampfgebiet den
    Gleichgewichtszustand.
    """
    T, v = float(T), float(v_molar)
    if not (np.isfinite([T, v]).all() and T > 0 and v > 0):
        raise ValueError("T und molares Volumen muessen endlich und positiv sein.")

    try:
        if phase is not None:
            state.specify_phase(phase)
        try:
            state.update(CP.DmolarT_INPUTS, 1.0 / v, T)
            p, g, f = float(state.p()), float(state.gibbsmass()), float(state.helmholtzmass())
        finally:
            if phase is not None:
                state.unspecify_phase()
    except (ValueError, RuntimeError) as exc:
        raise ValueError(f"HEOS bei T={T:.8g} K, v={v:.8g} m3/mol: {exc}") from exc

    if not np.isfinite([p, g, f]).all():
        raise ValueError(f"Nicht endliche Zustandsgroessen bei T={T}, v={v}.")

    return p, g, f


def critical_gf(stoff):
    """(Tc, pc, vc_molar, gc, fc) am kritischen HEOS-Zustand (DmolarT mit rhomolar_critical)."""
    state = _new_state(stoff)
    Tc, pc, rho_c = float(state.T_critical()), float(state.p_critical()), float(state.rhomolar_critical())

    if not (np.isfinite([Tc, pc, rho_c]).all() and min(Tc, pc, rho_c) > 0):
        raise ValueError("HEOS liefert ungueltige kritische Stoffdaten.")

    _, g, f = gibbs_helmholtz_tv(state, Tc, 1.0 / rho_c)

    return Tc, pc, 1.0 / rho_c, g, f


HOMOGENEOUS_LAYER_KEYS = ("metastable_liquid_surface", "metastable_vapor_surface", "unstable_surface")


def homogeneous_gf_layers(stoff, T_sub, sat_liquid, sat_vapor, mode, *, n_points=20, include_unstable=True,
                          n_scan=200):
    """Spinodalen und homogene (metastabile/instabile) Flaechen fuer g-p-T bzw. f-v-T.

    Wiederverwendet werden die validierten p-v-T-Verfahren (`homogeneous_pvt.solve_spinodal_branch` und
    `build_homogeneous_surfaces`: Spinodal-Solver, Stabilitaetsklassifikation, Randdruckpruefung, Luecken statt
    Interpolation). Die Knoten (T, v, p) werden anschliessend mit derselben HEOS-Zustandsgleichung
    (DmolarT, Phase vorgegeben) in g und f umgerechnet. Es wird nichts interpoliert oder geglaettet.

    mode='gpt': Spalten (T, p, g);  mode='fvt': Spalten (T, v, f).
    sat_liquid/sat_vapor: (n_T,3) exakte Binodalwerte in den Spalten des Modus (Randspalten der metastabilen
    Flaechen werden exakt auf diese Werte gesetzt).
    Rueckgabe: dict mit Spinodalen (inkl. kritischem Endpunkt), den drei Flaechen, `homogeneous_v`
    (molares Volumen je Knoten, fuer die Zustandsrekonstruktion) und `homogeneous_info`.
    """
    if mode not in ("gpt", "fvt"):
        raise ValueError("mode muss 'gpt' oder 'fvt' sein.")
    from homogeneous_pvt import build_homogeneous_surfaces, solve_spinodal_branch

    T_sub = np.asarray(T_sub, dtype=float)
    sat_l = np.asarray(sat_liquid, dtype=float)
    sat_g = np.asarray(sat_vapor, dtype=float)
    state = _new_state(stoff)
    Tc, pc, vc, gc, fc = critical_gf(stoff)

    def native(T, v, p, g, f):
        return (T, p, g) if mode == "gpt" else (T, v, f)

    # Binodalen im p-v-T-Format fuer den bewaehrten Builder
    tvp_l = np.empty((len(T_sub), 3))
    tvp_g = np.empty_like(tvp_l)
    for i, T in enumerate(T_sub):
        for row, q in ((tvp_l, 0.0), (tvp_g, 1.0)):
            state.update(CP.QT_INPUTS, q, float(T))
            row[i] = (float(T), 1.0 / float(state.rhomolar()), float(state.p()))

    spin_l, rep_l = solve_spinodal_branch(stoff, T_sub, "liquid", n_scan=n_scan)
    spin_g, rep_g = solve_spinodal_branch(stoff, T_sub, "vapor", n_scan=n_scan)
    surfaces_tvp, builder_report = build_homogeneous_surfaces(
        stoff, tvp_l, tvp_g, spin_l, spin_g, n_points=n_points, include_unstable=include_unstable)

    phases = {"metastable_liquid_surface": CP.iphase_liquid, "metastable_vapor_surface": CP.iphase_gas,
              "unstable_surface": CP.iphase_gas}
    info = {"spinodal_liquid": {k: rep_l[k] for k in ("requested_points", "valid_points", "missing_temperatures")},
            "spinodal_vapor": {k: rep_g[k] for k in ("requested_points", "valid_points", "missing_temperatures")},
            "surfaces": {}, "mode": mode, "n_points": int(n_points)}

    def convert(tvp, phase):
        out = np.full(tvp.shape, np.nan)
        f_nodes = np.full(tvp.shape[:-1], np.nan)
        for idx in zip(*np.nonzero(np.isfinite(tvp).all(axis=-1))):
            T, v, _p = tvp[idx]
            try:
                p, g, f = gibbs_helmholtz_tv(state, T, v, phase)
            except ValueError:
                continue
            out[idx] = native(T, v, p, g, f)
            f_nodes[idx] = f
        return out, f_nodes

    # Gleichgewichtsbedingung: bei gleichem (T, v) ist f_homogen >= f_Gleichgewicht (gemeinsame Tangente der Binodalpunkte).
    # Knoten darunter sind Artefakte der Mehrparameter-Zustandsgleichung im tiefen instabilen Bereich und werden
    # nicht dargestellt (nur Klassifikation, kein Eingriff in Werte).
    f_sat = np.empty((len(T_sub), 2))
    for i, T in enumerate(T_sub):
        for k, q in enumerate((0.0, 1.0)):
            state.update(CP.QT_INPUTS, q, float(T))
            f_sat[i, k] = float(state.helmholtzmass())
    gas_constant = 8.31446261815324 / float(state.molar_mass())

    result = {"homogeneous_v": {}}
    for key in HOMOGENEOUS_LAYER_KEYS:
        tvp = surfaces_tvp[key]
        arr, f_nodes = convert(tvp, phases[key])
        finite_before = np.isfinite(arr).all(axis=-1)
        v = tvp[..., 1]
        tangent = (f_sat[:, 0, None] + (v - tvp_l[:, 1, None]) / (tvp_g[:, 1, None] - tvp_l[:, 1, None])
                   * (f_sat[:, 1, None] - f_sat[:, 0, None]))
        tolerance = 1e-7 * np.maximum(np.abs(f_nodes), gas_constant * tvp[..., 0])
        below_hull = finite_before & (f_nodes < tangent - tolerance)
        arr[below_hull] = np.nan
        # Randspalten der metastabilen Flaechen = exakte Binodalknoten
        if key == "metastable_liquid_surface":
            same = np.isfinite(arr[:, 0]).all(axis=-1) & (tvp[:, 0, 1] == tvp_l[:, 1])
            arr[same, 0] = sat_l[same]
        elif key == "metastable_vapor_surface":
            same = np.isfinite(arr[:, -1]).all(axis=-1) & (tvp[:, -1, 1] == tvp_g[:, 1])
            arr[same, -1] = sat_g[same]
        if not include_unstable and key == "unstable_surface":
            continue
        finite = np.isfinite(arr).all(axis=-1)
        result[key] = arr
        result["homogeneous_v"][key] = np.where(finite, tvp[..., 1], np.nan)
        p_nodes = np.where(finite, tvp[..., 2], np.nan)
        info["surfaces"][key] = {
            "valid_nodes": int(finite.sum()), "nodes": int(finite.size),
            "below_equilibrium_hull_nodes": int(below_hull.sum()),
            "nonpositive_pressure_nodes": int(np.count_nonzero(finite & (p_nodes <= 0))),
            "gap_rows": len(builder_report[key].get("gap_rows", {})) if key == "unstable_surface" else 0}

    # Spinodalkurven mit exaktem kritischem Endpunkt (dort laufen Binodale und Spinodalen zusammen)
    for name, spin, phase in (("spinodal_liquid", spin_l, CP.iphase_liquid), ("spinodal_vapor", spin_g, CP.iphase_gas)):
        curve = np.full((len(spin) + 1, 3), np.nan)
        for i, (T, v, _p) in enumerate(spin):
            if np.isfinite([v, _p]).all():
                try:
                    p, g, f = gibbs_helmholtz_tv(state, T, v, phase)
                    curve[i] = native(T, v, p, g, f)
                except ValueError:
                    pass
        curve[-1] = native(Tc, vc, pc, gc, fc)
        result[name] = curve
    result["homogeneous_info"] = info
    return result

