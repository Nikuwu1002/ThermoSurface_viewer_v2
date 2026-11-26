
"""
Datei: spinodal_vapor_solver.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Numerische Bestimmung der Dampfspinodale auf Grundlage
der homogenen HEOS-Zustandsgleichung.

Aufgaben:
- Berechnung der Druckableitung nach der molaren Dichte
- Suche nach physikalisch plausiblen Druckmaxima
- Numerische Nullstellenbestimmung mittels Brent-Verfahren
- Prüfung der berechneten Spinodalpunkte
- Kennzeichnung nicht berechenbarer Temperaturstützstellen
- Bereitstellung von Diagnosen und Vergleichsdaten

Abhängigkeiten:
- coolprop_heos_common.py: HEOS-Zustandserzeugung
- CoolProp: Homogene Zustands- und Ableitungsberechnungen
- SciPy: Numerische Nullstellenbestimmung
- NumPy: Raster- und Datenverarbeitung

Verwendet von:
- coolprop_heos_pvt.py: Berechnung der Dampfspinodale
"""


from __future__ import annotations

import logging
from time import perf_counter

import numpy as np
from scipy.optimize import brentq

LOGGER = logging.getLogger("thermosurface.spinodal_solver")


def _check_options(T_min, T_max, n_T, n_scan, residual_tolerance):
    if not np.isfinite([T_min, T_max, residual_tolerance]).all() or not (
            0 < T_min < T_max and 0 < residual_tolerance <= 1e-2):
        raise ValueError("Spinodalsolver: Temperaturgrenzen bzw. Residualtoleranz ungueltig.")
    if (type(n_T) is not int or type(n_scan) is not int or
            not 2 <= n_T <= 10000 or not 16 <= n_scan <= 10000 or
            n_T * n_scan > 2_000_000):
        raise ValueError("Spinodalsolver: n_T 2..10000, n_scan 16..10000, "
                         "n_T*n_scan <= 2000000 erforderlich.")


def _homogeneous_derivative(state, cp, T, rho, counters):
    """P-T-Ableitung fuer die fortgesetzte Gasphase; Fehler werden sichtbar."""
    try:
        state.specify_phase(cp.iphase_gas)
        state.update(cp.DmolarT_INPUTS, float(rho), float(T))
        counters["derivative_evaluations"] += 1
        value = float(state.first_partial_deriv(cp.iP, cp.iDmolar, cp.iT))
        return value if np.isfinite(value) else float("nan")
    except (ValueError, RuntimeError, OverflowError):
        counters["derivative_errors"] += 1
        return float("nan")
    finally:
        state.unspecify_phase()


def _homogeneous_pressure(state, cp, T, rho):
    try:
        state.specify_phase(cp.iphase_gas)
        state.update(cp.DmolarT_INPUTS, float(rho), float(T))
        return float(state.p())
    finally:
        state.unspecify_phase()


def _brackets(densities, slopes):
    """Maxima und Minima unterscheiden, auch bei exakt getroffener Raster-Nullstelle.

    Ein ungueliger Gitterpunkt trennt Intervalle; wir ueberbruecken ihn nie.
    """
    found = []
    for j in range(len(densities) - 1):
        f0, f1 = slopes[j:j+2]
        if not np.isfinite([f0, f1]).all():
            continue
        if f0 * f1 < 0:  # Weder Interpolation ueber NaN noch Mehrfachbracket.
            found.append((float(densities[j]), float(densities[j+1]),
                          "maximum" if f0 > 0 else "minimum"))
    # Exact-root-Fall: nur mit aussagekraeftigen unmittelbaren Nachbarn.
    for j in range(1, len(densities)-1):
        a, z, b = slopes[j-1:j+2]
        if z == 0 and np.isfinite([a, b]).all() and a*b < 0:
            found.append((float(densities[j]), float(densities[j]),
                          "maximum" if a > 0 else "minimum"))
    return found


def calc_spinodale_vapor_temperature(
        stoff: str, T_min: float, T_max: float, n_T: int = 250,
        n_scan: int = 200, residual_tolerance: float = 1e-8, *,
        with_diagnostics: bool = False, state_factory=None, cp_module=None):
    """Berechnet eine Dampfspinodale auf der homogenen HEOS-Fortsetzung.

    Gesucht wird das Druckmaximum (dp/drho)_T: positiv -> negativ,
    rho_sat,v < rho_spin,v < rho_c, p_spin,v > p_sat.
    Weitere stationaere Stellen sind KEINE akzeptierten Dampfspinodalpunkte.
    Die Kriterien sind eine begrenzte, fluideigenschaftliche Plausibilitaetspruefung,
    kein allgemeiner Beweis fuer eine beliebige Fluideigenschaftsmodellierung.

    Gibt standardmaessig ein (n_T,3)-Array zurueck; wahlweise (Array, Bericht).
    """
    _check_options(T_min, T_max, n_T, n_scan, residual_tolerance)
    if cp_module is None:
        import CoolProp.CoolProp as cp_module
    if state_factory is None:
        from coolprop_heos_common import _new_state
        state_factory = lambda: _new_state(stoff)
    state = state_factory()
    sat = state_factory()
    Tc = float(state.T_critical())
    pc = float(state.p_critical())
    rho_c = float(state.rhomolar_critical())
    Ttr = float(sat.Ttriple())
    Tmin = float(sat.Tmin())
    Tmax = float(sat.Tmax())
    if not np.isfinite([Tc, pc, rho_c, Ttr, Tmin, Tmax]).all() or min(
            Tc, pc, rho_c, Ttr, Tmin) <= 0:
        raise ValueError(f"{stoff}: nichtendliche/ungueltige HEOS-Kennwerte")
    lower = max(float(T_min), Ttr + 1e-6 * Tc, Tmin + 1e-8 * Tc)
    # Keine fest vorgetaeuschte Spinodalberechnung am kritischen Punkt:
    upper = min(float(T_max), Tmax, Tc - max(1e-5, 1e-6 * Tc))
    if not lower < upper:
        raise ValueError(f"{stoff}: kein subkritischer HEOS-Suchbereich in [{T_min},{T_max}] K")
    temperatures = np.linspace(lower, upper, n_T)
    points = np.column_stack((temperatures, np.full(n_T, np.nan), np.full(n_T, np.nan)))
    diagnostics = {"requested_range": (float(T_min), float(T_max)),
                   "evaluated_range": (float(lower), float(upper)),
                   "requested_points": int(n_T), "valid_points": 0,
                   "derivative_evaluations": 0, "derivative_errors": 0,
                   "root_solver_calls": 0, "found_minima": 0,
                   "found_maxima": 0, "ambiguous_temperatures": [],
                   "missing_temperatures": [], "failures": {},
                   "max_residual": None, "median_residual": None,
                   "residual_tolerance": float(residual_tolerance),
                   "fluid": stoff, "n_scan": int(n_scan)}
    start = perf_counter()
    residuals = []

    for idx, T in enumerate(temperatures):
        T = float(T)
        reason = None
        try:
            sat.update(cp_module.QT_INPUTS, 1.0, T)
            rho_sat = float(sat.rhomolar())
            p_sat = float(sat.p())
            if not (0 < rho_sat < rho_c and p_sat > 0):
                raise ValueError("ungueltige Saettigungsdaten")
            low = rho_sat * (1.0 + 1e-10)
            high = rho_c * (1.0 - 1e-9)
            if low >= high:
                raise ValueError("kein nutzbares Dichteintervall")
            # Im Zweifel adaptiv verfeinern, statt fehlende Nullstellen
            # mathematisch/geometrisch zu erfinden.
            candidates = []
            for factor in (1, 2, 4):
                densities = np.geomspace(low, high, n_scan * factor)
                derivatives = np.array([_homogeneous_derivative(
                    state, cp_module, T, float(r), diagnostics) for r in densities])
                brackets = _brackets(densities, derivatives)
                maxima = [b for b in brackets if b[2] == "maximum"]
                if maxima:
                    candidates = brackets
                    break
            if not candidates:
                reason = "kein klassifizierbarer positiver->negativer Nullstellenwechsel"
                continue
            diagnostics["found_maxima"] += sum(b[2] == "maximum" for b in candidates)
            diagnostics["found_minima"] += sum(b[2] == "minimum" for b in candidates)
            qualified = []
            for a, b, kind in candidates:
                if kind != "maximum":
                    continue
                try:
                    diagnostics["root_solver_calls"] += int(a != b)
                    root = float(a if a == b else brentq(
                        lambda r: _homogeneous_derivative(
                            state, cp_module, T, r, diagnostics), a, b,
                        xtol=max(1e-13, 1e-13 * rho_c), rtol=1e-12,
                        maxiter=100))
                    res = abs(_homogeneous_derivative(
                        state, cp_module, T, root, diagnostics)) * rho_c / pc
                    pressure = _homogeneous_pressure(state, cp_module, T, root)
                    if (np.isfinite([res, pressure]).all() and
                            res <= residual_tolerance and
                            rho_sat < root < rho_c and
                            pressure > p_sat * (1.0 + 1e-11)):
                        qualified.append((root, pressure, res))
                except (ValueError, RuntimeError, OverflowError) as exc:
                    LOGGER.debug("%s: T=%.6f K: Brent-Fehler: %s", stoff, T, exc)
            if len(qualified) != 1:
                if len(qualified) > 1:
                    diagnostics["ambiguous_temperatures"].append(T)
                reason = "keine eindeutige physikalisch plausible Maximum-Nullstelle"
                continue
            root, pressure, residual = qualified[0]
            points[idx, 1:] = (1.0/root, pressure)
            residuals.append(float(residual))
        except (ValueError, RuntimeError, OverflowError) as exc:
            reason = f"HEOS-Auswertung fehlgeschlagen: {exc}"
        finally:
            if reason is not None:
                diagnostics["missing_temperatures"].append(T)
                diagnostics["failures"][str(idx)] = reason

    # Plausibilitaet der Fortsetzung: kein unbemerkter Wechsel auf anderen Ast.
    # Jede verworfene Stuetze bleibt explizit NaN und unterbricht das Rendering.
    for i in range(1, len(points)):
        if not (np.isfinite(points[i, 1:]).all() and
                np.isfinite(points[i-1, 1:]).all()):
            continue
        rho_prev, rho_now = 1/points[i-1, 1], 1/points[i, 1]
        if rho_now < rho_prev * (1 - 1e-3) or points[i, 2] < points[i-1, 2] * (1 - 1e-4):
            diagnostics["missing_temperatures"].append(float(temperatures[i]))
            diagnostics["failures"][str(i)] = "Monotonieverletzung/Astwechsel"
            points[i, 1:] = np.nan

    good = np.isfinite(points).all(axis=1) & (points[:, 1] > 0) & (points[:, 2] > 0)
    diagnostics["valid_points"] = int(np.sum(good))
    diagnostics["success_fraction"] = float(np.mean(good))
    diagnostics["first_valid_T"] = float(points[good, 0][0]) if np.any(good) else None
    diagnostics["last_valid_T"] = float(points[good, 0][-1]) if np.any(good) else None
    diagnostics["max_residual"] = float(max(residuals)) if residuals else None
    diagnostics["median_residual"] = float(np.median(residuals)) if residuals else None
    diagnostics["large_gaps"] = bool(np.any(~good[1:-1]) and np.any(good))
    diagnostics["elapsed_seconds"] = float(perf_counter() - start)
    diagnostics["covers_start"] = bool(good[0])
    diagnostics["covers_end"] = bool(good[-1])
    diagnostics["is_complete"] = bool(np.all(good))
    if not diagnostics["is_complete"]:
        LOGGER.warning("%s: Dampfspinodale %d/%d gueltig, Fehlstellen=%d",
                       stoff, int(sum(good)), n_T, int(sum(~good)))
    return (points, diagnostics) if with_diagnostics else points


def compare_spinodal_overlap(temperature_curve, build_curve, *, critical_T=None):
    """Nur Diagnose: HEOS-build-Ergebnisse zwischen benachbarten NEUEN Werten
    interpolieren. Die Interpolation fliesst NICHT in die Produktivkurve ein.
    """
    new = np.asarray(temperature_curve, dtype=float)
    old = np.asarray(build_curve, dtype=float)
    if (new.ndim != 2 or old.ndim != 2 or new.shape[1] != 3 or old.shape[1] != 3):
        raise ValueError("Vergleich benoetigt (n,3)-Kurven")
    good = np.isfinite(new).all(axis=1) & (new[:, 1] > 0) & (new[:, 2] > 0)
    new = new[good]
    good = np.isfinite(old).all(axis=1) & (old[:, 1] > 0) & (old[:, 2] > 0)
    old = old[good]
    if critical_T is not None:
        new = new[new[:, 0] < critical_T - 0.5]
        old = old[old[:, 0] < critical_T - 0.5]
    if len(new) < 2 or len(old) < 2:
        return {"compared": 0, "max_rel_volume": None, "max_rel_pressure": None}
    # Keine Interpolation ueber verworfene neue Gitterpunkte.
    exact_grid = np.asarray(temperature_curve, dtype=float)
    rel_v, rel_p = [], []
    for T, v, p in old:
        j = int(np.searchsorted(new[:, 0], T))
        if j <= 0 or j >= len(new):
            continue
        left, right = new[j-1], new[j]
        if np.count_nonzero((exact_grid[:, 0] > left[0]) & (exact_grid[:, 0] < right[0])):
            continue
        w = (T - left[0])/(right[0] - left[0])
        v_ref = (1-w)*left[1]+w*right[1]
        p_ref = (1-w)*left[2]+w*right[2]
        rel_v.append(abs(v_ref-v)/v)
        rel_p.append(abs(p_ref-p)/p)
    return {"compared": len(rel_v),
            "max_rel_volume": float(max(rel_v)) if rel_v else None,
            "max_rel_pressure": float(max(rel_p)) if rel_p else None}
