"""
Datei: homogeneous_pvt.py
Projekt: ThermoSurface
Kategorie: Thermodynamik / homogene p-v-T-Flächen

Beschreibung:
Berechnet und validiert homogene HEOS-Fortsetzungen um die Spinodalen.

Aufgaben:
- Sucht Spinodaläste und bestimmt metastabile und instabile Stützstellen.
- Prüft Flächenabdeckung, Randkontinuität und Netz-Topologie numerisch.
Einheiten: T [K], v [m³/mol], p [Pa].

Abhängigkeiten:
- Projektmodule (direkt): coolprop_heos_common.py, zustandsflaeche.py
- Bibliotheken: CoolProp, numpy, scipy

Verwendet von:
- Projektmodule: coolprop_heos_common.py, main.py, visualization.py, zustandsflaeche.py
"""

from __future__ import annotations
import numpy as np
from scipy.optimize import brentq
from time import perf_counter
from zustandsflaeche import PRESSURE_AXES


def _phase(cp, branch):
    return cp.iphase_liquid if branch == 'liquid' else cp.iphase_gas


def homogeneous_point(state, cp, T, v, branch):
    """p und normierte dp/drho fuer vorgeschriebene homogene Phase.

    Auch im metastabilen Bereich keinen Gleichgewichts-Flash ausloesen.
    """
    if not (np.isfinite([T, v]).all() and T > 0 and v > 0):
        raise ValueError('Ungueltige Temperatur/Volumen')
    state.specify_phase(_phase(cp, branch))
    try:
        state.update(cp.DmolarT_INPUTS, 1.0 / float(v), float(T))
        p = float(state.p())
        slope = float(state.first_partial_deriv(cp.iP, cp.iDmolar, cp.iT))
    finally:
        state.unspecify_phase()
    if not np.isfinite([p, slope]).all():
        raise ValueError('Homogene EOS nicht auswertbar')
    return p, slope


def solve_spinodal_branch(stoff, temperatures, branch, n_scan=200,
                          residual_tolerance=1e-8, *, state_factory=None, cp_module=None):
    """Saettigungsnaechstes Ableitungsextremum auf dem jeweiligen Dichteintervall.

    T_values wird ohne Temperaturinterpolation ausgewertet. Luecken bleiben
    als (T,NaN,NaN) erhalten. Der kritische Punkt wird NICHT vorgetaeuscht.
    Gibt es mehrere Extrema (EOS-Artefakt), gilt die von der Saettigung aus
    zuerst erreichte Nullstelle; die Temperaturen werden in
    report['ambiguous_temperatures'] protokolliert.
    """
    if branch not in ('liquid', 'vapor'):
        raise ValueError('branch muss liquid oder vapor sein')
    T = np.asarray(temperatures, dtype=float)
    if T.ndim != 1 or not 2 <= len(T) <= 10000 or not np.isfinite(T).all() or not np.all(np.diff(T) > 0):
        raise ValueError('Ungültiges monotones Temperaturarray')
    if type(n_scan) is not int or not 16 <= n_scan <= 10000 or len(T)*n_scan > 2_000_000:
        raise ValueError('Scan zu gross oder ungueltig')
    if not 0 < residual_tolerance <= 1e-2:
        raise ValueError('Ungültige Toleranz')
    if cp_module is None:
        import CoolProp.CoolProp as cp_module
    if state_factory is None:
        from coolprop_heos_common import _new_state
        state_factory = lambda: _new_state(stoff)
    state, sat = state_factory(), state_factory()
    Tc, pc, rhoc = float(state.T_critical()), float(state.p_critical()), float(state.rhomolar_critical())
    if min(Tc, pc, rhoc) <= 0:
        raise ValueError('Ungültige kritische Stoffdaten')
    out = np.column_stack((T, np.full((len(T), 2), np.nan)))
    report = {'branch': branch, 'requested_points': len(T), 'valid_points': 0,
              'missing_temperatures': [], 'failures': {}, 'ambiguous_temperatures': [],
              'nonpositive_pressures': 0, 'max_residual': None, 'elapsed_seconds': None}
    started = perf_counter()
    residuals = []
    for i, temp in enumerate(T):
        try:
            if not (max(float(sat.Ttriple()), float(sat.Tmin())) < temp < min(Tc, float(sat.Tmax()))):
                raise ValueError('Temperatur ausserhalb des unterkritischen HEOS-Bereichs')
            sat.update(cp_module.QT_INPUTS, 0. if branch == 'liquid' else 1., float(temp))
            rho_sat, p_sat = float(sat.rhomolar()), float(sat.p())
            if p_sat <= 0:
                raise ValueError('Ungültiger Sättigungsdruck')
            if branch == 'liquid':
                lo, hi = rhoc*(1. + 1e-9), rho_sat*(1. - 1e-10)
                expected = 'minimum'
            else:
                lo, hi = rho_sat*(1. + 1e-10), rhoc*(1. - 1e-9)
                expected = 'maximum'
            if not 0 < lo < hi:
                raise ValueError('Zu schmales Suchintervall')
            roots = []
            # Adaptive Verfeinerung nur ohne ausreichend klare Nullstelle.
            for scale in (1, 2, 4):
                rho_grid = np.geomspace(lo, hi, n_scan*scale)
                slopes = np.full(len(rho_grid), np.nan)
                for j, rho in enumerate(rho_grid):
                    try:
                        slopes[j] = homogeneous_point(state, cp_module, temp, 1./rho, branch)[1]
                    except (ValueError, RuntimeError, OverflowError):
                        pass
                brackets = []
                for j in range(len(rho_grid)-1):
                    a,b = slopes[j:j+2]
                    if not np.isfinite([a,b]).all():
                        continue
                    if (a < 0 < b if branch == 'liquid' else a > 0 > b):
                        brackets.append((rho_grid[j], rho_grid[j+1]))
                # Knapp exakt getroffene Gitterwurzel nicht uebergehen.
                for j in range(1,len(rho_grid)-1):
                    if slopes[j] == 0 and np.isfinite(slopes[j-1:j+2]).all() and (
                       slopes[j-1] < 0 < slopes[j+1] if branch == 'liquid'
                       else slopes[j-1] > 0 > slopes[j+1]):
                        brackets.append((rho_grid[j], rho_grid[j]))
                if brackets:
                    roots = brackets
                    break
            if not roots:
                raise ValueError(f'0 {expected}-Nullstellen gefunden')
            if len(roots) > 1:
                # Zusaetzliche EOS-Extrema (Mehrparameter-Artefakt): Die Spinodale ist die erste
                # Nullstelle von dp/drho, die von der Saettigung aus erreicht wird.
                report['ambiguous_temperatures'].append(float(temp))
            a,b = roots[-1] if branch == 'liquid' else roots[0]
            def deriv(rho):
                return homogeneous_point(state, cp_module, temp, 1./float(rho), branch)[1]
            rho_root = float(a if a == b else brentq(deriv, a,b, rtol=1e-12,
                              xtol=max(1e-13, 1e-13*rhoc),maxiter=100))
            p, dpdrho = homogeneous_point(state, cp_module, temp, 1./rho_root, branch)
            residual = abs(dpdrho)*rhoc/pc
            if not (np.isfinite(residual) and residual <= residual_tolerance):
                raise ValueError(f'Residuum {residual} ueber Toleranz')
            if branch == 'liquid':
                if not rhoc < rho_root < rho_sat or not p < p_sat:
                    raise ValueError('Falscher Fluessigkeits-Spinodalast')
            elif not rho_sat < rho_root < rhoc or not p > p_sat:
                raise ValueError('Falscher Dampf-Spinodalast')
            out[i,1:] = 1./rho_root, p
            residuals.append(residual)
        except (ValueError, RuntimeError, OverflowError) as exc:
            report['missing_temperatures'].append(float(temp))
            report['failures'][str(i)] = str(exc)
    # Keine Astwechsel: Monotoniepruefung fuer Volumen in dem Isothermenintervall.
    for i in range(1,len(T)):
        if not (np.isfinite(out[i,1:]).all() and np.isfinite(out[i-1,1:]).all()):
            continue
        if ((branch == 'vapor' and out[i,1] > out[i-1,1]*(1+1e-3)) or
            (branch == 'liquid' and out[i,1] < out[i-1,1]*(1-1e-3))):
            out[i,1:] = np.nan
            report['missing_temperatures'].append(float(T[i]))
            report['failures'][str(i)] = 'Unplausibler Astverlauf'
    good = np.isfinite(out).all(axis=1) & (out[:,1]>0)
    report['valid_points'] = int(good.sum())
    report['nonpositive_pressures'] = int(np.sum(good & (out[:,2]<=0)))
    report['max_residual'] = float(max(residuals)) if residuals else None
    report['success_fraction'] = float(np.mean(good))
    report['elapsed_seconds'] = perf_counter()-started
    return out, report


def _pressure_tolerance(*pressures, pc):
    return 1e-9 * max(pc, *(abs(x) for x in pressures))


def unstable_node_status(state, cp, T, v, window, rhoc, pc):
    """AP19: Klassifiziert einen Knoten der instabilen Fortsetzung.

    Rueckgabe (p, norm, grund); grund is None <=> gueltig. Gueltig heisst:
    (dp/dv)_T > 0 (aequivalent dp/drho < 0) UND p im Spinodaldruckfenster
    [p_spin,l, p_spin,g]. Der EOS-Wert wird nie veraendert, nur klassifiziert.
    """
    p, slope = homogeneous_point(state, cp, T, v, 'vapor')
    norm = slope * rhoc / pc
    if not norm < 0:
        return p, norm, 'mechanisch stabil (dp/dv<=0): zusaetzliches EOS-Extremum'
    if p < window[0]:
        return p, norm, 'p unter Fluessigkeitsspinodaldruck (EOS-Artefakt)'
    if p > window[1]:
        return p, norm, 'p ueber Dampfspinodaldruck (EOS-Artefakt)'
    return p, norm, None


def _bisect_edge(is_ok, v_ok, v_bad, rel=1e-13, iterations=100):
    """Grenze des gueltigen Intervalls; Rueckgabe liegt immer auf der gueltigen Seite."""
    for _ in range(iterations):
        mid = float(np.sqrt(v_ok * v_bad))
        if mid == v_ok or mid == v_bad:
            break
        if is_ok(mid):
            v_ok = mid
        else:
            v_bad = mid
        if abs(v_bad / v_ok - 1.) < rel:
            break
    return v_ok


def trace_unstable_row(state, cp, T, v_l, v_g, p_l, p_g, rhoc, pc, scan_points=256):
    """AP19: Zusammenhaengende instabile Aeste einer Isotherme zwischen den Spinodalen.

    Von der Fluessigkeitsspinodale aus wird der zusammenhaengende Bereich mit
    dp/dv>0 und p in [p_l,p_g] bestimmt, ebenso rueckwaerts von der
    Dampfspinodale. Ist die Isotherme regulaer, ueberdecken beide Aeste das
    gesamte Intervall. Sonst liegt dazwischen ein echtes EOS-Artefakt
    (zusaetzliche Extrema), das als Luecke dokumentiert wird.
    """
    tol = _pressure_tolerance(p_l, p_g, pc=pc)
    window = (p_l - tol, p_g + tol)
    if not (np.isfinite([v_l, v_g, p_l, p_g]).all() and 0 < v_l < v_g and p_l < p_g):
        raise ValueError('Spinodalgrenzen verletzen v_l<v_g bzw. p_l<p_g')

    def status(v):
        try:
            return unstable_node_status(state, cp, T, v, window, rhoc, pc)
        except (ValueError, RuntimeError, OverflowError) as exc:
            return np.nan, np.nan, f'EOS nicht auswertbar: {exc}'

    def is_ok(v):
        return status(v)[2] is None

    grid = np.geomspace(v_l, v_g, scan_points + 2)[1:-1]
    results = [status(v) for v in grid]
    ok = np.array([r[2] is None for r in results])
    p_scan = np.array([r[0] for r in results])
    norm_scan = np.array([r[1] for r in results])
    info = {'T': float(T), 'regular': bool(ok.all()), 'v_a': None, 'v_b': None,
            'left_stop': None, 'right_stop': None,
            'scan_p_min': float(np.nanmin(p_scan)), 'scan_p_max': float(np.nanmax(p_scan)),
            'scan_slope_sign_changes': int(np.count_nonzero(
                np.diff(np.sign(norm_scan[np.isfinite(norm_scan)])) != 0)),
            'window': (float(p_l), float(p_g))}
    if info['regular']:
        return info
    first_bad = int(np.argmax(~ok))
    last_bad = int(len(ok) - 1 - np.argmax((~ok)[::-1]))
    left_ok = grid[first_bad - 1] if first_bad > 0 else v_l
    info['v_a'] = _bisect_edge(is_ok, left_ok, grid[first_bad])
    info['left_stop'] = results[first_bad][2]
    right_ok = grid[last_bad + 1] if last_bad + 1 < len(grid) else v_g
    info['v_b'] = _bisect_edge(is_ok, right_ok, grid[last_bad])
    info['right_stop'] = results[last_bad][2]
    return info


def unstable_row_layout(n_points):
    """Spalten: linker Ast [0,nL), Trennspalte nL, rechter Ast (nL,n_points)."""
    n_left = (n_points - 1) // 2
    return n_left, n_points - 1 - n_left


def _unstable_row_volumes(n_points, v_l, v_g, trace):
    n_left, n_right = unstable_row_layout(n_points)
    if trace['regular']:
        return np.geomspace(v_l, v_g, n_points)
    vs = np.full(n_points, np.nan)
    vs[0], vs[-1] = v_l, v_g
    v_a, v_b = trace['v_a'], trace['v_b']
    if n_left >= 2 and v_a > v_l * (1. + 1e-12):
        vs[:n_left] = np.geomspace(v_l, v_a, n_left)
    if n_right >= 2 and v_b < v_g * (1. - 1e-12):
        vs[n_left + 1:] = np.geomspace(v_b, v_g, n_right)
    return vs


def _enforce_monotone_pressure(row, direction, rel_tol=1e-12):
    """direction=+1 (instabil, dp/dv>0) bzw. -1 (metastabil, dp/dv<0) entlang eines Astes."""
    failures = {}
    finite = np.isfinite(row[:, 2]) & np.isfinite(row[:, 1])
    for j in range(1, len(row)):
        if finite[j] and finite[j - 1]:
            tol = rel_tol * max(abs(row[j, 2]), abs(row[j - 1, 2]), 1.)
            if direction * (row[j, 2] - row[j - 1, 2]) < -tol:
                row[j, 2] = np.nan
                finite[j] = False
                failures[str(j)] = 'Druck nicht monoton in v (EOS-Artefakt zwischen Stuetzstellen)'
    return failures


def build_homogeneous_surfaces(stoff, binodal_liquid, binodal_vapor,
                               spinodal_liquid, spinodal_vapor, *,
                               n_points=20, include_unstable=False,
                               state_factory=None, cp_module=None,
                               slope_tolerance=2e-4, scan_points=256, verify_unstable=True):
    """Homogene EOS-Stuetzstellen: Punktfehler lokal behalten statt ganze Zeilen loeschen.

    Fehlende Randpunkte oder widerspruechliche EOS-Randdruecke sperren die
    gesamte Zeile; isolierte Fehler innerer Punkte nur den jeweiligen Knoten.
    KEINE Druckinterpolation, keine kuenstliche Fortsetzung ueber NaN.

    AP19: Die instabile Flaeche enthaelt nur Knoten mit dp/dv>0 und
    p_spin,l<=p<=p_spin,g, die vom jeweiligen Spinodalast aus zusammenhaengend
    erreichbar sind. Zwischenbereiche mit zusaetzlichen EOS-Extrema bleiben als
    dokumentierte Trennspalte (NaN) leer.
    """
    arrays = [np.asarray(x, dtype=float) for x in (
        binodal_liquid, binodal_vapor, spinodal_liquid, spinodal_vapor)]
    l, g, sl, sg = arrays
    if any(x.ndim != 2 or x.shape != l.shape or x.shape[1] != 3 for x in arrays):
        raise ValueError('Alle Randkurven muessen (n_T,3) mit gleicher Form sein')
    if not all(np.array_equal(x[:, 0], l[:, 0]) for x in arrays):
        raise ValueError('Randkurven benoetigen exakte gleiche Temperaturstuetzstellen')
    if not isinstance(n_points, int) or not 3 <= n_points <= 1000:
        raise ValueError('3..1000 Stuetzpunkte pro Temperaturzeile')
    if len(l) * n_points * 3 > 1_000_000:
        raise ValueError('Oberflaechenbudget von 1 Mio Punkten ueberschritten')
    if cp_module is None:
        import CoolProp.CoolProp as cp_module
    if state_factory is None:
        from coolprop_heos_common import _new_state
        state_factory = lambda: _new_state(stoff)
    state = state_factory()
    pc = float(state.p_critical())
    rhoc = float(state.rhomolar_critical())
    if not np.isfinite([pc, rhoc]).all() or min(pc, rhoc) <= 0:
        raise ValueError('Ungueltige kritische Zustandsdaten')

    keys = ('metastable_liquid_surface', 'metastable_vapor_surface', 'unstable_surface')
    result = {name: np.full((len(l), n_points, 3), np.nan, dtype=float) for name in keys}
    reports = {name: dict(valid_rows=0, boundary_rows=0, partial_rows=0,
                          missing_rows=[], reasons={}, point_reasons={},
                          valid_nodes=0, invalid_nodes=len(l)*n_points,
                          nonpositive_pressures=0, display_clipped_cells=0)
               for name in keys}
    specs = [('metastable_liquid_surface', l, sl, 'liquid', 1),
             ('metastable_vapor_surface', sg, g, 'vapor', 1)]
    if include_unstable:
        specs.append(('unstable_surface', sl, sg, 'vapor', -1))
    reports['unstable_surface'].update(regular_rows=0, gap_rows={})
    for name in keys:
        reports[name]['order_violations'] = {}
    # AP19.2/AP20.1: v_sat,l < v_spin,l < v_spin,g < v_sat,g (und p_spin,l < p_spin,g).
    blocked = {}
    for i in range(len(l)):
        if np.isfinite([l[i, 1], sl[i, 1], sg[i, 1], g[i, 1]]).all():
            if not l[i, 1] < sl[i, 1] < sg[i, 1] < g[i, 1]:
                blocked[i] = 'Grenzordnung v_sat,l<v_spin,l<v_spin,g<v_sat,g verletzt'
            elif np.isfinite([sl[i, 2], sg[i, 2]]).all() and not sl[i, 2] < sg[i, 2]:
                blocked[i] = 'Spinodaldruecke verletzen p_spin,l<p_spin,g'
    n_left, _ = unstable_row_layout(n_points)

    for name, left, right, phase, expected in specs:
        for i, (a, b) in enumerate(zip(left, right)):
            try:
                if i in blocked:
                    reports[name]['order_violations'][str(i)] = blocked[i]
                    raise ValueError(blocked[i])
                if not (np.isfinite(a).all() and np.isfinite(b).all()
                        and a[0] == b[0] == l[i, 0] and 0 < a[1] < b[1]):
                    raise ValueError('Fehlende oder falsch sortierte Randpunkte')
                T = float(a[0])
                # Randdruecke validieren, BEVOR innere Punkte teilweise gespeichert werden.
                for j, point in ((0, a), (-1, b)):
                    branch = ('liquid' if (name == 'metastable_liquid_surface' or
                              (name == 'unstable_surface' and j == 0)) else 'vapor')
                    computed, _ = homogeneous_point(state, cp_module, T, point[1], branch)
                    if abs(computed - point[2]) > max(
                            100., 2e-4 * max(abs(point[2]), abs(computed))):
                        raise ValueError(f'Randdruckabweichung {j}: {computed} vs {point[2]}')
                if name == 'unstable_surface':
                    trace = trace_unstable_row(state, cp_module, T, float(a[1]), float(b[1]),
                                               float(a[2]), float(b[2]), rhoc, pc, scan_points)
                    vs = _unstable_row_volumes(n_points, float(a[1]), float(b[1]), trace)
                else:
                    trace = None
                    vs = np.geomspace(float(a[1]), float(b[1]), n_points)
            except (ValueError, RuntimeError, OverflowError) as exc:
                reports[name]['missing_rows'].append(float(l[i, 0]))
                reports[name]['reasons'][str(i)] = str(exc)
                continue

            row = np.column_stack((np.full(n_points, T), vs,
                                   np.full(n_points, np.nan)))
            row[0] = a
            row[-1] = b
            failures = {}
            if name == 'unstable_surface':
                tol_p = _pressure_tolerance(a[2], b[2], pc=pc)
                window = (float(a[2]) - tol_p, float(b[2]) + tol_p)
            else:
                tol_p = _pressure_tolerance(a[2], b[2], pc=pc) * 1e3
                window = (min(a[2], b[2]) - tol_p, max(a[2], b[2]) + tol_p)
            for j in range(1, n_points - 1):
                if not np.isfinite(vs[j]):
                    continue
                try:
                    if name == 'unstable_surface':
                        p, norm, reason = unstable_node_status(
                            state, cp_module, T, vs[j], window, rhoc, pc)
                        if reason is not None:
                            raise ValueError(reason)
                    else:
                        p, slope = homogeneous_point(state, cp_module, T, vs[j], phase)
                        norm = slope * rhoc / pc
                        if (not np.isfinite(norm) or norm < -slope_tolerance):
                            raise ValueError(f'Abweichende Stabilitaetsklasse (norm={norm})')
                        if not window[0] <= p <= window[1]:
                            raise ValueError(f'Druck {p} ausserhalb Randdruckfenster')
                    row[j, 2] = p
                except (ValueError, RuntimeError, OverflowError) as exc:
                    failures[str(j)] = str(exc)
            failures.update(_enforce_monotone_pressure(
                row, 1 if name == 'unstable_surface' else -1))

            result[name][i] = row
            if trace is not None:
                if trace['regular']:
                    reports[name]['regular_rows'] += 1
                else:
                    reports[name]['gap_rows'][str(i)] = dict(
                        T=T, v_a=trace['v_a'], v_b=trace['v_b'],
                        left_stop=trace['left_stop'], right_stop=trace['right_stop'],
                        scan_p_min=trace['scan_p_min'], scan_p_max=trace['scan_p_max'],
                        slope_sign_changes=trace['scan_slope_sign_changes'],
                        window=trace['window'], gap_column=n_left)
            if failures:
                reports[name]['partial_rows'] += 1
                reports[name]['point_reasons'][str(i)] = failures
            elif trace is not None and not trace['regular']:
                pass
            else:
                reports[name]['valid_rows'] += 1
            reports[name]['boundary_rows'] += 1

        finite = np.isfinite(result[name]).all(axis=2)
        reports[name]['valid_nodes'] = int(np.count_nonzero(finite))
        reports[name]['invalid_nodes'] = int(finite.size - np.count_nonzero(finite))
        reports[name]['nonpositive_pressures'] = int(np.count_nonzero(
            finite & (result[name][..., 2] <= 0)))
        reports[name]['display_clipped_cells'] = reports[name]['nonpositive_pressures']
    if include_unstable and verify_unstable:
        check = validate_unstable_nodes(state, cp_module, result['unstable_surface'],
                                        sl, sg, l=l, g=g)
        reports['unstable_surface']['independent_check'] = dict(
            checked=check['checked'], failures=len(check['failures']),
            passed=check['passed'], examples=check['failures'][:5])
    return result, reports


def critical_refinement_rows(stoff, T_last, n_rows, *, state_factory=None, cp_module=None,
                             min_distance=1e-2):
    """AP20.6: Zusaetzliche Temperaturzeilen Tc-(Tc-T_last)/2^k, k=1..n_rows.

    Rueckgabe (T, liquid(n,3), vapor(n,3)) nur fuer Zeilen mit erfolgreicher
    Saettigungsrechnung (Spalten T,v,p); Tc-T muss mindestens min_distance Kelvin
    betragen, damit HEOS-Saettigungsrechnungen zuverlaessig bleiben.
    """
    if cp_module is None:
        import CoolProp.CoolProp as cp_module
    if state_factory is None:
        from coolprop_heos_common import _new_state
        state_factory = lambda: _new_state(stoff)
    state = state_factory()
    Tc = float(state.T_critical())
    rows_l, rows_g = [], []
    for k in range(1, int(n_rows) + 1):
        T = Tc - (Tc - float(T_last)) / 2.**k
        if not (T_last < T < Tc - min_distance):
            continue
        try:
            state.update(cp_module.QT_INPUTS, 0., T)
            liquid = (T, 1. / float(state.rhomolar()), float(state.p()))
            state.update(cp_module.QT_INPUTS, 1., T)
            vapor = (T, 1. / float(state.rhomolar()), float(state.p()))
        except (ValueError, RuntimeError):
            continue
        if np.isfinite(liquid + vapor).all() and liquid[1] < vapor[1]:
            rows_l.append(liquid)
            rows_g.append(vapor)
    return (np.array([r[0] for r in rows_l]), np.array(rows_l).reshape(-1, 3),
            np.array(rows_g).reshape(-1, 3))


def validate_homogeneous_boundaries(result, l, g, sl, sg, critical):
    """Prueft auch partiell gefuellte Zeilen und meldet fehlende Flaechen ehrlich."""
    l, g, sl, sg = [np.asarray(x, dtype=float) for x in (l, g, sl, sg)]
    checked = {}
    for name, a, b in (('metastable_liquid_surface', l, sl),
                       ('metastable_vapor_surface', sg, g),
                       ('unstable_surface', sl, sg)):
        x = np.asarray(result.get(name, []))
        if x.ndim != 3 or x.shape[-1] != 3 or x.shape[0] != len(a):
            raise ValueError(f'{name}: falsch geformt')
        finite = np.isfinite(x).all(axis=2)
        boundary = finite[:, 0] & finite[:, -1]
        complete = finite.all(axis=1)
        if np.any(boundary):
            if not np.array_equal(x[boundary, 0], a[boundary]) or \
                    not np.array_equal(x[boundary, -1], b[boundary]):
                raise AssertionError(f'{name}: Randknoten nicht exakt identisch')
            for row in x[boundary]:
                volumes = row[np.isfinite(row).all(axis=1), 1]
                if not np.all(np.diff(volumes) > 0):
                    raise AssertionError(f'{name}: Volumenraster nicht streng monoton')
        checked[name] = dict(valid_rows=int(np.count_nonzero(complete)),
                             partial_rows=int(np.count_nonzero(boundary & ~complete)),
                             boundary_rows=int(np.count_nonzero(boundary)),
                             missing_rows=int(np.count_nonzero(~boundary)),
                             valid_nodes=int(np.count_nonzero(finite)),
                             exact_boundary=bool(np.any(boundary)))
    Tc, pc, vc = critical
    if not (np.isfinite([Tc, pc, vc]).all() and min(Tc, pc, vc) > 0):
        raise ValueError('Kritischer Punkt ungueltig')
    return checked


def count_renderable_triangles(points, pressure_axis='log'):
    """Zaehlt moegliche Dreiecke aus 3/4 Gueltigknoten, ohne VTK.

    Die Positivdruckbedingung gilt ausschliesslich fuer die Log-Achse.
    """
    x = np.asarray(points, dtype=float)
    if x.ndim != 3 or x.shape[-1] != 3:
        raise ValueError('Erwarte (n_T,n_v,3)')
    if min(x.shape[:2]) < 2:
        return 0
    if pressure_axis not in PRESSURE_AXES:
        raise ValueError('pressure_axis muss eines von ' + ', '.join(PRESSURE_AXES) + ' sein')
    valid = (np.isfinite(x).all(axis=2) &
             (x[..., 0] > 0) & (x[..., 1] > 0))
    if pressure_axis == 'log':
        valid &= x[..., 2] > 0
    corners = (valid[:-1, :-1].astype(np.int8) +
               valid[1:, :-1] + valid[1:, 1:] + valid[:-1, 1:])
    return int(np.count_nonzero(corners == 3) +
               2 * np.count_nonzero(corners == 4))


# ---------------------------------------------------------------------------
# AP19.4 / AP20.5 / AP20.7: VTK-freie Netzlogik, Topologie- und Abdeckungspruefung
# ---------------------------------------------------------------------------

def node_validity(values, pressure_axis='log'):
    """Knoten, die ueberhaupt gerendert werden duerfen (identisch zum Renderer)."""
    x = np.asarray(values, dtype=float)
    valid = np.isfinite(x).all(axis=2) & (x[..., 0] > 0) & (x[..., 1] > 0)
    if pressure_axis == 'log':
        valid &= x[..., 2] > 0
    return valid


def grid_triangles(valid):
    """Dreiecke aus einem (n_T,n_v)-Knotengitter.

    4 gueltige Ecken -> 2 Dreiecke, 3 Ecken -> 1 Dreieck, sonst keines. Die
    zyklische Eckenreihenfolge bleibt erhalten; dadurch sind alle Dreiecke
    einheitlich orientiert. Rueckgabe: (n,3)-Knotenindizes und (n,2)-Zellindizes.
    """
    valid = np.asarray(valid, dtype=bool)
    if valid.ndim != 2 or min(valid.shape) < 2:
        return np.empty((0, 3), dtype=np.int64), np.empty((0, 2), dtype=np.int64)
    idx = np.arange(valid.size, dtype=np.int64).reshape(valid.shape)
    ids = np.stack((idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]), axis=-1)
    ok = np.stack((valid[:-1, :-1], valid[1:, :-1], valid[1:, 1:], valid[:-1, 1:]), axis=-1)
    count = ok.sum(axis=-1)
    tris, cells = [], []
    ci, cj = np.nonzero(count == 4)
    quad = ids[ci, cj]
    cell = np.column_stack((ci, cj))
    tris += [quad[:, [0, 1, 2]], quad[:, [0, 2, 3]]]
    cells += [cell, cell]
    for missing in range(4):
        ci, cj = np.nonzero((count == 3) & ~ok[..., missing])
        keep = [k for k in range(4) if k != missing]
        tris.append(ids[ci, cj][:, keep])
        cells.append(np.column_stack((ci, cj)))
    return (np.concatenate(tris).astype(np.int64),
            np.concatenate(cells).astype(np.int64))


def drop_degenerate(triangles, cells, xyz_flat):
    """Entfernt Dreiecke mit nicht endlicher oder verschwindender Flaeche."""
    if len(triangles) == 0:
        return triangles, cells, 0
    a, b, c = (xyz_flat[triangles[:, k]] for k in range(3))
    area2 = np.linalg.norm(np.cross(b - a, c - a), axis=1)
    keep = np.isfinite(area2) & (area2 > 0)
    return triangles[keep], cells[keep], int(np.count_nonzero(~keep))


def mesh_topology_report(triangles, cells, valid, xy=None):
    """Komponenten, offene Kanten, Orientierung und Projektionsfaltungen.

    Offene Kanten werden klassifiziert in
      frame        : liegen auf dem Gitterrand (Randkurven der Flaeche),
      gap          : grenzen an eine unvollstaendige/leere Zelle (dokumentierte Luecke),
      unexplained  : grenzen an eine vollstaendig gueltige Nachbarzelle (Nahtfehler).
    xy: (n_nodes,2) Projektion (log v, T); ein Vorzeichenwechsel der Projektionsflaeche
    wuerde eine Faltung/Selbstueberschneidung der Flaeche bedeuten.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    valid = np.asarray(valid, dtype=bool)
    n_nodes = valid.size
    report = dict(triangles=int(len(triangles)), nodes_used=0, components=0,
                  component_sizes=[], open_edges=0, open_edges_frame=0,
                  open_edges_gap=0, open_edges_unexplained=0,
                  nonmanifold_edges=0, orientation_conflicts=0,
                  projected_orientation_flips=0, degenerate_projected=0)
    if len(triangles) == 0:
        return report
    used = np.unique(triangles)
    report['nodes_used'] = int(len(used))
    edges = np.concatenate((triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]))
    owner = np.concatenate([np.arange(len(triangles))] * 3)
    graph = coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])),
                       shape=(n_nodes, n_nodes))
    _, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels[used])
    report['component_sizes'] = sorted((int(s) for s in sizes if s > 0), reverse=True)
    report['components'] = len(report['component_sizes'])
    undirected = np.sort(edges, axis=1)
    key = undirected[:, 0] * n_nodes + undirected[:, 1]
    uniq, first, counts = np.unique(key, return_index=True, return_counts=True)
    report['nonmanifold_edges'] = int(np.count_nonzero(counts > 2))
    directed_key = edges[:, 0] * n_nodes + edges[:, 1]
    _, dcounts = np.unique(directed_key, return_counts=True)
    report['orientation_conflicts'] = int(np.count_nonzero(dcounts > 1))
    nr, nc = valid.shape
    full = np.zeros((nr - 1, nc - 1), dtype=bool)
    if nr > 1 and nc > 1:
        full = (valid[:-1, :-1] & valid[1:, :-1] & valid[1:, 1:] & valid[:-1, 1:])
    for k in np.flatnonzero(counts == 1):
        a, b = edges[first[k]]
        ra, ca = divmod(int(a), nc)
        rb, cb = divmod(int(b), nc)
        ci, cj = cells[owner[first[k]]]
        report['open_edges'] += 1
        if ra == rb:
            nci, ncj = (ci - 1 if ra == ci else ci + 1), cj
        elif ca == cb:
            nci, ncj = ci, (cj - 1 if ca == cj else cj + 1)
        else:
            report['open_edges_gap'] += 1  # Diagonale einer Dreiecks-Teilzelle
            continue
        if not (0 <= nci < nr - 1 and 0 <= ncj < nc - 1):
            report['open_edges_frame'] += 1
        elif full[nci, ncj]:
            report['open_edges_unexplained'] += 1
        else:
            report['open_edges_gap'] += 1
    if xy is not None:
        p = [np.asarray(xy)[triangles[:, k]] for k in range(3)]
        signed = ((p[1][:, 0] - p[0][:, 0]) * (p[2][:, 1] - p[0][:, 1]) -
                  (p[1][:, 1] - p[0][:, 1]) * (p[2][:, 0] - p[0][:, 0]))
        report['degenerate_projected'] = int(np.count_nonzero(signed == 0))
        nonzero = signed[signed != 0]
        if len(nonzero):
            majority = 1. if np.count_nonzero(nonzero > 0) >= len(nonzero) / 2 else -1.
            report['projected_orientation_flips'] = int(np.count_nonzero(nonzero * majority < 0))
    return report


def gap_mask(shape, gap_rows):
    """Strukturelle Trennspalten (EOS-Artefakt) aus dem Report der instabilen Flaeche."""
    mask = np.zeros(shape[:2], dtype=bool)
    for i, info in (gap_rows or {}).items():
        mask[int(i), int(info['gap_column'])] = True
    return mask


def surface_coverage(surface, gap_rows=None, pressure_axis='log'):
    """AP20.7: Knoten- und Zellabdeckung einer Flaeche.

    Zuverlaessig berechenbar sind Zellen zwischen zwei Temperaturzeilen mit
    vorhandenen Randpunkten, ohne strukturelle EOS-Artefakt-Spalte. Davon sind
    'covered' die Zellen mit vier gueltigen Knoten.
    """
    x = np.asarray(surface, dtype=float)
    if x.ndim != 3 or x.shape[-1] != 3:
        raise ValueError('Erwarte (n_T,n_v,3)')
    finite = np.isfinite(x).all(axis=2)
    thermo = finite & (x[..., 0] > 0) & (x[..., 1] > 0)
    valid = node_validity(x, pressure_axis)
    rows = finite[:, 0] & finite[:, -1]
    gaps = gap_mask(x.shape, gap_rows)
    result = dict(nodes_total=int(finite.size), nodes_valid=int(np.count_nonzero(thermo)),
                  nodes_nan=int(np.count_nonzero(~finite)),
                  nodes_nonrenderable=int(np.count_nonzero(thermo & ~valid)),
                  nodes_structural_gap=int(np.count_nonzero(gaps)),
                  rows_available=int(np.count_nonzero(rows)),
                  rows_missing=int(np.count_nonzero(~rows)))
    if min(x.shape[:2]) < 2:
        result.update(cells_evaluable=0, cells_covered=0, cells_partial=0,
                      cells_excluded_gap=0, coverage=None, coverage_rendered=None)
        return result
    row_pair = rows[:-1] & rows[1:]
    corner_gap = (gaps[:-1, :-1] | gaps[1:, :-1] | gaps[1:, 1:] | gaps[:-1, 1:])
    evaluable = row_pair[:, None] & ~corner_gap

    def quads(mask):
        four = mask[:-1, :-1] & mask[1:, :-1] & mask[1:, 1:] & mask[:-1, 1:]
        count = (mask[:-1, :-1].astype(int) + mask[1:, :-1] + mask[1:, 1:] + mask[:-1, 1:])
        return four, count

    four, count = quads(thermo)
    covered = evaluable & four
    four_r, _ = quads(valid)
    result.update(cells_evaluable=int(np.count_nonzero(evaluable)),
                  cells_covered=int(np.count_nonzero(covered)),
                  cells_partial=int(np.count_nonzero(evaluable & (count == 3))),
                  cells_excluded_gap=int(np.count_nonzero(row_pair[:, None] & corner_gap)),
                  cells_covered_rendered=int(np.count_nonzero(evaluable & four_r)))
    n_eval = result['cells_evaluable']
    result['coverage'] = result['cells_covered'] / n_eval if n_eval else None
    result['coverage_rendered'] = result['cells_covered_rendered'] / n_eval if n_eval else None
    return result


def combined_surface_report(result, pressure_axis='asinh', unstable_gap_rows=None):
    """AP20.3/20.5: Die drei Teilflaechen als ein Netz mit exakt gemeinsamen Knoten.

    Knoten gelten als identisch, wenn (T,v,p) bitgleich sind. Geprueft werden
    Komponenten, offene Kanten nach dem Verschmelzen und die Naehte
    Fluessigkeit|instabil und instabil|Dampf.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    names = ('metastable_liquid_surface', 'unstable_surface', 'metastable_vapor_surface')
    per_surface, coords, tri_global, cursor = {}, [], [], 0
    for name in names:
        x = np.asarray(result[name], dtype=float)
        valid = node_validity(x, pressure_axis)
        tri, _cells = grid_triangles(valid)
        used = np.unique(tri) if len(tri) else np.empty(0, dtype=np.int64)
        remap = np.full(valid.size, -1, dtype=np.int64)
        remap[used] = np.arange(len(used)) + cursor
        coords.append(x.reshape(-1, 3)[used])
        per_surface[name] = (tri, remap, x.shape)
        tri_global.append(remap[tri] if len(tri) else tri)
        cursor += len(used)
    all_coords = np.concatenate(coords) if coords else np.empty((0, 3))
    summary = dict(triangles=int(sum(len(t) for t in tri_global)), nodes_unmerged=int(len(all_coords)))
    if not len(all_coords) or summary['triangles'] == 0:
        summary.update(nodes_merged=0, components=0, open_edges=0, seams={})
        return summary
    uniq, inverse = np.unique(all_coords, axis=0, return_inverse=True)
    inverse = inverse.reshape(-1)
    summary['nodes_merged'] = int(len(uniq))
    shared = np.bincount(inverse, minlength=len(uniq)) >= 2
    summary['shared_nodes'] = int(np.count_nonzero(shared))
    triangles = np.concatenate([t for t in tri_global if len(t)])
    triangles = inverse[triangles]
    edges = np.concatenate((triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]))
    n = len(uniq)
    graph = coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels[np.unique(triangles)])
    summary['component_sizes'] = sorted((int(s) for s in sizes if s > 0), reverse=True)
    summary['components'] = len(summary['component_sizes'])
    undirected = np.sort(edges, axis=1)
    keys = undirected[:, 0] * n + undirected[:, 1]
    uniq_keys, first, counts = np.unique(keys, return_index=True, return_counts=True)
    open_idx = first[counts == 1]
    summary['open_edges'] = int(len(open_idx))
    summary['open_edges_between_shared_nodes'] = int(np.count_nonzero(
        shared[undirected[open_idx, 0]] & shared[undirected[open_idx, 1]]))
    summary['nonmanifold_edges'] = int(np.count_nonzero(counts > 2))

    # Zusammenhaengende, durchgaengig berechenbare Temperaturintervalle (alle drei Randzeilen da)
    rows_ok = np.ones(np.asarray(result[names[0]]).shape[0], dtype=bool)
    for name in names:
        x = np.asarray(result[name], dtype=float)
        rows_ok &= np.isfinite(x[:, 0]).all(axis=1) & np.isfinite(x[:, -1]).all(axis=1)
    runs, start = [], None
    for i, ok in enumerate(list(rows_ok) + [False]):
        if ok and start is None:
            start = i
        elif not ok and start is not None:
            if i - start >= 2:
                runs.append((start, i - 1))
            start = None
    gap_set = {int(k) for k in (unstable_gap_rows or {})}
    expected = 0
    for a, b in runs:
        regular = any(i not in gap_set for i in range(a, b + 1))
        expected += 1 if regular else 2
    summary['temperature_intervals'] = [(int(a), int(b)) for a, b in runs]
    summary['components_expected'] = expected

    def column_edges(name, column):
        tri, _remap, shape = per_surface[name]
        nc = shape[1]
        found = set()
        for t in tri:
            nodes = [(int(k) // nc, int(k) % nc) for k in t]
            on = [r for r, c in nodes if c == column]
            if len(on) == 2 and abs(on[0] - on[1]) == 1:
                found.add(min(on))
        return found

    seams = {}
    for left, right in (('metastable_liquid_surface', 'unstable_surface'),
                        ('unstable_surface', 'metastable_vapor_surface')):
        xl, xr = np.asarray(result[left], dtype=float), np.asarray(result[right], dtype=float)
        cl = per_surface[left][2][1] - 1
        a, b = column_edges(left, cl), column_edges(right, 0)
        fin_l, fin_r = np.isfinite(xl[:, cl]).all(axis=1), np.isfinite(xr[:, 0]).all(axis=1)
        both = fin_l & fin_r
        expected = {i for i in range(len(both) - 1) if both[i] and both[i + 1]}
        identical = bool(np.array_equal(xl[both, cl], xr[both, 0]))
        seams[f'{left}|{right}'] = dict(
            shared_rows=int(both.sum()), expected_edges=len(expected), identical_nodes=identical,
            missing_left=sorted(expected - a), missing_right=sorted(expected - b),
            # Kanten nur einer Seite ohne gemeinsame Randknoten: andere Flaeche hat dort keine Zeile
            one_sided=sorted((a ^ b) - expected))
    summary['seams'] = seams
    return summary


def flag_triangles(values, triangles, windows):
    """AP19.4: Dreiecke mit Knoten ausserhalb des Spinodaldruckfensters protokollieren.

    windows[i] = (p_spin,l, p_spin,g) der Temperaturzeile i. Rueckgabe: Liste mit den
    urspruenglichen (T,v,p)-Koordinaten der betroffenen Knoten.
    """
    x = np.asarray(values, dtype=float)
    nc = x.shape[1]
    flat = x.reshape(-1, 3)
    flagged = []
    for t, tri in enumerate(np.asarray(triangles)):
        bad = []
        for k in tri:
            row = int(k) // nc
            lo, hi = windows[row]
            p = flat[k, 2]
            if not lo - 1e-9 * abs(lo) - 1. <= p <= hi + 1e-9 * abs(hi) + 1.:
                bad.append((int(k) // nc, int(k) % nc, *map(float, flat[k])))
        if bad:
            flagged.append(dict(triangle=t, nodes=[tuple(map(float, flat[k])) for k in tri],
                                nodes_outside=bad))
    return flagged


def validate_unstable_nodes(state, cp, surface, sl, sg, *, l=None, g=None, fd_rel=1e-5):
    """AP19.2/19.3: Unabhaengige Nachpruefung jedes akzeptierten instabilen Knotens.

    Je Knoten: v_spin,l<=v<=v_spin,g, p im Spinodaldruckfenster, (dp/dv)_T>0 ueber
    zentrale Differenzen der p(v)-Werte (nicht ueber first_partial_deriv) und
    identischer Druck bei expliziter Fluessig-/Gasphasenvorgabe. Optional
    v_sat,l<v_spin,l und v_spin,g<v_sat,g.
    """
    x = np.asarray(surface, dtype=float)
    sl, sg = np.asarray(sl, dtype=float), np.asarray(sg, dtype=float)
    pc = float(state.p_critical())
    failures, checked = [], 0
    for i in range(x.shape[0]):
        finite_row = np.isfinite(x[i]).all(axis=1)
        for j in range(x.shape[1]):
            T, v, p = x[i, j]
            if not np.isfinite([T, v, p]).all():
                continue
            checked += 1
            reasons = []
            if not (sl[i, 1] <= v <= sg[i, 1]):
                reasons.append('ausserhalb der Spinodalgrenzen')
            tol = _pressure_tolerance(sl[i, 2], sg[i, 2], pc=pc)
            if not (sl[i, 2] - tol <= p <= sg[i, 2] + tol):
                reasons.append('p ausserhalb [p_spin,l; p_spin,g]')
            if l is not None and not l[i, 1] < sl[i, 1]:
                reasons.append('v_sat,l<v_spin,l verletzt')
            if g is not None and not sg[i, 1] < g[i, 1]:
                reasons.append('v_spin,g<v_sat,g verletzt')
            has_prev = j > 0 and finite_row[j - 1]
            has_next = j + 1 < x.shape[1] and finite_row[j + 1]
            try:
                p_liq = homogeneous_point(state, cp, T, v, 'liquid')[0]
                p_gas = homogeneous_point(state, cp, T, v, 'vapor')[0]
                if abs(p_liq - p) > 1e-6 * max(abs(p), pc) or abs(p_gas - p) > 1e-6 * max(abs(p), pc):
                    reasons.append('Druck weicht von der EOS-Neuberechnung ab (Phasenvorgabe)')
                # Differenzenquotient, nach innen gerichtet an Astenden (dort ist dp/dv=0).
                h_f = min(fd_rel * v, 0.25 * (x[i, j + 1, 1] - v)) if has_next else fd_rel * v
                h_b = min(fd_rel * v, 0.25 * (v - x[i, j - 1, 1])) if has_prev else fd_rel * v
                if has_prev and has_next:
                    p_hi = homogeneous_point(state, cp, T, v + h_f, 'vapor')[0]
                    p_lo = homogeneous_point(state, cp, T, v - h_b, 'vapor')[0]
                    dp_dv = (p_hi - p_lo) / (h_f + h_b)
                elif has_next or not has_prev:
                    dp_dv = (homogeneous_point(state, cp, T, v + h_f, 'vapor')[0] - p) / h_f
                else:
                    dp_dv = (p - homogeneous_point(state, cp, T, v - h_b, 'vapor')[0]) / h_b
                if not dp_dv > 0:
                    reasons.append(f'(dp/dv)_T={dp_dv:.3e} <= 0')
            except (ValueError, RuntimeError, OverflowError) as exc:
                reasons.append(f'EOS-Pruefung fehlgeschlagen: {exc}')
            if reasons:
                failures.append(dict(i=i, j=j, T=float(T), v=float(v), p=float(p), reasons=reasons))
    return dict(checked=checked, failures=failures, passed=not failures)


def critical_approach_report(l, sl, sg, g, critical, last_rows=6):
    """AP20.6: Verhalten der vier Grenzkurven nahe Tc (nur Diagnose, keine Schliessung).

    Prueft, ob in den letzten gueltigen Zeilen alle Intervallbreiten monoton
    schrumpfen, die Grenzordnung gilt und die Zeilen auf die Kritische-Punkt-Naeherung
    zulaufen. closure_supported ist nur dann True, wenn alle Bedingungen gelten.
    """
    Tc, pc, vc = critical
    l, sl, sg, g = [np.asarray(a, dtype=float) for a in (l, sl, sg, g)]
    ok = np.isfinite(np.column_stack((l[:, 1], sl[:, 1], sg[:, 1], g[:, 1]))).all(axis=1)
    idx = np.flatnonzero(ok)
    idx = idx[-last_rows:]
    report = dict(rows=[int(i) for i in idx], closure_supported=False, reasons=[])
    if len(idx) < 3:
        report['reasons'].append('weniger als 3 vollstaendige Randzeilen nahe Tc')
        return report
    T = l[idx, 0]
    widths = {'binodal': g[idx, 1] - l[idx, 1], 'spinodal': sg[idx, 1] - sl[idx, 1],
              'metastable_liquid': sl[idx, 1] - l[idx, 1], 'metastable_vapor': g[idx, 1] - sg[idx, 1]}
    report['T_last'] = float(T[-1])
    report['reduced_distance_T'] = float((Tc - T[-1]) / Tc)
    report['widths_over_vc'] = {k: [float(w / vc) for w in v] for k, v in widths.items()}
    if not np.all(np.diff(T) > 0):
        report['reasons'].append('Temperaturen nicht monoton')
    for key, w in widths.items():
        if not (np.all(w > 0) and np.all(np.diff(w) < 0)):
            report['reasons'].append(f'{key}: Breite schrumpft nicht streng monoton')
    inside = (np.all(l[idx, 1] < vc) and np.all(g[idx, 1] > vc) and
              np.all(sl[idx, 1] < vc) and np.all(sg[idx, 1] > vc))
    if not inside:
        report['reasons'].append('v_c liegt nicht zwischen den Grenzkurven')
    pressures = np.column_stack((l[idx, 2], sl[idx, 2], sg[idx, 2]))
    if not (np.all(sl[idx, 2] < g[idx, 2]) and np.all(g[idx, 2] < sg[idx, 2])):
        report['reasons'].append('Druckordnung p_spin,l<p_sat<p_spin,g verletzt')
    if np.any(~np.isfinite(pressures)) or not np.all(np.abs(pressures[-1] - pc) < 0.1 * pc):
        report['reasons'].append('letzte Zeile nicht innerhalb 10 % von p_c')
    report['closure_supported'] = not report['reasons']
    return report
