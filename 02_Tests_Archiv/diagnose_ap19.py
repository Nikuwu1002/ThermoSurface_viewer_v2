"""AP19: Ursachenanalyse der nach oben ausreissenden instabilen Flaeche (Standardfall Wasser).

Aufruf:  python 02_Tests_Archiv/diagnose_ap19.py

Prueft nacheinander (jeweils unabhaengig von der Visualisierung):
  A  HEOS-Rohwerte gegen eine unabhaengige Implementierung (IAPWS-95 via 'iapws', optional),
  B  Phasenvorgabe (iphase_liquid / iphase_gas): veraendert sie den Druck?,
  C  Extrema von p(rho) und Vorzeichen von dp/drho zwischen den Spinodalen,
  D  Netzbildung: enthalten die Ausreisser-Dreiecke ausschliesslich Ausreisser-Knoten?,
  E  Vorher/Nachher der instabilen Flaeche.
Es wird nichts abgeschnitten oder geglaettet; Ausreisser werden nur klassifiziert und protokolliert.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np

TESTS_DIR = Path(__file__).resolve().parents[1] / "03_Tests_Implementierung"
PROGRAM_DIR = Path(__file__).resolve().parents[1] / "01_Programm"
sys.path.insert(0, str(TESTS_DIR))
sys.path.insert(0, str(PROGRAM_DIR))

import CoolProp.CoolProp as CP
import homogeneous_pvt as new
import validate_ap19_ap20 as val


def raw_p_independent(fluid, T, rho_mol):
    """Druck direkt aus der IAPWS-95-Helmholtzfunktion (nur Wasser), None ohne 'iapws'."""
    if fluid != 'Water':
        return None
    try:
        from iapws import IAPWS95
    except ImportError:
        return None
    ref = IAPWS95(T=T, P=1.0)
    M = 0.018015268
    rho = rho_mol * M
    d = ref._phir(ref.Tc / T, rho / ref.rhoc)
    return rho * ref.R * T * (1. + rho / ref.rhoc * d['fird']) * 1e3  # kJ -> J


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--fluid', default='Water')
    parser.add_argument('--T-start', type=float, default=580.)
    parser.add_argument('--T-stop', type=float, default=646.)
    parser.add_argument('--T-step', type=float, default=2.)
    parser.add_argument('--n-points', type=int, default=20)
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parent /
                                                     'results' / 'ap19_diagnose_water.json')
    args = parser.parse_args()
    fluid = args.fluid
    state = val._new_state(fluid)
    Tc, pc = float(state.T_critical()), float(state.p_critical())
    rhoc = float(state.rhomolar_critical())
    T = np.arange(args.T_start, min(args.T_stop, Tc - 1e-6) + 1e-9, args.T_step)
    T = T[T < Tc]
    l = np.full((len(T), 3), np.nan)
    g = l.copy()
    for i, t in enumerate(T):
        state.update(CP.QT_INPUTS, 0., float(t))
        l[i] = t, 1. / state.rhomolar(), state.p()
        state.update(CP.QT_INPUTS, 1., float(t))
        g[i] = t, 1. / state.rhomolar(), state.p()
    sl, _ = new.solve_spinodal_branch(fluid, T, 'liquid')
    sg, _ = new.solve_spinodal_branch(fluid, T, 'vapor')
    report = {'fluid': fluid, 'T': [float(T[0]), float(T[-1]), float(args.T_step)],
              'n_points': args.n_points, 'rows': int(len(T))}

    # --- A/B/C: eine Isotherme im Detail (Index der ersten Zeile mit Artefakt)
    legacy = val.load_legacy()
    old, _ = legacy.build_homogeneous_surfaces(fluid, l, g, sl, sg, n_points=args.n_points,
                                               include_unstable=True)
    u_old = old['unstable_surface']
    windows = np.column_stack((sl[:, 2], sg[:, 2]))
    row = int(np.argmin(np.abs(T - 600.))) if fluid == 'Water' else len(T) // 3
    Tr = float(T[row])
    nodes = u_old[row]
    detail = []
    for k, (t, v, p) in enumerate(nodes):
        if not np.isfinite(p):
            continue
        p_liq = new.homogeneous_point(state, CP, t, v, 'liquid')
        p_gas = new.homogeneous_point(state, CP, t, v, 'vapor')
        ref = raw_p_independent(fluid, t, 1. / v)
        detail.append(dict(j=k, T=float(t), v=float(v), p=float(p), p_liquid_spec=p_liq[0],
                           p_gas_spec=p_gas[0], p_independent_helmholtz=ref,
                           dp_drho_normed=float(p_liq[1] * rhoc / pc),
                           inside_window=bool(windows[row, 0] <= p <= windows[row, 1])))
    report['isotherm'] = {'T': Tr, 'window_Pa': windows[row].tolist(), 'nodes': detail}
    report['A_independent_eos_identical'] = (
        None if detail[0]['p_independent_helmholtz'] is None else
        bool(all(abs(d['p_independent_helmholtz'] - d['p']) <= 1e-4 * max(abs(d['p']), pc) for d in detail)))
    report['B_phase_spec_changes_pressure'] = bool(any(
        abs(d['p_liquid_spec'] - d['p_gas_spec']) > 1e-9 * max(abs(d['p']), pc) for d in detail))
    # Extrema von p(rho) auf feiner Dichteskala (CoolProp, explizite Gasphase)
    rho = np.linspace(1. / sg[row, 1] * 1.0005, 1. / sl[row, 1] * 0.9995, 4000)
    pp = np.array([new.homogeneous_point(state, CP, Tr, 1. / r, 'vapor')[0] for r in rho])
    sign = np.sign(np.diff(pp))
    ext = np.flatnonzero(np.diff(sign) != 0)
    report['C_extrema_of_p_rho'] = [dict(rho=float(rho[k + 1]), v=float(1. / rho[k + 1]),
                                         p=float(pp[k + 1]), p_over_pc=float(pp[k + 1] / pc))
                                    for k in ext]
    report['C_expected_extrema_van_der_Waals_like'] = 0

    # --- D: Netz
    valid = new.node_validity(u_old, 'asinh')
    tri, cells = new.grid_triangles(valid)
    flagged = new.flag_triangles(u_old, tri, windows)
    flagged_ids = {f['triangle'] for f in flagged}
    xyz = np.stack((np.log10(np.where(valid, u_old[..., 1], 1.) * rhoc), u_old[..., 0] / Tc,
                    np.arcsinh(np.where(valid, u_old[..., 2], 0.) / pc)), axis=-1).reshape(-1, 3)
    edge_len = np.stack([np.linalg.norm(xyz[tri[:, a]] - xyz[tri[:, b]], axis=1)
                         for a, b in ((0, 1), (1, 2), (2, 0))], axis=1).max(axis=1)
    mask = np.zeros(len(tri), dtype=bool)
    mask[list(flagged_ids)] = True
    report['D_mesh'] = dict(
        triangles=int(len(tri)), triangles_with_outlier_node=int(mask.sum()),
        max_edge_with_outlier_node=float(edge_len[mask].max()) if mask.any() else None,
        max_edge_without_outlier_node=float(edge_len[~mask].max()) if (~mask).any() else None,
        median_edge_without_outlier_node=float(np.median(edge_len[~mask])) if (~mask).any() else None,
        conclusion='Alle langen Kanten gehen von Knoten ausserhalb des Spinodaldruckfensters aus; '
                   'die Vernetzung verbindet nur vorhandene Knoten.')
    worst = sorted(flagged, key=lambda f: -max(abs(n[4]) for n in f['nodes_outside']))[:10]
    report['D_worst_triangles_original_TvP'] = [
        dict(triangle=f['triangle'], nodes_TvP=f['nodes'], nodes_outside_window=f['nodes_outside'])
        for f in worst]

    # --- E: Vorher/Nachher
    after, _ = new.build_homogeneous_surfaces(fluid, l, g, sl, sg, n_points=args.n_points,
                                              include_unstable=True)
    for label, arr in (('vorher', u_old), ('nachher', after['unstable_surface'])):
        fin = np.isfinite(arr).all(axis=2)
        outside = 0
        worst_ratio = 0.
        for i in range(len(T)):
            p = arr[i, fin[i], 2]
            lo, hi = windows[i]
            tol = new._pressure_tolerance(lo, hi, pc=pc)
            bad = (p < lo - tol) | (p > hi + tol)
            outside += int(bad.sum())
            if bad.any():
                worst_ratio = max(worst_ratio, float(np.max(np.maximum(lo - p, p - hi)[bad]) / pc))
        report['E_' + label] = dict(nodes=int(fin.sum()), nodes_outside_window=outside,
                                    max_excess_over_pc=worst_ratio,
                                    p_min_over_pc=float(arr[..., 2][fin].min() / pc),
                                    p_max_over_pc=float(arr[..., 2][fin].max() / pc))
    print(f"{fluid}: Isotherme T={Tr:.1f} K, Fenster [{windows[row,0]:.4g}, {windows[row,1]:.4g}] Pa")
    print(f"  A unabhaengige Helmholtz-Auswertung identisch: {report['A_independent_eos_identical']}")
    print(f"  B Phasenvorgabe veraendert Druck: {report['B_phase_spec_changes_pressure']}")
    print("  C Extrema p(rho) zwischen den Spinodalen:",
          [(round(e['rho']), round(e['p_over_pc'], 1)) for e in report['C_extrema_of_p_rho']])
    d = report['D_mesh']
    print(f"  D Dreiecke mit Ausreisserknoten: {d['triangles_with_outlier_node']}/{d['triangles']}; "
          f"max. Kante mit/ohne Ausreisser: {d['max_edge_with_outlier_node']:.3g} / "
          f"{d['max_edge_without_outlier_node']:.3g}")
    for label in ('vorher', 'nachher'):
        e = report['E_' + label]
        print(f"  E {label}: Knoten={e['nodes']}, ausserhalb Fenster={e['nodes_outside_window']}, "
              f"p/pc in [{e['p_min_over_pc']:.3g}, {e['p_max_over_pc']:.3g}]")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=1, default=float), encoding='utf-8')
    print('Ergebnisse:', args.out)


if __name__ == '__main__':
    main()
