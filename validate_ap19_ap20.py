"""AP19/AP20: reproduzierbare Vorher-Nachher-Validierung der homogenen Flaechen (ohne GUI/VTK).

Beispiele:
    python 03_Tests_Implementierung/validate_ap19_ap20.py --fluids Water CO2 Argon --out 03_Tests_Implementierung/results/ap19_ap20_validation.json
    python 03_Tests_Implementierung/validate_ap19_ap20.py --fluids Water --quick

Je Fluid und Aufloesung (n_T, n_points) werden die Spinodalen einmal berechnet und die
drei homogenen Flaechen mit dem alten (Legacy) und dem korrigierten Aufbau erzeugt.
Alle Kennzahlen stammen aus den Rohdaten; es wird nichts abgeschnitten oder geglaettet.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np

PROGRAM_DIR = Path(__file__).resolve().parents[1] / "01_Programm"
sys.path.insert(0, str(PROGRAM_DIR))

import CoolProp.CoolProp as CP
from coolprop_heos_common import _new_state
import homogeneous_pvt as new

ROOT = Path(__file__).resolve().parent
RESOLUTIONS = {'niedrig': (14, 10), 'mittel': (30, 20), 'hoch': (60, 40)}
SURFACES = ('metastable_liquid_surface', 'unstable_surface', 'metastable_vapor_surface')


def load_legacy():
    path = ROOT / 'regression' / 'ap19_vorher' / 'homogeneous_pvt_legacy.py'
    spec = importlib.util.spec_from_file_location('homogeneous_pvt_legacy', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def boundary_curves(fluid, n_T, low=0.80):
    """Gemeinsames Temperaturraster (nach Tc verdichtet) und beide Binodaeste."""
    state = _new_state(fluid)
    Tc = float(state.T_critical())
    T_low = max(low * Tc, 1.001 * float(state.Ttriple()), 1.001 * float(state.Tmin()))
    T = Tc - (Tc - T_low) * np.geomspace(1., 0.002, n_T)
    T = np.unique(T[(T > T_low) & (T < Tc)])
    liquid = np.full((len(T), 3), np.nan)
    vapor = liquid.copy()
    for i, t in enumerate(T):
        try:
            state.update(CP.QT_INPUTS, 0., float(t))
            liquid[i] = t, 1. / state.rhomolar(), state.p()
            state.update(CP.QT_INPUTS, 1., float(t))
            vapor[i] = t, 1. / state.rhomolar(), state.p()
        except (ValueError, RuntimeError):
            liquid[i] = vapor[i] = np.nan
            liquid[i, 0] = vapor[i, 0] = t
    return T, liquid, vapor, state


def surface_metrics(name, arr, critical, reports, windows, axis='log'):
    pc = critical[1]
    finite = np.isfinite(arr).all(axis=2)
    valid = new.node_validity(arr, axis)
    tri, cells = new.grid_triangles(valid)
    out = dict(nodes_valid=int(valid.sum()), nodes_total=int(valid.size), triangles=int(len(tri)))
    if finite.any():
        p = arr[..., 2][finite]
        out['p_min_over_pc'] = float(p.min() / pc)
        out['p_max_over_pc'] = float(p.max() / pc)
    if name == 'unstable_surface' and windows is not None:
        outside = 0
        worst = 0.
        for i in range(arr.shape[0]):
            ok = finite[i]
            if not ok.any() or not np.isfinite(windows[i]).all():
                continue
            lo, hi = windows[i]
            tol = 1e-9 * max(abs(lo), abs(hi), pc)
            pi = arr[i, ok, 2]
            bad = (pi < lo - tol) | (pi > hi + tol)
            outside += int(bad.sum())
            if bad.any():
                worst = max(worst, float(np.max(np.maximum(lo - pi, pi - hi)[bad])) / pc)
        out['nodes_outside_spinodal_pressure_window'] = outside
        out['worst_excess_over_pc'] = worst
        flagged = new.flag_triangles(arr, tri, windows) if len(tri) else []
        out['triangles_with_node_outside_window'] = len(flagged)
        out['flagged_examples'] = flagged[:3]
    return out


def evaluate(module, fluid, T, l, g, sl, sg, n_points, critical, label):
    result, reports = module.build_homogeneous_surfaces(
        fluid, l, g, sl, sg, n_points=n_points, include_unstable=True)
    windows = np.column_stack((sl[:, 2], sg[:, 2]))
    data = {'surfaces': {}}
    for name in SURFACES:
        entry = surface_metrics(name, result[name], critical, reports[name], windows)
        gaps = reports[name].get('gap_rows')
        entry['coverage'] = new.surface_coverage(result[name], gaps)
        valid = new.node_validity(result[name], 'asinh')
        tri, cells = new.grid_triangles(valid)
        from_xyz = np.stack((np.log10(np.where(valid, result[name][..., 1], 1) / critical[2]),
                             result[name][..., 0] / critical[0]), axis=-1).reshape(-1, 2)
        entry['topology'] = new.mesh_topology_report(tri, cells, valid, xy=from_xyz)
        entry['reports'] = {k: (v if not isinstance(v, (dict, list)) or len(str(v)) < 400 else '...')
                            for k, v in reports[name].items() if k not in ('reasons', 'point_reasons')}
        data['surfaces'][name] = entry
    gap_rows = reports['unstable_surface'].get('gap_rows')
    # 'asinh' = alle thermodynamisch gueltigen Knoten (auch p<=0); 'log' = Darstellung mit Log-Achse.
    data['combined'] = new.combined_surface_report(result, 'asinh', gap_rows)
    data['combined_log'] = new.combined_surface_report(result, 'log', gap_rows)
    data['_arrays'] = result
    data['_reports'] = reports
    return data


def run(fluid, resolution, legacy, quick=False):
    n_T, n_points = RESOLUTIONS[resolution]
    T, l, g, state = boundary_curves(fluid, n_T)
    critical = (float(state.T_critical()), float(state.p_critical()), 1. / float(state.rhomolar_critical()))
    sl, rl = new.solve_spinodal_branch(fluid, T, 'liquid')
    sg, rg = new.solve_spinodal_branch(fluid, T, 'vapor')
    available = np.isfinite(np.column_stack((l[:, 1], sl[:, 1], sg[:, 1], g[:, 1]))).all(axis=1)
    row = dict(fluid=fluid, resolution=resolution, n_T=int(len(T)), n_points=n_points,
               rows_all_four_boundaries=int(available.sum()),
               missing_boundary_temperatures=[float(t) for t in T[~available]],
               spinodal_valid=dict(liquid=rl['valid_points'], vapor=rg['valid_points']))
    after = evaluate(new, fluid, T, l, g, sl, sg, n_points, critical, 'nachher')
    row['nachher'] = {k: v for k, v in after.items() if not k.startswith('_')}
    row['nachher']['unstable_rows'] = dict(
        regular=after['_reports']['unstable_surface']['regular_rows'],
        eos_artifact=len(after['_reports']['unstable_surface']['gap_rows']))
    check = new.validate_unstable_nodes(state, CP, after['_arrays']['unstable_surface'], sl, sg, l=l, g=g)
    row['nachher']['unstable_independent_check'] = dict(
        checked=check['checked'], failures=len(check['failures']), first=check['failures'][:3])
    row['nachher']['boundary_check'] = new.validate_homogeneous_boundaries(
        after['_arrays'], l, g, sl, sg, critical)
    row['critical_approach'] = new.critical_approach_report(l, sl, sg, g, critical)
    if legacy is not None:
        before = evaluate(legacy, fluid, T, l, g, sl, sg, n_points, critical, 'vorher')
        row['vorher'] = {k: v for k, v in before.items() if not k.startswith('_')}
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--fluids', nargs='+', default=['Water', 'CO2', 'Argon'])
    parser.add_argument('--resolutions', nargs='+', default=list(RESOLUTIONS), choices=list(RESOLUTIONS))
    parser.add_argument('--no-legacy', action='store_true')
    parser.add_argument('--quick', action='store_true', help='nur Aufloesung mittel')
    parser.add_argument('--out', type=Path, default=None)
    args = parser.parse_args()
    resolutions = ['mittel'] if args.quick else args.resolutions
    legacy = None if args.no_legacy else load_legacy()
    rows = []
    for fluid in args.fluids:
        for resolution in resolutions:
            row = run(fluid, resolution, legacy)
            rows.append(row)
            u_new = row['nachher']['surfaces']['unstable_surface']
            line = (f"{fluid:8s} {resolution:7s} T-Zeilen={row['rows_all_four_boundaries']}/{row['n_T']} "
                    f"| instabil NACHHER: Knoten={u_new['nodes_valid']}, ausserhalb Fenster="
                    f"{u_new.get('nodes_outside_spinodal_pressure_window')}, "
                    f"regulaer={row['nachher']['unstable_rows']['regular']}, "
                    f"EOS-Artefaktzeilen={row['nachher']['unstable_rows']['eos_artifact']}, "
                    f"Pruefung={row['nachher']['unstable_independent_check']['failures']} Fehler")
            if 'vorher' in row:
                u_old = row['vorher']['surfaces']['unstable_surface']
                line += (f" | VORHER: Knoten={u_old['nodes_valid']}, ausserhalb Fenster="
                         f"{u_old.get('nodes_outside_spinodal_pressure_window')}, "
                         f"max|p|/pc={max(abs(u_old.get('p_min_over_pc', 0)), abs(u_old.get('p_max_over_pc', 0))):.3g}")
            print(line)
            for name in SURFACES:
                s = row['nachher']['surfaces'][name]
                cov = s['coverage']['coverage']
                cov_r = s['coverage']['coverage_rendered']
                top = s['topology']
                print(f"    {name:27s} Dreiecke={s['triangles']:5d} Abdeckung(thermo)="
                      f"{'n/a' if cov is None else f'{cov:.4f}'} (log-Achse "
                      f"{'n/a' if cov_r is None else f'{cov_r:.4f}'}) Komponenten={top['components']} "
                      f"offene Kanten(unerklaert)={top['open_edges_unexplained']} "
                      f"Orientierungskonflikte={top['orientation_conflicts']} "
                      f"Projektionsfaltungen={top['projected_orientation_flips']}")
            comb = row['nachher']['combined']
            print(f"    kombiniert (alle gueltigen Knoten): Komponenten={comb.get('components')} "
                  f"(erwartet {comb.get('components_expected')}), offene Kanten zwischen gemeinsamen Knoten="
                  f"{comb.get('open_edges_between_shared_nodes')}, Naehte (fehlende Kanten L/R, identisch)="
                  f"{ {k.split('|')[0][:12] + '|' + k.split('|')[1][:12]: (len(v['missing_left']), len(v['missing_right']), v['identical_nodes']) for k, v in comb.get('seams', {}).items()} }")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(rows, indent=1, default=lambda o: o.tolist()
                                       if hasattr(o, 'tolist') else str(o)), encoding='utf-8')
        print('Ergebnisse:', args.out)


if __name__ == '__main__':
    main()
