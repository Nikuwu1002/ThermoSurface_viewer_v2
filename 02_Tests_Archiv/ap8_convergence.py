"""AP8: HEOS-Rasterkonvergenz bei einphasigem ueberkritischem T-v-Gas.

Unabhaengige, direkt neu berechnete HEOS-Zwischenpunkte messen bilinearen
Interpolationsfehler (nicht die Differenz gemeinsamer Rasterknoten).
Laufzeit und Peak-RSS sind system-/lastabhaengig; keine GUI-Gesamtzeit.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
from time import perf_counter
import numpy as np

DENSITIES = ((21, 41), (41, 81), (81, 161))


def bilinear(t, v, T_grid, V_grid, P):
    """Interpolation innerhalb einer regulaeren Rasterzelle."""
    T_grid = np.asarray(T_grid)
    V_grid = np.asarray(V_grid)
    P = np.asarray(P)
    if P.shape != (len(T_grid), len(V_grid)) or len(T_grid) < 2 or len(V_grid) < 2:
        raise ValueError('Unpassende Rasterdimensionen')
    if (not np.all(np.diff(T_grid) > 0) or not np.all(np.diff(V_grid) > 0)
            or not T_grid[0] <= t <= T_grid[-1] or not V_grid[0] <= v <= V_grid[-1]):
        raise ValueError('Punkt ausserhalb des Interpolationsrasters')
    i = min(int(np.searchsorted(T_grid, t, side='right')-1), len(T_grid)-2)
    j = min(int(np.searchsorted(V_grid, v, side='right')-1), len(V_grid)-2)
    ti = (t - T_grid[i])/(T_grid[i+1] - T_grid[i])
    vj = (v - V_grid[j])/(V_grid[j+1] - V_grid[j])
    return ((1-ti)*(1-vj)*P[i,j] + ti*(1-vj)*P[i+1,j]
            + (1-ti)*vj*P[i,j+1] + ti*vj*P[i+1,j+1])


def _memory_rss():
    try:
        import psutil
        return psutil.Process(os.getpid()).memory_info().rss
    except ImportError:
        return None


def calculate(fluid='Water', counts=DENSITIES, n_tests=225, repetitions=2, seed=42):
    """Vergleicht drei Punktdichten mit denselben zufaelligen Testzuständen."""
    from CoolProp import AbstractState
    import CoolProp.CoolProp as CP
    state = AbstractState('HEOS', fluid)
    Tc = float(state.T_critical())
    vc = 1. / float(state.rhomolar_critical())
    limits = ((1.07*Tc, 1.25*Tc), (2*vc, 8*vc))
    if limits[0][1] > state.Tmax():
        raise ValueError('Gewaehlter Temperaturbereich liegt ausserhalb des HEOS-Modells')
    if not (isinstance(repetitions, int) and repetitions >= 1 and n_tests >= 20):
        raise ValueError('Ungueltige Benchmarkparameter')
    rng = np.random.default_rng(seed)
    # Testpunkte gleichverteilt im Gebiet; dieselben Zustaende pro Aufloesung.
    probe = np.column_stack((rng.uniform(*limits[0], n_tests), rng.uniform(*limits[1], n_tests)))
    reference = []
    for t, v in probe:
        state.update(CP.DmolarT_INPUTS, 1./float(v), float(t))
        p = float(state.p())
        if not math.isfinite(p) or p <= 0:
            raise ValueError('HEOS-Referenzdruck ausserhalb des Einphasenbereichs')
        reference.append(p)
    reference = np.asarray(reference)
    results = []
    for nt, nv in counts:
        if nt < 2 or nv < 2 or nt*nv > 1_000_000:
            raise ValueError('Ungueltiges Raster')
        T = np.linspace(*limits[0], nt)
        V = np.linspace(*limits[1], nv)
        times, peaks = [], []
        P = None
        for rep in range(repetitions):
            before = _memory_rss()
            start = perf_counter()
            grid = np.empty((nt, nv), dtype=float)
            for i, temp in enumerate(T):
                for j, volume in enumerate(V):
                    state.update(CP.DmolarT_INPUTS, 1./float(volume), float(temp))
                    grid[i, j] = float(state.p())
            elapsed = perf_counter()-start
            times.append(elapsed)
            after = _memory_rss()
            if after is not None:
                peaks.append(after)
            P = grid
        interp = np.array([bilinear(t,v,T,V,P) for t,v in probe])
        errors = np.abs(interp-reference)/reference
        results.append(dict(fluid=fluid, n_T=nt, n_v=nv, n_points=nt*nv,
                            duration_median_s=float(np.median(times)),
                            grid_bytes=int(P.nbytes),
                            process_rss_after_max_bytes=max(peaks) if peaks else None,
                            max_relative_error=float(np.max(errors)),
                            mean_relative_error=float(np.mean(errors)),
                            rms_relative_error=float(np.sqrt(np.mean(errors**2))),
                            test_states=n_tests, seed=seed, repetitions=repetitions))
    return results, limits


def choose_resolution(rows, threshold=1e-3):
    """Kleinste gepruefte Aufloesung mit max. relativer Interpolationsabweichung <= Ziel."""
    if threshold <= 0:
        raise ValueError('Zielschwelle muss positiv sein')
    feasible = [r for r in rows if r['max_relative_error'] <= threshold]
    if not feasible:
        return {'threshold': threshold, 'recommended': None,
                'reason': 'Keine der untersuchten Rasterdichten erfuellt die Fehlerschwelle.'}
    smallest = min(feasible, key=lambda r:r['n_points'])
    return {'threshold': threshold, 'recommended': [smallest['n_T'],smallest['n_v']],
            'reason': 'Kleinstes getestetes Raster im untersuchten einphasigen Teilgebiet.'}


def write_results(rows, limits, directory, tolerance):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    with (directory/'ap8_grid_benchmark.csv').open('w', encoding='utf-8', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    summary = dict(domain_T_K=limits[0], domain_v_m3mol=limits[1],
                   results=rows, choice=choose_resolution(rows, tolerance),
                   warning='Nur ueberkritischer homogener T-v-Rasterbereich, kein p-v-T-Gesamtlauf mit Phasenflaechen und Rendering.')
    (directory/'ap8_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fluid',default='Water')
    parser.add_argument('--output', default=str(Path(__file__).resolve().parent / 'results' / 'ap8'))
    parser.add_argument('--repetitions',type=int,default=2)
    parser.add_argument('--threshold',type=float,default=0.001)
    args=parser.parse_args()
    rows, limits=calculate(args.fluid, repetitions=args.repetitions)
    summary=write_results(rows,limits,args.output,args.threshold)
    print(json.dumps(summary['choice'],indent=2))


if __name__=='__main__':
    main()
