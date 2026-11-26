"""Ausfuehrbarer HEOS-Zielrechnervergleich und Konvergenzbenchmark.

Aufruf im Projektordner: python 03_Tests_Implementierung/validate_spinodal_heos.py
Erzeugt spinodal_heos_validation.json (kein Ueberschreiben vorhandener Baselines).
Ohne CoolProp bricht dieses Skript sichtbar ab; synthetische Tests separat.
"""
from __future__ import annotations

import json
import platform
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np

PROGRAM_DIR = Path(__file__).resolve().parents[1] / "01_Programm"
sys.path.insert(0, str(PROGRAM_DIR))

from coolprop_heos_common import critical_point
from coolprop_heos_pvt import (calc_spinodale, filter_spinodale,
                                calc_spinodale_vapor_temperature,
                                compare_spinodal_overlap, prepare_spinodal_branch)
from Debug import test_spinodale_ableitung


def _versions():
    values = {"python": platform.python_version()}
    for name in ("CoolProp", "numpy", "scipy", "PySide6", "pyvista", "pyvistaqt"):
        try:
            values[name] = version(name)
        except Exception:
            values[name] = "nicht installiert"
    return values


def _test_fluid(fluid, start, stop):
    critical = critical_point(fluid)
    Tc, pc, vc = critical
    raw_liquid, raw_vapor = calc_spinodale(fluid)
    old, old_geo = prepare_spinodal_branch(raw_vapor, critical, start, stop, "vapor")
    results = {"fluid": fluid, "critical": [float(x) for x in critical],
               "legacy": {"geometry": old_geo,
                          "derivative": test_spinodale_ableitung(
                              fluid, np.empty((0,3)), old, critical)["vapor"]},
               "runs": []}
    for n_T, n_scan in [(250,200), (500,400), (500,800)]:
        arr, metadata = calc_spinodale_vapor_temperature(
            fluid, start, stop, n_T=n_T, n_scan=n_scan, with_diagnostics=True)
        curve, geometry = prepare_spinodal_branch(arr, critical, start, stop, "vapor")
        AP3 = test_spinodale_ableitung(
            fluid, np.empty((0,3)), curve, critical)["vapor"]
        results["runs"].append({"n_T": n_T, "n_scan": n_scan,
                               "numerical": metadata,
                               "geometry": geometry,
                               "ap3_derivative": AP3,
                               "old_overlap": compare_spinodal_overlap(
                                   arr, old, critical_T=Tc)})
        print(f"{fluid} {n_T}x{n_scan}: {metadata['valid_points']}/{n_T}, "
              f"T={metadata['first_valid_T']}..{metadata['last_valid_T']} K "
              f"Rmax={metadata['max_residual']}, t={metadata['elapsed_seconds']:.2f}s")
    return results


def main():
    import CoolProp.CoolProp as CP
    from CoolProp import AbstractState
    results = {"versions": _versions(), "water_reference": [], "fluids": []}
    water = AbstractState("HEOS", "Water")
    Tc, pc, vc = critical_point("Water")
    rho_c = 1/vc
    reference = {
        580: (.274950, 3.637028, .515090, 1.202898),
        600: (.355207, 2.815258, .628857, 1.123961),
        620: (.464580, 2.152482, .762474, 1.058026),
        630: (.541137, 1.847960, .839532, 1.030857),
        640: (.658638, 1.518284, .926810, 1.009076),
    }
    for T, expected in reference.items():
        array, status = calc_spinodale_vapor_temperature(
            "Water", float(T), float(T)+0.0001, n_T=2, n_scan=300,
            with_diagnostics=True)
        water.update(CP.QT_INPUTS, 1., float(T))
        psat = float(water.p())
        value = array[0]
        actual = (float(1/value[1]/rho_c), float(value[1]/vc),
                  float(value[2]/pc), float(value[2]/psat))
        diff = [abs(a-b) for a,b in zip(actual, expected)]
        row = {"T_K": T, "measured_old": expected, "new": actual,
               "absolute_difference": diff, "max_residual": status["max_residual"],
               "passed": bool(np.allclose(actual, expected, rtol=1e-4, atol=2e-5))}
        print(f"Water {T} K: {actual}; match={row['passed']}")
        results["water_reference"].append(row)
    results["fluids"].append(_test_fluid("Water", 580.0, 647.0))
    state = AbstractState("HEOS", "R134a")
    results["fluids"].append(_test_fluid(
        "R134a", 0.9*state.T_critical(), 0.997*state.T_critical()))
    out = Path(__file__).with_name("spinodal_heos_validation.json")
    if out.exists():
        raise FileExistsError(f"Alten validierten Ergebnisdatensatz nicht ueberschreiben: {out}")
    out.write_text(json.dumps(results, ensure_ascii=False, indent=2, allow_nan=False),
                   encoding="utf-8")
    print(f"Messbericht erstellt: {out}")


if __name__ == "__main__":
    main()
