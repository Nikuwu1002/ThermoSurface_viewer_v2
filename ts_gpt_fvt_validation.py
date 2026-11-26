"""Validierungslauf g-p-T / f-v-T: Thermodynamik gegen unabhaengige CoolProp-Aufrufe, Laufzeit, Knotenzahlen.

Aufruf: python scripts/ts_gpt_fvt_validation.py   (schreibt results/gpt_fvt_validation.json)
Es werden je Stoff 300 zufaellige Stuetzknoten pro Teilflaeche mit PropsSI (neuer AbstractState, anderer Codepfad) geprueft:
  g = h - T s,  f = u - T s,  g = f + p v  (alle massenspezifisch, v_mass = v_molar / M)
Randknoten (Saettigung) und der kritische Zustand werden ausgenommen (dort ist (T,p) bzw. (T,rho) allein nicht eindeutig).
"""
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import CoolProp.CoolProp as CP  # noqa: E402
from CoolProp import AbstractState  # noqa: E402

import coolprop_heos_fvt as FVT  # noqa: E402
import coolprop_heos_gpt as GPT  # noqa: E402

FLUIDS = ["Water", "CO2", "Nitrogen", "Methane", "Argon", "R134a", "Ammonia", "Ethanol", "n-Propane", "n-Hexane"]
N_SAMPLES = 300


def rel(a, b, scale):
    return abs(a - b) / max(abs(b), scale)


def check_gpt(fluid, d, rng):
    M = float(AbstractState("HEOS", fluid).molar_mass())
    worst = {"g_vs_PropsSI": 0.0, "g_eq_h_minus_Ts": 0.0, "g_eq_f_plus_pv": 0.0}
    checked = 0
    for name in ("surface_liquid", "surface_vapor", "surface_high_temperature"):
        arr = d[name]
        for _ in range(N_SAMPLES):
            i = int(rng.integers(0, arr.shape[0]))
            j = int(rng.integers(0, arr.shape[1]))
            T, p, g = arr[i, j]
            on_edge = (name == "surface_liquid" and j == 0) or (name == "surface_vapor" and j == arr.shape[1] - 1)
            if on_edge or (name != "surface_high_temperature" and i == arr.shape[0] - 1) or (name == "surface_high_temperature" and i == 0):
                continue
            try:
                g_ref = CP.PropsSI("G", "T", T, "P", p, fluid)
                h = CP.PropsSI("H", "T", T, "P", p, fluid)
                s = CP.PropsSI("S", "T", T, "P", p, fluid)
                rho = CP.PropsSI("D", "T", T, "P", p, fluid)
                u = CP.PropsSI("U", "T", T, "P", p, fluid)
            except ValueError:
                continue
            scale = 8.31446261815324 / M * T
            f = u - T * s
            worst["g_vs_PropsSI"] = max(worst["g_vs_PropsSI"], rel(g, g_ref, scale))
            worst["g_eq_h_minus_Ts"] = max(worst["g_eq_h_minus_Ts"], rel(g, h - T * s, scale))
            worst["g_eq_f_plus_pv"] = max(worst["g_eq_f_plus_pv"], rel(g, f + p / rho, scale))
            checked += 1
    return worst, checked


def check_fvt(fluid, d, rng):
    M = float(AbstractState("HEOS", fluid).molar_mass())
    worst = {"f_vs_PropsSI": 0.0, "g_eq_f_plus_pv": 0.0, "two_phase_common_tangent_slope": 0.0}
    checked = 0
    for name in ("surface_liquid", "surface_vapor", "surface_high_temperature", "two_phase_surface"):
        arr = d[name]
        for _ in range(N_SAMPLES):
            i = int(rng.integers(0, arr.shape[0] - 1))
            j = int(rng.integers(0, arr.shape[1]))
            T, v, f = arr[i, j]
            rho = M / v
            try:
                props = {k: CP.PropsSI(k, "T", T, "D", rho, fluid) for k in ("P", "G", "S", "U")}
            except ValueError:
                continue
            scale = 8.31446261815324 / M * T
            worst["f_vs_PropsSI"] = max(worst["f_vs_PropsSI"], rel(f, props["U"] - T * props["S"], scale))
            worst["g_eq_f_plus_pv"] = max(worst["g_eq_f_plus_pv"], rel(props["G"], f + props["P"] / rho, scale))
            checked += 1
    mix = d["two_phase_surface"]
    for i in range(0, mix.shape[0] - 1, 5):
        psat = CP.PropsSI("P", "T", mix[i, 0, 0], "Q", 0, fluid)
        slopes = np.diff(mix[i, :, 2]) / np.diff(mix[i, :, 1])
        worst["two_phase_common_tangent_slope"] = max(worst["two_phase_common_tangent_slope"],
                                                      float(np.max(np.abs(slopes / (-psat / M) - 1.0))))
    return worst, checked


def main():
    rng = np.random.default_rng(20261111)
    report = {}
    for fluid in FLUIDS:
        entry = {}
        for kind, compute, check in (("gpt", GPT.berechne_gpt_daten, check_gpt), ("fvt", FVT.berechne_fvt_daten, check_fvt)):
            t0 = time.perf_counter()
            d = compute(fluid)
            elapsed = time.perf_counter() - t0
            worst, checked = check(fluid, d, rng)
            info = d["homogeneous_info"]["surfaces"]
            entry[kind] = {
                "compute_s": round(elapsed, 3), "checked_nodes": checked, "max_relative_deviation": worst,
                "homogeneous_valid_nodes": {k: [v["valid_nodes"], v["nodes"]] for k, v in info.items()},
                "below_equilibrium_hull_nodes": {k: v["below_equilibrium_hull_nodes"] for k, v in info.items()},
                "stable_nodes": int(sum(np.asarray(d[k]).reshape(-1, 3).shape[0] for k in
                                        ("surface_liquid", "surface_vapor", "surface_high_temperature"))),
            }
            print(f"{fluid:10s} {kind}: {elapsed:5.2f} s, geprueft {checked}, max. Abweichung "
                  f"{max(worst.values()):.2e}")
        report[fluid] = entry
    out = os.path.join(ROOT, "results", "gpt_fvt_validation.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=1)
    print("geschrieben:", out)


if __name__ == "__main__":
    main()
