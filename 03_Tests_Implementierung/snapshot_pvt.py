"""Erzeugt vollständige p-v-T-Snapshots; Baselines nur nach Freigabe."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import shutil
import sys
import tempfile
from importlib.metadata import version
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
PROGRAM_ROOT = ROOT.parent / "01_Programm"
REGRESSION = ROOT / "regression"
sys.path.insert(0, str(PROGRAM_ROOT))

import main as app
CASE_NAMES = ("water_config", "water_auto", "argon_auto", "water_auto_debug")
ARRAY_NAMES = (
    "surface", "isothermen", "isochoren", "binodal_liquid", "binodal_vapor",
    "spinodale_liquid", "spinodale_vapor", "critical_point", "phase_mask",
    "two_phase_surface", "Q_values",
    "metastable_liquid_surface", "metastable_vapor_surface", "unstable_surface",
    "homogeneous_spinodal_liquid", "homogeneous_spinodal_vapor",
)
SOURCE_NAMES = (
    "main.py", "Debug.py", "coolprop_heos_common.py", "coolprop_heos_pvt.py",
    "coolprop_heos_pht.py", "coolprop_heos_pst.py", "coolprop_heos_hsp.py",
    "coolprop_heos_uvt.py", "gui.py", "visualization.py", "config.yaml", "limits.py",
    "homogeneous_pvt.py", "spinodal_vapor_solver.py",
)
SCHEMA_VERSION = 2
CASE_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}\Z")


def metadata_versions() -> dict[str, str]:
    versions = {"python": platform.python_version()}
    for pkg in ("numpy", "CoolProp", "PyYAML", "pyvista", "pyvistaqt", "PySide6"):
        versions[pkg] = version(pkg)
    return versions


def code_hashes() -> dict[str, str]:
    return {name: hashlib.sha256((PROGRAM_ROOT / name).read_bytes()).hexdigest()
            for name in SOURCE_NAMES}


def save_case(directory: Path, name: str, stoff: str,
              raster: dict, visu: dict, debug: bool,
              versions: dict, hashes: dict) -> None:
    app.DEBUG_MODE = debug
    app.configure_debug(debug)
    daten = app.berechne_daten(stoff, raster, visu)

    arrays = {key: np.asarray(daten[key]) for key in ARRAY_NAMES}
    for i, curve in enumerate(daten["isobaren"]):
        arrays[f"isobare_{i:03d}"] = np.asarray(curve)
    if not daten["isobaren"] or arrays["surface"].ndim != 3:
        raise ValueError(f"{name}: Leere Isobaren oder ungültige Oberflächenform.")
    if any(not np.issubdtype(arr.dtype, np.number) and arr.dtype != np.bool_
           for arr in arrays.values()):
        raise TypeError(f"{name}: nichtnumerische Arraydaten.")

    np.savez_compressed(directory / f"{name}.npz", **arrays)
    settings = {
        "schema_version": SCHEMA_VERSION,
        "case": name,
        "stoff": stoff,
        "debug": debug,
        "raster": raster,
        "visu_raster": visu,
        "versions": versions,
        "platform": platform.platform(),
        "source_sha256": hashes,
        "pressure_axis": daten.get("pressure_axis", "log"),
        "homogeneous_state": daten.get("homogeneous_validation", {}).get("state"),
    }
    (directory / f"{name}.json").write_text(
        json.dumps(settings, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"{name}: gespeichert, surface={arrays['surface'].shape}")


def create_snapshot(run: str, approve_baseline: bool = False) -> Path:
    if not CASE_PATTERN.fullmatch(run):
        raise ValueError("Laufname: nur Buchstaben, Ziffern, _ und - verwenden.")
    if run == "baseline" and not approve_baseline:
        raise PermissionError("Baseline benötigt ausdrücklich --approve-baseline.")
    if run != "baseline" and approve_baseline:
        raise ValueError("--approve-baseline ist nur mit 'baseline' erlaubt.")

    target = REGRESSION / run
    if target.exists():
        raise FileExistsError(f"Snapshot existiert bereits: {target}. Nicht überschreiben!")
    REGRESSION.mkdir(exist_ok=True)
    temp = Path(tempfile.mkdtemp(prefix="._snapshot_", dir=REGRESSION))
    try:
        _stoff, raster, visu = app.lade_konfiguration()
        versions, hashes = metadata_versions(), code_hashes()
        save_case(temp, "water_config", "Water", raster, visu, False, versions, hashes)
        wr, wv, _ = app.generate_stoff_grid("Water", raster, visu)
        save_case(temp, "water_auto", "Water", wr, wv, False, versions, hashes)
        ar, av, _ = app.generate_stoff_grid("Argon", raster, visu)
        save_case(temp, "argon_auto", "Argon", ar, av, False, versions, hashes)
        save_case(temp, "water_auto_debug", "Water", wr, wv, True, versions, hashes)
        (temp / "manifest.json").write_text(
            json.dumps({"schema_version": SCHEMA_VERSION, "cases": CASE_NAMES}, indent=2),
            encoding="utf-8")
        temp.rename(target)
    finally:
        if temp.exists():
            shutil.rmtree(temp)
        app.DEBUG_MODE = False
        app.configure_debug(False)
    print(f"Snapshot vollständig: {target}")
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", help="Laufname, etwa 'baseline' oder 'nach_ap1'")
    parser.add_argument("--approve-baseline", action="store_true",
                        help="Einmalig neue Referenzdateien anlegen")
    args = parser.parse_args()
    create_snapshot(args.run, args.approve_baseline)


if __name__ == "__main__":
    main()
