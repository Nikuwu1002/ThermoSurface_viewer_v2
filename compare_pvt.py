"""Strenger p-v-T-Regressionsvergleich. Keine HEOS-/GUI-Imports erforderlich."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent / "regression"
SCHEMA_VERSION = 2
EXPECTED_CASES = ("water_config", "water_auto", "argon_auto", "water_auto_debug")
CASE_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}\Z")


def load_manifest(directory: Path) -> list[str]:
    if not directory.is_dir():
        raise FileNotFoundError(f"Snapshot-Verzeichnis fehlt: {directory}")

    manifest_path = directory / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Manifest fehlt: {manifest_path}")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cases = manifest.get("cases")
    if (manifest.get("schema_version") != SCHEMA_VERSION
            or not isinstance(cases, list) or not cases
            or any(not isinstance(name, str) or not CASE_PATTERN.fullmatch(name)
                   for name in cases)
            or len(set(cases)) != len(cases)
            or tuple(cases) != EXPECTED_CASES):
        raise ValueError(f"Ungültiges/leeres Regressionsmanifest: {manifest_path}")

    expected = {"manifest.json"}
    for name in cases:
        expected.update((f"{name}.npz", f"{name}.json"))

    actual = {p.name for p in directory.iterdir()
              if p.is_file() and p.suffix in {".npz", ".json"}}
    if expected != actual:
        raise ValueError(f"Unvollständiger/zusätzlicher Snapshot in {directory}: "
                         f"fehlend={sorted(expected - actual)}, "
                         f"zusätzlich={sorted(actual - expected)}")
    return cases


def compare_runs(baseline: Path, candidate: Path) -> int:
    """Vergleicht alle Fälle; jeder fehlende Vergleich ist ein harter Fehler."""
    baseline, candidate = Path(baseline), Path(candidate)
    base_cases = load_manifest(baseline)
    new_cases = load_manifest(candidate)
    if base_cases != new_cases:
        raise AssertionError(f"Falllisten unterschiedlich: {base_cases} != {new_cases}")

    for name in base_cases:
        old = json.loads((baseline / f"{name}.json").read_text(encoding="utf-8"))
        new = json.loads((candidate / f"{name}.json").read_text(encoding="utf-8"))
        required = {"schema_version", "case", "stoff", "debug",
                    "raster", "visu_raster", "versions", "source_sha256"}
        if not required.issubset(old) or not required.issubset(new):
            raise ValueError(f"{name}: Metadatenfelder fehlen.")
        if old["schema_version"] != SCHEMA_VERSION or new["schema_version"] != SCHEMA_VERSION:
            raise ValueError(f"{name}: Unbekannte Metadatenschema-Version.")
        if old["case"] != name or new["case"] != name:
            raise ValueError(f"{name}: Fallname im JSON inkonsistent.")

        # Versions- und Parametervergleich ist strikt. Quellcode-Hashes sind
        # Provenienzangaben, da Codeänderungen gerade geprüft werden sollen.
        for key in ("stoff", "debug", "raster", "visu_raster", "versions", "pressure_axis"):
            if key not in old or key not in new:
                raise ValueError(f"{name}: {key} in den Snapshot-Metadaten nicht vorhanden")
            if old[key] != new[key]:
                raise AssertionError(f"{name}: Abweichung in {key}: "
                                     f"baseline={old[key]!r}, neu={new[key]!r}")

        with np.load(baseline / f"{name}.npz", allow_pickle=False) as ref, \
             np.load(candidate / f"{name}.npz", allow_pickle=False) as cur:
            if not ref.files or set(ref.files) != set(cur.files):
                raise AssertionError(f"{name}: Array-Keys fehlen oder unterscheiden sich: "
                                     f"{ref.files} != {cur.files}")
            for key in sorted(ref.files):
                a, b = ref[key], cur[key]
                if a.shape != b.shape or a.dtype != b.dtype:
                    raise AssertionError(f"{name}/{key}: Shape oder dtype verändert: "
                                         f"{a.shape}/{a.dtype} -> {b.shape}/{b.dtype}")
                if not (np.issubdtype(a.dtype, np.floating)
                        or np.issubdtype(a.dtype, np.integer)
                        or np.issubdtype(a.dtype, np.bool_)):
                    raise TypeError(f"{name}/{key}: Nichtnumerischer Datentyp: {a.dtype}")
                try:
                    if a.dtype == np.bool_:
                        np.testing.assert_array_equal(a, b)
                    else:
                        np.testing.assert_array_equal(np.isnan(a), np.isnan(b))
                        np.testing.assert_allclose(a, b, rtol=1e-9,
                                                   atol=1e-12, equal_nan=True)
                except AssertionError as exc:
                    raise AssertionError(f"{name}/{key}: {exc}") from exc

        if old["source_sha256"] != new["source_sha256"]:
            print(f"{name}: INFO – Quellcodedateien wurden verändert.")
        print(f"{name}: OK (alle Arrays und Einstellungen verglichen)")

    print(f"REGRESSION BESTANDEN: {len(base_cases)} Fälle")
    return len(base_cases)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", help="Name des aktuellen Snapshot-Laufs")
    args = parser.parse_args()
    if not CASE_PATTERN.fullmatch(args.run) or args.run == "baseline":
        parser.error("Ungültiger Laufname (oder 'baseline' als Kandidat).")
    compare_runs(ROOT / "baseline", ROOT / args.run)


if __name__ == "__main__":
    main()
