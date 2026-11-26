# Tatsächliches Testprotokoll – AP5–AP10 (19.10.2026)

## Ausführbare Tests

Befehl: `python -m unittest Test.RegressionTests test_ap4.AP4SyntheticTests test_ap5_ap8.AP5IntegrationTests test_ap5_ap8.AP6ReferenceTests test_ap5_ap8.AP7AnalyticalTests test_ap5_ap8.AP8NumericalTests -v` (aus `ap10_verify.py`).

- **Anzahl:** 23 ausgeführte Tests, alle erfolgreich.
- **AP1-Regressionstestlogik:** 6 erfolgreich; hier werden ausschließlich synthetische `fixture()`-Snapshots verwendet, **keine echten Projekt-Baselines**.
- **AP4-Synthetik:** 7 erfolgreich; reine Struktur-/Transformations-/Mocktests.
- **AP5:** 4 erfolgreich; inklusive Original-Callback, Daten-Dictionary, originalem Renderer und NaN-Unterbrechungen im Mock.
- **AP6:** 2 erfolgreich; Referenz-Zahlenkonstanten und Fehlerrechnung.
- **AP7:** 2 erfolgreich; analytische kritische Beziehungen und Ausnahmeprüfungen.
- **AP8:** 2 erfolgreich; bilineare Interpolation und Rasterwahl.
- **Syntax:** `python -m compileall -q .` / `compileall.compile_dir(...)` erfolgreich.

## Nicht abgeschlossene Realtests

Befehl: `python ap10_verify.py`, **Exit-Code 2** (gewollt, Abnahmekriterien nicht erfüllt). `test_ap5_ap8.AP5HEOSTests`, `test_ap4.AP4HEOSTests`, `test_ap4.AP4RenderTests`, `test_ap4.AP4GUISmokeTests`:

- **9 Testmethoden**, **9 als SKIP** protokolliert, **0 tatsächlich durchgeführte Realprüfungen**.
- Abhängigkeiten fehlen: `CoolProp, pyvista, pyvistaqt, PySide6`.
- `baseline_available=false`: Originale Snapshot-Verzeichnisse wurden nicht hochgeladen.
- Das separate `python -m unittest Test.RecordedRegressionTests -v` wurde ausgeführt und schlug wegen `FileNotFoundError: regression/baseline` fehl. Dieser erwartete Fehlzustand wurde **nicht** als bestandener Test gezählt.
- Die getrennt versuchten Skripte `ap6_reference.py`, `ap7_vdw_compare.py` und `ap8_convergence.py` scheiterten jeweils mit `ModuleNotFoundError: No module named 'CoolProp'`. Keine HEOS-Vergleichstabellen, Zahlenwerte oder Performancewerte wurden generiert.

## Artefakte und weitere Prüfung

`results/ap10/ap10_audit.json` enthält bei tatsächlicher Prüfung den Paketstatus, die Hashsummen und Einzeltestausgaben. Das analytische Bild `figures/vdw_reduced_isotherms.png` wurde dagegen erfolgreich mit NumPy/Matplotlib erzeugt und ist **kein** Ergebnis des HEOS-Modellvergleichs.

**Schlussfolgerung:** Softwareteilbereiche sind real synthetisch/analytisch geprüft, die **wissenschaftliche HEOS-Abnahme und die endgültige Reproduzierbarkeit fehlen**. Es besteht insbesondere kein begründeter Grund, die Arbeitspakete AP6–AP10 als vollständig abgeschlossen zu kennzeichnen.
