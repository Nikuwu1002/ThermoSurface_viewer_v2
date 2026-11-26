# ThermoSurface-Projektstruktur

## Kategorien

1. [01_Programm](../01_Programm): Anwendungsmodule, `main.py`, `config.yaml` und `requirements.txt`.
2. [02_Tests_Archiv](../02_Tests_Archiv): einmalige Diagnose- und Untersuchungsdateien; sie werden nicht automatisch mit der Testsuite ausgeführt.
3. [03_Tests_Implementierung](../03_Tests_Implementierung): wiederverwendbare Tests und Validierungen sowie `scripts`, Regressionseingaben, Messergebnisse und Beispielbilder.
4. [04_Dokumentation](.): Projekt- und Arbeitspaketdokumentation, PDF-/ZIP-Archive und Literatur.

Git-Metadaten, `.gitignore`, VS-Code-Einstellungen und Test-/Python-Caches bleiben außerhalb dieser Kategorien im Projektstamm.

## Start und Tests

Die Befehle im Projektstamm ausführen:

```powershell
python -m pip install -r 01_Programm/requirements.txt
python 01_Programm/main.py
python -m pytest -q 03_Tests_Implementierung
```

Die GUI benötigt eine OpenGL-fähige Umgebung. Die dokumentierten historischen Arbeitspaketberichte behalten ihre damaligen Befehle und Testergebnisse als Aufzeichnungen; die Befehle oben entsprechen der aktuellen Ordnerstruktur.
