# AP18 – Abschluss, Regression und Performance

## Vorliegender Stand

Die bestehende Architektur mit fünf Diagrammtypen wurde soweit möglich unverändert übernommen. Die Erweiterung wird bei Auswahl von `pvt` berechnet, nicht bei `pht/pst/hsp/uvt`. Der bestehende Snapshot-Vertrag in `snapshot_pvt.py` enthält weiterhin die ursprünglichen Arrays; die neuen optionalen Ergebnisdaten sind dort **noch nicht als eigene freigegebene Baseline** enthalten. `compare_pvt.py` benötigt weiterhin eine echte vorher angelegte Baseline, die in den bereitgestellten Dateien fehlt. Die im Entwicklungsplan benannten Test-/Hilfsmodule `Test.py`, `test_ap4.py`, `test_ap5_ap8.py`, `ap10_verify.py`, `ap4_validation.py`, `ap6_reference.py` und `ap7_vdw_compare.py` waren im bereitgestellten Stand **nicht verfügbar**. Sie wurden nicht erfunden oder als ausgeführt ausgewiesen.

Die neue automatische Benchmarkmatrix `validate_ap17_ap18.py` benutzt drei Auflösungen `(nT,nscan,npoints)=(16,64,12),(32,128,20),(48,256,28)`. Berechnungszeit der beiden Spinodaläste und der zusätzlichen Flächen sowie Array-Speicherbedarf werden getrennt erfasst. Wegen fehlendem CoolProp wurden **keine HEOS-Performancezahlen gemessen**. Es ist daher wissenschaftlich nicht zulässig, bereits eine durch Messungen abgesicherte Standardrasterauflösung anzugeben; `homogeneous_points=20` ist eine vorläufige Bedienungs-Voreinstellung.

## Reproduktion am Windows-Zielrechner

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m compileall -q .
python -m unittest discover -v
python validate_ap17_ap18.py --output results/ap17_ap18_validation.json
python main.py
```

Die Test-JSON ist zu prüfen: `BLOCKED_NO_COOLPROP`/`HEOS_CORE_INCOMPLETE` sind **keine** Freigabe, und selbst `HEOS_CORE_PASS` deckt Qt-Interaktion sowie die fehlenden Alt-Baselines nicht ab. Für den grafischen Test die drei neuen Checkboxen, Oberflächengitter und jeden Diagrammwechsel kontrollieren; Bildschirmfotos verschiedener Kamerawinkel erstellen; kritische und negative Druckbereiche in Wasser und R134a prüfen. Die Befunde zusammen mit echten Versionen und Messwerten archivieren.

## Doppelprüfung und Ergebnis

**Prüfung A:** Ressourcenlimit, Trennung der Messstufen, Konfigurations- und Regressionseinfluss untersucht. **Prüfung B:** Syntax und 17 Tests ausgeführt, die fehlenden Realstoff- und GUI-Abnahmen explizit identifiziert. Der wissenschaftlich/numerische **Projektabschluss ist noch nicht möglich**, solange die vorgeschriebenen Stofftests, GUI-Tests, Baselines und Messungen fehlen. Die bereitgestellte ZIP ist eine implementierte, statisch und synthetisch geprüfte Entwicklungsfassung, **kein vollständig validiertes Release**.
