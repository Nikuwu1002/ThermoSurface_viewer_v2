# AP10 – Endprüfung, Dateisicherung und Reproduzierbarkeit

**Ziel:** Gesamtprüfung und Reproduktion auf einem passenden zweiten System, Testmatrix, Versions- und Ergebnisarchiv, Abschlussdokumentation.

## Dateien und Mechanismen

`ap10_verify.py` erstellt ein JSON-Protokoll (`results/ap10/ap10_audit.json`) mit UTC-Zeit, Python-/Betriebssysteminformationen, Paketversionen, fehlenden Abhängigkeiten, SHA-256-Summen aller Python-/YAML-Dateien, Baselineverfügbarkeit sowie Status von Syntax-, Unit- und (optional) GUI-/HEOS-Tests. Die Funktion `provenance()` greift auf die installierte Softwareumgebung zu, nicht auf angenommene Versionen.

Der Schalter `--static-only` startet nur Syntaxprüfung und ausführbare Unit-/Mocktestklassen. **`accepted` ist in diesem Modus immer `false`.** Ohne Schalter werden außerdem die bestehenden echten `AP4HEOSTests`, `AP4RenderTests` und `AP4GUISmokeTests` versucht. Eine vollständige Abnahme setzt erfolgreiche Tests **ohne SKIP**, verfügbare Bibliotheken, reale AP1-Snapshot-Baselines und explizites GUI-Opt-in voraus. Das Skript gibt dann und nur dann den Flag `accepted=true` zurück. Für die grafischen GUI-Tests unter Windows wird `THERMOSURFACE_AP4_GUI=1` benötigt.

Die neue `requirements.txt` enthält Versions**bereiche**, keinen in dieser Umgebung geprüften, exakten Lock. Bei der lokalen Installation die funktionierenden Versionen durch `python -m pip freeze > requirements-lock-local.txt` dokumentieren, für zweiten Rechner denselben Lock verwenden. Python 3.11+ empfohlen. Ein weiterer Rechner muss die Dateihashes, Testergebnisse und Beispielarrays erneut nachvollziehen.

## Verbindliche Schrittfolge auf dem Zielrechner

```powershell
python -m pip install -r requirements.txt
python -m pip freeze > requirements-lock-local.txt
python -m unittest Test.RegressionTests test_ap4.AP4SyntheticTests test_ap5_ap8 -v
python -m unittest test_ap5_ap8.AP5HEOSTests test_ap4.AP4HEOSTests test_ap4.AP4RenderTests -v
$env:THERMOSURFACE_AP4_GUI = "1"
python -m unittest test_ap4.AP4GUISmokeTests -v
python ap6_reference.py
python ap7_vdw_compare.py
python ap8_convergence.py
python ap10_verify.py
python main.py
```

Die **gespeicherten AP1-Baselines** des eigenen Originalprojekts sind einzuspielen (`regression/baseline` und ein gültiges Vergleichsverzeichnis). Die Original-Baselines fehlen im hier übermittelten Datensatz; sie werden *nicht* fingiert. Nach bewusstem, fachlich kontrolliertem Aufbau der neuen Referenz-/Vergleichssnapshots kann `python -m unittest Test.RecordedRegressionTests -v` laufen. Ein Ergebnisvergleich über Libraryversionen hinweg ist wegen der strikten Snapshot-Metadaten eingeschränkt.

Bei der GUI die fünf Diagrammtypen, Stoffwechsel Water/R134a, Phasengrenzen, p-h-T-Isothermen/Isobaren/Isochoren, Schalterverhalten, kritische Bereiche und Q-Fläche manuell oder über echte GUI-Tests prüfen. Die GUI-Tests prüfen nicht alle denkbaren Nutzerinteraktionen, insbesondere keine Power-Wall-Projektion.

## Tatsächlich ausgeführt und Status

- Gesamter Python-Quellcode: Syntax mit `compileall` bestanden.
- AP1-Regressionstestlogik + AP4-Synthetik + AP5–AP8-Unit-Tests: **23 Tests erfolgreich** (Details in `results/ap10/ap10_audit.json`).
- Vollständige HEOS-, Renderer- und GUI-Tests: **nicht als bestanden nachgewiesen**. Fehlende CoolProp, PyVista, PyVistaQt und PySide6.
- Echte AP1-Dateisnapshots fehlen aus den Uploads. `RecordedRegressionTests` kann deshalb keine grüne Endfreigabe liefern.
- Kein zweiter Rechner verfügbar; versionsidentische und Ergebnis-identische Reproduzierbarkeit extern **nicht nachgewiesen**.

**Status:** Audit- und Reproduktionsroutine implementiert und statisch erfolgreich geprüft; **vollständige Endabnahme ausdrücklich noch offen**.
