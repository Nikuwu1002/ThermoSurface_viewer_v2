# AP17 – Thermodynamische und geometrische Validierung

## Neu erstellte Prüfinfrastruktur

`test_ap13_ap18.py` prüft ohne HEOS eine unabhängige analytische Modell-EOS: identifizierbare Maxima und Minima, Temperaturentwicklung der Äste, NaN-Lücken, negatives Flüssigkeitsspinodal-p, mathematisch korrekt berechnete Zwischenzustände, exakte Randknoten, Ablehnungsverhalten bei inkonsistenten Binodalen und Netzunterbrechung vor `log(p<=0)`. Die Alt-Testsuite `test_spinodal_extension.py` wurde um einen **expliziten Skip** bei der fehlenden, bisher importierten Datei `Test.py` ergänzt; ihre synthetischen Einzeltests sind unverändert.

`validate_ap17_ap18.py` enthält einen reproduzierbaren echten HEOS-Validierungslauf für **Water, R134a, CO2 und Argon**; er ermittelt Zahl der validen Spinodal- und Flächenzeilen, normiertes Residuum, echte gemeinsame Randknoten und transformierte Fehler. Bei installierter PyVista-Bibliothek prüft er Flächennetzpunktzahlen und positive Zellflächen. Drei zusätzliche Rasterkonfigurationen bilden den AP18-Benchmark. Ergebnisse landen als JSON unter `results/`.

## Tatsächliche Ergebnisse

- `python -m compileall -q .`: **bestanden** (Python-Syntax der gesamten mitgelieferten Codebasis).
- `python -m unittest discover -v`: **17 Tests, 14 bestanden, 3 übersprungen**, 0 fehlgeschlagen. Zwei Skips: HEOS fehlt; ein Skip: bislang benötigte Datei `Test.py` fehlt.
- `python validate_ap17_ap18.py`: **Exitcode 2: BLOCKED_NO_COOLPROP**; Bericht `results/ap17_ap18_validation.json`. Die Realstoffprüfungen wurden nicht ausgeführt.

## Akzeptanzstatus / Doppelprüfung

**Prüfung A:** Testmatrix, Formel- und Datenschnittstellenprüfung gegen Analytik. **Prüfung B:** Unabhängige Mehrfachprüfung des generischen Solvers und des Meshes; dabei Monotoniefehler in der Flüssigkeitsspinodale gefunden und behoben. Die synthetischen und statischen Abnahmekriterien sind erfüllt. Echte HEOS-Residuen, Volumen-/Druckabweichungen, realer kritischer Anschluss, dreidimensionale Kameraprüfung und die vollständigen Regressionstests sind **offen**. Die Diagnosefunktion `Debug.diagnose_homogeneous` ergänzt das bisherige Datenlogging.
