# ThermoSurface – Entwicklungsstand AP11 bis AP18

**Ausgabe:** 21.10.2026 | **Status: implementierte, synthetisch geprüfte Entwicklungsversion; keine vollständige Realstoff-/GUI-Abnahme.**

## Fachliche Aussage

Die bisherige Zweiphasenfläche `p=p_sat(T)` gehört zum **thermodynamischen Gleichgewicht**, die Spinodalen gehören zu einer **homogenen, möglicherweise metastabilen oder mechanisch instabilen Fortsetzung**. Das neue Modul erzeugt deshalb zwei getrennte metastabile Flächen und auf Wunsch einen dritten instabilen Bereich. Die zusätzliche Fläche enthält niemals künstlich zwischen Sättigungs- und Spinodaldruck interpolierte Werte. Fehlen für einen Stoff tatsächliche homogene EOS-Werte, bleiben die betreffenden Flächenabschnitte leer. Die Referenzabbildung aus der Nualsri-Arbeit (Abbildung 71, gedruckte Seite 121) zeigt eine unkorrigierte van-der-Waals-Fläche, nicht die vorhandene HEOS-Gleichgewichtsfläche.

## Dateizuordnung

| Datei | Änderungen |
|---|---|
| `homogeneous_pvt.py` | **Neu:** zwei Spinodal-Wurzelzweige; homogene EOS-Auswertung; drei Flächenarrays; Randvalidierung |
| `main.py` | AP13-AP15-Berechnung hinzugefügt, neue optionale Ergebnis-Keys, alter p-v-T-Datenpfad bleibt erhalten |
| `visualization.py` | getrennte homogene PyVista-Actors, keine Netzbrücken über ungültige/logarithmisch undarstellbare Punkte |
| `gui.py` | drei zusätzliche p-v-T-Ebenen, Oberflächengitter-Einstellung erweitert |
| `Debug.py` | Diagnose zusätzlicher Flächen |
| `config.yaml` | `diagram_settings.pvt` für neue Flächenparameter |
| `test_ap13_ap18.py` | **Neu:** analytische numerische Tests und synthetischer Renderer-Test |
| `test_integration_mocked_main.py` | **Neu:** End-to-End-Verkabelung von `main` mit EOS-/GUI-Doubles |
| `validate_ap17_ap18.py` | **Neu:** echte HEOS-Stoffmatrix, PyVista-Mesh-Prüfung und 3 Auflösungen, sobald Abhängigkeiten installiert sind |
| `test_spinodal_extension.py` | Nicht vorhandene Legacy-Testabhängigkeit explizit als Skip statt Importfehler dokumentiert |
| `README_AP11.md` – `README_AP18.md` | separate Berichte je Arbeitspaket |
| `README_AP19.md`, `README_AP20.md` | **Neu (22.10.2026):** Korrektur der Ausreißerfläche (instabile Fortsetzung) und zusammenhängende Dreiflächen-Oberfläche inkl. Abdeckungs-/Netzkennzahlen; zugehörig `diagnose_ap19.py`, `validate_ap19_ap20.py`, `test_ap19_ap20.py`, `regression/ap19_vorher/` |
| `analyse_zustandsflaeche.md` | **Neu:** Ursachenanalyse H1–H6 und gekennzeichnete Ersatzfläche zwischen den Spinodalen (`zustandsflaeche.py`, `scripts/`, `test_zustandsflaeche.py`); neue config‑Keys `include_replacement`, `replacement_points`, `auto_pressure_axis` |

Die unveränderten Module für die anderen Diagrammtypen, Konvergenzskripte und ursprünglichen Vergleichsskripte sind im Archiv enthalten.

## Betriebsparameter

Unter `config.yaml` → `diagram_settings.pvt`:

```yaml
include_homogeneous: true
homogeneous_points: 20
include_unstable: false
# ab AP19/AP20 zusaetzlich:
unstable_scan_points: 256      # Scan je Isotherme fuer die zusammenhaengenden instabilen Aeste
critical_refinement_rows: 4    # Zusatz-T-Zeilen nahe Tc nur fuer homogene Flaechen (0 = aus)
close_critical_point: false    # Dreiecksfaecher zum kritischen Punkt, nur wenn validiert
```

Im p-v-T-Modus werden zusätzliche Flächen derzeit im Standardlauf mitberechnet, wenn `include_homogeneous` aktiv ist. Bei nicht unterstützten HEOS-Zuständen wird deren Berechnung mit Warnung eingeschränkt, statt das gesamte Gleichgewichtsdiagramm zu ersetzen. Die instabile Fläche ist standardmäßig ausgeschaltet und kann zur Diagnose aktiv geschaltet werden. Sie darf wissenschaftlich nur ausgewertet werden, wenn die Stabilitätsprüfungen bestanden sind.

## Teststand

| Testart | Stand |
|---|---|
| Vollständige Python-Syntax aller gelieferten `.py` | **Bestanden** |
| Unittest-Gesamtsuite | **18 Tests gelaufen: 15 bestanden, 3 übersprungen, 0 fehlgeschlagen** |
| Analytischer Test zweier Spinodaläste | **Bestanden** |
| Metastabile/instabile analytische Flächendaten | **Bestanden** |
| Exakte gemeinsame Randpunkte | **Bestanden an synthetischer EOS** |
| Triangulation ohne Sprung über NaN oder negative Drücke | **Bestanden mit PyVista-Double** |
| `main.py`-Datenverkabelung | **Bestanden mit Mocks**, kein echter HEOS-Lauf |
| Reale HEOS-Stoffmatrix (Water, R134a, CO2, Argon) | **BLOCKED:** CoolProp nicht installiert |
| Echte PyVista-/Qt-Oberfläche | **BLOCKED:** PyVista, PyVistaQt und PySide6 nicht installiert |
| Alte Regressionsbaselines und Legacy-Pakete | **BLOCKED:** nicht bereitgestellt |
| Reale Performance / mindestens drei Raster | **BLOCKED:** HEOS nicht installiert |

Maschinenlesbarer Laufbericht: `results/ap17_ap18_validation.json`, Status `BLOCKED_NO_COOLPROP` (Exitcode 2). Ausgeführte Tests in `results/unittest_output.txt`.

## Nachzuholende Schritte zur wissenschaftlichen Freigabe

1. Auf einem Computer mit CoolProp/PyVista/PySide6 `pip install -r requirements.txt` ausführen und die vier Reinstoffe über `python validate_ap17_ap18.py` prüfen. Die hierbei ausgewiesenen Lücken, Drücke und Residuen kritisch beurteilen.
2. Die tatsächliche GUI aller fünf Diagramme starten (`python main.py`), die drei neuen Ebenen sowie die Kurvenverläufe aus verschiedenen Perspektiven prüfen und Screenshots archivieren. Auf `T≈Tc` und nichtpositive Drücke besonders achten.
3. Die zuvor nicht bereitgestellten Referenz-Baselines und Tests wiederherstellen. Insbesondere die alte `Test.py` und die sonstigen im AP18-Plan benannten Dateien müssen für eine vollständige Regression vorliegen.
4. Drei Rasterauflösungen auf dem Zielrechner messen, geeignete Toleranzen mit Realstoffdaten festlegen und erst anschließend die finale Standardauflösung und die wissenschaftliche Vollabnahme bestimmen.

## Einschränkungen der Zielgeometrie

Die zusätzliche homogene Fläche besitzt absichtlich keinen künstlichen Dreiecksfächer zum kritischen Punkt. Erfordert der jeweilige Stoff nahe `Tc` nicht zuverlässig berechenbare HEOS-Zustände, kann dort eine echte Lücke bleiben. Ebenso kann ein Teil der Flüssigkeitsspinodale einen nichtpositiven Druck besitzen und deshalb im logarithmischen p-v-T-Renderer nicht angezeigt werden. Das sind **explizit ausgewiesene** physikalische bzw. numerische Grenzen, keine stillschweigend als valide dargestellten Geometrien.

Die im Entwicklungsplan geforderte vollständig abgeschlossene Umsetzung ist wegen dieser nicht ausführbaren Realstoff-, GUI- und Regressionstests **noch nicht erreicht**. Das ZIP liefert den geänderten Quellcode, die eigenständigen Tests, den Realstoff-Testläufer und nachvollziehbare Dokumentation zur Fertigstellung.
