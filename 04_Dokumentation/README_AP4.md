# ThermoSurface – AP4: Systematische Validierung aller fünf 3D-Diagrammtypen

**Stand:** 19. Oktober 2026  
**Projektgrundlage:** `ThermoSurface_Projektanalyse_und_Arbeitspakete_AP1-AP10.md`, AP4  
**Lieferform:** Update-Paket, **keine Neuinstallation des Gesamtprojekts**. Die neuen Dateien in denselben Ordner wie `main.py` legen. Bestehende AP1-/AP2-/AP3-Dateien, Konfigurationen und Regressions-Baselines nicht überschreiben.

> **Validierungsstatus:** AP4-Tests und Prüflogik sind im Python-Code implementiert. Synthetische Tests, statische Prüfungen und vorhandene AP1-/AP3-Regressionstests sind hier erfolgreich gelaufen. Die **verbindliche numerische Abnahme mit echten CoolProp-HEOS-Zuständen für zwei Fluide pro Diagramm sowie die GUI-/VTK-Abnahme stehen noch aus**, da CoolProp, PyVista, PyVistaQt und PySide6 in der Ausführungsumgebung nicht verfügbar sind. Eine vollständig nachgewiesene AP4-Abnahme wird ausdrücklich nicht behauptet.

## 1. Übersicht und Abgrenzung

AP4 heißt gemäß Projektanalyse **„Alle fünf Diagramme validieren“**. Verbindliches Abnahmekriterium: Jeder der fünf Diagrammtypen muss für **mindestens zwei repräsentative Fluide** numerisch überprüft werden; **GUI und Renderer müssen auf dem Zielrechner** praktisch getestet werden.

| Diagramm | Rohdatenkoordinaten (SI) | Zu prüfende Sachverhalte |
|---|---|---|
| p-v-T | `(T [K], v_m [m³/mol], p [Pa])` | Stabile Fläche, Phasenmaske, Binodalen, kritischer Punkt, getrennte Q-Fläche; AP2 bleibt erhalten |
| p-h-T | `(T [K], h [J/kg], p [Pa])` | Drei Einphasen-Teilnetze, kritischer Anschluss, Siede-/Taulinie, Q-Mischung |
| p-s-T | `(T [K], s [J/(kg K)], p [Pa])` | Wie p-h-T, mit Entropie; negative Entropiereferenzwerte sind grundsätzlich zulässig |
| h-s-p | `(h [J/kg], s [J/(kg K)], p [Pa])` | Konsistente Zusammensetzung aus p-h-T und p-s-T bei gleichem `(T,p)`; Mischung |
| u-v-T | `(T [K], v_m [m³/mol], u [J/kg])` | Molares Volumen und massenspezifische innere Energie; Einphasen- und Q-Flächen |

**Prüffluide:** `Water` und `R134a` mit dem bestehenden Backend `HEOS`. Beide sind in der GUI-Stoffliste vorhanden. Diese Kombination soll unterschiedliche kritische Größen und Dampf-Flüssig-Gleichgewichte abdecken. Sie ist **keine unabhängige Stoffdatenreferenz**.

**Bereits abgeschlossen:** AP1 (Testbasis), AP2 (p-v-T-Randknoten), AP3 (Spinodalen). Diese Programmdateien werden nicht verändert. Keine Umsetzung von AP5 (p-h-T-Isolinienintegration), AP6 (externe Referenzdaten), AP7 (vdW), AP8 (Performance), AP9 (Studienarbeit), AP10 (Endprüfung).

## 2. Ausgangszustand vor AP4

Die Codeinspektion des bereitgestellten Standes ergab:

- `main.berechne_auswahl(...)` dispatcht anhand des Diagrammtyps an `berechne_daten()` beziehungsweise `berechne_pht_daten()`, `berechne_pst_daten()`, `berechne_hsp_daten()` und `berechne_uvt_daten()`.
- Die Module `coolprop_heos_pht.py`, `coolprop_heos_pst.py`, `coolprop_heos_hsp.py` und `coolprop_heos_uvt.py` verfügen bereits über geometrische/strukturelle Validatoren (`validate_pht_data`, `validate_hsp_data`, `validate_uvt_data`) beziehungsweise wiederverwendete Prüfungen.
- `coolprop_heos_pvt.py` verfügt bereits über `validate_pvt_phase_boundaries()`; die Spinodalenprüfung stammt aus AP3.
- `gui.py` bietet fünf Diagrammtypen und die Darstellungsebenen; `visualization.py` stellt die fünf zugehörigen Renderer `render_pvt`, `render_pht`, `render_pst`, `render_hsp`, `render_uvt` bereit.
- Die vorhandene `Test.py` deckt insbesondere AP1/AP2/AP3 ab. Eine **einheitliche, explizite Zwei-Fluid-Validierung aller fünf Diagramme** einschließlich Prüfung physikalischer Achsen und realer Render-/GUI-Pfade existierte nicht.

### Ist-Soll-Vergleich

| AP4-Anforderung | Ist vor AP4 | Soll | Umsetzung / betroffene Dateien |
|---|---|---|---|
| Alle fünf Diagrammtypen strukturell prüfen | Unterschiedliche vorhandene Einzelvalidatoren | Einheitlich auf Arrayform, Randknoten, Q, Tc prüfbar | `ap4_validation.py` (neu) |
| Zwei reale Fluide je Diagramm | Keine gemeinsame nachweisbare Testmatrix | Water und R134a, zehn reale Prüffälle | `test_ap4.py` (neu) |
| Kritische Zustände und Phasengrenzen prüfen | Teilweise numerische Selbsttests | Kritische HEOS-Daten, Q=0/1, Sättigung und Einphasen-Stichproben | `test_ap4.py` |
| SI-Einheiten, Datentransformation | Formeln bereits in Renderer enthalten | Kontrollierte Größen- und Achsenzuordnung | `ap4_validation.py`, `test_ap4.py` |
| GUI-/Rendererpfad prüfen | Ausführbare Funktionen, noch kein Fünf-Diagramm-Smoke-Test | Fünf reale Renderer + GUI-Auswahl und Sichtbarkeit | `test_ap4.py` (Zielrechner) |
| AP1–AP3 nicht beeinträchtigen | Bestehende Regressionstests | Alte Module unverändert; Regressionen erneut ausführen | Keine Produktivcode-Änderung |
| Dokumentation und Abnahmestatus | Keine AP4-README | Reproduzierbare Testmatrix mit ehrlichem Ergebnis | `README_AP4.md` (neu) |

**Strategieentscheidung:** Statt jeden Berechnungspfad mit zusätzlicher Debug-Logik zu versehen, wird eine **unabhängige, wiederverwendbare Prüfschicht** eingesetzt. Das vermeidet neues Kopplungsverhalten im GUI-/Produktions-Dispatcher, Duplikation bestehender Validatoren und versehentliche Änderungen an AP1 bis AP3. Der Validator prüft *Datenkonsistenz*, nicht die thermodynamische Richtigkeit durch externe Quellen. Die echte HEOS-Nachrechnung bleibt ein separater Testpfad.

## 3. Tatsächlich geänderte beziehungsweise neue Dateien

| Datei | Funktion / Bereich | Vor AP4 | Durchgeführte Änderung | Grund |
|---|---|---|---|---|
| **Neu:** `ap4_validation.py` | `validate_diagram_data()`, Hilfsfunktionen, `AP4ValidationError` | Nicht vorhanden | Pure-NumPy-Validierung aller fünf Ergebnis-Dictionaries | Gemeinsamer prüfbarer Vertrag ohne GUI-/HEOS-Abhängigkeit |
| **Neu:** `test_ap4.py` | `AP4SyntheticTests` | Nicht vorhanden | Synthetische positive/negative Array-, Q-, Tc- und Achsentests | Sofort reproduzierbar, auch ohne CoolProp |
| **Neu:** `test_ap4.py` | `AP4HEOSTests` | Nicht vorhanden | Reale Fluidmatrix, HEOS-Flash-Kontrollen, Modellgrenzen | Zielrechner-Abnahme für zwei Fluide je Diagramm |
| **Neu:** `test_ap4.py` | `AP4RenderTests` | Nicht vorhanden | Fünf echte PyVista-Renderer im Offscreen-Plotter | VTK-/Actor-Schnittstellen überprüfen |
| **Neu:** `test_ap4.py` | `AP4GUISmokeTests` | Nicht vorhanden | Diagrammwechsel, Renderpfad und Checkboxen im echten Qt-Fenster | GUI-/Renderer-Integration am Zielrechner |
| **Neu:** `README_AP4.md` | Gesamtdokumentation | Nicht vorhanden | Implementierungs- und Testprotokoll | Reproduzierbarkeit und ehrliche Abnahme |

**Keine Änderungen** an `main.py`, `gui.py`, `visualization.py`, `coolprop_heos_*.py`, `Debug.py`, `Test.py`, `limits.py`, `config.yaml`, `snapshot_pvt.py` oder `compare_pvt.py`. Bereits implementierte Funktionen werden bewusst weitergenutzt. Neue Imports sind auf die beiden neuen Python-Dateien beschränkt. Keine öffentlichen Funktionssignaturen oder Rückgabeformate des Programms wurden geändert; kein bestehender Code wurde entfernt.

## 4. Durchführung und Implementierungsschritte

| Schritt | Datei und Codestelle | Tatsächlich umgesetzt | Prüfung |
|---|---|---|---|
| 1 | `ap4_validation.py` / `validate_diagram_data` | Feldverfügbarkeit, `(n,3)` bzw. `(n,m,3)`, endliche Werte, Punktzahlen, Q-Raster und kritische Koordinaten prüfen | Positives und negatives synthetisches Testset |
| 2 | `ap4_validation.py` / Diagrammzweige | PVT: `phase_mask`, fehlende kritische Q-Zeile, `v_l < v_g`; übrige: drei Teilflächen und Hochtemperatur-Anschluss | Fehlerfälle mit verschobenen Grenzen |
| 3 | `ap4_validation.py` / Q-Interpolation | Physikalische Q-Mischung: p-v-T Volumen, p-h-T Enthalpie, p-s-T Entropie, h-s-p beide Größen, u-v-T Volumen und Energie | Absichtlich verfälschte mittlere Q-Stützstelle wird abgelehnt |
| 4 | `test_ap4.py` / `AP4SyntheticTests` | Originale Achsentransformationen mittels AST aus `visualization.py` laden, ohne PyVista zu benötigen, kritisch normierte Koordinaten kontrollieren | Fünf Transformationsvergleiche in einem Test |
| 5 | `test_ap4.py` / `AP4HEOSTests._calculate` | Reale Daten für je zwei Fluide, verkleinerte reproduzierbare Raster für Laufzeitbegrenzung | Berechnung im Zielsystem vorgesehen |
| 6 | `test_ap4.py` / `_check_critical`, `_check_sat`, `_check_single_phase` | Pointwise-HEOS-Kontrolle von Tc, pc, h/s/u/v, Q=0/1 sowie repräsentativen Einphasenzuständen | HEOS-Tests im Zielsystem vorgesehen |
| 7 | `test_ap4.py` / `AP4RenderTests`, `AP4GUISmokeTests` | Offscreen-Renderer und opt-in GUI-Smoke mit realem PyVista-/Qt-Stack | Grafisches Zielsystem erforderlich |
| 8 | `README_AP4.md` | Ist/Soll, Schnittstellen, Tests, Einschränkungen und Abschlussstatus | Inhalts- und Syntaxprüfung |

### Warum keine neue Berechnungslogik?

Alle fünf Diagramm-Berechnungsfunktionen und die meisten Topologietests waren bereits implementiert. In AP4 ist **deren systematisches Testen** gefordert, nicht das erneute Entwickeln thermodynamischer Oberflächen. Eine darüber hinausgehende Produktivcodeänderung wäre ohne nachgewiesenen Defekt unnötig und würde das Regressionsrisiko erhöhen. Falls auf dem Zielrechner numerische Tests fehlschlagen, muss die konkrete Fehlerursache im betreffenden Modul behoben und anschließend die Matrix erneut ausgeführt werden.

## 5. Neue Funktionen, Parameter und Rückgaben

### `ap4_validation.validate_diagram_data(diagram_type: str, daten: dict) -> dict`

- **Eingaben:** Diagrammtyp (`pvt`, `pht`, `pst`, `hsp`, `uvt`), unverändertes Ergebnis-Dictionary aus dem zugehörigen Modul.
- **Ausgabe:** Dictionary `{"diagram": str, "binodal_points": int, "two_phase_points": int, "surface_points": int}`. Die Rückgabe ist eine Prüfstatistik, **keine Veränderung** des Berechnungsergebnisses.
- **Fehler:** `AP4ValidationError(ValueError)` bei fehlenden Feldern, nichtnumerischen/falschen Arraydimensionen, `NaN/Inf`, inkonsistenter Mischung, abweichenden Grenzknoten, falscher kritischer Spitze oder unplausibler Ordnung/positiven Druck- und Volumenbedingungen.
- **Abhängigkeiten:** Nur NumPy. Weder CoolProp noch Renderer oder GUI werden importiert.
- **Numerische Toleranz:** Exakte Arraygleichheit für explizit wiederverwendete gemeinsame Randpunkte. Für interpolierte mittlere Q-Stützstellen Standard `rtol=2e-8`, `atol=1e-10`; für Mischungsdruck `atol=1e-6 Pa`. Dies sind geometrisch-numerische Programmprüfungen, keine physikalischen Messunsicherheiten.
- **Hinweis:** Werte der Entropie, Enthalpie und inneren Energie werden *nicht* pauschal auf positive Zahlen eingeschränkt, da ihre Referenznullpunkte modellabhängig sind.

### `test_ap4.AP4SyntheticTests`

Sieben Tests ohne CoolProp/Qt/PyVista: gültige Minimalgeometrien, absichtlich verschobene Randknoten, kritische Fehler, ungültige Arrays, Q-Interpolationsfehler, tatsächliche Renderer-Transformationsformeln und der Qt-freie Import des originalen Berechnungsdispatchers mit Mock-Aufrufen. Die Fixtures sind mathematische Kunstwerte und werden **nicht** als HEOS-Stoffdaten ausgegeben.

### `test_ap4.AP4HEOSTests`

- `_calculate(diagram, fluid)` ruft für `pvt` die bestehende AP2/AP3-Berechnung in `main.berechne_daten()` mit verkleinertem, stoffbezogenem Raster auf; zum numerischen Test wird `main.py` mit einem isolierten Qt-/GUI-Importstub geladen, während die echten HEOS-Berechnungsmodule unverändert arbeiten; andere Diagramme werden direkt über ihre vorhandenen Berechnungsmodule mit `n_T_sub=8`, `n_T_super=4`, `n_p=9` beziehungsweise `n_v=13` und `n_quality=9` berechnet.
- `_check_critical(...)`: Vergleich mit einem frischen `AbstractState("HEOS", fluid)` bei kritischer molarer Dichte und Temperatur.
- `_check_sat(...)`: Drei unterkritische Temperaturen, jeweils Q=0 und Q=1, Vergleich mit frischen HEOS-Sättigungs-Flashes.
- `_check_single_phase(...)`: Direkte `PT_INPUTS`- bzw. `DmolarT_INPUTS`-Neuberechnung ausgewählter Einphasenzustände. Physikalische Einheiten: `p` in Pa, `h/u` in J/kg, `s` in J/(kg K), `v_m` in m³/mol.
- `_run_fluid_matrix(diagram)`: Führt Prüfbaugruppen für Water und R134a durch, inklusive vorhandener modulspezifischer Validatoren.
- `test_invalid_model_limits_two_fluids`: Verlangt `ValueError` bei Temperaturen außerhalb des zulässigen Zweiphasenbereichs und bei ungültigen Qualitätsrastern.
- **Ausgaben:** `unittest`-Ergebnisse mit `subTest` je Diagramm/Fluid. Keine Änderungen am gespeicherten Programmzustand. Reale Modellgrenzen müssen auf dem Zielsystem nachgewiesen werden.

### `test_ap4.AP4RenderTests` und `AP4GUISmokeTests`

- `AP4RenderTests` verwendet `pyvista.Plotter(off_screen=True)` und die tatsächlichen `render_*`-Funktionen; geprüft werden Oberfläche, Zweiphasenfläche, Binodalen und kritischer Marker.
- `AP4GUISmokeTests` erstellt das reale `ThermoSurfaceWindow`, wählt alle fünf Diagrammtypen, führt mit echten berechneten Daten jeweils `_calculation_ready()` aus, prüft `rendered_diagram_type` und schaltet den Oberflächen-Actor aus und wieder ein. Der direkte GUI-Renderpfad wird damit abgedeckt, nicht aber der vollständige asynchrone QThread-Lebenszyklus.
- Ohne geeignete Bibliotheken beziehungsweise ohne `THERMOSURFACE_AP4_GUI=1` werden diese Tests **als SKIP** protokolliert. Das ist keine bestandene GUI-Abnahme.

## 6. Architektur und Integration

```text
                Bestehendes ThermoSurface
        main.py  ───────>  coolprop_heos_* (5 Module)
           |                        |
           v                        v
         gui.py  <─────── Ergebnis-Dictionaries
           |
           v
    visualization.py (5 Renderer)

          AP4: separat ausführbare Tests
    test_ap4.py  ───────>  vorhandene Module
          |               (ohne Quellcodeänderung)
          +───────────>  ap4_validation.py (NumPy)
          +───────────>  bestehende Renderer / GUI
```

Im Produktionsbetrieb wird `ap4_validation.py` **nicht automatisch** bei jedem GUI-Aufruf ausgeführt. Das vermeidet redundante Prüfkosten und verändert das Verhalten der abgeschlossenen Arbeitspakete nicht. Es kann von Tests oder späteren Diagnosewerkzeugen explizit aufgerufen werden. Die bestehenden fünf Diagramm-Datenverträge, Actor-Namen und Oberflächen-/Binodalen-Arrays bleiben unverändert.

**Wichtig für die wissenschaftliche Einordnung:** Der Stichprobenvergleich nutzt auf beiden Seiten *dasselbe* CoolProp-HEOS-Modell mit verschiedenen Eingabepaaren und Datenpfaden. Er überprüft daher die **interne Konsistenz und Einheiten-/Schnittstellenkorrektheit**, nicht eine unabhängig gemessene Zustandsgleichungsgenauigkeit. Letzteres gehört zu AP6.

## 7. Tests und Validierung

### 7.1 Tatsächlich hier ausgeführt

Umgebung: Python **3.13.5**, NumPy **2.3.5**, PyYAML **6.0.3**. CoolProp, PyVista, PyVistaQt und PySide6 waren nicht installiert; ein Installationsversuch war wegen fehlender Netzwerk-Namensauflösung nicht erfolgreich.

| Test / Akzeptanzkriterium | Befehl / Methodik | Erwartetes Resultat | Tatsächlicher Befund | Status |
|---|---|---|---|---|
| Python-Syntax sämtlicher Projektmodule inkl. AP4 | `python -m compileall -q .` | Kein Syntaxfehler | Keine Syntaxfehler | **Bestanden** |
| AP4: fünf synthetische Diagrammformen und Fehlereinspeisung | `python -m unittest test_ap4.AP4SyntheticTests -v` | Sieben positive/negative Tests bestehen | **7 Tests bestanden** | **Bestanden (synthetisch)** |
| AP4: fünf Renderer-Achsentransformationen | In `AP4SyntheticTests` | Kritische Punkte in normierten Koordinaten und endliche Binodalen | Test enthalten in den sieben bestandenen Tests | **Bestanden (reine NumPy-Formeln)** |
| AP1: Snapshot-Comparator-Logik | `python -m unittest Test.RegressionTests -v` | Sechs isolierte Tests | **6 bestanden** | **Bestanden** |
| AP2/AP3: synthetische Rand-/Spinodalen-/Dispatcher-Prüfungen | `python -m unittest Test.AP3SyntheticTests -v` | Neun Tests | **9 bestanden** | **Bestanden** |
| AP4: 5 Diagramme × 2 reale Fluide | `python -m unittest test_ap4.AP4HEOSTests -v` | 10 Fälle, HEOS-Vergleich erfolgreich | CoolProp fehlt; 6 Klassen-Tests als **SKIP** | **Nicht ausführbar** |
| AP4: echter PyVista-Renderer | `python -m unittest test_ap4.AP4RenderTests -v` | Fünf Actor-Dictionaries fehlerfrei | PyVista/CoolProp fehlen; **SKIP** | **Nicht ausführbar** |
| AP4: echte Qt-Interaktion | `THERMOSURFACE_AP4_GUI=1` + GUI-Test | Fünf Diagrammwechsel/Sichtbarkeitszyklen | Qt-/VTK-Umgebung fehlt; **SKIP** | **Nicht ausführbar** |
| AP1: gespeicherte echte Vergleichs-Baselines | `Test.RecordedRegressionTests` | Vergleich echter Baseline und aktuellem Lauf | Verzeichnis `regression/` im bereitgestellten Projektstand nicht enthalten | **Nicht durchgeführt (bewusst fail-closed)** |

**Zusammenfassung des AP4-Testlaufs:** `python -m unittest test_ap4 -v` ergab **15 Tests: 7 erfolgreich, 8 ausdrücklich übersprungen, 0 fehlgeschlagen**. Ein grüner `unittest`-Exitcode bedeutet hier **nicht**, dass die acht übersprungenen echten HEOS-/Grafikprüfungen bestanden hätten.

### 7.2 Tests am Windows-Zielrechner durchführen

Im bestehenden Projektordner in PowerShell:

```powershell
# Ggf. bereits vorhandene Projekt-Virtualenv aktivieren.
python -m pip install numpy CoolProp PyYAML pyvista pyvistaqt PySide6

# Zunächst AP4-Synthetik, dann zehn echte HEOS-Fluidfälle.
python -m unittest test_ap4.AP4SyntheticTests -v
python -m unittest test_ap4.AP4HEOSTests -v

# Echte VTK-/Offscreen-Renderer (erfordert funktionsfähige Grafikumgebung).
python -m unittest test_ap4.AP4RenderTests -v

# Auf einem Rechner mit sichtbarer Desktop-Sitzung, keine Headless-Session:
$env:THERMOSURFACE_AP4_GUI = "1"
python -m unittest test_ap4.AP4GUISmokeTests -v
Remove-Item Env:\THERMOSURFACE_AP4_GUI

# Vorhandene AP1-/AP2-/AP3-Prüfungen zusätzlich:
python -m unittest Test.RegressionTests Test.AP3SyntheticTests -v
python -m unittest Test.ThermoSurfaceTests Test.AP3HEOSTests -v

# Manuelle grafische Endkontrolle mit beiden Prüffluiden:
python main.py
```

Auf dem Zielrechner die tatsächlichen Testausgaben speichern. **Für die Abnahme müssen zehn reale Diagramm-Fluidfälle einzeln als bestanden nachweisbar sein.** Fehlgeschlagene Tests führen nicht zum Abbruch der Analyse, sondern sind als konkrete Codefehler zu untersuchen und nach Korrektur erneut auszuführen. Ist die PyVista-Offscreen-Ausgabe unter einer bestimmten Grafiktreiber-/Remote-Desktop-Konfiguration nicht möglich, diese Einschränkung protokollieren und den Renderpfad zusätzlich direkt im sichtbaren GUI-Fenster testen.

### 7.3 Manuelle GUI-Abnahme (zusätzlich zum automatisierten Smoke-Test)

Für **Water und R134a** jeweils alle fünf Diagrammtypen auswählen und die Berechnung über **„Berechnung starten“** auslösen. Verifizieren, dass die QThread-Berechnung beendet wird, die Fortschrittsanzeige 100 % erreicht, kein Fehlerdialog erscheint und sich bei jedem Typ jeweils **Zustandsoberfläche**, **Zweiphasengebiet**, **Siede-/Taulinie** und **kritischer Punkt** getrennt schalten lassen. Kameradrehung, Gitternetz-Checkbox und Diagrammwechsel ebenfalls prüfen. Im p-v-T-Diagramm AP2-Randkontinuität und AP3-Warnstatus separat beachten. Die dokumentierten Tests ergänzen, ersetzen aber nicht diese visuelle Kontrolle.

### 7.4 Numerische Prüfkriterien

- Randkontinuität: **bitgenau gleiche Knoten**, da dieselben Datenarrays zur Vernetzung und zur Phasengrenzlinie benutzt werden.
- Q-Interpolation: Vergleich mit der linearen Massendampfgehaltsmischung der jeweiligen extensiven spezifischen Größe; `Q=0`, `Q=1`, innere Werte und kritische Kollapszeile.
- Reale HEOS-Stichproben: relative Abweichung höchstens `2e-6` in den fachlich vergleichbaren Größen, plus absolute Toleranzen entsprechend der Einheit (typisch `1e-6` bei Sättigungsdaten, `1e-3` bei h/s/u und `0.5 Pa` beim einphasigen p-v-T-Druck). Diese Grenzen sind bewusst **numerische Tests**, keine Messunsicherheiten.
- Kritischer Punkt: von `rhomolar_critical` und `T_critical` frisch berechnete Zustandsgrößen; `rtol=2e-7`.
- Raster-/Gültigkeitsgrenzen: unzulässige Temperaturen und `n_quality=1` müssen durch `ValueError` zurückgewiesen werden.

## 8. Einschränkungen und offene Risiken

1. **Entscheidend:** Der wissenschaftlich relevante AP4-Nachweis mit echtem CoolProp HEOS auf zwei Fluiden ist hier **offen**. Die implementierten numerischen Tests wurden nicht mit realen Stoffdaten ausgeführt.
2. Renderer und GUI wurden hier nicht grafisch ausgeführt. Die AST-Transformationsprüfung deckt nur die reinen Zahlenformeln ab, nicht VTK-Mesh-Erstellung, Triangulation, Grafikkartentreiber oder Qt-Signalverbindungen.
3. `AP4GUISmokeTests` prüft den Aufruf von `_calculation_ready()` auf einem echten Qt-Fenster, aber nicht den vollständigen QThread-Zeitablauf; dieser ist manuell zu prüfen.
4. Bei fehlenden optionalen Abhängigkeiten **skipt** die AP4-Suite bewusst. Für eine endgültige Abnahme zählt die Zahl erfolgreicher **real ausgeführter** Prüfungen, nicht nur der Prozess-Exitcode.
5. PVT kann fluidabhängig nicht verfügbare oder unvollständige Spinodaläste liefern. Dies wird durch AP3 kenntlich gemacht und darf nicht als bestandene Spinodalresidualvalidierung fehlinterpretiert werden.
6. Größere Grenzfälle und Rasterkonvergenz werden hier nicht systematisch profiliert (AP8); unabhängige Literaturdaten fehlen bewusst (AP6).
7. AP1-Baselines (`regression/`) sind nicht Bestandteil der bereitgestellten Projektdateien; eine echte Snapshot-Regression kann daher nicht nachgewiesen werden. Die **sechs isolierten Regressionstests** wurden davon unabhängig durchgeführt.
8. Die Zusatzdiagramm-Fluidergebnisse aus dem Zielsystem können bislang nicht auf fluidabhängige HEOS-Fehler hin korrigiert werden; etwaige Failures sind **offen zu behandeln, nicht als Erfolg zu deklarieren**.

## 9. Abschlussstatus AP4

| Verbindliche AP4-Anforderung | Implementiert | Validiert in dieser Umgebung | Bemerkung |
|---|---|---|---|
| Numerische Prüfung p-v-T mit ≥2 Fluiden | Ja (Testcode) | **Nein, SKIP** | Water/R134a stehen aus |
| Numerische Prüfung p-h-T mit ≥2 Fluiden | Ja (Testcode) | **Nein, SKIP** | Water/R134a stehen aus |
| Numerische Prüfung p-s-T mit ≥2 Fluiden | Ja (Testcode) | **Nein, SKIP** | Water/R134a stehen aus |
| Numerische Prüfung h-s-p mit ≥2 Fluiden | Ja (Testcode) | **Nein, SKIP** | Water/R134a stehen aus |
| Numerische Prüfung u-v-T mit ≥2 Fluiden | Ja (Testcode) | **Nein, SKIP** | Water/R134a stehen aus |
| Geometrie-, Q-, Tc- und Achsenprüfungen | Ja | **Ja, synthetisch** | Sieben Tests erfolgreich |
| Reale HEOS-Gültigkeits-/Phasenrandfälle | Ja (Testcode) | **Nein, SKIP** | Temperatur- und Qualitätsraster im Zielsystem prüfen |
| Echter Renderer jedes Typs | Ja (Smoke-Test) | **Nein, SKIP** | VTK-/PyVista-Zielsystem erforderlich |
| Echter GUI-Pfad | Ja (Smoke-Test + Checkliste) | **Nein, SKIP** | Desktop-Zielsystem erforderlich |
| AP1–AP3-Schutz | Keine Altdateien verändert | **15 isolierte Tests bestanden** | Echte Baseline-/HEOS-Regressionsläufe offen |
| Technische Dokumentation | Ja | Ja, anhand tatsächlich vorhandener Änderungen | Diese README |

**AP4 ist softwareseitig als Test- und Prüfpaket implementiert, aber bis zur erfolgreichen Ausführung der realen Fluid- und Grafiktests noch nicht fachlich endgültig abgenommen.** Diese Unterscheidung ist für die Studienarbeit zwingend.

## 10. Voraussetzungen für AP5 bis AP10 und Update-Archiv

- **AP5:** Die vorhandenen `render_pht`- und GUI-Actor-Schnittstellen bleiben unangetastet. Isolinien können später separat ergänzt werden.
- **AP6:** Unabhängige experimentelle/literaturbasierte Stoffdaten sind gesondert einzubinden. Die hier ausgeführten Selbstkonsistenzvergleiche dürfen nicht als solche ausgegeben werden.
- **AP7:** EOS-Vergleich erfordert keine Umstellung der HEOS-Produktivmodule.
- **AP8:** Für Rasterkonvergenz/Benchmarking später separate reproduzierbare Gittergrößen erfassen; die AP4-Prüfungen arbeiten zur Zeitbegrenzung bewusst auf reduzierten Netzen.
- **AP9:** Die numerischen Soll-/Ist-Protokolle bilden einen späteren dokumentierbaren Teil der Methoden-/Validierungskapitel.
- **AP10:** Später vollständige Versionen und echten Regression-Snapshot erfassen; keine Baseline aus dem Nichts erzeugen.

### Inhalt des AP4-Update-Pakets

```text
ThermoSurface_AP4.zip
├── ap4_validation.py   # Neu: reiner NumPy-Datenvalidator
├── test_ap4.py         # Neu: synthetische, HEOS-, Renderer- und GUI-Tests
└── README_AP4.md       # Neu: Ist/Soll, Implementierung, Test- und Abnahmeprotokoll
```

Die drei Dateien in den **bestehenden** ThermoSurface-Projektordner kopieren; vorhandene Module und deren Dateinamen beibehalten. Zum Start der normalen Anwendung weiterhin `python main.py` verwenden. Die Projektanalyse, die PDF-Masterarbeit und bestehende AP1-/AP2-/AP3-Dateien werden nicht ersetzt.
