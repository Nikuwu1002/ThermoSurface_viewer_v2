# ThermoSurface – AP3: Validierung der Spinodalen

**Stand:** 19. Oktober 2026  
**Paket:** AP3 („Spinodalen validieren“) gemäß `ThermoSurface_Projektanalyse_und_Arbeitspakete_AP1-AP10.md`, Abschnitt 5.3 und AP3-Abnahmekriterium  
**Lieferform:** Update-Paket, **kein vollständiger Projektstand**. Die Dateien aus dem ZIP in den bestehenden ThermoSurface-Projektordner kopieren. Andere Projektdateien einschließlich `config.yaml` unverändert weiterverwenden.

## 1. Ziel, Umfang und Abgrenzung

AP3 überprüft die beiden mit CoolProp HEOS ermittelten p-v-T-Spinodaläste quantitativ. Die Bedingungen sind:

- normiertes dimensionsloses Ableitungsresiduum der **homogenen** Zustandsgleichung,
- konsistente Phasenvorgabe der Dichte-Temperatur-Zustände,
- Eingrenzung auf gültige Temperaturstützstellen und Einhaltung des GUI-Temperaturfilters,
- begründeter statt erzwungener Anschluss am kritischen Punkt,
- Kennzeichnung fehlender, lückenhafter, nicht darstellbarer oder unzureichend geprüfter Äste.

**Unverändert:** AP1-Snapshot- und Vergleichsinfrastruktur; AP2-Raster, Binodal- und Mesh-Randkontinuität; Berechnungs- und Rendererlogik für `pht`, `pst`, `hsp`, `uvt`; die Trennung GUI ↔ `main.py` ↔ HEOS ↔ Visualisierung.

**Nicht Teil von AP3:** unabhängige Stoffdatenreferenzen (AP6), neue metastabile/instabile 3D-Flächen, Anpassung der anderen vier Diagrammtypen (AP4/AP5), Performance-Benchmark (AP8) oder weitere EOS (AP7).

## 2. Ausgangszustand und Anforderungsanalyse

Die vorherige `calc_spinodale(stoff)` nutzte `AbstractState.build_spinodal()` und `get_spinodal_data()` und rekonstruierte Punkte `(T, v, p)` aus reduzierter Temperatur und Dichte. Für die homogenen Druckauswertungen wurde eine Phase vorgeschrieben. Schwierigkeiten:

1. Ein Fehler beim homogenen Druck-Flash erzeugte `NaN`, wurde aber nicht quantitativ in einer Freigabe erfasst.
2. Der alte Test `Debug.test_spinodale_ableitung()` protokollierte ausschließlich Median und Maximum der normierten Ableitung, ohne Toleranz/Status/Rückgabe.
3. Sowohl `calc_spinodale()` als auch `main.py` ergänzten künstlich den kritischen Punkt. `main.py` tat dies **nach** dem Temperaturfilter – also gegebenenfalls außerhalb von `T_stop_spinodal`.
4. Ein Ast ohne reale Unterkritisch-Stützpunkte konnte durch das Anfügen eines kritischen Punkts fälschlich geschlossen erscheinen.
5. Spinodaldrücke `p <= 0` sind auf der verwendeten `log10(p/pc)`-Achse nicht darstellbar; die bisherige GUI hatte keine explizite AP3-Statusanzeige.
6. Phasenrichtung, unvollständige Temperaturabdeckung, größere Stützstellenlücken und fehlgeschlagene Ableitungsberechnungen hatten keine gemeinsame Diagnose.

| AP3-Anforderung | Vor AP3 | Änderung | Dateien | Priorität |
|---|---|---|---|---|
| Ableitungsresiduum mit Toleranz | Nur Logging | Explizites Prüfergebnis und Stichprobenzahlen | `Debug.py`, `main.py` | Hoch |
| Homogene Phasenvorgaben | Vorhanden, nicht testbar abgesichert | Beibehalten, Phasenzuordnung separat prüfen | `coolprop_heos_pvt.py`, `Test.py` | Hoch |
| Temperaturfilter | Vorhanden, nachträglich unterlaufen | Grenzen validieren, keine unzulässige Tc-Ergänzung | `coolprop_heos_pvt.py`, `main.py` | Hoch |
| Kritischer Endpunkt | Unbedingt angehängt | Nur innerhalb Filter und mit Annäherungsheuristik | `coolprop_heos_pvt.py`, `main.py` | Hoch |
| Fehlende/unvollständige Äste | Kaum kenntlich | Strukturbericht und Ableitungsstatus, Warnhinweis in GUI | `coolprop_heos_pvt.py`, `main.py`, `gui.py` | Hoch |
| AP1/AP2 erhalten | Abgeschlossen | Synthetische Erhaltungs- und Integrationstests | `Test.py` | Hoch |

## 3. Tatsächlich geänderte Dateien

| Datei | Funktion / Codebereich | Vorher | Nachher | Begründung |
|---|---|---|---|---|
| `coolprop_heos_pvt.py` | `calc_spinodale()` | Rohe Punkte plus immer angehängter kritischer Punkt | Rohe, sortierte unterkritische Punkte; gleichförmige `(n,3)`-Arrays; NaN bleibt sichtbar | Keine vorgetäuschte Kurvenvollständigkeit |
| `coolprop_heos_pvt.py` | `filter_spinodale()` | Stilles Filtern | Temperaturgrenzen und Arrayform geprüft | Filter kann nicht unbemerkt überschritten werden |
| `coolprop_heos_pvt.py` | neu: `prepare_spinodal_branch()` | Fehlte | Astgeometrie, Phase, Druck, Temperaturabdeckung und Tc-Anschluss beurteilen | Reproduzierbare AP3-Strukturdiagnose |
| `Debug.py` | `test_spinodale_ableitung()` | Nur Median/Max-Log | Normiertes Residuum, Grenzwert, Zähler, expliziter Pass/Fail; optionaler Mock-Zustand | Quantitative und testbare Thermodynamikdiagnose |
| `main.py` | `berechne_daten()` | Kritischer Punkt nach Filter forciert; Ableitung nur im Debug-Modus | Spinodalen getrennt vorbereiten und **immer** prüfen; `spinodal_validation` hinzufügen | Integrierte Freigabe im Normalbetrieb |
| `gui.py` | `_calculation_ready()`, `_thread_finished()` | Allgemeiner Abschlussstatus, durch Threadabschluss zurückgesetzt | Warnhinweis und Tooltip; Abschlussstatus beim Threadende bewahren | Nutzer sieht Einschränkungen dauerhaft |
| `Test.py` | `AP3SyntheticTests`, `AP3HEOSTests` | Fehlten | Reine NumPy-/Mocktests und separater echter HEOS-Zielrechnertest | Automatisierter Prüfnachweis |

Unveränderte Python-Module werden **nicht** im Update-Paket überschrieben.

## 4. Implementierungsplan und durchgeführte Reihenfolge

| Schritt | Datei | Funktion | Umsetzung | Abhängigkeit | Prüfung |
|---|---|---|---|---|---|
| 1 | `coolprop_heos_pvt.py` | `calc_spinodale` | Kritische Kunstpunkte aus Rohberechnung entfernen, zwei sortierte Äste erhalten | CoolProp HEOS | Mock-HEOS-Stützstellen |
| 2 | `coolprop_heos_pvt.py` | `filter_spinodale` | Bereich und Eingabeformat absichern | NumPy | Grenz-/Leertests |
| 3 | `coolprop_heos_pvt.py` | `prepare_spinodal_branch` | Kritischen Anschluss sowie Sichtbarkeits-/Plausibilitätsstatus bestimmen | `filter_spinodale` | Gültige, fehlende und unplottbare Äste |
| 4 | `Debug.py` | `test_spinodale_ableitung` | Residuen mit Grenzwert und Fehlerrückgabe bewerten | HEOS `first_partial_deriv` | Fake-State für Pass, Fail und Exceptions |
| 5 | `main.py` | `berechne_daten` | Beide Berichte zusammenführen und in Ergebnissen ausgeben | Schritte 1–4 | Mock-Integration des realen Dispatchers |
| 6 | `gui.py` | `_calculation_ready`, `_thread_finished` | Bericht in Warnstatus und Tooltip überführen; Status nach Threadende erhalten | `main.py` | Statische Integration; echter GUI-Test offen |
| 7 | `Test.py` | neue Tests | AP1/AP2-Schutz, Berechnungs- und Integrationstests | Schritte 1–6 | `unittest` |
| 8 | `README_AP3.md` | Dokumentation | Entscheidungen, Prüfergebnisse, Grenzen | Vorherige Schritte | Inhaltskontrolle |

### 4.1 Spinodalpunkte und Phasenzuordnung

`calc_spinodale(stoff)` verwendet weiterhin `build_spinodal()` und `get_spinodal_data()` und rechnet mit

- `T = T_reducing / tau` [K],
- `rho_molar = delta * rho_reducing` [mol/m³],
- `v = 1 / rho_molar` [m³/mol],
- `p = AbstractState.p()` [Pa] nach `specify_phase(iphase_liquid/iphase_gas)` und `DmolarT_INPUTS`.

`rho_molar > rho_c` bezeichnet den flüssigkeitsseitigen Ast, `rho_molar < rho_c` den dampfseitigen Ast. Für die Berechnung fehlgeschlagener Drücke bleibt `NaN` im Array: Die Darstellung trennt entsprechende Segmente, die AP3-Diagnostik meldet sie als ungültig. Bereits durch HEOS nicht gelieferte Bereiche können nicht automatisch rekonstruiert werden.

Die Funktionssignatur und das **Zwei-Array-Rückgabeformat** bleiben erhalten. **Bewusste Verhaltensänderung:** Die Rohfunktion hängt nicht mehr selbst `Tc` an. Die einzige bedingte Endpunktergänzung liegt in `prepare_spinodal_branch()`.

### 4.2 Temperatur-, Phasen- und Endpunktkontrolle

`prepare_spinodal_branch(kurve, critical, T_min, T_max, branch)` erwartet Rohpunkte `(T,v,p)`, `(Tc,pc,vc)` und den Astnamen `liquid` oder `vapor` und liefert `(points, report)`.

Geprüft werden mindestens:

- `0 < T_min < T_max`, endliche Temperaturen, streng aufsteigende Temperaturpunkte,
- mindestens drei unterkritische Stützstellen für vollständige geometrische Abdeckung,
- `v > 0`, `T > 0`, `p` endlich; für die logarithmische Grafik zusätzlich `p > 0`,
- richtige Seite des kritischen molaren Volumens (`v < vc` flüssig, `v > vc` dampfförmig),
- Anfang und Ende des angeforderten Temperaturintervalls innerhalb einer stützstellenabhängigen Toleranz,
- auffallend große Lücken in der Temperaturabfolge.

**Kritischer Anschluss:** `(Tc,vc,pc)` wird nur angefügt, wenn `T_min <= Tc <= T_max`, mindestens zwei tatsächliche Punkte vorhanden sind, der oberste unterkritische Punkt um höchstens `0.01*Tc` vom kritischen Temperaturwert abweicht und das relative Volumen weniger als 50 % vom kritischen Volumen entfernt liegt. Dies ist ausdrücklich eine **geometrische Heuristik**, kein EOS-Residualbeweis für den kritischen Punkt. Der kritische Punkt wird nicht in den Ableitungstest aufgenommen.

Ein negativer homogene-EOS-Druck ist nicht prinzipiell ein numerischer Fehler, aber in der bestehenden `log10(p/pc)`-Darstellung unsichtbar. Er wird gemeldet und verhindert die Einstufung **„geometrisch vollständig darstellbar“**. Für echte Brüche in der Kurve werden keine fiktiven Punkte interpoliert.

`report` enthält `points`, `subcritical_points`, `invalid_states`, `nonpositive_pressure`, `wrong_phase`, `covers_start`, `covers_end`, `critical_appended`, `temperature_order_ok`, `large_gaps`, `geometry_complete`.

### 4.3 Thermodynamisches Spinodalresiduum

Die physikalische Spinodalbedingung lautet:

\[
\left(\frac{\partial p}{\partial v}\right)_T=0
\quad\Longleftrightarrow\quad
\left(\frac{\partial p}{\partial\rho_m}\right)_T=0,
\qquad \rho_m=\frac1v.
\]

Das implementierte dimensionslose Residuum ist:

\[
\boxed{R_s=\left|\frac{\rho_c}{p_c}\left(\frac{\partial p}{\partial\rho_m}\right)_T\right|}
\]

Die Ableitung kommt aus `state.first_partial_deriv(CP.iP, CP.iDmolar, CP.iT)` nach `state.specify_phase(...)`, anschließend wird die Vorgabe immer zurückgesetzt. Der numerische Zustand ist **homogen**, nicht der stabile Zweiphasengleichgewichtszustand.

**Verwendete Akzeptanzparameter:**

- `tolerance = 1e-2` (maximaler dimensionsloser Residualbetrag höchstens 1 %), als pragmatischer Startwert zur quantitativen Prüfung. Der Wert ist **kein universeller Referenzgrenzwert** und muss für konkrete Fluidversionen in der echten HEOS-Umgebung hinterfragt werden.
- `critical_exclusion_rel = 1e-3`, zusätzlich mindestens `0.5 K` Ausschluss unmittelbar unter `Tc`, um numerisch empfindliche kritische Zustände und den expliziten Tc-Knoten auszuschließen.
- `max_samples = 128` pro Ast, gleichmäßig über die **für den Test geeigneten** Stützstellen verteilt. Die Prüfung ist ausdrücklich **stichprobenartig**, nicht flächendeckend.
- Mindestens zwei erfolgreich getestete Zustände pro Ast und **keine** fehlgeschlagenen Ableitungen in den gewählten Stichproben.

`Debug.test_spinodale_ableitung()` gibt für jeden Ast `candidates`, `sampled`, `excluded`, `evaluated`, `errors`, `max_residual`, `median_residual`, `tolerance`, `passed` zurück. `passed` ist bei zu wenigen Daten, zu hohem Residuum oder mindestens einem Berechnungsfehler `False`. Die Funktion akzeptiert für Tests optionale `state_factory`- und `cp_module`-Argumente; im Produktivbetrieb werden die echten CoolProp-Objekte verwendet.

**Abschlussregel je Ast:** `complete = geometry_complete and derivative.passed`. Die Bewertung gilt ausschließlich für den eingestellten sichtbaren T-Bereich und die getesteten Stützstellen. Sie ist keine Aussage über das Verhalten jedes thermodynamisch möglichen Zustands.

### 4.4 Anwendungsarchitektur und Resultate

```text
GUI (gui.py)
   └── main.berechne_daten() / main.berechne_auswahl('pvt')
         ├── vorhandene p-v-T-/AP2-Berechnungen (unverändert)
         ├── pvt.calc_spinodale(stoff) -> (raw_liquid, raw_vapor)
         ├── pvt.prepare_spinodal_branch(...) -> (gefilterte Kurve, Geometriebericht)
         ├── Debug.test_spinodale_ableitung(...) -> Residuen und Pass/Fail
         └── daten['spinodal_validation'] -> GUI-Status/Tooltip
                └── bestehendes PyVista render_pvt() (unverändert)
```

Das Ergebnis-Dictionary besitzt weiterhin die unveränderten Schlüssel `spinodale_liquid`, `spinodale_vapor`, `binodal_liquid`, `binodal_vapor`, `critical_point`, `surface` usw. **Hinzugekommen ist ein optionaler Metadatenschlüssel `spinodal_validation`**, mit den Unterwörterbüchern `liquid` und `vapor`. Für die bestehende PyVista-Darstellung oder die AP1-Snapshot-Arrayliste entsteht kein neues Arrayformat.

Auch bei `DEBUG_MODE = False` werden die Ableitungen überprüft; `DEBUG_MODE` steuert nur die zusätzliche Diagnostik. Eine fehlerhafte Spinodale blockiert nicht automatisch die gültigen AP2-Einphasen- und Binodalflächen. Die GUI gibt stattdessen ausdrücklich einen Warnstatus aus. Für die anderen Diagrammtypen wird die Spinodalprüfung weiterhin nicht gestartet.

## 5. Ausgeführte Tests und tatsächliche Ergebnisse

**Lokale Prüfungen in dieser Bearbeitung:** Python 3.13; NumPy und Standardbibliothek verfügbar. CoolProp, PyVista, PyVistaQt und PySide6 **nicht verfügbar**. Ein Versuch, CoolProp mittels `pip` zu installieren, schlug wegen fehlender Namensauflösung für den Paketindex fehl. Keine numerischen Real-HEOS-Ergebnisse werden daher behauptet.

Ausgeführter Befehl im Projekt-Testordner:

```powershell
python -m unittest Test.AP3SyntheticTests Test.RegressionTests Test.AP3HEOSTests -v
```

**Ergebnis:** 16 Testfälle insgesamt: **15 bestanden, 1 übersprungen**, 0 fehlgeschlagen. Darin enthalten sind 9 reine AP3/AP2-Mock- und Strukturtests sowie 6 bestehende AP1-Regressionstests. `test_real_spinodal_water` wird nur bei tatsächlich installiertem CoolProp ausgeführt.

| AP3-Akzeptanzkriterium | Test / Prüfmethode | Erwartung | Ergebnis hier | Status |
|---|---|---|---|---|
| Normierte Ableitung und begründete Toleranz | Fake-HEOS: Null, zu großes Residuum, geworfene Exception | Pass / Fail / Fail | Entsprechend klassifiziert, einschließlich Rücksetzen der Phasenvorgabe | Bestanden (Mock) |
| Kritischer Ausschluss | Kunstpunkte bei `Tc-0.2 K` und `Tc` | Nicht als Residualbeweis zählen | 0 geprüfte Punkte, `passed=False` | Bestanden (Mock) |
| Richtige Phasenseite | Falsches `v/vc` in Flüssigkeitsast | Kennzeichnung `wrong_phase` | 1 erkannt, geometrisch unvollständig | Bestanden (synthetisch) |
| Temperaturgrenzen | Kurve bei `T_max < Tc` | Kein Tc-Knoten außerhalb Fenster | Kurve bleibt innerhalb des Fensters | Bestanden (synthetisch) |
| Kritischer Endpunkt | Nahes und fernes Kurvenende | Nur nahes zulässig | Erwartete bedingte Ergänzung | Bestanden (synthetisch) |
| Fehlende/ungültige/negative Drücke | Leerer Ast, `NaN`, negativer Druck | Keine vollständige Freigabe | Entsprechende Diagnosen vorhanden | Bestanden (synthetisch) |
| HEOS-Rohdaten, Sortierung, Phase | Fake-`get_spinodal_data`, fake-Zustands-Flash | Zwei unterkritische sortierte Äste ohne Tc | Entspricht Erwartung | Bestanden (Mock) |
| Integration `main.berechne_daten` | Produktiver Dispatcher mit gemockten Berechnungsteilen | AP2-Tc bleibt, AP3-Filter respektiert Grenzen; Report enthalten | Alle Assertions bestanden | Bestanden (Mock) |
| AP1-Beibehaltung | `Test.RegressionTests` | 6 Tests grün | 6/6 bestanden | Bestanden |
| AP2-Raster und Q-Ränder | Gemeinsames Raster und Zweiphasenfunktion synthetisch | Identische Sat-Ränder | Bestanden | Bestanden (synthetisch) |
| Physikalisches Residuum für Wasser | `Test.AP3HEOSTests.test_real_spinodal_water` | Reale HEOS-Stützstellen und Derivate auswertbar | CoolProp fehlt | **Nicht ausführbar** |
| Vollständiger GUI-/PyVista-Pfad | GUI starten, Wasser/Argon auswählen, Warnstatus/Tooltip prüfen | Sinnvolle Linien und Status | Qt/VTK fehlt | **Nicht ausführbar** |
| Weitere Diagrammtypen real berechnen | `pht/pst/hsp/uvt` am Zielrechner | Bestehende Berechnungen unverändert | HEOS + grafische Laufzeit fehlt | **Nicht ausführbar** |
| Persistente AP1-Snapshot-Baselines | `Test.RecordedRegressionTests` | Frühere echte Baselines vergleichen | `regression/` nicht mitgeliefert | **Nicht ausführbar** (bewusst kein Fake) |

Zusätzlich wurde die Python-Syntax der geänderten Dateien per `compileall` und `ast.parse` erfolgreich geprüft. Statische Prüfung der Aufrufstellen zeigt: nur die p-v-T-Verarbeitung verwendet die geänderte Rückgabe von `Debug.test_spinodale_ableitung()`, deren ursprüngliche vier Positionsargumente unverändert bleiben.

**Hinweis zu echten Regressionen:** Da AP3 absichtlich Spinodallängen, Rohpunktdaten und kritische Endpunkte verändert, können alte `.npz`-Snapshots für die Spinodalen vom neuen Stand abweichen. Die bestehenden AP1-Baselines **nicht überschreiben**; Unterschiede fachlich prüfen und erst dann neue freigegebene Referenzen anlegen.

## 6. Einrichtung und Zielrechner-Prüfung

Im bestehenden Projektordner müssen CoolProp, NumPy, PyYAML, PySide6, PyVista und PyVistaQt passend zum verwendeten Python installiert sein. Anschließend:

```powershell
# Ohne Qt/VTK, mit NumPy: reine AP1-/AP3-Synthetiktests
python -m unittest Test.AP3SyntheticTests Test.RegressionTests -v

# Mit installiertem CoolProp, keine grafische Umgebung erforderlich:
python -m unittest Test.AP3HEOSTests -v

# Nur mit vollständig installierter ThermoSurface-Laufzeit:
python -m unittest Test.ThermoSurfaceTests -v
python main.py
```

Am Zielrechner insbesondere für **Wasser und Argon** je einen normal eingestellten und einen absichtlich verkürzten Spinodal-Temperaturbereich auswählen, `T_stop_spinodal` unter und anschließend über `Tc` setzen, im GUI-Tooltip/Status den AP3-Bericht prüfen und auf gebrochene Linien bei nicht-positivem Druck achten. Im Debug-Log den Betrag `max_residual` gegen `tolerance` prüfen. Keine Freigabe auf Grundlage einer rein visuellen Kurvenform.

Die numerische Residualtoleranz und die geometrischen Kriterien können durch gemessene Real-HEOS-Ergebnisse weiter kalibriert werden. Änderungen an solchen Schwellenwerten gehören in dokumentierte, reproduzierbare Tests.

## 7. Bekannte Grenzen, Risiken und Restarbeiten

1. **Kein lokaler Real-HEOS-Test möglich:** Die thermodynamische Eignung von 1 % Toleranz, die Abdeckung beider Spinodaläste sowie kritische HEOS-Fortsetzungen bleiben bis zur Zielrechnerprüfung offen.
2. `build_spinodal()` kann fluidabhängig unvollständige Zweige liefern; AP3 **kennzeichnet**, aber extrapoliert fehlende Punkte nicht. Ein fehlender Bereich darf nicht als EOS-Datenlage interpretiert werden.
3. Bei logarithmischem Druck sind `p <= 0` nicht darstellbar. Eine andere Projektion wäre eine eigenständige Darstellungsentscheidung, die AP3 nicht vorwegnimmt.
4. Die geometrische Tc-Heuristik (1 % Tc, 50 % vc) und die Prüfung großer Temperaturlücken sind **plausibilitätsorientiert**, nicht mathematische Beweise eines kritischen Grenzwerts.
5. Das Residuum wird auf maximal 128 geeigneten Stützstellen pro Ast geprüft. Nicht ausgewählte Punkte bleiben ungetestet; `candidates` und `sampled` dokumentieren dies.
6. Andere drei-Diagrammtypen, Powerwall-Ausgabe und echte GUI-Interaktionen wurden nicht an der hier verfügbaren Laufzeitumgebung geprüft.
7. Der bestehende AP1-Snapshot-Vergleich benötigt die echten Baselines; diese liegen nicht vor und werden nicht erfunden.

## 8. Abschlussstatus AP3

| Verbindliche Anforderung | Implementiert? | Validiert in dieser Umgebung? | Status |
|---|---|---|---|
| Spinodal-Ableitungsresiduum mit Grenzwert | Ja | Synthetisch/Mock, nicht reale HEOS-Zustände | Technisch umgesetzt, wissenschaftliche Realprüfung offen |
| Phasen-/Einheiten-/Plausibilitätsprüfung | Ja | Synthetisch/Mock | Technisch umgesetzt |
| Temperaturfilter und kritischer Endpunkt | Ja | Synthetisch + Integrationsmock | Technisch umgesetzt |
| Kennzeichnung unvollständiger Äste | Ja | Datenbericht getestet; GUI statisch geprüft | Echte GUI-Prüfung offen |
| Keine AP1/AP2-Regressionen | Bestehende Logik beibehalten | AP1 synthetisch und AP2 Teilprüfungen | Vollständige Zielrechner-Regression offen |
| Fünf Diagrammtypen weiterhin betriebsbereit | Module/Schnittstellen unverändert | Nicht real ausführbar | Echte Funktionsfreigabe offen |

**Fazit:** AP3 ist im Quellcode implementiert und durch **15 ausführbare lokale Tests** (einschließlich AP1-Tests) abgesichert. Der wissenschaftliche AP3-Abschluss ist **noch nicht uneingeschränkt nachgewiesen**, weil reale HEOS-Ableitungsresiduen und die GUI hier nicht getestet werden konnten.

## 9. Übergabe an die nächsten Arbeitspakete

- **AP4:** Die Zusatzdiagramme bleiben unverändert; beim Zielrechnertest die neue `spinodal_validation` für p-v-T separat auswerten.
- **AP6:** Referenzdaten und gegebenenfalls strengere fluidabhängige Toleranzen für die Spinodale ergänzen; dies erfolgt ausdrücklich später.
- **AP8/AP10:** Gezielte Messung des Overheads der zusätzlichen HEOS-Ableitungsstichproben; die dokumentierte Testmatrix reproduzieren und echte Baselines einbinden.

## 10. Inhalt und Einspielen des Update-Pakets

```text
ThermoSurface_AP3.zip
├── Debug.py
├── coolprop_heos_pvt.py
├── main.py
├── gui.py
├── Test.py
└── README_AP3.md
```

Alle übrigen Dateien einschließlich `coolprop_heos_common.py`, `visualization.py`, `coolprop_heos_pht.py`, `coolprop_heos_pst.py`, `coolprop_heos_hsp.py`, `coolprop_heos_uvt.py`, `limits.py`, `config.yaml`, `compare_pvt.py` und `snapshot_pvt.py` aus der **bestehenden** Projektversion weiterverwenden. Die Projektbeschreibung, PDF-Referenzarbeit und privaten lokalen Dateien sind nicht Teil des Update-ZIP.
