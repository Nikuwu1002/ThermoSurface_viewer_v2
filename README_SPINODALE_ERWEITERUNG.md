# ThermoSurface – Erweiterte HEOS-Dampfspinodale

**Projektstand:** 19.10.2026  
**Umfang:** eigenständige Ergänzung zu AP1–AP4; kein Ersatz für AP5–AP10.  
**Implementierung:** vorhanden, **numerische Abnahme mit echtem HEOS: offen**.  
**Ausgangsbasis:** beigefügte Python-Quelldateien, `README_AP3.md`, `README_AP4.md`, `ThermoSurface_Spinodalen_Gesamtanalyse_Tests_Ergebnisse_2026-10-19.md`.

## 1. Ziel und Ausgangsproblem

Die bestehende CoolProp-Routine `calc_spinodale(stoff)` benutzt `AbstractState.build_spinodal()` und `get_spinodal_data()`. Bei Wasser zeigt der frühere lokale Windows-Test für den Dampfast erst ab **625,128101 K** einen positiven Druck; das GUI-Intervall beginnt bei 580 K. Eine zusätzliche Oberflächenauflösung behebt fehlende Daten nicht. Auf der flüssigkeitsseitigen Spinodale können physikalisch berechnete **nichtpositive Drücke** auftreten; diese sind auf der logarithmischen Achse nicht darstellbar.

Diese Erweiterung berechnet den Dampfspinodalzweig selbst **aus der in CoolProp HEOS hinterlegten Zustandsgleichung**. Die Flüssigkeitsspinodale, Gleichgewichtsflächen und weiteren Diagrammtypen bleiben im bisherigen Rechenweg. Es werden keine Punkte durch Kurvenextrapolation oder graphisches Ergänzen erzeugt.

## 2. Bisherige Methode und Datenvertrag

`calc_spinodale(stoff)` liefert unverändert zwei sortierte `(n,3)`-Arrays für flüssigkeits- und dampfseitige Stabilitätsgrenzen. Aus den HEOS-Werten werden `T=T_reducing/tau`, `rho_molar=delta*rhomolar_reducing` und `v_molar=1/rho_molar` bestimmt. Für den homogenen Druck erfolgt ein `DmolarT_INPUTS`-Flash mit vorgegebener Phase. Die Schnittstelle bleibt **`(T [K], v_molar [m³/mol], p [Pa])`**.

AP3-`prepare_spinodal_branch()` und `Debug.test_spinodale_ableitung()` werden weiterhin nach Auswahl des Datenquellenzweigs aufgerufen. AP3-Anzeigekriterien bleiben erhalten; insbesondere beträgt die bestehende Prüftoleranz **`1e-2`** für normierte Ableitungsresiduen.

## 3. Bisherige Wasser-Messergebnisse (NICHT hier neu gemessen)

Quelle: vorherige nutzerseitig dokumentierte Tests vom 19.10.2026.

| T [K] | ρ/ρc | v/vc | p/pc | p/p_sat |
|---:|---:|---:|---:|---:|
| 580 | 0,274950 | 3,637028 | 0,515090 | 1,202898 |
| 600 | 0,355207 | 2,815258 | 0,628857 | 1,123961 |
| 620 | 0,464580 | 2,152482 | 0,762474 | 1,058026 |
| 630 | 0,541137 | 1,847960 | 0,839532 | 1,030857 |
| 640 | 0,658638 | 1,518284 | 0,926810 | 1,009076 |

Im gescannten Dichtebereich wurde jeweils eine **Druckmaximum-Nullstelle** (`dp/drho_m`: `+ -> -`) und eine weitere **Druckminimum-Nullstelle** (`- -> +`) beobachtet. Die zweite Nullstelle darf nicht als gesuchte Dampfspinodale verwendet werden, auch wenn sie mathematisch zur homogenen EOS-Fortsetzung gehört. Die interne Ursache des früheren CoolProp-Tracer-Endes wurde nicht nachgewiesen.

## 4. Thermodynamische Herleitung

Für die mechanische Stabilitätsgrenze einer homogenen Phase gilt:

\[
\left(\frac{\partial p}{\partial v_m}\right)_T=0,
\qquad v_m=\rho_m^{-1},
\qquad \left(\frac{\partial p}{\partial v_m}\right)_T
=-\rho_m^2\left(\frac{\partial p}{\partial\rho_m}\right)_T.
\]

Bei positiver Dichte haben daher beide Ableitungen dieselben Nullstellen. Der Solver verwendet ausschließlich

```python
state.specify_phase(CP.iphase_gas)
state.update(CP.DmolarT_INPUTS, rho, T)
state.first_partial_deriv(CP.iP, CP.iDmolar, CP.iT)
```

Eine numerische Nullstelle allein garantiert **keine** korrekte physikalische Astzuordnung.

## 5. Nullstellensuche und Mehrfachnullstellen

`spinodal_vapor_solver.py` erzeugt gleichmäßig verteilte unterkritische Temperaturen. Für jede Temperatur werden `rho_sat,v` und `p_sat` aus `QT_INPUTS(Q=1,T)` mit HEOS berechnet. Das Dichtesuchintervall lautet näherungsweise `(rho_sat,v, rho_c)` ohne die numerisch kritischen exakten Endpunkte. Darin werden logarithmisch verteilte Dichten untersucht.

Vorzeichenwechsel identifizieren Kandidaten; identisch getroffene Raster-Nullstellen werden mit den unmittelbaren Nachbarableitungen klassifiziert. Maxima (`+ -> -`) und Minima (`- -> +`) werden unterschieden. Die Verfeinerung der Maximum-Kandidaten erfolgt mit SciPy `brentq` (`rtol=1e-12`, `maxiter=100`). Eine Ableitung mit NaN unterbricht das betroffene Suchintervall, anstatt darüber hinweg zu interpolieren. Falls zunächst kein Maximum gefunden wird, wird das Dichtegitter kontrolliert auf das Zwei- bzw. Vierfache verdichtet. Ein Fehlschlag bleibt ein Fehlschlag.

**Freigabebedingungen eines Dampfpunkts:**

- Maximum nach beobachtetem Vorzeichenwechsel, nicht das Druckminimum;
- `rho_sat,v < rho_spin,v < rho_c`;
- `p_spin,v > p_sat` mit kleinem numerischem Vergleichsabstand;
- endlicher homogener Druck und internes dimensionsloses Ableitungsresiduum `R_s = |rho_c/p_c * dp/drho_m|` höchstens `1e-8` (konfigurierbar);
- genau **ein** physikalisch qualifizierbares Maximum in diesem Temperaturintervall.

Werden mehrere qualifizierte Maxima gefunden, wird **kein** Zweig geraten: Die Temperaturzeile wird als ungültig gekennzeichnet. Dadurch kann eine heiße EOS-Fortsetzung mit mehrfachen stationären Stellen keine scheinbar gültige Dampfspinodale erzeugen.

## 6. Kritischer Bereich, Kontinuität, Fehlstellen

Berechnet wird höchstens bis `Tc - max(1e-5 K, 1e-6*Tc)` und bis zur vom Benutzer verlangten Grenze. Das ersetzt keine kritische Skalierungsanalyse; bei konditionell problematischen Temperaturen können Punkte ausfallen. Ein unmotiviertes Anhängen an `Tc` erfolgt nicht. Der bisherige AP3-Kriterium für einen kritischen Endpunkt gilt unverändert.

Die Dampfspinodale wird als `(n_T,3)` ausgegeben. Erfolgreiche Zeilen enthalten `(T,v,p)`, fehlende Zeilen ausdrücklich `(T,NaN,NaN)`. Monotonieverletzungen bei aufeinanderfolgenden Stützstellen werden als möglicher Astwechsel gekennzeichnet und ebenfalls verworfen. `visualization._draw_curves()` zerlegt die Polyline an nichtendlichen Punkten; **keine Verbindung über Lücken**.

Der Bericht protokolliert unter anderem: genutztes T-Intervall, Gesamt-/Erfolgszahl, fehlende T-Werte, Fehlergründe, Anteil gültiger Temperaturen, erstes/letztes gültiges T, Zahl der Druckmaxima/-minima, Newton-unabhängige Brent-Aufrufe, HEOS-Ableitungsauswertungen, deren Exceptions, Median/Maximum der Residuen und echte Laufzeit.

## 7. Funktionsschnittstellen und geänderte Dateien

| Datei | Tatsächliche Änderung | Zweck |
|---|---|---|
| `spinodal_vapor_solver.py` **neu** | `calc_spinodale_vapor_temperature()`, `_brackets()`, `compare_spinodal_overlap()` | numerischer Solver, Klassifikation, nur vergleichende Interpolation |
| `coolprop_heos_pvt.py` | Import der neuen Solverfunktionen; zusätzliche `spinodale_grid`-Validierung | kompatibler öffentlicher Berechnungseinstieg |
| `main.py` | temperaturparametrisierter Dampfzweig mit Legacy-Fallback, neuer Bericht `spinodal_solver` | keine Änderung am Flüssigast, AP3 bleibt erhalten |
| `config.yaml` | zusätzliche `n_T`, `n_scan`, `residual_tolerance`; sinnvoller Wasser-Startbereich | Rechen- und Darstellungsbereich getrennt |
| `gui.py` | numerische Spinodalparameter, Tooltip/Status, Statusbewahrung nach Threadende | bestehende QThread-Berechnung weiterverwenden |
| `visualization.py` | Spinodalen **nicht** mehr an Untergrenze des Einphasen-Meshs clippen | valide Spinodalpunkte sichtbar halten |
| `test_spinodal_extension.py` **neu** | Fake-EOS, Integration/Legacy-Fallback, optional echte HEOS-Prüfung | automatisierte Tests |
| `validate_spinodal_heos.py` **neu** | Wasser-/R134a-Test und drei Rasterkonvergenzen, JSON-Bericht | Messprotokoll auf Zielrechner |
| `requirements.txt` **neu** | SciPy und vorhandene Bibliotheken aufführen | Installationshinweise |
| `README_SPINODALE_ERWEITERUNG.md` **neu** | diese Datei | Nachvollziehbarkeit |

Die bisherige AP1–AP4-Quelldokumentation, Siede- und Taulinienfunktion, übrigen vier Berechnungsmodule und ursprüngliche AP3-Diagnose werden beibehalten.

## 8. Neue Konfigurationsparameter

```yaml
grid_points_creation:
  spinodale_grid:
    T_start_spinodal: 580.0    # Berechnungs-T-Untergrenze [K]
    T_stop_spinodal: 647.0     # Berechnungs-T-Obergrenze [K]
    n_T: 250                   # Zahl der Temperaturzeilen
    n_scan: 200                # initiale Dichtestützstellen je Temperatur
    residual_tolerance: 1.0e-8

visu_points_creation:
  spinodale_visu:
    T_start_spinodal: 580.0    # reine Anzeige-/Filtergrenze [K]
    T_stop_spinodal: 647.0
```

Die Berechnung verwendet die **Schnittmenge** aus Berechnungs- und Darstellungsintervall und schränkt sie auf den zulässigen unterkritischen HEOS-Bereich ein. Der automatische Stoffwechsel aktualisiert beide Temperaturgrenzen über die schon vorhandene `generate_stoff_grid()`-Logik. Der neue Solver nutzt nur die zusätzliche Konfiguration und wird bei älteren Konfigurationen ohne `n_T` nicht aufgerufen (rückwärtskompatible Legacy-Funktion). `n_T` liegt zwischen 2 und 10000, `n_scan` zwischen 16 und 10000, deren Produkt höchstens bei 2 Mio. Prüfpaaren; 0 < Residualtoleranz <= 1e-2.

Die Grenze von 2 Mio. ist eine Lastschutzmaßnahme. Durch adaptive Wiederholung kann die tatsächliche Zahl der Auswertungen größer sein; sie ist daher **keine** maximale tatsächliche Evaluationszahl.

## 9. Auswahlstrategie und Rohdatenqualität

`main.py` berechnet zunächst die bisherige CoolProp-Kurve für beide Äste und danach – bei entsprechender Konfiguration – die neue Dampfkurve. Nur dann wird umgeschaltet, wenn:

1. mindestens 10 und mindestens so viele gültige neue Stützstellen wie vorhandene gültige Legacy-Punkte vorliegen;
2. das neue Temperaturgitter am unteren Rand beginnt und mindestens **90 %** gültige Zeilen besitzt;
3. keine mehrdeutigen physikalisch qualifizierten Maxima gemeldet wurden;
4. das obere Ende des alten, sichtbaren Dampfastes bis auf eine rasterbedingte Toleranz ebenfalls abgedeckt ist;
5. bei mindestens drei überlappenden Referenzpunkten die maximale Vergleichsabweichung höchstens **6 % für `v`** und **3 % für `p`** beträgt.

Diese letzten beiden Vergleichsschwellen sind **vorläufige konservative Heuristiken**, keine nachgewiesenen physikalischen Genauigkeitsgrenzen. Die zum Vergleich notwendige lineare Interpolation wird **niemals** zur Erzeugung von Produktivpunkten eingesetzt. Bestehen Zweifel, bleibt der Legacy-Ast erhalten und der Fallback wird dokumentiert. Nach der Auswahl werden die regulären AP3-Prüfungen **immer** auf den ausgewählten Ast angewendet. Auch ein ausgewählter Ast kann daher transparent als unvollständig gekennzeichnet sein.

Eine Ergebnis-Dictionary besitzt weiterhin `spinodale_liquid`, `spinodale_vapor`, `spinodal_validation`; zusätzlich enthält sie `spinodal_solver` mit `selected`, `temperature`, `overlap`, ggf. `fallback_reason`.

## 10. Visualisierung und GUI

PyVista bildet physikalische Rohdaten ab durch `x=log10(v/vc)`, `y=T/Tc`, `z=log10(p/pc)`. Die zuvor verwendete Beschneidung `min_z=single_phase_grid.bounds[4]` wurde für Spinodalen aufgehoben: das Mesh bestimmt **keine** thermodynamische Gültigkeitsgrenze. Drucke `p<=0` bleiben korrekt nicht logarithmisch darstellbar, werden nicht positiv umgerechnet. Nullstellen, fehlende Temperaturen und Clipping sind getrennte Fehlerursachen.

Im linken Bedienpanel unter **Berechnung: Dampfspinodale (HEOS)** lassen sich Temperaturgrenzen, Temperaturpunktzahl, Dichtepunktzahl und Residualtoleranz konfigurieren. Unter **Darstellung: Spinodale** erfolgt unverändert die Sichtbarkeitsauswahl. Die Berechnung startet ausschließlich per **Berechnung starten** im bestehenden `CalculationWorker`-`QThread`. Der Abschlussstatus und der Tooltip enthalten die ausgewählte Datenquelle und die berechnete Abdeckung. Ein eigener Fortschritt **innerhalb** jeder Brent-Iteration sowie ein Abbruchknopf wurden nicht hinzugefügt, da der vorhandene Worker keinen Cancel-Vertrag bietet.

## 11. Hier tatsächlich ausgeführte synthetische Tests

Ohne CoolProp können Fake-Zustände die Algorithmen testen, aber **keine realen Wasser-/R134a-Werte bestätigen**. Die folgenden Resultate sind **neue Messungen in der verfügbaren Testumgebung**:

- Analytisch konstruierte Isothermen mit je einem lokalen Maximum und Minimum: **30/30** Temperaturpunkte gültig, maximale Residuen unter `1e-9`; Maximum korrekt, Minimum nicht ausgewählt.
- Absichtlich ungültiges Temperaturfenster (418–424 K) in einer anderen synthetischen Rechnung: **27/30** gültig, drei fehlende Zeilen behalten `NaN`.
- Konstruierte EOS-Fortsetzung mit **drei** physikalisch qualifizierbaren Maxima: **0/8** frei gegebene Stützstellen; alle acht Zeilen als mehrdeutig zurückgewiesen.
- Tests mit exakt auf Rasterpunkten liegenden Nullstellen, NaN-Barrieren, ungültigen Parametern, kritischer T-Grenze und Vergleich ohne Dateneinmischung: erfolgreich.
- `main.py`-Dispatcher unter kontrollierten Stubs: qualifizierter neuer Ast wird ausgewählt; unzureichende Lösung führt zu dokumentiertem Legacy-Fallback.

Diese Tests zeigen die **Programmstruktur und Algorithmik**, nicht die physikalische Genauigkeit von HEOS.

## 12. Wasser-Referenzen: offene numerische Abnahme

`test_spinodal_extension.VaporSolverHEOSTests.test_water_reference_points` und `validate_spinodal_heos.py` enthalten die Referenzliste `T=[580,600,620,630,640] K`. Im echten HEOS-Test werden `rho/rho_c`, `v/v_c`, `p/p_c`, `p/p_sat` und interne Residuen verglichen (`rtol=1e-4`, `atol=2e-5` für die gerundeten reduzierten Größen; internes `R_s<=1e-8`). Damit ist der Solltest **implementiert**, aber **hier nicht bestanden/nachgewiesen**, da CoolProp nicht installiert ist. Zuverlässige Aussagen über das reproduzierte Intervall 580–647 K lassen sich aus synthetischen Tests nicht ableiten.

## 13. Legacy-Vergleich: offene numerische Abnahme

`compare_spinodal_overlap()` interpoliert ausschließlich den neuen Referenzzweig auf Temperaturpunkte der bisherigen `build_spinodal()`-Kurve (keine eingefügten Produktivpunkte; kein Vergleich über ausdrücklich fehlgeschlagene Temperaturzeilen). Sie ermittelt maximale relative Druck-/Volumenabweichungen. Der Zielrechner-Test vergleicht insbesondere den vorher dokumentierten Bereich ab ungefähr 625 K. Aktuelle echte Vergleichszahlen liegen **nicht** vor. Das erste gemessene CoolProp-Tracer-Ende darf nicht als unveränderliches Fluidgesetz behandelt werden.

## 14. Zweites Fluid und Konvergenz-/Laufzeittests

`validate_spinodal_heos.py` berechnet `Water` (580–647 K) und `R134a` (0,9 Tc bis 0,997 Tc) jeweils mit **250×200**, **500×400** und **500×800** (Temperaturpunkte × initiale Dichtestützstellen). Gemessen werden zusätzlich HEOS-Ableitungsanzahl, Brent-Aufrufe, Residuen, gültige Temperaturabdeckung und reale Dauer. Die zu erstellende Datei `spinodal_heos_validation.json` enthält Versionsangaben, numerische Resultate und AP3-Bewertung. **Keiner dieser echten HEOS-Läufe wurde hier durchgeführt**; eine Rechenzeit für echte Fluidmodelle wird nicht behauptet.

Der Regelwert 250×200 entspricht bereits 50.000 initialen Ableitungsauswertungen zuzüglich Root- und Sättigungsaufrufen, bei adaptiver Verdichtung deutlich mehr. Eine Optimierung durch Kurvenfortsetzung ist noch nicht implementiert, da sie eine zusätzliche Astwechselprüfung erfordern würde.

## 15. Regressionstests (lokal ohne CoolProp)

Testbefehle im entpackten Projektordner:

```powershell
python -m pip install -r requirements.txt
python -m unittest Test.RegressionTests Test.AP3SyntheticTests test_ap4.AP4SyntheticTests test_spinodal_extension.VaporSolverSyntheticTests test_spinodal_extension.VaporSolverIntegrationTests -v
python -m unittest test_spinodal_extension.VaporSolverHEOSTests Test.AP3HEOSTests test_ap4.AP4HEOSTests -v
python -m unittest test_ap4.AP4RenderTests -v
set THERMOSURFACE_AP4_GUI=1
python -m unittest test_ap4.AP4GUISmokeTests -v
python validate_spinodal_heos.py
python main.py
```

**Tatsächliche Umgebung:** NumPy und SciPy vorhanden, CoolProp/PyVista/PyVistaQt/PySide6 nicht vorhanden. Ein Installationsversuch für CoolProp scheiterte an fehlendem Netz-Zugriff. Der abschließende kombinierte Testlauf umfasste **42 Tests: 31 bestanden, 11 übersprungen, 0 unerwartet fehlgeschlagen** (siehe `TESTPROTOKOLL.txt`). Echte HEOS-/Qt-/VTK-Tests stehen aus und werden durch `unittest` als `SKIP` gekennzeichnet.

`Test.RecordedRegressionTests` schlägt **erwartungsgemäß** fehl, weil die ursprünglichen Snapshot-Dateien im bereitgestellten Projektstand nicht enthalten sind (`regression/baseline`). Diese Baselines dürfen nicht durch neue Ergebnisse ersetzt werden. Änderungen an `spinodale_vapor` verändern berechtigt bisherige p-v-T-Snapshots, weshalb Vergleiche separat zu dokumentieren sind.

## 16. Bibliotheksstände

Die tatsächlich geprüfte Container-Umgebung verwendet **Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0 und PyYAML 6.0.3**. **CoolProp, PyVista, PyVistaQt und PySide6 sind in diesem Container nicht installiert**. Versionen auf dem Windows-Zielrechner wurden nicht automatisch ermittelt und werden nicht erfunden. `validate_spinodal_heos.py` erfasst diese bei Ausführung in einem JSON-Bericht.

## 17. Bekannte Grenzen und endgültiger Stand

**Implementiert und synthetisch getestet:** neuer HEOS-basierter Berechnungspfad, Max-/Min-Klassifikation, adaptive Dichtescans, Residualkontrolle, NaN-Gap-Markierung, Legacy-Fallback, GUI-Konfiguration/-Status, freie Spinodal-Zeichnung ohne Mesh-Clipping sowie neue Testinfrastruktur.

**Nicht verifiziert:** numerische Wasser-Referenztreffer, Vergleich zur echten bisherigen CoolProp-Kurve, vollständige 580–647-K-Abdeckung, zweite reale Flüssigkeit, Leistung und Rasterkonvergenz mit HEOS sowie GUI/VTK auf Windows. Deshalb ist die wissenschaftlich-numerische Endabnahme des Arbeitspakets **offen** und es darf **nicht** als vollständig physikalisch validiert gelten.

**Fachlich offen:** Die einfachen Ast- und Monotoniekriterien sind nicht für beliebige EOS-Fluide beweisbar; nicht erkannte tangentiale Doppelpunkte ohne Vorzeichenwechsel sind möglich; extrem kritische Zustände bleiben numerisch schwierig. Der Vergleich gegen einen ggf. selbst unvollständigen Legacy-Tracer ist eine Vorsichtsmaßnahme, keine unabhängige thermodynamische Referenz.

**Übergabe:** Das ZIP enthält den kompletten vorliegenden Quellcode samt neuem Solver, Konfiguration, neuen Tests, Dokumentation und Abhängigkeitsliste. Den echten Zielrechnerlauf zunächst mit `test_spinodal_extension.py`, anschließend `validate_spinodal_heos.py`, zuletzt `main.py` und der vollständigen AP4-GUI-Testmatrix ausführen. Bei Abweichungen Fallback und Fehlerberichte untersuchen, statt Druckwerte oder Stützstellen optisch zu korrigieren.
