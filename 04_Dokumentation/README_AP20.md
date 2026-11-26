# AP20 – Vollständige und zusammenhängende homogene Zustandsoberfläche zwischen den Binodalen

**Stand:** 22.10.2026 · **Voraussetzung:** AP19 abgeschlossen und validiert (`README_AP19.md`)
**Leitgedanke (laut Aufgabenstellung):** AP20 deckt keine AP19‑Fehler geometrisch zu. Die zusammengeführte Fläche besteht ausschließlich aus Knoten, die das AP19‑Kriterium erfüllen; wo die EOS keine gültigen Zustände liefert, bleibt eine **dokumentierte Lücke**.

---

## 1. Umsetzung

Reihenfolge der Teilflächen je Temperatur:
`Siedelinie → Flüssigkeitsspinodale → Dampfspinodale → Taulinie` = **metastabile Flüssigkeit | instabile Fortsetzung | metastabiler Dampf**.

| AP | Umsetzung | Ort |
|---|---|---|
| 20.1 gemeinsames T‑Raster, vier Grenzpunkte | alle drei Flächen bekommen dieselben `sat_T`‑Zeilen; `build_homogeneous_surfaces()` verlangt exakt gleiche Temperaturen der vier Randkurven und prüft je Zeile `v_sat,l<v_spin,l<v_spin,g<v_sat,g`. Fehlende Spinodalpunkte bleiben NaN‑Zeile (`missing_rows`, `reasons`), nie interpoliert | `homogeneous_pvt.py` |
| 20.2 drei randkonforme Volumenraster | `[v_sat,l, v_spin,l]`, `[v_spin,l, v_spin,g]`, `[v_spin,g, v_sat,g]`; innere Drücke immer aus der homogenen HEOS (explizite Phasenvorgabe, kein Flash) | `homogeneous_pvt.py` |
| 20.3 exakte gemeinsame Randknoten | Randzeilen werden aus **denselben** Arrays kopiert; Test auf **bitgleiche** `(T,v,p)` | `test_ap19_ap20.py`, `validate_homogeneous_boundaries()` |
| 20.4 fehlende Punkte minimieren | siehe Abschnitt 3 | |
| 20.5 zusammenhängendes Netz | VTK‑freie Triangulation + Topologieprüfung; Renderer und Prüfung nutzen denselben Code | `grid_triangles`, `mesh_topology_report`, `combined_surface_report`, `visualization.homogeneous_*` |
| 20.6 kritischer Punkt | zusätzliche T‑Zeilen, Diagnosebericht, **optionale** validierte Schließung | Abschnitt 4 |
| 20.7 Abdeckung quantitativ | `surface_coverage`, `combined_surface_report`, Debug‑Log, GUI‑Tooltip | Abschnitt 2 |
| 20.8 Validierung | Wasser, CO₂, Argon, Methan, R134a × 3 Auflösungen | Abschnitt 5 |

**Zeilenlayout der instabilen Fläche** (nötig, damit alle Flächen rechteckige `(n_T, n_v, 3)`‑Arrays und die GUI‑Schnittstelle behalten): `n_v = n_L + 1 + n_R` mit `n_L=⌊(n_v−1)/2⌋`. *Reguläre* Zeilen sind ein durchgehendes Raster (Trennspalte = gewöhnlicher Knoten). *Artefakt*‑Zeilen (AP19) haben links `[v_l, v_a]`, rechts `[v_b, v_g]` und eine NaN‑Trennspalte; benachbarte Zeilen gehen am Übergang kontinuierlich ineinander über (`v_a→v_b`, wenn die Artefakte verschwinden). Der Mesher verbindet nur Zellen mit ≥3 gültigen Ecken, die Trennspalte wird also nie überbrückt.

**GUI:** Die Actor‑Schlüssel (`metastable_liquid_surface`, `unstable_surface`, `metastable_vapor_surface`) und alle Checkboxen sind unverändert; jede Fläche bleibt ein eigenes Mesh und einzeln ein‑/ausblendbar. Die separate Gleichgewichts‑Zweiphasenfläche (`zweiphasengebiet`) wurde nicht angefasst. Die Tooltips der drei Flächen nennen jetzt zusätzlich Zellabdeckung und Anzahl der Artefaktzeilen.

## 2. Mesh‑ und Abdeckungskennzahlen (AP20.5 / AP20.7)

`main.berechne_daten()` legt sie unter `daten["homogeneous_validation"]` ab (`surfaces[<Fläche>]["coverage"]`, `["combined"]`, `["critical_approach"]`); `Debug.diagnose_homogeneous()` protokolliert sie.

* **Knoten:** gültig / NaN / nicht renderbar (p≤0 auf Log‑Achse) / strukturelle Lücke.
* **Zellen:** *auswertbar* = Zelle zwischen zwei T‑Zeilen mit vorhandenen Randpunkten **ohne** Beteiligung der EOS‑Artefakt‑Trennspalte; *abgedeckt* = vier gültige Ecken; *teilabgedeckt* (3 Ecken); *ausgeschlossen* = Artefaktzellen (eigene Zahl, nicht stillschweigend aus der Statistik entfernt).
  `coverage = abgedeckt / auswertbar` (thermodynamisch) und `coverage_rendered` (mit gewählter p‑Achse).
* **Dreiecke**, **Komponenten** (Union‑Find über Dreieckskanten), **offene Kanten** (*Rand* = Gitterrand/Randkurve, *Lücke* = grenzt an dokumentierte Lücke/Teilzelle, **unerklärt** = grenzt an vollständig gültige Nachbarzelle → Nahtfehler), **Orientierungskonflikte** (gerichtete Kante doppelt), **Projektionsfaltungen** (Vorzeichenwechsel der Dreiecksfläche in der `(log v, T)`‑Ebene ⇒ Fläche wäre nicht mehr einwertig ⇒ Selbstüberschneidung), nicht‑mannigfaltige Kanten.
* **Zusammengeführte Fläche:** Knoten mit bitgleichen `(T,v,p)` werden verschmolzen; geprüft werden Komponenten, offene Kanten zwischen gemeinsamen Knoten und die beiden Nähte (erwartete Kanten, fehlende Kanten links/rechts, Knotenidentität).

## 3. Fehlende Zustandspunkte (AP20.4)

Analyse der bisherigen Behandlung und Ergebnis:

| Ursache für fehlende Punkte | Frühere Behandlung | Jetzt |
|---|---|---|
| Spinodalzeile verworfen wegen Mehrfachextrema (EOS‑Artefakt) | gesamte Zeile verworfen (Wasser 2/59, CO₂ 4/29, R134a 20/26 Zeilen) | sättigungsnächste Nullstelle (physikalisch definiert: erster Vorzeichenwechsel von dp/dρ ab Sättigung); Mehrfachextrema in `ambiguous_temperatures` protokolliert. Wasser, CO₂, R134a: alle zuvor verworfenen Zeilen wieder da (Methan war schon vollständig) |
| Spinodalzeile verworfen: „Unplausibler Astverlauf“ (v(T) nicht monoton) | NaN | unverändert (Argon: 2–3 von 29–59 Zeilen). Ursache ist EOS‑Welligkeit (Argon ≈131 K: erste Nullstelle springt auf ρ/ρ_c=1,24, p≈−57 MPa). Eine Zeile mit solch einem Sprung zu behalten oder die Nachbarzeilen zu verwerfen wäre willkürlich; sie bleiben als Lücke ausgewiesen |
| EOS‑Artefaktzone der instabilen Fläche | Ausreißerfläche | dokumentierte Trennspalte (AP19), Zellen separat gezählt |
| einzelne innere Knoten (HEOS‑Exception, Klassenabweichung, Druckfenster, Monotonie) | NaN‑Knoten | unverändert lokal, `point_reasons`; **in allen 15 Validierungsläufen keine einzige solche Lücke** (`partial_rows=0`) |
| Druck ≤ 0 bei `pressure_axis: log` | nicht renderbar | unverändert (Darstellungsgrenze, keine Datenlücke): `coverage_rendered` weist sie getrennt aus; `asinh` zeigt sie |

**Adaptive Stützstellen**, soweit thermodynamisch zulässig:
1. *Randgrenzen‑adaptiv:* in Artefaktzeilen werden die Äste per Bisektion exakt bis zur Gültigkeitsgrenze (`v_a`, `v_b`) bestimmt und mit voller Knotenzahl aufgelöst – bei Wasser/600 K ist der gültige Ast nur ≈0,1 % breit und mit einem festen Raster nicht erfassbar.
2. *Kritische Annäherung:* zusätzliche T‑Zeilen `T_c − (T_c − T_letzte)/2^k, k=1..4` (`critical_refinement_rows`, nur für die homogenen Flächen; Binodal‑/Zweiphasen‑Arrays bleiben unverändert).
3. Nicht umgesetzt: krümmungsgesteuerte Knoteneinfügung je Zelle und Ersatzknoten für gescheiterte innere Punkte. Letztere würden die in `test_corrections.py` festgelegte Regel („Punktfehler lokal behalten, nichts ersetzen“) verletzen, und es gab im Validierungsumfang keine solchen Fehler; erstere hätte die bitgleiche Rückwärtskompatibilität regulärer Zeilen (`geomspace`) aufgegeben. Fehlende Drücke werden **nirgends** zwischen den Spinodaldrücken interpoliert.

## 4. Kritischer Punkt (AP20.6)

Konfiguration (`config.yaml`, `diagram_settings.pvt`): `critical_refinement_rows: 4`, `close_critical_point: false`, `unstable_scan_points: 256`.

`critical_approach_report()` prüft in den letzten ≤6 Zeilen: T monoton; alle vier Breiten (Binodal, Spinodal, metastabile Flüssigkeit/Dampf) **streng monoton schrumpfend**; `v_c` zwischen allen Grenzkurven; `p_spin,l < p_sat < p_spin,g`; letzte Zeile innerhalb 10 % von p_c. Bei ΔT/T_c = 4·10⁻⁴ (Validierungsraster) ergibt sich `closure_supported=True` für Water, CO₂, Argon, Methan, R134a; die Breiten der letzten Zeile (in v_c) sind aber noch groß:

| Stoff | Binodal | Spinodal | metast. Flüssigkeit | metast. Dampf |
|---|---|---|---|---|
| Water | 0,311 | 0,211 | 0,041 | 0,060 |
| CO₂ | 0,264 | 0,152 | 0,050 | 0,062 |
| Argon | 0,224 | 0,141 | 0,034 | 0,048 |
| Methan | 0,233 | 0,154 | 0,032 | 0,046 |

Auf dem Standardraster der `config.yaml` ist die Schließung für **CO₂ nicht** unterstützt (Breite der metastabilen Flüssigkeit nicht streng monoton); für Wasser, Argon, R134a ja.

Entscheidung: **Schließung standardmäßig aus.** Ein Dreiecksfächer von der letzten Zeile zum kritischen Punkt (`homogeneous_surface_mesh(..., critical_cap=True)`) ist eine lineare Ausfüllung über ≈0,2–0,3 v_c, keine EOS‑Daten – das widerspräche „keine künstlichen Druckflächen“. Wer ihn dennoch möchte (z. B. für Abbildungen), setzt `close_critical_point: true`; er wird nur gezeichnet, wenn `closure_supported` **und** die letzte Zeile jeder angezeigten Fläche vollständig ist (`daten["critical_closure"]`). Unabhängig davon verfeinern die Zusatzzeilen die Annäherung mit echten HEOS‑Werten.

## 5. Validierung (AP20.8)

`python validate_ap19_ap20.py --fluids Water CO2 Argon Methane R134a --out results/ap19_ap20_validation.json` (je Fall ≈1–3 s). Volle Tabelle mit allen Kennzahlen: `results/ap19_ap20_validation.json`.

| Stoff | Aufl. | Dreiecke L / U / D | Zellabdeckung (thermo.) L / U / D | ausgeschl. Artefaktzellen (U) | Komponenten je Fläche L/U/D | Komponenten kombiniert (erwartet) | offene Kanten unerklärt | Nähte: fehlende Kanten, identische Knoten |
|---|---|---|---|---|---|---|---|---|
| Water | niedrig | 201 / 184 / 216 | 100 % / 100 % / 100 % | 14 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |
| Water | mittel | 962 / 942 / 1064 | 100 % / 100 % / 100 % | 32 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |
| Water | hoch | 4058 / 4094 / 4524 | 100 % / 100 % / 100 % | 68 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |
| CO2 | niedrig | 203 / 180 / 216 | 100 % / 100 % / 100 % | 18 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |
| CO2 | mittel | 987 / 944 / 1064 | 100 % / 100 % / 100 % | 40 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |
| CO2 | hoch | 4157 / 4162 / 4524 | 100 % / 100 % / 100 % | 78 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |
| Argon | niedrig | 180 / 170 / 216 | 100 % / 100 % / 100 % | 6 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |
| Argon | mittel | 877 / 870 / 1064 | 100 % / 100 % / 100 % | 14 | 2 / 3 / 1 | 2 (3) | 0 | 0/True, 0/True |
| Argon | hoch | 3912 / 3914 / 4524 | 100 % / 100 % / 100 % | 36 | 2 / 3 / 1 | 2 (3) | 0 | 0/True, 0/True |
| Methane | niedrig | 203 / 190 / 216 | 100 % / 100 % / 100 % | 12 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |
| Methane | mittel | 987 / 975 / 1064 | 100 % / 100 % / 100 % | 28 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |
| Methane | hoch | 4180 / 4230 / 4524 | 100 % / 100 % / 100 % | 56 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |
| R134a | niedrig | 201 / 202 / 216 | 100 % / 100 % / 100 % | 4 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |
| R134a | mittel | 950 / 993 / 1064 | 100 % / 100 % / 100 % | 8 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |
| R134a | hoch | 4022 / 4220 / 4524 | 100 % / 100 % / 100 % | 18 | 1 / 1 / 1 | 1 (1) | 0 | 0/True, 0/True |

Lesehilfe und ehrliche Einordnung:

* **100 % ist relativ zum „auswertbaren“ Bereich** (Abschnitt 2). Der Nenner schließt Artefaktzellen und T‑Zeilen ohne vier Randpunkte aus; beide werden getrennt ausgewiesen. Anteil der Artefaktzellen an allen Zellen der instabilen Fläche, z. B. Wasser/mittel: 32 von 560 (5,7 %); Wasser/hoch: 68 von 2 262 (3,0 %). Das 98 %‑Ziel bezieht sich auf das „zuverlässig berechenbare Untersuchungsgebiet“; die Artefaktzone ist nach AP19 nicht zuverlässig berechenbar und wird deshalb nicht mitgezählt.
* **Nicht auswertbare Zeilen:** Argon 1/13, 2/29, 3/59 (Astverlauf, s. Abschnitt 3); alle anderen Stoffe vollständig.
* **Komponenten kombiniert (erwartet):** „erwartet“ zählt je durchgängig berechenbarem T‑Intervall 1 Komponente (2, falls dort nur Artefaktzeilen liegen). Bei Argon sind es 2 statt 3, weil die metastabile Dampffläche über die Zeilen ohne Flüssigkeitsspinodale hinweg durchläuft und die Intervalle so verbindet – kein Netzfehler. In einem durchgehend berechenbaren Intervall mit mindestens einer regulären Zeile ist die zusammengeführte Fläche **eine** Komponente (Wasser, CO₂, Methan, R134a: ✔).
* **Rein artifaktbehaftete Intervalle** (niedrige T): linker und rechter Ast sind zwei Komponenten (keine Überbrückung – eine solche wäre eine thermodynamisch unbegründete Druckfläche).
* **Normalen/Orientierung:** 0 Orientierungskonflikte, 0 Projektionsfaltungen, 0 nicht‑mannigfaltige Kanten, 0 unerklärte offene Kanten in allen Läufen.
* **Selbstüberschneidungen:** nur über die Projektionsfaltung geprüft (Fläche bleibt einwertiger Graph über `(v,T)`); eine allgemeine 3D‑Schnittprüfung von Dreiecken wurde nicht implementiert.

Fünf Diagrammtypen: `berechne_auswahl()` für `pvt, pht, pst, hsp, uvt` (Wasser, Standardkonfiguration) läuft fehlerfrei; die Module `coolprop_heos_pht/pst/hsp/uvt` und deren Visualisierung wurden nicht verändert.

## 6. Dateien, Tests, Reproduktion

Geändert/neu: siehe `README_AP19.md` Abschnitt 6 (`config.yaml` bekommt `unstable_scan_points`, `critical_refinement_rows`, `close_critical_point`). Neu in AP20 zusätzlich im Code: `surface_coverage`, `combined_surface_report`, `critical_approach_report`, `critical_refinement_rows`, `mesh_topology_report`, `homogeneous_mesh_diagnostics`, `critical_cap`.

```powershell
python -m unittest test_ap19_ap20 test_ap13_ap18 test_corrections test_spinodal_extension test_integration_mocked_main -v   # 51 Tests, 1 vorbestehend übersprungen
python validate_ap19_ap20.py --fluids Water CO2 Argon Methane R134a --out results/ap19_ap20_validation.json
python diagnose_ap19.py
python validate_ap17_ap18.py                                   # HEOS_CORE_PASS
```

Kontinuitäts‑/Mesh‑Tests (`test_ap19_ap20.py`): bitgleiche Randknoten (Flüssigkeit|instabil|Dampf), Naht‑Kantenabgleich, eine Komponente ohne offene Nähte bei regulärer EOS, keine Überbrückung der Artefaktlücke, Abdeckung/Teilzellen/fehlende Zeile, Erkennung von Nahtfehler/Orientierungsfehler/Faltung, Dreiecksfächer zum kritischen Punkt, Zusatzzeilen nahe T_c, gesperrte Zeile bei verletzter Grenzordnung.

## 7. Grenzen

* Die Qt‑Oberfläche wurde importiert, aber nicht interaktiv bedient (Checkbox‑Klicktests stehen aus); Berechnung, Debug‑Pfad und PyVista‑Offscreen‑Darstellung wurden ausgeführt.
* Die kritische Zusatzverfeinerung ist ab Werk aktiv (4 Zeilen). Für reine Rückwärtskompatibilität der Arrayform `critical_refinement_rows: 0` setzen. Bereits gespeicherte Regressions‑Snapshots (`regression/…`) enthalten die alte Zeilenzahl; ein neuer Snapshot braucht eine bewusste Freigabe (`snapshot_pvt.py`), er wurde nicht überschrieben.
* Die Spinodal‑Wurzelwahl (sättigungsnächste Nullstelle) ist eine Konvention; sie ist mit der Definition der metastabilen Zone (zusammenhängend ab Sättigung) konsistent, ändert aber `solve_spinodal_branch` gegenüber AP13–AP18 (vorher: Zeile verworfen).
* Ein ZIP für AP19 und AP20 gemeinsam: die Änderungen an `homogeneous_pvt.py`/`visualization.py` sind ineinander verzahnt und nicht sinnvoll getrennt lieferbar (`ThermoSurface_AP19_AP20.zip`).
