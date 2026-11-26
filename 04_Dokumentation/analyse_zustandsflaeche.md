# Analyse und Umsetzung: zusammenhängende homogene Zustandsfläche p(T,v) zwischen den Spinodalen

**Stand:** 22.10.2026 · **Branch:** `feature/homogene-zustandsflaeche` (Basis `master` = Projektstand nach AP20) · **Stoff der config.yaml:** Water · CoolProp 7.2.0, PyVista 0.47.3, NumPy 2.3.4, Python 3.13
Kennzeichnung: **[Code: Datei, Funktion]** · **[Berechnet: …]** · **[Annahme]** · **[Unbestätigt]**. Reproduktion der Analyse: `python scripts/diagnose_zustandsflaeche.py` (→ `results/diagnose_zustandsflaeche.json`).

---

## 1. Verifizierte Ursache

Zuerst der Datenweg: EOS → Spinodalen → Ableitungen/Filter → Raster → Triangulierung → Achse → Mesh/GUI.
[Code: `homogeneous_pvt.homogeneous_point` (HEOS, explizite Phase) → `solve_spinodal_branch` → `unstable_node_status`/`trace_unstable_row` → `build_homogeneous_surfaces` (Raster) → `grid_triangles` → `visualization.to_xyz` → `render_pvt`].

| Nr. | Hypothese | Ergebnis | Beleg | Ursachenklasse |
|---|---|---|---|---|
| H1 | zusätzliche Extrema zwischen den Spinodalen | **bestätigt** | [Berechnet: Water, 4000 Dichten] Vorzeichenwechsel von (∂p/∂ρ)_T: **4** bei T/T_c = 0,6; 0,8; 0,95; 0,99 (2 Spinodalen + 1 Extremapaar); Argon **6** bei 0,8, 4 bei 0,6/0,95, 2 bei 0,99; CO₂ 4. Überkritisch (1,0; 1,05): 0. p/p_c, Water 0,8: −8,8·10⁴ … +1,6·10⁵; 0,95: −1,7 … +5,0 | Modellverhalten (EOS) |
| H2 | Filter entfernen Punkte | **bestätigt** | [Berechnet: Water, 20×40 / 60×120 / 150×300] gezeigt 16,0 / 13,3 / 14,6 %; entfernt durch dp/dv<0: 13,0 / 12,4 / 12,0 %; Druckfenster: 55,7 / 57,0 / 57,5 %; Nicht-Zusammenhang: 15,2 / 17,3 / 15,9 %. Argon: gezeigt 35,0 / 33,1 / 26,0 %. [Code: `homogeneous_pvt.unstable_node_status`, `zustandsflaeche.apply_display_filter`] | Filter (Folge von H1) |
| H3 | verworfene Punkte werden NaN, Rohwert geht verloren | **bestätigt (bis AP20), jetzt behoben** | [Code: `homogeneous_pvt.build_homogeneous_surfaces`: `row[j,2]=p` nur für bestandene Knoten, sonst bleibt NaN; Rohwert wird nicht gespeichert]. [Berechnet: Water] allein der Vorzeichenfilter verwarf 104 / 896 / 5383 endliche Rohwerte. Neu: `HomogeneousField` hält `p_raw` unverändert, Filter setzen nur Darstellungsstatus | Filter + Datenmodell |
| H4 | Vernetzung überbrückt Lücken nicht | **bestätigt; Vernetzung ist nicht die Ursache** | [Code: `homogeneous_pvt.grid_triangles`, `zustandsflaeche.grid_triangles_strict`] Dreieck nur bei gültigen Ecken; Test `test_ap19_ap20.test_all_rows_have_gap_and_are_not_bridged`. Die Zacken entstehen, weil die Anzahl gezeigter Knoten je T-Zeile springt (H6), nicht durch fehlerhafte Dreiecke | Filter / Modellverhalten |
| H5 | Log-Achse blendet p ≤ 0 aus | **bestätigt** | [Berechnet: Water] Anteil p ≤ 0 der Rohknoten 55–56 % (Argon 36–37 %); je Isotherme 64 % (0,6), 54 % (0,8), 18 % (0,95), 0 % (0,99). [Code: `visualization.to_xyz` (vorher `np.log10(p/pc)`), `homogeneous_pvt.node_validity`: log verlangt p>0] | Achsendarstellung |
| H6 | Raster beeinflusst nur Zacken, nicht Lücken | **bestätigt** | Anteile gezeigt/entfernt unabhängig vom Raster (Water 13–16 %); „Zackenmaß“ (mittlere Änderung gezeigter Knoten zwischen T-Zeilen) 6,7 / 18,6 / 39,9, also konstant ≈ 13–17 % der Spaltenzahl: feineres Raster macht die Zacken nur feiner | Filter / Modellverhalten |

**Pflichtprüfungen**

| Prüfung | Befund |
|---|---|
| Phasenvorgabe | [Berechnet: Water, 600 K, ρ = 21 372,7 mol/m³ (Zweiphasengebiet)] mit `specify_phase`: **p = 1,1639·10⁹ Pa** (rohe EOS); ohne Vorgabe: **p = 1,2345·10⁷ Pa** = p_sat (Flash; `PropsSI` identisch). [Code: `homogeneous_point` ruft `state.specify_phase(...)` vor `update(DmolarT_INPUTS,…)`] ⇒ die Auswertung liefert die rohe EOS, keinen Flash. Gegenprobe der Rohwerte mit unabhängiger IAPWS‑95‑Implementierung: identisch (siehe `README_AP19.md`) |
| Spinodalzuordnung | [Code: `solve_spinodal_branch`: bei mehreren Brackets wird die sättigungsnächste gewählt (`roots[-1]` Flüssigkeit, `roots[0]` Dampf)] = Regel „erste Nullstelle ab ρ_sat“. [Berechnet: Water, Argon, CO₂ bei T/T_c = 0,8; 0,95; 0,99] Übereinstimmung auf 10⁻⁶ (Test `test_zustandsflaeche.P2HeosSpinodalzuordnung`). Isothermen mit **mehr als zwei Nullstellen:** Water alle unterkritischen Testisothermen (4), Argon 0,6/0,8/0,95 (4/6/4), CO₂ (4); nur Argon 0,99 hat genau 2. Vor AP20 wurde bei >1 Treffer die Zeile verworfen („falsche Spinodalauswahl“ als Ursache fehlender Zeilen, behoben). CO₂ bei 0,6 T_c liegt unter dem Tripelpunkt: nicht berechenbar |
| Einheiten/Vorzeichen | [Code] molare Basis (`Dmolar`, v in m³/mol, p in Pa). [Berechnet: Water, 600 K, ρ = 12 000 mol/m³] (∂p/∂v)_T aus Differenzenquotient 8,026129987·10¹² Pa·mol/m³ = −ρ²(∂p/∂ρ)_T (rel. Abweichung 9·10⁻¹¹). Konvention im Code: `norm = slope·ρ_c/p_c`; `norm<0` ⇔ (∂p/∂v)_T>0 ⇔ instabil |
| Kenngrößen je Testisotherme | vollständige Tabelle (Extrema, p_min/p_max, Anteil p ≤ 0, Fehler/NaN) in `results/diagnose_zustandsflaeche.json` → `isotherms`; Fehler-/NaN‑Anteil der Rohauswertung: **0 %** für alle 18 Isothermen |

## 2. Unterschied zur Van-der-Waals-Darstellung und Aussagegrenzen

| Aspekt | vdW (generisch) | HEOS (diese Arbeit) | Ersatzfläche |
|---|---|---|---|
| p(v) zwischen Spinodalen | ein monotoner Anstieg p_L → p_V | 2–4 zusätzliche Extrema, p bis ±10⁴ p_c außerhalb [p_L, p_V] | monotoner Anstieg, kubisch (C¹) |
| Aussage | modellabhängige Fortsetzung | rohe EOS‑Extrapolation, in der instabilen Zone nicht an Messdaten angepasst | Modellfläche, **keine** EOS‑Auswertung |
| Randpunkte | Spinodalen des Modells | HEOS‑Spinodalen | **dieselben** HEOS‑Spinodalpunkte (identische Vertex‑Indizes) |

Grenzen: Die Ersatzfläche ist weder Gleichgewichtsfläche noch Maxwell‑Konstruktion; sie hat keine Fundamentalgleichung (keine s, h, u, g, c_p; Hover zeigt nur T, v, p und das Ersatz‑Label). Ihre Krümmung stimmt nicht mit der HEOS überein (s. 3). [Unbestätigt] Wie die Masterarbeit ihre vdW‑Fläche implementiert: die PDF liegt zwar in `literatur/`, wurde für diese Arbeit aber **nicht** ausgewertet.

## 3. Entscheidung

| Weg | Kontinuität | Anschluss an HEOS‑Spinodalen | Aussagekraft | Aufwand | Befund (Zahlen) |
|---|---|---|---|---|---|
| **A** ungefilterte HEOS + Klassifikation | **nein** | exakt | hoch (EOS), aber nur auf Ästen | gering | H2: nur 13–16 % (Water) / 26–35 % (Argon) der Knoten sind instabil, im Fenster und zusammenhängend; bei 600 K: p(ρ) hat 4 Nullstellen ⇒ geht nicht gleichzeitig mit (∂p/∂v)>0 überall und Kontinuität |
| **B** kubische EOS | ja | **nein** (eigene Spinodalen) | eigenes Modell | hoch | [Berechnet: CoolProp 7.2.0] Backend `PR`/`SRK` für **Argon verfügbar**, für **Water nicht** („type not set“). Nicht Teil des Auftrags, nicht umgesetzt |
| **C** gekennzeichnete Ersatzfläche | ja | exakt (Randwerte, ∂p/∂v=0) | Modell, gekennzeichnet | mittel | **Hauptansatz** |

**Datenbasis ist A:** `HomogeneousField.p_raw` bleibt unverändert, nichts wird gelöscht; Filter setzen nur Stabilitätsklasse, Darstellungsstatus und Filtercode.

**Konflikt „unveränderte HEOS‑Werte / überall (∂p/∂v)>0 / Kontinuität“:** nicht gleichzeitig erreichbar. [Berechnet: Water, T/T_c = 0,8] ∂p/∂ρ hat 4 Nullstellen bei ρ/ρ_c = 0,116; 0,933; 1,235; 2,144 (Dampf‑/Flüssigkeitsspinodale 0,116 und 2,144, Zusatzextremapaar 0,933/1,235); p/p_c läuft über die Isotherme von −8,8·10⁴ bis +1,6·10⁵. Zwischen den Zusatzextrema (ρ/ρ_c 0,933 … 1,235) ist (∂p/∂v)_T<0, die Fläche dort also nicht „instabil“; bei 0,95 liegen die Nullstellen bei 0,432; 0,931; 1,233; 1,691. Messgröße der Lücke: nur 13–16 % (Water) bzw. 26–35 % (Argon) der Rohknoten zwischen den Spinodalen bleiben als instabiler, zusammenhängender Ast im Druckfenster.

**[Annahme]-Entscheidungen**

1. Ersatzfläche kubisch (C¹) statt quintisch: [Berechnet: Water/Argon/CO₂ bei 0,8; 0,95; 0,99] die HEOS‑Randkrümmungen haben die richtigen Vorzeichen (∂²p/∂v² > 0 flüssig, < 0 dampf), aber andere Beträge (z. B. Water 0,95: HEOS 8,2·10¹⁷ vs. kubisch 3,7·10¹⁵ an der Flüssigkeitsspinodale); eine quintische Hermite‑Variante mit den HEOS‑Krümmungen ist in **9 von 9** Fällen nicht monoton. Die Krümmungsunstetigkeit ist damit dokumentiert.
2. ε nahe T_c: **EPS_WIDTH = 10⁻² v_c** und **EPS_PRESSURE = 10⁻⁸ p_c**. [Berechnet: Water, 1−T/T_c = 10⁻⁶: Breite 0,0151 v_c, Δp 2,9·10⁻⁸ p_c; 10⁻⁷: 0,0049 v_c, 9,3·10⁻¹⁰ p_c, also unter der Druckauflösung 10⁻⁹ p_c der Spinodalen.] Kollabierte Zeilen erhalten keine Punkte; der Fächer startet an der letzten gültigen Zeile.
3. Ersatzfläche wird (wenn T_c erreicht wird) mit einem Dreiecksfächer an **einem** Apex‑Vertex (kritischer Punkt) geschlossen; das Fächer‑Stück ist Teil der gekennzeichneten Modellfläche.
4. HEOS‑Ansicht der instabilen Zone (`unstable_surface`) nutzt weiter die AP19‑Knotenverteilung (Bisektion an den Astenden), das Rohfeld (`homogeneous_raw_field`) liegt auf dem festen Indexgitter. Beide stammen aus derselben EOS.
5. Druckachse: `log` wird automatisch auf **`symlog`** umgestellt, wenn homogene Flächen p ≤ 0 enthalten (Hinweis in Log, `daten["pressure_axis_notice"]` und GUI‑Status; abschaltbar mit `auto_pressure_axis: false`). `symlog`: z = sign(p)·log₁₀(1+|p|/(10⁻² p_c)); Ticks in echten Druckwerten.
6. Git ist nicht installiert und das Verzeichnis war kein Repository: [Annahme] `dulwich` (nur Werkzeug, nicht in requirements) legt ein Repository an; Basis-Commit = Stand nach AP20; `literatur/` per `.gitignore` ausgeschlossen. `pytest` wurde nur zum Ausführen der Suite installiert (die Tests sind `unittest`‑Klassen und laufen unter beiden).
7. Die Ersatzfläche wird auch ohne `include_unstable` berechnet (`include_replacement: true`), 40 Spalten je Zeile (`replacement_points`).

## 4. Umgesetzte Pakete

| Paket | Commit | Inhalt | Tests (Befehl → Ergebnis) |
|---|---|---|---|
| P0 Baseline | `303de60e` (Basis `b7002588`) | `scripts/p0_baseline.py`, `regression/p0_baseline.json` (Gleichgewichtsfläche 1361 Punkte / 2613 Zellen, Bounds, SHA‑256 aller Diagramm‑Arrays), `scripts/diagnose_zustandsflaeche.py` | `python -m unittest …` Baseline 51 Tests OK; `python scripts/p0_baseline.py --compare regression/p0_baseline.json` → „OK (bitgleich)“ |
| P1 Datenmodell | `510743f2` | `zustandsflaeche.HomogeneousField`, `build_raw_field`, `apply_display_filter` (a)–(d) | `test_zustandsflaeche.P1Datenmodell` (9) OK: p_raw bitgleich vor/nach Filter/Vernetzung |
| P2 Randpunkte | `63b6bc6c` | `assemble_shared_vertices` (Union‑Find, identische Vertex‑Indizes), `first_zero_from_saturation` | `P2…` (6) OK: Randabweichung 0 (≤ 1e‑12), Zuordnung Water/Argon |
| P3 Ersatzfläche | `94664428` | `substitute_isotherm`, `build_substitute_field`, Statuscodes | `P3Ersatzflaeche` (10) OK: p_V>p_L monoton, p_V≤p_L „nicht konstruierbar“, Randwerte/Randableitung = 0 |
| P4 Vernetzung | `d411d1cd` | `grid_triangles_strict` (feste Diagonale), `mesh_report` (Randschleifen, innere Löcher) | `P4Vernetzung` (6) OK, inkl. Abgleich mit `extract_feature_edges(boundary_edges=True)` |
| P5 Kritischer Punkt | `21293cfa` | ε‑Kollaps, `fan_applicable`, `critical_fan` (ein Apex‑Vertex) | `P5KritischerPunktHeos` (2) OK: 0 degenerierte Dreiecke bei T/T_c ∈ {0,99; 1,0; 1,05} |
| P6 Druckachse | `07d51350` | `pressure_to_z` (linear/log/asinh/symlog), `effective_pressure_axis`, Ticks in echten Druckwerten, Renderer auf gemeinsamer Vertextabelle | `P6Druckachse` (5) OK: p<0, p=0, p>0 bleiben im Netz |
| P7 Kennzeichnung | `0290dbc0` | GUI‑Ebene „Ersatzfläche (Modell, keine EOS‑Auswertung)“ (Voreinstellung an), HEOS‑Rohfortsetzung exklusiv, Legende, Hover (`pick_hover_text`), CSV‑Export mit `source` | `P7Kennzeichnung` (4) OK; GUI‑Rauchtest (reale Qt‑Plattform) s. 5 |
| P8 Validierung | `3f13cf84` | `scripts/validate_zustandsflaeche.py`, `P8Validierungsmatrix`, `P8EndeZuEnde`, Regressionstest | siehe unten |

**Gesamtsuite:** `python -m pytest -q -p no:cacheprovider test_zustandsflaeche.py test_ap19_ap20.py test_ap13_ap18.py test_corrections.py test_spinodal_extension.py test_integration_mocked_main.py` → **100 passed, 2 skipped** (vorbestehend: fehlende `Test.py`), 46 Subtests.

**Validierungsmatrix** (`python scripts/validate_zustandsflaeche.py` → `results/validate_zustandsflaeche.json`, 21 Fälle: Water/Argon/CO₂ × 20×40 / 60×120 / 150×300 × Achse linear/log, plus „echte“ Log‑Achse im p>0‑Bereich T/T_c ≥ 0,97): **0 Verletzungen** von: innere Löcher (Ersatz und kombinierte Fläche) = 0 · degenerierte Dreiecke = 0 · Randabweichung 0 ≤ 1e‑12 mit identischen Vertex‑Indizes · Randpolylinie monoton in T · Orientierungskonflikte, nicht‑mannigfaltige Kanten, Projektionsfaltungen = 0 · Herkunftssumme = Gesamtpunktzahl · HEOS‑Rohwerte bitgleich · Regression gegen P0 bitgleich (Vertex‑/Zellzahl und Bounds der Gleichgewichtsfläche, alle Diagramm‑Arrays). Bei Achse `log` wird in den Wide‑Range‑Fällen auf `symlog` umgestellt (p ≤ 0 vorhanden), im p>0‑Bereich bleibt `log`. Einzelisothermen: Water und Argon: 0,6/0,8/0,95/0,99 „ok“ (Breite/v_c 90,2/8,13/1,73/0,76 bzw. 8,69/2,86/1,16/0,53), 1,0 „kritische Temperatur“, 1,05 „überkritisch“; CO₂ bei 0,6 „Spinodale fehlt“ (unter Tripelpunkt).
Visuell: `figures/zustandsflaeche_water_heos_gefiltert.png` (A: HEOS gefiltert, zwei schmale Wände) vs. `figures/zustandsflaeche_water_ersatz.png` (C: zusammenhängende Ersatzfläche), jeweils Wasser, 580–647 K, `symlog`.

Offene oder blockierte Pakete: **keine**.

## 5. Offene Grenzen, ungetestete Teile, Unbestätigtes

* **Hover im laufenden Qt‑Fenster [ungetestet]:** Die Logik (`pick_hover_text`, `hover_label`) ist getestet, auch mit dem VTK‑Picker offscreen. Der Maus‑Observer in `gui.ThermoSurfaceWindow._install_hover` (`iren.add_observer("MouseMoveEvent", …)`) lief nicht mit echter Mausbewegung. Qt `offscreen` scheitert hier an OpenGL (Absturz), der Rauchtest lief deshalb auf der normalen Windows‑Plattform mit kurz sichtbarem Fenster.
* **GUI‑Rauchtest (real):** Standardzustand (Ersatz an, HEOS‑instabil aus), Exklusivität der beiden Ebenen und Actor‑Sichtbarkeit, Rendern der Wasserdaten, Export (24 320 Knoten: HEOS 22 876 + Ersatz 1 444 = Gesamt). **Nicht** per Klick getestet: der Export‑Button/Dateidialog (nur `export_homogeneous_csv` direkt) und das Legenden‑Layout (im Screenshot am rechten Rand abgeschnitten).
* **Auto‑Achse ändert das Standardverhalten:** Bei `pressure_axis: log` und p ≤ 0 wird `symlog` verwendet. Die config.yaml des Nutzers enthielt zuletzt `asinh`; sie wurde nicht verändert.
* **Fächer zum kritischen Punkt:** Bei der Standardkonfiguration (Water) entsteht kein Fächer: die letzte Zeile (1−T/T_c ≈ 1·10⁻⁴, Breite ≈ 0,12 v_c) liegt über ε. Er wird nur genutzt, wenn Zeilen unter ε kollabieren oder T = T_c im Raster liegt (Validierungsmatrix: ja).
* Rasterzeilen ohne Spinodale (z. B. unter dem Tripelpunkt, Argon‑Astverlauf) bleiben Lücken mit Statuscode; die Ersatzfläche füllt sie nie.
* Stabilität der Ersatzfläche in (v,T): sie ist je Isotherme monoton; zwischen Isothermen wird nur linear vernetzt. Eine allgemeine 3D‑Selbstüberschneidungsprüfung gibt es nicht (nur Projektionsfaltung in (log v, T)).
* [Unbestätigt] Implementierung der Masterarbeit (`literatur/`), nicht gelesen.
* Die HEOS‑instabile Ansicht bleibt (bewusst) lückenhaft; die Lücke ist EOS‑Eigenschaft (H1).

## 6. Geänderte und neue Dateien

**NEU:** `zustandsflaeche.py`, `test_zustandsflaeche.py`, `scripts/diagnose_zustandsflaeche.py`, `scripts/p0_baseline.py`, `scripts/validate_zustandsflaeche.py`, `regression/p0_baseline.json`, `results/diagnose_zustandsflaeche.json`, `results/validate_zustandsflaeche.json`, `figures/zustandsflaeche_water_*.png`, `.gitignore`, diese Datei.
**Geändert:** `homogeneous_pvt.py` (Achsenliste), `visualization.py` (`to_xyz` über `pressure_to_z`, gemeinsame Vertextabelle, Ersatzfläche, Ticks, Legende, Hover‑Picking), `main.py` (Rohfeld, Ersatzfläche, Achsenumschaltung, neue `daten`‑Keys `replacement_surface`, `replacement_field`, `replacement_status`, `homogeneous_raw_field`, `pressure_axis_notice`), `gui.py` (Ebene, Exklusivität, Export, Hover, Hinweise), `config.yaml` (`include_replacement`, `replacement_points`, `auto_pressure_axis`), `test_integration_mocked_main.py` (Mock für `build_raw_field`).
