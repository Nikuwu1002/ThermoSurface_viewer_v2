# AP19 – Korrektur der nach oben ausreißenden instabilen Zustandsfläche

**Stand:** 22.10.2026 · **Basis:** Projektstand nach AP18 (CoolProp 7.2.0, PyVista 0.47.3, NumPy 2.3.4, Python 3.13)
**Ergebnis in einem Satz:** Die Ausreißer sind **keine** Rechen-, Klassifikations-, Phasenvorgabe- oder Vernetzungsfehler, sondern eine **tatsächliche Eigenschaft der Mehrparameter-Zustandsgleichung** (zusätzliche Extrema von p(ρ) zwischen den Spinodalen). Die instabile Fläche enthält jetzt nur noch Zustände, die auf einem zusammenhängenden, mechanisch instabilen Ast zwischen den Spinodalen liegen; der Rest wird als „EOS-Artefaktlücke“ dokumentiert, nicht geglättet oder abgeschnitten.

---

## 1. AP19.1 – Ausgangszustand reproduzieren

```powershell
python diagnose_ap19.py                      # Wasser, 580–646 K, Schritt 2 K, 20 Volumenpunkte
python diagnose_ap19.py --fluid CO2 --T-start 250 --T-stop 303 --out results/ap19_diagnose_co2.json
```

| Größe (Wasser, Standardraster, alte Geometrie) | Wert |
|---|---|
| Knoten der instabilen Fläche | 562 |
| davon außerhalb `[p_spin,l ; p_spin,g]` der eigenen T‑Zeile | **229** |
| Druckbereich der Fläche | **−633 p_c … +1100 p_c** (bei 580 K bis 2,4·10¹⁰ Pa) |
| Dreiecke mit mindestens einem Ausreißerknoten | 451 von 954 |
| größte Dreieckskante (Plot‑Einheiten) mit / ohne Ausreißerknoten | 3,41 / 0,44 |

Bildvergleich (nur Metastabil‑ und Instabilflächen, `asinh`-Druckachse, Offscreen‑PyVista): [Wasser vorher](../../03_Tests_Implementierung/figures/ap19_water_vorher.png) · [Wasser nachher](../../03_Tests_Implementierung/figures/ap19_water_nachher.png) · [CO₂ vorher](../../03_Tests_Implementierung/figures/ap19_co2_vorher.png) · [CO₂ nachher](../../03_Tests_Implementierung/figures/ap19_co2_nachher.png).
Die rohen Knotenkoordinaten der schlimmsten Dreiecke (`(T,v,p)` je Knoten) stehen in `results/ap19_diagnose_water.json` unter `D_worst_triangles_original_TvP`.

## 2. Ursachenanalyse (AP19.3 / AP19.4)

| Hypothese | Prüfung (unabhängig von der Visualisierung) | Befund |
|---|---|---|
| **A** fehlerhafte HEOS‑Zustandswerte | Druck direkt aus der IAPWS‑95‑Helmholtzfunktion (`iapws`, eigene Implementierung, kein CoolProp, kein Flash) bei 600 K und allen 20 Knoten | **identisch** (Abweichung < 10⁻⁴ p_c; Wasser 600 K: p(ρ=21373 mol/m³)=1163,9 MPa in beiden) |
| **B** ungeeignete Phasenvorgabe / ρ‑T‑Update | `specify_phase(liquid)` vs. `specify_phase(gas)`; Wert gegen Helmholtz‑Rechnung | Druck **unabhängig** von der Phasenvorgabe; kein Flash |
| **C** falsch klassifizierte Zustände | Vorzeichen von dp/dρ und Extrema von p(ρ) auf 4000 Dichten zwischen den Spinodalen (600 K) | p(ρ) hat **zwei zusätzliche Extrema** zwischen den Spinodalen: Minimum −31,2 p_c bei ρ=16 729 mol/m³, Maximum +57,1 p_c bei ρ=22 152 mol/m³. Ein van‑der‑Waals‑artiger Verlauf hätte keines. Die Vorzeichenklassifikation war korrekt, aber sie sah nur das Vorzeichen, nicht den Druckbereich |
| **D** fehlerhafte PyVista‑Vernetzung | Kantenlängen aller Dreiecke, getrennt nach „enthält Ausreißerknoten“ / „enthält keinen“ | Alle langen Kanten gehen von Knoten außerhalb des Druckfensters aus (3,41 vs. 0,44); die Vernetzung verbindet nur vorhandene Knoten und erzeugt keine eigenen Ausreißer |

Gleiches gilt für CO₂, Argon und Methan (HEOS‑Mehrparameter‑Gleichungen; Stickstoff wurde nur mit der Legacy‑Geometrie reproduziert); bei R134a ist der Effekt schwächer (siehe Validierung). Die Gleichungen sind nur im stabilen/metastabilen Bereich an Messdaten angepasst; ihr Verlauf in der instabilen Zone ist nicht Teil der Anpassung und wird nicht durch Messdaten gestützt. Das ist eine tatsächliche EOS‑Eigenschaft und nach Aufgabenstellung (Abnahmekriterium 1) ausdrücklich zulässig einzuordnen – aber die Fläche darf diese Bereiche nicht als „gültige instabile Fortsetzung“ ausgeben.

Warum der alte Aufbau trotzdem Ausreißer akzeptierte: `build_homogeneous_surfaces()` verlangte nur das Vorzeichen `dp/dρ<0` je Knoten. Nach dem ersten zusätzlichen Extremum wechselt das Vorzeichen mehrfach; Knoten im hinteren „instabilen“ Teilbereich (z. B. ρ≈22 000–26 000) bestehen den Test, haben aber Druck bis +1100 p_c, also weit oberhalb von `p_spin,g`.

## 3. Korrektur (AP19.5) – `homogeneous_pvt.py`

**Akzeptanzkriterium für Knoten der instabilen Fläche** (`unstable_node_status`):
1. `(∂p/∂v)_T > 0` (äquivalent `∂p/∂ρ < 0`), streng,
2. `p_spin,l ≤ p ≤ p_spin,g` (Toleranz 10⁻⁹·max(p_c,|p_l|,|p_g|)) – exakt die Forderung aus AP19.2,
3. entlang jedes Astes steigt p mit v (Monotonieprüfung der akzeptierten Knoten),
4. der Knoten ist **vom jeweiligen Spinodalast aus zusammenhängend erreichbar**.

**Ast‑Bestimmung** (`trace_unstable_row`): je T‑Zeile werden 256 Punkte (`unstable_scan_points`) zwischen `v_spin,l` und `v_spin,g` klassifiziert. Vom Flüssigkeitsast aus wird das erste Ende des gültigen Bereichs (Vorzeichenwechsel *oder* Verlassen des Druckfensters) durch Bisektion auf 10⁻¹³ relativ genau bestimmt (`v_a`), ebenso rückwärts vom Dampfast (`v_b`).

* **Reguläre Isotherme** (Wasser ab ≈ 643,7 K ≈ 0,995 T_c, CO₂ ab ≈ 302,4 K, Argon ab ≈ 148,4 K, R134a ab ≈ 348,6 K ≈ 0,93 T_c; Raster `validate_ap19_ap20.py`, hoch): beide Äste decken `[v_l, v_g]` vollständig ab → Zeile ist **bitgleich** zum alten Raster `geomspace(v_l, v_g, n)`.
* **Isotherme mit EOS‑Artefakt**: `[v_l, v_a]` und `[v_b, v_g]` werden je mit eigener Knotenzahl vernetzt; dazwischen liegt eine **NaN‑Trennspalte** (Spalte `⌊(n−1)/2⌋`). Es wird **nicht** über die Lücke interpoliert oder trianguliert. Ursache, `v_a`, `v_b`, Anzahl der Vorzeichenwechsel und der p‑Bereich des Scans stehen in `reports["unstable_surface"]["gap_rows"][<T‑Index>]`.

**Weitere Änderungen:**

* `solve_spinodal_branch()` wählt bei mehreren Extrema die **sättigungsnächste Nullstelle** (erste von der Siedelinie/Taulinie aus erreichte), statt die Zeile zu verwerfen. Die Temperaturen mit Mehrfachextrema bleiben in `ambiguous_temperatures` protokolliert (jetzt rein informativ). Wirkung: Wasser 57/59 → 59/59, CO₂ 25/29 → 29/29, R134a 6/26 → 26/26 gültige Flüssigkeitsspinodal‑Zeilen (jeweils gleiches Testraster vor/nach der Änderung).
* `build_homogeneous_surfaces()`: Grenzordnung `v_sat,l < v_spin,l < v_spin,g < v_sat,g` und `p_spin,l < p_spin,g` je Zeile (AP19.2); Zeilen, die das verletzen, werden für alle Flächen gesperrt (`order_violations`). Metastabile Flächen: zusätzlich Druckfenster `[p_sat, p_spin]` und Monotonie (dp/dv<0).
* `validate_unstable_nodes()` (läuft automatisch im Build, `reports["unstable_surface"]["independent_check"]`): jeder akzeptierte Knoten wird unabhängig geprüft – Spinodalgrenzen, Druckfenster, Druck bei erzwungener Flüssig‑ **und** Gasphase, `dp/dv` aus Differenzenquotienten der p(v)‑Werte (nicht aus `first_partial_deriv`; an Astenden nach innen gerichtet), sowie `v_sat,l<v_spin,l`, `v_spin,g<v_sat,g`.
* `validate_homogeneous_boundaries()` und `Debug.diagnose_homogeneous()` vertragen jetzt die NaN‑Trennspalte (vorher hätte der Debug‑Modus mit „Volumen nicht monoton“ abgebrochen).
* `visualization.py`: Triangulation (`grid_triangles`, `drop_degenerate`) in das VTK‑freie Modul ausgelagert, damit Renderer und Tests dieselbe Logik nutzen; `homogeneous_mesh_diagnostics()` protokolliert Dreiecke, Komponenten, offene Kanten, Orientierung, Projektionsfaltungen. `flag_triangles()` listet Dreiecke mit Knoten außerhalb des Fensters samt originalen `(T,v,p)`.

**Nicht getan (ausdrücklich):** kein Abschneiden großer Drücke, keine Glättung, keine Interpolation, keine Entfernung gültiger Werte nur wegen ihrer Größe. Die Ausreißerknoten sind echte EOS‑Werte; sie liegen nur nicht auf dem zusammenhängenden Ast zwischen den Spinodalen und werden deshalb nicht Teil *dieser* Fläche. (Die Rohdaten‑Reproduktion bleibt über `regression/ap19_vorher/homogeneous_pvt_legacy.py` möglich.)

## 4. Validierung (AP19.6)

Drei Rasterauflösungen (niedrig 13×10, mittel 29×20, hoch 59×40 Knoten; T‑Raster nach T_c verdichtet, bis ΔT/T_c=4·10⁻⁴) und fünf Stoffe:
`python validate_ap19_ap20.py --fluids Water CO2 Argon Methane R134a --out results/ap19_ap20_validation.json`

| Stoff | Aufl. (n_T,n_v) | T-Zeilen vollst. | Vorher: Knoten außerhalb p-Fenster | max. Überschreitung vorher / p_c | **Nachher: außerhalb Fenster** | unabh. Prüfung (Fehler) | reguläre Zeilen / EOS-Artefaktzeilen |
|---|---|---|---|---|---|---|---|
| Water | niedrig (13,10) | 13/13 | 18 von 108 | 939 | **0** von 119 | 0 | 6 / 7 |
| Water | mittel (29,20) | 29/29 | 84 von 472 | 3,14e+04 | **0** von 531 | 0 | 13 / 16 |
| Water | hoch (59,40) | 59/59 | 361 von 1895 | 1e+05 | **0** von 2171 | 0 | 25 / 34 |
| CO2 | niedrig (13,10) | 13/13 | 23 von 107 | 8,61e+03 | **0** von 118 | 0 | 5 / 8 |
| CO2 | mittel (29,20) | 29/29 | 100 von 458 | 2,52e+04 | **0** von 535 | 0 | 10 / 19 |
| CO2 | hoch (59,40) | 59/59 | 432 von 1866 | 3,63e+04 | **0** von 2210 | 0 | 21 / 38 |
| Argon | niedrig (13,10) | 12/13 | 2 von 110 | 5,25 | **0** von 113 | 0 | 8 / 4 |
| Argon | mittel (29,20) | 27/29 | 15 von 481 | 13,4 | **0** von 514 | 0 | 18 / 9 |
| Argon | hoch (59,40) | 56/59 | 76 von 1991 | 16,7 | **0** von 2132 | 0 | 35 / 21 |
| Methane | niedrig (13,10) | 13/13 | 5 von 117 | 1,47 | **0** von 121 | 0 | 7 / 6 |
| Methane | mittel (29,20) | 29/29 | 25 von 515 | 2,4 | **0** von 546 | 0 | 15 / 14 |
| Methane | hoch (59,40) | 59/59 | 114 von 2082 | 2,86 | **0** von 2234 | 0 | 31 / 28 |
| R134a | niedrig (13,10) | 13/13 | 1 von 125 | 0,928 | **0** von 123 | 0 | 11 / 2 |
| R134a | mittel (29,20) | 29/29 | 6 von 549 | 10,8 | **0** von 545 | 0 | 25 / 4 |
| R134a | hoch (59,40) | 59/59 | 34 von 2226 | 44,4 | **0** von 2210 | 0 | 50 / 9 |

(„Vorher“ = alter Netzaufbau, aber mit den hier verwendeten Spinodalen, damit nur die Flächengeometrie verglichen wird.)
Zusätzlich: `python validate_ap17_ap18.py` → `HEOS_CORE_PASS` für Water, R134a, CO2, Argon (PyVista‑Zellflächen > 0 und endlich); Ende‑zu‑Ende‑Läufe von `main.berechne_daten()` für Water, CO2, Argon, R134a (auch mit `DEBUG_MODE=True`).

## 5. Abnahmekriterien

| Kriterium | Stand |
|---|---|
| Ausreißer ursachengerecht beseitigt oder als EOS‑Eigenschaft eingeordnet | ✔ eingeordnet (Abschnitt 2, IAPWS‑95‑Gegenprobe) und aus der Fläche entfernt (Abschnitt 3) |
| alle akzeptierten Knoten innerhalb der Spinodalgrenzen | ✔ 0 Verletzungen in 15 Fällen; je Knoten unabhängig geprüft |
| mechanische Klassifikation korrekt | ✔ strikt `dp/dv>0`; Gegenprobe über Differenzenquotienten |
| keine ungültigen/degenerierten Dreiecke | ✔ `drop_degenerate`, 0 Orientierungskonflikte, 0 Projektionsfaltungen |
| Vorher‑Nachher‑Vergleiche, reproduzierbare Tests | ✔ `diagnose_ap19.py`, `validate_ap19_ap20.py`, `test_ap19_ap20.py` |

## 6. Dateien

Geändert: `homogeneous_pvt.py`, `visualization.py`, `main.py`, `Debug.py`, `gui.py` (Tooltips), `config.yaml`, `test_integration_mocked_main.py` (Mock für Zusatzzeilen), `validate_ap17_ap18.py` (fehlender Import `count_renderable_triangles`, hatte alle Stoffe mit `NameError` scheitern lassen).
Neu: `diagnose_ap19.py`, `validate_ap19_ap20.py`, `test_ap19_ap20.py`, `regression/ap19_vorher/homogeneous_pvt_legacy.py`, `results/ap19_*.json`, `figures/ap19_*.png`, diese Datei.

## 7. Tests

```powershell
python -m unittest test_ap19_ap20 test_ap13_ap18 test_corrections test_spinodal_extension test_integration_mocked_main -v
```
51 Tests, davon 1 vorbestehend übersprungen. `test_ap19_ap20.py` enthält eine synthetische EOS mit `dp/dρ = K(ρ−1)(ρ−3)(ρ−6)(ρ−8)` (zusätzliches Extremapaar außerhalb des Druckfensters), reproduziert darin die Ausreißer mit der Legacy‑Geometrie, belegt ihre Beseitigung, prüft Trennspalte, Randknoten, Fensterkanten, die Manipulationserkennung der Knotenprüfung und – bei installiertem CoolProp – Wasser/CO₂ (3 Auflösungen) sowie Argon/Methan.

## 8. Grenzen und offene Punkte

* Bei tiefen Temperaturen besteht die instabile Fläche für Mehrparameter‑EOS nur aus schmalen Streifen entlang beider Spinodalen; dazwischen fehlt sie **bewusst**. Wer die EOS‑Rohwerte der Artefaktzone sehen will, nutzt die Legacy‑Routine (kein Produktivpfad).
* Die Scan‑Auflösung (256 Punkte) kann ein extrem schmales Extremapaar zwischen zwei Scanpunkten übersehen; jeder Knoten wird zusätzlich einzeln geprüft (Fenster, Vorzeichen, Monotonie, unabhängige Nachprüfung), eine so übersehene Zone würde dort als Punktfehler auftauchen. Ein Test dafür mit einer echten EOS liegt nicht vor.
* Die Spinodal‑Zweigprüfung „Unplausibler Astverlauf“ (Volumenmonotonie in T, 10⁻³) bleibt unverändert. Sie verwirft bei Argon 2–3 Flüssigkeitsspinodal‑Zeilen (EOS‑Welligkeit bei ≈131 K, die Wurzel springt auf ρ/ρ_c=1,24 bei p≈−57 MPa); diese Zeilen bleiben NaN, nicht interpoliert.
* `pressure_axis: log` blendet Knoten mit p≤0 aus (bestehendes Verhalten); für vollständige Darstellung `asinh` wählen.
* Die Qt‑Oberfläche wurde importiert, aber **nicht interaktiv bedient**; getestet wurden Berechnung, `Debug`‑Pfad und PyVista‑Offscreen‑Rendering.
