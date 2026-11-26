# TS-VP07 – Integration, Performance, Regression und Dokumentation

**Stand:** 16.11.2026 · **Umgebung:** Windows, Python 3.13, CoolProp 7.2.0, PyVista 0.47, PySide6 (native Qt-Plattform mit OpenGL für die GUI-Tests)

## 1. Testmatrix (Spezifikation Abschn. 4) – reproduzierbar

Aufruf je Datei (GUI-Dateien ohne `QT_QPA_PLATFORM=offscreen`, die übrigen mit): `python -m pytest 03_Tests_Implementierung/<datei> -q -p no:cacheprovider`

| Test-ID | Nachweis | Datei :: Test | Ergebnis |
|---|---|---|---|
| L01 | 10…1000, f = 3.14 → 10 / 31.4 / 98.596 / 309.59144 / 972.1171216 | `test_ts_vp02_levels.py::GeneratorFormula::test_L01…` | bestanden |
| L02 | start ≤ 0, f ≤ 1, stop < start, NaN/Inf | `…::GeneratorValidation::test_L02…` | bestanden |
| L03 | f ≈ 1, `max_levels`-Schutz, Überlauf | `…::test_L03…`, `test_max_levels…`, `test_overflow…` | bestanden |
| L04 | T-/v-/p-Residuen entlang jeder Kurve (Wasser, CO₂; 7 Diagramme) | `test_ts_vp03_isolines.py::ConstantAlongCurve` | bestanden (273 Subtests) |
| L05 | Isolinien-Actors in allen sieben Diagrammen | `…::RendererActors` | bestanden |
| V01 | Checkbox aus → kein Clipping, physikalische Labels | `test_ts_vp_gui.py::…test_controls_follow_the_axis_matrix_and_start_disabled` | bestanden |
| V02 | gültige sechs Werte + Anwenden → Clipping und manuelle Achsengrenzen | `…test_V02_apply_clips_original_geometry_without_any_recomputation` | bestanden |
| V03 | Grenzen umstellen, wieder erweitern/zurücksetzen → Original vollständig | `…test_V03_reset_restores_full_original_geometry`, `ActorRegistry::test_clip_restore_and_repeat…` | bestanden |
| V04 | leeres Fenster → Rahmen/Labels bleiben, keine Renderfehler | `…test_V04_empty_window_keeps_frame_and_axes`, `ActorRegistry::test_V04…` | bestanden |
| V05 | Layer-Toggle → Box und Kamera unverändert | `…test_V05_layer_toggles_keep_box_and_camera` | bestanden |
| V06 | Mehrfachansicht, zwei verschiedene Ausschnitte | `…test_V06_V07_independent_viewports_survive_maximize_and_mode_switch` | bestanden |
| V07 | Zoom/Maximieren/Moduswechsel | ebenda | bestanden |
| V08 | negative Energieachse (h, s, u, g, f) | `test_ts_vp04_05_06_viewport.py::BoundsValidation::test_V08…` | bestanden |
| V09 | p ≤ 0 mit Log/Asinh/Symlog | `…::test_V09…` | bestanden |
| V10 | Clipping bei NaN-Lücken → keine Brücken | `…::TriangleClipping::test_V10…` | bestanden |
| V11 | Hover/Picker/Vergleichsmarker an geclippten Rändern | `test_ts_vp_gui.py::…test_V11_picking…`, `…test_markers_outside_the_box…`, `ActorRegistry::test_new_points_are_flagged…`, `…test_pick_geometry…` | bestanden |
| C01 | Viewport mehrfach ändern → keine neue HEOS-Berechnung, gleiche Rohdaten/Cache-ID | `…test_V02…` (Aufrufzähler, `daten is daten`, `DiagramKey`, Arrays), `scripts/ts_vp07_performance.py` (0 Zustandserzeugungen, Hash gleich) | bestanden |
| C02 | Isolinienlevel ändern → neuer Cache-Schlüssel | `test_ts_vp03_isolines.py::CacheIdentity`, `test_ts_vp_gui.py::…test_C02…` | bestanden |

## 2. Prüfumfang (Aufgabenliste TS-VP07)

1. **Diagrammtypen/Stoffe:** alle sieben; Wasser und CO₂ (L04, L05, Levelvergleich), zusätzlich R134a und Argon (`FluidMatrix`, Standardisolinien); bestehende Matrix `test_fluid_diagram_matrix.py` (10 GUI-Stoffe).
2. **Isolinien:** exakte Formel, Endpunktverhalten, Größenordnungen (T, v, p), Werte nahe `factor = 1` (L02/L03).
3. **Viewports:** vollständig innen (`box ⊇ Daten`: alle Dreiecke erhalten, bitgleich), Teil-/Randschnitt (Flächeninhalt exakt, Kantenpunkte auf der Boxfläche), leeres Fenster, negativer Energiebereich, Druck-Sonderfälle.
4. **Cacheidentität:** unverändert bei Viewportwechsel (C01), neu bei geänderter Isolinienkonfiguration (C02), gleich bei Vorgabe = config.yaml.
5. **Picker/Zustandsvergleich:** Treffer auf beschnittenen Flächen (Ursprungsdreieck, `clipped`), Marker außerhalb der Box ausgeblendet und beim Zurücksetzen wieder sichtbar; Hover-Kennzeichnung der Schnittpunkte.
6. **Layer/Kamera/Mehrfachlayout:** Layer inkl. homogener/instabiler/Ersatzflächen (Actor-Schlüssel unverändert; `test_gui_multiview.py`, `test_gpt_fvt.py`), Show-Edges (`apply_grid` unverändert), Kamera (V05), Mehrfachlayout (V06/V07).
7. **Bitgleichheit:** Rohdatenarrays vor/nach mehrfachen Bereichsänderungen (`test_clip_restore_and_repeat_leave_raw_data_untouched`, GUI-V02, Performance-Skript: SHA-1 aller Arrays identisch).
8. **Reaktionszeit/Speicher:** Abschn. 3.
9. **Dokumentation:** `README_TS-VP01.md` … `README_TS-VP07.md`, `README_TS-VP_GESAMT.md`; Beispielbilder `03_Tests_Implementierung/figures/ts_vp07_*`.

## 3. Reaktionszeit und Speicher des Clippings (ohne CoolProp-Neuberechnung)

`python 03_Tests_Implementierung/scripts/ts_vp07_performance.py` → `03_Tests_Implementierung/results/ts_vp07_performance.json` (Wasser, Standardraster, Offscreen-Rendering 900×700; 20 zufällige Bereichswechsel je Diagramm inkl. Achsen und Render, danach Zurücksetzen). Während der Messung ist `coolprop_heos_common._new_state` instrumentiert: **0 HEOS-Zustandserzeugungen**.

| Diagramm | Dreiecke/Zellen (Original) | Bereichswechsel Median [ms] | max [ms] | Zurücksetzen [ms] | Python-Speicherspitze [MB] | Geometrie Original → letztes Fenster [MB] | einmalige Berechnung [s] |
|---|---|---|---|---|---|---|---|
| pvt | 180 364 | 315 | 456 | 23 | 45.7 | 8.61 → 1.32 | 2.9 |
| pht | 15 783 | 130 | 152 | 20 | 2.1 | 0.44 → 0.18 | 0.15 |
| pst | 15 783 | 107 | 152 | 20 | 2.0 | 0.44 → 0.10 | 0.14 |
| hsp | 15 783 | 128 | 156 | 21 | 2.1 | 0.44 → 0.36 | 0.40 |
| uvt | 15 580 | 118 | 137 | 21 | 2.1 | 0.43 → 0.14 | 0.14 |
| gpt | 15 585 | 127 | 150 | 21 | 2.1 | 0.44 → 0.00 | 0.58 |
| fvt | 20 465 | 121 | 166 | 21 | 2.0 | 0.57 → 0.01 | 0.52 |

Prozess-Spitzenspeicher des gesamten Messlaufs (alle sieben Berechnungen, Renderer, Clipping): ca. 430 MB. Das reine Clipping (NumPy, ~0.25 s für 351 000 Dreiecke) ist kleiner als der Render-/Achsenanteil der gemessenen Zeit. Das Original bleibt zusätzlich einmal im Speicher (Referenz im `ClipRegistry`, geometrische Kopie nur als normalisiertes `PolyData`).
**Begrenzung:** Clipping wirkt nur auf bereits berechnete Daten. Limits außerhalb des berechneten Bereichs ergeben Leerräume; ein feineres Raster im Ausschnitt oder eine Neuberechnung findet nicht statt.

## 4. Gesamtergebnis der Testläufe

| Datei | Ergebnis |
|---|---|
| `test_ts_vp02_levels.py` | 17 bestanden (+12 Subtests) |
| `test_ts_vp03_isolines.py` | 9 bestanden (+273 Subtests) |
| `test_ts_vp04_05_06_viewport.py` | 28 bestanden |
| `test_ts_vp_gui.py` | 12 bestanden |
| `test_gui_multiview.py` | 19 bestanden (Erwartung Layer `isothermen` in p-h-T an TS-VP03 angepasst) |
| `test_state_comparison_gui.py` | 4 bestanden |
| `test_state_comparison.py` | 41 bestanden (+58 Subtests) |
| `test_gpt_fvt.py` | 47 bestanden (+5989 Subtests; GUI-Anteil nativ) |
| `test_fluid_diagram_matrix.py` | 3 bestanden (+43 Subtests) |
| `test_ap13_ap18.py`, `test_ap19_ap20.py`, `test_corrections.py`, `test_two_phase_limits.py`, `test_integration_mocked_main.py`, `test_spinodal_extension.py` | 8 / 29 / 4 / 7 / 1 / 9 bestanden (spinodal: 2 bewusste Skips) |
| `test_zustandsflaeche.py` | 48 bestanden, **1 vorbestehender Fehlschlag** (`…test_regression_against_p0_baseline_is_bitwise_equal`, `single_phase_surface n_points 9733 != 9614`) – mit den gesicherten Originaldateien identisch, unabhängig von TS-VP01…07 |

## 5. Anwenderhinweise, Einheiten, Formeln, Grenzen

* **Isolinien:** `x[n+1] = f·x[n]`, Start/Stop in der GUI (T [K], v [m³/mol], p [MPa]; 0 = automatisch), Faktor f > 1 (Standard 3.14). Bei `T_stop/T_start < 3.14` entsteht nur eine Isotherme – für dichte Scharen f senken (z. B. 1.1). Levels außerhalb des Diagrammbereichs werden ausgewiesen und nicht gezeichnet.
* **Ansichtsbereich:** Checkbox → sechs Felder (Anzeigeeinheiten) → „Ansicht anwenden“; „Zurücksetzen“ stellt die Originaldarstellung her. Rein darstellungsbezogen, keine Neuberechnung.
* **Achsen:** physikalische Werte/Einheiten, Skalenart im Titel (log, asinh, symlog). p ≤ 0 nur mit `pressure_axis` ≠ `log` (`config.yaml: diagram_settings.pvt/gpt.pressure_axis`).
* **Schnittpunkte** sind geometrische Interpolationen im Renderer (Hover: „Schnittpunkt … keine HEOS-Auswertung“); Zustände werden aus dem Ursprungsdreieck rekonstruiert.
* **Bekannte Einschränkungen:** (1) Negative/Null-Druckwerte homogener Zustände bleiben Oberflächendaten, werden aber nie als geometrisch verteilte Isobaren erzeugt. (2) Für h, s, u, f, g existiert keine Progression (nicht definiert). (3) Tick-Kanten werden nach jeder Mausinteraktion, nicht während des Ziehens, neu gesetzt. (4) Die Orientierungswidget-Animation löst kein Tick-Update aus (Labels bleiben lesbar, liegen aber ggf. an der zuvor gewählten Kante bis zur nächsten Interaktion). (5) p-v-T-Rasterfeld „Anzahl Isothermen“ ist wirkungslos, solange `isolines` konfiguriert ist (Standard).

## 6. Beispielbilder

`03_Tests_Implementierung/figures/ts_vp07_<pvt|pht|pst|hsp|uvt|gpt|fvt>_voll.png` und `…_ausschnitt.png` (Wasser; Ausschnitt = mittlere 60 % je Achse der automatischen Grenzen, Isolinienfaktoren 1.12/2/2), erzeugt mit `python 03_Tests_Implementierung/scripts/ts_vp07_beispielbilder.py`.

## 7. Abnahme

Reproduzierbare Testmatrix ohne Regression gegenüber dem Ausgangszustand (einziger Fehlschlag vorbestehend), dokumentierte Laufzeiten, Beispielbilder und vollständige Übergabe (`README_TS-VP01…07`, Gesamtnote).
