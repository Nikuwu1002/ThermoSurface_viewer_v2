# TS-A / TS-B – Änderungs- und Abnahmeprotokoll

Grundlage: `ThermoSurface_Implementierungsplan_final.pdf` und `Plan(1).pdf`.
Umgebung der Prüfung: Windows, Python 3.13.9, CoolProp 7.2.0, pyvista 0.47.3, PySide6 6.12.0, pyvistaqt, h5py 3.16.0 (nur für den Plattencache-Test nachinstalliert).
Ausgangszustand: `Doku/TS_A01_Ausgangszustand.md`.

## 1. Freigabestatus

| Paket | Status | Bemerkung |
|---|---|---|
| TS-A01 … A06 | abgenommen (Entwicklungsrechner) | Kein ausdrücklicher Freigabeschritt zwischen Gate A und B (Auftrag: A und B in einem Lauf). Interaktives Bedienen auf dem Zielrechner steht aus. |
| TS-B01 … B06 | abgenommen (Entwicklungsrechner) | s. u. |
| TS-B07 | abgenommen mit Einschränkungen | Vorhandene Fehlschläge siehe Abschnitt 5. Echte Mausbedienung durch einen Menschen steht aus. |

## 2. Neue und geänderte Dateien

**Neu (A):** `diagram_registry.py` (Diagramm-Registry, Layer-Definitionen), `diagram_cache.py` (DiagramKey, LRU-Cache mit Budget), `calculation_manager.py` (sequentielle Berechnung im QThread, Cache-Abfrage im GUI-Thread, Teilfehler), `diagram_viewport.py` (ein Plotter je Diagramm, Layer/Gitter/Hover/Klick), `diagram_workspace.py` (Kachel-Layout, Maximieren/Zurück), `tile_layout.py`, `scripts/ts_a01_referenzlauf.py`, `scripts/ts_a06_performance.py`, `test_gui_multiview.py`.

**Neu (B):** `diagram_coords.py` (Hin- und Rücktransformation der Darstellungskoordinaten), `thermo_state.py` (ThermoState, Statuscodes), `state_store.py` (Quellflächen-Adapter, KD-Index, RAM-LRU, optionaler HDF5-Cache), `state_picker.py` (vtkCellPicker, Baryzentrik nur zur Lokalisierung), `state_resolver.py` (Rekonstruktion über HEOS), `selection_manager.py`, `comparison_view.py`, `test_state_comparison.py`, `test_state_comparison_gui.py`.

**Geändert:**
- `gui.py`: Modus Einzel/Mehrfach, Diagramm-Checkboxen, Workspace statt Einzel-Plotter, Layer-Panel je Viewport, Gruppe „Zustandsvergleich“.
- `main.py`: Berechnung über `BERECHNUNGEN`, Cache-Parameter je Diagrammtyp (nur tatsächlich genutzte Parameter).
- `visualization.py`: die fünf `to_xyz*`-Transformationen nach `diagram_coords.py` verschoben und wieder importiert (Namen bleiben exportiert). Renderer unverändert.
- `limits.py`, `config.yaml`: `limits.cache_max_mb`, `limits.cache_warn_mb`, Abschnitt `state_cache` (Plattencache standardmäßig **aus**).
- `Doku/requirements.txt`: `h5py` als optional vermerkt.

Nicht verändert: EOS- und Thermodynamikmodule (`coolprop_heos_pvt.py`, `zustandsflaeche.py` usw.), bestehende Tests.

## 3. Schnittstellen-Kompatibilität

- `main.berechne_auswahl(stoff, raster, visu_raster, diagram_type, progress_callback)` und die `render_*`-Funktionen unverändert (Signatur und Rückgabe).
- Einzelmodus ist Standard und liefert denselben Bericht, dieselben Tooltips und Layer wie zuvor (`test_gui_multiview.py::test_single_mode_is_default_and_unchanged`).
- Die HEOS-/Ersatzflächen-Herkunft bleibt bis in den ThermoState erhalten (`source`-Codes der Ersatzfläche, Status `synthetic_substitute`).

## 4. Abnahmetabelle

Testlauf: `python -m pytest test_gui_multiview.py test_state_comparison.py test_state_comparison_gui.py` ohne `QT_QPA_PLATFORM=offscreen`: **64 bestanden**, 56 Subtests.

| ID | Prüfszenario / Eingaben | Erwartet | Test | Ergebnis |
|---|---|---|---|---|
| A-T01 | Einzelmodus, pvt, pht, pst, hsp, uvt nacheinander (Water, CO2) | Jeder Renderer und die Layer-Checkboxen funktionieren | `test_single_mode_all_types_and_fluid_change`, `scripts/ts_a01_referenzlauf.py` | bestanden |
| A-T02 | 2, 3, 4, 5 Diagramme | Raster 1×2, 2×2, 2×2, 2×3; Ansichten unabhängig | `test_layouts_independence_maximize_cache_and_errors`, `GridShapeTests.test_fixed_automatic_arrangement` | bestanden |
| A-T03 | Doppelklick, Zurück, anderes Diagramm | Kein `berechne_auswahl`-Aufruf | gleicher Test (Zähler `compute_calls`) | bestanden |
| A-T04 | Stoff, Raster, Layer, Kamera ändern | Nur Stoff/Raster → Cache-Miss | `CacheKeyRelevanceTests`, Stoffwechsel-Test | bestanden |
| A-T05 | Ein Diagramm schlägt fehl | Übrige Ansichten bleiben nutzbar | `test_failed_diagram_does_not_block_others`, `test_layouts_independence_maximize_cache_and_errors` | bestanden |
| B-T01 | Gitterknoten in allen 5 Diagrammen, Water + CO2 | Exakte Quelldaten, molare Einheit, Herkunft, Fluid | `test_exact_nodes_all_diagrams_reproduce_source_data` | bestanden |
| B-T02 | Punkte zwischen Stützstellen (einphasig, alle Gebiete) | Rekonstruktion gegen direkten HEOS-Aufruf, Plausibilitätsprüfung, Status | `test_single_phase_interior_points_all_regions`, `test_pvt_single_phase_interior_matches_direct_heos`, `test_interior_hit_has_valid_weights_and_phase_information` | bestanden |
| B-T03 | Zweiphasenpunkt Q = 0, 0,5, 1 | Konsistente molare Größen, Phasenkennung | `test_saturated_liquid_vapor_and_mixture_with_given_Q`, `test_T_and_p_alone_at_saturation_do_not_define_Q` | bestanden |
| B-T04 | Ersatzfläche | Als synthetisch markiert, h/s/u/f/g/ρ/Q = None | `test_substitute_surface_gets_no_invented_heos_values`, `test_synthetic_state_cannot_carry_heos_energies` | bestanden |
| B-T05 | Drei Punkte aus zwei Diagrammtypen, Achsen wechseln | Verbindung/Beschriftung korrekt, keine Neuberechnung | `test_three_points_two_diagram_types_switch_axes`, `test_axis_change_does_not_touch_states`, GUI: `test_click_selects_marks_and_double_click_only_maximizes` | bestanden |
| B-T06 | HDF5-Datei nach Neustart | Gleiche Datei, Gültigkeit nach Schema/CoolProp-Version/Einheiten/Raster | `test_single_hdf5_file_is_updated_and_validated`, `test_grid_hash_changes_with_grid_and_fluid` | bestanden (h5py 3.16.0) |
| B-T07 | Anwendung schließen | Auswahl nicht persistent, Threads/Plotter freigegeben | `test_selection_is_not_persisted`, GUI: `test_exit_releases_threads_views_and_persists_no_selection` | bestanden |
| B07 | Ersatzfläche, negativer Druck, kritische Nähe, nicht berechenbarer Punkt, fehlerhafter Datensatz, T außerhalb | Eindeutiger Status, keine Scheinpräzision | `test_negative_pressure_is_kept_but_flagged`, `test_critical_nearness_is_flagged`, `test_not_computable_point_gets_status_and_no_values`, `test_corrupted_dataset_is_reported_not_silently_fixed`, `test_out_of_range_temperature_is_rejected` | bestanden |
| B07 | Vergleichsmodus aus / an | Aus: keine Auswertung. An: Wiederholung aus Cache | `test_inactive_comparison_does_nothing_expensive`, `test_repeated_query_costs_at_most_one_heos_evaluation`, GUI `test_picking_inactive_until_mode_enabled` | bestanden |
| B07 | Mindestens drei Punkte sortieren/löschen/zurücksetzen, Stoffwechsel | Nummern eindeutig, Marker nur im Ursprungsdiagramm | GUI-Tests, `test_numbering_origin_color_add_delete_reset`, `test_reorder_changes_only_connection_order` | bestanden |
| B07 | Speicherbudget | LRU, Warnung | `test_lru_budget`, `CacheTests` | bestanden |

## 5. Gesamtlauf bestehender Tests

`QT_QPA_PLATFORM=offscreen python -m pytest`: 160 bestanden, 11 übersprungen (GUI-Tests brauchen ein echtes OpenGL-Fenster), **4 fehlgeschlagen** – dieselben vier wie vor den Änderungen (Ausgangslage 106 bestanden / 4 fehlgeschlagen):

- `test_zustandsflaeche.py::P8EndeZuEnde::test_regression_against_p0_baseline_is_bitwise_equal` (9733 ≠ 9614) – bereits vor TS-A vorhanden, nicht verursacht und nicht behoben.
- `P6Druckachse::test_mesh_keeps_all_signed_pressure_states…`, `P7Kennzeichnung::test_all_substitute_vertices_carry_label…`, `P7Kennzeichnung::test_pick_hover_text_offscreen`: schlagen nur im Gesamtlauf fehl, weil `test_ap13_ap18` pyvista durch einen Stub ersetzt (Testverschmutzung); einzeln laufen sie.

## 6. Messwerte (`results/ts_a06_performance.json`, Water)

| | Standardraster | Feines Raster |
|---|---|---|
| pvt Berechnung / Render | 3,1 s / 1,2 s | 4,6 s / 1,3 s |
| pht, pst, hsp, uvt je | 0,1 – 0,6 s | 0,2 – 1,2 s |
| RAM-Spitze bei 5 Viewports | ≈ 660 MB | ≈ 694 MB |
| Cache-Größe (5 Diagramme) | 4,2 MB | 7,8 MB |
| Wiederholung | 0 neue Berechnungen | 0 neue Berechnungen |

Jeder zusätzliche Viewport kostet etwa 40 MB.

## 7. Bekannte Einschränkungen

- Der Fächer zum kritischen Punkt der Zweiphasenfläche ist `not_unique` (Q dort nicht definiert); der kritische Knoten selbst ist gültig.
- Zustände der Ersatzfläche enthalten nur T, v, p (Modellwert); h, s, u, f, g, ρ, Q sind None.
- Auf der pvt-Gleichgewichtsfläche liegen Plateau-Knoten im Zweiphasengebiet; dort wird Q aus HEOS bestimmt.
- Die Plausibilitätsprüfung der Rekonstruktion toleriert 25 % der Dreiecksausdehnung; darüber hinaus Status `partially_valid`.
- Energiegrößen stehen im CoolProp-Referenzzustand; h, s, u werden aus J/kg über die Molmasse in J/mol umgerechnet.
- Metastabile/instabile Zustände sind homogene HEOS-Fortsetzungen (Phase vorgegeben, nicht als Gleichgewicht zu lesen).
- Der Plattencache ist standardmäßig aus (`state_cache.disk_enabled: false`) und braucht `h5py`; ohne h5py läuft alles im RAM.
- `QT_QPA_PLATFORM=offscreen` unterstützt QtInteractor nicht; GUI-Tests werden dort übersprungen.
- `Doku/manifest_sha256.json` ist nicht aktualisiert.

## 8. Auf Zielrechner offen

- Interaktive Bedienung durch den Anwender (Mausgefühl beim Picken, Marker-Sichtbarkeit bei Rotation, Maximieren mit der echten GPU-Konfiguration).
- Weitere Reinstoffe außer Water und CO2 und große Raster mit Plattencache bei > 1000 Zuständen.
- Freigabeentscheidung: [ ] abgenommen   [ ] Nacharbeit erforderlich   [x] auf Zielrechner ausstehend (interaktive Punkte)
