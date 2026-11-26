# TS-A01 – Ausgangszustand und API-Verträge (Stand vor Änderungen)

Erstellt am 11.11.2026 mit `scripts/ts_a01_referenzlauf.py` (Rohdaten: `results/ts_a01_referenz.json`).
Umgebung: Python 3.13.9, CoolProp 7.2.0, pyvista 0.47.3, PySide6 6.12.0, pyvistaqt installiert, pytest vorhanden, h5py **nicht** installiert.

## 1. Berechnungs- und Render-Schnittstellen

Einstieg: `main.berechne_auswahl(stoff, raster, visu_raster, diagram_type, progress_callback)` → `dict` (`daten`).
Renderer: `visualization.render_<typ>(plotter, daten)` → Actor-Dictionary. Jeder Renderer ruft `plotter.clear()`, `show_bounds(...)` und `reset_camera()` auf und arbeitet ausschließlich auf dem übergebenen Plotter (geeignet für mehrere Viewports).

| Typ | Punktspalten der Flächen | `critical_point` | Einheiten | Rechenparameter (tatsächlich verwendet) |
|---|---|---|---|---|
| `pvt` | (T, v, p) | (Tc, pc, vc) | K, m³/mol, Pa | komplettes `raster` + `visu_raster` |
| `pht` | (T, h, p) | (Tc, pc, hc) | K, J/kg, Pa | `two_phase_visu.Q_points` |
| `pst` | (T, s, p) | (Tc, pc, sc) | K, J/(kg K), Pa | `two_phase_visu.Q_points` |
| `hsp` | (h, s, p) | (hc, sc, pc) | J/kg, J/(kg K), Pa | `Q_points` + `config.yaml: diagram_settings.hsp` |
| `uvt` | (T, v, u) | (Tc, vc, uc) | K, m³/mol, J/kg | `Q_points` + `config.yaml: diagram_settings.uvt` |

### Rückgabeschlüssel (Referenzlauf Water/CO2, Shapes Water)
- **pvt**: `surface (99,430,3)`, `isothermen`, `isochoren`, `isobaren` (Liste), `binodal_liquid/vapor (n,3)`, `spinodale_liquid/vapor`, `spinodal_validation`, `spinodal_solver`, `critical_point`, `phase_mask`, `two_phase_surface (n_T,n_Q,3)`, `two_phase_reaches_critical`, `Q_values`, `pressure_axis`, `pressure_axis_notice`, `critical_closure`, `replacement_surface`, `replacement_field`, `homogeneous_raw_field`, `replacement_status`, `metastable_liquid_surface`, `metastable_vapor_surface`, `unstable_surface`, `homogeneous_spinodal_liquid/vapor`, `homogeneous_validation`, `homogeneous_mesh_info` (wird vom Renderer ergänzt!).
- **pht**: `surface_liquid`, `surface_vapor`, `surface_high_temperature`, `binodal_liquid`, `binodal_vapor`, `two_phase_surface`, `Q_values`, `critical_point`, `raster`, `meta` (Tc, pc, hc).
- **pst**: wie pht plus `isothermen`, `isobaren`, `isochoren`; `meta` enthält sc.
- **hsp / uvt**: wie pst; `meta` enthält hc/sc bzw. vc/uc.

### Actor-/Layer-Keys der Renderer
Gemeinsam: `surface`, `surface_single_phase` (immer `None`), `zweiphasengebiet`, `binodale`, `critical_point`, `isothermen`, `isobaren`, `isochoren`, `spinodale`.
Zusätzlich nur `pvt`: `metastable_liquid_surface`, `metastable_vapor_surface`, `unstable_surface`, `replacement_surface`, `pressure_ticks`.
Referenzlauf Actor-Anzahl (Water): pvt 1 Fläche/17 Isothermen/27 Isochoren/11 Isobaren/2 Binodalen/2 Spinodalen; pht 3 Flächen, keine Isolinien; pst 3/7/6/5; hsp 3/7/6/5; uvt 3/10/3/4.

### Referenzlauf (Zeiten, Standardraster aus config.yaml)
| Stoff | pvt calc/render [s] | pht | pst | hsp | uvt |
|---|---|---|---|---|---|
| Water | 2.91 / 0.39 | 0.11 / 0.05 | 0.15 / 0.10 | 0.50 / 0.10 | 0.11 / 0.10 |
| CO2 | 1.64 / 0.29 | 0.11 / 0.05 | 0.14 / 0.11 | 0.46 / 0.11 | 0.10 / 0.11 |

Alle 10 Läufe (5 Typen × 2 Stoffe) ohne Fehler.

## 2. Kritische Randfälle / Verträge
- p-v-T verlangt vorher `generate_stoff_grid` („Stoff festlegen“); `pvt_grid_stoff` muss zum bestätigten Stoff passen.
- p-v-T: Druckachse `log/asinh/symlog/linear`; bei p ≤ 0 automatische Umstellung (`pressure_axis_notice`). `render_pvt` wirft `ValueError`, wenn keine gültigen Einphasenzellen existieren.
- `unstable_surface` (HEOS-Rohfortsetzung) und `replacement_surface` (Ersatz, keine EOS) sind exklusiv (`zustandsflaeche.HOMOGENEOUS_CHOICE`). Provenienz bleibt in `source`-Punktdaten (HEOS | Ersatz) erhalten.
- `render_pvt` verändert `daten` (`homogeneous_mesh_info`) und setzt `plotter._thermosurface_meshes`; Wiederverwendung desselben `daten` in mehreren Plottern ist unkritisch (Renderer idempotent).
- Hover (`pick_hover_text`) liest nur Punktdaten `T, v, p_raw, source` (nur p-v-T).
- Zusatzdiagramme haben eigene Raster (Q_points); Spinodalen/Isolinien in `pht` nicht dargestellt.
- Qt: Berechnung im `QThread` (`CalculationWorker`), Rendering im GUI-Thread.

## 3. Baseline der vorhandenen Tests (vor Änderungen)
`python -m pytest -q` (QT_QPA_PLATFORM=offscreen): **106 bestanden, 4 fehlgeschlagen, 2 übersprungen** (34 s).
- `test_zustandsflaeche.py::P8EndeZuEnde::test_regression_against_p0_baseline_is_bitwise_equal`: schlägt auch isoliert fehl (`9733 != 9614 single_phase_surface`); bereits vor TS-A vorhanden (vermutlich CoolProp-/pyvista-Version gegenüber `regression/p0_baseline.json`).
- 3 weitere Fehler (`P6Druckachse`, `P7Kennzeichnung` ×2) treten nur im Gesamtlauf auf, weil `test_ap13_ap18.py` `pyvista` durch Stubs ersetzt (Test-Pollution); isoliert grün.
- Einzelne Testdateien isoliert: alle grün außer der oben genannten Regression.

Diese Befunde sind **nicht** Teil von TS-A; Vergleichsmaßstab für alle weiteren Läufe: keine zusätzlichen Fehler.

## 4. Betroffene Schnittstellen für TS-A02–A06
`main.berechne_auswahl` (Dispatch → Registry), `gui.CalculationWorker`/`ThermoSurfaceWindow` (globaler `plotter`, `actors`, `layer_preferences`, `_last_daten`), `visualization.render_*`, `visualization.pick_hover_text`, `visualization.update_pvt_axes`.
