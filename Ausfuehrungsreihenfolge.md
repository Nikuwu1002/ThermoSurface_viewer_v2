# ThermoSurface – Dokumente nach Wichtigkeit für die Ausführung

Reihenfolge abgeleitet aus der Importkette ab `01_Programm/main.py`. Pfade relativ zum Projektstamm.

## 1. Start und Konfiguration (ohne diese läuft nichts)

| # | Datei | Rolle |
|---|---|---|
| 1 | `01_Programm/main.py` | Einstiegspunkt: lädt Konfiguration und HEOS-Module, startet die GUI |
| 2 | `01_Programm/config.yaml` | Konfiguration (gelesen über `limits.py` und `main.py`) |
| 3 | `01_Programm/requirements.txt` | Abhängigkeiten (CoolProp, PySide6, pyvista, scipy, h5py …) |
| 4 | `01_Programm/limits.py` | Lädt `config.yaml`, von fast allen Rechenmodulen genutzt |

## 2. Rechenkern (Thermodynamik)

| # | Datei | Rolle |
|---|---|---|
| 5 | `01_Programm/coolprop_heos_common.py` | Gemeinsame Basis aller HEOS-Module |
| 6 | `01_Programm/homogeneous_pvt.py` | Homogene p-v-T-Berechnung |
| 7 | `01_Programm/zustandsflaeche.py` | Zustandsfläche (Kern, Achsen, Visualisierung) |
| 8 | `01_Programm/spinodal_vapor_solver.py` | Spinodal-/Dampf-Solver |
| 9 | `01_Programm/coolprop_heos_pvt.py` | p-v-T-Fläche |
| 10 | `01_Programm/isoline_levels.py` | Isolinien-Niveaus |
| 11 | `01_Programm/coolprop_heos_pht.py` | p-h-T |
| 12 | `01_Programm/coolprop_heos_pst.py` | p-s-T |
| 13 | `01_Programm/coolprop_heos_hsp.py` | h-s-p |
| 14 | `01_Programm/coolprop_heos_uvt.py` | u-v-T |
| 15 | `01_Programm/coolprop_heos_fvt.py` | f-v-T |
| 16 | `01_Programm/coolprop_heos_gpt.py` | g-p-T |

## 3. GUI und Darstellung

| # | Datei | Rolle |
|---|---|---|
| 17 | `01_Programm/gui.py` | Hauptfenster, bindet fast alles ein |
| 18 | `01_Programm/visualization.py` | 3D-Rendering |
| 19 | `01_Programm/diagram_registry.py` | Registrierung der Diagrammtypen |
| 20 | `01_Programm/diagram_coords.py` | Koordinaten |
| 21 | `01_Programm/diagram_axes.py` | Achsen |
| 22 | `01_Programm/axes_view.py` | Achsendarstellung |
| 23 | `01_Programm/diagram_viewport.py` | Viewport |
| 24 | `01_Programm/viewport_bounds.py` | Viewport-Grenzen |
| 25 | `01_Programm/viewport_clipping.py` | Clipping |
| 26 | `01_Programm/geometry_clip.py` | Geometrie-Clipping |
| 27 | `01_Programm/diagram_workspace.py` | Arbeitsbereich (Multi-View) |
| 28 | `01_Programm/tile_layout.py` | Kachel-Layout |
| 29 | `01_Programm/diagram_cache.py` | Cache |
| 30 | `01_Programm/calculation_manager.py` | Berechnungsverwaltung |

## 4. Zustandsauswahl und Vergleich (Point-Picker)

| # | Datei | Rolle |
|---|---|---|
| 31 | `01_Programm/thermo_state.py` | Zustandsobjekt |
| 32 | `01_Programm/state_store.py` | Zustandsspeicher |
| 33 | `01_Programm/state_picker.py` | Picker |
| 34 | `01_Programm/state_resolver.py` | Zustandsauflösung |
| 35 | `01_Programm/selection_manager.py` | Auswahlverwaltung |
| 36 | `01_Programm/comparison_view.py` | Zustandsvergleich |

## 5. Optional zur Laufzeit

| # | Datei | Rolle |
|---|---|---|
| 37 | `01_Programm/Debug.py` | Diagnosewerkzeug, von `main.py` und `gui.py` importiert |

## 6. Nicht für die Ausführung nötig – Doku nach Lesepriorität

| # | Datei | Rolle |
|---|---|---|
| 38 | `04_Dokumentation/README_Projektstruktur.md` | Projektstruktur |
| 39 | `04_Dokumentation/README_Implementierung.md` | Implementierungsüberblick |
| 40 | `04_Dokumentation/Doku/README_START.md` | Einstieg |
| 41 | `04_Dokumentation/Doku/README_GESAMT.md` | Gesamtübersicht |
| 42 | `04_Dokumentation/Doku/STUDIENARBEIT_ThermoSurface.md` | Studienarbeit |
| 43 | `04_Dokumentation/Doku/return_of_functions.md` | Rückgabewerte der Funktionen |
| 44 | `Implementierung_1_fvT_tangente_APT1-8.md` | Implementierung 1 |
| 45 | `Implementierung_2_QOL.md` | Implementierung 2 |
| 46 | `Implementierung_3_PointPicker_3D_Zustandsflaeche.md` | Implementierung 3 |
| 47 | `04_Dokumentation/ThermoSurface_Implementierungsplan_final.pdf`, `04_Dokumentation/Plan(1).pdf` | Pläne |
| 48 | `04_Dokumentation/Doku/README_AP2.md` … `README_AP20.md`, `README_AP5-AP10_Gesamtuebersicht.md` | Arbeitspakete |
| 49 | `04_Dokumentation/Doku/README_TS-VP01.md` … `README_TS-VP07.md`, `README_TS-VP_GESAMT.md` | Viewport-Arbeitspakete |
| 50 | `04_Dokumentation/Doku/README_GPT_FVT.md`, `README_SPINODALE_ERWEITERUNG.md` | Erweiterungen |
| 51 | `04_Dokumentation/Doku/analyse_zustandsflaeche.md`, `ThermoSurface_Projektanalyse_und_Arbeitspakete_AP1-AP10.md`, `ThermoSurface_Entwicklungsplan_AP11-AP18.md`, `ThermoSurface_Spinodalen_Gesamtanalyse_Tests_Ergebnisse_2026-10-19.md` | Analysen und Pläne |
| 52 | `04_Dokumentation/Doku/TESTPROTOKOLL_AP5-AP10.md`, `TS_A01_Ausgangszustand.md`, `TS_A_B_Abnahmeprotokoll.md`, `manifest_sha256.json` | Protokolle und Prüfsummen |

## 7. Tests und Archiv

| # | Pfad | Rolle |
|---|---|---|
| 53 | `03_Tests_Implementierung/conftest.py` | Gemeinsame pytest-Fixtures |
| 54 | `03_Tests_Implementierung/test_integration_mocked_main.py`, `test_gui_multiview.py`, `test_zustandsflaeche.py`, `test_state_comparison*.py` | Kerntests |
| 55 | `03_Tests_Implementierung/test_*.py` (übrige) | Weitere Tests |
| 56 | `03_Tests_Implementierung/validate_*.py`, `ts_gpt_fvt_validation.py`, `compare_pvt.py`, `snapshot_pvt.py`, `scripts/` | Validierungs- und Hilfsskripte |
| 57 | `03_Tests_Implementierung/regression/`, `results/`, `figures/` | Referenz- und Ausgabedaten |
| 58 | `02_Tests_Archiv/` | Archivierte Diagnosen, nicht eingebunden |
