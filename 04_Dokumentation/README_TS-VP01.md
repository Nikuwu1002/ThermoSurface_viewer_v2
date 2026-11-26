# TS-VP01 – Bestandsaufnahme, gemeinsame Verträge und Testbasis

**Stand:** 16.11.2026 · **Bezug:** `ThermoSurface_Technische_Spezifikation_TS-VP01-VP07.md` (Abschn. 3, TS-VP01)
**Ergebnis:** Matrix der sieben Diagrammtypen, schriftliche Verträge (Isolinien, Viewport, Cache), Baseline, dokumentierte Vorab-Defizite.
Die Spezifikation hatte `diagram_coords.py`, `diagram_registry.py` sowie `coolprop_heos_pht/pst/hsp/uvt.py` als „nicht vorliegend“ geführt. Alle waren im aktuellen Projektstand vorhanden und wurden vor jedem Eingriff gelesen; die Matrix unten stützt sich auf den tatsächlichen Code.

## 1. Datenfluss je Diagrammtyp (Ausgangszustand vor TS-VP02…06)

Rechenweg: `main.berechne_auswahl(stoff, raster, visu_raster, diagram_type)` → `BERECHNUNGEN[id]` → HEOS-Modul → Ergebnis-Dictionary → `DiagramRegistry.render` → `visualization.render_*` → Actor-Dictionary → `DiagramViewport`.

| ID | HEOS-Modul / Funktion | Native Spalten (Datenarrays) | Darstellung x / y / z (Transformation, `diagram_coords.py`) | Isolinien vorher | Renderer | Actor-Schlüssel (unverändert) |
|---|---|---|---|---|---|---|
| `pvt` | `coolprop_heos_pvt` via `main.berechne_daten` | (T [K], v [m³/mol], p [Pa]) | `log10(v/vc)` / `T/Tc` / `pressure_to_z(p)` (`log`, `asinh`, `symlog`, `linear`) | T_start/T_stop/`n_isotherms`, v_start/v_stop/v_step, p_start/p_stop/p_step (lineare Raster) | `render_pvt` | `surface`, `surface_single_phase`(leer), `metastable_*`, `unstable_surface`, `replacement_surface`, `zweiphasengebiet`, `isothermen`, `isochoren`, `isobaren`, `binodale`, `critical_point`, `spinodale`, `pressure_ticks` |
| `pht` | `coolprop_heos_pht.berechne_pht_daten` | (T, h [J/kg], p) | `(h-hc)/1e6` / `T/Tc` / `log10(p/pc)` | **keine** (Funktionen `calc_isotherme/isobare/isochore` existierten, wurden nicht aufgerufen) | `render_pht` | Basis + leere Isolinienlisten |
| `pst` | `coolprop_heos_pst.berechne_pst_daten` | (T, s [J/(kg K)], p) | `(s-sc)/1e3` / `T/Tc` / `log10(p/pc)` | fest: 7 Isothermen (T aus 7 Rasterindizes), 6 Isobaren (geometrisch zwischen p_min und p_max), 5 Isochoren (0.55…3 vc) | `render_pst` | Basis + Isolinien |
| `hsp` | `coolprop_heos_hsp.berechne_hsp_daten` (pht + pst, PS-Flash-Umrechnung) | (h, s, p) | `(h-hc)/1e6` / `(s-sc)/1e3` / `log10(p/pc)` | wie pst | `render_hsp` | wie pst |
| `uvt` | `coolprop_heos_uvt.berechne_uvt_daten` | (T, v, u [J/kg]) | `log10(v/vc)` / `T/Tc` / `(u-uc)/1e6` | fest: 6+4 Isothermen aus Netzzeilen, 4 Isochoren (0.55…5 vc), 3 Isobaren (0.2…1.2 pc) | `render_uvt` | wie pst |
| `gpt` | `coolprop_heos_gpt.berechne_gpt_daten` | (T, p, g [J/kg]) | `pressure_to_z(p)` (+ Druckfenster bei `log`) / `T/Tc` / `(g-gc)/1e6` | fest: 6+4 Isothermen aus Netzzeilen, 3 Isobaren, 4 Isochoren | `render_gpt` → `_render_free_energy` | wie pst + homogene Ebenen, Spinodalen |
| `fvt` | `coolprop_heos_fvt.berechne_fvt_daten` | (T, v, f [J/kg]) | `log10(v/vc)` / `T/Tc` / `(f-fc)/1e6` | fest wie uvt | `render_fvt` → `_render_free_energy` | wie gpt |

Koordinatenvertrag (Spezifikation 1.3) bestätigt: pvt = (v, T, p), pht = (h, T, p), pst = (s, T, p), hsp = (h, s, p), uvt = (v, T, u), gpt = (p, T, g), fvt = (v, T, f). Alle Vorwärtstransformationen sind je Achse streng monoton steigend (Voraussetzung für das Clipping an einem Quader). Maschinenlesbar abgelegt in `diagram_axes.DIAGRAM_AXES`; Rundlauftest in `test_ts_vp04_05_06_viewport.py::AxisMatrix`.

Layer: `LAYER_LABELS` in `diagram_registry.py`; pvt alle zwölf, pht/pst/hsp/uvt `BASE_LAYERS (+ ISOLINE_LAYERS)`, gpt/fvt zusätzlich homogene Ebenen/Spinodalen. **Änderung in TS-VP03:** `pht` führt jetzt ebenfalls die Isolinien-Layer.

Cache (`diagram_cache.DiagramKey`): Stoff, Diagramm-ID, SHA-1 der tatsächlich verwendeten Parameter, Modell, CoolProp-Version. pvt: `raster`, `visu_raster`, `diagram_settings`; Zusatzdiagramme: `Q_points` (+ `diagram_settings`). Kamera, Layer, Fenstergröße gehören nicht in den Schlüssel.

## 2. Festgeschriebene Verträge (spätere Pakete ändern sie nicht eigenmächtig)

### 2.1 Isolinienkonfiguration (TS-VP02/03)

```yaml
diagram_settings:
  <id>:                      # pvt | pht | pst | hsp | uvt | gpt | fvt
    isolines:
      isothermen: {start: auto|<K>,      stop: auto|<K>,      factor: 3.14}
      isochoren:  {start: auto|<m3/mol>, stop: auto|<m3/mol>, factor: 3.14}
      isobaren:   {start: auto|<Pa>,     stop: auto|<Pa>,     factor: 3.14}   # optional: max_levels (Standard 200)
```

* Intern ausschließlich SI; `auto` = Bereich des Diagramms (Details `README_TS-VP03.md`).
* Wirksame Vorgabe: GUI-/Aufgabenvorgabe (`DiagramTask.isolines`) vor `config.yaml`; fehlt beides, gilt für p-v-T das bisherige Raster, für alle anderen der Standard (Faktor 3.14, automatischer Bereich).
* Ergebnis-Dictionary: unveränderte Schlüssel `isothermen/isochoren/isobaren` (Listen von `(n,3)`-Arrays) plus `isoline_report` (Levels, Anzahl, verworfene Levels, Abtastpunkte je Kurve).

### 2.2 Viewport (TS-VP04…06)

* `viewport_bounds.ViewportState` (`diagram_id`, `enabled`, `bounds` {Größe: (min_si, max_si)}, `fields` letzte gültige Eingaben, `camera_state`, `user_defined`); je Diagramminstanz, in der GUI je Diagrammtyp gespeichert.
* Grenzen intern SI; GUI-Anzeigeeinheiten: T [K], p [MPa], v [m³/mol], h/u/g/f [kJ/kg], s [kJ/(kg K)].
* Änderungen betreffen nur Darstellung: kein `berechne_*`, kein Cache-Zugriff, kein Eingriff in das Ergebnis-Dictionary.

### 2.3 Cache (TS-VP03/07)

* Neu im Schlüssel: kanonisierte Isolinienparameter (`isoline_levels.config_signature`) – nur Berechnungsparameter.
* Nicht im Schlüssel: Viewportgrenzen, Kamera, Layer, Gitter.

## 3. Baseline

* Vor dem ersten Eingriff wurden alle `*.py` und `config.yaml` gesichert.
* Numerische Referenz- und Strukturtests (bestehend, verbindlich): `test_fluid_diagram_matrix.py` (10 Stoffe × Diagramme), `test_gpt_fvt.py`, `test_state_comparison.py`, `test_zustandsflaeche.py`, `test_ap13_ap18.py`, `test_ap19_ap20.py`, `test_spinodal_extension.py`, `test_integration_mocked_main.py`, `test_gui_multiview.py`, `test_state_comparison_gui.py`. Ergebnisse nach TS-VP07 siehe `README_TS-VP07.md`.
* **Vorbestehender Fehlschlag (nicht Folge dieser Arbeit):** `test_zustandsflaeche.py::P8EndeZuEnde::test_regression_against_p0_baseline_is_bitwise_equal` meldet `single_phase_surface n_points 9733 != 9614` (Referenz `regression/p0_baseline.json`). Identisches Ergebnis mit den gesicherten Originaldateien – die gespeicherte Baseline stammt offenbar von einer anderen CoolProp-/Rasterversion (installiert: CoolProp 7.2.0).

## 4. Dokumentierte geometrische/funktionale Defizite des Ausgangszustands (nicht Folge der Erweiterung)

1. `render_pht()` zeichnete keine Isolinien; die Registry bot dafür keine Layer (behoben in TS-VP03).
2. uvt/gpt/fvt entnahmen Isothermen aus **Netzzeilen** (nur Gitter-T möglich) und nutzten fest codierte Faktoren (0.55, 1, 2, 5 · vc; 0.2, 0.5, 1.2 · pc) – nicht als geometrische Progression steuerbar (behoben in TS-VP03).
3. Achsenbeschriftung zeigte transformierte Größen (`log10(v_r)`, `T_r`, `h - h_c [MJ/kg]` …) statt physikalischer Werte; pvt hatte als einziges Diagramm eigene Druck-Ticks (behoben in TS-VP06).
4. Achsen wurden diagrammspezifisch an verschiedenen Stellen gesetzt (`show_bounds` in jedem Renderer, `update_pvt_axes` nur bei pvt).
5. g-p-T mit `pressure_axis: log` blendet homogene Knoten außerhalb des berechneten Druckfensters aus (dokumentiert, Daten bleiben unverändert; Statuszeilenhinweis vorhanden).
6. Hover-Text nur für pvt (`spec.hover`); Punktpicker für alle sieben Diagramme.
7. Die gespeicherte `p0_baseline.json` passt nicht zur installierten CoolProp-Version (siehe 3).
8. Bei `factor = 3.14` und engen Temperaturfenstern (Verhältnis `T_stop/T_start` ≈ 1.6 in allen Zusatzdiagrammen) ergibt sich nur eine Isotherme – mathematische Folge der gewünschten Progression, kein Fehler (siehe `README_TS-VP03.md`).

## 5. Abnahme TS-VP01

Matrix konsistent (Abschn. 1), Verträge schriftlich (Abschn. 2), Baseline gesichert (Abschn. 3), Defizite dokumentiert (Abschn. 4); die Matrix ist über `test_ts_vp04_05_06_viewport.py::AxisMatrix` und `test_ts_vp03_isolines.py` automatisiert geprüft.
