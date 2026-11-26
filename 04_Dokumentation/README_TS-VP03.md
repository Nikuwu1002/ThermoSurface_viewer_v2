# TS-VP03 – Isolinien in alle sieben Berechnungspfade integriert

**Module:** `coolprop_heos_pvt.py`, `coolprop_heos_pht.py`, `coolprop_heos_pst.py`, `coolprop_heos_hsp.py`, `coolprop_heos_uvt.py`, `coolprop_heos_gpt.py`, `coolprop_heos_fvt.py`, `main.py`, `calculation_manager.py`, `gui.py`, `visualization.py`, `diagram_registry.py`, `config.yaml`
**Tests:** `test_ts_vp03_isolines.py` (L04, L05, C02, Stoffmatrix) · Verwendung des Generators: `README_TS-VP02.md`

## 1. Was sich ändert

* **Alle sieben Diagramme** setzen die konstanten T-, v-, p-Werte aus `isoline_levels.resolve_levels`; die hartcodierten Arrays in g-p-T, f-v-T (und u-v-T, p-s-T) sind entfernt. Die Kurvenberechnung selbst (Phasenregeln, Sättigungsknoten, NaN-Trennung, `HmassP`/`PSmass`-Flashs) bleibt unverändert.
* **p-v-T:** `multi_isochore(..., v_values=…)` und `multi_isobare(..., p_values=…)` nehmen explizite Levelarrays; `multi_isotherm(..., T_values=…)` konnte das bereits. Der separate Oberflächenrasterpfad `multi_isotherm(..., T_values=T_all)` bleibt unberührt.
* **p-h-T:** `berechne_pht_daten(..., include_isolines=…, isolines=…)` liefert erstmals Isolinien (`calc_isotherme/isobare/isochore`), `render_pht()` zeichnet sie (NaN-getrennte Abschnitte), `pht` besitzt jetzt die Layer *Isothermen/Isobaren/Isochoren*. Die Voreinstellung der Funktion ist `include_isolines=False` (unverändertes Rückgabeformat bei direkten Aufrufen); `config.yaml` schaltet sie ein.
* **u-v-T/f-v-T/g-p-T:** Isothermen für beliebige T (nicht mehr nur Netzzeilen): `isotherm_over_volume` (Flüssig-Ast → lineare Q-Strecke → Dampf-Ast bzw. ein Ast oberhalb Tc) bzw. `_isotherme` in g-p-T (Dampfast → gemeinsamer Sättigungsknoten → Flüssigast).
* **Bestehendes bleibt:** `include_isolines`-Optionen, Layer-Checkboxen (`isothermen`, `isochoren`, `isobaren`), Stoffauswahl, Ergebnis-Schlüssel und Rückgabeformat (Listen von `(n,3)`-Arrays; pvt: Arrays bzw. Liste, wie bisher).

## 2. Bereiche („auto“) und Gültigkeitsfenster

`start`/`stop` = `auto` setzt die Grenzen des Diagrammbereichs; Levels außerhalb des Fensters (dort existieren keine berechneten Flächendaten) werden **nicht gezeichnet** und im Bericht als `dropped_outside_domain` ausgewiesen. Es wird nicht nachgerechnet.

| Diagramm | Isothermen T | Isochoren v | Isobaren p | Abtastung entlang der Linie (unverändert getrennt) |
|---|---|---|---|---|
| pvt | `isotherm_visu`: T_start…T_stop | `isochore_visu`: v_start…v_stop | `isobare_visu`: p_start…p_stop | Isotherme: v-Raster von `isotherm_visu`; Isochore/Isobare: T-Raster (T_step) |
| pht, pst, hsp | `T_sub[0]`…`T_high[-1]` (Standard 0.75…1.2 Tc) | 0.5…30 vc | `p_min`…`p_max` der Fläche (Standard bis 2 pc) | 41+41 p-Punkte + n_quality; Isobare 81 + n_quality; Isochore 121 |
| uvt, fvt | wie oben | `v_min`…`v_max` der Fläche | 0.1…10 pc | Isotherme 2·n_v−1 + Q-Strecke; Isobare 61 + n_quality; Isochore = T-Raster der Fläche |
| gpt | wie oben | 0.5…30 vc | `p_min`…`p_max` (Druckfenster der Fläche) | Isotherme 2·n_p−1; Isobare 81; Isochore = T-Raster |

Isochoren/Isobaren, die den Flächenbereich verlassen (z. B. `p` außerhalb `[p_min, p_max]`), erhalten an diesen Punkten NaN; die Renderer trennen dort die Linie (keine Brücke). Levels ohne mindestens zwei gültige Punkte erscheinen als `empty_levels`, HEOS-Fehler als `failed_levels`.

## 3. Konfiguration und Bericht

```yaml
diagram_settings:
  gpt:
    include_isolines: true
    isolines:
      isothermen: {start: auto, stop: auto, factor: 3.14}
      isochoren:  {start: auto, stop: auto, factor: 3.14}
      isobaren:   {start: 2.0e6, stop: 2.0e7, factor: 2.0}     # Pa
```

`daten["isoline_report"][familie]`: `levels`, `n_levels`, `start`, `stop`, `factor`, `domain`, `dropped_outside_domain`, `computed_levels`, `failed_levels`, `empty_levels`, `samples_per_curve`, `n_curves` – Levelanzahl und Abtastpunktzahl sind getrennt protokolliert. Die GUI-Statuszeile zeigt `Isolinien T/v/p: n/n/n`.

**Wichtiger Hinweis zum Standardfaktor 3.14:** Mit `T_stop/T_start ≈ 1.6` (alle Zusatzdiagramme) und 580…1200 K (p-v-T/Wasser) ergibt der Faktor nur **eine** Isotherme; Druck- und Volumenfamilien über mehrere Dekaden liefern 3…10 Linien. Das ist die mathematische Folge der gewählten Progression (`x[n+1] = f·x[n]`). Für dichtere Scharen den Faktor in der GUI oder in `config.yaml` verkleinern (z. B. `1.1`).

## 4. GUI (Gruppe „Isolinien (geometrische Progression)“)

* Je Familie `Start`, `Stop` (0 = auto) und `Faktor f`; Einheiten: T [K], v [m³/mol], p [**MPa**] (Umrechnung nach SI in `ThermoSurfaceWindow.isoline_config_si`, Test vorhanden).
* Gilt für das gewählte Diagramm (Einzelansicht) bzw. das aktive/erste gewählte Diagramm (Mehrfachansicht); Eingaben bleiben je Diagramm erhalten.
* Ungültige Eingaben (`stop < start`, `f ≤ 1`, zu viele Levels) werden vor dem Berechnungsstart gemeldet.
* Die alten Felder des p-v-T-Rasters werden nicht zweckentfremdet: „Anzahl Isothermen“ ist deaktiviert (Tooltip: ersetzt durch die Progression); Temperatur-/Volumen-/Druckbereiche der `*_visu`-Gruppen sind Automatikbereiche, die Schrittweiten `T_step`, `v_step` bleiben Abtastung entlang der Linien.

## 5. Cache (C02)

`cache_params(raster, visu_raster, isolines=None)` ergänzt `isolines = config_signature(wirksame Vorgabe)`. Gleiche Parameter (auch GUI-Vorgabe = config.yaml-Vorgabe) → gleicher Schlüssel (Cache-Treffer); geänderter Faktor/Start/Stop → neuer Schlüssel und Neuberechnung nur dieses Diagramms. Viewportgrenzen sind nicht Teil des Schlüssels.

## 6. Nachweis

* **L04** (`ConstantAlongCurve`): für Wasser und CO₂, alle sieben Diagramme, alle Familien: T (|ΔT| ≤ 2·10⁻³ K), v und p (relativ ≤ 10⁻⁵) konstant entlang jeder Kurve – direkt aus den Daten, sonst per HEOS rekonstruiert (`HmassP`, `PSmass`, `PT`, `DmolarT`). Ausnahme dokumentiert: g-p-T-Isochorenpunkte im Nassdampfgebiet liegen auf der Sättigungslinie (g_l = g_g) und sind aus (T, p) nicht rekonstruierbar.
* **L05**: Isolinien-Actors in allen sieben Renderern, Layer in allen Registry-Einträgen (`RendererActors`).
* Gleiche Levels in allen Diagrammen: dieselbe Formel `x_n = start·f^n` im jeweiligen Gültigkeitsfenster (`LevelsAcrossDiagrams`).
* Flächen-/Binodalen-Arrays sind bei unterschiedlichen Isolinienparametern bitgleich (`SeparationOfLevelsAndGrids`).
* Stoffmatrix Wasser, CO₂, R134a, Argon × 7 Diagramme mit Standardisolinien (`FluidMatrix`).
