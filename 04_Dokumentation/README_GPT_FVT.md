# g-p-T und f-v-T (3D) – Implementierung, Thermodynamik, Validierung

Stand: 11.11.2026 · Umgebung der Prüfung: Windows, Python 3.13.9, CoolProp 7.2.0, pyvista 0.47.3, PySide6 6.12.0, pyvistaqt, h5py 3.16.0 (optional).
Aufbauend auf dem Stand nach TS-A/TS-B (Registry, Mehrfachansicht, Zustandsvergleich).

## 1. Architektur und verwendete Schnittstellen

**Analyseergebnis.** Die Anwendung ist registry-getrieben: `diagram_registry.DIAGRAMS` (Dropdown, Checkboxen, Layer, Renderer-Name, Cache-Schlüssel) → `main.BERECHNUNGEN` (Berechnung je ID) → `coolprop_heos_<typ>.berechne_<typ>_daten(stoff, …)` → `visualization.render_<typ>(plotter, daten)` → `DiagramViewport` (Layer, Raster, Picking, Marker) → `state_store`-Adapter + `state_resolver` (Punktauswahl/Zustandsvergleich). `gui.py`, `diagram_viewport.py`, `diagram_workspace.py`, `calculation_manager.py` und `diagram_cache.py` enthalten **keine** diagrammspezifische Logik und blieben unverändert.

Ein neues Diagramm benötigt daher genau diese Bausteine (alle umgesetzt):

| Baustein | g-p-T | f-v-T |
|---|---|---|
| Berechnungsmodul (neu) | `coolprop_heos_gpt.py` | `coolprop_heos_fvt.py` |
| Registry-Eintrag, Layer, Info-Text | `diagram_registry.py` | `diagram_registry.py` |
| Berechnungsaufruf + Cache-Parameter | `main.BERECHNUNGEN["gpt"]`, `DIAGRAMS.attach` | `main.BERECHNUNGEN["fvt"]`, `DIAGRAMS.attach` |
| Darstellungskoordinaten (PyVista-frei) | `diagram_coords.to_xyz_gpt/from_xyz_gpt` | `to_xyz_fvt/from_xyz_fvt` |
| Renderer | `visualization.render_gpt` | `visualization.render_fvt` |
| Picker/Zustandsrekonstruktion | Adapter `_GptAdapter`, Flash `PT` | Adapter `_FvtAdapter`, Flash `DmolarT` |
| Konfiguration | `config.yaml: diagram_settings.gpt` | `diagram_settings.fvt` |

**Wiederverwendet (keine zweite Logik):** Temperaturraster, Druckfenster und `MAX_POINTS`-Prüfung aus `coolprop_heos_pht.generate_pht_grid` (kritisch verdichtetes T-Raster), `pt_update_robust` (PT-Flash am kritischen Rand), die validierten p-v-T-Verfahren `homogeneous_pvt.solve_spinodal_branch` und `build_homogeneous_surfaces` (Spinodalen, metastabile/instabile Knoten, Lücken statt Interpolation), die Mesh-Bausteine `_new_surface_mesh`, `_new_phase_mesh`, `_new_curve`, `grid_triangles`, `drop_degenerate` sowie die komplette Viewport-/Workspace-/Picker-/Vergleichsinfrastruktur.

Neue gemeinsame Funktionen in `coolprop_heos_common.py` (rein additiv): `gibbs_helmholtz_tv`, `critical_gf`, `homogeneous_gf_layers`.

## 2. Thermodynamische Grundlagen

Massenspezifische Größen (J/kg), direkt aus CoolProp: `gibbsmass` (g = h − T·s) und `helmholtzmass` (f = u − T·s). Beide Identitäten und g = f + p·v sind in `results/gpt_fvt_validation.json` und in den Tests nachgewiesen (Abweichung ≤ 1,4·10⁻⁸ relativ, s. Abschnitt 7).

### g-p-T – Achsen: x = log₁₀(p/pc), y = T/Tc, z = (g − gc)/10⁶ [MJ/kg]
* **Natürliche Variablen (T,p).** (∂g/∂p)_T = v > 0, (∂g/∂T)_p = −s. Die Gleichgewichtsfläche ist stetig mit **Knick an der Sättigungslinie** (dg/dp springt von v_g auf v_l) und ist die untere Hülle aus Flüssigkeits- und Dampfast.
* **Das Zweiphasengebiet entartet zu einer Kurve:** Gleichgewicht verlangt g_l = g_g bei p = p_sat(T). Es gibt **keine Zweiphasenfläche** und keinen Dampfgehalt Q aus (T,p,g). Siede- und Taulinie liegen aufeinander (`binodal_liquid`, `binodal_vapor`: gleiches (T,p), g nur im Rundungsbereich verschieden; die Gleichheit wird geprüft, nicht erzwungen). Die Ebene „Zweiphasengebiet“ ist für g-p-T bewusst nicht vorhanden; „Siede- und Taulinie“ zeigt die Sättigungslinie.
* **Metastabile/instabile Zustände** bilden mit dem stabilen Ast den **Schwalbenschwanz** (in (T,p) mehrwertig). Sie liegen energetisch **über** der Gleichgewichtsfläche (Test: g_homogen ≥ g_stabil bei gleichem (T,p)). Die Flächen werden über (T,v) parametrisiert.

### f-v-T – Achsen: x = log₁₀(v/vc), y = T/Tc, z = (f − fc)/10⁶ [MJ/kg]
* **Natürliche Variablen (T,v).** p = −(∂f/∂v)_T, s = −(∂f/∂T)_v. f(T,v) ist auch im metastabilen und instabilen Bereich einwertig.
* **Zweiphasengebiet = Regelfläche (gemeinsame Tangente / Maxwell-Konstruktion):** f = (1−Q)·f_l + Q·f_g, v = (1−Q)·v_l + Q·v_g, df/dv = −p_sat/M. Anders als in g-p-T gilt f_l ≠ f_g. g = f + p_sat·v ist im ganzen Zweiphasengebiet konstant (= g_sat).
* **Homogene Fortsetzung** liegt über der Tangentenfläche (f_homogen ≥ f_Tangente).

### Volumeneinheit (bewusste Abweichung von „konsequent massenspezifisch“)
Der Auftrag fordert massenspezifische Größen; die Codebasis führt das Volumen in **allen** Diagrammen mit v-Achse (p-v-T, u-v-T, Zustandsdatenbank) als **molares Volumen v [m³/mol]**. Um die Achse mit p-v-T/u-v-T und dem Zustandsvergleich (`ThermoState.v_molar`) kompatibel zu halten, steht v in f-v-T ebenfalls molar in den Daten (wie in u-v-T), f und g sind massenspezifisch. Die Darstellungsachse log₁₀(v/vc) ist davon unabhängig (Verhältnis). Das massenspezifische Volumen ist `v/M` (`meta["molar_mass"]` in kg/mol); alle Konsistenzprüfungen (g = f + p·v) rechnen mit v_mass = v/M.

### Referenzzustand
Energien, Enthalpien, Entropien stehen im CoolProp-Standardreferenzzustand des Fluids (unverändert). Physikalisch aussagekräftig sind Differenzen; deshalb die Darstellung relativ zum kritischen Zustand (g − gc, f − fc).

## 3. Implementierte Funktionen

| Funktion | g-p-T | f-v-T |
|---|---|---|
| Zustandsfläche (Flüssigkeit, Dampf, Hochtemperatur; randkonform, identische Randknoten) | ✔ | ✔ |
| Zweiphasengebiet | Sättigungslinie (keine Fläche) | ✔ Regelfläche |
| Siede-/Taulinie, kritischer Punkt | ✔ | ✔ |
| Spinodalen (Flüssigkeit/Dampf) | ✔ | ✔ |
| Metastabile Flüssigkeit/Dampf, instabile Fortsetzung (homogene HEOS) | ✔ | ✔ |
| Isothermen, Isobaren, Isochoren | ✔ | ✔ |
| Rasteranzeige (Gitter ein/aus) | ✔ | ✔ |
| Einzeln ein-/ausblendbare Ebenen je Viewport | ✔ | ✔ |
| Mehrfachdiagramm (mit allen 7 Typen, Maximieren, Cache) | ✔ | ✔ |
| Point Picker + Zustandsvergleich (Koordinaten, g, f, h, s, u, v, p, Einheiten, Status) | ✔ | ✔ |
| Stoffauswahl mit automatischer Bereichsbestimmung (T: 0,75 Tc … 1,2 Tc, p: ≈ 0,02 pc … 2 pc (bzw. 0,25·p_sat(T_min), falls kleiner), v: aus den Binodalen) | ✔ | ✔ |
| Rasterauflösung (`n_T_sub`, `n_T_super`, `n_p` bzw. `n_v`, `n_quality`=`Q_points`) | ✔ | ✔ |

Standardmäßig sind in der GUI nur Fläche, Siede-/Taulinie, (Zweiphasenfläche bei f-v-T) und kritischer Punkt an; Spinodalen, homogene Ebenen und Isolinien sind **aus**, damit die erste Ansicht die reine Gleichgewichtsfläche zeigt (keine Vermischung).

**Bewusst nicht vorhanden:** „Ersatzfläche“ (pvt-spezifisches Modell zwischen den Spinodalpunkten, keine EOS-Auswertung) und die Dreiecksfächer zum kritischen Punkt der homogenen Flächen (pvt-AP20.6): sie sind für g-p-T/f-v-T weder erforderlich noch thermodynamisch begründet. Die Option ist in der Layer-Liste dieser Diagramme nicht enthalten.

## 4. Geänderte Dateien und relevante Funktionen

**Neu:** `coolprop_heos_gpt.py` (`berechne_gpt_daten`, `multi_binodale_gpt`, `validate_gpt_data`, `g_eos`, `critical_point_gpt`), `coolprop_heos_fvt.py` (`berechne_fvt_daten`, `validate_fvt_data`, `f_eos`, `critical_point_fvt`), `test_gpt_fvt.py`, `scripts/ts_gpt_fvt_validation.py`, `results/gpt_fvt_validation.json`, `Doku/README_GPT_FVT.md`.

**Geändert (nur Ergänzungen, bestehende Funktionen unverändert):**
* `coolprop_heos_common.py` – `gibbs_helmholtz_tv`, `critical_gf`, `homogeneous_gf_layers`, `HOMOGENEOUS_LAYER_KEYS`.
* `diagram_coords.py` – `to_xyz_gpt/fvt`, `from_xyz_gpt/fvt` (g-p-T mit Druckachse log/asinh/symlog/linear und Druckfenster).
* `diagram_registry.py` – Einträge `gpt`, `fvt`, `GPT_LAYERS`, `FVT_LAYERS`, `_gpt_info`.
* `main.py` – `_berechne_zusatzdiagramm` akzeptiert Modulnamen (Import erst beim Aufruf, damit die isolierten Ladetests von `main` mit ersetzten HEOS-Modulen unverändert bleiben), `BERECHNUNGEN`/`attach` für `gpt`, `fvt`.
* `visualization.py` – `render_gpt`, `render_fvt`, `_render_free_energy`, `_masked_grid_mesh`; `render_*` der alten Diagramme unverändert.
* `state_store.py` – `SourceSurface.aux`/`node_aux`, `_GptAdapter`, `_FvtAdapter`, `_FreeEnergyMixin`.
* `state_resolver.py` – Flash-Plan `gpt: PT`, `fvt: DmolarT`; `_update("PT")`; Rückfall auf die Phase der Quellfläche auf der Sättigungslinie; Plausibilität auch für g und f; Volumen der homogenen g-p-T-Knoten aus `aux` (`_aux_volume`).
* `state_picker.py` – Hinweis `vertex_aux_v`. `thermo_state.py` – Einheiten `g`, `f` [J/kg]. `selection_manager.py` – Farben `gpt`, `fvt`.
* `config.yaml` – `diagram_settings.gpt` und `diagram_settings.fvt` (Rückwärtskompatibel: fehlen die Abschnitte, gelten die Funktionsvorgaben).
* Tests: `test_gui_multiview.py` (7 registrierte Typen, Zähler `len(DIAGRAMS)`), `test_state_comparison.py` (Zweiphasen-Test überspringt g-p-T, das keine Zweiphasenfläche hat).

## 5. Unterschiede zur bisherigen Diagrammberechnung

* **g-p-T hat keine Zweiphasenfläche** (Kurve), p-h-T/p-s-T/h-s-p/u-v-T besitzen eine.
* **Geometrische Transformation reicht nicht:** g und f werden aus HEOS-Zuständen (PT-Flash bzw. DmolarT) neu berechnet; sie sind keine Umskalierung der vorhandenen Flächen (g ≠ f ≠ h ≠ u; verschiedene Hüllen- und Knickstruktur).
* **Homogene Ebenen in zwei Diagrammen, die dieselben (T,v)-Knoten teilen:** f-v-T zeigt sie einwertig, g-p-T als Schwalbenschwanz (mehrwertig, daher Darstellung über (T,v)).
* **Picker g-p-T:** Flash (T,p); auf der Sättigungslinie liefert HEOS ohne Phase keinen Zustand – es wird der gesättigte Zustand der Seite der angeklickten Fläche rekonstruiert (Q = `None`, Hinweis im Zustand). Ein Klick knapp jenseits der Sättigungslinie ergibt den **Gleichgewichtszustand** (z. B. Dampf) mit Status `partially_valid` und Hinweis „Phasenrandnähe“.

## 6. Phasengrenzen, metastabile und instabile Zustände

* **Gleichgewichtsflächen:** stabile Phase wird bei T < Tc vorgegeben (`iphase_liquid/gas`), bei T ≥ Tc ohne Vorgabe (robuster PT-Flash). Randknoten sind exakt die Sättigungswerte (bitgleich, geprüft in `validate_*`).
* **Homogene Fortsetzung:** HEOS mit vorgegebener Phase ohne Gleichgewichts-Flash (`DmolarT`). Metastabile Äste von der Binodale bis zur Spinodale, instabiler Zweig zwischen den Spinodalen. Die Knotenauswahl stammt aus den validierten p-v-T-Verfahren (Stabilitätsklassifikation (dp/dv), Druckfenster, Lücken statt Interpolation).
* **Zusätzliches Kriterium (neu, thermodynamisch begründet):** bei gleichem (T,v) muss f_homogen ≥ f_Gleichgewicht (gemeinsame Tangente) gelten. Knoten darunter sind Artefakte der Mehrparameter-Zustandsgleichung im tiefen instabilen Bereich und werden **nicht dargestellt** (Wasser: 29 von 900 instabilen Knoten; Methan: 15; CO₂, N₂, Ethanol: 0). Werte bleiben unverändert, es wird nur klassifiziert (`homogeneous_info["surfaces"][…]["below_equilibrium_hull_nodes"]`).
* **Instabile Fläche mit Lücken:** wie in p-v-T enthält die Zeilenmitte eine Trennspalte (NaN); sie wird nicht überbrückt.
* **Druckachse in g-p-T:** `log` (Standard) kann p ≤ 0 und p → 0⁺ nicht abbilden. Homogene Knoten außerhalb des berechneten Druckfensters [p_min, p_max] werden deshalb nicht dargestellt (Daten unverändert; Statuszeile nennt die Anzahl). `pressure_axis: symlog|asinh|linear` oder `auto_pressure_axis: true` zeigt alle (Wasser: 634 von 2629 homogenen Knoten liegen außerhalb des Fensters, 623 davon p ≤ 0 – typisch für gestreckte Flüssigkeit).
* **Kritische Nähe:** das kritisch verdichtete Temperaturraster aus p-h-T wird wiederverwendet; der kritische Punkt ist exakter gemeinsamer Endpunkt aller Kurven (geprüft). Spinodalkurven enden mit dem exakten kritischen Zustand (kein Fächer).
* **Nicht dargestellt / nicht unterstützt:** Zustände jenseits des HEOS-Gültigkeitsbereichs führen zu `ValueError` bzw. im Picker zu Statuscodes `out_of_range`/`numerical_failure` (keine Scheinwerte).

## 7. Durchgeführte Tests und Ergebnisse (tatsächlich ausgeführt)

**Neuer Testsatz `test_gpt_fvt.py`:** 47 Tests (40 ohne GL, 7 mit nativem OpenGL-Fenster). Referenzen aus unabhängigen CoolProp-Aufrufen (`PropsSI`, eigene `AbstractState`), Fluide Water, CO₂, Nitrogen, Methane (Matrixtest: zehn GUI-Stoffe).

| Prüfung | Ergebnis |
|---|---|
| g = h − T·s, f = u − T·s, g = f + p·v (massenspezifisch) auf Stützknoten aller Teilflächen | bestanden |
| (∂g/∂p)_T = v, (∂g/∂T)_p = −s; p = −(∂f/∂v)_T, s = −(∂f/∂T)_v (Differenzenquotienten) | bestanden (rel. 10⁻⁵ / 10⁻⁴) |
| g-p-T: g_l = g_g auf der Sättigungslinie, Knick (Steigung v_g links, v_l rechts), keine Zweiphasenfläche | bestanden |
| g-p-T: Gleichgewichtsfläche = untere Hülle (g_metastabil ≥ g_stabil, >200 Knoten je Stoff) | bestanden |
| f-v-T: Regelfläche = gemeinsame Tangente (Steigung −p_sat/M, f gegen unabhängigen (T,ρ)-Zustand, g_sat konstant) | bestanden |
| f-v-T: homogene f ≥ Tangente; (d²f/dv²) > 0 metastabil, < 0 instabil | bestanden |
| Spinodalen: (dp/dρ)_T = 0 (unabhängig), g-p-T und f-v-T beschreiben denselben Zustand | bestanden |
| Homogene Knoten = direkte HEOS-Werte (T,v, Phase vorgegeben); Monotonie p(v) je Klasse | bestanden |
| Kritischer Punkt gemeinsamer Endpunkt aller Kurven; Randkontinuität bitgleich | bestanden |
| Absichtlich fehlerhafte Datensätze (NaN, Randabweichung, dg/dp ≤ 0, f-Rand, Binodale) werden abgelehnt | bestanden |
| Ungültige Eingaben (`n_p=1`, `n_v=1`, falsche Druckachse, v-Bereich, unzulässiges p_min, unbekannter Stoff, T < 0) | bestanden |
| Rasterauflösungen, Optionen (`include_homogeneous`, `include_isolines`) | bestanden |
| Registry (7 Typen, Reihenfolge, Layer, Standard-Ebenen), `berechne_auswahl`, Cache-Parameter, Koordinaten-Umkehrung (4 Druckachsen), Druckfenster | bestanden |
| Zehn GUI-Stoffe × g-p-T/f-v-T berechenbar, alle Stützflächen endlich | bestanden |
| Picker-Zustände: exakte Knoten (Einheit J/kg → J/mol), Innenpunkte gegen direkte HEOS-Aufrufe, Sättigungsrand (Q = None, Hinweis), Punkt jenseits der Sättigung, homogene Punkte (Volumen aus Knoten), Zweiphasen-Q, kritischer Knoten | bestanden |
| Zustandsvergleich: g-, f-Achsen, gemischte Diagrammtypen ohne Neuberechnung | bestanden |
| GUI nativ: 7 Typen in Auswahl/Checkboxen, Layer ein/aus je Viewport, Diagrammwechsel ohne Neustart (Cache-Treffer), Mehrfachansicht mit allen 7 (nur neue berechnet), echter Mausklick-Picker auf g-p-T/f-v-T und auf homogenen Ebenen, Marker, Doppelklick maximiert ohne Auswahl, Picker aus bei deaktiviertem Vergleichsmodus, Renderer mit allen Druckachsen und ohne optionale Ebenen | bestanden |

**Validierungslauf** (`python scripts/ts_gpt_fvt_validation.py`, 10 Stoffe, je 300 Zufallsknoten je Teilfläche, über 20 000 Prüfknoten): größte relative Abweichung g-p-T **1,4·10⁻⁸** (Ethanol; Genauigkeit des HEOS-PT-Dichteflashs), f-v-T **≤ 5,5·10⁻¹³**; Laufzeit je Stoff und Diagramm 0,4–1,1 s (Standardraster, inkl. homogener Ebenen und Isolinien).

**Regression bestehender Funktionen:**
* `QT_QPA_PLATFORM=offscreen python -m pytest` (gesamtes Projekt): 4 Fehlschläge, **dieselben wie vor der Änderung**: `P8EndeZuEnde::test_regression_against_p0_baseline_is_bitwise_equal` (9733 ≠ 9614; bereits vor TS-A vorhanden, nicht behoben) sowie drei Tests (`P6…symlog_and_drops_them_on_log`, `P7…substitute_vertices…`, `P7…pick_hover_text_offscreen`), die nur im Gesamtlauf scheitern, weil `test_ap13_ap18` pyvista durch einen Stub ersetzt (einzeln ausgeführt bestehen sie). Alle übrigen Tests grün.
* Nativ (ohne `offscreen`): `test_gui_multiview.py`, `test_state_comparison.py`, `test_state_comparison_gui.py`, `test_gpt_fvt.py`: **111 bestanden**.
* Die fünf bisherigen Diagramme liefern für Water und CO₂ unveränderte Ergebnis-Schlüssel und Array-Formen gegenüber `results/ts_a01_referenz.json` (einzige Abweichung: `homogeneous_mesh_info` wird bei p-v-T erst vom Renderer ergänzt, so schon in der Baseline vermerkt).

**Ausschließlich statisch geprüft bzw. nicht ausgeführt:** händisches Bedienen mit der Maus durch einen Anwender (die Klicks wurden mit `QTest` synthetisch ausgelöst), andere Betriebssysteme/Grafikkarten, Stoffe außerhalb der zehn getesteten, die optische Beurteilung der Pressure-Achsen `asinh`/`linear` (nur Rendern ohne Fehler und Koordinaten-Umkehrung getestet).

## 8. Bekannte Einschränkungen und offene Probleme

1. **v molar statt massenspezifisch** in f-v-T (Begründung Abschnitt 2); eine Umstellung der Datenspalte auf m³/kg würde auch p-v-T/u-v-T-Konventionen, `ThermoState.v_molar` und den Zustandsvergleich berühren und wurde nicht vorgenommen.
2. **g-p-T-Auflösung der Genauigkeit:** g stammt aus dem HEOS-PT-Flash; Abweichungen gegen unabhängige Aufrufe bis ≈ 10⁻⁸ relativ (Dichteiteration). f-v-T (DmolarT) ist auf ≈ 10⁻¹³ reproduzierbar.
3. **Instabile Ebene:** mit Trennspalte und Lücken (EOS-Artefakte im tiefen instabilen Bereich, besonders Wasser bei T ≈ 0,75 Tc); 29 Wasserknoten wegen f < f_Gleichgewicht entfernt. Die instabile Fortsetzung ist keine Gleichgewichts- oder Stabilitätsaussage.
4. **g-p-T Log-Achse** blendet homogene Knoten außerhalb des Druckfensters aus (Option `symlog` zeigt sie); Gleichgewichtsflächen werden nie abgeschnitten.
5. **Z-Skalierung** ist wie in p-h-T/u-v-T fest in MJ/kg: bei leichten Fluiden mit kleiner g/f-Spanne (CO₂, N₂, Methan) erscheint die Fläche flach; die Achsenbeschriftung zeigt die echten Werte, Zoom/Kamera wie gewohnt.
6. **Sättigungsrandknoten in g-p-T** sind einseitige Zustände (Q nicht bestimmt); Klicks exakt auf die Linie liefern den Zustand der geklickten Seite.
7. **Zielrechner:** interaktive Bedienung und große Raster (`n_p`, `n_v`, hohe `homogeneous_points`) sind dort zu bestätigen; Cache-Speicher wie in TS-A (`limits.cache_max_mb`).
8. **Offene Punkte außerhalb des Auftrags:** der Fehlschlag `P8EndeZuEnde` und die Stub-Testverschmutzung durch `test_ap13_ap18` bestanden vor dieser Arbeit und wurden nicht angefasst; `Doku/manifest_sha256.json` ist nicht aktualisiert.

## 9. Interaktive gemeinsame Tangente

Die nachträglich ergänzte Temperatursteuerung, HEOS-Gleichgewichtstangente, nicht-modale 2D-Isothermenansicht und Validierung sind in [README_FVT_TANGENTE.md](README_FVT_TANGENTE.md) dokumentiert. Die Begleitfunktionen für Viewport-Recovery, Hover, Isolinien und diagrammübergreifende Auswahlmarker stehen in [README_QOL_UND_ZUSTANDSPROJEKTION.md](README_QOL_UND_ZUSTANDSPROJEKTION.md).
