# AP11 – Bestandsaufnahme und Ursachenanalyse

## Istzustand und Datenfluss

`gui.CalculationWorker.run()` ruft `main.berechne_auswahl()` und im Modus `pvt` `main.berechne_daten()` auf. In `coolprop_heos_pvt.py` werden zunächst die gemeinsamen Temperaturstützstellen (`build_pvt_temperature_grid`), die einphasig beziehungsweise als Gleichgewicht ausgewertete Oberfläche (`multi_isotherm`), die unterkritischen Binodalen (`multi_binodale`) und die Maske (`berechne_zweiphasenmaske`) berechnet. Die zweiphasige Gleichgewichtsfläche (`calc_two_phase_surface`) entsteht **vor** der späteren Ergänzung des kritischen Punktes an den Binodalkurven. Die Spinodalen entstehen separat über `calc_spinodale` und `calc_spinodale_vapor_temperature`, mit der bisherigen Auswahl- und Validierungslogik in `main.py`.

`visualization.render_pvt` erzeugt aus `create_single_phase_surface` den stabilen Flächen-Actor `surface` (Randknoten aus den Binodalen), aus `add_two_phase_surface` den Actor `zweiphasengebiet` und die separat berechneten Spinodal-Actors. Die GUI schaltet Actors um, ohne die Zustandspunkte neu zu berechnen. Die anderen Diagrammtypen benutzen die eigenständigen Module `coolprop_heos_pht.py`, `coolprop_heos_pst.py`, `coolprop_heos_hsp.py` und `coolprop_heos_uvt.py`.

## Datenvertrag im Istzustand

| Schlüssel | Arrayform | Erzeuger | Bedeutung |
|---|---|---|---|
| `surface` | `(n_T,n_v,3)` | `multi_isotherm()` | T-v-p-Daten des Gleichgewichts-/Einphasenpfades |
| `phase_mask` | `(n_T,n_v)` Boolean | `berechne_zweiphasenmaske()` | Zweiphasenmaske |
| `binodal_liquid/vapor` | `(n_sat+1,3)` | `multi_binodale` + `main` | Unterkritische Sättigungszustände und Tc |
| `two_phase_surface` | `(n_sat,n_Q,3)` | `calc_two_phase_surface()` | Gleichgewichtsmischung aus Siedelinie und Taulinie |
| `spinodale_liquid/vapor` | `(n_s,3)` | HEOS-Tracer/Temperatursolver | mechanische Grenzen homogener Fortsetzung |
| `critical_point` | Tupel `(Tc,pc,vc)` | `coolprop_heos_common` | Kritische Referenzgrößen |

Koordinaten aller p-v-T-Zustandsarrays: `(T [K], v [m³/mol], p [Pa])`. Die Darstellung transformiert in `(log10(v/vc), T/Tc, log10(p/pc))`.

## Ursachenklassifikation

**Thermodynamisch:** Die Gleichgewichtsfläche hat `p=p_sat(T)` und `v=(1-Q)v_l+Qv_g`. Die Spinodalen besitzen dagegen `p=p_hom(T,v_spin)`; sie sind keine Ränder dieser Gleichgewichtsfläche. Eine Abweichung in `p` ist somit physikalisch begründet. Quantifizierbar bei gleichem `T` sind `Δp = p_spin-p_sat` und für positive Drücke `Δz=log10(p_spin/p_sat)`. Ein geometrischer Druckabgleich wäre falsch.

**Numerisch:** Der Dampf-Temperatursolver liefert explizite NaN für fehlende Wurzeln, der Flüssigkeitszweig war bislang vorwiegend auf den Tracer angewiesen. Die Temperaturabdeckung beider Spinodaläste ist nicht automatisch die der Binodalen.

**Geometrisch:** Die alte Zweiphasenfläche erhält ihren kritischen Abschluss über einen separaten Dreiecksfächer. `_mesh_from_rows` trianguliert die einphasigen Teilnetze mit gemeinsamen Binodalknoten. Ein Gleichgewichtsnetz hat aber keine gemeinsamen Spinodalknoten. Die Referenz-Abbildung 71 zeigt eine **unkorrigierte van-der-Waals-Zustandsfläche** und damit eine andere thermodynamische Teilfläche. Unterschiedliche Blickwinkel können die Lücke optisch verbergen, nicht physikalisch beheben.

## Prüfung und Grenzen

**Prüfung A:** Datenfluss und Einheiten durch vollständige Modulinspektion und Analyse der vorliegenden Projektdateien zugeordnet. **Prüfung B:** Renderpfade und bestehende AP2-Randknotenlogik unabhängig nachverfolgt; keine Produktivänderung in AP11. Für einen quantitativen **Wasser**-Datensatz und eine echte VTK-Dreiecksprüfung fehlen hier CoolProp und PyVista. Das quantitative Akzeptanzkriterium ist daher noch **offen**, nicht bestanden. Die identifizierten Eingriffspunkte sind `main.berechne_daten`, `homogeneous_pvt` (neu ab AP13), `visualization.render_pvt` und `gui.layer_checks`.
