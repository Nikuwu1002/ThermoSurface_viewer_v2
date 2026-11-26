# E4 – Chemisches Potenzial im 2D-Schnitt und Isolinien im 3D-Modell beschriften

**Projekt und Grundlage:** ThermoSurface. **Module:** `fvt_slice_view.py`, `fvt_tangent.py`, `gui.py`, `visualization.py`, `diagram_viewport.py`, `isoline_levels.py`. 2D unabhängig. Für 3D die tatsächlich erzeugten Linienwerte aus E3 übernehmen. Die Abnahmetests sind vorgesehen und noch nicht durchgeführt.

## 1. Ziel und Umfang

Im f-v-Schnitt den Zahlenwert von μ in J/mol am Punkt vₘ = 0 anzeigen. In den sieben 3D-Diagrammen sollen Nutzer einzelne berechnete Isothermen und Isobaren auswählen können, deren Wert direkt an der dargestellten Kurve erscheint. allerdings auch ausgeblendet werden kann mit einer checkbox

## 2. Ist-Zustand im Code

`calculate_fvt_tangent()` in `fvt_tangent.py` bestimmt das molare chemische Potenzial über $\mu=f_{m,l}+p_\mathrm{sat}v_{m,l}$. `tangent_plot_data()` liefert `intercept_volume=0.0` und `intercept=tangent.mu`. `FvtSliceView.update_tangent()` zeichnet das Kreuz bereits an diesem Punkt, annotiert es aber nur mit `μ = gₘ`. Der 2D-Plot verwendet molare Energie. `tangent_display_points()` dividiert die molare Energie für den **3D**-Renderer durch die molare Masse, da dessen f-Achse massenspezifisch ist.

In 3D erzeugen `_draw_curves()`, `_pht_add_curve()`, `_pst_add_curve()` und `_new_curve()` die Linien-Actors. Die Renderer sammeln sie nach Familie. Eine Kurve kann durch ungültige Werte oder Clipping in mehrere Segmente zerfallen. Bei Zusatzdiagrammen führt `build_isolines()` die erfolgreich berechneten Niveaus als `computed_levels`. Bei p-v-T wird in `berechne_daten()` die Abtastung separat über `record_sampling()` protokolliert. Eine flache Actor-Liste genügt deshalb nicht zur sicheren Zuordnung „Actor Nummer i = Niveau i“.

## 3. Umsetzungsschritte

1. **2D-Annotation ändern.** In `FvtSliceView.update_tangent()` für Marker und Annotation dieselben Koordinaten `data["intercept_volume"]`, `data["intercept"]` verwenden. Den Text beispielsweise als `μ = gₘ = 1234,56 J/mol` formatieren. Zahlenwert direkt aus `data["intercept"]` übernehmen. Eine Umrechnung über die molare Masse wäre hier falsch.
2. **Aktualisierung prüfen.** Den Weg von `_update_fvt_tangent()` über direkte Anzeige bzw. `FvtIsothermSampler` und `_on_fvt_isotherm_sampled()` verfolgen. Auch bei schnellem Verschieben des Temperaturschiebers müssen Kurve, Titel und Zahlenwert zum selben `FvtTangent` gehören. `show_error()` muss eine vorherige gültige Annotation mit entfernen.
3. **3D-Linienidentität beim Erzeugen festhalten.** Als neue Datenstruktur je erzeugter Linie Familie, exakten SI-Level und alle zugehörigen Segment-Actors speichern. Beim p-v-T-Pfad den Level zusammen mit der Kurve weiterreichen. Bei den Zusatzdiagrammen `computed_levels` verwenden und bei der h-s-p-Transformation die Zuordnung mitführen. So bleibt auch eine Isotherme im h-s-p-Diagramm identifizierbar, dessen drei Achsen keine Temperaturkoordinate enthalten.
4. **Auswahl im aktiven Viewport anbieten.** Im GUI-Block eine Mehrfachauswahl der tatsächlich erzeugten Isothermen und Isobaren anzeigen. Schlüssel ist `(Familie, SI-Level)` innerhalb des aktiven Diagramms. Für den Text Temperatur in K und Druck in der für das Diagramm festgelegten Anzeigeeinheit formatieren. Die GUI-Auswahl soll keine HEOS-Neuberechnung auslösen.
5. **Anker auf dem sichtbaren Segment bestimmen.** Nach Anwendung der manuellen Box die aktuellen Linien-Datasets der zugeordneten Actors auswerten. Das längste zusammenhängende gültige Teilstück nach Bogenlänge wählen und den Anker entlang seiner Liniensegmente auf halber Länge interpolieren. Ein Mittelwert aller Kurvenpunkte kann neben der Kurve oder innerhalb einer Lücke liegen.
6. **Labels im Viewport verwalten.** Kameraorientierte Punktlabels nach dem im Projekt bereits verwendeten Muster `add_point_labels()` erzeugen, getrennt von Clip- und Flächen-Actors. Nach `show_result()`, `apply_layers()` und `_apply_view()` Sichtbarkeit und Anker aktualisieren. Bei leerem Segment das Label entfernen. Bei Neuberechnung alte Actors vollständig freigeben und nur noch vorhandene ausgewählte Niveaus wiederherstellen.

## 4. Abnahmetests

| Ausgangszustand und Bedienhandlung | Soll-Ergebnis |
|---|---|
| f-v-T-Schnitt bei zwei gültigen subkritischen Temperaturen öffnen | Annotation liegt jeweils bei (0, μ). Zahlenwert stimmt innerhalb der letzten angezeigten Dezimalstelle mit `tangent_plot_data(tangent)["intercept"]` überein |
| Temperatur rasch mehrfach wechseln, homogene Schnittkurve aktiviert | Nach Abschluss zeigt der Schnitt für Titel, Tangente, homogene Kurve und Annotation dieselbe letzte Temperatur |
| Zwei Isothermen und eine Isobare auswählen, Kamera drehen und zoomen | Drei lesbare Labels mit korrekten Werten. Jeder Anker liegt auf einem Segment seiner eigenen Kurve |
| Kurve teilweise und anschließend vollständig clippen, Familie ausblenden, neu berechnen, Diagramm wechseln | Anker wandert auf ein verbleibendes Segment oder verschwindet. Keine veralteten oder doppelten Labels. Die Auswahl eines anderen Viewports bleibt eigenständig |

## 5. Abhängigkeiten und offene Entscheidungen

Die 2D-Änderung lässt sich unabhängig durchführen. Für 3D sollte die stabile Zuordnung von Level und Kurve mit E3 geteilt werden. Die Labelauswahl gehört zum Darstellungszustand des Viewports und damit nicht in den Berechnungs-Cache.
