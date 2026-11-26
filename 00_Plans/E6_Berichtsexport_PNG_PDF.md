# E6 – Bericht aus dem aktiven Viewport als PNG und PDF exportieren

**Projekt und Grundlage:** ThermoSurface. **Module:** `gui.py`, `diagram_workspace.py`, `diagram_viewport.py`, `calculation_manager.py`, `diagram_cache.py`, `main.py`, `selection_manager.py`, `thermo_state.py`. Eigenständig umsetzbar. Beschriftungen aus E4 werden mit dem Viewport aufgenommen. Die Abnahmetests sind vorgesehen und noch nicht durchgeführt.

## 1. Ziel und Umfang

Ein Button „Bericht exportieren“ erzeugt ein hochauflösendes PNG des aktiven 3D-Viewports und ein PDF mit demselben Bild, Stoffdaten, wirksamen Parametern und ausgewählten Zustandspunkten. PNG und PDF müssen auch bei Mehrfachansicht, Cache-Treffer und noch nicht übernommenen GUI-Eingaben eindeutig zur dargestellten Berechnung gehören.

## 2. Ist-Zustand im Code

`DiagramWorkspace.active` identifiziert den aktiven Viewport. `selected_diagram_ids()` liefert dagegen die zur Berechnung ausgewählten Diagramme und kann mehrere IDs enthalten. `DiagramViewport.show_result(daten, key)` bindet `daten`, `result_key`, Renderer-Actors und Plotter an genau einen Viewport.

`TaskOutcome` verbindet den tatsächlich bearbeiteten `DiagramTask`, Ergebnisdaten, Cache-Status und `DiagramKey`. `DiagramKey` enthält Stoff, Diagramm-ID, Modell, Solverversion und `params_hash`. Aus dem Hash lassen sich die verwendeten Parameter nicht zurückgewinnen. Das ist die wesentliche Lücke für einen reproduzierbaren Bericht.

| PDF-Inhalt | Verbindliche Datenquelle |
|---|---|
| Stoff, Diagrammtyp, Modell, Solverversion | `viewport.result_key`, `viewport.diagram_id`, Bezeichnung aus `DIAGRAMS` |
| Kritische Stoffdaten | Diagrammspezifische `viewport.daten["meta"]` und `critical_point`. Bei p-v-T Reihenfolge `(Tc, pc, vc)` mit K, Pa, m³/mol. Bei Zusatzdiagrammen die jeweilige Konvention verwenden, etwa h_c, g_c, f_c, u_c in J/kg oder s_c in J/(kg·K) |
| Berechnungsparameter | Neu zu sichernde Kopie der effektiven Eingaben des tatsächlich ausgeführten Tasks einschließlich aufgelöster Diagrammoptionen. Erzeugte Raster und `isoline_report` aus genau diesem Ergebnis ergänzen |
| Darstellungsparameter | `effective_layers()`, tatsächliche Achsentransformation, `view_state`, `manual_box()`, aktuelle Kamera und Gridzustand dieses Viewports. Bei aktivierter Tangente bzw. Labels deren dargestellten Stand erfassen |
| Vergleichspunkte | `SelectionManager.points` mit `selection_index`, `selection_order`, `source_diagram`, zugehörigem `ThermoState`, Gültigkeit und verfügbaren Größen. Herkunft und Darstellbarkeit im exportierten Viewport getrennt ausweisen |

## 3. Umsetzungsschritte

1. **Angewendete Parameter am Ergebnis sichern.** Beim Auflösen der Berechnungsoptionen in `main.py` eine unabhängige Kopie der tatsächlich verwendeten Werte erzeugen und gemeinsam mit dem Ergebnis cachen. `TaskOutcome` bzw. `_after_render()` ordnet diese Daten dem angezeigten `result_key` zu. Den Snapshot auch bei Cache-Treffern und beim Wiederaufbau eines Viewports verfügbar halten. Ein späteres erneutes Lesen der YAML-Datei könnte bereits andere Einstellungen liefern.
2. **Exportziel eindeutig festlegen.** Beim Betätigen des Buttons `workspace.active` erfassen und über `has_result`, `render_state`, `scene_is_ready()` und `result_key` prüfen. Nach Rückkehr aus dem Dateidialog erneut prüfen, ob Viewport und Ergebnis noch identisch sind. Falls inzwischen ein neues Ergebnis dargestellt wird, den Export mit verständlicher Meldung abbrechen oder ausdrücklich den neuen Stand wählen lassen.
3. **Einen gemeinsamen Ansichtsstand aufnehmen.** In einem zusammenhängenden GUI-Schritt Kamera, wirksame Grenzen, Ebenen, Marker, Beschriftungen, Zustandspunkte und Ergebniskennung erfassen. Den bereits aufgebauten `viewport.plotter` rendern und mit einem festen Auflösungsfaktor, beispielsweise 3, aufnehmen. Kein zweites Diagramm aus den Rohdaten rendern. Seitenverhältnis, Kameraprojektion und aktuelles Clipping müssen erhalten bleiben. Die tatsächlichen PNG-Pixelmaße aus dem erzeugten Bild bestimmen.
4. **PDF aus diesem Snapshot erzeugen.** Exakt das erzeugte PNG einbetten. Auf der ersten Seite Stoff, Diagramm, Ergebniskennung und Bildabmessungen nennen. Anschließend Stoffdaten, Berechnungsparameter, Ansichtsparameter und Punkte mit Einheiten aufführen. Für die Punktliste empfiehlt sich die gesamte Vergleichsauswahl des zum Ergebnis gehörenden Stoffs in `selection_order`, jeweils gekennzeichnet als Quellenpunkt, sichtbare Projektion oder im aktiven Viewport nicht dargestellt. Fehlende Zustandsgrößen als „nicht verfügbar“ ausgeben. Ohne Punkte erscheint „Keine Vergleichspunkte ausgewählt“.
5. **Offene Eingaben sauber behandeln.** Ein noch nicht berechneter neuer Isolinienfaktor und nur eingetippte Ansichtsgrenzen dürfen nicht im Abschnitt „Angewendete Parameter“ landen. Ein bereits ausgeführtes Ein- oder Ausblenden einer Ebene gehört dagegen zum aufgenommenen Ansichtsstand. `get_settings()` ist deshalb keine geeignete alleinige Berichtsquelle.
6. **PNG und PDF als Paar schreiben.** Einmal Zielordner und gemeinsamen Basisnamen wählen, beispielsweise `Argon_pvt_20261126_143012`. Beide Dateien zunächst unter temporären Namen im Zielverzeichnis erzeugen und auf Lesbarkeit prüfen. Erst danach die endgültigen Namen vergeben. Vorhandene Namenspaare durch eindeutigen Suffix umgehen. Schlägt die zweite Umbenennung fehl, die erste neu angelegte Zieldatei zurücknehmen. Bei Abbruch des Dialogs entstehen keine Dateien. CSV bleibt eine optionale spätere Ergänzung.

## 4. Abnahmetests

| Ausgangszustand und Bedienhandlung | Soll-Ergebnis |
|---|---|
| p-v-T und f-v-T desselben Stoffs mit verschiedenen Kameras, Ebenen und Ausschnitten öffnen. f-v-T aktivieren und exportieren | PNG zeigt ausschließlich den f-v-T-Viewport. PDF nennt f-v-T, dessen Stoffdaten und dieselbe Ergebniskennung. Das eingebettete Bild stammt aus dem erzeugten PNG |
| Viewport mit bekanntem Seitenverhältnis, einer sichtbaren Markierung, Tangente und Label aufnehmen | PNG besitzt die vorgesehenen Pixelmaße. Ausschnitt, Kamerarichtung, projizierte Lage der Markierung sowie sichtbare Tangente und Label stimmen mit dem gewählten Viewport überein |
| Vor Export einen Isolinienfaktor verändern, ohne neu zu berechnen. Eine andere Ebene direkt ausblenden | PDF nennt den zuletzt **berechneten** Faktor. PNG und PDF zeigen die inzwischen ausgeblendete Ebene als unsichtbar |
| Mehrere Punkte aus zwei Diagrammen auswählen, einen projizieren und einen weiteren durch die manuelle Box ausblenden | PDF nennt korrekte IDs, Reihenfolge, Ursprung, Werte und Einheiten. Die Kennzeichnung der Sichtbarkeit entspricht dem PNG |
| Export nach Cache-Treffer sowie ohne Vergleichspunkte wiederholen | Angewendete Parameter sind auch beim Cache-Treffer vollständig. Ohne Punkte entstehen ein gültiges PNG und ein PDF mit entsprechender Angabe |
| Dateidialog abbrechen bzw. Schreibfehler beim zweiten Ausgabeschritt simulieren | Bei Abbruch entsteht keine Ausgabe. Bei Schreibfehler bleibt kein als erfolgreich gemeldetes, unvollständiges Dateipaar zurück. Vorhandene Berichte bleiben lesbar |

## 5. Abhängigkeiten und offene Entscheidungen

Empfohlen ist ein PDF mit eingebettetem PNG und einem gemeinsamen Basisnamen für beide Dateien. Dadurch lässt sich die Bildzuordnung auch nach dem Weitergeben des PDFs nachvollziehen. Die entscheidende Voraussetzung ist der **am Berechnungsergebnis gespeicherte Parametersatz**. Eine nachträgliche Rekonstruktion aus aktuellen GUI-Feldern oder allein aus `params_hash` kann diese Zuordnung nicht gewährleisten.
