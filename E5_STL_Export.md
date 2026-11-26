# E5 – STL-Export des reinen Oberflächennetzes

**Projekt und Grundlage:** ThermoSurface. **Module:** `gui.py`, `diagram_registry.py`, `diagram_viewport.py`, `visualization.py`, `viewport_clipping.py`, `diagram_coords.py`. Vorhandene Renderer und geometrisches Clipping reichen als Grundlage. Die Abnahmetests sind vorgesehen und noch nicht durchgeführt.

## 1. Ziel und Umfang

Eine Aktion „Oberflächen als STL exportieren“ exportiert die im aktiven Viewport sichtbaren und ausgewählten Flächen. Die STL soll die dargestellte dreidimensionale Form einschließlich eines angewendeten geometrischen Beschnitts reproduzieren. Kamerazoom, Drehung und Verschiebung ändern die exportierte Geometrie nicht.

## 2. Ist-Zustand im Code

`DIAGRAMS` in `diagram_registry.py` registriert sieben Diagrammtypen. **Alle sieben besitzen exportierbare Flächen.** g-p-T besitzt keine Zweiphasenfläche, aber reguläre Einphasenflächen. Die Berechnungsergebnisse sind überwiegend NumPy-Arrays. Die Renderer bauen daraus Netze und geben Actors zurück. `DiagramViewport._surface_actors()` erfasst nur den Schlüssel `surface` und wäre als alleinige Exportquelle unvollständig.

| Diagramm | Exportierbare Flächen und Datenquelle | Vorliegender Typ und Aufbereitung |
|---|---|---|
| p-v-T | Einphasenfläche aus `create_single_phase_surface()` mit binodalkonformen Rändern. `two_phase_surface` über `add_two_phase_surface()`. Metastabile Flüssigkeit, metastabiler Dampf, instabile HEOS-Fläche und `replacement_surface` über `homogeneous_shared_meshes()` | Rohdaten als Arrays bzw. homogene Felddaten. Renderer erzeugt Teilnetze und zusammengeführte Datasets. Aktuellen Mapper-Datensatz jeder zulässigen Ebene in Dreiecks-`PolyData` überführen. `surface_single_phase` ist hier ein leerer Kompatibilitätseintrag |
| p-h-T | `surface_liquid`, `surface_vapor`, `surface_high_temperature`, `two_phase_surface` | NumPy-Arrays → `_pht_surface_mesh()` bzw. `_pht_two_phase_mesh()` → triangulierte Oberflächen mit eigener kritischer Spitze |
| p-s-T | Dieselben vier Ergebnisschlüssel im p-s-T-Ergebnis | NumPy-Arrays → `_pst_surface_mesh()` bzw. `_pst_two_phase_mesh()` → triangulierte Oberflächen |
| h-s-p | Drei Einphasenflächen und `two_phase_surface` aus `berechne_hsp_daten()` | NumPy-Arrays → `_new_surface_mesh()` bzw. `_new_phase_mesh()` mit h-s-p-Transformation → triangulierte Oberflächen |
| u-v-T | Drei Einphasenflächen und `two_phase_surface` aus `berechne_uvt_daten()` | NumPy-Arrays → `_new_surface_mesh()` bzw. `_new_phase_mesh()` mit u-v-T-Transformation → triangulierte Oberflächen |
| g-p-T | Drei Einphasenflächen, optional `metastable_liquid_surface`, `metastable_vapor_surface`, `unstable_surface`. Das Sättigungsgleichgewicht erscheint als Linie | Reguläre Arrays über `_new_surface_mesh()`, optionale Arrays mit ungültigen Bereichen über `_masked_grid_mesh()`. Nur die daraus erzeugten Flächen übernehmen |
| f-v-T | Drei Einphasenflächen, die Zweiphasen-Regelfläche und die drei optionalen homogenen Flächen | Reguläre Arrays über `_new_surface_mesh()` bzw. `_new_phase_mesh()`, maskierte homogene Arrays über `_masked_grid_mesh()`. Tangente und Schnittebene gesondert ausschließen |

## 3. Umsetzungsschritte

1. In `gui.py` den Export an `workspace.active` und ein gültig dargestelltes Ergebnis binden. Im Dialog die verfügbaren, aktuell sichtbaren Flächen dieses Viewports mit ihren Ebenennamen anbieten. Sichtbare Flächen sind vorgewählt. Die synthetische p-v-T-Ersatzfläche ausdrücklich als Modellfläche kennzeichnen.
2. In `diagram_viewport.py` eine **neue** Funktion zum Einsammeln der freigegebenen Flächen vorsehen. Zulässige Schlüssel sind `surface`, `zweiphasengebiet`, `metastable_liquid_surface`, `metastable_vapor_surface`, `unstable_surface`, `replacement_surface`. Einzel-Actors und Actor-Listen behandeln. Vorhandensein, tatsächliche Sichtbarkeit und nichtleere Flächenzellen prüfen. Eine Positivliste schließt auch die f-v-T-Schnittebene zuverlässig aus, obwohl sie geometrisch selbst ein Polygon ist.
3. Den aktuellen Mapper-Eingang auslesen. `ClipRegistry.apply()` setzt diesen bei einer manuellen Box auf das beschnittene Dataset und beim Reset auf `entry.original`. Die aktuelle Kamera beeinflusst den Mapper-Eingang nicht. So stimmen Export und geometrisch sichtbarer Ausschnitt überein.
4. In **VTK-Darstellungskoordinaten** exportieren. Die Renderer in `visualization.py` verwenden die Transformationen aus `diagram_coords.py`, einschließlich der wirksamen Druckskalierung. Eine Rücktransformation in SI würde bei logarithmischen oder anderen nichtlinearen Achsen eine andere Form erzeugen. Im Exportdialog die Bedeutung als „dargestellte Form, dimensionslose Modellkoordinaten“ angeben. STL speichert keine thermodynamischen Achseneinheiten.
5. Jeden Datensatz als Oberfläche extrahieren, triangulieren und ausschließlich gültige, endliche Dreiecksflächen übernehmen. Punkte, Linien und entartete Dreiecke verwerfen. Teilnetze ohne künstliche Verbindungsflächen zusammenführen. Das vorbereitete `pyvista.PolyData` speichern. Bei leerer Auswahl oder vollständig weggeclippten Flächen keine Datei erzeugen.

## 4. Abnahmetests

| Ausgangszustand und Bedienhandlung | Soll-Ergebnis |
|---|---|
| Je registriertem Diagramm eine Einphasenfläche exportieren und STL erneut öffnen | Reimport enthält ausschließlich Dreiecke. Dreiecksanzahl stimmt mit dem vorbereiteten Exportnetz überein. Die sechs räumlichen Grenzen stimmen innerhalb einer zur STL-Koordinatenpräzision passenden Toleranz überein |
| p-v-T mit eingeblendeter Zweiphasenfläche und Modellfläche, zusätzlich Isolinien, Marker und Labels sichtbar | STL enthält genau die gewählten Flächen. An- und Ausschalten von Linien, Markern oder Labels ändert die exportierte Dreiecksmenge nicht |
| Erst nur Kamera drehen und zoomen, danach manuelle Box anwenden | Nach Kamerabewegung bleiben Dreiecksgeometrie und Grenzen gleich. Nach geometrischem Clipping liegen alle exportierten Punkte innerhalb der Box zuzüglich Rundungstoleranz |
| g-p-T öffnen bzw. nur Linien sichtbar lassen | g-p-T bietet seine Einphasenflächen und vorhandene homogene Flächen an. Bei leerer Flächenauswahl erhält der Nutzer eine konkrete Meldung und es entsteht keine leere STL |

## 5. Abhängigkeiten und offene Entscheidungen

Empfohlen ist der Export der dargestellten Oberfläche einschließlich offener Ränder. Ein druckbarer, geschlossener Volumenkörper würde zusätzliche Wandstärken, Deckflächen und eine Modellskalierung erfordern. Das wäre ein eigener Auftrag. Der vorhandene Export in `on_export_homogeneous()` hat einen anderen Datenumfang und sollte als eigenständige Funktion erhalten bleiben.
