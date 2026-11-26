# TS-VP05 – Universelles 3D-Clipping der vorhandenen Actors

**Module:** `geometry_clip.py` (NumPy), `viewport_clipping.py` (PyVista), `diagram_viewport.py`, `state_picker.py`, `state_resolver.py`, `zustandsflaeche.py` · **Tests:** `test_ts_vp04_05_06_viewport.py`, `test_ts_vp_gui.py`

## 1. Verfahren

1. **Box:** Die sechs SI-Grenzen werden mit **derselben Vorwärtstransformation wie die Meshpunkte** (`diagram_axes.to_display`, Grundlage `diagram_coords.py`) in eine VTK-Box umgerechnet (alle Achsentransformationen sind monoton steigend).
2. **Originale:** Beim Rendern registriert `ClipRegistry.register` das Eingabedatenset jedes Geometrie-Actors. **Jedes Clipping geht vom Original aus** – kein kumulatives Wegschneiden; Zurücksetzen setzt die Originale per `SetInputData` wieder ein (V03, bitgleich geprüft).
3. **Flächen:** Dreiecke werden an den sechs Boxebenen exakt geschnitten (Sutherland–Hodgman je Ebene, vektorisiert); Schnittkanten liegen auf den Boxflächen, Wicklungssinn bleibt erhalten. **Kurven:** Liang–Barsky je Segment; zusammenhängende, nicht geschnittene Segmente bleiben eine Polylinie. **Punkte** (kritischer Punkt): In-Box-Prüfung.
4. **Darstellung bleibt:** Actor-Objekte, Farbe, Transparenz, `show_edges`, Layer-Sichtbarkeit bleiben erhalten (nur das Mapper-Eingabedatenset wird getauscht); Normalen für `smooth_shading` werden aus dem geclippten Netz neu berechnet.
5. **Punktdaten:** Arrays (`T`, `v`, `p_raw`, `source`, Zelldaten) werden mitgeführt. Neue Schnittpunkte erhalten `source = SRC_CLIPPED (2)`; ihre `T/v/p_raw` stammen aus der Rückwärtstransformation der Darstellungskoordinaten – **keine HEOS-Auswertung**.
6. **NaN-Lücken** werden nie geschlossen: Zellen mit nichtendlichen Koordinaten werden verworfen; es entstehen keine künstlichen Brücken, keine neuen Zellen außerhalb der Ursprungsdreiecke.
7. **Leere Fenster:** Actors erhalten ein leeres Datenset (kein Renderfehler); Rahmen, Ticks und Titel bleiben sichtbar; die Statuszeile meldet „Ansichtsbereich enthält keine Geometrie“.
8. **Nicht geclippt:** Achsen, Ticklabels, Orientierungswidget.

## 2. Zustände, Picker und Marker (V11)

* Geclippte Netze tragen das Zellarray `ts_orig_cell` (Index des Ursprungsdreiecks). `pick_from_viewport` rekonstruiert bei einem Treffer auf ein beschnittenes Netz **aus dem unveränderten Ursprungsdreieck** (Pickposition liegt darin); `MeshPick.clipped = True`, `state.provenance["viewport_clipped"] = True`.
* Dadurch sind Schnittpunkte nie „fälschlich originale HEOS-Rasterknoten“: ein exakter Knoten wird nur erkannt, wenn die Pickposition in der Originaldreieckskoordinate tatsächlich an einem Knoten liegt; sonst Interpolation wie bisher.
* Hover (p-v-T): Schnittpunkte werden als „Schnittpunkt (Renderer-Interpolation) – keine HEOS-Auswertung“ gekennzeichnet (`zustandsflaeche.hover_label`).
* Zustandsvergleichsmarker außerhalb der Box werden ausgeblendet, aber nicht aus der Auswahl gelöscht; beim Zurücksetzen erscheinen sie wieder.

## 3. Grenzen und Nichtumfang

* Es wird nur auf bereits berechneten Daten geschnitten; Limits außerhalb des Datenbereichs lassen Leerräume (kein Nachrechnen, keine Verdichtung).
* Die Schnittgeometrie ist eine lineare Interpolation im Darstellungsraum. Bei starker Krümmung der Transformation (z. B. Log-Achsen) ist sie nicht dieselbe Kurve wie eine neu berechnete HEOS-Fläche – genau deshalb ist sie als Renderer-Geometrie gekennzeichnet.
* Dreiecke der Netze sind unverändert; es gibt keine Neuvernetzung.

## 4. Abnahme

Alle gerenderten Geometriepunkte liegen innerhalb der Box (Toleranz 10⁻⁹ relativ; GUI-Test `test_V02…`), Flächeninhalte stimmen exakt (`TriangleClipping::test_area_is_exact…`, `test_z_clipping_of_an_inclined_plane`), kein kumulatives Wegschneiden (`test_repeated_clipping_always_starts_from_the_original`, GUI `test_V03…`), kein neuer HEOS-Aufruf (GUI-Test und `scripts/ts_vp07_performance.py`: 0 Zustandserzeugungen), Picker auf beschnittenen Flächen (`test_V11…`).
