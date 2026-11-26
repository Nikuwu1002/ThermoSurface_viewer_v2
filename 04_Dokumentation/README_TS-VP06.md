# TS-VP06 – Physikalische Achsen, benutzerdefinierte Ticklabels und Kameraführung

**Module:** `01_Programm/axes_view.py` (PyVista), `01_Programm/diagram_axes.py` (VTK-frei), `01_Programm/visualization.py`, `01_Programm/diagram_viewport.py` · **Tests:** `03_Tests_Implementierung/test_ts_vp04_05_06_viewport.py`, `03_Tests_Implementierung/test_ts_vp_gui.py`

## 1. Zentrale Achsenaktualisierung

`axes_view.apply_axes(plotter, diagram_id, daten, actors, box=None)` ersetzt alle diagrammspezifischen `show_bounds()`-Aufrufe und `update_pvt_axes()` (letzteres ist ein dünner Wrapper). Die Geometrie-Transformationen (`log10(v/vc)`, `T/Tc`, `asinh(p/pc)`, `(h-hc)/1e6` …) bleiben unverändert.

* **Rahmen:** Quader aus den **sechs manuellen Grenzen** (Viewport aktiv) – unabhängig von Layer-Sichtbarkeit und vorhandener Geometrie – sonst aus den automatischen Datengrenzen (p-v-T: sichtbare Ebenen wie bisher, alle anderen: gesamte Geometrie, stabil bei Layerwechsel).
* **Ticks:** physikalische Zahlen, Positionen über **dieselbe Vorwärtstransformation wie die Meshpunkte** (`diagram_axes.to_display`). Lineare Abstände physikalischer Tickwerte erscheinen auf Log-Achsen ungleich – das ist erwartet.
* **Titel:** `v [m³/mol] (log)`, `T [K]`, `p [MPa] (asinh)`, `h [kJ/kg]`, `s [kJ/(kg K)]`, `u/g/f [kJ/kg]`; der Zusatz zeigt die Skalenart (log, asinh, symlog; linear ohne Zusatz). Doppelte VTK-Beschriftungen entfallen (kein `vtkCubeAxesActor` mehr).
* **Lesbarkeit:** Ticks/Titel sind VTK-2D-Labels (`always_visible`, immer zur Kamera gerichtet), außerhalb der Box versetzt. Die Kanten (vordere Bodenkante für x/y, linke senkrechte Kante für z) folgen der Kamera und werden nach jeder Interaktion (`EndInteractionEvent`) neu gesetzt – die Kamera wird dabei nicht verändert.

## 2. Tickregeln (`diagram_axes.physical_ticks`)

| Skala | Regel | Beispiel |
|---|---|---|
| linear (T, h, s, u, g, f, p linear) | 1-2-2.5-5-Raster, ca. 6 Ticks | p 4…26 MPa → 5, 10, 15, 20, 25 |
| log (v, p-log) | Dekaden; bei < 3 Dekaden zusätzlich 2 und 5 | v 10⁻⁵…10⁻² → 1e-5, 1e-4, 1e-3, 1e-2 |
| asinh / symlog (p) | 0 und symmetrische Dekaden bzw. 1-2-5-Werte; Ticks im nahezu linearen Bereich (|p| < 0.01 pc) werden ausgeblendet | p −5…30 MPa → −5, −2, −1, 0, 1, 2, 5, 10, 20 … |

Formatierung: `format_tick` (z. B. `1e-5`, `972.117`, `0`). Druckachsenmodi `log`, `linear`, `symlog`, `asinh` und die Vorzeichengültigkeit sind in `BoundsValidation` und `Ticks` geprüft.

## 3. Kamera

* Beim ersten Anwenden bzw. „Ansicht anwenden“: `plotter.reset_camera(bounds=box)` – Mittelpunkt und Ausdehnung folgen der Box; **Blickrichtung, `view_up` und Projektionsart bleiben** (`axes_view.fit_camera`).
* **Layer-Toggle, Gitter-Toggle, Auswahlmodus** setzen weder Kamera noch Box zurück (GUI-Test V05 vergleicht `camera_state()` und Rahmen-Bounds nach jedem Umschalten).
* Die Standardkamera des Renderers wird separat gespeichert (`DiagramViewport.default_camera`) und über „Standardkamera“ wiederhergestellt; Orientierung ist nicht Teil des Darstellungsfensters (`ViewportState.camera_state`).
* „Zurücksetzen“ passt die Kamera auf die automatischen Grenzen an (Blickrichtung bleibt).

## 4. Abnahme

Alle sieben Diagramme zeigen physikalische Achsengrenzen/Ticks/Titel mit Einheiten (Beispielbilder: `03_Tests_Implementierung/figures/ts_vp07_<id>_voll.png`, `…_ausschnitt.png`, erzeugt mit `03_Tests_Implementierung/scripts/ts_vp07_beispielbilder.py`); Grenzen ändern sich nicht beim Layerwechsel; Achsen bleiben bei leerem Fenster sichtbar (V04).
