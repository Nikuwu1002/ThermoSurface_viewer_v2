# TS-VP04 – Manuelle GUI und persistenter Viewportzustand

**Module:** `gui.py`, `diagram_viewport.py`, `diagram_workspace.py`, `viewport_bounds.py`, `diagram_axes.py` · **Tests:** `test_ts_vp_gui.py`, `test_ts_vp04_05_06_viewport.py`

## 1. Bedienung (Gruppe „Ansichtsbereich (nur Darstellung)“)

1. **Checkbox „Ansichtsbereich manuell begrenzen“** – zunächst deaktiviert (und ausgegraut, solange das aktive Diagramm kein Ergebnis hat).
2. Nach dem Aktivieren werden sechs Felder (Start/Stop je Achse) freigeschaltet; die Beschriftung folgt dem aktiven Diagramm, z. B. u-v-T: `x: v von [m³/mol]`, `y: T von [K]`, `z: u von [kJ/kg]`; h-s-p: `h [kJ/kg]`, `s [kJ/(kg K)]`, `p [MPa]`. Als Vorschlag stehen die aktuellen automatischen Grenzen in den Feldern.
3. Erst **„Ansicht anwenden“** wirkt auf das aktive Diagramm. Eingabe allein ändert nichts.
4. **„Zurücksetzen“** schaltet die manuelle Begrenzung aus: vollständige Originalgeometrie, automatische Grenzen, Kamera auf die gesamte Szene. Die zuletzt gültigen Eingaben bleiben als Vorschlag erhalten. Wird die Checkbox abgewählt, gilt dasselbe.
5. **„Standardkamera“** stellt die Standardkamera des Renderers wieder her (unabhängig vom Darstellungsfenster).

Fehler (Text, Start ≥ Stop, Start ≤ 0 bei T/v, p ≤ 0 auf Log-Achse, nicht darstellbare Transformation) erscheinen rot unter den Feldern und in der Statuszeile; der **zuletzt gültige Ausschnitt bleibt unverändert** (`ViewportState.apply_fields` übernimmt atomar, `DiagramViewport.apply_manual_view` rollt bei Darstellungsfehlern zurück). Es wird nie eine HEOS-Berechnung ausgelöst.

## 2. Koordinaten- und Einheitenvertrag (`diagram_axes.DIAGRAM_AXES`)

| Diagramm | x | y | z |
|---|---|---|---|
| pvt | v [m³/mol] | T [K] | p [MPa] |
| pht | h [kJ/kg] | T [K] | p [MPa] |
| pst | s [kJ/(kg K)] | T [K] | p [MPa] |
| hsp | h [kJ/kg] | s [kJ/(kg K)] | p [MPa] |
| uvt | v [m³/mol] | T [K] | u [kJ/kg] |
| gpt | p [MPa] | T [K] | g [kJ/kg] |
| fvt | v [m³/mol] | T [K] | f [kJ/kg] |

Eingaben werden in diese Anzeigeeinheiten gelesen und intern in SI gespeichert (`ViewportState.bounds`, z. B. `{"v": (1e-5, 1e-3), "T": (300, 800), "p": (1e5, 6e7)}`).

## 3. Prüfregeln (`diagram_axes.validate_bounds`)

* alle sechs Werte endlich, Start < Stop;
* T und v > 0; p > 0, wenn die Druckachse `log` ist (pht/pst/hsp: immer `log`; pvt/gpt: `pressure_axis` des Ergebnisses). Mit `asinh`, `symlog`, `linear` sind p ≤ 0 zulässig (V09);
* h, s, u, g, f dürfen negativ sein (V08);
* die Transformation der Grenzen muss endlich und streng wachsend sein.

## 4. Zustand je Diagramm und Persistenz

* Jeder `DiagramViewport` besitzt einen eigenen `ViewportState`; die GUI bearbeitet den **aktiven** Viewport (Klick auf den Kopf in der Mehrfachansicht).
* Die GUI speichert je Diagrammtyp eine Kopie (`_view_prefs`); beim Neuerzeugen eines Viewports desselben Typs (Wechsel Einzel ↔ Mehrfach, anderes Diagramm und zurück) wird sie übernommen, und die Box wird für die neuen/gleichen Daten wieder angewendet.
* Layerwechsel, erneutes Rendern desselben Ergebnisses (Cache-Treffer), Maximieren/Entmaximieren behalten den Zustand (Test V06/V07).
* Mehrfachansicht: zwei Diagramme mit unterschiedlichen Ausschnitten sind unabhängig.

## 5. Abnahme

Pro Diagramm kann ein eigener Grenzquader vorgegeben und auf dasselbe berechnete Ergebnis angewendet werden. Geprüft in `test_ts_vp_gui.py` (V02, V03, V04, V06, V07, ungültige Eingaben) und `test_ts_vp04_05_06_viewport.py` (V08, V09, Atomarität, Unabhängigkeit der Zustände). Cache-/Rohdaten-Identität (C01): `test_V02_apply_clips_original_geometry_without_any_recomputation` zählt die `berechnen`-Aufrufe, prüft `viewport.daten is daten`, gleichen `DiagramKey` und bitgleiche Rohdatenarrays.
