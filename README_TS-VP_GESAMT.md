# TS-VP01 … TS-VP07 – Gesamtnote

**Stand:** 16.11.2026 · Umsetzung der „Verbindlichen technischen Spezifikation TS-VP01 bis TS-VP07“ (alle sieben Arbeitspakete).

| Paket | Ergebnis | Hauptdateien |
|---|---|---|
| TS-VP01 | Matrix der sieben Diagramme, Verträge, Baseline, Defizite | `README_TS-VP01.md` |
| TS-VP02 | zentraler Levelgenerator `x_n = start·f^n` | `isoline_levels.py`, `test_ts_vp02_levels.py` |
| TS-VP03 | Isolinien in allen sieben Berechnungspfaden, p-h-T-Renderer, GUI-Gruppe, Cache | `coolprop_heos_*.py`, `main.py`, `gui.py`, `visualization.py`, `config.yaml` |
| TS-VP04 | Checkbox + sechs Felder + Anwenden/Zurücksetzen je Diagramm, persistenter Zustand | `viewport_bounds.py`, `diagram_axes.py`, `gui.py`, `diagram_viewport.py` |
| TS-VP05 | exaktes 3D-Clipping aller Layer aus dem unveränderten Original | `geometry_clip.py`, `viewport_clipping.py` |
| TS-VP06 | physikalische Achsen mit Einheiten, zentrale Achsenaktualisierung, Kameraführung | `axes_view.py`, `diagram_axes.py` |
| TS-VP07 | Testmatrix, Performance, Beispielbilder, Dokumentation | `test_ts_vp_*.py`, `scripts/ts_vp07_*.py`, `results/ts_vp07_performance.json`, `figures/ts_vp07_*` |

Die Laufzeitmodule liegen in `01_Programm/`; Tests, Skripte und zugehörige Ergebnisse in `03_Tests_Implementierung/`.

**Kernaussagen.** Der Ansichtsbereich begrenzt nur die gerenderte Darstellung: kein CoolProp-Lauf, kein Eingriff in Ergebnis-Dictionary oder Cache-Schlüssel (Rohdaten bitgleich). Neue Schnittpunkte sind Renderer-Interpolationen und werden so gekennzeichnet; Zustände werden aus dem Ursprungsdreieck rekonstruiert. Isolinienlevels (konstante T/v/p) sind von der Abtastung entlang der Linie und vom Flächenraster getrennt; Änderungen der Level-Parameter verändern den Cache-Schlüssel der betroffenen Diagramme.

**Zu beachten.** Mit dem Standardfaktor 3.14 ergibt sich in engen Temperaturfenstern nur eine Isotherme (mathematische Folge der gewünschten Progression); der Faktor ist in der GUI und in `config.yaml` einstellbar. Ein vorbestehender Fehlschlag in `test_zustandsflaeche.py` (p0-Baseline, andere CoolProp-Version) ist unverändert gegenüber dem Ausgangszustand.
