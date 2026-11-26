# AP5 – Einbindung von p-h-T-Isolinien

**Ziel und Abnahme:** Die in `coolprop_heos_pht.py` schon vorhandenen `calc_isotherme`, `calc_isobare` und `calc_isochore` sollen in die reguläre Ergebnisübergabe gelangen und in PyVista/GUI separat sichtbar sein. Keine Neubestimmung der grundlegenden thermodynamischen Methoden.

## Umgesetzte Änderungen

- `coolprop_heos_pht.py`: `berechne_pht_daten(..., include_isolines=True)` ergänzt die drei Listen `isothermen`, `isobaren`, `isochoren` im Ergebnis-Dictionary. Mit `False` bleiben die Listen leer; Flächen-/Binodalberechnung bleibt unverändert. Neue `_pht_isolinien(stoff, raster, critical)` berechnet standardmäßig sieben ausgewählte Isothermen, sechs Isobaren, fünf Isochoren durch Wiederverwendung vorhandener Einzelfunktionen. Die 3D-Punkte bleiben im Format `(T [K], h [J/kg], p [Pa])`.
- `visualization.py`: `render_pht()` konstruiert die drei separaten Actor-Listen. Nichtendliche Zustands- oder Transformationskoordinaten unterbrechen einen Linienzug, statt als unzulässige Verbindungen dargestellt zu werden. Die ursprüngliche Flächen- und Phasengrenzvisualisierung bleibt erhalten.
- `gui.py`: Die p-h-T-Ebenen `isothermen`, `isobaren`, `isochoren` sind separat schaltbar und im Standardprofil aktiviert. `spinodale` bleibt in p-h-T deaktiviert.
- `coolprop_heos_hsp.py`: HSP ruft `berechne_pht_daten(..., include_isolines=False)` auf, um zusätzliche, für die HSP-Linien nicht benötigte p-h-T-Auswertungen zu vermeiden. Der ursprüngliche HSP-Rechenweg über p-s-T bleibt erhalten.
- `main.py`: Keine Signaturänderung notwendig. Dispatcher ruft die p-h-T-Funktion wie bisher auf; deren standardmäßig aktive optionale Linien werden jetzt mitgeliefert.

## Tests und Grenzen

Der Test `test_ap5_ap8.AP5IntegrationTests` prüft per Original-AST die p-h-T-Signatur und die HSP-Auslassung, führt die originale Linien-Sammelfunktion mit drei Mock-Einzelberechnern aus und rendert den Original-`render_pht` mit einem Plotter-/PyVista-Mock auf getrennte Actor-Gruppen (auch bei NaN-Unterbrechung). **4 von 4 ausführbaren AP5-Tests bestanden.** Zusätzlich sind bisherige synthetische AP4-Tests und AP1-Testlogik gelaufen.

Der zusätzlich bereitgestellte `python -m unittest test_ap5_ap8.AP5HEOSTests -v` prüft die tatsächliche HEOS-Erzeugung für Water und R134a sowie die Linienkonstanz und den optionalen Ausschalter. In dieser Umgebung wurde dieser eine Test **als SKIP** protokolliert.

**Offen:** Der reale HEOS-Durchlauf einschließlich geometrischer Lage aller drei Linien auf der p-h-T-Fläche, reale PyVista-Ausgabe und der GUI-Schaltertest am Zielrechner. Ursache: CoolProp, PyVista, PyVistaQt und PySide6 stehen hier nicht zur Verfügung. Für eine physikalische und visuelle Freigabe `python -m unittest test_ap4.AP4HEOSTests test_ap4.AP4RenderTests -v` und den GUI-Smoke ausführen, zusätzlich die p-h-T-Linien prüfen.

**Status:** Softwareseitige Integration implementiert und in Teilpfaden synthetisch getestet; **nicht vollständig HEOS-/GUI-validiert**.
