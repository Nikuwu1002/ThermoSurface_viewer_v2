# AP8 – Rasterkonvergenz und Performance

**Ziel:** Mindestens drei Rasterdichten auf Rechenzeit, Speicherbedarf und geometrische Interpolationsgenauigkeit hin untersuchen und anhand eines definierten Kriteriums eine Auflösung wählen.

## Technische Umsetzung

`ap8_convergence.py` implementiert eine Auswertung direkt über `AbstractState('HEOS',fluid)` mit `DmolarT_INPUTS`. Untersucht wird absichtlich nur der überkritische, gasähnliche und stabile Teilbereich `T/Tc ∈ [1,07;1,25]`, `v/vc ∈ [2;8]`. Die Raster sind `21×41 (861)`, `41×81 (3321)` und `81×161 (13041)` Punkte. Für den Rasterdruck wird der HEOS-Zustand an sämtlichen Knoten ausgewertet.

225 deterministisch gewählte Zusatzpunkte (`seed=42`) werden getrennt mit HEOS berechnet; `bilinear()` liefert zugehörige bilineare Interpolationen aus jedem diskreten Raster. Fehlerkennzahl `|p_interpoliert-p_HEOS|/p_HEOS`: Mittel, RMS und Maximum. Aus der `perf_counter`-Laufzeit mehrerer Rasterdurchläufe wird der Median gebildet. `grid_bytes` erfasst den Datenspeicher der Druckmatrix; die optionale Prozess-RSS von `psutil` wird nach Berechnungsdurchläufen aufgenommen. Dies ist keine exakte maximale Speicherbelastung und kein vollständiges GUI-Profiling.

`choose_resolution(rows,threshold=1e-3)` wählt die kleinste Rastergröße, deren gemessener *maximaler* relativer Interpolationsfehler die angegebene, selbstgewählte Grenze unterschreitet. Erfüllt kein getestetes Raster das Kriterium, lautet die Empfehlung ausdrücklich **keine**. Bericht und Messwerte werden nach `results/ap8` als CSV/JSON geschrieben. Der Vergleich betrifft eine homogene Teilfläche und kann daher nicht ohne weitere Messungen die optimale Auflösung der gesamten Anwendung mit Phasengeometrie und kritischem Punkt begründen.

## Aufruf und Verifikation

```powershell
python -m unittest test_ap5_ap8.AP8NumericalTests -v
python ap8_convergence.py --fluid Water --repetitions 2 --threshold 0.001 --output results/ap8
```

**Bestanden (2 Tests):** Der originale bilineare Interpolationscode reproduziert eine analytische Funktion `a+bT+cv+dTv` exakt an mehreren Kontrollpunkten; außerhalb des Gebiets wird eine Ausnahme ausgelöst. Die Rasterwahl wählt bei synthetischen Vorgabefehlern korrekt `41×81` beziehungsweise keine Auflösung.

**Noch offen:** Keine echte Laufzeit-/Speichermessung oder HEOS-Interpolationsfehlerreihe wegen fehlender CoolProp-Installation. Damit fehlt eine **gemessene Standardrasterempfehlung**. Die Interpolationskonvergenz nahe der Zweiphasengrenze und der kritischen Umgebung ist nicht untersucht.

**Status:** Vollständig ausführbarer Benchmark-Algorithmus entworfen und synthetisch geprüft; **quantitative Messungen/Abnahme noch offen**.
