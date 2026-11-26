# AP7 – Vergleich HEOS und eigene van-der-Waals-Zustandsgleichung

**Ziel:** Eigenständige begrenzte vdW-EOS, identische homogene Einphasenzustände und reproduzierbare Druckfehler; keine Modellumschaltung in der GUI, keine Gleichsetzung von vdW-Schleife und Phasengleichgewicht.

## Implementierung

`ap7_vdw_compare.py` implementiert `vdw_parameters(Tc,pc)`, `vdw_pressure(T,v,a,b)`, `vdw_derivative(T,v,a,b)`, `vdw_critical(Tc,pc)`, `compare(fluid,n_v)` und `write_results(...)`. Molares SI-System: `R=8,31446261815324 J/(mol K)`, `a=27R²Tc²/(64pc)`, `b=RTc/(8pc)`, `p_vdW=RT/(v-b)-a/v²`. Bereiche `T>0`, `v>b`, `Tc>0`, `pc>0` werden explizit durch `ValueError` abgesichert. Für vdW gilt analytisch `vc=3b` und `Zc=0,375`.

Bei `Tr=0,85` und `0,95` beginnt das Volumenraster bei `1,2*v_g,HEOS(T)` und endet bei `5*v_g,HEOS(T)`; die Grenze wird durch HEOS-Sättigungsdampfwerte ermittelt. Der Vergleich findet nur auf der stabilen Gasseite statt. Bei `Tr=1,10` werden gasähnliche überkritische Zustände zwischen etwa `1,5*vc,HEOS` und `8*vc,HEOS` gewählt. Für jedes identische `(T,v)` werden beide Drücke ausgewertet und gespeichert: `results/ap7/ap7_vdw_vs_heos.csv`, `ap7_summary.json`, optional drei PNG-Diagramme (Matplotlib).

Eine Maxwell-Konstruktion oder die Berechnung von vdW-Sättigungsdrücken ist **absichtlich nicht Teil dieses eingeschränkten Vergleichs**. Das echte kritische HEOS-Volumen ist nicht mit `3b` identisch. Abweichende kritische Kompressibilitätsfaktoren sind bereits analytisch nachweisbar.

## Tests und Status

```powershell
python -m unittest test_ap5_ap8.AP7AnalyticalTests -v
python ap7_vdw_compare.py --fluid Water --n-v 16 --output results/ap7
```

**Bestanden (2 Tests):** Kritische vdW-Identitäten, normierte Ableitungsbedingung, Fehlerfälle `v<=b` und ungültiges `T`. Eine anfänglich ungeeignete absolute Nulltoleranz für die Differenz großer Ableitungsterme wurde durch eine dimensionslose Prüfung ersetzt. Reiner analytischer Wasservergleich anhand veröffentlichter IAPWS-Konstanten ergibt `vc,vdW≈9,1443e-5 m³/mol`, `vc,Wasser≈5,5948e-5 m³/mol` und `Zc,Wasser≈0,22944`; dies sind keine HEOS-Durchlaufmessungen.

**Offen:** Die eigentlichen punktweisen Druckabweichungen und der automatische Export einer HEOS-vdW-Vergleichsreihe konnten wegen fehlendem CoolProp nicht ausgeführt werden.

**Status:** Vergleichscode implementiert, analytische Grundlagen getestet; **quantitative Modellvergleichsmatrix auf Zielrechner erforderlich**.
