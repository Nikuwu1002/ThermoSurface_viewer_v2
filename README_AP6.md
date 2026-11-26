# AP6 – Externe thermodynamische Referenzvalidierung

**Ziel:** Ausgewählte HEOS-Zustandsgrößen, Sättigungsdrücke, kritische Werte und die Clapeyron-Gleichung quantitativ und mit Quellennachweis prüfen.

## Datenquelle und Berechnung

`ap6_reference.py` enthält originalgetreu in SI-Einheiten übertragene IAPWS-IF97-Testwerte aus der offiziellen [IAPWS R7-97(2012)](https://iapws.org/documents/release/IF97-Rev): Tab. 5 (flüssige Einphasenpunkte, `300 K/3 MPa`, `300 K/80 MPa`, `500 K/3 MPa`), Tab. 15 (Dampf-Einphasenpunkte, `700 K/3,5 kPa`, `700 K/30 MPa`) und Tab. 35 (Sättigung bei `300`, `500`, `600 K`). Dazu kritische Wasserkonstanten `Tc=647,096 K`, `pc=22,064 MPa`, `rho_c=322 kg/m³`.

Ein separat initialisierter `AbstractState('HEOS','Water')` bewertet für identische `(T,p)`-Zustände `v_mass`, `h_mass`, `s_mass`, `u_mass`; für die Sättigung wird `QT_INPUTS` verwendet. Ausgabe: `results/ap6/ap6_references.csv` mit Quelle, Einheiten, Soll, Ist, absolutem/relativem Fehler und Pass/Fail. Für Tabellengrößen gilt als **eigene vorläufige Prüfung** maximal 0,5 % relative Abweichung; kritische Werte werden mit `1e-5` relativer Toleranz geprüft. Diese Grenzen sind keine experimentelle Fehlerangabe. Ein fehlgeschlagener Vergleich beendet das Skript mit Exit-Code 1.

Für `T=373,15`, `500` und `600 K` wird der numerische Sättigungsdruckgradient (zentral, Schritt 0,05 K) gegen `(h_g-h_l)/(T(v_g,mass-v_l,mass))` verglichen. Beide Seiten stammen hier aus HEOS: eine **Konsistenzprüfung**, ausdrücklich kein weiterer unabhängiger Datensatz. Das Vorzeichen und die massenspezifischen Einheiten werden berücksichtigt. Auswertungsstatus und Anzahl fehlgeschlagener Vergleiche stehen in `ap6_summary.json`.

## Aufruf

```powershell
python ap6_reference.py --output results/ap6
python -m unittest test_ap5_ap8.AP6ReferenceTests -v
```

## Prüfung und Status

**Bestanden:** Zwei rein numerische Tests zur unverfälschten Referenztabellendatenlage und relativen Fehlerlogik einschließlich erzwungenem Fehlfall/NaN. Die Tabellenwerte wurden vor Aufnahme gegen das offizielle IAPWS-Dokument (Tabellen 5, 15 und 35) kontrolliert.

**Nicht ausgeführt:** Die eigentliche quantitative HEOS-vs.-IF97-Vergleichstabelle und Clapeyron-Zahlen, da CoolProp nicht installiert und nicht installierbar war. Es werden keine beispielhaft erfundenen Ist-Werte angegeben. Außerdem ist IF97 eine andere, aber auf IAPWS-95 abgestimmte Referenz**formulierung**, keine unabhängige experimentelle Datensammlung. Andere Fluide sind außerhalb der Wasserreferenzmatrix.

**Status:** Vollständiges validierendes Auswerteskript implementiert, Teiltests erfolgreich; **Abnahmekriterium (quantitative Vergleichsergebnisse) noch offen**.
