# AP13 – Temperaturparametrisierte Spinodalen beider Phasen

## Änderungen

Neues Modul `homogeneous_pvt.py`, Funktion `solve_spinodal_branch(fluid, temperatures, branch, n_scan, residual_tolerance)`. Die bestehende Funktion `calc_spinodale_vapor_temperature` und der HEOS-Spinodaltracer bleiben erhalten. Der neue Solver ist bei `liquid` und `vapor` symmetrisch aufgebaut; der kritische Punkt wird nicht in die Daten hineingezwungen. Die Aufruf- und Ausgabeschnittstelle liefert `(n_T,3)` mit `(T,v,p)` und ein Diagnosedictionary. Er verwendet für beide Äste `T_values` direkt statt einer interpolierten Temperaturachse.

Je Temperatur: Sättigungsdichte des ausgewählten Astes bestimmen; Intervalle `(ρ_sat,v,ρc)` bzw. `(ρc,ρ_sat,l)` ableiten; `dp/dρ_m` mit `specify_phase()` auswerten; gesuchte Nullstelle mit `brentq` finden und mit `|ρc/pc * dp/dρ_m| ≤ residual_tolerance` nachprüfen. Der Dampfzweig ist ein lokales Druckmaximum bezüglich ρ, der Flüssigkeitszweig ein lokales Minimum. Zusätzliche Plausibilität: `p_spin,v > p_sat`, `p_spin,l < p_sat`, korrekte Dichteordnung. Nichtpositive Flüssigkeitsspinodaldrücke bleiben als reale EOS-Daten erhalten. Multiple Wurzeln und fehlende Auswertungen führen zu NaN statt Kurvenbrücken. Bei fehlendem Bracket werden bis zu zwei gezielte Verfeinerungen vorgenommen.

## Doppelprüfung und Tests

**Prüfung A:** Vorzeichenanalyse und Grenzintervalle; unabhängige Test-EOS mit bekannten stationären Dichten `ρ=3` (Maximum) und `ρ=7` (Minimum). **Prüfung B:** Fehlerfallprüfung (NaN-Zeilen, negativer physikalischer Druck, ungültige Temperaturen). Ein zunächst falsches Monotoniekriterium für die Flüssigkeitsspinodale wurde korrigiert: Mit steigender Temperatur kann `v_spin,l` **steigen**, während `v_spin,v` sinkt. Ein temperaturabhängiger analytischer Test bestätigt beide Trends.

`python -m unittest discover -v` hat 17 Tests ausgeführt, davon 14 erfolgreich und 3 explizit übersprungen (2 HEOS, 1 alte Testintegration wegen fehlender `Test.py`). **Nicht nachgewiesen:** beide HEOS-Zweige bei Wasser/R134a, CO₂, Argon und Residuen unter realen HEOS-Daten. Die Algorithmik ist getestet; das stoffübergreifende AP13-Akzeptanzkriterium ist weiterhin offen.
