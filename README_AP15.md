# AP15 – Rekonstruktion zusätzlicher homogenen Zustandsflächen

## Änderungen

`build_homogeneous_surfaces()` in `homogeneous_pvt.py` erzeugt die zusätzlichen Ergebnisarrays `metastable_liquid_surface`, `metastable_vapor_surface` und `unstable_surface`. Diese werden in `main.berechne_daten()` zum bestehenden Ergebnis-Dictionary hinzugefügt; die vorhandenen Schlüssel und die Gleichgewichtsberechnung bleiben erhalten. Die neuen Flächen werden ausschließlich aus homogen vorgeschriebenen HEOS-Zuständen berechnet: `p=state.p()` nach `DmolarT_INPUTS`, **keine Druckinterpolation**. Lediglich die zu berechnenden **Volumenstützstellen** werden logarithmisch zwischen den Randvolumina verteilt. Unter den neuen Arrays sind Fehler als NaN-Zeilen sichtbar.

Die beiden metastabilen Bereiche müssen `(dp/dρ_m)_T ≥ -ε` erfüllen, die optionale instabile Fortsetzung `(dp/dρ_m)_T ≤ ε`. Es werden weder Sättigungs-Flashdrücke im Innenraum eingesetzt noch negative Drücke auf Null gesetzt. Jeder erste und letzte Punkt einer Zeile wird exakt aus den in AP14 gespeicherten Randkurven übernommen; ein erneuter Homogen-EOS-Vergleich kontrolliert die Randdrücke. Bei unstimmigen Rändern wird die gesamte Temperaturzeile ausgelassen.

`config.yaml` steuert über `diagram_settings.pvt.include_homogeneous`, `homogeneous_points` (Standard 20) und `include_unstable` (Standard `false`) die neuen Berechnungen. Das Punktlimit von einer Million wird bei der neuen Berechnung separat geprüft.

## Doppelprüfung

**A:** Ableitungsvorzeichen und homogene Druckberechnung anhand analytischer Test-EOS. **B:** alle Innenpunkte gegen die analytische EOS und alle Randpunkte bitweise gegen ihre Quellarrays verglichen. Fehlerhafte Grenzdrücke führen nachweislich zu abgewiesenen Zeilen. Das Testverfahren wurde ausgeführt und bestanden.

**Noch offen:** tatsächliche Auswertbarkeit aller entsprechenden HEOS-Zustände und damit Freigabe bei realen Fluiden. Die instabile Fläche bleibt standardmäßig deaktiviert, da eine stoffabhängige Prüfung unverzichtbar ist. Auch positive / negative Druckabschnitte sind getrennt zu bewerten (nichtpositive Drücke sind in `log(p)` nicht zeigbar).
