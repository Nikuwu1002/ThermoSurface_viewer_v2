# AP14 – Gemeinsame Randkurven und Temperaturstützstellen

## Änderungen

`main.berechne_daten` verwendet die durch AP2 bereits vorhandene Temperaturbasis `T_sat = binodal_liquid[:-1,0]` für **beide** neuen Spinodalsuchen. Dadurch besitzen `homogeneous_spinodal_liquid` und `homogeneous_spinodal_vapor` exakt dieselben Temperaturstützstellen wie die Unterkritik-Binodalen. Fehlende Wurzeln verbleiben als `(T,NaN,NaN)` im entsprechenden Index; es erfolgt weder eine Interpolation zwischen erfolgreichen Zeilen noch eine künstliche Ergänzung fehlender Spinodalpunkte.

Das Teilnetz einer metastabilen Flüssigkeit nutzt direkt `binodal_liquid[:-1]` und `homogeneous_spinodal_liquid` als Spalte 0 bzw. letzte Spalte. Für den Dampf werden `homogeneous_spinodal_vapor` und `binodal_vapor[:-1]` entsprechend übernommen. Für die optionale instabile Fläche bilden die beiden neuen Spinodaläste den Rand. Der kritische Punkt bleibt als Singularität im alten Gleichgewichtsnetz angebunden, nicht jedoch in unvollständige homogene Datensätze hineininterpoliert.

## Prüfungen

**Prüfung A:** Form `(n_T,n_v,3)`, monotone Volumenkoordinaten, bitweise identische Grenzen. **Prüfung B:** `validate_homogeneous_boundaries` vergleicht tatsächliche gespeicherte Endspalten mit den Randquellen per `np.array_equal`; das synthetische Testmodell erfüllt diese Gleichheit. Der Renderpfad verwendet dieselbe deterministische logarithmische Transformation für Linien und Netzknoten. `validate_ap17_ap18.py` prüft die transformierten Randfehler bei echten HEOS-Daten auf dem Zielrechner.

**Einschränkung:** Bei einzelnen fehlenden Zwischenzeilen sind Netzstreifen absichtlich offen. Der exakte geometrische Anschluss einer homogenen Spinodalfläche bei `Tc` ist ohne nachgewiesene HEOS-Grenzwertkonvergenz **nicht** zugesichert. Realer PyVista-/HEOS-Test steht aus.
