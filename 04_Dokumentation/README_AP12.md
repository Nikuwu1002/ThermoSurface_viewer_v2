# AP12 – Thermodynamische Zielgeometrie und Machbarkeit

## Physikalischer Vertrag

Für gewöhnliche Flüssigkeit-Dampf-Koexistenz bei `T<Tc` gilt in der **homogenen** EOS die Volumenordnung `v_sat,l < v_spin,l < v_spin,v < v_sat,v`. Damit bestehen drei zusätzliche, vom Gleichgewicht strikt zu trennende Bereiche: metastabile Flüssigkeit `[v_sat,l,v_spin,l]`, mechanisch instabile homogene Fortsetzung `[v_spin,l,v_spin,v]` und metastabiler Dampf `[v_spin,v,v_sat,v]`. Positive isotherme Kompressibilität bedeutet `(dp/dv)_T<0`, also `(dp/dρ_m)_T>0`. Die instabile Fortsetzung hat umgekehrtes Vorzeichen. Beide Spinodalen erfüllen `(dp/dρ_m)_T=0`.

Der kritische Punkt `(Tc,vc,pc)` ist singulär. Er ist ein physikalischer Grenzwert, darf aber nur nach nachgewiesener Konvergenz an neue homogene Teilflächen angeschlossen werden. Die bestehenden Binodalen besitzen bereits den exakten kritischen Endpunkt. Ein künstlicher zusätzlicher Spinodalendpunkt ist ausdrücklich nicht Bestandteil der Umsetzung.

Die Oberfläche der Gleichgewichtsmischung bleibt unverändert mit `p=p_sat(T)` bestehen. Die Drücke der neuen Flächen werden aus der homogenen Zustandsgleichung `p_hom(T,v)` berechnet, nicht zwischen Spinodaldruck und Sättigungsdruck interpoliert. Das reguläre Gleichgewichtsnetz und die homogenen Fortsetzungen werden separat gerendert.

## Gewählte Implementierungsstrategie

**Strategie A (HEOS mit Phasenvorgabe)** ist implementiert: `AbstractState('HEOS', fluid)`, `specify_phase(iphase_liquid/iphase_gas)`, `update(DmolarT_INPUTS,1/v,T)` und direkte Ableitungen über `first_partial_deriv`. Die Phasenvorgabe wird in `finally` zurückgesetzt. Ungültige Zustände werden als NaN/fehlende Zeilen markiert. Eine zusätzliche vdW-EOS wurde bewusst **nicht** mit HEOS-Binodalen vermischt. **Strategie B** wurde nicht implementiert.

## Qualität und offener Entscheidungspunkt

**Prüfung A:** Dimensionsanalyse, Ableitungsvorzeichen, Dichteintervalle und Trennung der physikalischen Datensätze. **Prüfung B:** Mit einer eigenständigen analytischen, dichtekubischen Test-EOS wurden alle drei Bereiche auf Existenz und konsistente Druck-/Ableitungswerte geprüft (`test_ap13_ap18.py`). Das ist ausdrücklich **kein** HEOS-Machbarkeitsnachweis. CoolProp ist in der Ausführungsumgebung nicht installiert und lässt sich hier nicht beziehen. Der vorgeschriebene Versuch mit Wasser und einem zweiten realen Fluid ist daher **nicht durchgeführt**.

**Freigabestatus:** nur algorithmische Vorprüfung; die Freigabe für reale Stoffe und den kritischen Bereich ist bis `validate_ap17_ap18.py` auf dem Zielrechner offen. Die späteren Module sind deshalb fehlertolerant implementiert, aber **nicht als stoffübergreifend validiert** auszuweisen.
