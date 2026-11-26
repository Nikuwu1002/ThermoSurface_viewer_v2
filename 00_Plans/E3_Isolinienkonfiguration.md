# E3 – Isolinienkonfiguration zusammenführen und GUI bereinigen

**Projekt und Grundlage:** ThermoSurface. **Module:** `gui.py`, `isoline_levels.py`, `main.py`, `calculation_manager.py`, die sieben diagrammspezifischen HEOS-Module, `config.yaml`. Zwischenstand A vor B abnehmen. Gemeinsame Rasterfelder mit E2 abstimmen. Die Abnahmetests sind vorgesehen und noch nicht durchgeführt.

## 1. Ziel und Umfang

**Zwischenstand A:** Im Block „Isolinien – Linienwerte und Auflösung“ je Diagramm und Familie die konstanten Linienwerte und die Abtastung entlang einer Linie unabhängig einstellen. Die alten GUI-Blöcke bleiben für den Vergleich zunächst vorhanden.

**Zwischenstand B:** Nach erfolgreichem Vergleich die alten Blöcke entfernen. Weiterhin wirksame Phasen-, Solver- und Abtastparameter müssen auch beim Laden einer bestehenden `config.yaml` erhalten bleiben. Standard einsmtellung soll hierbei default geometrisch sein

## 2. Ist-Zustand im Code

Die drei alten Liniengruppen gehören zu `visu_points_creation`. Die zentralen Vorgaben liegen unter `diagram_settings.<id>.isolines`.

| Bisheriges GUI-Feld oder Vorgabe | Konfigurationsschlüssel | Tatsächlicher Verbraucher | Künftiger Ersatz oder Standard |
|---|---|---|---|
| Zentraler Start, Stop, Faktor je Familie | `isolines.<familie>.start/stop/factor` | `isoline_config_si()`, `normalise_isoline_config()`, `resolve_levels()`, `generate_levels()` | Beibehalten, um Modus und additiven Abstand ergänzen. Bestehende Konfiguration ohne Modus bleibt geometrisch |
| Isothermen: T-Bereich und deaktivierte Anzahl | `isotherm_visu.T_start_isotherm`, `T_stop_isotherm`, `T_step_isotherm` und GUI-Hilfswert `n_isotherms` | `main.berechne_daten()` nutzt die Grenzen als Leveldomäne. Die Anzahl wirkt im alten Aufrufpfad ohne explizite Levels | Auto-Domäne aus den bestehenden Grenzen erhalten. `n_isotherms` aus der Oberfläche entfernen |
| Isothermen: Volumenbereich, Schritt, Verteilung, Punktzahl | `isotherm_visu.v_start_isotherm`, `v_stop_isotherm`, `v_step_isotherm`, `v_spacing`, `v_points` | `multi_isotherm()` → `calc_isotherme()` → `_volume_samples()` | Zentrale Linienabtastung über das Molvolumen. Bestehende Grenzen und Verteilung übernehmen |
| Isochoren: Volumenbereich und Volumenfolge | `isochore_visu.v_start_isochore`, `v_stop_isochore`, `v_step_isochore`, `v_spacing`, `v_points` | `resolve_levels()` verwendet die Domäne. `multi_isochore()` erzeugt aus den alten Folgeneinstellungen nur dann Werte, wenn `v_values` fehlt | Zentraler Levelmodus für das Molvolumen. Alte Werte im Kompatibilitätspfad erhalten |
| Isochoren: Temperaturbereich und Schritt | `isochore_visu.T_start_isochore`, `T_stop_isochore`, `T_step_isochore` | `calc_isochore()` durchläuft die Temperaturen einer einzelnen Isochore | Zentrale Linienabtastung über T |
| Isobaren: Druckbereich und Druckabstand | `isobare_visu.p_start_isobare`, `p_stop_isobare`, `p_step_isobare` | Domäne für `resolve_levels()`. `multi_isobare()` verwendet den alten Druckabstand nur ohne `p_values` | Zentraler Levelmodus für p, intern Pa |
| Isobaren: Temperaturbereich und Schritt | `isobare_visu.T_start_isobare`, `T_stop_isobare`, `T_step_isobare` | `calc_isobare()` tastet T ab und kann zusätzlich 30 Qualitätsstützstellen am Phasenübergang einfügen | Zentrale T-Abtastung, zusätzliche Phasenpunkte separat ausweisen |
| Binodale/Zweiphasenfläche: T-Bereich und Schritt | `binodale_visu.T_start_sat`, `T_stop_sat`, `T_step_sat` | `build_pvt_temperature_grid()` erzeugt daraus zusätzliche Flächenzeilen. `multi_binodale(..., T_values=T_sat)` verwendet das gemeinsame Raster | Aus wirksamer Konfiguration bzw. Stoffautomatik beziehen. Auch bei ausgeblendeter Binodale erhalten |
| Spinodale: Temperaturbereich | `spinodale_visu.T_start_spinodal`, `T_stop_spinodal` | `berechne_daten()` begrenzt die Zweige und schneidet das Fenster mit `raster["spinodale_grid"]` für den Dampfsolver | Stoffbezogene Grenzen weiterführen. `n_T`, `n_scan` und `residual_tolerance` des Solvers separat erhalten |
| Zweiphasenauflösung und Druckfenster | `two_phase_visu.Q_points`, `p_start_twophase`, `p_stop_twophase` | `calc_two_phase_surface()`, `clip_two_phase_surface()` und `_berechne_zusatzdiagramm()` | `Q_points` als gemeinsame Konfigurationsquelle behalten. Druckfenster nur für p-v-T, aus wirksamer Konfiguration bzw. Stoffautomatik |
| Zusatzdiagramme: implizite Linienauflösung | Teilweise `diagram_settings.<id>.n_p/n_v/n_T_sub/n_T_super`, teilweise lokale Konstanten | Diagrammspezifische Linienerzeuger, siehe unten | Eigene optionale Auflösung je Familie. Ohne neue Vorgabe das jeweilige bisherige Verhalten reproduzieren |

Die Zusatzdiagramme lesen die Linienabtastung der alten p-v-T-Blöcke nicht einheitlich. Für Zwischenstand A sind folgende Verbraucher anzuschließen:

| Diagramm | Isotherme | Isochore | Isobare |
|---|---|---|---|
| p-v-T | `calc_isotherme()`, Volumenvektor | `calc_isochore()`, T-Vektor | `calc_isobare()`, T-Vektor und zusätzliche Q-Punkte |
| p-h-T | `coolprop_heos_pht._isolinien()` → `calc_isotherme()`, derzeit `n_p=41` und Q-Abtastung | `calc_isochore()`, `n_T=121`, anschließend Druckfilter | `calc_isobare()`, `n_T=81` und Q-Abtastung |
| p-s-T | `coolprop_heos_pst._isolinien()` → `calc_isotherme()`, `n_p=41` und Q-Abtastung | `calc_isochore()`, `n_T=121`, Druckfenster | `calc_isobare()`, `n_T=81` und Q-Abtastung |
| h-s-p | `berechne_hsp_daten()` übernimmt die p-s-T-Isothermen | Übernimmt die p-s-T-Isochoren | Übernimmt die p-s-T-Isobaren. `_convert_pst_curve()` transformiert die Kurven. Eigene h-s-p-Vorgaben bereits an den p-s-T-Aufruf weitergeben |
| u-v-T | `_calculate_isolines()` → `isotherm_over_volume()`, bisher `n_v` je Ast und Q-Abtastung | Lokale `isochore()` verwendet vorhandene Temperaturstützstellen | `_isobare_uvt()` verwendet p-h-T-Isobare mit `n_T=61` und Q-Abtastung |
| g-p-T | `_calculate_isolines()` → `_isotherme()` mit `n_p` | Lokale `isochore()` verwendet vorhandene Temperaturen und Druckfilter | `_isobare()` mit 81 Temperaturstützstellen |
| f-v-T | `_calculate_isolines()` → `isotherm_over_volume()`, bisher `n_v` und Q-Abtastung | Lokale `isochore()` verwendet vorhandene Temperaturen | `_isobare()` mit 61 Temperaturstützstellen und Q-Abtastung, anschließend Volumenfilter |

## 3. Umsetzungsschritte

1. **A – Konfiguration erweitern.** In `isoline_levels.py` je Familie die neuen Schlüssel `mode`, `step` und einen separaten Abtastabschnitt vorsehen. Aktiv sind für die lineare Folge Start, Stop und positiver additiver Abstand, für die geometrische Folge Start, Stop und Faktor > 1. Endliche positive physikalische Werte verlangen, `auto` weiterhin über die Diagrammdomäne auflösen. Gleiche explizite Grenzen können genau eine Linie bezeichnen. Bei vertauschten Grenzen abbrechen.
2. **A – Endwert eindeutig behandeln.** Linear $x_k=x_0+k\Delta x$, geometrisch $x_k=x_0 f^k$. Den Endwert nur aufnehmen, wenn die Folge ihn innerhalb der numerischen Toleranz trifft. `max_levels`, Domänenfilter sowie Meldungen für verworfene oder leere Linien erhalten. Bei einem Druckabstand dieselbe Umrechnung von der GUI-Einheit nach Pa wie bei den Druckgrenzen anwenden.
3. **A – GUI und Datenfluss gemeinsam umstellen.** `_build_isoline_group()`, `_isoline_gui_values()`, `_refresh_isoline_panel()`, `_on_isoline_field_changed()`, `isoline_config_si()` und `_validate_isolines()` erweitern. `_isoline_target()` bestimmt das bearbeitete Diagramm. Die effektive Konfiguration über `DiagramTask.isolines`, `main._effektive_isolinien()` und die jeweiligen Rechenfunktionen weiterreichen. `config_signature()` und die Cacheparameter müssen alle rechnerisch wirksamen Felder berücksichtigen.
4. **A – Linienauflösung direkt an den Verbraucher anschließen.** Drei Isothermen bei 300, 320 und 340 K mit je 100 Volumenstützstellen bedeuten drei HEOS-Kurven. Der Abstand 20 K steuert die Anzahl der Linien, die 100 steuert den Volumenvektor je Linie. Bei phasenweise erzeugten Kurven die Eingabe präzise als „Stützstellen je Ast“, „Temperaturstützstellen“ oder „Q-Stützstellen“ beschriften. `record_sampling()` erfasst die ausgegebene Arraylänge, einschließlich möglicher `NaN`-Zeilen. Zusätzlich die Zahl endlicher Punkte ausweisen, wenn der Bericht die tatsächlich darstellbare Auflösung angibt.
5. **A – Alte Abtastung reproduzierbar halten.** Bestehende Schrittweiten und Verteilungen beim Laden unverändert auf die neuen Felder abbilden. Eine alte halboffene Schrittfolge darf beim Vergleich nicht unbemerkt zu einer endpunktinklusiven `linspace`-Folge werden. Neue explizite Punktzahlen können einen eigenen, dokumentierten Modus erhalten. Während A pro Parameter genau eine wirksame Eingabe zulassen, damit die beiden Panels keine widersprüchlichen Werte liefern.
6. **B – Erst nach Referenzvergleich Widgets entfernen.** In `setup_controls()` die sechs alten Darstellungsgruppen entfernen. `_read_grids()` darf danach keine verborgenen Werte aus fehlenden Widgets rekonstruieren oder `n_isotherms` zwingend voraussetzen. Normalisierte Basiswerte und `generate_stoff_grid()` versorgen weiterhin die in der Zuordnungstabelle genannten Verbraucher. Insbesondere muss `two_phase_visu.Q_points` auch bei ausschließlicher Auswahl eines Zusatzdiagramms verfügbar bleiben.
7. **B – Ebenensteuerung erhalten.** Die Checkboxen im Block „Darstellungsebenen“, `_on_layer_toggled()` und `DiagramViewport.apply_layers()` weiterhin für die Sichtbarkeit verwenden. Ein Klick auf „Binodale“ darf insbesondere keine Änderung des gemeinsamen p-v-T-Temperaturrasters auslösen.

## 4. Abnahmetests

| Zwischenstand | Ausgangszustand und Bedienhandlung | Soll-Ergebnis |
|---|---|---|
| A | Geeignete p-v-T-Domäne, linear 300 bis 340 K in 20 K, je 100 Volumenpunkte. Danach nur auf 150 Punkte erhöhen | In beiden Läufen exakt 300, 320, 340 K. Die drei Rohkurven haben zunächst je 100, danach je 150 Zeilen. Das Flächenraster bleibt gleich |
| A | Geometrische Folge 100 bis 450 mit Faktor 2, innerhalb gültiger Domäne. Danach Stop auf 400 setzen | Beide Läufe liefern 100, 200, 400. Ein Faktor 1 oder linearer Abstand 0 wird zurückgewiesen |
| A | In der Mehrfachansicht nur die Isobarenauflösung von h-s-p ändern | Die an p-s-T delegierte h-s-p-Linienberechnung erhält den neuen Wert. Ein separat geöffnetes p-s-T-Diagramm behält seine eigene Konfiguration |
| A → B | Vor Entfernen der Blöcke für alle sieben Diagramme wirksame Parameter und Ergebnisarrays sichern. Dieselben Eingaben nach B erneut berechnen | Gleiche Level und Kurvenreihenfolge, gleiche Arrayformen, gleiche `NaN`-Positionen und numerisch gleiche endliche Werte innerhalb vorab festgelegter Toleranz. Binodalen, verfügbare Spinodalen und Zweiphasenflächen ebenfalls vergleichen |
| B | Alte Konfiguration mit abweichendem `Q_points` laden, nur ein Zusatzdiagramm berechnen. Anschließend Ebenen aus- und einblenden | Alte Q-Auflösung bleibt wirksam. Umschalten der Sichtbarkeit verändert weder Arrays noch Cache-Schlüssel |

## 5. Abhängigkeiten und offene Entscheidungen

Die wichtigste Designentscheidung ist die Bedeutung von „100 Stützstellen“. Bei einer einfachen p-v-T-Isotherme kann das exakt die Länge der Rohkurve sein. Bei einer subkritischen Zusatzdiagramm-Isotherme existieren Flüssigkeitsast, Dampfast und Qualitätsstützstellen. Empfohlen ist dort eine explizite Angabe **je Ast** mit separater Q-Auflösung. Damit bleiben Sättigungsrandpunkte und bestehende Phasenübergänge reproduzierbar.
