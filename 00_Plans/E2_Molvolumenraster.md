# E2 – Molvolumenraster der p-v-T-Zustandsfläche konfigurieren

**Projekt und Grundlage:** ThermoSurface. **Module:** `gui.py`, `coolprop_heos_pvt.py`, `main.py`, `config.yaml`, `diagram_cache.py`. Unabhängig von den Anzeigegrenzen aus E1. Die Abnahmetests sind vorgesehen und noch nicht durchgeführt.

## 1. Ziel und Umfang

Im Block „Berechnung: Zustandsfläche“ Start, obere Grenze und Verteilung des Molvolumens vₘ in m³/mol zugänglich machen. Linear wird Δv eingegeben, logarithmisch die Stützstellenzahl nᵥ. Die Einstellung steuert die HEOS-Auswertung jeder Flächenisotherme.

## 2. Ist-Zustand im Code

`config.yaml` enthält unter `grid_points_creation.isotherm_grid` bereits `v_start_isotherm`, `v_stop_isotherm` und `v_step_isotherm`. `ThermoSurfaceWindow.__init__()` ergänzt bei alten Konfigurationen `v_spacing="linear"` und `v_points=1000`. `setup_controls()` zeigt in dieser Gruppe jedoch nur die drei Temperaturschlüssel. `_read_grids()` behandelt Berechnungsfelder außerdem durchgehend mit `field.value()`. Eine dort ergänzte Modus-Combobox benötigt `currentData()`.

Der Rechenweg existiert bereits. `main.berechne_daten()` liest `cfg = raster["isotherm_grid"]`, erzeugt mit `build_pvt_temperature_grid()` die gemeinsamen Temperaturzeilen und ruft `multi_isotherm(**cfg, ..., T_values=T_all)` auf. `multi_isotherm()` und `calc_isotherme()` verwenden `_volume_samples()`. Jeder Eintrag der Rohfläche enthält `(T, v, p)`. Bei einer fehlgeschlagenen Druckauswertung bleibt der Volumenstützpunkt erhalten und der Druck wird `NaN`.

**Randkonvention:** `_volume_samples()` verwendet linear `np.arange(start, stop, step)`, logarithmisch `np.geomspace(start, stop, num=v_points)`. Das Temperaturraster hängt seinen Stopwert ausdrücklich an. Das lineare Volumenraster hat bereits eine andere Konvention.

## 3. Umsetzungsschritte

1. `setup_controls()` um die fünf vorhandenen Volumenschlüssel erweitern. Bei linearer Verteilung nur die Schrittweite aktivieren, bei logarithmischer nur die Punktzahl. Feldtitel für linear als „Obere Volumengrenze (exklusiv)“, für logarithmisch als „Endvolumen (inklusiv)“ formulieren. Kleine Molvolumina müssen mit ausreichenden Nachkommastellen editierbar sein.
2. `_read_grids()` und `_write_grids()` symmetrisch für Zahlenfelder und `v_spacing` behandeln. Die Werte ausschließlich nach `raster["isotherm_grid"]` schreiben. `visu_raster["isotherm_visu"]` bleibt für die separat gezeichneten Isothermen zuständig.
3. Die bestehende lineare Konvention beibehalten. Die Stützstellen liegen bei $v_{\min}+i\Delta v<v_{\max}$. Auch bei nicht teilbarem Intervall bleibt der Abstand Δv. Bei Gleitkomma-Grenzfällen sicherstellen, dass kein numerisch mit `stop` identischer Endpunkt hineinrutscht. Für logarithmische Verteilung gilt $q=(v_{\max}/v_{\min})^{1/(n_v-1)}$.
4. Die bestehenden Stoffvorgaben aus `generate_stoff_grid()` über `on_confirm_stoff()` und `on_reset_auto()` in allen fünf Feldern anzeigen. Manuelle Änderungen gelten bis zum nächsten ausdrücklichen Übernehmen der Stoffvorgaben oder Auto-Reset. Eine reine Diagrammaktivierung darf sie nicht überschreiben. Die Automatik darf die Auflösung an das Punktbudget anpassen. Eine zu große **manuelle** Eingabe mit konkreter Punktzahl zurückweisen.
5. `validate_rasters()` und die abschließende Grenze in `multi_isotherm()` auf denselben tatsächlich erzeugten Vektor beziehen. Für die Fläche zählt `len(T_all) * len(volumes)`, einschließlich zusätzlicher Sättigungszeilen und T_c. `MAX_POINTS` stammt aus der Konfiguration und beträgt im Anhang 1.000.000. `_cache_params_pvt()` enthält das Raster bereits. Prüfen, dass Modus, Grenzen und aktive Auflösung über `DiagramKey.build()` einen anderen Schlüssel erhalten, sobald sich die Rechengeometrie ändert.

## 4. Abnahmetests

Bei einem bestätigten Stoff und gültigem Temperaturraster direkt die Volumenspalte aller Zeilen des von `multi_isotherm()` gelieferten Arrays prüfen, **vor** Phasenmaske, Einfügen der Binodalränder und Rendering.

| Eingabe in m³/mol | Erwartetes Array beziehungsweise Verhalten |
|---|---|
| Linear, Start 1·10⁻⁵, Stop 5·10⁻⁵, Schritt 1·10⁻⁵ | Vier Werte: [1, 2, 3, 4]·10⁻⁵. Alle Abstände 1·10⁻⁵ |
| Linear, Start 1·10⁻⁵, Stop 6·10⁻⁵, Schritt 2·10⁻⁵ | Drei Werte: [1, 3, 5]·10⁻⁵. Letzter Stützpunkt 5·10⁻⁵, beide Abstände 2·10⁻⁵ |
| Logarithmisch, Start 1·10⁻⁵, Stop 1·10⁻³, nᵥ = 3 | [10⁻⁵, 10⁻⁴, 10⁻³]. Beide Endwerte enthalten, beide Quotienten 10 |
| Start ≥ Stop, nᵥ = 1 oder tatsächliches Raster oberhalb des Punktbudgets | Eingabe wird vor der Flächenberechnung abgelehnt. Das zuletzt gültige Diagramm bleibt verwendbar |
| Manuelles Raster berechnet, danach anderen Stoff bestätigen und Auto-Reset ausführen | Felder und nächstes Flächenarray entsprechen den Vorgaben des neuen Stoffs. Eine Änderung des reinen Ansichtsausschnitts verändert das Array nicht |

## 5. Abhängigkeiten und offene Entscheidungen

Empfohlen ist das halboffene lineare Raster, weil es die bestehende Volumenberechnung reproduziert. Ein obligatorischer Endpunkt wäre eine zusätzliche fachliche Änderung: Bei nicht teilbarem Intervall müsste entweder der letzte Abstand verkürzt oder die Schrittweite aller Intervalle angepasst werden.
