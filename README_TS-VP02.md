# TS-VP02 – Zentraler geometrischer Isolinien-Levelgenerator

**Modul:** `isoline_levels.py` (Qt-, VTK- und CoolProp-frei) · **Tests:** `test_ts_vp02_levels.py`

## 1. Formel und Definition

Für jede Isolinienfamilie und jeden Diagrammtyp: `x_n = start · factor^n`, solange `x_n ≤ stop` (Grenztoleranz `rtol = 1e-12`), mit `0 < start ≤ stop` und `factor > 1`. Beispiel `start = 10, stop = 1000, f = 3.14`:

`10 · 31.4 · 98.596 · 309.59144 · 972.1171216` (der nächste Wert, 3052.4, liegt außerhalb; **1000 wird nicht ergänzt**).

* Die Werte entstehen aus der Potenz (`np.power`), nicht durch kumulatives Multiplizieren – identische SI-Eingaben liefern bitgleiche Ergebnisse.
* Das Raster ist geometrisch, nicht äquidistant; auf einer Log-Achse sind die Abstände konstant (`log f`). Es wird keine zusätzliche räumliche Dichteanpassung behauptet.
* Der Generator berechnet nur die Konstanten. Kurven, Abtastpunkte entlang einer Linie (`n_T`, `n_v`, `n_p`, `n_quality`, 121/81/61 Punkte …) und das 2D-Flächenraster bleiben davon getrennt.

## 2. API

| Funktion | Zweck |
|---|---|
| `generate_levels(start, stop, factor=3.14, max_levels=200)` | streng wachsendes `float64`-Array |
| `validate_levels_input(...)` | Prüfung der vier Eingaben |
| `normalise_isoline_config(raw)` / `config_signature(raw)` | kanonische Konfiguration je Familie bzw. JSON-fähige Cache-Signatur |
| `resolve_levels(config, domains)` | Levels je Familie; `domains[familie] = (lo, hi)` ist Automatikbereich und Gültigkeitsfenster; verworfene Levels im Bericht |
| `build_isolines(levels, report, makers)` | Kurven je Level aus Rechen-Callbacks, nicht berechenbare Levels werden ausgewiesen (nie erfunden) |
| `record_sampling(report, familie, kurven)` | Abtastpunktzahl je Kurve getrennt von der Levelanzahl |
| `isotherm_over_volume(...)` | gemeinsame (T,v,x)-Isotherme für u-v-T/f-v-T |

Einheiten (SI): T [K], v [m³/mol], p [Pa]. Familien: `isothermen`, `isochoren`, `isobaren`.

## 3. Validierung und Schutz (Testfälle L01–L03)

| Fall | Verhalten |
|---|---|
| `start ≤ 0`, `stop < start`, `factor ≤ 1 + 1e-9`, NaN/Inf, Nicht-Zahl | `IsolineLevelError` (Unterklasse von `ValueError`) mit lesbarer Meldung |
| `factor ≈ 1` | vor jeder Schleife über `ln(stop/start)/log1p(f-1)` abgefangen; kein Hänger |
| mehr als `max_levels` (Standard 200, Obergrenze 10 000) | `IsolineLevelError` |
| Überlauf (`stop` nahe `float.max`) / fehlender Rundungsfortschritt | `IsolineLevelError` |
| Grenzgleichheit (`x_n == stop` bis auf Rundung) | Wert bleibt erhalten |
| Familien mit vorzeichenbehafteter Konstante (h, s, u, f, g) | nicht definiert, keine stillschweigende Betrags-/Offsetumrechnung (nicht Teil der Familien) |

Konfiguration: unbekannte Schlüssel/Familien → Fehler; `start`/`stop` = `auto`, leer, `null` oder `0` → automatisch (GUI-Feld 0 = „auto“).

## 4. Abnahme

`test_ts_vp02_levels.py`: 17 Tests (Beispiel L01, exakte Formeltreue für mehrere Größenordnungen, strenge Monotonie, konstantes Verhältnis, keine Linearinterpolation/kein erzwungener Endpunkt, L02, L03, `max_levels`-Grenzwerte, Überlauf, Konfigurationsauflösung, Cache-Signatur) – alle bestanden.
