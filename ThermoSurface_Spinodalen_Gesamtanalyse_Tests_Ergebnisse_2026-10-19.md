# ThermoSurface: Vollständige Analyse der HEOS-Spinodalen, Testskripte und Messergebnisse

**Stand:** 19. Oktober 2026  
**Zweck:** Vollständige, eigenständig lesbare Übergabedokumentation dieses Chatverlaufs als Daten- und Entscheidungsgrundlage für ein **später gesondert zu definierendes Arbeitspaket**.  
**Projekt:** ThermoSurface, interaktive thermodynamische 3D-Zustandsdiagramme mit CoolProp HEOS, PyVista, PySide6 und Python.  
**Geltungsbereich:** p-v-T-Spinodalen, insbesondere der verkürzte dampfseitige Spinodalast von Wasser.  
**Wichtig:** Diese Dokumentation ist **keine Aufforderung zur sofortigen Implementierung**. In diesem Gespräch wurden keine Projekt-Quelldateien durch den Assistenten verändert. Die Testskripte wurden dem Nutzer zur lokalen Ausführung vorgeschlagen; die hier aufgeführten tatsächlich ausgeführten Tests sind durch die vom Nutzer zurückgesendeten Konsolen-Screenshots belegt.

## 0. Leseschlüssel: belegte Fakten vs. Vorschläge

Die folgenden Bezeichnungen verhindern, dass bei der späteren Arbeitspaketdefinition Vermutungen als Messergebnisse übernommen werden.

- **[CODE]**: In den am 19.10.2026 bereitgestellten Projektdateien nachgelesen.
- **[MESSUNG]**: Lokales Ergebnis, das der Nutzer als Screenshot mitgeteilt hat; Ablesegenauigkeit entsprechend Konsolenausgabe.
- **[AUSWERTUNG]**: Aus Code und Messung hergeleitete fachliche Interpretation; keine unabhängige neue HEOS-Messung.
- **[VORSCHLAG]**: In diesem Chat empfohlene, **nicht als umgesetzt oder getestet bestätigte** Codeänderung.
- **[OFFEN]**: Noch nicht experimentell oder numerisch nachgewiesen.
- **[FRÜHERE ANGABE]**: Im Gespräch genannte bibliotheksinterne Informationen, für die in diesem Dokument keine erneute eigenständige Überprüfung vorgenommen wurde.

**Provenienz:** Die Messergebnisse stammen aus einer lokalen Windows-PowerShell/VS-Code-Umgebung des Nutzers unter Python 3.13. Konkrete installierte Versionen von CoolProp, NumPy, SciPy, PyVista und den übrigen Abhängigkeiten wurden nicht mitgeteilt. Die vier Bildquellen in diesem Chat heißen `image(20261019-074618).png`, `image(20261019-075314).png`, `image(20261019-075930).png` und `image(20261019-080315).png`. Die relevanten Messinhalte sind unten vollständig als Text und Tabellen übertragen, damit diese Markdown-Datei auch ohne die Bilder nutzbar bleibt. Ein gesondertes ZIP mit den Original-Screenshots kann zusätzlich archiviert werden.

---

## 1. Fragestellung und Verlauf des gesamten Chats

### 1.1 Ausgangsfrage: Funktionsweise von CoolProp HEOS

Der Nutzer fragte, **wie das HEOS-Backend Spinodalen erzeugt und welche Parameter definierbar sind**. Geklärt wurde:

- Die Spinodale ist eine mechanische Stabilitätsgrenze eines **homogenen** Zustands, nicht die Binodale bzw. Sättigungslinie.
- CoolProp HEOS verwendet eine stoffspezifische Helmholtz-Energie-Zustandsgleichung und stellt `build_spinodal()`/`get_spinodal_data()` bereit.
- Interne Kurventracer-Parameter sind über die reguläre Python-Funktion `build_spinodal()` nicht unmittelbar einstellbar.
- ThermoSurface kann zwar einen Temperaturbereich zum **nachträglichen Filtern und Anzeigen** festlegen, dies veranlasst CoolProp bisher aber nicht, beliebige neue Spinodalpunkte in diesem Bereich zu berechnen.
- Eine eigene temperaturparametrisierte Nullstellensuche könnte die HEOS-Ableitung verwenden, ohne eine neue Zustandsgleichung zu implementieren.

### 1.2 Beobachtung: ungleich lange Spinodaläste

Der Nutzer zeigte eine Abbildung der Spinodalen im reduzierten logarithmischen p-v-Plot:

- Violetter Ast: Flüssigkeitsseite; reicht in der Darstellung zu deutlich niedrigeren Drücken.
- Orangefarbener Ast: Dampfseite; wirkt stark verkürzt und endet ungefähr bei `log10(v/vc) = 0,4` auf der gezeigten Kurve.
- Schwarzer Punkt: kritischer Punkt/markierter gemeinsamer Endpunkt.
- Abbildungsachsen: horizontal `log10(v_r)` bzw. `log10(v/vc)`; vertikal `log10(p_r)` bzw. `log10(p/pc)`.
- Nutzerziel: Der Dampfast soll nicht frühzeitig aufhören, sondern mindestens über den gewünschten thermodynamischen Temperaturbereich untersucht und, soweit gültig, dargestellt werden.

**Wichtige fachliche Einschränkung:** Gleiche **Temperaturabdeckung** bedeutet nicht gleiche geometrische Kurvenlänge, gleichen Druckbereich oder Ausdehnung bis zu den äußeren Volumengrenzen der Zustandsoberfläche. Eine Spinodale wird durch die Ableitungs-Nullstelle bestimmt, nicht durch die Mesh-Grenze.

### 1.3 Diagnose und wiederholte Tests

Danach wurden nacheinander drei Arten lokaler Tests vorgeschlagen und vom Nutzer ausgeführt:

1. **Rohdaten/Filter/AP3-Status prüfen**: Wie viele Spinodalpunkte liefert die bisherige Methode, welcher Temperaturbereich ist wirklich abgedeckt, liegen nichtpositive Drücke vor?
2. **Ableitungs-Scan im Dampf-Dichteintervall**: Existieren zwischen 580 und 640 K bereits HEOS-Nullstellen, obwohl `build_spinodal()` dort keine Dampfpunkte liefert?
3. **Brent-Nullstellensuche und Klassifizierung**: Welche von zwei gefundenen Nullstellen gehört wahrscheinlich zur gewünschten dampfseitigen Stabilitätsgrenze? Wie groß sind Druck, Dichte, Volumen und numerische Residuen?

Die folgenden Abschnitte enthalten **alle im Chat gezeigten Skripte, Einstellungen und gemeldeten Ergebnisse**.

---

## 2. Thermodynamische und numerische Grundlagen

### 2.1 Definition der Spinodale

Für eine homogene Phase lautet die mechanische Spinodalbedingung

\[
\left(\frac{\partial p}{\partial v_m}\right)_T=0
\quad\Longleftrightarrow\quad
\left(\frac{\partial p}{\partial \rho_m}\right)_T=0,
\qquad v_m=\frac1{\rho_m}.
\]

Dabei ist \(\rho_m\) die **molare Dichte** in `mol/m³`, \(v_m\) das **molare Volumen** in `m³/mol`. Für den kritischen Punkt einer gewöhnlichen Flüssig-Dampf-Kritikalität gilt zusätzlich die verschwindende zweite Ableitung nach dem Volumen.

Die für diese Untersuchung benötigte HEOS-Ableitung wird mit

```python
state.first_partial_deriv(CP.iP, CP.iDmolar, CP.iT)
```

abgefragt, nachdem ein homogener Zustand durch `state.specify_phase(...)` und `state.update(CP.DmolarT_INPUTS, rho_molar, T)` gesetzt wurde.

### 2.2 Was `build_spinodal()` liefert

```python
from CoolProp import AbstractState
state = AbstractState("HEOS", "Water")
state.build_spinodal()
data = state.get_spinodal_data()
```

Die relevanten Datenfelder sind:

| Feld | Bedeutung | Umrechnung |
|---|---|---|
| `data.tau` | Inverse reduzierte Temperatur | `T = T_reducing / tau` |
| `data.delta` | Reduzierte molare Dichte | `rho_molar = delta * rhomolar_reducing` |
| `data.M1` | Zusätzliche CoolProp-Kritikalitäts-/Spinodalinformation | In der bestehenden `calc_spinodale()` nicht verwendet |

**Nicht verwechseln:** `T_reducing`/`rhomolar_reducing` sind Reduktionsgrößen des Fluidmodells; sie müssen nicht für jede mögliche Stoffmodellkonfiguration exakt `T_critical`/`rhomolar_critical` entsprechen. In deinem Quellcode werden die API-Reduktionswerte korrekt separat abgefragt.

Die Daten werden in `(T, v_molar, p)` umgerechnet, und nach der kritischen Dichte aufgeteilt: `rho_molar > rho_c` = flüssigkeitsseitig, `rho_molar < rho_c` = dampfseitig. Anschließend wird der **homogene** Druck für die Stützstelle mit vorgeschriebener Phase nachgerechnet.

### 2.3 Interner CoolProp-Kurventracer: bereits im Gespräch genannte Daten

**[FRÜHERE ANGABE; vor einer wissenschaftlichen Zitierung anhand genau der verwendeten CoolProp-Version überprüfen]** Das frühere Chat-Ergebnis nannte für eine eingesehene C++-Implementierung unter anderem eine numerische Kurvenverfolgung mit Startwertsuche und Newton-/Brent-/Halley-Verfahren sowie diese internen Einstellungen:

| Im Chat genannter Wert | Zahl | Status |
|---|---:|---|
| Maximale Kurvenverfolgungsschritte | 300 | Bibliotheksintern, nicht per normaler Python-API einstellbar |
| Anfangsschritt in reduzierten Variablen `tau` und `delta` | 0,001 | Bibliotheksintern |
| Regulärer Suchradius im `delta`-Raum | 0,025 | Bibliotheksintern |
| Regulärer Suchradius im `tau`-Raum | 0,1 | Bibliotheksintern |
| Genannte Newton-Toleranz | `1e-10` | Bibliotheksintern |
| Genanntes Newton-Iterationsmaximum | 100 | Bibliotheksintern |

Diese Werte dürfen **nicht** ungeprüft als garantiert gültige Einstellungen einer unbekannten lokal installierten CoolProp-Version dargestellt werden. Der bisherige Verlauf liefert insbesondere **keinen direkten Tracer-Debugnachweis**, welcher konkrete interne Abbruchmechanismus zum fehlenden Dampfast führte.

### 2.4 Unterscheidung der einstellbaren Parameter

| Einstellung | Bisher direkt über `build_spinodal()` möglich? | ThermoSurface-Zuständigkeit |
|---|---|---|
| Fluidname `Water` usw. | Ja, bei `AbstractState` | GUI/`config.yaml` |
| Backend `HEOS` | Ja, bei `AbstractState` | Modellfestlegung |
| Beliebige Anzahl Spinodal-Temperaturpunkte | Nein | Nur mit zusätzlichem eigenem Solver |
| Beliebiges `T_step` der internen Spinodale | Nein | Nur mit zusätzlichem eigenem Solver |
| Dichte-/Volumen-Suchraster der internen Spinodale | Nein | Nur mit zusätzlichem eigenem Solver |
| Interne Tracer-Startwerte/-Schritte/-Toleranzen | Nicht über die verwendete Standard-Python-Schnittstelle | Eventuell Bibliotheksmodifikation nötig |
| Sichtbares `T_min` und `T_max` | Nicht als Neuerzeugung von HEOS-Punkten | Bereits als nachträglicher Filter umgesetzt |
| AP3-Ableitungs-Residualtoleranz | Separater Python-Code | Bereits definierbar |
| Astdarstellung/Clipping | Separater PyVista-Code | Bereits vorhanden |

---

## 3. Tatsächlich vorhandene ThermoSurface-Architektur vor einer neuen Implementierung

**[CODE]** Überprüfte Projektdateien: `coolprop_heos_pvt.py`, `coolprop_heos_common.py`, `main.py`, `Debug.py`, `visualization.py`, `gui.py`, `config.yaml`, `README_AP3.md`; ergänzend AP2-/AP4-README und Projektübersicht.

### 3.1 `coolprop_heos_pvt.py`: `calc_spinodale(stoff)`

Der für die Ursache maßgebliche bestehende Ablauf lautet gekürzt, aber semantisch unverändert:

```python
spin_state = _new_state(stoff)
spin_state.build_spinodal()
data = spin_state.get_spinodal_data()
tau = np.asarray(data.tau, dtype=float)
delta = np.asarray(data.delta, dtype=float)

T_reducing = float(spin_state.T_reducing())
rho_reducing = float(spin_state.rhomolar_reducing())
Tc = float(eos_state.T_critical())
rho_c = float(eos_state.rhomolar_critical())

T = T_reducing / tau_i
rho_molar = delta_i * rho_reducing
liquid = rho_molar > rho_c
phase = CP.iphase_liquid if liquid else CP.iphase_gas

eos_state.specify_phase(phase)
eos_state.update(CP.DmolarT_INPUTS, rho_molar, T)
pressure = float(eos_state.p())
eos_state.unspecify_phase()
```

Zusätzliche Details des real vorhandenen Codes:

- Verwirft ungültige `tau`/`delta`, nichtpositive Dichten/Temperaturen, `T >= Tc` und Punkte praktisch exakt bei `rho_c`.
- Sortiert beide Äste aufsteigend nach Temperatur, verwirft doppelte Temperaturstützstellen.
- Fehlgeschlagene Druckberechnungen bleiben als `NaN` statt unbemerkt gelöscht zu werden.
- Hängt **keinen künstlichen kritischen Punkt** in dieser Rohfunktion an (Änderung aus AP3).
- Gibt zwei `np.ndarray` mit der Form `(n,3)` im Format `(T [K], v [m³/mol], p [Pa])` zurück.
- Die Rohfunktion besitzt bislang **keine** Parameter `n_T`, `n_rho_scan`, `T_min` oder `T_max` für die interne Stützstellensuche.

### 3.2 `filter_spinodale` und `prepare_spinodal_branch`

```python
filter_spinodale(kurve, T_min, T_max)
prepare_spinodal_branch(kurve, critical, T_min, T_max, branch)
```

`filter_spinodale` erhält einen bestehenden Rohdatensatz und behält nur Temperaturen innerhalb `[T_min, T_max]`. Es berechnet **keine neuen Spinodalpunkte**.

`prepare_spinodal_branch` prüft Form, Phasenseite, Druck-Vorzeichen, Sortierung, Lücken, minimale Punktzahl, Temperaturabdeckung sowie den optionalen kritischen Anschluss. Das zugehörige `report` enthält:

```text
branch, points, subcritical_points, invalid_states,
nonpositive_pressure, wrong_phase, covers_start, covers_end,
critical_appended, temperature_order_ok, large_gaps,
geometry_complete
```

Kritischer Punkt wird nur bei passendem Nutzerintervall (`T_min <= Tc <= T_max`), mindestens zwei vorhandenen Punkten und ausreichender geometrischer Annäherung angefügt. Die bestehende Annäherungsheuristik aus AP3 prüft u. a. `Tc - T_last <= 0.01*Tc` und `abs(v_last/vc - 1) <= 0.5`. Eine ausgedehnte Kurve wird nicht durch einen künstlichen Endpunkt vorgetäuscht.

Das aktuelle `geometry_complete` verlangt unter anderem, dass **keine** unplottbaren (insbesondere nichtpositiven Druck-)Punkte auftreten. Eine mathematisch vorhandene Spinodale kann daher bei logarithmischer Darstellung berechtigt als `geometry_complete=False` gemeldet werden.

### 3.3 `main.py`: tatsächlicher Berechnungs- und Validierungspfad

```python
spin_cfg = visu_raster["spinodale_visu"]
T_min = float(spin_cfg["T_start_spinodal"])
T_max = float(spin_cfg["T_stop_spinodal"])

raw_liquid, raw_vapor = _berechnen(
    "Spinodale", pvt.calc_spinodale, stoff
)

spinodale_liquid, liquid_report = pvt.prepare_spinodal_branch(
    raw_liquid, kritisch, T_min, T_max, "liquid"
)
spinodale_vapor, vapor_report = pvt.prepare_spinodal_branch(
    raw_vapor, kritisch, T_min, T_max, "vapor"
)

spinodal_derivatives = test_spinodale_ableitung(
    stoff, spinodale_liquid, spinodale_vapor, kritisch
)
```

Der Originalcode besitzt hierfür Fehlerbehandlung um `calc_spinodale`, führt die Prüfung auch ohne aktiven `DEBUG_MODE` aus und erzeugt `daten["spinodal_validation"]` mit Status je Ast. Der Gesamtstatus folgt der Regel

```python
complete = geometry_complete and derivative.passed
```

Die GUI stellt Warnstatus und Tooltip für nicht vollständige Äste bereit. AP1/AP2/AP3 sollen bei einer Erweiterung nicht beschädigt werden.

### 3.4 `Debug.py`: AP3-Ableitungsprüfung

Auszug der vorhandenen Schnittstelle:

```python
def test_spinodale_ableitung(
    stoff, liquid, vapor, critical, *,
    tolerance=1e-2,
    critical_exclusion_rel=1e-3,
    max_samples=128,
    state_factory=None,
    cp_module=None
):
    ...
```

Das geprüfte dimensionslose Residuum ist

\[
R_s = \left|\frac{\rho_c}{p_c}
     \left(\frac{\partial p}{\partial\rho_m}\right)_T\right|.
\]

- Standard-Residualtoleranz `1e-2` (kein universell normativer Referenzwert).
- `critical_exclusion_rel=1e-3` und mindestens `0.5 K` Schutzabstand zum kritischen Punkt (gemäß AP3-Dokumentation).
- `max_samples=128` je Ast.
- Bei zu wenigen Punkten oder fehlgeschlagenen Ableitungen kein `passed=True`.
- Statuswerte: `candidates`, `sampled`, `excluded`, `evaluated`, `errors`, `max_residual`, `median_residual`, `tolerance`, `passed`.
- Phasenfestlegung und `state.unspecify_phase()` nach jeder Berechnung.

### 3.5 `config.yaml`: tatsächliche Ausgangswerte

```yaml
grid_points_creation:
  spinodale_grid:
    T_start_spinodal: 100.0
    T_stop_spinodal: 1000.0

visu_points_creation:
  spinodale_visu:
    T_start_spinodal: 580.0
    T_stop_spinodal: 647.0
```

`spinodale_grid` war in der GUI als **Reserve** beschriftet. Für die tatsächliche Spinodalberechnung in `main.py` werden die Werte von `spinodale_visu` verwendet. Der automatische Stoffrastergenerator kann diese Werte beim Wechsel des Fluids überschreiben. Beim künftigen Implementieren deshalb **nicht davon ausgehen**, dass die ursprünglichen YAML-Werte immer den aktuellen GUI-Werten entsprechen.

### 3.6 `visualization.py`: Projektion und Clipping

Die p-v-T-Rohkoordinaten `(T, v, p)` werden nach

```python
x = np.log10(v / vc)
y = T / Tc
z = np.log10(p / pc)
```

transformiert. `p <= 0` führt zu keinem endlichen logarithmischen `z`. Aktuell:

```python
actors["spinodale"] = add_spinodals(
    plotter, daten, min_z=single_phase_grid.bounds[4]
)
```

Das kann Teile der Kurve unterhalb der unteren sichtbaren Flächengrenze zusätzlich abschneiden. Ein vorgeschlagener **reiner Diagnosetest**, bislang nicht mit separatem Resultat nachgewiesen, war:

```python
actors["spinodale"] = add_spinodals(
    plotter, daten, min_z=None
)
```

**Das Umstellen von `min_z` erzeugt keine fehlenden HEOS-Daten und macht nichtpositive Drücke auf einer Logarithmusachse nicht plottbar.**

---

## 4. TEST 1: Bestehende Roh-Spinodaldaten, Druckvorzeichen und AP3-Validierung

### 4.1 Testzweck und Durchführung

**[MESSUNG]** Der Nutzer startete seine reguläre ThermoSurface-Anwendung und ergänzte nach dem vorhandenen `calc_spinodale`-Aufruf eine Diagnose zur Anzahl und Temperaturspanne der Rohpunkte sowie zu positiven/nichtpositiven Drücken.

Im Chat vorgeschlagener Diagnosecode:

```python
for name, kurve in (
    ("Flüssigkeit", raw_liquid),
    ("Dampf", raw_vapor),
):
    if len(kurve) == 0:
        print(f"{name}: Keine Spinodalpunkte")
        continue

    finite = np.isfinite(kurve).all(axis=1)
    positive = finite & (kurve[:, 2] > 0)

    print(f"\n{name}")
    print("Rohpunkte:", len(kurve))
    print("T-Bereich:", kurve[:, 0].min(), kurve[:, 0].max())
    print("Endliche Punkte:", np.count_nonzero(finite))
    print("Positive Drücke:", np.count_nonzero(positive))
    print("Nichtpositive Drücke:",
          np.count_nonzero(finite & (kurve[:, 2] <= 0)))

    if positive.any():
        print("Darstellbarer T-Bereich:",
              kurve[positive, 0].min(),
              kurve[positive, 0].max())
```

### 4.2 Screenshot-Ergebnis: direkt sichtbare Diagnose

Auf dem Screenshot `image(20261019-075314).png` sind unter anderem folgende Zeilen direkt zu erkennen (der Anfang der Rohdatenausgabe ist im Screenshot abgeschnitten):

```text
Positive Drücke: 20
Nichtpositive Drücke: 0
Darstellbarer T-Bereich: 625.1281013277346 647.0886875644494
```

Diese sichtbaren Zeilen gehören zum zuletzt ausgegebenen, in der Auswertung als **Dampf** zugeordneten Ast. Da der Screenshot die davor ausgegebenen Rohpunkt- und Flüssigkeitszeilen nicht vollständig zeigt, werden **keine fehlenden Konsolenwerte erfunden**. Der Screenshot stützt eine rohdatenseitige positive Druckabdeckung des Dampfzweigs erst ab ungefähr **625,1281 K** bei einer oberen berechneten Rohpunkttemperatur von ungefähr **647,0887 K**.

Im selben Bild erscheint `AP3` bei beiden Ästen als Warnung (`nicht vollstaendig validiert`), während die numerischen Ableitungen bestanden sind.

### 4.3 AP3-Report aus demselben Screenshot

| Merkmal | Flüssigkeitsast (`liquid`) | Dampfast (`vapor`) |
|---|---:|---:|
| `points` | 43 | 13 |
| `subcritical_points` | 43 | 13 |
| `invalid_states` | 0 | 0 |
| `nonpositive_pressure` | 17 | 0 |
| `wrong_phase` | 0 | 0 |
| `covers_start` | `True` | `False` |
| `covers_end` | `True` | `True` |
| `critical_appended` | `False` | `False` |
| `temperature_order_ok` | `True` | `True` |
| `large_gaps` | `False` | `False` |
| `geometry_complete` | `False` | `False` |
| Ableitung `candidates` | 43 | 13 |
| Ableitung `sampled` | 43 | 13 |
| Ableitung `excluded` | 0 | 0 |
| Ableitung `evaluated` | 43 | 13 |
| Ableitung `errors` | 0 | 0 |
| Ableitung `passed` | `True` | `True` |
| Residualgrenze `tolerance` | `0.01` | `0.01` |
| `max_residual` | `6.645380120632806e-14` | `1.3205596303343482e-14` |
| `median_residual` | `9.771269006556085e-15` | `5.7231510792322905e-15` |
| `complete` | `False` | `False` |

Die vollständigen numerischen AP3-Maxima und -Mediane sind in der obigen Tabelle originalgetreu aus der Konsolenausgabe übertragen. Diese AP3-Residuen dürfen nicht mit den separat berechneten Residuen der Tests 2 und 3 verwechselt werden. Die für die spätere Entscheidung relevante Aussage ist eindeutig: **Beide AP3-Derivatprüfungen bestanden mit sehr kleinen Residuen; beide Geometrieprüfungen scheiterten aus unterschiedlichen Gründen.**

**Exakte Fehlermeldungslogik:**

- Flüssigkeit: `geometry_complete=False`, obwohl `covers_start=True` und `covers_end=True`, maßgeblich wegen **17 nichtpositiver Drücke** auf einer logarithmischen Druckachse.
- Dampf: `geometry_complete=False` aufgrund von **`covers_start=False`**; es fehlen Starttemperaturen im gewünschten Bereich `580...647 K`.
- Beide: `critical_appended=False`, da `T_stop_spinodal=647 K` **unter** `Tc(Water)` von etwa `647,096 K` liegt. Der kritische Punkt wird nicht unzulässig außerhalb des Nutzerfilters ergänzt. Das ist korrektes AP3-Verhalten, kein Beleg für einen Solverfehler.
- Die Rohdatenausgabe mit `20` positiven Dampfpunkten und der AP3-Report mit `13` Dampfpunkten widersprechen sich nicht: Sie stammen aus **verschiedenen Phasen des Verarbeitungspfads (Rohdaten gegenüber auf `T <= 647 K` gefilterten Daten)**. Die im Screenshot abgeschnittene restliche Rohdatenausgabe lässt keine zusätzliche exakte Gesamtpunktzahl erkennen.

### 4.4 Zusätzliche Konsolenmeldung

```text
PyVistaFutureWarning: The default value of 'algorithm' for the filter
'StructuredGrid.extract_surface' will change in the future.
...
mesh = grid.extract_surface().triangulate()
```

Fundort laut Screenshot: `visualization.py`, ungefähr Zeile 260. Dies ist eine **FutureWarning zur PyVista-API**, nicht die Ursache des verkürzten Dampfastes. Eine explizite `algorithm`-Angabe wäre ggf. eine **separate** Pflegeaufgabe; sie ist nicht die gesuchte numerische Spinodalkorrektur.

### 4.5 Test-1-Schlussfolgerung

**[AUSWERTUNG]** Das bloße Abschalten von Clipping oder eine höhere PyVista-Meshauflösung reicht nicht. Der Dampfast besitzt im geprüften Rohdatensatz bereits keine geeigneten Stützstellen für den Bereich **580 bis etwa 625 K**. Bei der Flüssigkeit besteht zusätzlich ein separates **Druckvorzeichen-/Logarithmusproblem**.

---

## 5. TEST 2: Existieren HEOS-Ableitungs-Nullstellen bei 580 bis 640 K?

### 5.1 Frage

Ist die fehlende Dampfkurve eine Beschränkung der internen Spinodal-Kurvenverfolgung oder kann die zugrunde liegende HEOS-Zustandsgleichung in den fehlenden Temperaturbereichen tatsächlich keinen Spinodalpunkt erzeugen?

### 5.2 Vollständiges Testskript aus dem Chat

Dem Nutzer wurde empfohlen, ein eigenständiges Skript `test_spinodal_bereich.py` im Projektordner anzulegen (der Screenshot zeigt die Ausführung einer lokalen Datei `Untitled-1.py`). **Das Skript scannt die Ableitung, löst aber die Nullstellen noch nicht genau auf.**

```python
import numpy as np
import CoolProp.CoolProp as CP
from coolprop_heos_common import _new_state

state = _new_state("Water")
rho_c = state.rhomolar_critical()

for T in (580, 590, 600, 610, 620, 625, 630, 640):
    state.update(CP.QT_INPUTS, 1.0, float(T))
    rho_g = state.rhomolar()

    # Dampfast zwischen Sättigungsdichte und kritischer Dichte
    if not 0 < rho_g < rho_c:
        print(f"{T} K: Kein gültiges Suchintervall")
        continue

    rho_grid = np.geomspace(
        rho_g * (1 + 1e-6),
        rho_c * (1 - 1e-6),
        500
    )

    deriv = np.full(len(rho_grid), np.nan)

    state.specify_phase(CP.iphase_gas)
    try:
        for i, rho in enumerate(rho_grid):
            try:
                state.update(CP.DmolarT_INPUTS, float(rho), float(T))
                deriv[i] = state.first_partial_deriv(
                    CP.iP, CP.iDmolar, CP.iT
                )
            except (ValueError, RuntimeError):
                continue
    finally:
        state.unspecify_phase()

    valid = np.isfinite(deriv)
    crossings = np.where(
        valid[:-1] & valid[1:] &
        (np.signbit(deriv[:-1]) != np.signbit(deriv[1:]))
    )[0]

    print(
        f"{T} K: "
        f"{len(crossings)} Nullstellenintervalle, "
        f"{valid.sum()}/{len(valid)} gültige Werte"
    )
```

**Testparameter:**

- Fluid: `Water`, Backend: `HEOS`.
- Temperaturwerte: `580`, `590`, `600`, `610`, `620`, `625`, `630`, `640 K` (acht Temperaturen).
- Dichteintervall pro Temperatur: vom minimal oberhalb der gesättigten Dampfdichte gelegenen Wert `rho_g*(1+1e-6)` bis minimal unter `rho_c*(1-1e-6)`.
- Raster: **500 logarithmisch äquidistante Dichtestützstellen** pro Temperatur (`np.geomspace`).
- Phasenvorgabe: `CP.iphase_gas`.
- Kriterium: Vorzeichenwechsel von `first_partial_deriv(CP.iP, CP.iDmolar, CP.iT)` zwischen benachbarten endlichen Rasterwerten.
- Es handelt sich um Vorzeichenwechselintervalle, nicht schon um durch einen Solver validierte exakte Nullstellen.

### 5.3 Vollständige vom Nutzer zurückgegebene Ausgabe

Screenshot: `image(20261019-075930).png`.

```text
580 K: 2 Nullstellenintervalle, 500/500 gültige Werte
590 K: 2 Nullstellenintervalle, 500/500 gültige Werte
600 K: 2 Nullstellenintervalle, 500/500 gültige Werte
610 K: 2 Nullstellenintervalle, 500/500 gültige Werte
620 K: 2 Nullstellenintervalle, 500/500 gültige Werte
625 K: 2 Nullstellenintervalle, 500/500 gültige Werte
630 K: 2 Nullstellenintervalle, 500/500 gültige Werte
640 K: 2 Nullstellenintervalle, 500/500 gültige Werte
```

Die am oberen Screenshotrand erscheinende Meldung `KeyboardInterrupt` gehört erkennbar zu einer vorherigen/unterbrochenen Sitzung und **nicht** als Fehler zu den hier vollständig ausgegebenen acht Testzeilen.

### 5.4 Auswertung Test 2

- **[MESSUNG]** Bei allen acht Temperaturen findet der Dichtescan **zwei** Vorzeichenwechselintervalle.
- **[MESSUNG]** Alle `8 x 500 = 4.000` geprüften Ableitungswerte waren endlich bzw. laut Skript gültig.
- **[AUSWERTUNG]** Die fehlende Rohdatenabdeckung ab `580 K` ist **nicht** dadurch erklärbar, dass HEOS dort prinzipiell keine auswertbare Ableitung bereitstellt.
- **[OFFEN]** Der Test belegt noch **nicht**, dass beide gefundenen Nullstellen physikalisch zur gewünschten Dampfspinodale gehören.
- **[OFFEN]** Die konkrete Abbruchursache des eingebauten `build_spinodal()`-Kurventracers ist nicht gemessen worden.

---

## 6. TEST 3: Nullstellen exakt bestimmen und beide Kandidaten thermodynamisch klassifizieren

### 6.1 Ziel

Jeder der zwei Vorzeichenwechsel aus Test 2 wird mit `scipy.optimize.brentq` als Nullstelle bestimmt. Ausgegeben werden reduzierte Dichte, reduziertes molares Volumen, reduzierter Druck, Druck relativ zum Sättigungsdruck, normiertes Ableitungsresiduum und Vorzeichenrichtung. Damit soll die **lokale Druckmaximum-Nullstelle** von der weiteren Nullstelle unterschieden werden.

### 6.2 Vollständiges Testskript aus dem Chat

```python
import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import brentq
from coolprop_heos_common import _new_state

state = _new_state("Water")
rho_c = state.rhomolar_critical()
pc = state.p_critical()

for T in [580, 600, 620, 630, 640]:
    state.update(CP.QT_INPUTS, 1.0, float(T))
    rho_g = state.rhomolar()
    p_sat = state.p()

    def deriv(rho):
        state.specify_phase(CP.iphase_gas)
        try:
            state.update(CP.DmolarT_INPUTS, float(rho), float(T))
            return float(state.first_partial_deriv(
                CP.iP, CP.iDmolar, CP.iT
            ))
        finally:
            state.unspecify_phase()

    densities = np.geomspace(
        rho_g * (1 + 1e-6),
        rho_c * (1 - 1e-6),
        500
    )
    values = np.array([deriv(rho) for rho in densities])

    print(f"\nTemperatur: {T} K, p_sat = {p_sat/1e6:.5f} MPa")

    for i in range(len(densities) - 1):
        a, b = values[i], values[i+1]

        if not np.isfinite([a, b]).all() or a * b >= 0:
            continue

        rho_root = brentq(
            deriv, densities[i], densities[i+1],
            rtol=1e-12
        )

        state.specify_phase(CP.iphase_gas)
        try:
            state.update(CP.DmolarT_INPUTS, rho_root, float(T))
            p_root = state.p()
        finally:
            state.unspecify_phase()

        residual = abs(deriv(rho_root)) * rho_c / pc
        crossing = "+ nach -" if a > 0 else "- nach +"

        print(
            f"  rho/rho_c = {rho_root/rho_c:.6f}, "
            f"v/v_c = {rho_c/rho_root:.6f}, "
            f"p/p_c = {p_root/pc:.6f}, "
            f"p/p_sat = {p_root/p_sat:.6f}, "
            f"R = {residual:.3e}, "
            f"Vorzeichen: {crossing}"
        )
```

**Testparameter:** `T = [580, 600, 620, 630, 640] K`; `500` logarithmische Dichtestützstellen je T; `brentq(rtol=1e-12)`; Sättigung aus `QT_INPUTS` bei `Q=1`; `iphase_gas`; dimensionsloses Residuum `R=abs(dp/drho)_T * rho_c / pc`.

### 6.3 Vollständig transkribierte Testausgabe

Screenshot: `image(20261019-080315).png`. Die ausgegebenen Dezimalzahlen sind **reduzierte** Größen; negative `p/p_c`-Werte der zweiten Nullstelle sind echte im Skript ausgegebene homogene EOS-Fortsetzungswerte, keine negativen Sättigungsdrücke.

```text
Temperatur: 580 K, p_sat = 9.44797 MPa
  rho/rho_c = 0.274950, v/v_c = 3.637028, p/p_c = 0.515090, p/p_sat = 1.202898, R = 2.819e-15, Vorzeichen: + nach -
  rho/rho_c = 0.938031, v/v_c = 1.066062, p/p_c = -634.769867, p/p_sat = -1482.389091, R = 6.296e-09, Vorzeichen: - nach +

Temperatur: 600 K, p_sat = 12.34482 MPa
  rho/rho_c = 0.355207, v/v_c = 2.815258, p/p_c = 0.628857, p/p_sat = 1.123961, R = 3.141e-15, Vorzeichen: + nach -
  rho/rho_c = 0.935994, v/v_c = 1.068383, p/p_c = -31.239440, p/p_sat = -55.834493, R = 9.275e-12, Vorzeichen: - nach +

Temperatur: 620 K, p_sat = 15.90058 MPa
  rho/rho_c = 0.464580, v/v_c = 2.152482, p/p_c = 0.762474, p/p_sat = 1.058026, R = 8.925e-15, Vorzeichen: + nach -
  rho/rho_c = 0.928674, v/v_c = 1.076804, p/p_c = -0.159955, p/p_sat = -0.221957, R = 1.586e-13, Vorzeichen: - nach +

Temperatur: 630 K, p_sat = 17.96897 MPa
  rho/rho_c = 0.541137, v/v_c = 1.847960, p/p_c = 0.839532, p/p_sat = 1.030857, R = 4.876e-14, Vorzeichen: + nach -
  rho/rho_c = 0.924121, v/v_c = 1.082110, p/p_c = 0.707422, p/p_sat = 0.868640, R = 4.711e-15, Vorzeichen: - nach +

Temperatur: 640 K, p_sat = 20.26521 MPa
  rho/rho_c = 0.658638, v/v_c = 1.518284, p/p_c = 0.926810, p/p_sat = 1.009076, R = 1.675e-15, Vorzeichen: + nach -
  rho/rho_c = 0.952539, v/v_c = 1.049826, p/p_c = 0.911964, p/p_sat = 0.992912, R = 1.005e-14, Vorzeichen: - nach +
```

### 6.4 Nullstelle A: `+ nach -`, lokales Druckmaximum, gewünschter Kandidat

| T (K) | `p_sat` (MPa) | `rho/rho_c` | `v/v_c` | `p/p_c` | `p/p_sat` | R |
|---:|---:|---:|---:|---:|---:|---:|
| 580 | 9,44797 | 0,274950 | 3,637028 | 0,515090 | 1,202898 | 2,819e-15 |
| 600 | 12,34482 | 0,355207 | 2,815258 | 0,628857 | 1,123961 | 3,141e-15 |
| 620 | 15,90058 | 0,464580 | 2,152482 | 0,762474 | 1,058026 | 8,925e-15 |
| 630 | 17,96897 | 0,541137 | 1,847960 | 0,839532 | 1,030857 | 4,876e-14 |
| 640 | 20,26521 | 0,658638 | 1,518284 | 0,926810 | 1,009076 | 1,675e-15 |

**[MESSUNG]** Der Vorzeichenwechsel ist bei allen fünf Punkten `+ nach -`, also eine lokale Druckmaximum-Stelle entlang der untersuchten Isotherme. Stets `p_spin > p_sat`, `rho_spin < rho_c`, `v_spin > vc`. Die Werte nähern sich beim Anstieg von T sichtbar dem kritischen Zustand. Die normierten Residuen liegen im Bereich ungefähr `1.7e-15` bis `4.9e-14`.

**[AUSWERTUNG]** Die Reihe ist ein **sehr plausibler Kandidat für die bislang fehlende Dampfspinodale von Wasser**. Im gemeinsam mit der eingebauten HEOS-Kurve abgedeckten Bereich ist aber eine echte **punktweise Vergleichsmessung** noch offen. Die fünf Stichproben allein beweisen keine lückenlose Gültigkeit für jeden Temperaturwert und nicht für andere Fluide.

### 6.5 Nullstelle B: `- nach +`, lokales Druckminimum, nicht mit Dampfspinodale verwechseln

| T (K) | `rho/rho_c` | `v/v_c` | `p/p_c` | `p/p_sat` | R |
|---:|---:|---:|---:|---:|---:|
| 580 | 0,938031 | 1,066062 | -634,769867 | -1482,389091 | 6,296e-09 |
| 600 | 0,935994 | 1,068383 | -31,239440 | -55,834493 | 9,275e-12 |
| 620 | 0,928674 | 1,076804 | -0,159955 | -0,221957 | 1,586e-13 |
| 630 | 0,924121 | 1,082110 | 0,707422 | 0,868640 | 4,711e-15 |
| 640 | 0,952539 | 1,049826 | 0,911964 | 0,992912 | 1,005e-14 |

**[MESSUNG]** Der Vorzeichenwechsel lautet durchweg `- nach +`, also ein **lokales Druckminimum**, kein Dampf-Druckmaximum. Bei 580, 600 und 620 K gibt die homogene EOS an dieser zweiten Nullstelle **negative Drücke** aus, teils extrem negativ. Bei 630 und 640 K sind die Drücke positiv, aber **unter dem Sättigungsdruck**.

**[AUSWERTUNG]** Diese mathematische Nullstelle darf **nicht allein wegen des kleinen Residuums** als gewünschter dampfseitiger Spinodalast ausgewählt werden. Eine besondere Vorsicht gilt deshalb bei einem allgemeinen Nullstellensolver: Er kann abhängig vom Suchintervall und Startwert die unerwünschte zweite Nullstelle finden.

**Nuance:** Die bloße Bezeichnung der zweiten Nullstelle als allgemein „unphysikalisch“ wäre zu pauschal. Ihr Auftreten ist als mathematischer Extrempunkt der fortgesetzten verwendeten HEOS-EOS gezeigt. Welche modelltechnische Bedeutung, Extrapolationsqualität und Stabilitätszweigzuordnung sie besitzt, wurde nicht unabhängig untersucht. Für die **hier gesuchte metastabile Dampf-Stabilitätsgrenze** ist sie nach den beobachteten Kennzeichen nicht der passende Kandidat.

### 6.6 Reduzierte Diagrammkoordinaten eines ergänzbaren Punktes

Aus den **gemessenen** Werten für 580 K, Nullstelle A:

\[
\log_{10}(v/v_c) = \log_{10}(3{,}637028) \approx 0{,}561,
\]

\[
\log_{10}(p/p_c) = \log_{10}(0{,}515090) \approx -0{,}288.
\]

Für den 3D-Plot kommt `y = T/Tc` hinzu. Auch bei gleicher Temperaturabdeckung wird der orange Ast deshalb **nicht dieselbe räumliche Form oder Ausdehnung wie der violette** besitzen.

---

## 7. Gesicherte Ergebnisse aus allen drei Tests zusammen

| Frage | Antwort / Evidenz | Einstufung |
|---|---|---|
| Ist der kurze Dampfast im GUI allein ein Zeichen für zu wenige PyVista-Punkte? | Nein. Der relevante untere Temperaturbereich fehlt bereits in den ursprünglichen `build_spinodal()`-Rohdaten. | [MESSUNG] + [CODE] |
| Was ist die Rohdatengrenze der positiven Dampfpunkte? | `625.1281013277346 K` bis `647.0886875644494 K` (sichtbarer Ausschnitt der Konsolenausgabe). | [MESSUNG] |
| Was passiert bei Filterung auf 580 bis 647 K? | AP3-Dampfast enthält nur `13` Punkte, `covers_start=False`, `covers_end=True`. | [MESSUNG] |
| Können HEOS-Ableitungen unterhalb 625 K ausgewertet werden? | Ja, bei acht getesteten T-Werten jeweils `500/500` gültige Ableitungen. | [MESSUNG] |
| Werden dort Nullstellen gefunden? | Ja, jeweils zwei Vorzeichenwechselintervalle. | [MESSUNG] |
| Welche Nullstelle ist Dampfspinodal-Kandidat? | Der erste `+ nach -`-Wechsel mit `p > p_sat`. | [MESSUNG] + [AUSWERTUNG] |
| Warum nicht die zweite? | `- nach +` = lokales Minimum; bei 580–620 K stark negative und bei 630–640 K untersättigte Drücke. | [MESSUNG] + [AUSWERTUNG] |
| Sind die neuen Nullstellen präzise? | Sehr kleine Residuen, siehe Originaltabellen, `brentq(rtol=1e-12)`. | [MESSUNG] |
| Ist damit der interne CoolProp-Abbruchmechanismus bekannt? | Nein. Nur gezeigt: der bestehende `build_spinodal()`-Datensatz deckt den gewünschten Dampftemperaturbereich nicht ab. | [OFFEN] |
| Ist die Verlängerung bis 580 K plausibel? | Ja, anhand von fünf direkten HEOS-Wasser-Stichproben; Serienrechnung, Vergleich/Validierung noch nötig. | [AUSWERTUNG] |
| Muss der Flüssigkeitsast genauso lang werden? | Nein. Thermodynamische Asymmetrie erlaubt unterschiedliche Druck- und Volumenverläufe. | [AUSWERTUNG] |
| Können negative Drücke auf bestehender `log10(p/pc)`-Achse gezeigt werden? | Nein. 17 Flüssigkeitsstützstellen waren im AP3-Bericht nichtpositiv. | [MESSUNG] + [CODE] |
| Warum fehlt der kritische Endpunkt? | `T_stop_spinodal=647 K < Tc(Water)≈647.096 K`, AP3 hängt außerhalb des Nutzerfensters nichts an. | [CODE] + [AUSWERTUNG] |
| Wurde der Produktivcode bereits umgeschrieben? | Nein. Es wurden nur Entwürfe und separate lokale Diagnosetests besprochen. | [STATUS] |

**Wichtigste technische Folgerung:** HEOS kann für Wasser bei `580...640 K` auf Basis derselben Helmholtz-Zustandsgleichung Druckmaxima mit verschwindender Ableitung auswerten. Eine **eigene temperaturparametrisierte Nullstellensuche** ist daher eine plausible Ergänzung zum eingebauten Kurventracer. Diese Erkenntnis ist **kein Beweis**, dass eine konkrete neue Implementierung bereits robust, korrekt integriert oder fluidübergreifend validiert wurde.

---

## 8. Im Chat diskutierte Lösungsansätze, aber noch **nicht implementiert**

### 8.1 Variante A: Originale CoolProp-Spinodale unverändert übernehmen

- `build_spinodal()`/`get_spinodal_data()` nur einmal aufrufen.
- Beide geliefert Äste filtern, nachrechnen und via AP3 validieren.
- Vorteil: Geringer Wartungs- und Implementierungsaufwand, existierender Workflow erhalten.
- Nachteil: Bei Wasser zeigt die Nutzermessung einen fehlenden Dampfast unter etwa 625 K.
- Eine reine grafische Interpolation glättet nur vorhandene Stützstellen, ersetzt aber **keine fehlenden physikalisch berechneten Nullstellen**.

### 8.2 Variante B: Zusätzliche isotherme Nullstellensuche mit HEOS

**[VORSCHLAG]** Für festgelegte Temperaturen im subkritischen Intervall:

1. Für jedes `T` die gesättigte Dampfdichte `rho_g(T)` und `p_sat(T)` ermitteln.
2. Die kritische molare Dichte `rho_c` verwenden.
3. Homogene Dichtepunkte zwischen `rho_g(T)` und `rho_c` logarithmisch abtasten.
4. Die analytische HEOS-Ableitung \((\partial p/\partial\rho_m)_T\) auswerten.
5. Kandidaten-Intervalle mit Vorzeichenwechsel **`+ nach -`** suchen.
6. Die ausgewählte Nullstelle mit `scipy.optimize.brentq` bestimmen.
7. Dort den **homogenen** Druck nachrechnen und `p > p_sat` prüfen.
8. `R_s` quantitativ prüfen, Kontinuität und Monotonie/Branchenidentität über benachbarte Temperaturen kontrollieren.
9. `T`, `v=1/rho`, `p` als neues Punktarray `(n,3)` liefern; fehlende Temperaturen und Solverfehler **nicht stillschweigend interpolieren**.
10. Alte und neue HEOS-Spinodalen im überlappenden Temperaturbereich vergleichen; bestehende AP3-Validierung beibehalten.

### 8.3 Vorgeschlagene Signaturen und Parameter

Im Chat tauchten zwei Varianten auf:

```python
def calc_spinodale_isotherm(stoff, T_min, T_max, n_T=250, n_rho_scan=160):
    ...
```

und konkret für den Dampfzweig:

```python
def calc_spinodale_vapor_temperature(
    stoff, T_min, T_max, n_T=250, n_scan=200
):
    ...
```

Die zweite Variante würde den **Flüssigkeitsast zunächst unverändert** bei `build_spinodal()` belassen und nur den lückenhaften Dampfast ergänzen. Das senkt das Änderungsrisiko im Vergleich zu einem vollständigen Austausch der gesamten Spinodalroutine. Ein langfristig einheitlicher Solver für beide Äste wäre eine mögliche Folgeentscheidung, **aber nicht durch diese Tests notwendig geworden**.

### 8.4 Bereits vorgeschlagener grundlegender Solver-Baustein

**[VORSCHLAG, nicht selbst getestet]** Der frühere Chat enthielt als generischen Kern:

```python
from scipy.optimize import brentq

def finde_spinodal_dichte(state, T, rho_min, rho_max,
                          phase, n_scan=160):

    def ableitung(rho):
        state.specify_phase(phase)
        try:
            state.update(CP.DmolarT_INPUTS, float(rho), float(T))
            return float(state.first_partial_deriv(
                CP.iP, CP.iDmolar, CP.iT
            ))
        finally:
            state.unspecify_phase()

    rho_grid = np.geomspace(rho_min, rho_max, n_scan)
    roots = []

    for a, b in zip(rho_grid[:-1], rho_grid[1:]):
        try:
            fa, fb = ableitung(a), ableitung(b)
            if not np.isfinite([fa, fb]).all():
                continue
            if fa * fb < 0:
                roots.append(brentq(
                    ableitung, a, b, rtol=1e-10
                ))
        except (ValueError, RuntimeError, OverflowError):
            continue

    return roots
```

In der ursprünglichen Erklärung wurde vorgeschlagen, die Flüssigkeits-/Dampf-Suchintervalle mit

\[
\rho_g(T)<\rho_{\mathrm{spin},g}<\rho_c,
\qquad
\rho_c<\rho_{\mathrm{spin},l}<\rho_l(T)
\]

zu strukturieren. **Diese Intervalle sind Arbeitsannahmen für die betrachteten üblichen unterkritischen Flüssig-Dampf-Zustände und sind nicht ohne Validierung auf jedes Fluid bzw. alle Temperaturen verallgemeinerbar.**

Der generische Baustein liefert **beide** Vorzeichenwechsel; genau wegen der Testergebnisse ist eine fachlich geprüfte Auswahl des gewünschten Astes erforderlich.

### 8.5 Im Chat vollständig vorgeschlagener Prototyp zur Dampfspinodale

**[VORSCHLAG; nicht in Projekt integriert, nicht vom Nutzer als fertig getestet gemeldet]** Dieser Code wird unverändert als Referenz für das spätere Arbeitspaket aufbewahrt:

```python
def calc_spinodale_vapor_temperature(
    stoff, T_min, T_max, n_T=250, n_scan=200
):
    from scipy.optimize import brentq

    state = _new_state(stoff)
    Tc = float(state.T_critical())
    rho_c = float(state.rhomolar_critical())

    T_lower = max(T_min, state.Ttriple() + 0.001)
    T_upper = min(T_max, Tc - 0.05)

    if T_lower >= T_upper or n_T < 2 or n_scan < 20:
        raise ValueError("Ungueltiges Spinodalraster")

    temperatures = np.linspace(T_lower, T_upper, n_T)
    results = []

    def evaluate(T, rho, derivative=True):
        state.specify_phase(CP.iphase_gas)
        try:
            state.update(CP.DmolarT_INPUTS, float(rho), float(T))
            if derivative:
                return float(state.first_partial_deriv(
                    CP.iP, CP.iDmolar, CP.iT
                ))
            return float(state.p())
        finally:
            state.unspecify_phase()

    for T in temperatures:
        try:
            state.update(CP.QT_INPUTS, 1.0, float(T))
            rho_g = float(state.rhomolar())
            p_sat = float(state.p())

            if not 0 < rho_g < rho_c:
                continue

            densities = np.geomspace(
                rho_g * (1 + 1e-6),
                rho_c * (1 - 1e-6),
                n_scan
            )

            deriv = np.full(n_scan, np.nan)
            for i, rho in enumerate(densities):
                try:
                    deriv[i] = evaluate(T, rho)
                except (ValueError, RuntimeError):
                    pass

            for i in range(n_scan - 1):
                fa, fb = deriv[i], deriv[i + 1]

                # Dampfspinodale: lokales Druckmaximum
                if not (np.isfinite(fa) and np.isfinite(fb)):
                    continue
                if not (fa > 0 and fb < 0):
                    continue

                root = brentq(
                    lambda rho: evaluate(T, rho),
                    densities[i],
                    densities[i + 1],
                    rtol=1e-11
                )

                p = evaluate(T, root, derivative=False)

                if np.isfinite(p) and p > p_sat:
                    results.append((T, 1.0 / root, p))
                    break

        except (ValueError, RuntimeError):
            LOGGER.warning(
                "Dampfspinodale bei T=%g K nicht berechenbar", T
            )

    if not results:
        return np.empty((0, 3), dtype=float)

    return np.asarray(results, dtype=float)
```

**Wesentliche Eigenschaften des Prototyps:**

- Intern wird `state = _new_state(stoff)` nur einmal erzeugt.
- Temperatur wird durch `np.linspace(T_lower, T_upper, n_T)` vorgegeben.
- Obergrenze `Tc - 0.05 K`; Untergrenze `Ttriple + 0.001 K`.
- Dichtescan `np.geomspace` mit `n_scan` Punkten, Standardwert `200`.
- Rückgabewerte `(T, 1/rho, p)`; gleiche Datenkonvention wie bestehender Code.
- Akzeptiert nur `+ nach -` und `p > p_sat`.
- Bei fehlenden Zuständen wird kein Punkt erzeugt; es entstehen dann unter Umständen ungleichmäßige/gap-behaftete Temperaturkurven.

**Achtung für spätere Umsetzung:** Der Prototyp setzt voraus, dass im Zielmodul bereits `_new_state`, `CP`, `np`, `LOGGER` importiert sind (im jetzigen `coolprop_heos_pvt.py` grundsätzlich vorhanden). SciPy muss vorhanden und in `requirements`/Abhängigkeiten berücksichtigt werden. **Das ist ein Prototyp, kein robuster, fluidübergreifend akzeptierter Produktivcode.** Insbesondere müssen Exceptions, kritisch nahe Zweige, Grenzen, Mehrfachnullstellen, Abtastauflösung und Datenlücken weiter validiert werden.

### 8.6 Im Chat vorgeschlagene Anbindung in `main.py`

Nur als früherer Entwurf festgehalten:

```python
raw_vapor_new = _berechnen(
    "Dampfspinodale: temperaturgesteuert",
    pvt.calc_spinodale_vapor_temperature,
    stoff,
    T_min,
    T_max,
    n_T=250,
    n_scan=200
)

if len(raw_vapor_new) >= 2:
    raw_vapor = raw_vapor_new
```

**Noch nicht umgesetzt.** Fachlich wäre es vorzuziehen, die Auswahlregeln zu schärfen, statt `len >= 2` als ausreichendes Kriterium zu akzeptieren: Eine fast vollständig fehlende Kurve dürfte sonst bereits mit zwei Punkten den bisherigen Datensatz ersetzen. Das vorherige `prepare_spinodal_branch` plus die AP3-Ableitungsprüfung müssen auch für den neuen Ast weiterhin gelten.

### 8.7 Vorgeschlagene neue `config.yaml`-Parameter

**[VORSCHLAG, nicht aktueller Dateistand]**

```yaml
spinodale_grid:
  T_start_spinodal: 580.0
  T_stop_spinodal: 647.0
  n_T_spinodal: 250
  n_rho_scan: 160
```

Der Chat enthielt für die eigentliche Solverfunktion teils `n_scan=200` statt `160` sowie `n_T=250`. Es handelt sich **um vorgeschlagene Startauflösungen, nicht um empirisch konvergenzvalidierte optimale Werte**.

Wichtig für die spätere Architekturplanung:

- Aktuell ist `spinodale_grid` eine Reserve; **`main.py` verwendet `spinodale_visu`** für den Filter. Das zukünftige Verhalten beider Einstellgruppen muss bewusst spezifiziert werden.
- Neue Konfigurationsschlüssel benötigen Validierung, Übergabe durch `main.py`, Eingabefelder/Defaultwerte in `gui.py` und sinnvolle Grenzen.
- Die maximale globale Punktzahl aus `config.yaml` ist kein Ersatz für eine Abschätzung der erforderlichen **HEOS-Ableitungsauswertungen**.
- Bei `n_T=250` und `n_scan=200` fallen in einem naiven Vollscan rund **50.000 Dichteableitungsauswertungen** plus Root-Solver- und Sättigungsaufrufe an; die reale Laufzeit wurde bislang **nicht gemessen**.

---

## 9. Wissenschaftlich-technische Risiken, die bei der späteren Arbeitspaketdefinition nicht verloren gehen dürfen

### 9.1 Mehrfachnullstellen der Helmholtz-EOS

**Nachgewiesen:** Zwei Ableitungs-Nullstellen treten im untersuchten Intervall zwischen Sättigungsdampf- und kritischer Dichte auf. Ein generisches `brentq` mit beliebigem Intervall bzw. ein bloßes `abs(residual) < tolerance` reicht deshalb **nicht** aus. Der Dampfzweig muss nach *Druckmaximum*, Phasen-/Dichteseite, relativer Sättigungsdrucklage und Kurvenkontinuität klassifiziert werden.

### 9.2 Numerische vs. thermodynamische Gültigkeit

Ein kleines Ableitungsresiduum beweist nur, dass eine numerische stationäre Stelle der gewählten homogenen EOS gefunden wurde. Für die Interpretation als metastabile Dampfspinodale sind weitere physikalische Kriterien und die Anschlussfähigkeit an den vorhandenen Ast nötig. Die zweite Nullstelle illustriert genau diese Unterscheidung.

### 9.3 Nahe dem kritischen Punkt

Die Spinodaläste konvergieren dort, und ein festes Dichtescanning kann den Wechsel oder zwei sehr nahe Nullstellen übersehen. Der bisher vorgeschlagene Schutzabstand `Tc - 0.05 K` ist eine **Heuristik**, keine garantierte Gültigkeitsgrenze. AP3 schließt seinerseits einen kritischen Bereich aus der Residualstichprobe aus. Die Endpunktlogik darf nicht durch erzwungene kritische Kunstpunkte ersetzt werden.

### 9.4 Temperaturfilter, Abdeckung, Lücken

Bei streng vorgegebenem `T_values`-Raster sollte ein fehlgeschlagener Punkt nicht zu einem falsch als „durchgehend“ interpretierten Segment führen. Ausgabe und `spinodal_validation` müssen erkennbar machen, welche Temperaturen berechnet, verworfen oder mangels gültiger Nullstelle nicht abgedeckt sind. Interpolation ist keine Ersatzberechnung.

### 9.5 Nichtpositive Drücke und logarithmische Darstellung

Der Flüssigkeitsast hat bei der dokumentierten Messung **17 nichtpositive Punkte**, obwohl die Ableitungsresiduen bestanden haben. `p <= 0` darf nicht als EOS-Fehler entfernt werden, wenn es eine mathematisch gültige homogene Fortsetzung ist; die `log10(p/pc)`-Achse ist dafür schlicht ungeeignet. Optionen wie lineare/symmetrisch transformierte Achsen oder getrennte wissenschaftliche Diagnosedarstellung wären ein anderes, ausdrücklich zu entscheidendes Thema. Eine stillschweigende Druckverschiebung zum Erzwingen positiver Logwerte wäre falsch.

### 9.6 Clipping und Mesh-Grenzen

`min_z=single_phase_grid.bounds[4]` ist ein zusätzlicher Renderer-Clip; `min_z=None` ist für eine isolierte Sichtbarkeitsdiagnose sinnvoll, löst aber die im Rohdatensatz fehlenden Dampfpunkte nicht. Ein Verfahren soll nur bis zur Spinodal-Stabilitätsgrenze rechnen, **nicht** bis zur willkürlichen äußeren Zustandsoberflächenkante.

### 9.7 Andere Fluide

Der dokumentierte Realtest gilt **nur für Wasser**. Auf andere HEOS-Fluide darf insbesondere die Regel „erste Nullstelle mit `+ nach -` und `p > p_sat`“ nicht ohne konkrete Tests übertragen werden. Die spätere Abnahme sollte wenigstens ein zweites repräsentatives Fluid, fluidabhängige Gültigkeitsgrenzen sowie mögliche EOS-Artefakte umfassen.

### 9.8 Regression und Schnittstellen

- Bestehende Arrays bleiben `(T,v,p)` und nicht `(p,v,T)`.
- `calc_spinodale(stoff)` und neue Hilfsfunktionen dürfen keine ungewollten Seiteneffekte in anderen Diagrammtypen (`pht`, `pst`, `hsp`, `uvt`) erzeugen.
- Die vorhandene `spinodal_validation`-/GUI-Statuslogik und AP3-Ableitungsprüfung sind zu erhalten.
- AP1-Snapshots können sich durch geänderte Spinodallängen berechtigt ändern. Alte Baselines **nicht ungeprüft überschreiben**.
- AP2-Randkontinuität der Siede-/Taulinien und Zweiphasenflächen darf nicht verändert werden.

### 9.9 Dokumentations- und Performance-Ansprüche

Die Aussage „CoolProp spart X Iterationen gegenüber einem selbstgeschriebenen van-der-Waals-Solver“ wurde **nicht** durch diesen Chat nachgewiesen. Für eine solche Behauptung braucht es instrumentierte Laufzeit- und Solverzählungen unter identischen Bedingungen. Auch ein weiterer Vergleich `build_spinodal()` vs. eigener HEOS-Solver benötigt echte Messungen.

---

## 10. Nachweislich noch offene Punkte / mögliche spätere Prüfgegenstände

**Dieser Abschnitt ist eine vollständige Merkliste für die spätere Definition eines Arbeitspakets, aber absichtlich noch kein nummeriertes Arbeitspaket, kein Implementierungsauftrag und keine Aufwandsschätzung.**

1. **[OFFEN]** `calc_spinodale()`-Rohdaten für Wasser und `calc_spinodale_vapor_temperature()` punktweise im überlappenden Bereich `T≈625...647 K` vergleichen (z. B. in `T`, `v`, `p`, \(R_s\)); dabei berücksichtigen, dass die beiden Algorithmen unterschiedliche Temperaturstützstellen liefern.
2. **[OFFEN]** Vollständigen Dampfzweig für `580...647 K` mit dokumentiertem Temperaturgrid berechnen und AP3 (`covers_start`, `covers_end`, `large_gaps`, `wrong_phase`, `invalid_states`, `derivative.passed`) darauf anwenden.
3. **[OFFEN]** Tatsächlich mit der neuen Methode erreichte untere T-Grenze und effektive Punktzahl messen; nicht automatisch eine bis 580 K lückenlose Kurve annehmen.
4. **[OFFEN]** Beide Nullstellen bei weiteren Temperaturen, insbesondere 590, 610, 625 und kritisch nahen Temperaturen, quantitativ verfolgen; Vorzeichenwechsel, `p/p_sat`, Zweigidentität und numerische Präzision prüfen.
5. **[OFFEN]** Root-Solver-Fehlerbehandlung, Nichtkonvergenz, Klammerbildung, sehr enge Nullstellen, exakte Nullstellensamples, Ausnahmen und physikalische Fluidgrenzen formal spezifizieren.
6. **[OFFEN]** Bei Bedarf statt jedesmal vollständigem Dichtescanning eine kontrollierte Kurvenfortsetzung vom zuvor validierten Spinodalpunkt untersuchen, ohne fragwürdige Astsprünge zu erlauben.
7. **[OFFEN]** Laufzeit/HEOS-Funktionsaufrufe für verschiedene `n_T` und `n_scan` messen und einen begründeten Default wählen.
8. **[OFFEN]** Neues Verfahren an einem weiteren repräsentativen Fluid testen; andere EOS-Modelle/Fluide können andere Mehrfachnullstellen-/Kritikalitätsverläufe besitzen.
9. **[OFFEN]** Prüfen, ob der Flüssigkeitsast weiterhin mit `build_spinodal()` bleiben soll oder ob eine einheitliche getrennte Nullstellensuche für beide Äste wissenschaftlich begründet ist; keine voreilige Ausweitung des Umfangs.
10. **[OFFEN]** Entkoppelte Optionen für **Rechen-** und **Darstellungstemperaturraster** (aktuell `spinodale_grid` Reserve, `spinodale_visu` aktiv) definieren.
11. **[OFFEN]** Renderer-Clipping mit `min_z=None` isoliert auf dem Zielrechner prüfen; Ergebnis wurde noch nicht als separater Test übermittelt.
12. **[OFFEN]** `p<=0`-Druckabschnitte und die Wahl der logarithmischen Achse wissenschaftlich getrennt dokumentieren.
13. **[OFFEN]** GUI-Integration, Fortschrittsanzeige, Threading, Abbruchverhalten, Fehlerstatus und unveränderte Sichtbarkeit der anderen Ebenen auf dem Zielrechner prüfen.
14. **[OFFEN]** Verfügbarkeit von SciPy, genaue Paketversionen und deren Reproduzierbarkeit prüfen.
15. **[OFFEN]** Bibliotheksinterne Kurventracerwerte und die genaue Version der lokalen CoolProp-Implementierung nachprüfen, bevor sie in der Studienarbeit als überprüfte Quellcodefakten zitiert werden.

---

## 11. Projektbezug zu AP1 bis AP10 (nur Einordnung, keine neue AP-Definition)

Die vorhandene Referenz `ThermoSurface_Projektanalyse_und_Arbeitspakete_AP1-AP10.md` ordnet die bisherigen Entwicklungsschritte ein:

| Bestehendes AP | Status/Relevanz zur vorliegenden Diagnose |
|---|---|
| AP1 – Tests/Regression | Bereits als Infrastruktur vorhanden, bei neuer Spinodalberechnung schützen und Baselines prüfen. |
| AP2 – p-v-T-Raster/Phasengrenzen | Bereits realisiert; unberührt lassen. |
| AP3 – Spinodalen validieren | Prüft und kennzeichnet aktuell unvollständige Äste. Die **nachträgliche Erzeugung fehlender Punkte** war ausdrücklich nicht Teil von AP3. |
| AP4 – alle fünf Diagramme validieren | Tests/Validierungslogik vorhanden; reale Zielrechnerabnahme war laut README_AP4 noch gesondert auszuführen. |
| AP5–AP10 | Nicht Gegenstand der vorliegenden Implementierung; spätere Schnittstellen-/Performance-/Referenzbetrachtungen gegebenenfalls berücksichtigen. |

Das nächste, gesondert zu formulierende Arbeitspaket wäre damit fachlich eine **Funktionsverbesserung/Ergänzung der Dampfspinodalberechnung**, nicht eine nachträgliche Behauptung, AP3 habe fehlende Daten schon erzeugt.

---

## 12. Dateiverweise, Screenshotbelege und Reproduzierbarkeit

### 12.1 Relevante bestehende Projektdateien

| Datei | Bedeutung für spätere Fortsetzung |
|---|---|
| `coolprop_heos_common.py` | `_new_state()`, Stoff-/kritische Referenzdaten, HEOS-Zustandsobjekt |
| `coolprop_heos_pvt.py` | `calc_spinodale`, `filter_spinodale`, `prepare_spinodal_branch`, p-v-T-Datenvertrag |
| `main.py` | Berechnungsdispatcher, Aufrufe, AP3-Prüfstatus, Übergabe an GUI |
| `Debug.py` | `test_spinodale_ableitung` und weitere numerische Prüfungen |
| `config.yaml` | `spinodale_grid` (Reserve), `spinodale_visu` (aktiv), Grenz-/Standardwerte |
| `gui.py` | Spinodalen-Checkbox/Status-Tooltip und Parametereingaben |
| `visualization.py` | `log10(p/pc)`-Darstellung, `_draw_curves`, `add_spinodals`, `min_z`-Clipping |
| `Test.py` | Bestehende AP1/AP2/AP3-Tests, mögliche Erweiterung |
| `README_AP3.md` | AP3-Architektur, Residualschwellen, bekannte Grenzen |
| `README_AP2.md` | Unveränderte randkonforme Phasenflächen |
| `README_AP4.md`, `test_ap4.py` | Diagrammübergreifende Tests und Zielrechnerabnahme |
| `ThermoSurface_Projektanalyse_und_Arbeitspakete_AP1-AP10.md` | Bisherige Arbeitspaketdefinition, Rahmen der Studienarbeit |

### 12.2 Original-Screenshots dieses Chats

| Originalbild | Inhalt | Wesentliche Daten hier |
|---|---|---|
| `image(20261019-074618).png` | Ausgangsplot violetter/orangefarbener Spinodalast, kritischer Punkt | Abschnitt 1.2 |
| `image(20261019-075314).png` | Rohdatendiagnose, AP3-Reports, PyVista-Warnung | Abschnitt 4 |
| `image(20261019-075930).png` | Acht Temperaturen, je zwei Ableitungs-Nullstellenintervalle | Abschnitt 5 |
| `image(20261019-080315).png` | Fünf Temperaturen, je zwei mit Brent gefundene Nullstellen samt reduzierten Drücken und Residuen | Abschnitt 6 |

Die Skripte der **beiden separaten numerischen Testläufe** sind vollständig in Abschnitt 5.2 und 6.2 enthalten. Der vorangegangene Diagnose-Schnipsel steht vollständig in Abschnitt 4.1. Der später vorgeschlagene, **noch nicht real ausgeführte** Produktionsprototyp steht vollständig in Abschnitt 8.5. Es werden somit Testcode und Vorschlags-Code bewusst getrennt überliefert.

### 12.3 Umgebungsstatus

- Die Tests wurden vom **Nutzer auf seinem Windows-System** unter `Python313/python.exe` ausgeführt; die Screenshots stammen aus PowerShell/VS Code.
- Der Realtest zeigt, dass die verwendete lokale Umgebung CoolProp ausführen und beim Nullstellentest SciPy/`brentq` aufrufen konnte.
- Versionsnummern von CoolProp/SciPy/NumPy/PyVista fehlen. Keine weitere Laufzeit-/Speicher- oder API-Versionserhebung erfolgt.
- Keine lokalen Messergebnisse für Argon, CO2, R134a oder andere Fluide liegen für diese Spinodalverlängerung vor.
- Keine realen GUI-Screenshots mit bereits **verlängertem** dampfseitigen Ast wurden eingereicht.
- Keine Änderung an `coolprop_heos_pvt.py`, `main.py`, `config.yaml`, `gui.py` oder `visualization.py` wurde durch diesen Chat als **erfolgreich umgesetzt** nachgewiesen.

---

## 13. Kurzfazit für den Start des nächsten Chats

**Problem:** Die in ThermoSurface für Wasser aus `build_spinodal()` verwendeten Rohdaten beginnen für den dampfseitigen Spinodalast erst bei rund `625,128 K`. Im eingestellten Temperaturfenster `580...647 K` verbleiben nur `13` Punkte; der AP3-Geometriecheck meldet `covers_start=False`. Der Flüssigkeitsast enthält `43` AP3-Punkte, darunter `17` mit nichtpositivem Druck, die die logarithmische Darstellung behindern. Die AP3-Druckableitungsprüfungen für beide Äste sind bestanden.

**Experimenteller Befund:** Der eigene Dichtescan mit HEOS findet bei `580`, `590`, `600`, `610`, `620`, `625`, `630` und `640 K` je **zwei** Vorzeichenwechselintervalle bei `500/500` gültigen Ableitungswerten. Eine anschließende `brentq`-Bestimmung bei `580`, `600`, `620`, `630` und `640 K` identifiziert jeweils eine Nullstelle `+ nach -` mit `p_spin > p_sat` und sehr kleinen Residuen. Diese stellt einen plausiblen Ansatz zur Berechnung des fehlenden Dampfastes dar. Eine zweite Nullstelle `- nach +` mit problematischen Druckwerten darf nicht versehentlich verwendet werden.

**Potenzielle Lösung:** Eine zusätzliche temperaturgesteuerte Dampfspinodalroutine in `coolprop_heos_pvt.py`, die dieselbe HEOS-Zustandsgleichung mittels Ableitung und `brentq` verwendet und vorläufig den bestehenden Flüssigkeitsast sowie die AP3-Test-/GUI-/Renderer-Schnittstellen unverändert lässt. Alle in diesem Chat vorgeschlagenen Quellcodeänderungen sind **Entwürfe**, keine implementierten Ergebnisse.

**Nicht nachgewiesen:** Punktweiser Vergleich mit dem vorhandenen HEOS-Dampfast im gemeinsamen Temperaturbereich, komplette Konvergenz/Temperaturabdeckung des neuen Solvers, Laufzeit, Robustheit für andere Fluide und reale GUI-Darstellung des verlängerten Astes.

**Übergabeanforderung:** Im nächsten Chat erst den dort tatsächlich verfügbaren vollständigen Quellcode und diese Dokumentation prüfen; auf dieser Grundlage **ein neues, spezifisches Arbeitspaket mit Scope, Akzeptanzkriterien, Tests und Dateiauswirkungen definieren**, ohne Testergebnisse zu erfinden oder bereits abgeschlossene AP1–AP4-Arbeiten zu verlieren.
