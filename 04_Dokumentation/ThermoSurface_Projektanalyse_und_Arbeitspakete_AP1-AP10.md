# ThermoSurface – Vollständige Projektanalyse und Arbeitspakete AP1 bis AP10

**Stand der Analyse:** 18. Oktober 2026  
**Projekt:** Interaktive 3D-Visualisierung thermodynamischer Zustandsdiagramme in Python  
**Technologie:** CoolProp HEOS, PyVista/VTK, PyVistaQt, PySide6, NumPy  
**Zweck dieser Datei:** Dauerhafte Referenz für die Studienarbeit und künftige ChatGPT-Unterhaltungen. Die Datei dokumentiert die vollständige vorherige technische und wissenschaftliche Projektbewertung und definiert insbesondere die Bezeichnungen **AP1 bis AP10** verbindlich für die Projektkommunikation.

> **Hinweis für künftige Chats:** Wenn der Nutzer „AP1“, „AP2“ usw. schreibt, sind die unten beschriebenen Arbeitspakete dieses Dokuments gemeint. Prüfe **vor jeder Implementierung** den dann aktuellen Stand sämtlicher betroffener Projektdateien. Diese Analyse bildet den Stand vom 18.10.2026 ab; Funktionen, Fehler und Dateizeilen können sich später ändern. Schlage bereits erledigte Arbeiten nicht erneut vor. Ohne ausdrückliche Aufforderung zunächst nur analysieren oder eine Schritt-für-Schritt-Anleitung liefern, statt die Projektdateien eigenmächtig umzuschreiben. Die ursprüngliche Priorität liegt auf wissenschaftlicher Korrektheit, Validierung und einem begrenzten Arbeitsumfang, nicht auf zusätzlichen kosmetischen Funktionen.

## Schnellreferenz: Bedeutung von AP1 bis AP10

| AP | Bezeichnung | Geplante Aufgabe | Betroffene Module (Schwerpunkt) | Budget | Verbindliches Abnahmekriterium |
|---|---|---|---|---:|---|
| **AP1** | **Testbasis und Regression** | Vorhandene Testskripte zu ausführbaren, aussagekräftigen Tests ausbauen; Referenzdaten, Versionen und fehlende Baselines absichern | `Test.py`, `compare_pvt.py`, `snapshot_pvt.py`, `Debug.py` | **6 h** | Testlauf gibt nachvollziehbare Ergebnisse aus und schlägt bei fehlenden Referenzdaten sicher fehl. |
| **AP2** | **p-v-T-Raster und Phasengrenzen** | Temperatur- und Sättigungsraster vereinheitlichen; Sättigungsvolumina als exakte Randpunkte in Teilflächen integrieren; geometrische Kontinuität absichern | `main.py`, `coolprop_heos_pvt.py`, `visualization.py` | **11 h** | Siede-/Taulinien stimmen mit Flächenrändern überein; zulässige Rasteränderungen erzeugen keine fehlenden Sättigungszeilen. |
| **AP3** | **Spinodalen validieren** | Ableitungsresiduen, Temperaturfilter, Endpunkte, Phasenvorgaben und Plausibilität quantitativ prüfen | `coolprop_heos_pvt.py`, `Debug.py`, ggf. `main.py` | **7 h** | Normierte Spinodalresiduen gegen begründete Toleranzen geprüft; unvollständige Äste werden als solche kenntlich. |
| **AP4** | **Alle fünf Diagramme validieren** | p-v-T, p-h-T, p-s-T, h-s-p, u-v-T mit mehreren Fluiden und kritischen/Phasenrandfällen testen | `coolprop_heos_pvt.py`, `_pht.py`, `_pst.py`, `_hsp.py`, `_uvt.py`, GUI/Renderer | **7 h** | Jeder Diagrammtyp wird für mindestens zwei repräsentative Fluide numerisch überprüft; GUI- und Renderpfade am Zielrechner getestet. |
| **AP5** | **p-h-T-Isolinien integrieren** | Vorhandene Isothermen-, Isobaren- und Isochorenfunktionen in Ergebnis-Dictionary, Renderer und GUI einbinden | `coolprop_heos_pht.py`, `visualization.py`, `gui.py`, ggf. `main.py` | **5 h** | Die drei Isolinienklassen sind im p-h-T-Diagramm separat einblendbar. |
| **AP6** | **Thermodynamische Referenzvalidierung** | Ausgewählte Zustände, Sättigungsdaten, kritische Größen, Clapeyron-Gleichung mit unabhängigen Referenzdaten und Gleichungen prüfen | Berechnungs-/Debugmodule, Auswerteskripte, Dokumentation | **10 h** | Quantitative Vergleichstabelle mit Fehlergrößen, Quellen und Diskussion liegt vor. |
| **AP7** | **HEOS-vdW-Modellvergleich** | Begrenzte eigene van-der-Waals-EOS zur Gegenüberstellung homogener einphasiger Zustände implementieren und auswerten | Separates Vergleichsskript/-modul, bestehende HEOS-Auswertung | **12 h** | Druckabweichungen und kritisches Verhalten für identische einphasige Zustände reproduzierbar dargestellt; Zweiphasengrenzen korrekt eingeordnet. |
| **AP8** | **Konvergenz und Performance** | Drei oder mehr Rasterdichten im Hinblick auf Rechenzeit, Speicher und geometrischen bzw. Interpolationsfehler untersuchen | Berechnungsaufrufe, Profiling/Auswerteskripte | **6 h** | Messreihe mit begründeter Empfehlung für die Standardauflösung. |
| **AP9** | **Schriftliche Studienarbeit** | Theorie, numerische Methoden, Architektur, Ergebnisse, Diskussion, Grenzen und Fazit wissenschaftlich dokumentieren | Studienarbeit/Abbildungen/Dokumentation | **10 h** | Vollständiger Text mit messbaren Ergebnissen, Methoden und kritischer Diskussion. |
| **AP10** | **Endprüfung und Reproduzierbarkeit** | Letzter Gesamtfunktionstest, Testmatrix, Ergebnis- und Versionssicherung, Anleitung zur Reproduktion | Gesamtprojekt, Dokumentation, Testartefakte | **2 h** | Ergebnisse sind auf einem passend eingerichteten zweiten System reproduzierbar. |

**Geplanter Umfang:** Woche 1 = AP1 bis AP5 = **36 h**; Woche 2 = AP6 bis AP10 = **40 h**; zusammen **76 h**. Die Pakete sind als Plan **ab dem analysierten Stand** zu verstehen, nicht als automatisch zusätzlich zu bereits geleisteter Projektarbeit erforderliche Zeit. **Minimaler Abschluss: ca. 59 h**, erweiterter Abschluss einschließlich AP5 und AP7: **76 h**. AP5 und AP7 können bei begrenztem Budget entfallen. Die übrigen Arbeiten sind im Grundsatz der wissenschaftlichen und technischen Absicherung zugeordnet.

---

# Vollständige technische und wissenschaftliche Gesamtbewertung – ThermoSurface

## 1. Zusammenfassung der Ergebnisse

**Meine Gesamtbewertung: Dein Projekt besitzt bereits einen erheblichen technischen Umfang und ist als Grundlage einer zweiwöchigen Master-Studienarbeit sehr gut geeignet. Der größte zusätzliche Nutzen entsteht nicht durch weitere Diagrammtypen, sondern durch die Überprüfung der thermodynamischen Korrektheit, die Beseitigung einzelner struktureller Schwachstellen und eine systematische wissenschaftliche Auswertung.**

Die Codebasis ist bereits deutlich über einen einfachen Demonstrator hinaus entwickelt. Die Anwendung enthält:

- Fünf thermodynamische 3D-Diagrammtypen mit eigenen Berechnungsmodulen.
- Eine modulare Programmarchitektur mit PySide6-GUI, PyVista-Visualisierung und CoolProp-HEOS-Berechnung.
- Siede- und Taulinien, kritische Punkte, Zweiphasenflächen und teilweise Isolinien.
- Die Berechnung von Spinodalen im p-v-T-Diagramm.
- Automatische stoffabhängige Raster, Berechnung im Hintergrundthread und Fortschrittsanzeige.
- Mehrere bereits integrierte numerische Konsistenzprüfungen.

Allerdings sind nicht alle Funktionen im gleichen Entwicklungszustand. Es gibt insbesondere **zwei reproduzierbare Defizite beim p-v-T-Diagramm**, Unterschiede in der Integration der Isolinien und eine noch unzureichende unabhängige Validierung der thermodynamischen Ergebnisse.

### Bewertungsübersicht

| Bewertungsbereich | Qualitative Einschätzung | Begründung |
|---|---:|---|
| Programmarchitektur | **4/5** | Modular aufgebaut, grundsätzlich gut erweiterbar |
| Funktionsumfang | **4/5** | Fünf Diagrammtypen, zahlreiche Darstellungsebenen |
| Numerische Robustheit | **3/5** | Viele Prüfungen vorhanden, kritische Fälle offen |
| Thermodynamische Validierung | **2/5** | Konsistenztests vorhanden, unabhängige Validierung fehlt |
| Automatisierte Tests | **1/5** | Keine ausführbare vollständige Testsuite nachgewiesen |
| Wissenschaftliche Dokumentation | **1/5** | In den vorliegenden Projektdateien nicht ausgearbeitet |

*Qualitative Einschätzung anhand der Codeinspektion und der möglichen Funktionstests. Die Skala ist keine formale Benotung oder Messung der Laufzeitzuverlässigkeit.*

**Wesentliche Schlussfolgerung:** Für einen überzeugenden Studienarbeitsabschluss würde ich die vorhandenen fünf Diagramme nicht um weitere Typen ergänzen. Stattdessen würde ich einen wissenschaftlichen Schwerpunkt festlegen: die thermodynamisch konsistente 3D-Darstellung von Phasengleichgewichten und Stabilitätsgrenzen auf Basis von HEOS, ergänzt durch Validierung und gegebenenfalls einen begrenzten Vergleich mit van der Waals.

## 2. Umfang und Aussagekraft der durchgeführten Codeprüfung

Die damalige Analyse hat alle 14 vorhandenen Python-Dateien sowie `config.yaml` auf Syntax, Funktionsstruktur und modulübergreifende Zusammenhänge untersucht. Der Python-Code wurde vollständig eingelesen; für die fachliche Prüfung wurden die wesentlichen Berechnungs- und Visualisierungspfade gesondert nachvollzogen.

### Prüfergebnisse zum damaligen Stand

- **Bestanden:** Python-Syntaxprüfung aller Module.
- **Bestanden:** Konfiguration nach der in `main.py` vorgesehenen Typkonvertierung erfolgreich validiert.
- **Bestanden:** Fünf synthetische Einzelprüfungen der Raster-, Mischungs- und Datenumwandlungsfunktionen.
- **Reproduziert:** Zwei Fehlerfälle in der p-v-T-Verarbeitung mit synthetischen Daten.
- **Nicht ausführbar:** `Test.py`, da CoolProp in der damaligen Testumgebung fehlte.
- **Nicht nachgewiesen:** Von pytest automatisch erkannte Tests.

Die Installation von CoolProp war in der damaligen Umgebung nicht möglich. Auch PyVista, PyVistaQt und PySide6 standen dort nicht zur Verfügung. Deshalb konnten keine vollständigen HEOS-Berechnungsläufe, grafischen Interaktionstests oder Performance-Benchmarks durchgeführt werden.

Die folgenden Aussagen unterscheiden daher zwischen **statisch nachgewiesener Implementierung**, **synthetisch überprüfter Logik**, **reproduzierbaren Fehlern** und **noch nicht praktisch verifizierter Funktionalität**.

Wichtig ist diese Trennung besonders für die wissenschaftliche Bewertung: Dass eine Funktion die passende CoolProp-Schnittstelle verwendet, beweist noch nicht, dass sie für alle auswählbaren Fluide oder im kritischen Bereich zuverlässig funktioniert.

## 3. Analyse der bestehenden Programmarchitektur

### 3.1 Aufbau und Zuständigkeiten

Vereinfachte funktionale Architektur:

```text
                              main.py
                   (Start, Konfiguration, Dispatcher)
                             /       \
                            /         \
               gui.py                     visualization.py
          (PySide6, QThread)               (PyVista / VTK)
                            \         /
                             \       /
                  Berechnungsmodule nach Diagrammtyp
            coolprop_heos_pvt / pht / pst / hsp / uvt
                            /        \
                           /          \
              coolprop_heos_common     config.yaml / Debug.py
                 (CoolProp HEOS)          (Parameter, Diagnosen)
```

*Die tatsächlichen Modulaufrufe verlaufen insbesondere über die Berechnungssteuerung in `main.py`; der Renderer erhält bereits berechnete Zustandsdaten.*

| Modul | Bestehende Aufgabe | Bewertung |
|---|---|---|
| `main.py` | Programmstart, Konfiguration, Auswahl und Aufruf der Berechnungen | Gut |
| `gui.py` | Stoffauswahl, Bedienung, QThread, Darstellungssteuerung | Gut |
| `visualization.py` | Transformation und Vernetzung von Zustandsdaten, Darstellung | Gut, teilweise inkonsistent |
| `coolprop_heos_common.py` | Gemeinsame HEOS-Zugriffe und kritische Stoffdaten | Gut |
| `coolprop_heos_pvt.py` | p-v-T, Isolinien, Sättigung, Spinodalen, Zweiphasenmaske | Umfangreich, Schwachstellen |
| `coolprop_heos_pht.py` | p-h-T-Flächen und Phasengrenzen | Gut strukturiert, Isolinien nicht integriert |
| `coolprop_heos_pst.py` | p-s-T-Flächen, Phasengrenzen, Isolinien | Implementiert |
| `coolprop_heos_hsp.py` | h-s-p-Daten durch Zusammenführung von p-h-T und p-s-T | Sinnvoll wiederverwendet |
| `coolprop_heos_uvt.py` | u-v-T-Zustandsflächen und Isolinien | Implementiert |
| `Debug.py` | Diagnostik, Plausibilitätsprüfungen, Laufzeitmessung | Sinnvoll, enthält Fehler |
| `limits.py` | Zentrale Punktzahlbegrenzung | Zweckmäßig |
| `config.yaml` | Stoff, Raster, Diagrammeinstellungen | Teilweise uneinheitlich |
| `Test.py`, `snapshot_pvt.py`, `compare_pvt.py` | Importprüfungen und Regressionstest-Werkzeuge | Nicht ausreichend abgesichert |

### Architektururteil

Die bestehende Trennung zwischen Berechnung, GUI und Visualisierung ist zweckmäßig. Besonders positiv ist, dass die Berechnungsfunktionen der zusätzlichen Diagrammtypen überwiegend ohne Abhängigkeiten von PyVista oder PySide6 auskommen.

`main.berechne_auswahl()` (Zeilen 333-392 im damaligen Stand) fungiert bereits als zentraler Dispatcher, der nur den ausgewählten Diagrammtyp berechnet. Eine vollständige Umstellung auf eine neue Klassenarchitektur, einen Plugin-Manager oder ein anderes GUI-Framework wäre für die Studienarbeit nicht gerechtfertigt.

Verbesserungsbedarf besteht eher bei den **Schnittstellen und gemeinsamen Datenkonventionen**. Beispielsweise verwenden die Diagrammtypen unterschiedliche Reihenfolgen für ihre kritischen Punkte:

| Diagramm | Punktdaten | Kritischer Punkt |
|---|---|---|
| p-v-T | `(T, v, p)` | `(Tc, pc, vc)` |
| p-h-T | `(T, h, p)` | `(Tc, pc, hc)` |
| p-s-T | `(T, s, p)` | `(Tc, pc, sc)` |
| h-s-p | `(h, s, p)` | `(hc, sc, pc)` |
| u-v-T | `(T, v, u)` | `(Tc, vc, uc)` |

Die Implementierungen berücksichtigen diese Unterschiede überwiegend korrekt. Trotzdem entstehen dadurch Risiken bei generischen Tests, Debugging und zukünftigen Erweiterungen.

**Empfehlung:** Die Formate verbindlich dokumentieren und automatisiert prüfen. Eine umfangreiche Umstrukturierung ist dafür nicht erforderlich.

---

## 4. Tatsächlicher Funktionsstand der fünf Diagrammtypen

| Funktion | p-v-T | p-h-T | p-s-T | h-s-p | u-v-T |
|---|:---:|:---:|:---:|:---:|:---:|
| Zustandsfläche | ✓ | ✓ | ✓ | ✓ | ✓ |
| Siede-/Taulinie | ✓ | ✓ | ✓ | ✓ | ✓ |
| Zweiphasenfläche | ✓ | ✓ | ✓ | ✓ | ✓ |
| Kritischer Punkt | ✓ | ✓ | ✓ | ✓ | ✓ |
| Isothermen | ✓ | – | ✓ | ✓ | ✓ |
| Isobaren | ✓ | – | ✓ | ✓ | ✓ |
| Isochoren | ✓ | – | ✓ | ✓ | ✓ |
| Spinodalen | ✓ | – | – | – | – |

*✓ = Berechnung und Integration im Code vorhanden, nicht notwendigerweise praktisch funktionsverifiziert. – = Nicht in den regulären Berechnungs- und Darstellungsablauf eingebunden.*

### 4.1 p-v-T: funktional am weitesten ausgebaut

Das bestehende Modul bietet mit `multi_isotherm`, `multi_isochore`, `multi_isobare`, `multi_binodale`, `calc_spinodale` und `calc_two_phase_surface` bereits den größten Umfang.

Die GUI bindet sämtliche wesentlichen Ebenen ein, einschließlich einer Checkbox für das Oberflächengitter. Ein erneuter Ausbau dieser allgemeinen Bedienelemente ist deshalb nicht notwendig.

Allerdings ist die geometrisch randkonforme Oberflächenkonstruktion noch nicht für beliebige zulässige Parameterkonfigurationen abgesichert.

### 4.2 p-h-T: Flächen vorhanden, Isolinien noch nicht durchgängig integriert

In `coolprop_heos_pht.py` sind die drei Funktionen

- `calc_isotherme()` – Zeilen 364-388,
- `calc_isobare()` – Zeilen 391-427,
- `calc_isochore()` – Zeilen 430-456

bereits implementiert.

`berechne_pht_daten()` gibt allerdings nur die Flächen, Phasengrenzen und Metadaten zurück. Die Isolinien sind dort nicht Bestandteil des regulären Ergebnis-Dictionaries. `render_pht()` zeichnet diese ebenfalls nicht. Die Datenrückgabe befindet sich im damaligen Stand im Abschnitt von `coolprop_heos_pht.py` um Zeilen 333-359.

Das ist eine tatsächliche, überschaubare Integrationslücke. Die Berechnungslogik muss dafür nicht neu entwickelt werden.

### 4.3 p-s-T: weitgehend integriert

`coolprop_heos_pst.py` erstellt randkonforme Flüssigkeits-, Dampf- und Hochtemperaturflächen. Die Isolinien werden über `_isolinien()` erzeugt und an `render_pst()` übergeben.

Das Modul verwendet vorhandene Funktionen aus p-h-T zur Rastererzeugung und Geometrieprüfung. Das vermeidet Dopplungen, bedeutet aber auch, dass Änderungen am p-h-T-Raster auf p-s-T zurückwirken können.

### 4.4 h-s-p: wissenschaftlich interessantes Diagramm mit mehrfacher Berechnung

`coolprop_heos_hsp.py` führt p-h-T- und p-s-T-Datensätze über `_merge()` zu gemeinsamen `(h,s,p)`-Punkten zusammen. Die Übereinstimmung der Temperaturen und Drücke wird dabei überprüft.

Die Isolinien werden durch `_convert_pst_curve()` unter Verwendung eines Druck-Entropie-Flashes in Enthalpiekoordinaten umgerechnet. Die grundlegende Zusammenführung der Arrays wurde in der damaligen Prüfung synthetisch erfolgreich getestet.

Ein Nachteil: `berechne_hsp_daten()` ruft sowohl die vollständige p-h-T- als auch die vollständige p-s-T-Berechnung auf. Dadurch entstehen doppelte HEOS-Auswertungen gemeinsamer Zustände (`coolprop_heos_hsp.py`, `berechne_hsp_daten()`, im damaligen Stand um Zeilen 111-145).

Für die derzeitigen Rastergrößen ist das nicht zwingend ein Problem, sollte aber in der Performanceuntersuchung berücksichtigt werden.

### 4.5 u-v-T: weitgehend integriert

`coolprop_heos_uvt.py` erstellt die Zustandsflächen über `(T,v)`-Zustände und berechnet die innere Energie aus HEOS.

Besonders sinnvoll ist die explizite Berechnung der Zweiphasenfläche mit dem Dampfgehalt `Q`:

\[
v=(1-Q)v_l+Qv_g
\]

\[
u=(1-Q)u_l+Qu_g
\]

Diese Beziehungen sind für das Zweiphasengleichgewicht eines Reinstoffs bei einer gegebenen Sättigungstemperatur korrekt.

`validate_uvt_data()` überprüft außerdem, ob Flüssigkeits-, Dampf- und Zweiphasenflächen gemeinsame Randpunkte besitzen.

**Einordnung:** p-s-T, h-s-p und u-v-T sind nicht mehr bloß vorbereitete Dateien, sondern besitzen vollständige Berechnungs- und Renderpfade. Was fehlt, sind belastbare praktische Nachweise über mehrere Fluide hinweg.

---

## 5. Konkrete technische und thermodynamische Schwachstellen

### 5.1 Hohe Priorität: Inkonsistente p-v-T-Stützstellen

**Status in der damaligen Analyse: Reproduzierbarer Fehlerfall.**

**Betroffene Stellen:**

- `main.py`, Zeilen 175-224,
- `visualization.py`, Funktion `create_single_phase_surface()`, Zeilen 433-572.

In `main.berechne_daten()` werden zunächst die Temperaturstützstellen der Zustandsoberfläche mit denen der Binodalen zusammengeführt:

```python
T_all = np.unique(np.concatenate(
    (T_surface, T_sat_visu, T_critical_row)
))
```

Danach wird eine Binodale auf genau diesem Raster erzeugt und für die Zweiphasenmaske verwendet.

Später wird das Ergebnis aber durch einen zweiten Aufruf von `multi_binodale()` mit dem separaten `binodale_visu`-Raster überschrieben (`main.py`, damals Zeilen 199-224).

Die Funktion `create_single_phase_surface()` verlangt für jede unterkritische Temperaturzeile einen zugehörigen Sättigungspunkt.

Fehlt dieser im später übergebenen Binodalarray, wird folgende Ausnahme ausgelöst:

```text
ValueError: Kein Sättigungszustand für T=... K.
```

Dieser Fehler konnte mit einem kleinen synthetischen Datensatz reproduziert werden.

Er tritt nicht zwangsläufig bei den aktuellen Standardparametern auf, da dort viele Temperaturstützstellen zusammenfallen. Bei veränderten Rasterabständen ist die Übereinstimmung jedoch nicht garantiert.

**Notwendige Verbesserung:** Ein eindeutig definiertes gemeinsames Temperatur- und Sättigungsraster als Quelle für Phasengrenzen und Teilflächen verwenden, oder eine explizit validierte Zuordnung schaffen.

Zusätzlich sollten die exakten Sättigungsvolumina `v_l(T)` und `v_g(T)` Randknoten der einphasigen Flächen bilden.

Damit lässt sich nicht nur ein Softwarefehler vermeiden, sondern auch die geometrische Übereinstimmung der Phasengrenzen mathematisch nachweisen.

### 5.2 Hohe Priorität: Fehler in der Debug-Visualisierung

**Status in der damaligen Analyse: Reproduzierbarer Fehlerfall.**

**Betroffene Stellen:**

- `main.py`, Zeilen 240-297,
- `Debug.py`, Zeilen 303-335,
- `gui.py`, Zeilen 539-545.

Im p-v-T-Ablauf wird die Zweiphasenfläche zunächst ohne kritische Abschlusszeile erzeugt. Erst anschließend ergänzt `ergaenze_kritischen_punkt()` den kritischen Punkt zu beiden Binodalarrays.

`Debug.diagnose_visualisierung()` vergleicht jedoch die erste beziehungsweise letzte Spalte der Zweiphasenfläche direkt mit den gesamten Binodalarrays (`Debug.py`, damals Zeilen 315-322).

Damit unterscheiden sich ihre Dimensionen:

\[
N_{\mathrm{Binodale}}=N_{\mathrm{Zweiphasenfläche}}+1
\]

Diese Struktur wurde im synthetischen Test nachgebildet. Der Vergleich schlägt mit einem `AssertionError` fehl.

Bei aktiviertem `DEBUG_MODE` kann dies den ansonsten erfolgreich berechneten Darstellungsablauf unterbrechen.

**Notwendige Verbesserung:** Die regulären Randzeilen ohne kritischen Endpunkt vergleichen und den kritischen Abschluss getrennt prüfen.

### 5.3 Hohe wissenschaftliche Priorität: Spinodalen sind nicht ausreichend validiert

**Betroffene Stellen:**

- `coolprop_heos_pvt.py`, `calc_spinodale()`, Zeilen 586-691,
- `Debug.py`, `test_spinodale_ableitung()`, Zeilen 269-300,
- `visualization.py`, `add_spinodals()`, Zeilen 351-358.

Die Methode verwendet bereits die CoolProp-Funktionen `build_spinodal()` und `get_spinodal_data()`.

Die Spinodalpunkte werden anschließend über reduzierte Temperatur- und Dichteparameter rekonstruiert. Durch die Phasenvorgabe versucht der Code, den Druck der homogenen, gegebenenfalls metastabilen Zustandsfortsetzung zu bestimmen.

Dieser Ansatz ist grundsätzlich nachvollziehbar. Die wissenschaftliche Prüfung ist allerdings noch unvollständig.

Für einen Spinodalpunkt muss gelten:

\[
\boxed{\left(\frac{\partial p}{\partial v}\right)_T=0}
\]

Mit der molaren Dichte `rho_m = 1/v` ist dies äquivalent zu:

\[
\left(\frac{\partial p}{\partial \rho_m}\right)_T=0.
\]

Dein Debug-Modul berechnet bereits ein entsprechendes normiertes Ableitungsresiduum. Es protokolliert jedoch hauptsächlich Kennwerte, statt die Ergebnisse mit klaren numerischen Akzeptanzgrenzen zu bewerten.

Weitere Schwachstellen sind die fehlende Garantie auf durchgehend gültige Spinodaläste und die Tatsache, dass die grafische Transformation ausschließlich positive Druckwerte erlaubt.

**Wichtig für die thermodynamische Interpretation:** Die Spinodale muss nicht auf der stabilen Gleichgewichtsfläche liegen. Innerhalb der Binodalen zeigt HEOS bei der stabilen Gleichgewichtsauswertung das Zweiphasenplateau, während die Spinodale den homogenen Stabilitätsgrenzzustand beschreibt.

Das ist kein Visualisierungsfehler, sondern ein Unterschied der dargestellten thermodynamischen Zustände.

Außerdem fügt `main.py` den kritischen Punkt nach dem Temperaturfilter erneut an. Dadurch kann die angezeigte Spinodale die eingestellte obere Temperaturgrenze überschreiten. Auch das sollte explizit geprüft werden.

### 5.4 Mittlere Priorität: Punktzahlbegrenzung ohne systematisches Ressourcenmanagement

`limits.py` übernimmt aus `config.yaml` die Grenze:

\[
N_{\max}=1\,000\,000
\]

Diese Grenze wird in mehreren Modulen kontrolliert. Die Implementierung ist also nicht mehr grundsätzlich auf 10.000 Punkte beschränkt.

In `coolprop_heos_pvt.validate_rasters()` existieren noch lokale Zuweisungen `max_points = 10000`, die jedoch nicht verwendet werden. Sie stellen in der aktuellen Version keine wirksame Begrenzung dar.

Das eigentliche Problem ist die Speicher- und Rechenlast bei hohen Auflösungen.

**Speicherbedarf eines einzelnen NumPy-Arrays `(N,3)` mit float64:**

| Zustandspunkte | Rohdatenspeicher |
|---:|---:|
| 10.000 | ca. 0,229 MiB |
| 100.000 | ca. 2,289 MiB |
| 500.000 | ca. 11,444 MiB |
| 1.000.000 | ca. 22,888 MiB |

*Nur Rohdaten, ohne zusätzliche Arrays und VTK-Overhead.*

Eine Million Zustandspunkte benötigen als reines `(N,3)`-Array etwa 22,9 MiB. In der Anwendung können jedoch mehrere Arrays, transformierte Koordinaten, Kopien und die PyVista-/VTK-Geometrie gleichzeitig existieren.

Deshalb ist das Speichermaximum des gesamten Programms deutlich höher als das eines einzelnen Arrays. Eine belastbare Gesamtschätzung oder Rechenzeit kann ohne Profiling nicht angegeben werden.

Die fachlich sinnvollere Maßnahme ist eine Rasterkonvergenzanalyse anstelle einer weiteren Erhöhung der Punktzahlgrenze.

### 5.5 Mittlere Priorität: Vorhandene Tests sichern noch keine vollständige Reproduzierbarkeit

**Betroffene Dateien:** `Test.py`, `snapshot_pvt.py`, `compare_pvt.py`.

`Test.py` überprüft im Wesentlichen Modulzuordnungen und die Existenz bestimmter Funktionen. Es berechnet keine thermodynamischen Referenzzustände.

Die Snapshot-Dateien sind als Regressionstest-Werkzeuge sinnvoll angelegt. Allerdings befinden sich im bereitgestellten Projekt keine gespeicherten Baseline-Datensätze.

`compare_pvt.py` wurde zusätzlich ohne vorhandene Baseline ausgeführt: Das Skript beendete sich mit Exit-Code 0, obwohl tatsächlich kein einziger Datensatz verglichen wurde. Die Ursache ist die Iteration über einen möglicherweise leeren Dateibestand (damals `compare_pvt.py`, Zeilen 6-15).

**Notwendige Verbesserung:** Der Test muss einen fehlenden oder leeren Referenzdatenbestand ausdrücklich als Fehler behandeln. Außerdem sollten neben den Vergleichsdaten die verwendeten Bibliotheksversionen und Berechnungsparameter nachvollziehbar gespeichert werden.

### 5.6 Weitere festgestellte Einschränkungen

| Beobachtung | Konkrete Stelle | Einordnung |
|---|---|---|
| Zweifache kritische Stoffdatenbestimmung im p-v-T-Ablauf | `main.py`, `berechne_daten()` | Kleine Performanceineffizienz |
| Teilweise doppelte Sättigungsberechnungen | `main.py`, `berechne_daten()` | Optimierung sinnvoll |
| Erneute Berechnung gemeinsamer Zustände in h-s-p | `coolprop_heos_hsp.py`, `berechne_hsp_daten()` | Messbarer Optimierungskandidat |
| Keine allgemeine Ergebnisexportfunktion | Berechnungssteuerung und GUI | Optional, nützlich für Validierung |
| Keine Hover- oder Punktabfrage im 3D-Viewport | `gui.py`, `visualization.py` | Optional |
| Diagrammspezifische Einstellungen nur teilweise über GUI zugänglich | `gui.py`, `config.yaml` | Komforteinschränkung |
| Unterschiedliche Schattierungseinstellungen in den Renderern | `visualization.py` | Kosmetisch |
| Keine eigenständige kubische Zustandsgleichung in den vorliegenden Python-Dateien | Projektbestand | Nicht implementiert |

**Prioritätensetzung:** Weder Hover-Anzeigen noch EXE-Erstellung, zusätzliche Schattierung oder weitere Diagrammtypen sind derzeit entscheidend für den wissenschaftlichen Projektabschluss.

---

## 6. Wissenschaftliche Bewertung

### 6.1 Softwareentwicklung gegenüber Ingenieurarbeit

Hier liegt die entscheidende Frage für die Studienarbeit.

Die Verwendung von CoolProp, PyVista und PySide6 ist grundsätzlich kein Problem. Wissenschaftliche Software darf und soll etablierte Bibliotheken verwenden. Entscheidend ist, welche eigenständige Leistung darüber hinaus erbracht wird.

| Bestandteil | Softwaretechnische Leistung | Wissenschaftlicher Beitrag |
|---|---|---|
| GUI mit Dropdown und Checkboxen | Hoch | Gering |
| Erstellung mehrerer Zustandsdiagramme | Hoch | Mittel |
| Automatische Rastererzeugung um kritische Zustände | Mittel | Mittel bis hoch |
| Geometrisch konsistente Phasengrenzen | Hoch | Hoch |
| Berechnung und Interpretation der Spinodalen | Mittel | Hoch |
| Numerische Konvergenzuntersuchung | Mittel | Sehr hoch |
| Thermodynamische Validierung | Mittel | Sehr hoch |
| Vergleich HEOS gegen van der Waals | Mittel | Sehr hoch |

Der bisherige Code demonstriert insbesondere eine umfangreiche **eigenständige Integration thermodynamischer Berechnungsverfahren in eine interaktive Visualisierungsanwendung**.

Er beweist bislang noch nicht ausreichend, mit welcher Genauigkeit, Robustheit und numerischen Effizienz die Ergebnisse erzeugt werden.

Genau darin liegt die sinnvollste wissenschaftliche Ergänzung.

### 6.2 Darstellung stabiler, metastabiler und instabiler Zustände

Hier muss zwischen den bereits implementierten Zustandsflächen und den gewünschten Zustandsbereichen unterschieden werden.

Die Anwendung berechnet primär thermodynamische Gleichgewichtszustände.

Schematisch unterscheidet sich für eine unterkritische Isotherme die homogene EOS-Kurve mit eventueller van-der-Waals-artiger Schleife von einem Gleichgewichtsplateau bei Sättigungsdruck. Die konkrete Form ist modellabhängig. Diese schematische Erläuterung ist **keine numerische Ausgabe** des Programms.

Die Unterscheidung lautet:

- **Stabiler Bereich:** Thermodynamisch bevorzugter Gleichgewichtszustand.
- **Metastabiler Bereich:** Lokal stabiler homogener Zustand, aber nicht das globale Gleichgewicht.
- **Instabiler Bereich:** Homogene Zustände mit negativer isothermer Kompressibilität beziehungsweise `(∂p/∂v)_T > 0`.

Die Spinodale markiert die Grenze zwischen metastabilen und instabilen homogenen Zuständen. Die Binodale beschreibt dagegen das Phasengleichgewicht.

Die aktuelle p-v-T-Darstellung besitzt die Binodale, die Gleichgewichts-Zweiphasenfläche und gesondert berechnete Spinodalen. **Eine vollständige metastabile und instabile Zustandsfläche ist dagegen nicht implementiert.**

Diese Ergänzung sollte für eine zweiwöchige Studienarbeit nicht verpflichtend aufgenommen werden. Für eine nachvollziehbare wissenschaftliche Darstellung reicht zunächst:

1. Die stabile Gleichgewichtsfläche zeigen.
2. Binodale und Spinodale separat berechnen.
3. Den Unterschied thermodynamisch erklären.
4. Die Stabilitätsbedingung anhand ausgewählter Spinodalpunkte überprüfen.

Das ist inhaltlich fundiert und begrenzt den Aufwand.

### 6.3 Wissenschaftlicher Mehrwert eines HEOS-vdW-Vergleichs

**Ein begrenzter Vergleich mit einer selbst implementierten van-der-Waals-Zustandsgleichung wäre die bevorzugte optionale wissenschaftliche Erweiterung.**

Dafür müssen weder alle fünf Diagramme doppelt berechnet noch die bestehende GUI um ein zweites Modell erweitert werden.

Die van-der-Waals-Gleichung lautet:

\[
p(T,v)=\frac{RT}{v-b}-\frac{a}{v^2}
\]

Mit kritischer Temperatur und kritischem Druck ergeben sich:

\[
a=\frac{27R^2T_c^2}{64p_c},\qquad b=\frac{RT_c}{8p_c}
\]

Damit lässt sich eine eigene Zustandsgleichung mit verhältnismäßig geringem Implementierungsaufwand erzeugen.

Ein geeigneter Vergleich untersucht beispielsweise:

| Untersuchung | Erkenntnis |
|---|---|
| Isothermen bei `T_r = 0,85`, `0,95` und `1,10` | Unterschiede unterhalb und oberhalb des kritischen Punkts |
| Druckabweichungen bei gleichen `T,v` | Modellabweichung gegenüber HEOS |
| Kritischer Kompressibilitätsfaktor `Z_c` | Grenzen der universellen vdW-Modellstruktur |
| Spinodalbedingungen | Unterschiede im Stabilitätsverhalten |
| Berechnungszeit bei gleichem Raster | Kosten der Zustandsauswertung |

Für den Vergleich kann definiert werden:

\[
\varepsilon_p=\frac{|p_{\mathrm{vdW}}-p_{\mathrm{HEOS}}|}{|p_{\mathrm{HEOS}}|}
\]

Die Auswertung ist auf sinnvolle, physikalisch vergleichbare Zustände zu begrenzen.

Ein wesentlicher Punkt: Die unmodifizierte van-der-Waals-Schleife unterhalb von `T_c` beschreibt kein korrektes Gleichgewichts-Zweiphasenplateau. Wenn die Sättigungszustände beider Modelle verglichen werden sollen, muss für van der Waals zunächst das Phasengleichgewicht bestimmt werden, beispielsweise mit dem Maxwell-Flächenkriterium.

Das kann den Aufwand merklich erhöhen.

**Deshalb empfiehlt sich ein wissenschaftlich sauberer, aber begrenzter Vergleich zunächst in den homogenen einphasigen Bereichen.** Ein zusätzlicher Maxwell-Vergleich ist sinnvoll, wenn ausreichend Zeit verfügbar bleibt.

Wichtig für die Bewertung: HEOS ist dabei das höherwertige Vergleichsmodell, aber nicht automatisch eine unabhängige experimentelle Referenz. Für eine wirkliche Modellvalidierung werden zusätzlich Literatur- oder Referenzdaten benötigt.

### 6.4 Welche numerischen Untersuchungen fehlen?

Den höchsten Mehrwert besitzt eine Konvergenzuntersuchung.

Beispielsweise könnte die p-v-T-Fläche mit drei Rasterauflösungen berechnet werden:

| Raster | Stützstellen | Gesamtpunkte |
|---|---:|---:|
| A | 21 × 41 | 861 |
| B | 41 × 81 | 3.321 |
| C | 81 × 161 | 13.041 |

*Vorgeschlagene Versuchsgrößen, keine gemessenen Ergebnisse.*

Die Auswertung sollte dabei nicht nur die Differenz der direkt von HEOS berechneten Druckwerte betrachten: An identischen Zuständen können diese auch bei unterschiedlichen Rastern nahezu gleich sein.

Interessanter ist die Abweichung der **diskretisierten beziehungsweise interpolierten Zustandsfläche** gegenüber zusätzlichen HEOS-Kontrollpunkten.

Damit können die Beziehungen zwischen Rasterauflösung, geometrischem Fehler, Rechenzeit und Speicherbedarf untersucht werden.

Ein gutes Ergebnis wäre beispielsweise die Feststellung, ab welcher Rasterdichte eine weitere Verdoppelung der Stützpunkte keinen relevanten Vorteil mehr für die geometrische Darstellung bringt.

Das ist eine klare numerisch-ingenieurwissenschaftliche Leistung.

---

## 7. Priorisierte Übersicht der verbleibenden Arbeiten

Die folgende Tabelle berücksichtigt, dass wichtige Funktionen bereits vorhanden sind.

| Arbeit | Technischer Nutzen | Wissenschaftlicher Mehrwert | Aufwand | Schwierigkeit | Priorität |
|---|---|---|---|---|---|
| p-v-T-Phasengrenzen und Raster synchronisieren | Sehr hoch | Hoch | 8-12 h | Hoch | P0 |
| Debug- und Regressionstests korrigieren | Hoch | Mittel | 5-7 h | Niedrig bis mittel | P0 |
| Thermodynamische Referenzvalidierung | Hoch | Sehr hoch | 8-12 h | Mittel | P0 |
| Alle fünf Diagramme systematisch testen | Hoch | Hoch | 6-9 h | Mittel | P0 |
| Ergebnisse und Methodik dokumentieren | Hoch | Sehr hoch | 10-14 h | Mittel | P0 |
| Spinodalen quantitativ validieren | Hoch | Sehr hoch | 6-9 h | Mittel bis hoch | P1 |
| Rasterkonvergenz und Performance untersuchen | Hoch | Sehr hoch | 4-7 h | Mittel | P1 |
| p-h-T-Isolinien in GUI/Renderer integrieren | Mittel | Mittel | 3-6 h | Mittel | P1 |
| Eigene van-der-Waals-Vergleichsrechnung | Mittel | Sehr hoch | 10-16 h | Mittel bis hoch | P1 |
| Hover-Anzeige und Datenabfrage | Mittel | Gering | 4-8 h | Mittel | P2 |
| EXE-Erstellung und GUI-Feinschliff | Mittel | Gering | 4-10 h | Mittel | P2 |

Die Aufwände sind Schätzwerte für die Arbeit mit dem bestehenden Code. Sie beruhen nicht auf gemessenen Entwicklungszeiten und sind nicht unabhängig voneinander zu addieren, da sich einzelne Prüf- und Dokumentationsaufgaben überschneiden können.

P0 bezeichnet Arbeiten, die für einen belastbaren Projektabschluss als notwendig angesehen werden. P1 sind sinnvolle Vertiefungen, wobei die Spinodalprüfung verpflichtend wird, sobald quantitative Aussagen zur Stabilitätsgrenze veröffentlicht werden. P2 sollte innerhalb des engen Zeitrahmens entfallen.

---

## 8. Konkreter Zweiwochen-Arbeitsplan

Für den Fall, dass ab dem damals untersuchten Projektstand tatsächlich noch etwa 75 Stunden verfügbar sind, empfiehlt sich folgender Plan. Er baut auf der vorhandenen Architektur auf und vermeidet eine grundlegende Neuentwicklung.

**Zeitbudget: 76 Stunden / 10 Arbeitstage à durchschnittlich 7,6 Stunden.**

- Woche 1: **36 Stunden**.
- Woche 2: **40 Stunden**.

### Woche 1: Technische und thermodynamische Absicherung

| Paket | Konkrete Umsetzung | Aufwand |
|---|---|---:|
| **AP1 – Testbasis** | `Test.py`, `compare_pvt.py`, `snapshot_pvt.py`, `Debug.py`: ausführbare Tests, Versionserfassung, fehlende Baselines als Fehler behandeln | **6 h** |
| **AP2 – p-v-T-Raster** | `main.py`, `coolprop_heos_pvt.py`, `visualization.py`: gemeinsames Sättigungsraster, exakte Randpunkte, lückenfreie Teilflächen | **11 h** |
| **AP3 – Spinodalen** | `coolprop_heos_pvt.py`, `Debug.py`: Ableitungsresiduen, Temperaturfilter, Endpunkt- und Plausibilitätsprüfungen | **7 h** |
| **AP4 – Diagrammvalidierung** | `pht`, `pst`, `hsp`, `uvt`: unterschiedliche Fluide, Randkontinuität, Einheiten und kritische Zustände testen | **7 h** |
| **AP5 – p-h-T-Isolinien** | `coolprop_heos_pht.py`, `visualization.py`, `gui.py`: vorhandene Linienfunktionen in Berechnung und Renderer einbinden | **5 h** |
| **Gesamt** | | **36 h** |

Die überprüfbaren Ergebnisse dieser Woche sollten sein:

**AP1:** Ein Testlauf liefert nachvollziehbare Ergebnisse und schlägt bei fehlenden Referenzdaten zuverlässig fehl.

**AP2:** Siede- und Taulinien stimmen für alle getesteten Temperaturzeilen mit den Rändern der einphasigen und zweiphasigen Netze überein. Änderungen der zulässigen Rasterabstände führen nicht mehr zu fehlenden Sättigungszeilen.

**AP3:** Das dimensionslose Spinodalresiduum wird für geeignete Punkte außerhalb der unmittelbaren kritischen Nähe bestimmt und gegen eine begründete Toleranz geprüft. Unvollständige Äste werden nicht als vollständig validiert ausgegeben.

**AP4:** Jeder der fünf Diagrammtypen kann mit mindestens zwei repräsentativen Fluiden numerisch geprüft werden. Für GUI und Renderer ist danach auf dem Zielrechner ein praktischer Funktionstest erforderlich.

**AP5:** p-h-T-Isothermen, Isobaren und Isochoren können zusätzlich zu den vorhandenen Phasengrenzen eingeblendet werden.

AP5 ist optional, wenn der verbleibende Gesamtaufwand begrenzt werden muss.

### Woche 2: Wissenschaftliche Auswertung und Dokumentation

| Paket | Konkrete Umsetzung | Aufwand |
|---|---|---:|
| **AP6 – Referenzvalidierung** | HEOS-Ergebnisse mit geeigneten unabhängigen Referenzwerten vergleichen; Clapeyron-Gleichung, kritische Zustände und Phasengrenzen auswerten | **10 h** |
| **AP7 – vdW-Vergleich** | Kleine eigenständige vdW-Berechnung, identische einphasige Vergleichszustände, Druckabweichungen und kritisches Verhalten | **12 h** |
| **AP8 – Konvergenz und Performance** | Rasterauflösungen variieren, Rechenzeiten messen, Interpolationsfehler und Ressourcenbedarf auswerten | **6 h** |
| **AP9 – Schriftliche Studienarbeit** | Theorie, Methoden, Architektur, Ergebnisse, Diskussion, Grenzen und Fazit mit reproduzierbaren Abbildungen | **10 h** |
| **AP10 – Endprüfung** | Abschließende Testmatrix, Ergebnisse sichern, Reproduktionsanleitung und Abgabeprüfung | **2 h** |
| **Gesamt** | | **40 h** |

Die Abnahmekriterien sind:

**AP6:** Eine quantitative Vergleichstabelle enthält mindestens Druck- und Sättigungseigenschaften für ausgewählte Zustände. Abweichungen und Quellen sind angegeben. Die bereits vorhandene Clapeyron-Diagnostik wird nachvollziehbar ausgewertet.

**AP7:** Die selbst implementierte vdW-Gleichung liefert für identische homogene Zustände einen reproduzierbaren Vergleich mit HEOS. Die Grenzen des Vergleichs im Zweiphasengebiet werden ausdrücklich erläutert. Eine Maxwell-Konstruktion gehört nicht zwingend zum 12-Stunden-Paket.

**AP8:** Mindestens drei Rasterauflösungen sind hinsichtlich Rechenzeit und geometrischer beziehungsweise interpolationsbedingter Abweichung untersucht. Aus den Messungen wird eine begründete Standardauflösung abgeleitet.

**AP9:** Die schriftliche Arbeit enthält nicht nur Screenshots, sondern messbare Ergebnisse, nachvollziehbare Berechnungsmethoden und eine kritische Diskussion.

**AP10:** Ein anderer Rechner mit einer dokumentierten Python-Umgebung kann die ausgewählten Berechnungsergebnisse reproduzieren.

### Welche Arbeitspakete sind wirklich zwingend?

**Minimaler Abschluss: ca. 59 Stunden.** Technische Fehlerkorrektur, Testbasis, Spinodalprüfung, Diagrammvalidierung, unabhängige Referenzprüfung, Konvergenzuntersuchung und Dokumentation.

**Erweiterter Abschluss: 76 Stunden.** Zusätzlich p-h-T-Isolinienintegration und ein eigener, wissenschaftlich ausgewerteter van-der-Waals-Vergleich.

Die 59 Stunden sind nicht als nachträglich zwingend nachzuholender Aufwand zu verstehen. Sie sind eine Schätzung für die systematische Absicherung und Dokumentation ab dem untersuchten Stand.

**Entscheidender Hinweis zur Zeitplanung:** Wenn die zwei Wochen den gesamten Studienarbeitsaufwand einschließlich der bereits geleisteten Entwicklung bezeichnen, darf das 76-Stunden-Budget nicht einfach zusätzlich angesetzt werden. Dann sollten die bisherigen Entwicklungsstunden dokumentiert und zunächst die P0-Aufgaben mit dem tatsächlich verbliebenen Budget durchgeführt werden. Ein künstliches Ausweiten auf 76 weitere Stunden wäre fachlich nicht sinnvoll.

---

## 9. Wie sollte die wissenschaftliche Studienarbeit aufgebaut werden?

Eine geeignete Gliederung könnte folgendermaßen aussehen:

| Kapitel | Inhalt | Schwerpunkt |
|---|---|---|
| 1. Einleitung und Zielsetzung | Problemstellung, Anforderungen an interaktive 3D-Zustandsdiagramme | Motivation |
| 2. Thermodynamische Grundlagen | HEOS, Zustandsgrößen, Phasengleichgewicht, kritischer Punkt, Spinodale | Theorie |
| 3. Numerische Methodik | Zustandsberechnung, Rastererzeugung, Dampfgehaltsinterpolation, Netzgeometrie | Eigenleistung |
| 4. Softwarearchitektur | Module, GUI, PyVista, Datenfluss und Interaktion | Implementierung |
| 5. Ergebnisse und Validierung | Diagramme, Phasengrenzen, Referenzvergleich, Konvergenz, gegebenenfalls vdW | Wissenschaftlicher Kern |
| 6. Diskussion | Modellgrenzen, numerische Stabilität, Berechnungsaufwand, Fehlerquellen | Kritische Einordnung |
| 7. Fazit und Ausblick | Erfüllte Anforderungen, weitere Entwicklung und mögliche Power-Wall-Nutzung | Abschluss |

Für den wissenschaftlichen Ergebnisteil bieten sich insbesondere drei eigene Untersuchungen an.

**Erstens: Phasengrenzkontinuität.** Quantitativ zeigen, dass die von HEOS berechneten Sättigungspunkte identisch in die einzelnen geometrischen Flächen übernommen werden. Hier lässt sich eine relevante eigene Implementierungsentscheidung begründen.

**Zweitens: Thermodynamische Validierung.** Übereinstimmung von Sättigungseigenschaften, kritischen Daten und Clapeyron-Beziehung untersuchen. Numerische Konsistenz ausdrücklich von unabhängiger Stoffdatenvalidierung trennen.

**Drittens: Numerische Effizienz.** Darstellen, wie sich die Rastergröße auf Berechnungszeit und geometrische Genauigkeit auswirkt. Dadurch belegen, weshalb eine bestimmte Auflösung für den praktischen Einsatz ausreichend ist.

Der optionale vdW-Vergleich kann diese drei Untersuchungen sinnvoll ergänzen, ist aber keine Voraussetzung für eine wissenschaftlich fundierte Arbeit.

## 10. Abschließende Beantwortung der sechs Kernfragen

### 1. Wie ausgereift ist das Projekt?

Es handelt sich um einen fortgeschrittenen, modularen Softwareprototyp. Der technische Funktionsumfang ist hoch. Eine Einstufung als vollständig verifizierte Anwendung ist aufgrund fehlender End-to-End-Tests und der identifizierten Fehler nicht gerechtfertigt.

### 2. Welche wesentlichen Schwachstellen bestehen?

Die wichtigsten sind die nicht durchgängig abgesicherte p-v-T-Rasterkonsistenz, der Fehler in der Debug-Visualisierungsprüfung, die begrenzte quantitative Überprüfung der Spinodalen, fehlende unabhängige Referenzvalidierung und die unzureichend abgesicherte Regressionstestinfrastruktur.

### 3. Ist der bisherige Stand bereits ausreichend für eine zweiwöchige Studienarbeit?

Vom technischen Entwicklungsumfang her: **ja, durchaus**. Die fünf Diagrammtypen, GUI, Phasengrenzen und Datenverarbeitung repräsentieren bereits substanzielle Entwicklungsleistung.

Als wissenschaftlich abgeschlossene Master-Studienarbeit sollte der Stand aber erst nach belastbarer Validierung, Fehlerkorrektur und Dokumentation beurteilt werden. Allein aus den Projektdateien lässt sich außerdem nicht feststellen, wie viele Arbeitsstunden bisher investiert wurden.

### 4. Welche Ergänzungen bringen den größten wissenschaftlichen Mehrwert?

Die unabhängige Validierung, Konvergenzuntersuchung und quantitative Stabilitätsanalyse sind besonders wertvoll. Ein begrenzter vdW-HEOS-Vergleich wäre eine sehr geeignete zusätzliche Modelluntersuchung.

### 5. Was ist zwingend notwendig?

Die beiden reproduzierbaren p-v-T-Fehler müssen behoben werden. Die beanspruchten Funktionen müssen unter realen Bedingungen getestet werden. Thermodynamische Ergebnisse brauchen überprüfbare Vergleichswerte, und das Vorgehen muss reproduzierbar dokumentiert sein.

### 6. Wie wird die Eigenleistung nachweisbar?

Durch explizite Beschreibung und Prüfung der eigenen Implementierungsentscheidungen: Wahl und Konstruktion der Raster, geometrische Verbindung von Phasengrenzen und Teilflächen, Modellierung der Zweiphasenfläche, Validierung der Spinodalen, Behandlung numerischer Grenzfälle und quantitative Analyse der Berechnungseffizienz.

Diese Entscheidungen sind eigenständige Ingenieurarbeit, auch wenn CoolProp die zugrunde liegende Stoffdatenmodellierung übernimmt.

---

## Gesamturteil

**Das Projekt sollte nicht mehr wesentlich in die Breite entwickelt werden.** Fünf vollständig auswählbare Diagrammtypen sind für einen Bearbeitungszeitraum von zwei Wochen bereits ein sehr umfangreicher Ausgangspunkt.

Für einen überzeugenden Abschluss bietet sich folgende Schwerpunktsetzung an:

1. **Technisch:** Bestehende Implementierungen stabilisieren und die nachgewiesenen Fehler korrigieren.
2. **Wissenschaftlich:** Phasengleichgewicht, kritische Zustände und numerische Genauigkeit quantitativ untersuchen.
3. **Optional vertiefend:** Eine kleine eigene van-der-Waals-Implementierung mit HEOS vergleichen.
4. **Dokumentarisch:** Den Zusammenhang zwischen thermodynamischer Theorie, numerischer Berechnung und erzeugter 3D-Geometrie nachvollziehbar darstellen.

Damit wird aus einer bereits umfangreichen Visualisierungssoftware eine wissenschaftlich begründete Studienarbeit, ohne dass dafür eine neue Softwarearchitektur, weitere Diagrammtypen oder aufwendige GUI-Erweiterungen erforderlich sind.

**Stärkste Empfehlung:** **Den vorhandenen p-v-T-Code als zentralen wissenschaftlichen Untersuchungsgegenstand verwenden und die anderen vier Diagrammtypen als Nachweis der Erweiterbarkeit und Übertragbarkeit der Architektur behandeln.** So erhält die Arbeit einen klaren roten Faden und bleibt trotz des bereits erheblichen Funktionsumfangs in einem realistischen zeitlichen Rahmen.

---

## Verwendung dieser Datei in zukünftigen Chats

Folgende Bezugnahmen sind eindeutig:

- **„Beginne mit AP1“:** Analysiere den aktuellen Zustand von `Test.py`, `compare_pvt.py`, `snapshot_pvt.py` und `Debug.py`; anschließend schrittweise die Test- und Regressionstest-Infrastruktur verbessern. Prüfe explizit, ob fehlende Baselines fälschlich als Erfolg gelten.
- **„AP2 prüfen/umsetzen“:** Gemeinsames p-v-T-Raster, Sättigungsgrenzen, Binodale/Flächenränder und den reproduzierten Fehler mit fehlenden Sättigungszeilen untersuchen und korrigieren.
- **„AP3“:** Physikalisch/numerische Spinodalvalidierung einschließlich Ableitungsresiduum und kritischer Grenzfälle.
- **„AP4“:** Systematische Tests für alle fünf Diagrammtypen mit repräsentativen Fluiden.
- **„AP5“:** Vorhandene p-h-T-Isolinien in den Berechnungs-/GUI-/Renderpfad integrieren; nicht die Berechnungsfunktionen unnötig neu schreiben.
- **„AP6“:** Thermodynamische Validierung gegen unabhängige Referenzen und thermodynamische Identitäten.
- **„AP7“:** Begrenzter, wissenschaftlicher Vergleich eigener van-der-Waals-EOS mit CoolProp HEOS.
- **„AP8“:** Numerische Rasterkonvergenz und Performancebenchmark.
- **„AP9“:** Schriftliche Studienarbeit auf Grundlage des Codes, der Tests und der Validierungsergebnisse verfassen.
- **„AP10“:** Gesamtprüfung und reproduzierbares Abschluss-/Abgabepaket.

**Maßgebend sind die ausführlichen Definitionen und Abnahmekriterien in Abschnitt 8, nicht nur die Kurzformeln dieses Indexes.**
