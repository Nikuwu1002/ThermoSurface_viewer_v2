# ThermoSurface – Gesamtergebnis AP5 bis AP10 (19.10.2026)

## Ausgangslage und Architektur

Als Ausgangsbasis wurden die vom Nutzer gelieferten Originaldateien nach AP4 übernommen. Unter anderem existieren: fünf HEOS-Berechnungsmodule, `main.py` als Dispatcher, `gui.py`, `visualization.py`, `Test.py` (AP1–AP3) und `test_ap4.py`. Die Funktion von AP4 ist die Prüfung aller fünf Diagramme; die ursprüngliche AP4-README bezeichnet die reale HEOS-/Qt-/VTK-Abnahme bereits als offen. Die hier vorliegende Laufzeitumgebung kann CoolProp, PyVista, PyVistaQt und PySide6 nicht installieren. Der Quellcode wurde in einem neuen Verzeichnis bearbeitet; die vom Nutzer hochgeladenen Dateien wurden nicht überschrieben.

Der stabile Berechnungs-/GUI-/Renderer-Aufbau blieb bestehen: GUI → Dispatcher in `main.py` → gewähltes HEOS-Modul → Ergebnis-Dictionary → PyVista/Qt-Renderer. Separate Auswerteskripte für Wissenschaft und Reproduzierbarkeit wurden hinzugefügt, damit vdW nicht als zusätzliche Modelloption in den GUI-Produktionspfad eindringt.

## Implementierungs- und Prüfstatus

| Paket | Tatsächlich umgesetzt | Tatsächlich geprüft | Noch ausstehend |
|---|---|---|---|
| AP5 | p-h-T-Isothermen, Isobaren und Isochoren im Ergebnis-Dictionary, PyVista-Actor-Gruppen und GUI; HSP vermeidet p-h-T-Doppelarbeit | 4 AST-/Mock-Integrationstests bestanden | Realer HEOS-Lauf und interaktiver Qt-/VTK-Schaltertest |
| AP6 | Externe Wasserreferenzen aus IAPWS-IF97, quantitative Fehlertabellen- und Clapeyron-Prüfroutine | 2 Testfälle für Daten-/Fehlerlogik bestanden, Quellen manuell kontrolliert | HEOS-Fehlertabelle noch nicht gerechnet |
| AP7 | Eigenständige vdW-EOS, kritische Identität, Einphasenvergleichs-Skripte, CSV/JSON/PNG-Export | 2 analytische Tests bestanden | Tatsächliche HEOS-vdW-Druckabweichungen noch nicht gerechnet |
| AP8 | Drei Auflösungen, Zeit-, Speicher- und Interpolationsfehlerverfahren, Datenexport/Standardraster-Entscheidung | 2 algorithmische Tests bestanden | Reale Messwerte und empirische Rasterempfehlung fehlen |
| AP9 | Wissenschaftlich ausgearbeiteter Markdown-Text mit Abbildung, Theorie, Methoden, Architektur, Diskussion und Literatur | Abbildung generiert, fachlich/technisch mit vorliegenden Befunden abgeglichen | Numerische Ergebnisse von AP6/AP7/AP8, formale Abgabeprüfung |
| AP10 | Auditprogramm mit Versionen, Hashes, Tests und Reproduktionsworkflow | Statische Audits und 23 Unit-/Mocktests bestanden | Vollständige Endabnahme, Original-Baselines und Zweitrechnernachweis |

Die Summe **23** ergibt sich aus 6 AP1-Tests + 7 AP4-Synthetiktests + 10 neuen AP5–AP8-Tests. Alle positiven Ergebnisse sind reine Software-/Analytiktests; kein einziges HEOS-/GUI-Ergebnis wird als gemessen behauptet.

## Dateien

**Verändert:** `coolprop_heos_pht.py`, `coolprop_heos_hsp.py`, `visualization.py`, `gui.py`.

**Unverändert aus der Eingabe übernommen:** `main.py`, `coolprop_heos_common.py`, `coolprop_heos_pvt.py`, `coolprop_heos_pst.py`, `coolprop_heos_uvt.py`, `config.yaml`, `limits.py`, `Debug.py`, `Test.py`, `test_ap4.py`, `ap4_validation.py`, `compare_pvt.py`, `snapshot_pvt.py` und vorhandene READMEs.

**Neu:** `ap6_reference.py`, `ap7_vdw_compare.py`, `ap8_convergence.py`, `ap9_figure.py`, `ap10_verify.py`, `test_ap5_ap8.py`, `requirements.txt`, `README_AP5.md` bis `README_AP10.md`, die vorliegende Gesamtübersicht und `STUDIENARBEIT_ThermoSurface.md` mit analytischer Abbildung sowie `results/ap10/ap10_audit.json` aus dem hier ausgeführten statischen Audit.

## Ausführung und vollständige Abnahme

Vollständige Windows-Schrittfolge siehe `README_AP10.md`, Beschreibung wissenschaftlicher Auswertungen siehe `README_AP6.md` bis `README_AP9.md`.

```powershell
python -m pip install -r requirements.txt
python -m unittest Test.RegressionTests test_ap4.AP4SyntheticTests test_ap5_ap8 -v
python ap6_reference.py
python ap7_vdw_compare.py
python ap8_convergence.py
python ap10_verify.py
python main.py
```

Bei fehlenden Bibliotheken schlagen die wissenschaftlichen HEOS-Skripte **offen** fehl und geben keine erfundenen Werte aus. `ap10_verify.py` ohne `--static-only` darf bei SKIP-Tests/fehlenden Baselines keinen Erfolg melden. Im ZIP fehlen die originalen AP1-`regression/`-Referenzdateien, da sie vom Nutzer nicht bereitgestellt wurden. Ein Baseline-Neuaufbau darf erst nach fachlicher Prüfung erfolgen.

## Wissenschaftliche Schlussfolgerung

Die technische Integration der p-h-T-Isolinien lässt sich synthetisch nachvollziehen; ein Echtbetrieb ist noch zu belegen. Der vdW-Modellvergleich besitzt analytische Aussagekraft, insbesondere `Zc,vdW=0,375`, liefert aber ohne HEOS-Läufe keine quantitativen Druckabweichungsdiagramme. Für eine wissenschaftliche Studienarbeit sind die **externen Vergleichszahlen**, die **Konvergenzmessungen** und eine **echte GUI-Abnahme** entscheidende offene Nachweise.

**Gesamturteil:** Konsistenter, syntaktisch und in wesentlichen numerischen Hilfsfunktionen getesteter **Zwischenstand** mit umfangreich implementierten AP5–AP10-Werkzeugen. **Keine vollständige Abnahme aller sechs Arbeitspakete** in der verfügbaren Umgebung.
