# ThermoSurface – Installations- und Prüfhinweise

Die aktuelle Ordnerstruktur sowie ausführbare Start- und Testbefehle stehen in [README_Projektstruktur.md](../README_Projektstruktur.md). Der folgende Inhalt dokumentiert einen früheren AP5–AP10-Arbeitsstand.

Dieses Archiv enthält den vom Nutzer gelieferten aktuellen Python-Projektstand nach AP4, gezielte Erweiterungen für AP5 und separate wissenschaftliche Auswertungs-/Testwerkzeuge für AP6–AP10. **Ein vollständig abgenommener Stand der Arbeitspakete AP5–AP10 liegt noch nicht vor.** Insbesondere fehlt die Laufzeitprüfung mit CoolProp und die Grafik-/GUI-Abnahme; weitere Einzelheiten stehen in der Gesamtübersicht und den sechs AP-READMEs.

## Start unter Windows (PowerShell)

```powershell
cd .\ThermoSurface_AP5-AP10
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m unittest Test.RegressionTests test_ap4.AP4SyntheticTests test_ap5_ap8 -v
python main.py
```

Der Programmstart kann ein grafikfähiges System mit Qt-/OpenGL-/VTK-Treibern erfordern. Python 3.11 ist eine Empfehlung, kein in dieser Umgebung getestetes Garantie-Setup.

## Numerische Untersuchungen und Endprüfung

```powershell
python ap6_reference.py
python ap7_vdw_compare.py
python ap8_convergence.py
python ap9_figure.py
$env:THERMOSURFACE_AP4_GUI = "1"
python ap10_verify.py
```

Die Reproduktionsanleitung in `README_AP10.md` enthält die nötigen Detailtests einschließlich der ursprünglichen AP1-Referenzsnapshots, die dem verfügbaren Upload fehlten. Die aktuell mitgelieferte Datei `results/ap10/ap10_audit.json` zeigt die hier ausgeführten **23 bestandenen Unit-/Mocktests und neun übersprungene Realtests**. Ohne Baselines und ohne erfolgreich ausgeführte Realtests wird `accepted=false` ausgegeben.

**Wichtige Dokumente:** `README_AP5-AP10_Gesamtuebersicht.md` (Fortschritt), `TESTPROTOKOLL_AP5-AP10.md` (tatsächliche Prüfergebnisse), `STUDIENARBEIT_ThermoSurface.md` (wissenschaftliche Ausarbeitung).
