# ThermoSurface – AP2 integriert

Dieses Archiv enthält den bereitgestellten Quellcode mit gezielten Änderungen an
`coolprop_heos_pvt.py`, `main.py`, `visualization.py`, `gui.py` und `Test.py`.
Die übrigen Python-Module sowie `config.yaml` wurden aus der bereitgestellten
Projektversion unverändert übernommen. AP1 wird nicht zurückgesetzt.

## Installation und Start (Windows, im entpackten Projektordner)

```powershell
python -m pip install -r requirements.txt
python main.py
```

Die GUI setzt einen grafikfähigen Rechner mit funktionierender Qt-/VTK-Umgebung
voraus. Auf einem Rechner mit mehreren Python-Installationen muss der richtige
Interpreter verwendet werden.

## AP2-Änderungen

- Einheitliches Temperatur- und Sättigungsraster für p-v-T; kritischer Knoten
  ist in der Oberfläche enthalten. Der dargestellte Bereich muss die kritische
  Temperatur von unten bis oberhalb des kritischen Punkts umfassen.
- Nur eine Binodalberechnung; die Sättigungsdaten werden für alle p-v-T-Teilflächen
  konsistent wiederverwendet.
- Exakt gleiche Randknoten bei `Q=0` und `Q=1` wie auf Siede-/Taulinie.
- AP2-Datenvalidator gegen fehlende Sättigungszeilen, unpassende Arrayformen,
  verschobene Phasengrenzen und zu enge Volumengrenzen.
- Das schon vorhandene randkonforme Einphasenmesh ist jetzt die Standardfläche.
  Die Zweiphasenfläche wird unabhängig davon ein-/ausgeblendet.
- Zusätzliche AP2-Tests in `Test.py`.

## Tests

```powershell
# AP1: rein numerische Regressionstest-Logik, ohne CoolProp/GUI
python -m unittest Test.RegressionTests -v

# AP2: erfordern CoolProp / ggf. PyVista
python -m unittest Test.ThermoSurfaceTests.test_ap2_shared_temperature_grid -v
python -m unittest Test.ThermoSurfaceTests.test_ap2_exact_two_phase_edges -v
python -m unittest Test.ThermoSurfaceTests.test_ap2_heos_mesh_boundaries -v
```

Bei der GUI sind besonders die getrennte Sichtbarkeit der Zustandsoberfläche
und des Zweiphasengebiets sowie der Gitterschalter zu prüfen. Die anderen
Diagrammtypen (`pht`, `pst`, `hsp`, `uvt`) sind unverändert geblieben und sollten
auf dem Zielrechner stichprobenartig aufrufbar sein.

## Wichtige Hinweise zur Regression

Das ursprüngliche Verzeichnis `regression/` mit gespeicherten AP1-Snapshots war
in den bereitgestellten Quelldateien **nicht enthalten** und wird nicht künstlich
erzeugt. Deshalb schlägt `Test.RecordedRegressionTests` ohne mitgelieferte
Snapshots absichtlich fehl. Bitte die eigenen AP1-Baselines bei Bedarf separat
aus dem ursprünglichen Projekt übernehmen und **nicht überschreiben**.

AP2 verändert gezielt die Sättigungszeilen und daher möglicherweise Arrayformen.
Ein direkter numerischer Vergleich zu vor-AP2-Snapshots kann aus diesem Grund
fehlschlagen. Neue Baselines nur nach fachlicher Kontrolle der AP2-Daten anlegen.

## Grenzen der hier durchgeführten Validierung

- Python-Syntax aller Module, AP1-Regressionstest-Logik und isolierte
  synthetische AP2-Raster-, Grenz-, Integrations-, Mesh- und GUI-Mocktests wurden
  erfolgreich geprüft.
- CoolProp, PyVista, PyVistaQt und PySide6 standen in der Testumgebung nicht
  zur Verfügung. Ein realer HEOS-Lauf und die echte Qt-/VTK-Oberfläche konnten
  deshalb **nicht** ausgeführt werden; diese Tests sind lokal erforderlich.
- Die vorhandenen Spinodalen wurden nicht im Sinne des späteren AP3 neu bewertet.
- Die p-v-T-Teilflächen benötigen Volumengrenzen, die die Sättigungsvolumina
  sämtlicher unterkritischer Temperaturzeilen vollständig einschließen.
  Ungültige Kombinationen werden ausdrücklich abgelehnt.
