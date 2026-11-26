# Viewport-Wiederherstellung und diagrammübergreifende Zustandsmarker

## Viewport-Wiederherstellung

Berechnete HEOS-Daten und PyVista-Szenen werden getrennt behandelt. Ein Ergebnis-Cache-Treffer überspringt das Rendern nur, wenn der Viewport denselben Schlüssel besitzt und die Szene als `READY` geprüft wurde. Die Prüfung kontrolliert den Renderer, die vorhandenen Oberflächen-Actors und ihre ursprüngliche Geometrie mit Zellen und endlichen Bounds. Ein absichtlich leerer Ansichtsbereich bleibt gültig; ein schwarzes Bild durch einen GPU-/Treiberfehler kann die Prüfung nicht sicher erkennen.

Die Gruppe **Ansichtsbereich** enthält **Darstellung neu aufbauen**. Die Aktion ersetzt nur den aktiven Viewport und rendert dessen vorhandene Daten erneut. Layer, manuelle Grenzen, Gittereinstellung und Kamera werden übernommen. Wenn ein Renderer fehlschlägt, wird der Fehler protokolliert und in der Statusleiste angezeigt; es wird nicht endlos erneut versucht.

Für Diagnosezwecke protokolliert das Modul `thermosurface.viewport` unter anderem `render_begin`, `render_end`, `render_error` und `viewport_dispose`; `thermosurface.gui` protokolliert Cache-Übersprünge und manuelle Wiederherstellungen. Vergleichsfenster verwenden einen unabhängigen Plotter. Schließen blendet das Fenster zur Wiederverwendung aus; beim Beenden der Hauptanwendung werden Plotter und Selection-Listener freigegeben.

## Hover und Isolinien

Der Hover zeigt den nächsten gespeicherten **Gitterpunkt** des getroffenen Oberflächen-Meshes mit den drei nativen Diagrammgrößen und Einheiten. Durch manuelles Clipping erzeugte Schnittpunkte werden als interpoliert gekennzeichnet. Mausbewegungen lösen keinen CoolProp-Flash aus. Für strukturierte Oberflächen wird die native Punktmetadaten-Reihenfolge an die VTK-Punktordnung angepasst.

Die Isolinien heißen nach der konstanten physikalischen Größe: Isothermen `T = konstant`, Isobaren `p = konstant` und Isochoren `v = konstant`. Checkboxen bleiben deaktiviert, wenn für das aktive Diagramm kein entsprechender Actor vorhanden ist.

## Diagrammübergreifende Zustandsprojektion

In **Zustandsvergleich** das gewünschte Zieldiagramm auswählen und **Ausgewählte Punkte auf 3D-Fläche anzeigen** aktivieren. Vorhandene Zustände werden über ihre gespeicherten molaren SI-Werte und den normalen Ziel-Adapter in native Zielgrößen und anschließend in die Darstellungskoordinaten des tatsächlichen Ziel-Renderers überführt. Massenspezifische Achsen verwenden die im Zustand gespeicherte molare Masse. Es wird weder ein neuer Zustand geflasht noch an die Zieloberfläche „angeschnappt“.

Die Zieloberfläche kommt aus dem normalen Berechnungs-/Cacheweg. Jede darstellbare Auswahl wird mit ihrer ursprünglichen Farbe und stabilen Nummer als eigener, nicht auswählbarer Marker gezeichnet; zwischen Markern entstehen keine Linien. Nicht darstellbare Zustände werden in der Statusanzeige gezählt und ihre Gründe im Tooltip aufgeführt. Ausgeschlossen werden unter anderem mehrdeutige g-p-T-Sättigungszustände, homogene Fortsetzungen ohne gesicherte Zuordnung, Zustände außerhalb des sichtbaren Oberflächenbereichs und Zweiphasenzustände bei ausgeblendetem Zweiphasenlayer. Die Prüfung des Koordinatenbereichs weist nicht nach, dass die Oberfläche jeden Punkt geometrisch abdeckt. Marker außerhalb eines manuell begrenzten Ansichtsbereichs werden ausgeblendet und entsprechend gemeldet; nach dem Zurücksetzen erscheinen sie wieder. Für die Drucktransformation gilt die gewählte Zielachse: logarithmisch nur für positiven Druck, `asinh`, `symlog` und linear auch für nichtpositive Werte. Die bisherigen Ursprungsmarker und das Vergleichsdiagramm bleiben erhalten.

## Point Picker

Der erste VTK-Pick verwendet einen exakten Sichtstrahl. Die vorherige VTK-Toleranz konnte benachbarte Isolinien oder falsche Dreiecke anstelle der tatsächlich angeklickten Fläche treffen. Wenn dieser Strahl ein Dreieck trifft, die Position aber außerhalb dieses Dreiecks liegt, werden versetzte Sichtstrahlen bis maximal 7 physische Pixel Abstand versucht. Dabei werden nur Treffer auf demselben Actor und derselben Zelle angenommen. Bei einem vollständig fehlenden Treffer ist ein Fallback nur über leerem Hintergrund **außerhalb der konvexen Bildschirmhülle aller sichtbaren Zustandsflächen** zulässig. Der Kandidat muss sichtbar sein und zugeordnet werden können; sichtbare Vordergrundflächen werden auch dann berücksichtigt, wenn sie nicht pickbar sind. Andere Actor-Treffer, ungültige Zellen und Fehler werden nicht umgedeutet. Ein Fallback wird in der Statusleiste als solcher markiert und im Debug-Log protokolliert. Die baryzentrische Toleranz für Rundungsfehler beträgt `1e-6`; echte Extrapolation außerhalb des Dreiecks ist ausgeschlossen.

Die Hüllenprüfung ist bewusst konservativ: Randklicks innerhalb der konvexen Hülle, etwa an einem konkaven Außenrand, werden nicht angeschnappt. So kann der Picker nicht allein aufgrund einer nahegelegenen Fläche über Netzlöcher oder Lücken zwischen Zweigen hinweg wählen.

### Kontrollierte Picker-Messung

Auf Windows wurden 50 mit `QTest` ausgelöste Qt-Mausklicks auf unabhängig mit einem VTK-Sichtstrahl bestätigte sichtbare Zellen gemessen: je drei Blickrichtungen in p-v-T, p-h-T und f-v-T, insgesamt 18 nahe der Dreiecksmitte, 18 nahe einer Kante und 14 nahe einer Ecke. Ein Treffer zählt nur, wenn der gewählte Actor und die Zell-ID der Referenz entsprechen und alle Quellknoten zugeordnet werden konnten. Die gleiche Stichprobe wurde mit dem Picker aus Git-Revision `778a70b8` und mit der Änderung ausgewertet; die Messung ersetzt keine Prüfung mit physischer Maus:

| Diagramm | Vorher | Nachher |
|---|---:|---:|
| p-v-T | 6/18 | 18/18 |
| p-h-T | 3/15 | 15/15 |
| f-v-T | 11/17 | 17/17 |
| **Gesamt** | **20/50 (40 %)** | **50/50 (100 %)** |

Ein geometrischer Test prüft außerdem einen äußeren Randtreffer sowie die Ablehnung eines Netzlochs, fernen Hintergrunds, einer ausgeblendeten Fläche und verdeckender Vordergrund-Actors. Dies ist keine vollständige Abnahme aller Netzlücken, Phasenränder, Schichten oder DPI-Konfigurationen. Bei 50 aufeinanderfolgenden Auswahländerungen blieb die Zahl der Renderer-Actors und die Kamera des Projektionsziels unverändert.

## Reproduzierbare manuelle GPU-Abnahme

Die automatischen Tests decken Cache-Recovery, expliziten Szenen-Neuaufbau, Marker-Overlay, Hovermetadaten und den begrenzten Picker-Fallback ab. Eine manuelle Prüfung mit physischer Windows-GPU, mehreren DPI-Skalierungen und realem Maus-/Fensterverhalten wurde in dieser Implementierung **nicht ausgeführt**. Für eine Freigabe zusätzlich:

1. Einzel- und Mehrfachansicht mit allen Diagrammen berechnen, rotieren und zoomen.
2. Vergleichsfenster mindestens 20-mal öffnen und schließen; Hauptviewport danach bedienen.
3. Identische Berechnungen aus dem Cache ausführen und nach Bedarf **Darstellung neu aufbauen**.
4. Auswahlpunkte aus mehreren Diagrammen auf `p-v-T` projizieren, löschen und umsortieren; auf fehlende Linien und unveränderte Oberfläche/Kamera achten.
5. Point Picker bei 100 %, 125 % und 150 % Windows-Skalierung über Flächenränder, Lücken, verdeckte Bereiche und Nicht-Oberflächen-Actors testen.
6. Sidebar-Beschriftungen und Bedienbarkeit bei üblichen Fenstergrößen prüfen.
7. Die Picker-Stichprobe um die übrigen Diagrammtypen, reale Flächenränder, Datenlücken und Phasenbereiche erweitern und die Negativ-Testreihe auf der Zielhardware wiederholen.

