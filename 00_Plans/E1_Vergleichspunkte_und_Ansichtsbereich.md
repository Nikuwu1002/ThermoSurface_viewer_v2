# E1 – Ausgewählte Vergleichspunkte löschen und automatischen Ansichtsbereich wiederherstellen

**Projekt und Grundlage:** ThermoSurface. **Module:** `gui.py`, `selection_manager.py`, `comparison_view.py`, `diagram_viewport.py`, `viewport_bounds.py`, `axes_view.py`, `diagram_axes.py`. Beide Korrekturen sind unabhängig voneinander umsetzbar. Die Abnahmetests sind vorgesehen und noch nicht durchgeführt.

## 1. Ziel und Umfang

Zwei getrennte Fehler beheben. Erstens soll „Löschen“ genau den gewählten Vergleichspunkt einschließlich seiner Quellenmarker, Projektionen und Darstellung im Vergleichsfenster entfernen. Zweitens soll „Zurücksetzen“ im Block „Ansichtsbereich“ die vollständige Geometrie wiederherstellen und die sechs dazugehörigen automatischen Zahlenwerte anzeigen.

## 2. Ist-Zustand im Code

**Vergleichspunkte.** Das auswählbare Steuerelement in `ThermoSurfaceWindow._build_compare_group()` ist ein `QListWidget` namens `compare_list`. Das separate Fenster `ComparisonView` besitzt eine `QLabel`-Legende. Ein Klick auf einen Legendeneintrag wählt keine Zeile der Hauptliste aus.

`_on_point_clicked()` führt zur Aufnahme eines `SelectedPoint` in `SelectionManager`. `selection_index` ist seine dauerhafte Nummer, `selection_order` seine veränderliche Verbindungsreihenfolge. `_refresh_compare_list()` speichert `selection_index` in `Qt.ItemDataRole.UserRole`. `on_compare_delete()` liest diesen Wert über `_selected_compare_index()` und ruft `self.selection.remove(index)` auf. `SelectionManager.remove()` sucht mit `get()` nach genau dieser Nummer und meldet `"removed"`. `reset()` meldet `"reset"`. Beide Ereignisse erreichen `_on_selection_event()`, das Liste und Marker aktualisiert. Das geöffnete Vergleichsfenster aktualisiert sich über einen eigenen Listener.

**Konkrete Schwachstelle:** Fehlt `compare_list.currentItem()`, endet `on_compare_delete()` ohne Rückmeldung. Beim Neuaufbau stellt `_refresh_compare_list()` nur eine zuvor vorhandene, weiterhin gültige ID wieder her. Nach dem Hinzufügen ist der neue Punkt daher nicht ausdrücklich zur Löschung ausgewählt. Für den Fall „B wurde in der Hauptliste angeklickt und bleibt trotzdem erhalten“ zeigt der Code keinen eindeutigen Zuordnungsfehler.

**Ansichtsbereich.** `on_view_reset()` ruft `DiagramViewport.reset_manual_view()` auf. `ViewportState.reset()` setzt lediglich `enabled = False`. `fields`, `bounds` und `user_defined` bleiben erhalten. `_apply_view()` stellt über `clip_registry.apply(None, ...)` die ursprünglichen Datensätze wieder her und passt die Kamera an. Anschließend verwendet `_refresh_view_panel()` aber `state.fields or self._suggest_view_fields(viewport)`. Damit gewinnen die alten manuellen Texte gegen die bereits wiederhergestellte automatische Ansicht.

## 3. Umsetzungsschritte

1. **Löschpfad reproduzieren.** A, B und C anlegen. B in `compare_list` auswählen. Direkt beim Eintritt in `on_compare_delete()` die aktuelle Zeile, die IDs aus `selectedItems()`, den `UserRole`-Wert und die IDs in `selection.points` vergleichen. Danach `"removed"`, die verbleibenden IDs und den Aufruf von `_sync_markers()` verfolgen. Damit lässt sich zwischen fehlender Auswahl, fehlgeschlagener Modelländerung und veralteter Darstellung unterscheiden.
2. **Auswahl eindeutig machen.** Den Löschbutton an eine tatsächlich ausgewählte Listenzeile koppeln und bei leerer Auswahl deaktivieren. Die stabile ID aus dem gewählten Item verwenden. Beim Listenaufbau eine weiterhin vorhandene, ausdrücklich ausgewählte ID wiederherstellen. Nach dem Löschen kann der nächste vorhandene Listeneintrag ausgewählt werden. Keine Neunummerierung vornehmen.
3. **Marker über den vorhandenen Ereignispfad aktualisieren.** `_sync_markers()` erzeugt Quellenmarker über `points_of(viewport.diagram_id)` und Projektionen aus den verbliebenen Punkten. `set_selection_markers()` und `_draw_markers()` entfernen alte Actors und bauen die gültige Menge auf. Falls das Modell B korrekt entfernt, aber ein Marker bleibt, anhand seiner ID und seines Viewports prüfen, ob Quellenmarker, Projektion oder Text-Actor zurückbleibt. Nur die dabei identifizierte Stelle ändern.
4. **Expliziten Reset vom vorübergehenden Ausschalten trennen.** Für `on_view_reset()` den manuellen Zustand vollständig verwerfen, etwa durch einen neuen `ViewportState` für dieselbe Diagramm-ID. Das einfache Abwählen der Checkbox darf weiterhin die letzten Eingaben für eine spätere Wiederverwendung behalten. Beide Bedienhandlungen benutzen derzeit `reset_manual_view()` und benötigen deshalb eine unterscheidbare Semantik.
5. **Automatische Grenzen einmal aus der wiederhergestellten Szene ableiten.** Nach `ClipRegistry.apply(None, ...)` denselben Rückgabewert von `axes_view.auto_bounds()` für Kamera, Achsenrahmen und Feldanzeige verwenden. Bei p-v-T berücksichtigt die Funktion sichtbare Actors, bei Zusatzdiagrammen auch ausgeblendete Geometrie. Für die sechs Felder die Box mit `diagram_axes.display_bounds_to_physical()` in die jeweiligen Anzeigeeinheiten umrechnen. `_view_edit` und `_view_prefs` für genau dieses Diagramm aktualisieren. Fehlt eine gültige Box, den Grund anzeigen.

## 4. Abnahmetests

| Ausgangszustand | Bedienhandlung | Soll-Ergebnis |
|---|---|---|
| A, B, C mit verschiedenen Koordinaten, Vergleichsfenster und Zielprojektion geöffnet | B in der Hauptliste anklicken und löschen | Nur B verschwindet aus Liste, Quellenansicht, Zielprojektion, Legende und Vergleichsgeometrie. A und C behalten Nummer und Zustand |
| A und C verbleiben | C auswählen und löschen | Genau A bleibt erhalten. Der geänderte Zeilenindex von C beeinflusst die Zuordnung nicht |
| Liste enthält Punkte, keine Zeile ist ausgewählt | Zustand des Löschbuttons prüfen, anschließend „Alle zurücksetzen“ betätigen | Gezieltes Löschen ist deaktiviert. „Alle zurücksetzen“ leert Modell und betroffene Ansichten |
| Aktives Diagramm mit sechs angewendeten manuellen Grenzen | „Zurücksetzen“ klicken | Alle sechs Felder entsprechen innerhalb der Anzeigepräzision der rücktransformierten automatischen Box. Originalgeometrie, Achsenrahmen und Kameraausschnitt passen zusammen |
| p-v-T und ein Zusatzdiagramm mit unterschiedlichen Grenzen geöffnet | Nur p-v-T zurücksetzen, hin- und zurückwechseln | Jedes Diagramm zeigt seine eigenen Werte. Der manuelle Ausschnitt des Zusatzdiagramms bleibt erhalten |

## 5. Abhängigkeiten und offene Entscheidungen

Für die Löschkorrektur ist der Laufzeitbefund aus Schritt 1 entscheidend. Eine Umstellung auf Zeilennummern würde die bereits richtige Trennung von Identität und Reihenfolge beschädigen. Beim Ansichtsreset ist die Ursache hinreichend klar, um die Korrektur direkt umzusetzen.
