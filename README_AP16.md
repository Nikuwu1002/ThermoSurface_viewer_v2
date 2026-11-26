# AP16 – PyVista-Actors und GUI-Ebenen

## Änderungen

`visualization.py` erhielt `homogeneous_surface_mesh()` und `add_homogeneous_surface()`, die drei homogenen Actor-Keys werden in `render_pvt()` als eigenständige Meshes erzeugt. Jeder Netzstreifen verwendet zwei Dreiecke pro vollständiger benachbarter Viererknoten-Zelle. Nur gültige physikalische Koordinaten mit `p>0` werden gerendert: Über NaN-Zeilen und negative Druckbereiche werden **keine** Flächendreiecke gezogen. Der eigentliche thermodynamische Datensatz bleibt ungekürzt.

Neue GUI-Checkboxen in `gui.py`:

- `metastable_liquid_surface` – metastabile Flüssigkeit
- `metastable_vapor_surface` – metastabiler Dampf
- `unstable_surface` – instabile homogene Fortsetzung (anfangs aus)

Die bisherige `layer_checks`-, `layer_preferences`-, `_apply_render_layers`- und `set_layer_visibility`-Logik wird genutzt. Das Oberflächengitter wirkt auch auf diese drei Actors. Für `pht`, `pst`, `hsp` und `uvt` existieren die neuen Actors nicht und ihre Checkboxen sind deaktiviert. Die alten Siede-, Tau-, Spinodal- und Gleichgewichts-Actors bleiben bestehen. Für eine bei AP13 exakt berechnete Spinodalkurve priorisiert der Renderer die gemeinsamen Randdaten; andernfalls wird die vorhandene Kurve verwendet.

## Doppelprüfung

**A:** Struktur- und Renderpfadprüfung anhand der vorhandenen PyVista-API. **B:** Ein PyVista-Double überprüft ausgeführte Triangulationslogik mit vollständigen Zeilen, NaN-Unterbrechungen und negativen Drücken. Der analytische Mesh-Test besteht. Die echte PyVista-Bibliothek und Qt sind hier nicht installiert; deshalb wurden weder Kameratests noch GUI-Klicktests oder Screenshots angefertigt. Diese Akzeptanzkriterien bleiben **offen** und müssen am Zielrechner nachgeholt werden.
