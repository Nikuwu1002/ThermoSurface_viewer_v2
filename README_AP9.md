# AP9 – Schriftliche wissenschaftliche Studienarbeit

**Ziel:** Dokumentation der Thermodynamik, numerischen Methodik, Architektur, Ergebnisse, Bewertung und Grenzen in einem eigenständigen wissenschaftlichen Text.

## Erstelltes Artefakt

- `STUDIENARBEIT_ThermoSurface.md`: ausgearbeiteter Bericht mit Zusammenfassung, Einleitung, Grundlagen (HEOS, vdW, Zweiphasenzustände, kritischer Punkt, Clapeyron), numerischen Methoden (gemeinsame Sättigungsrandknoten, Interpolation, Validierung), Modularchitektur, Ergebnissen, Diskussion, Fazit und Literaturverzeichnis.
- `ap9_figure.py`: reproduzierbare **analytische** vdW-Isothermengrafik auf Grundlage von `p_r=8T_r/(3v_r−1)−3/v_r²`. Diese Gleichung ist äquivalent zur verwendeten vdW-EOS in reduzierten Koordinaten. Das Bild `figures/vdw_reduced_isotherms.png` wurde tatsächlich generiert und beigefügt. Darstellung der homogenen EOS und nicht der Phasengleichgewichtsschleife.
- Quellen: die vorhandene Masterarbeit von Tirapong Nualsri (2020) als fachlicher Hintergrund sowie das offizielle IAPWS-IF97-Release, CoolProp-Publikation und Projektquellen. Es wurden keine erfundenen Literaturangaben in die numerische Referenzvalidierung eingebracht.

## Tatsächliche Verifikation und offene Anforderungen

Die Abbildung wurde mit der installierten Python-/NumPy-/Matplotlib-Umgebung erzeugt. Die analytischen kritischen vdW-Eigenschaften sind durch AP7-Unit-Tests abgesichert. In der Studienarbeit werden die 23 tatsächlich ausgeführten softwaretechnischen Tests und sämtliche **ausstehenden** HEOS-, IAPWS-, Performance- sowie realen GUI-Prüfungen ausdrücklich getrennt. Keine HEOS-Messwerte wurden erfunden oder durch synthetische Zahlen ersetzt.

**Status:** Vollständige Textstruktur und fachlich ausgearbeiteter Entwurf mit nachgewiesener mathematischer Abbildung vorliegend; **Abgabeversion noch nicht vollständig**, weil die messbaren Ergebnis- und Diskussionsabschnitte aus AP6 bis AP8 erst nach realem Lauf der Skripte numerisch vervollständigt werden können. Zusätzlich sind formale hochschulspezifische Zitier-/Layoutvorgaben nicht bekannt und nicht implizit erfüllt.
