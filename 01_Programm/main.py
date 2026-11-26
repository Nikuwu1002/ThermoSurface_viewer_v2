# main.py
"""
Datei: main.py
Projekt: ThermoSurface
Kategorie: Programmeinstieg / Berechnungskoordination
Autor: Niklas King
Erstellt am: 04.10.2026
Version: 1.0

Beschreibung:
Startet ThermoSurface und vermittelt zwischen Konfiguration, HEOS und GUI.

Aufgaben:
- Lädt config.yaml und definiert die Berechnungs-Dispatcher der Diagramme.
- Initialisiert QApplication, ThermoSurfaceWindow und Ereignisschleife.
Einheiten der p-v-T-Daten: T [K], v [m³/mol], p [Pa].

Abhängigkeiten:
- Dynamisch: coolprop_heos_fvt.py, coolprop_heos_gpt.py
- Projektmodule (direkt): coolprop_heos_common.py, coolprop_heos_pvt.py,
  diagram_registry.py, gui.py, homogeneous_pvt.py, isoline_levels.py, zustandsflaeche.py
- Weitere Imports: Debug.py, coolprop_heos_hsp.py, coolprop_heos_pht.py,
  coolprop_heos_pst.py, coolprop_heos_uvt.py, limits.py
- Bibliotheken: PySide6, numpy, yaml

Verwendet von:
- Direkter Programmstart über main.py
"""

import os
os.environ["QT_API"] = "pyside6"

from homogeneous_pvt import (solve_spinodal_branch, build_homogeneous_surfaces, validate_homogeneous_boundaries, count_renderable_triangles, surface_coverage, combined_surface_report, critical_approach_report, critical_refinement_rows)
from zustandsflaeche import (PRESSURE_AXES, build_raw_field, apply_display_filter, build_substitute_field, effective_pressure_axis, ROW_OK, ROW_CRITICAL, ROW_SUPERCRITICAL, ROW_COLLAPSED)
from Debug import (configure_debug, diagnose_daten, diagnose_raster, diagnose_homogeneous, test_clapeyron, test_saturation_plateau, test_spinodale_ableitung, timed_call)
from isoline_levels import resolve_levels, record_sampling, config_signature, normalise_isoline_config
from coolprop_heos_common import critical_point, fluid_limits
from PySide6.QtWidgets import QApplication
from diagram_registry import DIAGRAMS
from gui import ThermoSurfaceWindow
from limits import MAX_POINTS
from pathlib import Path

import coolprop_heos_pvt as pvt
import coolprop_heos_pht as pht
import coolprop_heos_pst as pst
import coolprop_heos_hsp as hsp
import coolprop_heos_uvt as uvt
import importlib
import numpy as np
import logging
import yaml
import sys


# ==============================================
# 1. Konfiguration laden
# ==============================================

CONFIG_PATH = Path(__file__).resolve().parent / "config.yaml"
DEBUG_MODE = False  # Einziger Schalter: True = Diagnose, False = Normalbetrieb
LOGGER = logging.getLogger("thermosurface.main")


def lade_konfiguration():
    with CONFIG_PATH.open("r", encoding="utf-8") as datei:
        config = yaml.safe_load(datei)

    modell = str(config["thermodynamics"].get("model", "heos")).strip().lower()

    if modell != "heos":
        raise ValueError("Diese Version unterstützt nur CoolProp HEOS.")

    stoff = config["thermodynamics"]["stoff"]

    raster = config["grid_points_creation"]
    visu_raster = config["visu_points_creation"]

    # Numerische Parameter konvertieren, Auswahl des Rastertyps erhalten.
    def variablen(parameter):
        return {name: (str(wert) if name == "v_spacing" else float(wert))
                for name, wert in parameter.items()}

    raster = {name: variablen(param) for name, param in raster.items()}
    visu_raster = {name: variablen(param) for name, param in visu_raster.items()}
    pvt.validate_rasters(raster, visu_raster)

    return stoff, raster, visu_raster

def lade_diagramm_optionen(diagram_type):
    """Optionale Aufloesungen der Zusatzdiagramme; p-v-T-Config bleibt unveraendert."""
    with CONFIG_PATH.open("r", encoding="utf-8") as datei:
        config = yaml.safe_load(datei)

    options = config.get("diagram_settings", {}).get(diagram_type, {})

    if not isinstance(options, dict):
        raise ValueError(f"diagram_settings.{diagram_type} muss ein Mapping sein.")
    
    if "n_quality" in options:
        raise ValueError("n_quality ausschliesslich ueber two_phase_visu.Q_points festlegen.")
    
    return dict(options)


def lade_isolinien_vorgaben():
    """Isolinienvorgaben (diagram_settings.<id>.isolines) aller Diagramme fuer die GUI; fehlend = Standard."""
    vorgaben = {}
    for diagram_id in DIAGRAMS.ids():
        try:
            vorgaben[diagram_id] = normalise_isoline_config(lade_diagramm_optionen(diagram_id).get("isolines"))
        except ValueError as exc:
            LOGGER.warning("Isolinienvorgabe %s ungueltig (Standard verwendet): %s", diagram_id, exc)
            vorgaben[diagram_id] = normalise_isoline_config(None)
    return vorgaben


def generate_stoff_grid(stoff, raster, visu_raster):

    new_raster, new_visu, meta = pvt.generate_stoff_grid(stoff, raster, visu_raster)

    if DEBUG_MODE:
        diagnose_raster(stoff, new_raster, new_visu, meta)

    return new_raster, new_visu, meta


def _berechnen(name, function, *args, **kwargs):
    if DEBUG_MODE:
        return timed_call(name, function, *args, **kwargs)
    
    return function(*args, **kwargs)


# ==============================================
# 2. Thermodynamische Daten berechnen
# ==============================================

def berechne_daten(stoff, raster, visu_raster, progress_callback=None, isolines=None):
    """Berechnet wie bisher; optional phasenbasierte Statusmeldungen.

    isolines (TS-VP03): Mapping {isothermen|isochoren|isobaren: {start, stop, factor}} fuer die geometrisch
    verteilten Isolinienwerte; None = bisheriges Raster (T_step/v_step/p_step). `berechne_auswahl` setzt dafuer
    die Vorgabe aus diagram_settings.pvt.isolines (config.yaml) ein. Der Oberflaechenrasterpfad bleibt unberuehrt.
    """

    def progress(value, phase):
        if progress_callback is not None:
            progress_callback(int(value), str(phase))

    PROGRESS_STEPS = {
        "validierung": (2, "Prüfe Eingabeparameter ..."),
        "raster": (6, "Lese Stoffdaten und erstelle Raster ..."),
        "surface": (10, "Berechne Zustandsoberfläche (HEOS) ..."),
        "saturation": (55, "Berechne Sättigungszustände ..."),
        "phase_mask": (60, "Berechne Phasenmaske ..."),
        "isothermen": (65, "Berechne Isothermen ..."),
        "isochoren": (72, "Berechne Isochoren ..."),
        "isobaren": (78, "Berechne Isobaren ..."),
        "binodale": (84, "Berechne Binodalen ..."),
        "two_phase": (87, "Berechne Zweiphasenfläche ..."),
        "spinodale": (90, "Berechne Spinodalen ..."),
        "results": (93, "Bereite Ergebnisdaten auf ..."),
        "completed": (95, "Thermodynamische Berechnung abgeschlossen"),
    }

    progress(*PROGRESS_STEPS["validierung"])
    pvt.validate_rasters(raster, visu_raster)
    pvt.validate_heos_temperature_ranges(stoff, raster, visu_raster)
    progress(*PROGRESS_STEPS["raster"])

    # Zustandsfläche: dichteres Berechnungsraster
    # Kritische Daten
    kritisch = _berechnen("Kritischer Punkt", critical_point, stoff)
    Tc, pc, vc = kritisch
    Ttr, T_backend_min, T_backend_max = fluid_limits(stoff)

    # AP2: dieselben Temperaturstuetzstellen fuer Flaeche und Saettigung.
    cfg = raster["isotherm_grid"]
    sat_cfg = visu_raster["binodale_visu"]
    T_all, T_sat = pvt.build_pvt_temperature_grid(cfg, sat_cfg, Tc, Ttr, T_backend_min, T_backend_max)

    progress(*PROGRESS_STEPS["surface"])
    # Zustandsfläche auf dem gemeinsamen Raster berechnen
    surface = _berechnen("Zustandsoberfläche", pvt.multi_isotherm, **cfg, stoff=stoff, T_values=T_all)

    progress(*PROGRESS_STEPS["saturation"])

    # Saettigung exakt auf den unterkritischen Flaechentemperaturen.
    binodal_liquid, binodal_vapor = _berechnen("Binodale", pvt.multi_binodale, **sat_cfg, stoff=stoff, T_values=T_sat)

    if not np.array_equal(binodal_liquid[:, 0], T_sat):
        raise ValueError("HEOS konnte nicht sämtliche erforderlichen " "Sättigungszustände berechnen.")

    progress(*PROGRESS_STEPS["phase_mask"])
    
    phase_mask = _berechnen("Zweiphasenmaske", pvt.berechne_zweiphasenmaske, stoff, surface, binodal_liquid, binodal_vapor)

    
    # Isolinien: separate Visualisierungsparameter
    progress(*PROGRESS_STEPS["isothermen"])
    isotherm_cfg = dict(visu_raster["isotherm_visu"])
    isotherm_count = isotherm_cfg.pop("n_isotherms", None)
    isotherm_temperatures = None
    isoline_config = isolines
    isoline_levels, isoline_report = None, None

    if isoline_config is not None:
        # TS-VP03: konstante T-/v-/p-Werte aus der geometrischen Progression; Bereichsgrenzen der bisherigen
        # Visualisierungsparameter bilden Automatikbereich und Gueltigkeitsfenster.
        isoline_levels, isoline_report = resolve_levels(isoline_config, {
            "isothermen": (float(isotherm_cfg["T_start_isotherm"]), float(isotherm_cfg["T_stop_isotherm"])),
            "isochoren": (float(visu_raster["isochore_visu"]["v_start_isochore"]), float(visu_raster["isochore_visu"]["v_stop_isochore"])),
            "isobaren": (float(visu_raster["isobare_visu"]["p_start_isobare"]), float(visu_raster["isobare_visu"]["p_stop_isobare"]))})
        isotherm_temperatures = isoline_levels["isothermen"]

    elif isotherm_count is not None:
        isotherm_temperatures = np.linspace(float(isotherm_cfg["T_start_isotherm"]), float(isotherm_cfg["T_stop_isotherm"]), int(isotherm_count))

    if isoline_levels is not None and not len(isotherm_temperatures):
        isothermen = np.empty((0, 0, 3))
    else:
        isothermen = _berechnen("Isothermen", pvt.multi_isotherm, **isotherm_cfg, stoff=stoff, T_values=isotherm_temperatures)

    progress(*PROGRESS_STEPS["isochoren"])

    if isoline_levels is None:
        isochoren = _berechnen("Isochoren", pvt.multi_isochore, **visu_raster["isochore_visu"], stoff=stoff,)
    elif not len(isoline_levels["isochoren"]):
        isochoren = np.empty((0, 0, 3))
    else:
        isochoren = _berechnen("Isochoren", pvt.multi_isochore, **visu_raster["isochore_visu"], stoff=stoff, v_values=isoline_levels["isochoren"])

    progress(*PROGRESS_STEPS["isobaren"])

    if isoline_levels is None:
        isobaren = _berechnen("Isobaren", pvt.multi_isobare, **visu_raster["isobare_visu"], stoff=stoff,)
    else:
        isobaren = _berechnen("Isobaren", pvt.multi_isobare, **visu_raster["isobare_visu"], stoff=stoff, p_values=isoline_levels["isobaren"])
    if isoline_report is not None:
        for family, curves in (("isothermen", isothermen), ("isochoren", isochoren), ("isobaren", isobaren)):
            record_sampling(isoline_report, family, curves)

    # AP2: Die bereits auf T_sat berechneten Binodalen wiederverwenden.
    progress(*PROGRESS_STEPS["binodale"])

    # Gleichgewichts-Zweiphasenfläche
    progress(*PROGRESS_STEPS["two_phase"])
    Q_raw = visu_raster["two_phase_visu"]["Q_points"]

    if (not np.isfinite(Q_raw) or Q_raw < 2 or not float(Q_raw).is_integer()):
        raise ValueError("Q_points muss eine ganze Zahl >= 2 sein.")

    n_quality = int(Q_raw)

    two_phase_reaches_critical = True
    if len(binodal_liquid) >= 2 and len(binodal_vapor) >= 2:
        two_phase_surface, Q_values = _berechnen("Zweiphasenfläche", pvt.calc_two_phase_surface, binodal_liquid, binodal_vapor, n_quality=n_quality)
    else:
        LOGGER.warning("Keine ausreichenden Sättigungspunkte für %s; Zweiphasenfläche deaktiviert", stoff)
        two_phase_surface = np.empty((0, 0, 3))
        Q_values = np.empty(0)

    progress(*PROGRESS_STEPS["spinodale"])

    # AP3: Filter vor Endpunktergaenzung; kritischen Punkt nicht ausserhalb
    # der benutzerseitigen Temperaturgrenzen oder ohne Annaherungsnachweis
    # an eine unvollstaendige Spinodale anhaengen.
    spin_cfg = visu_raster["spinodale_visu"]
    T_min = float(spin_cfg["T_start_spinodal"])
    T_max = float(spin_cfg["T_stop_spinodal"])
    raw_liquid = np.empty((0, 3), dtype=float)
    raw_vapor = np.empty((0, 3), dtype=float)
    try:
        raw_liquid, raw_vapor = _berechnen("Spinodale", pvt.calc_spinodale, stoff)
    except (ValueError, RuntimeError, NotImplementedError) as exc:
        LOGGER.warning("HEOS-Spinodale fuer %s nicht verfuegbar: %s", stoff, exc)

    # Der bestehende CoolProp-Tracer bleibt als Referenz und Fallback bestehen.
    # Das neue HEOS-Dichtescanning ist getrennt konfiguriert und aendert
    # nur den Dampfzweig; fehlende Temperaturzeilen bleiben NaN.
    spinodal_solver = {"selected": "legacy", "fallback_reason": None}
    spin_grid = raster.get("spinodale_grid", {})

    if "n_T" in spin_grid:  # alte Konfiguration/Mock-Integration bleibt kompatibel
        requested_low = max(float(spin_grid["T_start_spinodal"]), T_min)
        requested_high = min(float(spin_grid["T_stop_spinodal"]), T_max)
        if requested_low < requested_high:
            try:
                new_vapor, numerical_report = _berechnen("Dampfspinodale (Temperatur-HEOS)", pvt.calc_spinodale_vapor_temperature, stoff, requested_low, requested_high,
                    n_T=int(spin_grid["n_T"]), n_scan=int(spin_grid.get("n_scan", 200)), residual_tolerance=float(spin_grid.get("residual_tolerance", 1e-8)), with_diagnostics=True)
                
                spinodal_solver["temperature"] = numerical_report
                overlap = pvt.compare_spinodal_overlap(new_vapor, raw_vapor, critical_T=Tc)
                spinodal_solver["overlap"] = overlap

                # Beide Methoden muessen auf identischem Ast liegen.
                # Interpolation ist AUSSCHLIESSLICH ein Vergleich, keine Kurvenerzeugung.
                match = (overlap["compared"] < 3 or (overlap["max_rel_volume"] <= 0.06 and overlap["max_rel_pressure"] <= 0.03))

                old_filtered = pvt.filter_spinodale(raw_vapor, T_min, T_max)
                old_valid = old_filtered[np.isfinite(old_filtered).all(axis=1) & (old_filtered[:, 1] > 0) & (old_filtered[:, 2] > 0)]

                valid_new = new_vapor[np.isfinite(new_vapor).all(axis=1) & (new_vapor[:, 1] > 0) & (new_vapor[:, 2] > 0)]

                # Nicht schon wegen weniger neuer Punkte den besser
                # aufgeloesten alten Ast ersetzen; der untere Temperaturbereich
                # muss tatsaechlich erschlossen und der alte abgedeckt sein.
                old_end_covered = (not len(old_valid) or (len(valid_new) and valid_new[-1, 0] >= old_valid[-1, 0] - 1.5 * (requested_high-requested_low) / (int(spin_grid["n_T"])-1)))
                
                accepts_new = (len(valid_new) >= max(3, len(old_valid), 10) and numerical_report["covers_start"] and numerical_report["success_fraction"] >= 0.90 and not numerical_report["ambiguous_temperatures"] and old_end_covered and match)
                
                if accepts_new:
                    raw_vapor = new_vapor
                    spinodal_solver["selected"] = "temperature_heos"
                else:
                    spinodal_solver["fallback_reason"] = ("Neuer Ast nicht freigegeben: Abdeckung/Mehrfachmaxima/" "Legacy-Uebereinstimmung nicht ausreichend")

            except (ValueError, RuntimeError, NotImplementedError) as exc:
                spinodal_solver["fallback_reason"] = f"Neuer Solver fehlgeschlagen: {exc}"
                LOGGER.warning("Dampfspinodale: %s", exc)
        else:
            spinodal_solver["fallback_reason"] = "Berechnungs-/Sichtbarkeitsbereich ohne Schnittmenge"

    spinodale_liquid, liquid_report = pvt.prepare_spinodal_branch(raw_liquid, kritisch, T_min, T_max, "liquid")
    spinodale_vapor, vapor_report = pvt.prepare_spinodal_branch(raw_vapor, kritisch, T_min, T_max, "vapor")

    vapor_report["source"] = spinodal_solver["selected"]

    if spinodal_solver.get("temperature"):
        vapor_report["grid_missing"] = int(spinodal_solver["temperature"]["requested_points"] - spinodal_solver["temperature"]["valid_points"])

        if spinodal_solver["selected"] == "temperature_heos" and vapor_report["grid_missing"]:
            vapor_report["geometry_complete"] = False

    if spinodal_solver.get("fallback_reason"):
        LOGGER.warning("Dampfspinodale: %s", spinodal_solver["fallback_reason"])

    # Die Residualdiagnostik ist AP3-relevant und laeuft auch ohne DEBUG_MODE;
    # sie darf mangelhafte Spinodalen nie als gueltig ausweisen.
    spinodal_derivatives = test_spinodale_ableitung(stoff, spinodale_liquid, spinodale_vapor, kritisch)
    spinodal_validation = {"liquid": liquid_report, "vapor": vapor_report}

    for branch in ("liquid", "vapor"):
        spinodal_validation[branch]["derivative"] = spinodal_derivatives[branch]
        spinodal_validation[branch]["complete"] = bool(spinodal_validation[branch]["geometry_complete"] and spinodal_derivatives[branch]["passed"])

        if not spinodal_validation[branch]["complete"]:
            LOGGER.warning("AP3 Spinodale %s nicht vollstaendig validiert: %s", branch, spinodal_validation[branch])

    # Binodale: AP2 setzt einen geometrisch exakt gemeinsamen Tc-Endpunkt.
    # Dieser Vorgang ist von der AP3-Spinodalbewertung getrennt.
    progress(*PROGRESS_STEPS["results"])
    kritisch_tvp = np.array([[Tc, vc, pc]])

    def ergaenze_kritischen_punkt(kurve):
        kurve = np.asarray(kurve, dtype=float).reshape(-1, 3)
        if len(kurve) < 2:
            return kurve
        
        kurve = kurve[kurve[:, 0] < Tc]

        return np.vstack((kurve, kritisch_tvp))

    binodal_liquid = ergaenze_kritischen_punkt(binodal_liquid)
    binodal_vapor = ergaenze_kritischen_punkt(binodal_vapor)

    daten = {
        "surface": surface,
        "isothermen": isothermen,
        "isochoren": isochoren,
        "isobaren": isobaren,
        "binodal_liquid": binodal_liquid,
        "binodal_vapor": binodal_vapor,
        "spinodale_liquid": spinodale_liquid,
        "spinodale_vapor": spinodale_vapor,
        "spinodal_validation": spinodal_validation,
        "spinodal_solver": spinodal_solver,
        "critical_point": kritisch,
        "phase_mask": phase_mask,
        "two_phase_surface": two_phase_surface,
        "two_phase_reaches_critical": two_phase_reaches_critical,
        "Q_values": Q_values,
        "isoline_report": isoline_report,}

    # AP13--AP15: homogene HEOS-Fortsetzung, strikt unabhaengig vom
    # Gleichgewichts-Dampfgehaltsnetz. Jede T-Zeile stammt direkt aus T_sat;
    # kein Interpolieren von Spinodalen ueber fehlende Zustaende.
    options = lade_diagramm_optionen("pvt")
    pressure_axis = options.get("pressure_axis", "log")

    if pressure_axis not in PRESSURE_AXES:
        raise ValueError("diagram_settings.pvt.pressure_axis: " + ", ".join(PRESSURE_AXES) + " erforderlich")
    
    requested_axis = pressure_axis
    auto_axis = bool(options.get("auto_pressure_axis", True))
    daten["pressure_axis"] = pressure_axis
    daten["pressure_axis_notice"] = None

    enabled = bool(options.get("include_homogeneous", True))
    n_hom = int(options.get("homogeneous_points", 20))
    unstable = bool(options.get("include_unstable", False))

    include_replacement = bool(options.get("include_replacement", True))
    n_sub = int(options.get("replacement_points", 40))

    unstable_scan = int(options.get("unstable_scan_points", 256))
    n_refine = int(options.get("critical_refinement_rows", 0))
    close_critical = bool(options.get("close_critical_point", False))

    hom_liquid, hom_vapor = binodal_liquid[:-1], binodal_vapor[:-1]
    refinement_note = None

    if enabled and n_refine > 0 and len(hom_liquid):
        # AP20.6: zusaetzliche Temperaturzeilen nur fuer die homogenen Flaechen.
        try:
            extra_T, extra_l, extra_g = critical_refinement_rows(stoff, float(hom_liquid[-1, 0]), n_refine)

            if len(extra_T):
                hom_liquid = np.vstack((hom_liquid, extra_l))
                hom_vapor = np.vstack((hom_vapor, extra_g))

        except (ValueError, RuntimeError, ImportError) as exc:
            refinement_note = str(exc)
            LOGGER.warning("Kritische Verfeinerung nicht moeglich: %s", exc)

    sat_T = hom_liquid[:, 0].copy()

    if not 3 <= n_hom <= 1000 or len(sat_T) * n_hom * 3 > MAX_POINTS:
        raise ValueError("Homogene Zusatzflaechen: ungueltiges Punktbudget")
    
    missing = np.column_stack((sat_T, np.full((len(sat_T), 2), np.nan)))
    boundaries = {"liquid": missing.copy(), "vapor": missing.copy()}

    hom_reports = {"enabled": enabled, "available": False, "state": "disabled",
                   "spinodals": {}, "surfaces": {}, "error": None,
                   "pressure_axis": pressure_axis, "temperature_rows": int(len(sat_T)),
                   "critical_refinement_rows": int(len(sat_T) - len(binodal_liquid[:-1])),
                   "critical_refinement_note": refinement_note}
    
    empty = np.full((len(sat_T), n_hom, 3), np.nan)
    new_surfaces = {name: empty.copy() for name in ("metastable_liquid_surface", "metastable_vapor_surface", "unstable_surface")}
    daten["critical_closure"] = {"requested": close_critical, "supported": False, "applied": False}
    n_sub = max(3, min(n_sub, 1000))
    daten["replacement_surface"] = np.full((len(sat_T), n_sub, 3), np.nan)
    daten["replacement_field"] = None
    daten["homogeneous_raw_field"] = None
    daten["replacement_status"] = {"state": "disabled", "rows": {}, "not_constructible": [], "notice": None}

    if enabled:
        try:
            if not 3 <= n_hom <= 1000:
                raise ValueError("homogeneous_points muss 3..1000 sein")
            
            scan = int(raster.get("spinodale_grid", {}).get("n_scan", 200))
            tol = float(raster.get("spinodale_grid", {}).get("residual_tolerance", 1e-8))

            # Beide unabhängigen Aeste und die kritische Umgebung pruefen.
            for branch in ("liquid", "vapor"):
                curve, diagnostics = _berechnen("AP13: " + branch + " Spinodale", solve_spinodal_branch, stoff, sat_T, branch, n_scan=scan, residual_tolerance=tol)
                boundaries[branch] = curve
                hom_reports["spinodals"][branch] = diagnostics

            surfaces, reports = _berechnen("AP15: homogene Flaechen", build_homogeneous_surfaces, stoff, hom_liquid, hom_vapor, boundaries["liquid"], boundaries["vapor"], n_points=n_hom, include_unstable=unstable, scan_points=unstable_scan)
            new_surfaces.update(surfaces)

            hom_reports["surfaces"] = reports
            hom_reports["boundary_check"] = validate_homogeneous_boundaries(new_surfaces, hom_liquid, hom_vapor, boundaries["liquid"], boundaries["vapor"], kritisch)

            # Zustandsflaeche: Rohfeld (HEOS unveraendert, Filter setzen nur den Darstellungsstatus)
            # und gekennzeichnete Ersatzflaeche zwischen den verifizierten HEOS-Spinodalen.
            if unstable:
                raw_field = build_raw_field(stoff, sat_T, boundaries["liquid"], boundaries["vapor"], n_hom)
                apply_display_filter(raw_field, boundaries["liquid"], boundaries["vapor"], float(kritisch[1]))
                daten["homogeneous_raw_field"] = raw_field
                hom_reports["raw_field"] = raw_field.counts()

            if include_replacement:
                sub = build_substitute_field(sat_T, boundaries["liquid"], boundaries["vapor"], kritisch, n_sub)
                daten["replacement_field"] = sub
                daten["replacement_surface"] = sub.raw_points()
                bad = [(float(sat_T[i]), sub.row_reason.get(i, ""))
                       
                       for i in range(len(sat_T)) if sub.row_status[i] not in (
                           ROW_OK, ROW_CRITICAL, ROW_SUPERCRITICAL, ROW_COLLAPSED)]
                
                state = ("not_constructible" if not np.any(sub.row_status == ROW_OK)
                         else "partial" if bad else "ok")
                
                notice = None

                if state != "ok":
                    notice = (f"Ersatzflaeche nicht konstruierbar fuer {len(bad)} von {len(sat_T)} " f"Temperaturzeilen (z. B. T={bad[0][0]:.4g} K: {bad[0][1]})" if bad else "Ersatzflaeche nicht konstruierbar")
                    LOGGER.warning(notice)

                daten["replacement_status"] = {"state": state, "rows": sub.counts()["rows"], "not_constructible": bad, "notice": notice}
                hom_reports["replacement"] = daten["replacement_status"]

            # Druckachse: p<=0 darf nicht ausgeblendet werden; Log nur fuer p>0, sonst Umschaltung.
            pressure_axis, axis_notice = effective_pressure_axis(pressure_axis, [new_surfaces[k] for k in ("metastable_liquid_surface", "metastable_vapor_surface", "unstable_surface")] + [daten["replacement_surface"]], auto=auto_axis)
            
            if axis_notice:
                daten["pressure_axis"] = pressure_axis
                daten["pressure_axis_notice"] = axis_notice
                LOGGER.warning(axis_notice)

            hom_reports["pressure_axis"] = pressure_axis
            hom_reports["pressure_axis_requested"] = requested_axis

            # Es zaehlen nur Gitterzellen mit >=3 benachbarten gueltigen Knoten.
            # Diese Anzahl ist VTK-unabhaengig und vermeidet falsch-positive
            # Verfuegbarkeit bei voneinander isolierten Temperaturzeilen.
            for name, info in reports.items():
                info["candidate_triangles"] = count_renderable_triangles( new_surfaces[name], pressure_axis)
                info["coverage"] = surface_coverage(new_surfaces[name], info.get("gap_rows"), pressure_axis)

            hom_reports["available"] = any(info["candidate_triangles"] > 0 for info in reports.values())
            hom_reports["state"] = ("available" if hom_reports["available"] else "no_renderable_cells")

            # AP20.7: Netz-/Naht-/Komponentenstatistik der zusammengefuehrten Flaeche.
            hom_reports["combined"] = combined_surface_report(new_surfaces, "asinh", reports["unstable_surface"].get("gap_rows"))

            # AP20.6: Schliessung am kritischen Punkt nur bei validierter Annaeherung.
            approach = critical_approach_report(hom_liquid, boundaries["liquid"], boundaries["vapor"], hom_vapor, kritisch)
            hom_reports["critical_approach"] = approach

            gap_rows = reports["unstable_surface"].get("gap_rows", {})

            shown = [n for n in new_surfaces if n != "unstable_surface" or unstable]
            last_row_complete = bool(all(np.isfinite(new_surfaces[n][-1]).all() for n in shown))
            closure_ok = bool(approach["closure_supported"] and last_row_complete)
            daten["critical_closure"] = {"requested": close_critical, "supported": bool(approach["closure_supported"]),"last_row_complete": last_row_complete, "applied": bool(close_critical and closure_ok)}
            
            # AP3: Die bisherigen sichtbaren, unabhaengig validierten Linien
            # werden NICHT ohne separate Freigabe durch die neuen ersetzt.
        except (ValueError, RuntimeError, OverflowError, ImportError) as exc:
            hom_reports["error"] = str(exc)
            hom_reports["state"] = "failed"

            if DEBUG_MODE:
                LOGGER.exception("Homogene Erweiterung fehlgeschlagen")
            else:
                LOGGER.warning("Homogene Erweiterung fehlgeschlagen: %s", exc)

    daten.update(new_surfaces)
    daten["homogeneous_spinodal_liquid"] = boundaries["liquid"]
    daten["homogeneous_spinodal_vapor"] = boundaries["vapor"]
    daten["homogeneous_validation"] = hom_reports

    # AP2: auch ohne Debug-Modus inkonsistente Phasengrenzen ablehnen.
    pvt.validate_pvt_phase_boundaries(daten)

    if DEBUG_MODE:
        diagnose_daten(daten)
        diagnose_homogeneous(daten)
        # Physikalischer Test nur an einem vorhandenen Sättigungspunkt.
        if len(binodal_liquid) >= 2:
            T_test = float(binodal_liquid[len(binodal_liquid) // 2, 0])
            sat_test = binodal_liquid[binodal_liquid[:, 0] < Tc]

            if len(sat_test) > 0:

                T_test = float(sat_test[len(sat_test) // 2, 0])

                test_clapeyron(stoff, T_test)

                test_saturation_plateau(stoff, T_test, pvt.calc_saturation_point, pvt.p_eos)
        # AP3-Ableitungsdiagnostik laeuft bereits im normalen Berechnungspfad.

    # Anzeigegrenzen der Zweiphasenfläche (Auto-Werte: p_sat(T_start_sat) bis pc);
    # erst nach der Validierung, die die vollständige Fläche voraussetzt.
    two_phase_cfg = visu_raster["two_phase_visu"]
    p_lo = two_phase_cfg.get("p_start_twophase")
    p_hi = two_phase_cfg.get("p_stop_twophase")
    if (p_lo is not None or p_hi is not None) and daten["two_phase_surface"].size:
        daten["two_phase_surface"], daten["two_phase_reaches_critical"] = pvt.clip_two_phase_surface(daten["two_phase_surface"], p_lo, p_hi, pc)

        if not daten["two_phase_surface"].size:
            LOGGER.warning("Zweiphasenfläche liegt vollständig außerhalb der Druckgrenzen [%s, %s] Pa", p_lo, p_hi)

    progress(95, "Thermodynamische Berechnung abgeschlossen")

    return daten

def _effektive_isolinien(optionen, isolines):
    """Wirksame Isolinienvorgabe: Aufgabenvorgabe (GUI) vor config.yaml; None = keine (bisheriges Raster)."""
    return isolines if isolines is not None else optionen.get("isolines")


def _berechne_pvt(stoff, raster, visu_raster, progress_callback=None, isolines=None):
    cfg = _effektive_isolinien(lade_diagramm_optionen("pvt"), isolines)
    return berechne_daten(stoff, raster, visu_raster, progress_callback=progress_callback, isolines=cfg)


def _berechne_zusatzdiagramm(label, modul, funktion, optionen_id=None):
    """Gemeinsame Huelle der Zusatzdiagramme (eigenes Raster, keine p-v-T-Berechnung).

    `modul` ist ein Modul oder ein Modulname; Namen werden erst beim Aufruf importiert (g-p-T, f-v-T), damit
    die isolierten Ladetests von main (mit ersetzten HEOS-Modulen) unveraendert bleiben.
    """

    def berechne(stoff, raster, visu_raster, progress_callback=None, isolines=None):
        if progress_callback is not None:
            progress_callback(5, f"Berechne {label}-Zustandsflächen ...")

        quelle = importlib.import_module(modul) if isinstance(modul, str) else modul
        optionen = lade_diagramm_optionen(optionen_id) if optionen_id else {}

        if isolines is not None:
            optionen["isolines"] = isolines

        daten = _berechnen(f"{label}-Zustandsflächen", getattr(quelle, funktion), stoff=stoff, n_quality=int(visu_raster["two_phase_visu"]["Q_points"]), **optionen)

        if progress_callback is not None:
            progress_callback(95, f"{label}-Berechnung abgeschlossen")

        return daten

    return berechne


def _cache_params_pvt(raster, visu_raster, isolines=None):
    optionen = lade_diagramm_optionen("pvt")
    out = {"raster": raster, "visu_raster": visu_raster, "diagram_settings": optionen}
    cfg = _effektive_isolinien(optionen, isolines)

    if cfg is not None:
        out["isolines"] = config_signature(cfg)

    return out


def _cache_params_zusatz(optionen_id=None):
    """Nur tatsaechlich verwendete Parameter: Q_points (+ diagram_settings des Typs, Isolinienlevels)."""

    def params(raster, visu_raster, isolines=None):
        out = {"Q_points": int(visu_raster["two_phase_visu"]["Q_points"])}
        if optionen_id:
            optionen = lade_diagramm_optionen(optionen_id)
            out["diagram_settings"] = optionen
            cfg = _effektive_isolinien(optionen, isolines)
            if cfg is not None:
                out["isolines"] = config_signature(cfg)

        return out

    return params


# Berechnungsfunktion je Registry-ID. Ein neues Diagramm braucht einen Registry-Eintrag
# (diagram_registry.py), einen Eintrag hier und eine cache_params-Funktion.
BERECHNUNGEN = {
    "pvt": _berechne_pvt,
    "pht": _berechne_zusatzdiagramm("p-h-T", pht, "berechne_pht_daten", "pht"),
    "pst": _berechne_zusatzdiagramm("p-s-T", pst, "berechne_pst_daten", "pst"),
    "hsp": _berechne_zusatzdiagramm("h-s-p", hsp, "berechne_hsp_daten", "hsp"),
    "uvt": _berechne_zusatzdiagramm("u-v-T", uvt, "berechne_uvt_daten", "uvt"),
    "gpt": _berechne_zusatzdiagramm("g-p-T", "coolprop_heos_gpt", "berechne_gpt_daten", "gpt"),
    "fvt": _berechne_zusatzdiagramm("f-v-T", "coolprop_heos_fvt", "berechne_fvt_daten", "fvt"),
}
DIAGRAMS.attach("pvt", cache_params=_cache_params_pvt)
DIAGRAMS.attach("pht", cache_params=_cache_params_zusatz("pht"))
DIAGRAMS.attach("pst", cache_params=_cache_params_zusatz("pst"))
DIAGRAMS.attach("hsp", cache_params=_cache_params_zusatz("hsp"))
DIAGRAMS.attach("uvt", cache_params=_cache_params_zusatz("uvt"))
DIAGRAMS.attach("gpt", cache_params=_cache_params_zusatz("gpt"))
DIAGRAMS.attach("fvt", cache_params=_cache_params_zusatz("fvt"))
assert set(BERECHNUNGEN) == set(DIAGRAMS.ids()), "Registry und BERECHNUNGEN muessen uebereinstimmen"


def berechne_auswahl(stoff, raster, visu_raster, diagram_type="pvt", progress_callback=None, isolines=None):
    """
    Wählt ausschließlich die Berechnung
    des angeforderten Diagrammtyps aus (Dispatch über die Registry-IDs).
    `isolines`: optionale Isolinienvorgabe (start/stop/factor je Familie); None = config.yaml.
    """
    DIAGRAMS.get(diagram_type)  # ValueError bei unbekanntem Typ
    return BERECHNUNGEN[diagram_type](stoff, raster, visu_raster, progress_callback=progress_callback, isolines=isolines)

# ==============================================
# 3. Hauptprogramm
# ==============================================

def main():
    configure_debug(DEBUG_MODE)
    app = QApplication(sys.argv)

    stoff, raster, visu_raster = lade_konfiguration()
    window = ThermoSurfaceWindow(
        start_stoff=stoff,
        surface_cfg=raster["isotherm_grid"],
        raster=raster,
        visu_raster=visu_raster,
        berechnen=berechne_auswahl,
        grid_generator=generate_stoff_grid,
        debug_mode=DEBUG_MODE,
        isoline_defaults=lade_isolinien_vorgaben())
    
    window.show()

    return app.exec()

if __name__ == "__main__":
    sys.exit(main())