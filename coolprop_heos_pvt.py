"""
coolprop_heos.py

HEOS-basierte Berechnungen für p-v-T-Zustandsdiagramme.

Verwendet wird das Helmholtz-Energie-Referenzmodell "HEOS" von CoolProp
verwendet (AbstractState("HEOS", stoff)). HEOS ist für sehr viele Reinstoffe
(z.B. Water, CO2, Nitrogen, Methane, ...) eine hochgenaue Referenz-
Zustandsgleichung und übernimmt intern bereits:

- die Auswertung der Helmholtz-Energie-Zustandsgleichung für p, v, rho
- die Bestimmung der Phasengrenze (Sättigungslinie)
- die Berechnung des Phasengleichgewichts
- die Bestimmung des kritischen Punkts

Die bestehenden Funktionssignaturen und Rückgabeformate bleiben unverändert.

Einheitenkonvention (wie in cubic_eos.py, SI-Basis):
    T            Temperatur                 [K]
    v            molares Volumen            [m^3/mol]
    rho_molar    molare Dichte = 1/v         [mol/m^3]
    p            Druck                      [Pa]
"""

import CoolProp.CoolProp as CP
import numpy as np
import logging

from copy import deepcopy
from limits import MAX_POINTS
from coolprop_heos_common import _new_state, fluid_limits, critical_point
from spinodal_vapor_solver import calc_spinodale_vapor_temperature, compare_spinodal_overlap

LOGGER = logging.getLogger("thermosurface.coolprop_heos")

# ============================================================
# p-v-T: Rasterverwaltung und Eingabevalidierung
# ============================================================
def validate_rasters(raster, visu_raster):
    """Gemeinsame Eingabevalidierung für automatische und manuelle Raster.

    Führt keine thermodynamischen HEOS-Rechnungen durch.
    """
    for group_name, groups in (("Berechnung", raster), ("Visualisierung", visu_raster)):
        for name, params in groups.items():
            for key, value in params.items():
                if key == "v_spacing":
                    if value not in ("linear", "log"):
                        raise ValueError(f"{group_name}.{name}: v_spacing muss linear oder log sein.")
                    continue

                if not np.isfinite(value):
                    raise ValueError(f"{group_name}.{name}.{key} muss endlich sein.")
                if key.startswith(("T_", "v_", "p_")) and key != "v_points" and value <= 0:
                    raise ValueError(f"{group_name}.{name}.{key} muss positiv sein.")

                if key == "Q_points":
                    if (value < 2 or value > MAX_POINTS or not float(value).is_integer()):
                        raise ValueError(f"Q_points muss ganzzahlig zwischen " f"2 und {MAX_POINTS} liegen.")

                if key == "v_points":
                    if (value < 2 or value > MAX_POINTS or not float(value).is_integer()):
                        raise ValueError( f"{name}.v_points muss ganzzahlig zwischen " f"2 und {MAX_POINTS} liegen.")

                if key == "n_isotherms":
                    if (value < 2 or value > MAX_POINTS or not float(value).is_integer()):
                        raise ValueError(
                            f"{name}.n_isotherms muss ganzzahlig zwischen 2 und {MAX_POINTS} liegen.")

                if "_start_" in key:
                    stop_key = key.replace("_start_", "_stop_", 1)
                    if stop_key in params and value >= params[stop_key]:
                        raise ValueError(f"{group_name}.{name}: {key} muss kleiner als {stop_key} sein.")

            # Punktzahl ohne np.arange() berechnen, um Speicherexzesse zu verhindern.
            if name in ("isotherm_grid", "isotherm_visu"):
                n_T = (int(params["n_isotherms"]) if name == "isotherm_visu" and "n_isotherms" in params
                       else int(np.ceil((params["T_stop_isotherm"] - params["T_start_isotherm"]) / params["T_step_isotherm"]))
                       + int(name == "isotherm_grid"))
                n_v = (int(params["v_points"]) if params.get("v_spacing") == "log"
                       else int(np.ceil((params["v_stop_isotherm"] - params["v_start_isotherm"]) / params["v_step_isotherm"])))

            elif name in ("isochore_grid", "isochore_visu"):
                n_T = int(np.ceil((params["T_stop_isochore"] - params["T_start_isochore"]) / params["T_step_isochore"]))
                n_v = (int(params["v_points"]) if params.get("v_spacing") == "log"
                       else int(np.ceil((params["v_stop_isochore"] - params["v_start_isochore"]) / params["v_step_isochore"])))
                
                max_points = 10000
            elif name in ("isobare_grid", "isobare_visu"):
                n_v = int(np.ceil((params["p_stop_isobare"] - params["p_start_isobare"]) / params["p_step_isobare"]))
                n_T = int(np.ceil((params["T_stop_isobare"] - params["T_start_isobare"]) / params["T_step_isobare"]))

                max_points = 10000
            elif name in ("binodale_grid", "binodale_visu"):
                n_T = int(np.ceil((params["T_stop_sat"] - params["T_start_sat"]) / params["T_step_sat"]))
                n_v = 1

            elif name == "two_phase_visu":
                continue

            elif name in ("spinodale_grid", "spinodale_visu"):
                if name == "spinodale_grid" and "n_T" in params:
                    n_T_spin, n_scan = params["n_T"], params.get("n_scan", 200)
                    tol = params.get("residual_tolerance", 1e-8)
                    if (not float(n_T_spin).is_integer() or not 2 <= n_T_spin <= 10000
                            or not float(n_scan).is_integer() or not 16 <= n_scan <= 10000
                            or n_T_spin * n_scan > 2000000 or not 0 < tol <= 1e-2):
                        raise ValueError("Spinodalraster: n_T=2..10000, n_scan=16..10000, "
                                         "n_T*n_scan<=2000000, 0<residual_tolerance<=1e-2")
                continue

            else:
                raise ValueError(f"Unbekannte Rastergruppe: {name}")
            
            if n_T < 2 or (n_v < 2 and name not in ("binodale_grid", "binodale_visu")):
                raise ValueError(f"{name}: mindestens zwei Werte je variabler " f"Achse erforderlich.")

            if n_T * n_v > MAX_POINTS:
                raise ValueError(f"{name}: {n_T * n_v:,} Rasterpunkte " f"überschreiten {MAX_POINTS:,}.")
            
            sat = visu_raster["binodale_visu"]
            n_sat = int(np.ceil((sat["T_stop_sat"] - sat["T_start_sat"]) / sat["T_step_sat"]))
            n_Q = int(visu_raster["two_phase_visu"]["Q_points"])

            if n_sat * n_Q > MAX_POINTS:
                raise ValueError("Zweiphasenraster überschreitet " f"{MAX_POINTS:,} Punkte.")

def validate_heos_temperature_ranges(stoff, raster, visu_raster):
    """Prüft berechnete Einzelphasenraster gegen die HEOS-Temperaturgrenzen."""
    _triple, Tmin, Tmax = fluid_limits(stoff)

    checks = (
        ("isotherm_grid", raster["isotherm_grid"], "T_start_isotherm", "T_stop_isotherm"),
        ("isotherm_visu", visu_raster["isotherm_visu"], "T_start_isotherm", "T_stop_isotherm"),
        ("isochore_visu", visu_raster["isochore_visu"], "T_start_isochore", "T_stop_isochore"),
        ("isobare_visu", visu_raster["isobare_visu"], "T_start_isobare", "T_stop_isobare"),)
    
    for name, cfg, lower, upper in checks:
        if cfg[lower] < Tmin or cfg[upper] > Tmax:
            raise ValueError(f"{name}: [{cfg[lower]:.4g}, {cfg[upper]:.4g}] K liegt " f"außerhalb des HEOS-Gültigkeitsbereichs [{Tmin:.4g}, {Tmax:.4g}] K für {stoff}.")


def generate_stoff_grid(stoff, raster, visu_raster):
    """HEOS-basierter, stoffabhängiger Rastervorschlag (ohne Berechnung der Fläche).

    Parameter und Rückgabe:
      stoff: CoolProp-Reinstoffname
      raster / visu_raster: bestehende Konfigurations-Dictionaries
      -> (neues_raster, neues_visu_raster, metadaten)

    Bestehende Dictionaries werden NICHT mutiert. Einheit: K, Pa, m³/mol.
    """
    Tc, pc, vc = critical_point(stoff)
    Ttr, T_backend_min, T_backend_max = fluid_limits(stoff)


    if not (np.isfinite([Tc, pc, vc, Ttr, T_backend_min, T_backend_max]).all()
            and min(Tc, pc, vc, Ttr) > 0 and Ttr < Tc):
        raise ValueError(f"HEOS liefert ungültige kritische/triple Stoffdaten für {stoff}.")

    T_min = max(0.75 * Tc, 1.01 * Ttr, 1.001 * T_backend_min)
    T_max = min(1.25 * Tc, 0.995 * T_backend_max)

    if not T_min < Tc < T_max:
        raise ValueError(f"Für {stoff} kann kein Bereich beidseits der kritischen Temperatur gewählt werden.")

    sat_start = T_min
    sat_stop = min(0.998 * Tc, T_max)
    
    if sat_start >= sat_stop:
        raise ValueError(f"Kein auswertbarer subkritischer Sättigungsbereich für {stoff}.")

    # Sättigungsvolumina über den gesamten unterkritischen Temperaturbereich
    # ableiten, damit die Oberflächen-/Darstellungsgrenzen für jeden Stoff
    # konsistent mit der HEOS-Sättigungslinie sind.
    sat_grid = np.linspace(sat_start, sat_stop, num=64, endpoint=True)
    sat_volumes = []
    sat_pressures = []
    for T in sat_grid:
        sat = calc_saturation_point(T, stoff)
        if sat is None:
            continue

        liquid, vapor = sat
        v_l = float(liquid[1])
        v_g = float(vapor[1])

        if np.isfinite([v_l, v_g]).all() and 0 < v_l < v_g:
            sat_volumes.append((v_l, v_g))
            if np.isfinite(liquid[2]) and liquid[2] > 0:
                sat_pressures.append(float(liquid[2]))

    if sat_volumes:
        v_l_min = min(v_l for v_l, _ in sat_volumes)
        v_g_max = max(v_g for _, v_g in sat_volumes)

    else:
        LOGGER.warning("Keine HEOS-Sättigungsvolumina für %s; kritisches Volumen als Fallback", stoff)
        v_l_min, v_g_max = 0.25 * vc, 30.0 * vc

    v_min = min(max(0.2 * vc, 0.75 * v_l_min), 0.95 * vc)
    v_max = max(30.0 * vc, 1.2 * v_g_max)

    if not (0 < v_min < vc < v_max and np.isfinite(v_max)):
        raise ValueError(f"Ungültiger Volumenbereich für {stoff}.")

    new_raster = deepcopy(raster)
    new_visu_raster = deepcopy(visu_raster)

    for params in (new_raster["isotherm_grid"], new_visu_raster["isotherm_visu"]):
        params.update(T_start_isotherm=T_min, T_stop_isotherm=T_max, v_start_isotherm=v_min, v_stop_isotherm=v_max,)

    for params in (new_raster["isochore_grid"], new_visu_raster["isochore_visu"]):
        params.update(T_start_isochore=T_min, T_stop_isochore=T_max, v_start_isochore=v_min, v_stop_isochore=v_max,)

    for params in (new_raster["isobare_grid"], new_visu_raster["isobare_visu"]):
        params.update(T_start_isobare=T_min, T_stop_isobare=T_max, p_start_isobare=0.1 * pc, p_stop_isobare=2.0 * pc,)


    for params in (new_raster["binodale_grid"], new_visu_raster["binodale_visu"]):
        params.update(T_start_sat=sat_start, T_stop_sat=sat_stop)
        params["T_step_sat"] = min(params["T_step_sat"], (sat_stop - sat_start) / 2.0)

    for params in (new_raster["spinodale_grid"], new_visu_raster["spinodale_visu"]):
        params.update(T_start_spinodal=sat_start, T_stop_spinodal=sat_stop)

    # Zweiphasenfläche: von der tiefsten Sättigungsdruck-Stützstelle bis pc.
    p_two_phase_min = min(sat_pressures) if sat_pressures else 1e-3 * pc
    new_visu_raster.setdefault("two_phase_visu", {}).update(
        p_start_twophase=p_two_phase_min, p_stop_twophase=pc)

    volume_grids = (
        (new_raster["isotherm_grid"], "T_start_isotherm", "T_stop_isotherm","T_step_isotherm", "v_start_isotherm", "v_stop_isotherm", "v_step_isotherm"),
        (new_visu_raster["isotherm_visu"], "T_start_isotherm", "T_stop_isotherm","T_step_isotherm", "v_start_isotherm", "v_stop_isotherm", "v_step_isotherm"),
        (new_raster["isochore_grid"], "T_start_isochore", "T_stop_isochore", "T_step_isochore", "v_start_isochore", "v_stop_isochore", "v_step_isochore"),
        (new_visu_raster["isochore_visu"], "T_start_isochore", "T_stop_isochore", "T_step_isochore", "v_start_isochore", "v_stop_isochore", "v_step_isochore"),)
    
    for params, T_start, T_stop, T_step, v_start, v_stop, v_step in volume_grids:
        n_T = (int(np.ceil((params[T_stop] - params[T_start]) / params[T_step]))
               + int(params is new_raster["isotherm_grid"]))
        max_v_points = MAX_POINTS // n_T

        if max_v_points < 2:
            raise ValueError("Automatisches Temperaturraster überschreitet ""die zulässige Punktzahl.")
        
        if params.get("v_spacing") == "log":
            params["v_points"] = min(int(params["v_points"]), max_v_points)

        else:
            min_v_step = (params[v_stop] - params[v_start]) / max_v_points
            params[v_step] = max(params[v_step], min_v_step * (1.0 + 1e-12))

    validate_rasters(new_raster, new_visu_raster)

    meta = {
        "Tc": Tc,
        "pc": pc,
        "vc": vc,
        "Ttr": Ttr,
        "T_min": T_min,
        "T_max": T_max,
        "v_min": v_min,
        "v_max": v_max,}

    return new_raster, new_visu_raster, meta


def build_pvt_temperature_grid(cfg, sat_cfg, Tc, Ttr, T_backend_min, T_backend_max):
    """Gemeinsames p-v-T-Temperaturraster mit vollstaendigen Saettigungszeilen.

    Jede unterkritische Oberflaechenzeile erscheint exakt einmal in T_sat.
    Der kritische Punkt ist fuer die Verbindung aller Teilnetze erforderlich.
    Rueckgabe: (T_all, T_sat), jeweils 1D-Arrays in Kelvin.
    """
    T_surface = np.arange(cfg["T_start_isotherm"], cfg["T_stop_isotherm"],
                          cfg["T_step_isotherm"], dtype=np.float64)
    if len(T_surface) and not np.isclose(
            T_surface[-1], cfg["T_stop_isotherm"], rtol=1e-12, atol=1e-12):
        T_surface = np.append(T_surface, cfg["T_stop_isotherm"])
    T_visu = np.arange(sat_cfg["T_start_sat"], sat_cfg["T_stop_sat"],
                       sat_cfg["T_step_sat"], dtype=np.float64)

    if (len(T_surface) < 2 or not np.isfinite(T_surface).all()
            or not np.isfinite(T_visu).all()
            or not np.isfinite([Tc, Ttr, T_backend_min, T_backend_max]).all()):
        raise ValueError("AP2: Ungueltiges p-v-T-Temperaturraster.")
    if not (T_surface[0] < Tc < cfg["T_stop_isotherm"]):
        raise ValueError("AP2: Das p-v-T-Oberflaechenraster muss die kritische "
                         "Temperatur von unten bis oberhalb Tc abdecken.")
    if T_surface[0] < max(Ttr, T_backend_min):
        raise ValueError("AP2: Das p-v-T-Raster beginnt unterhalb des Tripel- "
                         "bzw. HEOS-Temperaturminimums.")

    T_visu = T_visu[(T_visu >= max(Ttr, T_backend_min)) & (T_visu < Tc)]
    candidates = np.unique(np.r_[T_surface, T_visu, Tc])

    # Numerische Fast-Duplikate koennen sich bei inkommensurablen Raster-
    # schritten ergeben. Tc besitzt gegenueber solchen Werten Vorrang.
    tol = max(1e-10, 1e-10 * abs(Tc))
    candidates = candidates[np.abs(candidates - Tc) > tol]
    candidates = np.unique(np.r_[candidates, Tc])
    keep = np.r_[True, np.diff(candidates) > tol]
    T_all = candidates[keep]
    # Auch nach der Entdopplung muss der kritische Knoten exakt erhalten bleiben.
    if not np.any(T_all == Tc):
        T_all = np.unique(np.r_[T_all, Tc])
        near = (np.abs(T_all - Tc) <= tol) & (T_all != Tc)
        T_all = T_all[~near]

    if (len(T_all) < 3 or T_all[0] < T_backend_min
            or T_all[-1] > T_backend_max):
        raise ValueError("AP2: Gemeinsames Temperaturraster ausserhalb HEOS.")
    T_sat = T_all[T_all < Tc]
    if len(T_sat) < 2 or not np.any(T_all > Tc):
        raise ValueError("AP2: Mindestens zwei unterkritische und eine "
                         "ueberkritische Temperaturzeile erforderlich.")
    return T_all, T_sat


def _volume_samples(v_start, v_stop, v_step, v_spacing="linear", v_points=None):
    """Lineares Bestandsraster oder logarithmisches Raster mit fixer Punktzahl.

    v_step bleibt bei logarithmischer Abtastung als manuell nutzbarer
    linearer Parameter gespeichert, wird dort aber nicht ausgewertet.
    """
    if v_spacing not in ("linear", "log"):
        raise ValueError(f"v_spacing muss 'linear' oder 'log' sein, nicht {v_spacing!r}.")

    if v_spacing == "linear":
        if not np.isfinite(v_step) or v_step <= 0:
            raise ValueError("Lineares Volumenraster: v_step muss positiv sein.")
        if not 0 < v_start < v_stop:
            raise ValueError("Ungültiger linearer Volumenbereich.")
        
        return np.arange(v_start, v_stop, v_step)

    if (v_points is None or not float(v_points).is_integer() or not 2 <= v_points <= MAX_POINTS):
        raise ValueError(f"Logarithmisches Volumenraster: " f"v_points muss zwischen 2 und " f"{MAX_POINTS} liegen.")

    if not 0 < v_start < v_stop:
        raise ValueError("Ungültiger logarithmischer Volumenbereich.")

    return np.geomspace(v_start, v_stop, num=int(v_points))


# ============================================================
# 2. Einzelner p-v-T-Zustand (HEOS-Gegenstück zu p_eos(...))
# ============================================================

def p_eos(T, v, stoff):
    """Berechnet den Druck p zu gegebener Temperatur T und molarem Volumen v.

    HEOS bekommt die molare Dichte (rho_molar = 1/v) und die
    Temperatur (DmolarT_INPUTS) und berechnet den Druck.

    Ungültige/instabile Zustände (z.B. numerisch exakt auf der
    Sättigungslinie) werden abgefangen und als np.nan zurückgegeben,
    damit eine Punktreihe nicht an einem einzelnen Fehler abbricht.
    """
    if v <= 0:
        return np.nan
    
    if (not np.isfinite(T) or not np.isfinite(v) or T <= 0 or v <= 0):
        return np.nan

    state = _new_state(stoff)
    try:
        state.update(CP.DmolarT_INPUTS, 1.0 / v, T)
        return state.p()


    except (ValueError, RuntimeError):
        return np.nan


# ============================================================
# 3. Isothermen: T = konst., v variabel -> p
# ============================================================

def calc_isotherme(T_const, stoff, v_start_isotherm, v_stop_isotherm, v_step_isotherm, v_spacing="linear", v_points=None):
    """Isothermenberechnung mit dem HEOS-Backend.

    Für jeden Volumenpunkt v wird die molare Dichte rho_molar = 1/v
    gebildet und zusammen mit T_const an HEOS übergeben (DmolarT_INPUTS).
    HEOS liefert direkt den dazugehörigen Druck p - inklusive des
    physikalisch korrekten horizontalen Verlaufs im Zweiphasengebiet
    (p = p_sat(T) = konstant).
    """
    state = _new_state(stoff)
    isotherm_array = []

    for v in _volume_samples(v_start_isotherm, v_stop_isotherm, v_step_isotherm, v_spacing, v_points):
        if v <= 0:
            continue

        try:
            state.update(CP.DmolarT_INPUTS, 1.0 / v, T_const)
            p = state.p()
        except (ValueError, RuntimeError):
            # z.B. numerisch instabiler Punkt nahe der Phasengrenze
            p = np.nan

        isotherm_array.append((T_const, v, p))

    return isotherm_array


def _explicit_levels(values, name):
    """Explizite Levelarrays (TS-VP03): endlich, positiv, streng steigend."""
    levels = np.asarray(values, dtype=float)
    if levels.ndim != 1 or not np.isfinite(levels).all() or np.any(levels <= 0) or np.any(np.diff(levels) <= 0):
        raise ValueError(f"{name} muss ein endliches, positives, streng steigendes 1D-Array sein.")
    return levels


def multi_isotherm(T_start_isotherm, T_stop_isotherm, T_step_isotherm, v_start_isotherm, v_stop_isotherm, v_step_isotherm, stoff, v_spacing="linear", v_points=None, T_values=None):

    if T_values is None:
        temperatures = np.arange(T_start_isotherm, T_stop_isotherm, T_step_isotherm)

    else:
        temperatures = np.asarray(T_values, dtype=float)

    if (temperatures.ndim != 1 or not np.isfinite(temperatures).all() or np.any(np.diff(temperatures) <= 0)):
        raise ValueError("Ungültiges Temperaturraster.")

    volumes = _volume_samples(v_start_isotherm, v_stop_isotherm, v_step_isotherm, v_spacing, v_points)

    if len(temperatures) * len(volumes) > MAX_POINTS:
        raise ValueError("Erweitertes Oberflächenraster überschreitet " f"{MAX_POINTS:,} Punkte.")

    isotherm_array = []

    for T_const in temperatures:
        kurve = calc_isotherme(float(T_const), stoff, v_start_isotherm, v_stop_isotherm, v_step_isotherm, v_spacing, v_points)
        isotherm_array.append(kurve)

    return np.asarray(isotherm_array, dtype=float)



# ============================================================
# 4. Isochoren: v = konst., T variabel -> p
# ============================================================

def calc_isochore(v_const, stoff, T_start_isochore, T_stop_isochore, T_step_isochore):
    """Isochorenberechnung mit dem HEOS-Backend.

    Die molare Dichte rho_molar = 1/v_const bleibt über die gesamte
    Temperaturschleife konstant, nur T ändert sich - ebenfalls über
    DmolarT_INPUTS gelöst.
    """
    if v_const <= 0:
        return []

    state = _new_state(stoff)
    rho_molar = 1.0 / v_const

    isochore_array = []

    for T in np.arange(T_start_isochore, T_stop_isochore, T_step_isochore):
        try:
            state.update(CP.DmolarT_INPUTS, rho_molar, T)
            p = state.p()
        except (ValueError, RuntimeError):
            p = np.nan

        isochore_array.append((T, v_const, p))

    return isochore_array


def multi_isochore(v_start_isochore, v_stop_isochore, v_step_isochore, T_start_isochore, T_stop_isochore, T_step_isochore, stoff, v_spacing="linear", v_points=None, v_values=None):
    """Isochoren; `v_values` (TS-VP03) ersetzt das Volumenraster durch explizite konstante Volumina [m3/mol]."""
    isochore_array = []

    if v_values is None:
        volumes = _volume_samples(v_start_isochore, v_stop_isochore, v_step_isochore, v_spacing, v_points)
    else:
        volumes = _explicit_levels(v_values, "v_values")

    for v_const in volumes:
        isochore_array.append(calc_isochore(v_const, stoff, T_start_isochore, T_stop_isochore, T_step_isochore))

    return np.asarray(isochore_array, dtype=float)


# ============================================================
# 5. Isobaren: p = konst., T variabel -> v
# ============================================================

def calc_isobare(p_const, stoff, T_start_isobare, T_stop_isobare, T_step_isobare):
    """Isobarenberechnung mit dem HEOS-Backend.

    Hier werden T und p vorgegeben (PT_INPUTS) und HEOS bestimmt daraus
    die molare Dichte bzw. das molare Volumen. Intern löst HEOS dafür
    iterativ seine Helmholtz-Zustandsgleichung nach der Dichte auf.

    Achtung: Im Nassdampfgebiet sind T und p nicht unabhängig
    (p = p_sat(T)), d.h. PT_INPUTS ist dort nicht eindeutig lösbar.
    HEOS wirft in diesem Fall einen ValueError, der hier abgefangen
    und als np.nan markiert wird.
    """
    state = _new_state(stoff)
    isobar_array = []

    for T in np.arange(T_start_isobare, T_stop_isobare, T_step_isobare):
        try:
            state.update(CP.PT_INPUTS, p_const, T)
            v = 1.0 / state.rhomolar()
        except (ValueError, RuntimeError):
            # z.B. T,p exakt (oder numerisch zu nah) auf der Sättigungslinie
            v = np.nan

        isobar_array.append((T, v, p_const))

    # Zweiphasenabschnitt bei unterkritischem Druck ergänzen
    if 0 < p_const < state.p_critical():

        try:
            state.update(CP.PQ_INPUTS, p_const, 0.0)
            T_sat = state.T()
            v_l = 1.0 / state.rhomolar()

            state.update(CP.PQ_INPUTS, p_const, 1.0)
            v_g = 1.0 / state.rhomolar()

            if T_start_isobare <= T_sat < T_stop_isobare:
                for Q in np.linspace(0.0, 1.0, 30):
                    v = (1.0 - Q) * v_l + Q * v_g
                    isobar_array.append((T_sat, v, p_const))

        except (ValueError, RuntimeError):
            pass

    # Zustandspunkte nach Temperatur und Volumen sortieren
    isobar_array.sort(key=lambda punkt: (punkt[0], punkt[1] if np.isfinite(punkt[1]) else float("inf")))

    return isobar_array

def multi_isobare(p_start_isobare, p_stop_isobare, p_step_isobare, T_start_isobare, T_stop_isobare, T_step_isobare, stoff, p_values=None):
    """Isobaren; `p_values` (TS-VP03) ersetzt das Druckraster durch explizite konstante Druecke [Pa]."""
    isobar_array = []

    if p_values is None:
        pressures = np.arange(p_start_isobare, p_stop_isobare, p_step_isobare)
    else:
        pressures = _explicit_levels(p_values, "p_values")

    for p_const in pressures:
        isobar_array.append(calc_isobare(p_const, stoff, T_start_isobare, T_stop_isobare, T_step_isobare))

    return [np.asarray(kurve, dtype=float).reshape(-1, 3) for kurve in isobar_array]


# ============================================================
# 6. Phasengrenze / Zweiphasengebiet
# ============================================================


def calc_saturation_point(T, stoff):
    """
    Berechnet Siede- und Taupunkt bei einer Temperatur T.

    Rückgabe:
        liquid_point = (T, v_liquid, p_sat)
        vapor_point  = (T, v_vapor, p_sat)

    T [K], v [m³/mol], p [Pa]
    """

    state = _new_state(stoff)

    # Gültigen Temperaturbereich bestimmen
    T_crit = state.T_critical()
    T_triple = state.Ttriple()

    if not np.isfinite(T):
        return None

    if not (T_triple <= T < T_crit):
        return None

    try:
        # 1. Gesättigte Flüssigkeit (Q = 0)
        state.update(CP.QT_INPUTS, 0.0, T)
        p_sat = state.p()
        v_liquid = 1.0 / state.rhomolar()

        # 2. Gesättigter Dampf (Q = 1)
        state.update(CP.QT_INPUTS, 1.0, T)
        v_vapor = 1.0 / state.rhomolar()

    except (ValueError, RuntimeError):
        return None

    liquid_point = (T, v_liquid, p_sat)
    vapor_point = (T, v_vapor, p_sat)

    return liquid_point, vapor_point


def multi_binodale(T_start_sat, T_stop_sat, T_step_sat, stoff, T_values=None):
    """
    Berechnet die Binodale über einen Temperaturbereich.

    Rückgabe:
        binodal_liquid = [(T, v_liquid, p), ...]
        binodal_vapor  = [(T, v_vapor, p), ...]

    T [K], v [m³/mol], p [Pa]
    """

    if T_values is None:
        temperatures = np.arange(T_start_sat, T_stop_sat, T_step_sat)
    else:
        temperatures = np.asarray(T_values, dtype=float)

    if temperatures.ndim != 1 or not np.isfinite(temperatures).all():
        raise ValueError("Ungültiges Sättigungstemperaturraster.")

    binodal_liquid = []
    binodal_vapor = []

    for T in temperatures:
        points = calc_saturation_point(float(T), stoff)

        if points is None:
            raise ValueError(f"Sättigungszustand bei T={T:.8f} K " f"für {stoff} nicht berechenbar.")

        liquid_point, vapor_point = points
        binodal_liquid.append(liquid_point)
        binodal_vapor.append(vapor_point)

    return (np.asarray(binodal_liquid, dtype=float).reshape(-1, 3), np.asarray(binodal_vapor, dtype=float).reshape(-1, 3))


def berechne_zweiphasenmaske(stoff, surface, binodal_liquid=None, binodal_vapor=None):

    phase_mask = np.zeros(surface.shape[:2], dtype=bool)

    sat_map = {}

    if binodal_liquid is not None and binodal_vapor is not None:
        for liquid, vapor in zip(binodal_liquid, binodal_vapor):
            sat_map[float(liquid[0])] = (liquid, vapor)

    for i in range(surface.shape[0]):
        T = float(surface[i, 0, 0])

        sat = sat_map.get(T)

        if sat is None:
            sat = calc_saturation_point(T, stoff)

        if sat is None:
            continue

        liquid, vapor = sat

        v_l = liquid[1]
        v_g = vapor[1]
        p_sat = liquid[2]

        v = surface[i, :, 1]
        p = surface[i, :, 2]

        phase_mask[i, :] = (
            (v >= v_l)
            & (v <= v_g)
            & np.isfinite(p)
            & np.isclose(p, p_sat, rtol=1e-3)
        )

    return phase_mask


# ============================================================
# 7. Spinodale
# ============================================================


def calc_spinodale(stoff):
    """HEOS-Spinodalstuetzstellen als zwei rohe (T,v,p)-Arrays (SI).

    Die Arrays enthalten KEINEN synthetisch angehaengten kritischen Punkt.
    Fehlgeschlagene homogene Druckauswertungen verbleiben als NaN, damit
    unvollstaendige Bereiche bei der Validierung erkennbar bleiben.
    """
    spin_state = _new_state(stoff)
    spin_state.build_spinodal()
    data = spin_state.get_spinodal_data()
    tau = np.asarray(data.tau, dtype=float)
    delta = np.asarray(data.delta, dtype=float)
    if tau.ndim != 1 or delta.shape != tau.shape:
        raise ValueError("HEOS-Spinodaldaten: tau und delta inkompatibel.")

    T_reducing = float(spin_state.T_reducing())
    rho_reducing = float(spin_state.rhomolar_reducing())
    eos_state = _new_state(stoff)
    Tc = float(eos_state.T_critical())
    rho_c = float(eos_state.rhomolar_critical())
    if not np.isfinite([Tc, rho_c, T_reducing, rho_reducing]).all() or min(
            Tc, rho_c, T_reducing, rho_reducing) <= 0:
        raise ValueError("HEOS-Spinodale: ungueltige Reduktions- oder kritische Daten.")

    branches = {"liquid": [], "vapor": []}
    for tau_i, delta_i in zip(tau, delta):
        if not np.isfinite([tau_i, delta_i]).all() or min(tau_i, delta_i) <= 0:
            continue
        T = T_reducing / tau_i
        rho_molar = delta_i * rho_reducing
        if not np.isfinite([T, rho_molar]).all() or T <= 0 or rho_molar <= 0:
            continue
        # Der kritische Punkt ist kein regulaerer Unterkritisch-Stuetzpunkt.
        if T >= Tc or np.isclose(rho_molar, rho_c, rtol=1e-9, atol=0.0):
            continue
        liquid = rho_molar > rho_c
        phase = CP.iphase_liquid if liquid else CP.iphase_gas
        pressure = float("nan")
        try:
            eos_state.specify_phase(phase)
            eos_state.update(CP.DmolarT_INPUTS, float(rho_molar), float(T))
            pressure = float(eos_state.p())
        except (ValueError, RuntimeError) as exc:
            LOGGER.debug("Spinodaldruck bei T=%g K nicht auswertbar: %s", T, exc)
        finally:
            eos_state.unspecify_phase()
        branches["liquid" if liquid else "vapor"].append((T, 1.0/rho_molar, pressure))

    def sorted_branch(points):
        if not points:
            return np.empty((0, 3), dtype=float)
        out = np.asarray(points, dtype=float)
        out = out[np.argsort(out[:, 0], kind="stable")]
        # Bei doppelten T-Stuetzstellen keinen kuenstlichen Vertikalsprung erzeugen.
        return out[np.r_[True, np.diff(out[:, 0]) > 1e-10]]

    liquid = sorted_branch(branches["liquid"])
    vapor = sorted_branch(branches["vapor"])
    LOGGER.info("HEOS-Spinodale %s: Rohpunkte liquid=%d, vapor=%d", stoff,
                len(liquid), len(vapor))
    return liquid, vapor


def filter_spinodale(kurve, T_min, T_max):
    """Temperaturfilter [T_min,T_max] ohne erneutes Ergaenzen von Endpunkten."""
    if not np.isfinite([T_min, T_max]).all() or not 0 < T_min < T_max:
        raise ValueError("Spinodalenfilter: 0 < T_min < T_max erforderlich.")
    kurve = np.asarray(kurve, dtype=float)
    if kurve.size == 0:
        return np.empty((0, 3), dtype=float)
    if kurve.ndim != 2 or kurve.shape[1] != 3:
        raise ValueError("Spinodalast muss die Form (n,3) haben.")
    maske = np.isfinite(kurve[:, 0]) & (kurve[:, 0] >= T_min) & (kurve[:, 0] <= T_max)
    return kurve[maske]


def prepare_spinodal_branch(kurve, critical, T_min, T_max, branch):
    """AP3: Filtern, Endpunkt begruendet ergaenzen und Astqualitaet bewerten.

    Rueckgabe: (T-v-p-Array, Diagnosedict). Der kritische Punkt wird nur bei
    expliziter Bereichsabdeckung und nachgewiesener kritischer Annaeherung
    hinzugefuegt; andernfalls wird der Ast als unvollstaendig gekennzeichnet.
    """
    if branch not in ("liquid", "vapor"):
        raise ValueError("Spinodalast muss 'liquid' oder 'vapor' sein.")
    Tc, pc, vc = map(float, critical)  # Kritisches Tupel: (Tc,pc,vc).
    if not np.isfinite([Tc, pc, vc]).all() or min(Tc, pc, vc) <= 0:
        raise ValueError("Ungueltige kritische Stoffdaten fuer AP3.")
    points = filter_spinodale(kurve, T_min, T_max)
    if len(points):
        points = points[np.argsort(points[:, 0], kind="stable")]
    # Rohdrucke <= 0 koennen zur homogenen EOS gehoeren, sind aber auf log(p)
    # nicht darstellbar; diese Daten NICHT heimlich loeschen.
    valid = (np.isfinite(points).all(axis=1) & (points[:, 1] > 0)
             & (points[:, 0] > 0)) if len(points) else np.zeros(0, dtype=bool)
    right_phase = (points[:, 1] < vc if branch == "liquid" else
                   points[:, 1] > vc) if len(points) else np.zeros(0, dtype=bool)
    n_bad_phase = int(np.count_nonzero(valid & ~right_phase))
    n_invalid = int(np.count_nonzero(~valid))
    n_nonpositive_p = int(np.count_nonzero(valid & (points[:, 2] <= 0)))
    # Schwellwert ist explizit eine konservative geometrische Heuristik,
    # kein Beweis fuer die Konvergenz gegen den kritischen Punkt.
    last = points[-1] if len(points) else None
    near_critical = bool(
        last is not None and valid[-1] and right_phase[-1]
        and 0 < Tc - last[0] <= 0.01 * Tc
        and abs(last[1] / vc - 1.0) <= 0.5
    )
    # Der Benutzerbereich darf durch das Anhaengen niemals erweitert werden.
    add_critical = bool(len(points) >= 2 and T_min <= Tc <= T_max
                        and near_critical and points[-1, 0] < Tc)
    if add_critical:
        points = np.vstack((points, np.array([[Tc, vc, pc]])))
    sub = points[points[:, 0] < Tc]
    if len(sub) >= 2:
        spacings = np.diff(sub[:, 0])
        tolerance_T = min(0.01 * Tc, max(0.5, 1.5 * float(np.median(spacings))))
    else:
        tolerance_T = 0.5
    covers_start = bool(len(sub) >= 2 and sub[0, 0] <= T_min + tolerance_T)
    covers_end = bool(add_critical if T_max >= Tc else
                      len(sub) >= 2 and sub[-1, 0] >= T_max - tolerance_T)
    temperature_order_ok = bool(not len(points) or np.all(np.diff(points[:, 0]) > 0))
    # Grobe Luecken im temperaturparametrisierten Ast als Indiz fuer
    # fehlende Spinodalabschnitte klassifizieren (nicht als EOS-Nachweis).
    gaps = np.diff(sub[:, 0]) if len(sub) >= 2 else np.empty(0)
    large_gaps = bool(len(gaps) >= 2 and np.max(gaps) > max(
        10.0 * float(np.median(gaps)), 0.10 * (T_max - T_min)))
    # Volle Sichtbarkeit erfordert >0 Pa; alle numerisch unvollstaendigen
    # oder phasenfalschen Punkte setzen den Status 'unvollstaendig'.
    n_unplottable = n_invalid + n_nonpositive_p
    geometry_ok = bool(len(sub) >= 3 and covers_start and covers_end
                       and n_bad_phase == 0 and n_unplottable == 0
                       and temperature_order_ok and not large_gaps)
    report = {
        "branch": branch, "points": int(len(points)),
        "subcritical_points": int(len(sub)),
        "invalid_states": n_invalid, "nonpositive_pressure": n_nonpositive_p,
        "wrong_phase": n_bad_phase, "covers_start": covers_start,
        "covers_end": covers_end, "critical_appended": add_critical,
        "temperature_order_ok": temperature_order_ok,
        "large_gaps": large_gaps,
        "geometry_complete": geometry_ok,
    }
    return points, report


# ============================================================
# 8. Gleichgewichts-Zweiphasenfläche
# ============================================================

def clip_two_phase_surface(surface, p_min=None, p_max=None, pc=None, rel_tol=1e-6):
    """Beschneidet die Zweiphasenfläche auf T-Zeilen mit p_min <= p_sat <= p_max.

    Rückgabe: (Fläche, erreicht_kritischen_Punkt). Weniger als zwei verbleibende
    Zeilen ergeben eine leere Fläche (0, 0, 3). Der kritische Abschluss bleibt nur
    erhalten, wenn die obere Grenze pc erreicht und keine Zeile abgeschnitten wird.
    """
    surface = np.asarray(surface, dtype=float)
    if surface.ndim != 3 or surface.shape[0] == 0:
        return surface, False

    p_row = surface[:, :, 2].mean(axis=1)
    keep = np.ones(len(p_row), dtype=bool)
    if p_min is not None:
        keep &= p_row >= p_min * (1.0 - rel_tol)
    if p_max is not None:
        keep &= p_row <= p_max * (1.0 + rel_tol)

    if np.count_nonzero(keep) < 2:
        return np.empty((0, 0, 3)), False

    reaches_critical = bool(keep[-1] and (p_max is None or pc is None or p_max >= pc * (1.0 - rel_tol)))
    return surface[keep], reaches_critical


def calc_two_phase_surface(binodal_liquid, binodal_vapor, n_quality=40):
    """
    Berechnet die p-v-T-Zweiphasenfläche aus Binodaldaten.

    Eingabe:
        binodal_liquid: (n_T, 3)
        binodal_vapor:  (n_T, 3)
        n_quality: Anzahl der Dampfgehaltswerte

    Rückgabe:
        two_phase_surface: (n_T, n_Q, 3)
        Q_values:         (n_Q,)

    Reihenfolge: (T [K], v [m³/mol], p [Pa])
    """

    liquid = np.asarray(binodal_liquid, dtype=float)
    vapor = np.asarray(binodal_vapor, dtype=float)

    # 1. Arraystruktur überprüfen
    if (liquid.ndim != 2 or liquid.shape[1] != 3 or liquid.shape != vapor.shape or len(liquid) < 2):
        raise ValueError("Ungültige Binodalarrays: Erwartet (n_T, 3).")

    if not isinstance(n_quality, (int, np.integer)) or n_quality < 2:
        raise ValueError("n_quality muss ganzzahlig und >= 2 sein.")

    # 2. Thermodynamische Plausibilitätsprüfung
    if not (np.isfinite(liquid).all() and np.isfinite(vapor).all()):
        raise ValueError("Binodale enthält ungültige Werte.")

    if not np.allclose(liquid[:, 0], vapor[:, 0], rtol=0, atol=1e-8):
        raise ValueError("Temperaturraster stimmen nicht überein.")

    if not np.allclose(liquid[:, 2], vapor[:, 2], rtol=1e-8, atol=1e-4):
        raise ValueError("Sättigungsdrücke stimmen nicht überein.")

    if np.any(liquid[:, 1] >= vapor[:, 1]):
        raise ValueError("Flüssigkeitsvolumen muss kleiner sein.")

    if np.any(liquid[:, 2] <= 0):
        raise ValueError("Sättigungsdruck muss positiv sein.")

    n_T = len(liquid)

    if n_T * n_quality > MAX_POINTS:
        raise ValueError("Zweiphasenfläche würde " f"{n_T * n_quality:,} Punkte erzeugen. " f"Maximum: {MAX_POINTS:,}.")

    # 3. Dampfgehaltsraster erstellen
    Q_values = np.linspace(0.0, 1.0, n_quality)

    n_T = len(liquid)

    # 4. Temperatur und Druck auf das Q-Raster erweitern
    T_grid = np.broadcast_to(liquid[:, 0, None], (n_T, n_quality))
    p_grid = np.broadcast_to( liquid[:, 2, None], (n_T, n_quality))

    # 5. Molares Mischungsvolumen berechnen
    v_grid = ( (1.0 - Q_values[None, :]) * liquid[:, 1, None] + Q_values[None, :] * vapor[:, 1, None])

    # AP2: exakt dieselben Randknoten wie die Binodalen, ohne Rundungsdrift.
    v_grid[:, 0] = liquid[:, 1]
    v_grid[:, -1] = vapor[:, 1]

    # 6. Gemeinsames 3D-Datenarray
    two_phase_surface = np.stack((T_grid, v_grid, p_grid), axis=-1)

    return two_phase_surface, Q_values


def validate_pvt_phase_boundaries(daten):
    """AP2: Exakte Temperatur-, Binodal- und Zweiphasenrandkontinuitaet.

    Erwartet das bestehende p-v-T-Ergebnis-Dictionary nach Anhaengen des
    kritischen Punkts an beide Binodalzweige. Keine GUI/VTK-Abhaengigkeit.
    """
    surface = np.asarray(daten["surface"])
    liquid = np.asarray(daten["binodal_liquid"])
    vapor = np.asarray(daten["binodal_vapor"])
    mixture = np.asarray(daten["two_phase_surface"])
    q = np.asarray(daten["Q_values"])
    mask = np.asarray(daten["phase_mask"])
    Tc, pc, vc = daten["critical_point"]

    if (surface.ndim != 3 or surface.shape[-1] != 3
            or min(surface.shape[:2]) < 2
            or liquid.ndim != 2 or liquid.shape[1:] != (3,)
            or vapor.shape != liquid.shape or len(liquid) < 3
            or mixture.ndim != 3 or mixture.shape[-1] != 3
            or q.ndim != 1 or len(q) < 2
            or mixture.shape != (len(liquid) - 1, len(q), 3)
            or mask.shape != surface.shape[:2] or mask.dtype != np.bool_):
        raise ValueError("AP2: Inkompatible Oberflaechen- oder Phasenarrays.")
    if (not np.isfinite(liquid).all() or not np.isfinite(vapor).all()
            or not np.isfinite(mixture).all()
            or not np.isfinite([Tc, pc, vc]).all()
            or min(Tc, pc, vc) <= 0):
        raise ValueError("AP2: Nichtendliche oder ungueltige Phasengrenzen.")

    T_all = surface[:, 0, 0]
    if (not np.isfinite(T_all).all() or not np.all(np.diff(T_all) > 0)
            or not np.all(surface[..., 0] == T_all[:, None])
            or not np.any(T_all == Tc)
            or not np.any(T_all > Tc)):
        raise ValueError("AP2: Ungueltige Temperaturachse oder Tc-Zeile.")

    sub = T_all[T_all < Tc]
    if (not np.array_equal(sub, liquid[:-1, 0])
            or not np.array_equal(sub, vapor[:-1, 0])
            or not np.array_equal(liquid[:-1, 2], vapor[:-1, 2])):
        raise ValueError("AP2: Fehlende Saettigungszeilen oder verschiedene psat.")
    if (np.any(liquid[:-1, 1] >= vapor[:-1, 1])
            or np.any(liquid[:-1, 1] <= 0)
            or np.any(liquid[:-1, 2] <= 0)):
        raise ValueError("AP2: Saettigungsvolumina oder Druck unplausibel.")
    critical_tvp = np.asarray((Tc, vc, pc), dtype=float)
    if (not np.array_equal(liquid[-1], critical_tvp)
            or not np.array_equal(vapor[-1], critical_tvp)):
        raise ValueError("AP2: Binodalen enden nicht exakt bei (Tc,vc,pc).")
    if (not np.array_equal(q[[0, -1]], [0., 1.])
            or not np.array_equal(mixture[:, 0, :], liquid[:-1])
            or not np.array_equal(mixture[:, -1, :], vapor[:-1])):
        raise ValueError("AP2: Zweiphasenraender stimmen nicht mit Binodalen ueberein.")

    # Die randkonforme Netzbildung benoetigt beiderseits der Binodalen
    # wenigstens einen regulaeren Volumenstuetzpunkt.
    if (np.any(surface[T_all < Tc, 0, 1] >= liquid[:-1, 1])
            or np.any(surface[T_all < Tc, -1, 1] <= vapor[:-1, 1])
            or not np.all(np.diff(surface[..., 1], axis=1) > 0)):
        raise ValueError("AP2: Saettigungsvolumen ausserhalb des Volumenrasters.")

    # Geometrische HEOS-Daten sind nach diesem Test fuer die Renderer eindeutig.
    return None
