"""
Datei: isoline_levels.py
Projekt: ThermoSurface
Kategorie: Thermodynamik / Isolinien

Beschreibung:
Berechnet und validiert die Niveaus der Isothermen, Isochoren und Isobaren.

Aufgaben:
- Erzeugt Isolinienwerte über eine geometrische Progression.
- Normalisiert Konfigurationen und unterstützt die Kurvenabtastung.

Abhängigkeiten:
- Projektmodule (direkt): keine
- Bibliotheken: numpy

Verwendet von:
- Projektmodule: coolprop_heos_fvt.py, gui.py, main.py
"""

import math
from typing import Optional

import numpy as np

FACTOR_DEFAULT = 3.14
MAX_LEVELS_DEFAULT = 200      # empfohlene technische Obergrenze je Familie und Diagramm (TS-VP01, Abschn. 1.2)
MAX_LEVELS_CEILING = 10000    # harte Obergrenze fuer eine benutzerdefinierte max_levels-Vorgabe
FACTOR_MIN_STEP = 1e-9        # factor muss mindestens 1 + FACTOR_MIN_STEP betragen (sonst kein Fortschritt)
BOUNDARY_RTOL = 1e-12         # Grenztoleranz fuer x_n <= stop

FAMILIES = ("isothermen", "isochoren", "isobaren")
FAMILY_QUANTITY = {"isothermen": ("T", "K"), "isochoren": ("v", "m3/mol"), "isobaren": ("p", "Pa")}
_CONFIG_KEYS = ("start", "stop", "factor", "max_levels")


class IsolineLevelError(ValueError):
    """Kontrollierte Validierungsmeldung fuer ungueltige Isolinienvorgaben (keine Hänger, keine Ueberlaeufe)."""


def _finite_float(value, name):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise IsolineLevelError(f"{name} muss eine Zahl sein (erhalten: {value!r}).") from None
    if not math.isfinite(number):
        raise IsolineLevelError(f"{name} muss endlich sein (erhalten: {value!r}).")
    return number


def _is_number(value):
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


def validate_levels_input(start, stop, factor, max_levels=MAX_LEVELS_DEFAULT):
    """Validiert die vier Eingaben und gibt (start, stop, factor, max_levels) als float/float/float/int zurueck."""
    start = _finite_float(start, "start")
    stop = _finite_float(stop, "stop")
    factor = _finite_float(factor, "factor")
    if isinstance(max_levels, bool) or not _is_number(max_levels) or not float(max_levels).is_integer():
        raise IsolineLevelError(f"max_levels muss eine ganze Zahl sein (erhalten: {max_levels!r}).")
    max_levels = int(max_levels)
    if not 1 <= max_levels <= MAX_LEVELS_CEILING:
        raise IsolineLevelError(f"max_levels muss zwischen 1 und {MAX_LEVELS_CEILING} liegen (erhalten: {max_levels}).")
    if start <= 0.0:
        raise IsolineLevelError(f"start muss > 0 sein (erhalten: {start:g}); die geometrische Progression ist fuer "
                                "start <= 0 nicht definiert.")
    if stop < start:
        raise IsolineLevelError(f"stop ({stop:g}) darf nicht kleiner als start ({start:g}) sein.")
    if factor <= 1.0 + FACTOR_MIN_STEP:
        raise IsolineLevelError(f"factor muss > 1 sein (mindestens 1 + {FACTOR_MIN_STEP:g}; erhalten: {factor!r}).")
    return start, stop, factor, max_levels


def generate_levels(start, stop, factor=FACTOR_DEFAULT, max_levels=MAX_LEVELS_DEFAULT):
    """Streng wachsendes Array x_n = start * factor**n mit x_n <= stop (inklusive Grenztoleranz).

    Beispiel start=10, stop=1000, factor=3.14: 10, 31.4, 98.596, 309.59144, 972.1171216 (1000 wird nicht ergaenzt).
    Wirft IsolineLevelError bei ungueltigen Eingaben, Ueberlauf, fehlendem Rundungsfortschritt oder mehr als
    `max_levels` Werten (Schutz gegen factor ~ 1).
    """
    start, stop, factor, max_levels = validate_levels_input(start, stop, factor, max_levels)
    limit = stop * (1.0 + BOUNDARY_RTOL)

    # log1p(factor-1) bleibt auch fuer factor nahe 1 genau; die Schaetzung wird unten exakt nachgeprueft.
    ratio = math.log(stop / start) if stop > start else 0.0
    n_estimate = int(math.floor(ratio / math.log1p(factor - 1.0) + 1e-9))
    if n_estimate > max_levels:
        raise IsolineLevelError(
            f"Mehr als max_levels={max_levels} Isolinien fuer start={start:g}, stop={stop:g}, factor={factor:g}; "
            "Faktor vergroessern, Bereich verkleinern oder max_levels erhoehen.")

    def level(n):
        with np.errstate(over="ignore"):
            return float(start * np.power(factor, float(n)))

    n = max(n_estimate, 0)
    while n > 0 and level(n) > limit:
        n -= 1
    while level(n + 1) <= limit and n + 1 <= max_levels:
        n += 1
    count = n + 1
    if count > max_levels:
        raise IsolineLevelError(
            f"Mehr als max_levels={max_levels} Isolinien fuer start={start:g}, stop={stop:g}, factor={factor:g}.")

    with np.errstate(over="ignore"):
        values = start * np.power(factor, np.arange(count, dtype=float))
    if not np.isfinite(values).all():
        raise IsolineLevelError("Ueberlauf bei der Berechnung der Isolinienwerte.")
    if count > 1 and not np.all(np.diff(values) > 0.0):
        raise IsolineLevelError("Kein numerischer Fortschritt der Isolinienwerte (factor zu nahe an 1).")
    return values


# ---------------------------------------------------------------- Konfiguration je Familie

def normalise_family_config(raw, family):
    """Kanonische Familienkonfiguration {'start','stop','factor','max_levels'}; start/stop None = automatisch."""
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise IsolineLevelError(f"isolines.{family} muss ein Mapping sein.")
    unknown = sorted(set(raw) - set(_CONFIG_KEYS))
    if unknown:
        raise IsolineLevelError(f"isolines.{family}: unbekannte Schluessel {unknown}.")

    def optional_bound(key):
        value = raw.get(key)
        if value is None or (isinstance(value, str) and value.strip().lower() in ("", "auto")):
            return None
        number = _finite_float(value, f"isolines.{family}.{key}")
        # 0 steht in der GUI fuer "automatisch"; echte Werte sind positiv.
        return None if number == 0.0 else number

    factor = _finite_float(raw.get("factor", FACTOR_DEFAULT), f"isolines.{family}.factor")
    max_levels = raw.get("max_levels", MAX_LEVELS_DEFAULT)
    out = {"start": optional_bound("start"), "stop": optional_bound("stop"), "factor": factor,
           "max_levels": int(_finite_float(max_levels, f"isolines.{family}.max_levels"))}
    if out["start"] is not None and out["start"] < 0 or out["stop"] is not None and out["stop"] < 0:
        raise IsolineLevelError(f"isolines.{family}: start/stop muessen positiv sein.")
    return out


def normalise_isoline_config(raw):
    """Kanonische Konfiguration aller drei Familien (fehlende Familien = Standardfaktor, automatischer Bereich)."""
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise IsolineLevelError("isolines muss ein Mapping sein.")
    unknown = sorted(set(raw) - set(FAMILIES))
    if unknown:
        raise IsolineLevelError(f"isolines: unbekannte Familien {unknown}; erlaubt: {list(FAMILIES)}.")
    return {family: normalise_family_config(raw.get(family), family) for family in FAMILIES}


def default_isoline_config():
    return normalise_isoline_config({})


def resolve_levels(config, domains):
    """Levels je Familie aus der Konfiguration; `domains[family] = (lo, hi)` ist Auto-Bereich UND Gueltigkeitsfenster.

    - start/stop = None -> Grenzen der Domaene (die geometrische Gitterstruktur beginnt dann an `lo`).
    - Levels ausserhalb der Domaene (dort liegen keine berechneten Flaechendaten) werden nicht gezeichnet und im
      Bericht als `dropped_outside_domain` ausgewiesen.
    Rueckgabe (levels, report): levels[family] = ndarray (evtl. leer), report[family] = protokollierbares Dict.
    """
    config = normalise_isoline_config(config)
    levels, report = {}, {}
    for family in FAMILIES:
        cfg = config[family]
        quantity, unit = FAMILY_QUANTITY[family]
        domain = domains.get(family)
        entry = {"quantity": quantity, "unit": unit, "factor": cfg["factor"], "max_levels": cfg["max_levels"],
                 "domain": None if domain is None else (float(domain[0]), float(domain[1]))}
        if domain is None:
            levels[family] = np.empty(0)
            entry.update(start=cfg["start"], stop=cfg["stop"], levels=[], n_levels=0, dropped_outside_domain=[],
                         samples_per_curve=[], note="Familie in diesem Diagramm nicht verfuegbar")
            report[family] = entry
            continue
        lo, hi = float(domain[0]), float(domain[1])
        start = cfg["start"] if cfg["start"] is not None else lo
        stop = cfg["stop"] if cfg["stop"] is not None else hi
        try:
            generated = generate_levels(start, stop, cfg["factor"], cfg["max_levels"])
        except IsolineLevelError as exc:
            raise IsolineLevelError(f"{family}: {exc}") from None
        inside = (generated >= lo * (1.0 - 1e-9)) & (generated <= hi * (1.0 + 1e-9))
        levels[family] = generated[inside]
        entry.update(start=float(start), stop=float(stop), levels=[float(x) for x in levels[family]],
                     n_levels=int(inside.sum()), dropped_outside_domain=[float(x) for x in generated[~inside]],
                     samples_per_curve=[])
        report[family] = entry
    return levels, report


def record_sampling(report, family, curves):
    """Protokolliert die Abtastpunktzahl je Kurve getrennt von der Levelanzahl (Zahl der Punkte, nicht der Levels)."""
    report[family]["samples_per_curve"] = [int(len(np.asarray(curve))) for curve in curves]
    report[family]["n_curves"] = len(curves)


def config_signature(config):
    """Kanonische, JSON-faehige Darstellung fuer Cache-Schluessel (nur Berechnungsparameter)."""
    config = normalise_isoline_config(config)
    return {family: {key: cfg[key] for key in _CONFIG_KEYS} for family, cfg in config.items()}


# Automatikbereiche fuer Diagramme ohne eigenes v-/p-Fenster (relativ zum kritischen Zustand).
V_DOMAIN_FACTORS = (0.5, 30.0)
P_DOMAIN_FACTORS = (0.1, 10.0)


def isotherm_over_volume(T, Tc, v_min, v_max, n_v, n_quality, evaluate, saturation):
    """Isotherme (T, v, x) ueber das molare Volumen: Fluessigkeit -> Nassdampf -> Dampf (T < Tc), sonst ein Ast.

    evaluate(T, v, phase) -> (T, v, x) wirft ValueError bei nicht auswertbarem Zustand (Zeile wird NaN, keine
    Bruecke); phase in {"liquid", "gas", None}. saturation(T) -> ((v_l, x_l), (v_g, x_g)) oder None.
    Der Nassdampfabschnitt ist die lineare Mischungsgerade in Q (wie die Zweiphasenflaeche).
    """
    def points(volumes, phase):
        rows = []
        for v in volumes:
            try:
                rows.append(evaluate(T, float(v), phase))
            except ValueError:
                rows.append((T, float(v), np.nan))
        return rows

    sat = saturation(T) if T < Tc else None
    if sat is None or not (v_min < sat[0][0] < sat[1][0] < v_max):
        return np.asarray(points(np.geomspace(v_min, v_max, 2 * n_v - 1), None), dtype=float)

    (v_l, x_l), (v_g, x_g) = sat
    q = np.linspace(0.0, 1.0, n_quality)[1:-1]
    liquid = points(np.geomspace(v_min, v_l, n_v)[:-1], "liquid") + [(T, v_l, x_l)]
    mix = [(T, (1.0 - qq) * v_l + qq * v_g, (1.0 - qq) * x_l + qq * x_g) for qq in q]
    vapor = [(T, v_g, x_g)] + points(np.geomspace(v_g, v_max, n_v)[1:], "gas")
    return np.asarray(liquid + mix + vapor, dtype=float)


def build_isolines(levels, report, makers, log=None):
    """Berechnet die Kurven entlang der konstanten Levels mit `makers[family](level) -> (n,3)`-Array.

    Fehlgeschlagene oder nirgends auswertbare Levels werden nachvollziehbar ausgelassen (Bericht: failed_levels,
    empty_levels); es werden keine Ersatzzustaende erfunden. Rueckgabe: {family: [Kurve, ...]}. Die Abtastpunktzahl
    je Kurve wird getrennt von der Levelanzahl protokolliert.
    """
    curves = {}
    for family in FAMILIES:
        kept, failed, empty, out = [], [], [], []
        for level in levels[family]:
            try:
                curve = np.asarray(makers[family](float(level)), dtype=float)
            except (ValueError, RuntimeError, ArithmeticError) as exc:
                failed.append(float(level))
                if log is not None:
                    log("%s: Level %.9g nicht berechenbar: %s", family, float(level), exc)
                continue
            if curve.ndim != 2 or curve.shape[1] != 3 or np.count_nonzero(np.isfinite(curve).all(axis=1)) < 2:
                empty.append(float(level))
                continue
            kept.append(float(level))
            out.append(curve)
        curves[family] = out
        report[family].update(computed_levels=kept, failed_levels=failed, empty_levels=empty)
        record_sampling(report, family, out)
    return curves
