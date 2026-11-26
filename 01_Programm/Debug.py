"""Zentrale, optionale Diagnosen der p-v-T-Anwendung.

Keine Thermodynamik im Produktivpfad, keine Berechnungen beim Import.
Datenkonvention: (T [K], v [m³/mol], p [Pa]).
Importiert weder main noch gui noch visualization noch coolprop_heos.
"""

import logging
import sys
from time import perf_counter

import numpy as np


LOGGER = logging.getLogger("thermosurface.Debug")
ROOT_LOGGER = logging.getLogger("thermosurface")


def configure_debug(enabled):
    """Einziger zentraler Schalter für das anwendungseigene Logging."""
    ROOT_LOGGER.setLevel(logging.DEBUG if enabled else logging.WARNING)
    ROOT_LOGGER.propagate = False
    if not any(getattr(h, "_thermosurface_handler", False) for h in ROOT_LOGGER.handlers):
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter(
            "%(asctime)s %(levelname)-7s [%(name)s] %(message)s",
            datefmt="%H:%M:%S",
        ))
        handler._thermosurface_handler = True
        ROOT_LOGGER.addHandler(handler)
    LOGGER.debug("Debug-Modus aktiviert")


def diagnose_raster(stoff, raster, visu_raster, meta):
    """Überprüft automatisch ermittelte Raster ohne Debug-Importzyklen."""
    LOGGER.debug("Auto-Raster für %s: Tc=%.8g pc=%.8g vc=%.8g Ttr=%.8g",
                 stoff, meta["Tc"], meta["pc"], meta["vc"], meta["Ttr"])
    LOGGER.debug("Bereich: T=[%.8g, %.8g] K v=[%.8g, %.8g] m³/mol",
                 meta["T_min"], meta["T_max"], meta["v_min"], meta["v_max"])
    if not (meta["T_min"] < meta["Tc"] < meta["T_max"]
            and meta["v_min"] < meta["vc"] < meta["v_max"]):
        raise AssertionError("Kritischer Punkt nicht vollständig innerhalb des Auto-Rasters")
    for groupname, groups in (("Berechnung", raster), ("Darstellung", visu_raster)):
        for name, values in groups.items():
            LOGGER.debug("%s.%s = %s", groupname, name, values)


def timed_call(name, function, *args, **kwargs):
    """Ruft eine bestehende Berechnungsfunktion genau einmal auf."""
    started = perf_counter()
    LOGGER.debug("Berechnung gestartet: %s", name)
    try:
        result = function(*args, **kwargs)
    except Exception:
        LOGGER.exception("Berechnung fehlgeschlagen: %s", name)
        raise
    finally:
        LOGGER.debug("Laufzeit %s: %.4f s", name, perf_counter() - started)
    return result


def _as_points(name, value, ndim=2):
    a = np.asarray(value)
    if a.ndim != ndim or a.shape[-1] != 3:
        raise AssertionError(f"{name}: Erwartet {ndim}D-Array mit letzter Achse 3, erhalten {a.shape}")
    if not np.issubdtype(a.dtype, np.number):
        raise AssertionError(f"{name}: Kein numerischer Datentyp: {a.dtype}")
    finite = np.isfinite(a)
    LOGGER.debug("%s: shape=%s dtype=%s NaN=%d Inf=%d", name, a.shape, a.dtype,
                 int(np.isnan(a).sum()), int(np.isinf(a).sum()))
    if a.size:
        rows_ok = finite.all(axis=-1)
        LOGGER.debug("%s: vollständige Zustände=%d/%d", name,
                     int(rows_ok.sum()), rows_ok.size)
        for j, symbol in enumerate(("T [K]", "v [m³/mol]", "p [Pa]")):
            vals = a[..., j]
            valid = vals[np.isfinite(vals)]
            if valid.size:
                LOGGER.debug("%s %s: min=%.8g max=%.8g", name, symbol, valid.min(), valid.max())
                if np.any(valid <= 0):
                    LOGGER.warning("%s: %s enthält %d Werte <= 0 (logarithmisch ggf. nicht darstellbar)",
                                   name, symbol, int((valid <= 0).sum()))
    return a


def pruefe_isolinien(daten):
    """Überführt die Isolinien-Prüfungen aus der bisherigen main.py."""
    for name, spalte in (("isothermen", 0), ("isochoren", 1)):
        kurven = _as_points(name, daten[name], ndim=3)
        if not kurven.shape[0] or not kurven.shape[1]:
            LOGGER.warning("%s: kein auswertbares Raster", name)
            continue
        konstant = kurven[..., spalte]
        ok = np.allclose(konstant, konstant[:, :1], rtol=1e-12, atol=1e-14)
        LOGGER.debug("%s: konstante Koordinate korrekt=%s", name, ok)
        if not ok:
            raise AssertionError(f"{name}: Isolinien-Koordinate variiert")

    for index, kurve in enumerate(daten["isobaren"]):
        points = _as_points(f"isobaren[{index}]", kurve)
        if not len(points):
            LOGGER.warning("Isobare %d enthält keine Punkte", index)
            continue
        ok = np.allclose(points[:, 2], points[0, 2], rtol=1e-12, atol=1e-5)
        LOGGER.debug("Isobare %d: p konstant=%s", index, ok)
        if not ok:
            raise AssertionError(f"Isobare {index}: Druck nicht konstant")

    liquid = _as_points("binodal_liquid", daten["binodal_liquid"])
    vapor = _as_points("binodal_vapor", daten["binodal_vapor"])
    if liquid.shape != vapor.shape:
        raise AssertionError("Binodale: Äste besitzen keine kompatible Form")
    if len(liquid) < 2:
        LOGGER.warning("Binodale: weniger als zwei gültige Punkte; Zweiphasenfläche ausgelassen")
        return
    np.testing.assert_allclose(liquid[:, 0], vapor[:, 0], rtol=0, atol=1e-8)
    np.testing.assert_allclose(liquid[:, 2], vapor[:, 2], rtol=1e-8, atol=1e-4)
    # Bei einem angehaengten kritischen Punkt ist v_l == v_g zulassig.
    critical_is_last = bool(len(liquid) and np.array_equal(liquid[-1], vapor[-1]))
    last = -1 if critical_is_last else len(liquid)
    if not np.all(liquid[:last, 1] < vapor[:last, 1]):
        raise AssertionError("Unterkritisch muss v_flüssig kleiner als v_dampf sein")
    LOGGER.debug("Binodale: %d passende Flüssig-/Dampfpunkte", len(liquid))


def pruefe_uebereinstimmung(daten):
    """Prüft gemeinsame (T,v)-Stützstellen von Fläche und Isothermen."""
    surface = _as_points("surface", daten["surface"], ndim=3)
    iso = _as_points("isothermen", daten["isothermen"], ndim=3)
    if not len(surface) or not len(iso):
        return
    gemeinsame_T = np.intersect1d(surface[:, 0, 0], iso[:, 0, 0])
    if not gemeinsame_T.size:
        LOGGER.info("Oberfläche/Isothermen: keine gemeinsame Testtemperatur")
        return
    T_test = gemeinsame_T[0]
    surf = surface[np.flatnonzero(np.isclose(surface[:, 0, 0], T_test))[0]]
    line = iso[np.flatnonzero(np.isclose(iso[:, 0, 0], T_test))[0]]
    indices = np.searchsorted(surf[:, 1], line[:, 1])
    contained = indices < len(surf)
    match = np.zeros(len(indices), dtype=bool)
    match[contained] = np.isclose(surf[indices[contained], 1], line[contained, 1], rtol=0, atol=1e-12)
    if not np.any(match):
        LOGGER.info("Oberfläche/Isothermen: keine gemeinsamen Volumenpunkte")
        return
    actual, expected = surf[indices[match], 2], line[match, 2]
    ok = np.allclose(actual, expected, rtol=1e-6, atol=1e-2, equal_nan=True)
    LOGGER.debug("Oberfläche/Isothermen: T=%.5g K gemeinsame v=%d übereinstimmend=%s",
                 T_test, int(match.sum()), ok)
    if not ok:
        raise AssertionError("Isotherme und Fläche liefern unterschiedliche Drücke")


def diagnose_daten(daten):
    """Struktur, Datentypen, physikalische Plausibilität und Phasenmaske."""
    _as_points("surface", daten["surface"], ndim=3)
    _as_points("isothermen", daten["isothermen"], ndim=3)
    _as_points("isochoren", daten["isochoren"], ndim=3)
    for i, kurve in enumerate(daten["isobaren"]):
        _as_points(f"isobaren[{i}]", kurve)
    for name in ("binodal_liquid", "binodal_vapor", "spinodale_liquid", "spinodale_vapor"):
        _as_points(name, daten[name])
    phase_mask = np.asarray(daten["phase_mask"])
    if phase_mask.shape != daten["surface"].shape[:2] or phase_mask.dtype != np.bool_:
        raise AssertionError("phase_mask hat unerwartete Shape oder dtype")
    LOGGER.debug("phase_mask: shape=%s Zweiphasenpunkte=%d", phase_mask.shape, int(phase_mask.sum()))
    liquid = np.asarray(daten["binodal_liquid"])
    vapor = np.asarray(daten["binodal_vapor"])
    two_phase = _as_points("two_phase_surface", daten["two_phase_surface"], ndim=3)

    Q = np.asarray(daten["Q_values"])

    n_sat = two_phase.shape[0]

    if (
        Q.ndim != 1
        or len(Q) != two_phase.shape[1]
        or len(liquid) not in (n_sat, n_sat + 1)
    ):
        raise AssertionError(
            "Zweiphasenfläche und Binodale "
            "besitzen inkompatible Dimensionen."
        )

    if n_sat >= 2:

        np.testing.assert_allclose(
            Q[[0, -1]],
            [0.0, 1.0],
            rtol=0,
            atol=1e-12
        )

        np.testing.assert_allclose(
            two_phase[:, 0, :],
            liquid[:n_sat],
            rtol=1e-10,
            atol=1e-6
        )

        np.testing.assert_allclose(
            two_phase[:, -1, :],
            vapor[:n_sat],
            rtol=1e-10,
            atol=1e-6
        )

        if not (
            np.isfinite(two_phase).all()
            and np.all(two_phase > 0)
        ):
            raise AssertionError(
                "Ungültige Zweiphasenkoordinaten."
            )
        
        # Rand-/Meshpruefungen erfolgen separat in diagnose_visualisierung().
        # Die Variablen surface, transform und critical existieren hier nicht.
        
    Tc, pc, vc = daten["critical_point"]
    if not np.isfinite([Tc, pc, vc]).all() or min(Tc, pc, vc) <= 0:
        raise AssertionError("Kritischer Punkt ist ungültig")
    LOGGER.debug("Kritischer Punkt: Tc=%.8g K pc=%.8g Pa vc=%.8g m³/mol", Tc, pc, vc)
    pruefe_isolinien(daten)
    pruefe_uebereinstimmung(daten)


def test_saturation_plateau(stoff, T, calc_saturation_point, p_eos):
    """Alter Test_1: Sättigungsdruck beim mittleren Zweiphasenvolumen prüfen.

    Thermodynamische Aufrufe werden als bestehende Funktionen übergeben,
    deshalb benötigt Debug.py keinen Import aus coolprop_heos.py.
    """
    sat = calc_saturation_point(T, stoff)
    if sat is None:
        LOGGER.warning("Plateau-Test: kein Sättigungspunkt bei T=%g K", T)
        return None
    liquid, vapor = sat
    p_sat = liquid[2]
    v_mid = 0.5 * (liquid[1] + vapor[1])
    p_test = p_eos(T, v_mid, stoff)
    rel = abs(p_test - p_sat) / p_sat
    LOGGER.debug("Plateau-Test: stoff=%s T=%g K p_sat=%g Pa p(T,v_mid)=%g Pa rel. Fehler=%.3g",
                 stoff, T, p_sat, p_test, rel)
    if not np.isfinite(p_test) or not np.isclose(p_test, p_sat, rtol=1e-3):
        raise AssertionError("HEOS-Sättigungsdruck und DmolarT-Druck stimmen nicht überein")
    return rel


def test_clapeyron(stoff, T):
    """HEOS-Sättigungsableitung gegen Clapeyron (alter Test aus coolprop_heos.py)."""
    import CoolProp.CoolProp as CP
    from CoolProp import AbstractState

    state = AbstractState("HEOS", stoff)
    state.update(CP.QT_INPUTS, 0.0, float(T))
    h_l, v_l = state.hmolar(), 1.0 / state.rhomolar()
    dp_dT_heos = state.first_saturation_deriv(CP.iP, CP.iT)
    state.update(CP.QT_INPUTS, 1.0, float(T))
    h_g, v_g = state.hmolar(), 1.0 / state.rhomolar()
    dp_dT_clapeyron = (h_g - h_l) / (float(T) * (v_g - v_l))
    rel_error = abs(dp_dT_heos - dp_dT_clapeyron) / max(abs(dp_dT_heos), np.finfo(float).tiny)
    LOGGER.debug("Clapeyron: stoff=%s T=%g K HEOS=%g Pa/K Gleichung=%g Pa/K rel. Fehler=%.3g",
                 stoff, T, dp_dT_heos, dp_dT_clapeyron, rel_error)
    if not np.isclose(dp_dT_heos, dp_dT_clapeyron, rtol=1e-5):
        raise AssertionError("Clapeyron-Diagnose fehlgeschlagen")
    return rel_error


def test_spinodale_ableitung(stoff, liquid, vapor, critical, *,
                             tolerance=1e-2, critical_exclusion_rel=1e-3,
                             max_samples=128, state_factory=None, cp_module=None):
    """AP3: Prueft R=|(rho_c/p_c)*(dp/drho_m)_T| dimensionslos.

    Nur unterkritische, endliche homogene Zustaende ausserhalb eines
    relativen kritischen Fensters werden ausgewertet. Die Druckableitung
    wird bei PHASENVORGABE aus HEOS ermittelt (kein Gleichgewichts-Flash).
    Nicht auswertbare Punkte sind keine bestandenen Tests.

    Liefert fuer beide Aeste Zaehler und pass/fail, ohne die GUI zu importieren.
    Durch optionale state_factory/cp_module sind reine Mocktests moeglich.
    """
    if not np.isfinite([tolerance, critical_exclusion_rel]).all() or not (
            tolerance > 0 and 0 < critical_exclusion_rel < 0.1):
        raise ValueError("AP3: Ungueltige Residualtoleranz / kritischer Ausschluss.")
    if type(max_samples) is not int or max_samples < 2:
        raise ValueError("AP3: max_samples muss >=2 sein.")
    Tc, pc, vc = map(float, critical)
    if not np.isfinite([Tc, pc, vc]).all() or min(Tc, pc, vc) <= 0:
        raise ValueError("AP3: Ungueltiger kritischer Punkt.")

    if cp_module is None:
        import CoolProp.CoolProp as cp_module
    if state_factory is None:
        from CoolProp import AbstractState
        state_factory = lambda: AbstractState("HEOS", stoff)
    state = state_factory()
    rho_c = 1.0 / vc
    cutoff = Tc - max(0.5, Tc * critical_exclusion_rel)
    results = {}

    for name, points, phase in (("liquid", liquid, cp_module.iphase_liquid),
                                ("vapor", vapor, cp_module.iphase_gas)):
        arr = np.asarray(points, dtype=float)
        if arr.size == 0:
            arr = np.empty((0, 3), dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 3:
            raise ValueError(f"AP3: {name}-Ast benoetigt (n,3)-Punkte.")
        # Die kritische Endzeile ist absichtlich NICHT testbar: dort liegen
        # numerisch singulaere Grenzverhaeltnisse vor.
        selected = np.flatnonzero(
            np.isfinite(arr).all(axis=1) & (arr[:, 0] > 0)
            & (arr[:, 1] > 0) & (arr[:, 0] < cutoff)
        )
        excluded = int(len(arr) - len(selected))
        n_candidates = int(len(selected))
        if len(selected) > max_samples:
            selected = selected[np.unique(np.linspace(0, len(selected)-1,
                                                        max_samples, dtype=int))]
        residuals = []
        errors = 0
        for index in selected:
            T, v, _p = arr[index]
            try:
                state.specify_phase(phase)
                state.update(cp_module.DmolarT_INPUTS, float(1.0/v), float(T))
                deriv = float(state.first_partial_deriv(
                    cp_module.iP, cp_module.iDmolar, cp_module.iT))
                residual = abs(deriv) * rho_c / pc
                if not np.isfinite(residual):
                    raise ValueError("Nichtendliches Ableitungsresiduum")
                residuals.append(float(residual))
            except (ValueError, RuntimeError, OverflowError) as exc:
                errors += 1
                LOGGER.debug("AP3 %s Ableitung Index %d: %s", name, index, exc)
            finally:
                state.unspecify_phase()
        maximum = float(max(residuals)) if residuals else None
        passed = bool(len(residuals) >= 2 and errors == 0
                      and maximum is not None and maximum <= tolerance)
        results[name] = {
            "candidates": n_candidates, "sampled": int(len(selected)),
            "excluded": excluded, "evaluated": len(residuals), "errors": errors,
            "max_residual": maximum,
            "median_residual": float(np.median(residuals)) if residuals else None,
            "tolerance": float(tolerance), "passed": passed,
        }
        log = LOGGER.info if passed else LOGGER.warning
        log("AP3 Spinodale %s: %s, R_max=%s <= %.3g; %d/%d Stichproben, Fehler=%d, ausgeschlossen=%d",
            name, "BESTANDEN" if passed else "NICHT VALIDIERT",
            f"{maximum:.4g}" if maximum is not None else "n/a", tolerance,
            len(residuals), n_candidates, errors, excluded)
    return results


def diagnose_homogeneous(daten):
    """Numerische Verluste und Oberflaechenzellen ohne VTK protokollieren."""
    from homogeneous_pvt import (count_renderable_triangles, grid_triangles, node_validity,
                                 mesh_topology_report)
    meta = daten.get('homogeneous_validation', {})
    axis = daten.get('pressure_axis', 'log')
    critical = daten.get('critical_point')
    for key in ('metastable_liquid_surface', 'metastable_vapor_surface', 'unstable_surface'):
        surface = np.asarray(daten.get(key, np.empty((0, 0, 3))))
        if surface.ndim != 3 or surface.shape[-1] != 3:
            raise AssertionError(f'{key}: ungueltige Arraystruktur')
        finite = np.isfinite(surface).all(axis=2)
        full = finite.all(axis=1)
        borders = finite[:, 0] & finite[:, -1] if surface.shape[1] else np.zeros(len(surface), dtype=bool)
        n_tri = count_renderable_triangles(surface, axis)
        n_negative = int(np.count_nonzero(finite & (surface[..., 2] <= 0)))
        # NaN-Trennspalten (EOS-Artefakt) sind erlaubt; nur endliche Knoten muessen steigen.
        for row in surface[borders]:
            if not np.all(np.diff(row[np.isfinite(row).all(axis=1), 1]) > 0):
                raise AssertionError(f'{key}: Volumen nicht monoton')
        info = meta.get('surfaces', {}).get(key, {})
        gaps = info.get('gap_rows', {})
        LOGGER.info('%s: vollstaendige=%d partielle=%d EOS-Artefaktzeilen=%d ohne-Rand=%d '
                    'gueltige-Knoten=%d/%d p<=0=%d Dreiecke=%d',
                    key, int(full.sum()), int(np.count_nonzero(borders & ~full)) - len(gaps),
                    len(gaps), int(np.count_nonzero(~borders)), int(finite.sum()),
                    int(finite.size), n_negative, n_tri)
        coverage = info.get('coverage')
        if coverage:
            LOGGER.info('%s: Zellabdeckung=%s (%d/%d auswertbare Zellen), ausgeschlossene '
                        'EOS-Artefaktzellen=%d', key,
                        'n/a' if coverage['coverage'] is None else f"{coverage['coverage']:.2%}",
                        coverage['cells_covered'], coverage['cells_evaluable'],
                        coverage['cells_excluded_gap'])
        if critical is not None and surface.shape[0] > 1:
            Tc, pc, vc = critical
            valid = node_validity(surface, axis)
            xyz = np.stack((np.log10(np.where(valid, surface[..., 1], 1.) / vc),
                            surface[..., 0] / Tc), axis=-1).reshape(-1, 2)
            tri, cells = grid_triangles(valid)
            top = mesh_topology_report(tri, cells, valid, xy=xyz)
            LOGGER.info('%s: Dreiecke=%d Komponenten=%d offene Kanten (Rand/Luecke/unerklaert)='
                        '%d/%d/%d Orientierungskonflikte=%d Projektionsfaltungen=%d', key,
                        top['triangles'], top['components'], top['open_edges_frame'],
                        top['open_edges_gap'], top['open_edges_unexplained'],
                        top['orientation_conflicts'], top['projected_orientation_flips'])
            if top['open_edges_unexplained'] or top['orientation_conflicts'] \
                    or top['projected_orientation_flips']:
                LOGGER.warning('%s: Netzfehler erkannt: %s', key, top)
        for i, errors in info.get('point_reasons', {}).items():
            LOGGER.debug('%s Temperaturindex %s: verworfene Punkte %s', key, i, errors)
        for i, reason in info.get('reasons', {}).items():
            LOGGER.warning('%s Temperaturindex %s: Randproblem %s', key, i, reason)
        for i, gap in gaps.items():
            LOGGER.debug('%s Temperaturindex %s: EOS-Artefaktluecke v in (%s, %s), p-Scan '
                         '[%s, %s], Vorzeichenwechsel=%s', key, i, gap.get('v_a'), gap.get('v_b'),
                         gap.get('scan_p_min'), gap.get('scan_p_max'),
                         gap.get('slope_sign_changes'))
        check = info.get('independent_check')
        if check and not check.get('passed', True):
            LOGGER.warning('%s: unabhaengige Knotenpruefung fehlgeschlagen: %s', key, check)
    combined = meta.get('combined')
    if combined:
        LOGGER.info('Zusammengefuehrte Flaeche: Komponenten=%s (erwartet %s), offene Kanten zwischen '
                    'gemeinsamen Knoten=%s', combined.get('components'),
                    combined.get('components_expected'),
                    combined.get('open_edges_between_shared_nodes'))
    if daten.get('critical_closure'):
        LOGGER.info('Kritischer Punkt: %s', daten['critical_closure'])
    if meta.get('error'):
        LOGGER.warning('Homogene Berechnung fehlgeschlagen: %s', meta['error'])
    return meta


def diagnose_visualisierung(daten, actors, transform):
    """Pruft zweiphasige Mesh-Geometrie bei optionalem kritischem Binodal-Endpunkt."""
    import pyvista as pv

    required = {"surface", "isothermen", "isochoren", "isobaren", "binodale",
                "spinodale", "zweiphasengebiet", "critical_point"}
    missing = required - set(actors)
    if missing:
        raise AssertionError(f"Darstellungs-Ebenen fehlen: {missing}")
    for key, actor in actors.items():
        amount = len(actor) if isinstance(actor, (list, tuple)) else int(actor is not None)
        LOGGER.debug("Visualisierung %s: Actor-Anzahl=%d", key, amount)

    liquid = np.asarray(daten["binodal_liquid"], dtype=float)
    vapor = np.asarray(daten["binodal_vapor"], dtype=float)
    surface = np.asarray(daten["two_phase_surface"], dtype=float)
    critical = daten["critical_point"]  # (Tc, pc, vc)!

    if surface.shape == (0, 0, 3):
        LOGGER.debug("Visualisierung: keine Zweiphasenfläche verfügbar")
        return
    if surface.ndim != 3 or surface.shape[-1] != 3 or min(surface.shape[:2]) < 2:
        raise AssertionError("Zweiphasenfläche braucht (n_T>=2, n_Q>=2, 3).")
    n_T, n_Q = surface.shape[:2]
    if liquid.shape != vapor.shape or liquid.shape not in ((n_T, 3), (n_T + 1, 3)):
        raise AssertionError("Binodalarrays benötigen n_T oder n_T+1 Zeilen.")

    # Die letzte Binodalzeile kann der kritische Punkt sein, während die
    # Zweiphasenfläche nur die unterkritischen Profile enthält.
    np.testing.assert_allclose(surface[:, 0, :], liquid[:n_T], rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(surface[:, -1, :], vapor[:n_T], rtol=1e-10, atol=1e-12)

    if len(liquid) == n_T + 1:
        Tc, pc, vc = critical
        endpoint = np.array([Tc, vc, pc], dtype=float)
        np.testing.assert_allclose(liquid[-1], endpoint, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(vapor[-1], endpoint, rtol=1e-12, atol=1e-12)
        if not surface[-1, 0, 0] < Tc:
            raise AssertionError("Letzte Q-Zeile muss unterhalb Tc liegen.")

    xyz = transform(surface, critical)
    if not np.isfinite(xyz).all():
        raise AssertionError("Nichtendliche Zweiphasenkoordinaten")
    mesh = pv.StructuredGrid(xyz[..., 0], xyz[..., 1], xyz[..., 2])
    if mesh.n_points != n_T * n_Q or mesh.n_cells != (n_T - 1) * (n_Q - 1):
        raise AssertionError("Geometrischer Mesh-Test fehlgeschlagen")
    areas = np.asarray(mesh.compute_cell_sizes().cell_data["Area"])
    if not np.isfinite(areas).all() or not np.all(areas > 0):
        raise AssertionError("Mesh enthält ungültige oder entartete Flächen")
    LOGGER.debug("Zweiphasen-Mesh: Punkte=%d Zellen=%d Flächen min/max=%g/%g",
                 mesh.n_points, mesh.n_cells, areas.min(), areas.max())


def create_test_data():
    """Synthetischer Raster-Datensatz für unabhängige Koordinatentests (kein HEOS)."""

    T = np.linspace(400.0, 600.0, 9)
    v = np.geomspace(2e-5, 5e-4, 15)

    T_grid, V_grid = np.meshgrid(T, v, indexing="ij")

    P = 5e6 * (T_grid / 500.0) * (1e-4 / V_grid)

    surface = np.stack(
        (T_grid, V_grid, P),
        axis=-1
    )

    return {
        "surface": surface,
        "isothermen": surface[::4, :, :],
        "isochoren": np.swapaxes(surface[:, ::5, :], 0, 1),
        "isobaren": [],
        "binodal_liquid": np.empty((0, 3)),
        "binodal_vapor": np.empty((0, 3)),
        "critical_point": (500.0, 5e6, 1e-4)
    }

