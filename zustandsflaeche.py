"""
Datei: zustandsflaeche.py
Projekt: ThermoSurface
Kategorie: Thermodynamik / Flächenmodell

Beschreibung:
Verwaltet rohe homogene p-v-T-Daten, Stabilitätsklassen und Ersatzflächen.

Aufgaben:
- Trennt HEOS-Rohwerte von Anzeige- und Stabilitätsfiltern.
- Konstruiert Ersatzflächen und gemeinsame Netzpunkte bis zum kritischen Rand.
Einheiten: T [K], v [m³/mol], p [Pa].

Abhängigkeiten:
- Projektmodule (direkt): coolprop_heos_common.py, homogeneous_pvt.py
- Bibliotheken: CoolProp, numpy, scipy

Verwendet von:
- Projektmodule: diagram_axes.py, diagram_coords.py, gui.py, homogeneous_pvt.py,
  main.py, visualization.py
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

# (a) Auswertungsstatus des Rohwerts
EVAL_OK, EVAL_ERROR, EVAL_NONFINITE = 0, 1, 2
EVAL_LABELS = {EVAL_OK: 'ok', EVAL_ERROR: 'Fehler', EVAL_NONFINITE: 'nicht endlich'}
# (b) Stabilitaetsklasse
STAB_STABLE, STAB_METASTABLE, STAB_UNSTABLE, STAB_UNDETERMINED = 0, 1, 2, 3
STAB_LABELS = {STAB_STABLE: 'stabil', STAB_METASTABLE: 'metastabil',
               STAB_UNSTABLE: 'instabil', STAB_UNDETERMINED: 'unbestimmt'}
# (c) Filter-Grundcode (0 = gezeigt)
FILTER_NONE, FILTER_NOT_EVALUATED, FILTER_STABLE_SIGN, FILTER_PRESSURE_WINDOW, \
    FILTER_DISCONNECTED, FILTER_NO_ROW = 0, 1, 2, 3, 4, 5
FILTER_LABELS = {FILTER_NONE: 'gezeigt', FILTER_NOT_EVALUATED: 'Rohwert nicht auswertbar',
                 FILTER_STABLE_SIGN: 'dp/dv<0 (mechanisch stabil, EOS-Zusatzextremum)',
                 FILTER_PRESSURE_WINDOW: 'p ausserhalb [p_spin,l; p_spin,g]',
                 FILTER_DISCONNECTED: 'nicht zusammenhaengend ab Spinodale',
                 FILTER_NO_ROW: 'Temperaturzeile nicht konstruierbar'}
# (d) Herkunft
SRC_HEOS, SRC_SUBSTITUTE = 0, 1
# TS-VP05: Schnittpunkt eines viewportbegrenzten Netzes = geometrische Interpolation im Renderer, KEIN HEOS-Knoten.
SRC_CLIPPED = 2
SRC_LABELS = {SRC_HEOS: 'HEOS', SRC_SUBSTITUTE: 'Ersatz'}
SUBSTITUTE_HOVER_LABEL = 'Ersatzfläche – keine EOS-Auswertung'
CLIPPED_HOVER_LABEL = 'Schnittpunkt (Renderer-Interpolation) – keine HEOS-Auswertung'
# Zeilenstatus
ROW_OK, ROW_NO_SPINODAL, ROW_ORDER, ROW_COLLAPSED, ROW_CRITICAL, ROW_SUPERCRITICAL, ROW_NOT_FINITE = range(7)
ROW_LABELS = {ROW_OK: 'ok', ROW_NO_SPINODAL: 'Spinodale fehlt', ROW_ORDER: 'Ersatz nicht konstruierbar (Ordnung)',
              ROW_COLLAPSED: 'kollabiert nahe T_c', ROW_CRITICAL: 'kritische Temperatur',
              ROW_SUPERCRITICAL: 'ueberkritisch (keine Spinodale)', ROW_NOT_FINITE: 'nicht endlich'}

# Standardschwellen nahe T_c (aus scripts/diagnose_zustandsflaeche.py / Doku begruendet)
EPS_WIDTH = 1e-2      # (v_V - v_L) / v_c
EPS_PRESSURE = 1e-8   # (p_V - p_L) / p_c; 10x ueber der Druckaufloesung 1e-9 der Spinodalen
SYMLOG_LINTHRESH = 1e-2  # Linearschwelle p0/p_c der symlog-Achse


@dataclass
class HomogeneousField:
    """Alle Felder (n_T, n_s) ausser T und row_*; Rohwerte p_raw/dpdv bleiben unveraendert."""
    name: str
    T: np.ndarray
    v: np.ndarray
    p_raw: np.ndarray
    dpdv: np.ndarray
    eval_status: np.ndarray
    eval_reason: dict
    stability: np.ndarray
    display: np.ndarray
    filter_code: np.ndarray
    source: np.ndarray
    row_status: np.ndarray
    row_reason: dict = field(default_factory=dict)

    @property
    def shape(self):
        return self.v.shape

    def raw_points(self):
        """(n_T,n_s,3) mit allen Rohwerten (T,v,p_raw); NaN nur wo die Auswertung scheiterte."""
        out = np.empty(self.v.shape + (3,))
        out[..., 0] = self.T[:, None]
        out[..., 1] = self.v
        out[..., 2] = self.p_raw
        return out

    def display_points(self):
        """Abgeleitete Ansicht: ausgeblendete Knoten NaN. Aendert p_raw nicht (Kopie)."""
        out = self.raw_points()
        out[~self.display] = np.nan
        return out

    def counts(self):
        n = self.display.size
        by_source = {SRC_LABELS[k]: int(np.count_nonzero(self.source == k)) for k in SRC_LABELS}
        return dict(nodes=int(n), shown=int(self.display.sum()),
                    hidden=int(n - self.display.sum()),
                    eval_ok=int(np.count_nonzero(self.eval_status == EVAL_OK)),
                    eval_error=int(np.count_nonzero(self.eval_status == EVAL_ERROR)),
                    eval_nonfinite=int(np.count_nonzero(self.eval_status == EVAL_NONFINITE)),
                    by_filter={FILTER_LABELS[k]: int(np.count_nonzero(self.filter_code == k))
                               for k in FILTER_LABELS},
                    by_stability={STAB_LABELS[k]: int(np.count_nonzero(self.stability == k))
                                  for k in STAB_LABELS},
                    by_source=by_source, source_sum=int(sum(by_source.values())),
                    rows={ROW_LABELS[k]: int(np.count_nonzero(self.row_status == k)) for k in ROW_LABELS})


def raw_digest(field_or_array):
    """Bitgleicher Fingerabdruck der Rohwerte (Test: vor/nach Filter und Vernetzung identisch)."""
    import hashlib
    arr = field_or_array.p_raw if isinstance(field_or_array, HomogeneousField) else field_or_array
    return hashlib.sha256(np.ascontiguousarray(arr, dtype=np.float64).tobytes()).hexdigest()


def _empty_field(name, T, n_s, source):
    n_T = len(T)
    return HomogeneousField(
        name=name, T=np.asarray(T, dtype=float), v=np.full((n_T, n_s), np.nan),
        p_raw=np.full((n_T, n_s), np.nan), dpdv=np.full((n_T, n_s), np.nan),
        eval_status=np.full((n_T, n_s), EVAL_ERROR, dtype=np.int8), eval_reason={},
        stability=np.full((n_T, n_s), STAB_UNDETERMINED, dtype=np.int8),
        display=np.zeros((n_T, n_s), dtype=bool),
        filter_code=np.full((n_T, n_s), FILTER_NO_ROW, dtype=np.int8),
        source=np.full((n_T, n_s), source, dtype=np.int8),
        row_status=np.full(n_T, ROW_NO_SPINODAL, dtype=np.int8))


def _row_boundaries_ok(sl_row, sg_row):
    return bool(np.isfinite(sl_row).all() and np.isfinite(sg_row).all())


# ---------------------------------------------------------------------------
# (a) Rohfeld der HEOS-Fortsetzung zwischen den Spinodalen + (b),(c) Filterschritt
# ---------------------------------------------------------------------------

def build_raw_field(stoff, T, spinodal_liquid, spinodal_vapor, n_s, *, critical_temperature=None,
                    state_factory=None, cp_module=None):
    """Rohwerte p_raw auf dem festen Indexgitter (T_i, j), j=0..n_s-1, v=geomspace(v_L, v_V).

    Randknoten j=0 / n_s-1 sind die Spinodalpunkte selbst (identische Werte, einmal berechnet).
    Es wird nichts gefiltert: Fehlschlaege setzen nur den Auswertungsstatus (p_raw=NaN, Grund).
    """
    from homogeneous_pvt import homogeneous_point
    sl = np.asarray(spinodal_liquid, dtype=float)
    sg = np.asarray(spinodal_vapor, dtype=float)
    T = np.asarray(T, dtype=float)
    if sl.shape != (len(T), 3) or sg.shape != (len(T), 3):
        raise ValueError('Spinodalen muessen die Form (n_T,3) haben')
    if not isinstance(n_s, int) or n_s < 3:
        raise ValueError('n_s muss eine ganze Zahl >= 3 sein')
    if cp_module is None:
        import CoolProp.CoolProp as cp_module
    if state_factory is None:
        from coolprop_heos_common import _new_state
        state_factory = lambda: _new_state(stoff)
    state = state_factory()
    Tc = float(state.T_critical()) if critical_temperature is None else float(critical_temperature)
    fld = _empty_field('unstable_heos_raw', T, n_s, SRC_HEOS)
    for i in range(len(T)):
        if T[i] >= Tc:
            fld.row_status[i] = ROW_SUPERCRITICAL
            fld.row_reason[i] = 'T >= T_c: keine Spinodale'
            continue
        if not _row_boundaries_ok(sl[i], sg[i]):
            fld.row_reason[i] = 'Spinodalpunkt nicht berechenbar (NaN)'
            continue
        if not (0 < sl[i, 1] < sg[i, 1] and sl[i, 2] < sg[i, 2]):
            fld.row_status[i] = ROW_ORDER
            fld.row_reason[i] = 'Spinodalordnung v_l<v_g bzw. p_l<p_g verletzt'
            continue
        fld.row_status[i] = ROW_OK
        v = np.geomspace(sl[i, 1], sg[i, 1], n_s)
        v[0], v[-1] = sl[i, 1], sg[i, 1]
        fld.v[i] = v
        fld.p_raw[i, 0], fld.p_raw[i, -1] = sl[i, 2], sg[i, 2]
        fld.dpdv[i, 0] = fld.dpdv[i, -1] = 0.
        fld.eval_status[i, [0, -1]] = EVAL_OK
        for j in range(1, n_s - 1):
            try:
                p, slope = homogeneous_point(state, cp_module, float(T[i]), float(v[j]), 'vapor')
                fld.p_raw[i, j] = p
                fld.dpdv[i, j] = -(1. / v[j]) ** 2 * slope
                fld.eval_status[i, j] = EVAL_OK
            except (ValueError, RuntimeError, OverflowError) as exc:
                fld.eval_status[i, j] = EVAL_NONFINITE if 'nicht auswertbar' in str(exc) else EVAL_ERROR
                fld.eval_reason[(i, j)] = str(exc)
    return fld


def apply_display_filter(fld, spinodal_liquid, spinodal_vapor, pc, *, tol_rel=1e-9):
    """Setzt NUR (b) Stabilitaetsklasse und (c) Darstellungsstatus/Filtercode; p_raw bleibt unberuehrt.

    gezeigt <=> Rohwert ok, (dp/dv)_T>0, p in [p_spin,l; p_spin,g] und ueber eine ununterbrochene
    Kette gezeigter Knoten mit einer der beiden Spinodalen verbunden.
    """
    sl = np.asarray(spinodal_liquid, dtype=float)
    sg = np.asarray(spinodal_vapor, dtype=float)
    n_T, n_s = fld.shape
    fld.display[:] = False
    fld.filter_code[:] = FILTER_NO_ROW
    fld.stability[:] = STAB_UNDETERMINED
    for i in range(n_T):
        if fld.row_status[i] != ROW_OK:
            continue
        ok = fld.eval_status[i] == EVAL_OK
        slope = fld.dpdv[i]
        p = fld.p_raw[i]
        lo, hi = sl[i, 2], sg[i, 2]
        tol = tol_rel * max(abs(lo), abs(hi), pc)
        boundary = np.zeros(n_s, dtype=bool)
        boundary[[0, -1]] = True
        unstable = ok & (slope > 0)
        meta = ok & (slope < 0)
        fld.stability[i, unstable] = STAB_UNSTABLE
        fld.stability[i, meta] = STAB_METASTABLE
        in_window = ok & (p >= lo - tol) & (p <= hi + tol)
        shown = (unstable & in_window) | (boundary & ok)
        code = np.full(n_s, FILTER_NONE, dtype=np.int8)
        code[~ok] = FILTER_NOT_EVALUATED
        code[ok & ~boundary & ~unstable] = FILTER_STABLE_SIGN
        code[ok & ~boundary & unstable & ~in_window] = FILTER_PRESSURE_WINDOW
        run = np.zeros(n_s, dtype=bool)
        k = 0
        while k < n_s and shown[k]:
            run[k] = True
            k += 1
        k = n_s - 1
        while k >= 0 and shown[k]:
            run[k] = True
            k -= 1
        code[shown & ~run] = FILTER_DISCONNECTED
        fld.display[i] = run
        fld.filter_code[i] = np.where(run, FILTER_NONE, code)
    return fld


# ---------------------------------------------------------------------------
# (d) Ersatzflaeche (Modell): kubische Hermite-Verbindung der HEOS-Spinodalpunkte
# ---------------------------------------------------------------------------

def hermite_h(s):
    """h(s)=3s^2-2s^3 (C1, h(0)=0, h(1)=1, h'(0)=h'(1)=0, h'>=0 auf [0,1])."""
    s = np.asarray(s, dtype=float)
    return 3. * s * s - 2. * s ** 3


def hermite_dh(s):
    s = np.asarray(s, dtype=float)
    return 6. * s - 6. * s * s


def substitute_isotherm(v_l, p_l, v_v, p_v, n_s, *, v_c=None, p_c=None,
                        eps_width=EPS_WIDTH, eps_pressure=EPS_PRESSURE):
    """Ersatzisotherme p(s)=p_L+(p_V-p_L) h(s), v=v_L+s (v_V-v_L), s in [0,1].

    Rueckgabe (row_status, grund, v, p, dpdv). Bei Status != ROW_OK sind v,p,dpdv None: keine Punkte.
    """
    if not np.isfinite([v_l, p_l, v_v, p_v]).all():
        return ROW_NO_SPINODAL, 'Spinodalpunkt fehlt oder nicht endlich', None, None, None
    if not (0 < v_l < v_v):
        return ROW_ORDER, f'Ersatz nicht konstruierbar: v_V={v_v} <= v_L={v_l}', None, None, None
    if not p_v > p_l:
        return ROW_ORDER, f'Ersatz nicht konstruierbar: p_V={p_v} <= p_L={p_l}', None, None, None
    if v_c is not None and (v_v - v_l) < eps_width * v_c:
        return ROW_COLLAPSED, f'Breite {(v_v - v_l) / v_c:.3g} v_c < eps={eps_width}', None, None, None
    if p_c is not None and (p_v - p_l) < eps_pressure * p_c:
        return ROW_COLLAPSED, f'Druckdifferenz {(p_v - p_l) / p_c:.3g} p_c < eps={eps_pressure}', None, None, None
    s = np.linspace(0., 1., n_s)
    v = v_l + s * (v_v - v_l)
    p = p_l + (p_v - p_l) * hermite_h(s)
    v[0], v[-1], p[0], p[-1] = v_l, v_v, p_l, p_v
    dpdv = (p_v - p_l) / (v_v - v_l) * hermite_dh(s)
    if not (np.isfinite(p).all() and np.all(np.diff(p) >= 0) and np.all(np.diff(v) > 0)):
        return ROW_NOT_FINITE, 'Ersatzverlauf nicht monoton/endlich', None, None, None
    return ROW_OK, '', v, p, dpdv


def build_substitute_field(T, spinodal_liquid, spinodal_vapor, critical, n_s, *,
                           eps_width=EPS_WIDTH, eps_pressure=EPS_PRESSURE):
    """Ersatzflaeche auf dem Indexgitter (T_i, j); Randknoten = HEOS-Spinodalpunkte (identische Werte).

    critical=(T_c, p_c, v_c). Zeilen mit T>=T_c, fehlender/unplausibler Spinodale oder kollabierter
    Breite bekommen KEINE Punkte (Statuscode + Grund), niemals stillschweigend gefuellte Werte.
    """
    T = np.asarray(T, dtype=float)
    sl = np.asarray(spinodal_liquid, dtype=float)
    sg = np.asarray(spinodal_vapor, dtype=float)
    Tc, pc, vc = critical
    fld = _empty_field('replacement_surface', T, n_s, SRC_SUBSTITUTE)
    for i in range(len(T)):
        if T[i] > Tc * (1 + 1e-12):
            fld.row_status[i], fld.row_reason[i] = ROW_SUPERCRITICAL, 'T > T_c: keine instabile Zone'
            continue
        if T[i] >= Tc * (1 - 1e-12):
            fld.row_status[i], fld.row_reason[i] = ROW_CRITICAL, 'T = T_c: Zeile entspricht dem kritischen Punkt'
            continue
        status, reason, v, p, dpdv = substitute_isotherm(
            sl[i, 1], sl[i, 2], sg[i, 1], sg[i, 2], n_s, v_c=vc, p_c=pc,
            eps_width=eps_width, eps_pressure=eps_pressure)
        fld.row_status[i] = status
        if status != ROW_OK:
            fld.row_reason[i] = reason
            continue
        fld.v[i], fld.p_raw[i], fld.dpdv[i] = v, p, dpdv
        fld.eval_status[i] = EVAL_OK
        fld.stability[i] = np.where(dpdv > 0, STAB_UNSTABLE, STAB_UNDETERMINED)
        fld.display[i] = True
        fld.filter_code[i] = FILTER_NONE
        # Randknoten sind die unveraenderten HEOS-Spinodalpunkte, nur das Innere ist Modell.
        fld.source[i, [0, -1]] = SRC_HEOS
    return fld


# ---------------------------------------------------------------------------
# P2: gemeinsame Randpunkte (identische Vertex-Indizes)
# ---------------------------------------------------------------------------

@dataclass
class SharedVertices:
    vertices: np.ndarray            # (N,3) (T,v,p), jeder Randpunkt genau einmal referenziert
    ids: dict                       # Flaechenname -> (n_T,n_s) Vertex-Indizes
    apex_id: int | None             # Index des kritischen Punktes (oder None)
    mismatches: dict                # Verknuepfung -> Zeilenindizes mit ungleichen Randwerten


def assemble_shared_vertices(grids, links, critical_point=None):
    """grids: Name -> (n_T,n_s,3). links: [(a, b)] mit 'letzte Spalte von a == erste Spalte von b'.

    Beide Flaechen referenzieren dann DENSELBEN Vertex-Index (Wertgleichheit wird bitweise
    geprueft; bei Abweichung bleiben getrennte Vertices und die Zeile wird gemeldet).
    """
    ids, cursor, blocks = {}, 0, []
    for name, g in grids.items():
        g = np.asarray(g, dtype=float)
        ids[name] = np.arange(cursor, cursor + g.shape[0] * g.shape[1]).reshape(g.shape[:2])
        blocks.append(g.reshape(-1, 3))
        cursor += g.shape[0] * g.shape[1]
    vertices = np.concatenate(blocks) if blocks else np.empty((0, 3))
    parent = np.arange(len(vertices))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    mismatches = {}
    for a, b in links:
        ga, gb = np.asarray(grids[a]), np.asarray(grids[b])
        left, right = ga[:, -1], gb[:, 0]
        both = np.isfinite(left).all(axis=1) & np.isfinite(right).all(axis=1)
        same = both & np.all(left == right, axis=1)
        for i in np.flatnonzero(same):
            ra, rb = find(int(ids[a][i, -1])), find(int(ids[b][i, 0]))
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
        mismatches[f'{a}|{b}'] = [int(i) for i in np.flatnonzero(both & ~same)]
    for name in ids:
        ids[name] = np.vectorize(find, otypes=[np.int64])(ids[name]) if ids[name].size else ids[name]
    apex = None
    if critical_point is not None:
        Tc, pc, vc = critical_point
        vertices = np.vstack((vertices, [[Tc, vc, pc]]))
        apex = len(vertices) - 1
    return SharedVertices(vertices, ids, apex, mismatches)


# ---------------------------------------------------------------------------
# P2: Spinodalzuordnung (Pruefung gegen die Regel "erste Nullstelle ab Saettigung")
# ---------------------------------------------------------------------------

def first_zero_from_saturation(state, cp, T, rho_sat, rho_far, *, n=4000):
    """Erste Nullstelle von (dp/drho)_T, von rho_sat in Richtung rho_far gesucht (unabhaengig vom Solver).

    Fluessigkeitsspinodale: rho_sat=rho_sat,liq, rho_far=rho_c (erste Nullstelle zu kleineren rho).
    Dampfspinodale:         rho_sat=rho_sat,vap, rho_far=rho_c (erste Nullstelle zu groesseren rho).
    Rueckgabe rho der Nullstelle (brentq) oder None.
    """
    from homogeneous_pvt import homogeneous_point
    from scipy.optimize import brentq
    rho = np.geomspace(rho_sat * (1 + (1e-9 if rho_far > rho_sat else -1e-9)), rho_far, n)

    def slope(r):
        return homogeneous_point(state, cp, T, 1. / float(r), 'vapor')[1]

    prev_r, prev_s = None, None
    for r in rho:
        try:
            s = slope(r)
        except (ValueError, RuntimeError, OverflowError):
            prev_r, prev_s = None, None
            continue
        if prev_s is not None and np.sign(prev_s) != np.sign(s):
            return float(brentq(slope, prev_r, r, xtol=1e-12 * abs(r), rtol=1e-13, maxiter=200))
        prev_r, prev_s = r, s
    return None


# ---------------------------------------------------------------------------
# P4: Vernetzung mit festem Indexgitter und fester Diagonale
# ---------------------------------------------------------------------------

def grid_triangles_strict(shown):
    """Zwei Dreiecke je Zelle mit FESTER Diagonale (c0,c2); Dreieck nur, wenn alle drei Ecken gezeigt.

    Zellecken (T_i,j),(T_i+1,j),(T_i+1,j+1),(T_i,j+1) = c0..c3. Dreieck A=(c0,c1,c2), B=(c0,c2,c3).
    Rueckgabe (n,3) lokale Zeilen-Spalten-Indizes i*n_s+j. Fehlende Punkte erzeugen definierte
    Randkanten, keine Ueberbrueckung.
    """
    shown = np.asarray(shown, dtype=bool)
    if shown.ndim != 2 or min(shown.shape) < 2:
        return np.empty((0, 3), dtype=np.int64)
    idx = np.arange(shown.size, dtype=np.int64).reshape(shown.shape)
    c0, c1, c2, c3 = idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]
    s0, s1, s2, s3 = shown[:-1, :-1], shown[1:, :-1], shown[1:, 1:], shown[:-1, 1:]
    a = s0 & s1 & s2
    b = s0 & s2 & s3
    tri_a = np.column_stack((c0[a], c1[a], c2[a]))
    tri_b = np.column_stack((c0[b], c2[b], c3[b]))
    return np.concatenate((tri_a, tri_b)) if len(tri_a) + len(tri_b) else np.empty((0, 3), dtype=np.int64)


def surface_triangles(shown, ids):
    """Globale Vertex-Indizes der Dreiecke einer Flaeche (ids: (n_T,n_s) aus SharedVertices)."""
    local = grid_triangles_strict(shown)
    return np.asarray(ids).ravel()[local] if len(local) else local


def critical_fan(shown_last_row, ids_last_row, apex_id):
    """Dreiecksfaecher (j, apex, j+1) von der letzten Zeile zum kritischen Punkt (gleiche Orientierung
    wie die Gitterdreiecke B der letzten Zellzeile)."""
    shown_last_row = np.asarray(shown_last_row, dtype=bool)
    both = shown_last_row[:-1] & shown_last_row[1:]
    ids = np.asarray(ids_last_row)
    return np.column_stack((ids[:-1][both], np.full(int(both.sum()), apex_id, dtype=np.int64), ids[1:][both]))


def fan_applicable(row_status):
    """Faecher nur, wenn alle Zeilen nach der letzten OK-Zeile kollabiert/kritisch/ueberkritisch sind und
    mindestens eine kollabierte oder kritische Zeile folgt (T_c wird erreicht)."""
    row_status = np.asarray(row_status)
    ok = np.flatnonzero(row_status == ROW_OK)
    if not len(ok):
        return False, None
    tail = row_status[ok[-1] + 1:]
    reaches = np.any((tail == ROW_COLLAPSED) | (tail == ROW_CRITICAL))
    clean = np.all(np.isin(tail, (ROW_COLLAPSED, ROW_CRITICAL, ROW_SUPERCRITICAL)))
    return bool(reaches and clean), int(ok[-1])


def vertices_to_xyz(vertices, critical, axis):
    """Darstellungskoordinaten (log10(v/v_c), T/T_c, z(p)); p_raw bleibt in `vertices` unveraendert."""
    Tc, pc, vc = critical
    v = np.asarray(vertices, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.stack((np.log10(v[:, 1] / vc), v[:, 0] / Tc, pressure_to_z(v[:, 2], pc, axis)), axis=-1)


def mesh_report(triangles, vertices_xyz, *, xy=None):
    """Topologie ueber globale Vertex-Indizes: Komponenten, Randschleifen, innere Loecher, Entartung.

    innere Loecher = Randschleifen - Komponenten (Flaeche ist ein Graph ueber (v,T), also je Komponente
    ein Gebiet ohne Henkel; eine Randschleife ist die aeussere, jede weitere ein Loch).
    Entartet: Flaeche < 1e-12 x mittlere Flaeche (xyz in Darstellungskoordinaten).
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    tri = np.asarray(triangles, dtype=np.int64)
    rep = dict(triangles=int(len(tri)), vertices_used=0, components=0, boundary_loops=0, inner_holes=0,
               open_edges=0, nonmanifold_edges=0, orientation_conflicts=0, degenerate=0,
               projected_orientation_flips=0)
    if not len(tri):
        return rep
    n = int(max(tri.max() + 1, len(vertices_xyz)))
    rep['vertices_used'] = int(len(np.unique(tri)))
    edges = np.concatenate((tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]))
    g = coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(n, n))
    _, labels = connected_components(g, directed=False)
    used = np.unique(tri)
    comp_ids = np.unique(labels[used])
    rep['components'] = int(len(comp_ids))
    und = np.sort(edges, axis=1)
    key = und[:, 0] * n + und[:, 1]
    uniq, first, counts = np.unique(key, return_index=True, return_counts=True)
    rep['nonmanifold_edges'] = int(np.count_nonzero(counts > 2))
    _, dcounts = np.unique(edges[:, 0] * n + edges[:, 1], return_counts=True)
    rep['orientation_conflicts'] = int(np.count_nonzero(dcounts > 1))
    open_e = und[first[counts == 1]]
    rep['open_edges'] = int(len(open_e))
    if len(open_e):
        gb = coo_matrix((np.ones(len(open_e)), (open_e[:, 0], open_e[:, 1])), shape=(n, n))
        _, bl = connected_components(gb, directed=False)
        rep['boundary_loops'] = int(len(np.unique(bl[np.unique(open_e)])))
    rep['inner_holes'] = max(0, rep['boundary_loops'] - rep['components'])
    xyz = np.asarray(vertices_xyz, dtype=float)
    a, b, c = (xyz[tri[:, k]] for k in range(3))
    area2 = np.linalg.norm(np.cross(b - a, c - a), axis=1)
    mean = float(np.mean(area2)) if len(area2) else 0.
    rep['degenerate'] = int(np.count_nonzero(~np.isfinite(area2) | (area2 < 1e-12 * mean)))
    if xy is not None:
        p = [np.asarray(xy)[tri[:, k]] for k in range(3)]
        signed = ((p[1][:, 0] - p[0][:, 0]) * (p[2][:, 1] - p[0][:, 1]) -
                  (p[1][:, 1] - p[0][:, 1]) * (p[2][:, 0] - p[0][:, 0]))
        nz = signed[signed != 0]
        if len(nz):
            maj = 1. if np.count_nonzero(nz > 0) >= len(nz) / 2 else -1.
            rep['projected_orientation_flips'] = int(np.count_nonzero(nz * maj < 0))
    return rep


# ---------------------------------------------------------------------------
# P6: Druckachse (Rohwert bleibt unveraendert, nur z_display wird berechnet)
# ---------------------------------------------------------------------------

PRESSURE_AXES = ('log', 'asinh', 'linear', 'symlog')


def pressure_to_z(p, pc, axis):
    """Darstellungskoordinate z(p). log: nur p>0 (sonst NaN, kein Clipping); symlog/asinh/linear: alle p."""
    p = np.asarray(p, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        if axis == 'log':
            return np.log10(p / pc)
        if axis == 'asinh':
            return np.arcsinh(p / pc)
        if axis == 'linear':
            return p / pc
        if axis == 'symlog':
            return np.sign(p) * np.log10(1. + np.abs(p) / (SYMLOG_LINTHRESH * pc))
    raise ValueError(f'pressure_axis muss eines von {PRESSURE_AXES} sein')


def z_to_pressure(z, pc, axis):
    """Umkehrung von pressure_to_z: Darstellungskoordinate z -> Rohdruck [Pa]."""
    z = np.asarray(z, dtype=float)
    with np.errstate(over='ignore', invalid='ignore'):
        if axis == 'log':
            return pc * 10. ** z
        if axis == 'asinh':
            return pc * np.sinh(z)
        if axis == 'linear':
            return pc * z
        if axis == 'symlog':
            return np.sign(z) * (10. ** np.abs(z) - 1.) * SYMLOG_LINTHRESH * pc
    raise ValueError(f'pressure_axis muss eines von {PRESSURE_AXES} sein')


def effective_pressure_axis(requested, arrays, *, auto=True):
    """Wirksame Druckachse: 'log' wird auf 'symlog' umgestellt, wenn eine Flaeche p<=0 enthaelt.

    arrays: Iterable von (...,3)-Arrays (T,v,p). Rueckgabe (achse, hinweis oder None).
    """
    if requested not in PRESSURE_AXES:
        raise ValueError(f'pressure_axis muss eines von {PRESSURE_AXES} sein')
    if requested != 'log' or not auto:
        return requested, None
    for arr in arrays:
        a = np.asarray(arr, dtype=float)
        if a.ndim and a.shape[-1] == 3:
            fin = a[np.isfinite(a).all(axis=-1)]
            if fin.size and np.any(fin[:, 2] <= 0):
                return 'symlog', ("Druckachse automatisch von 'log' auf 'symlog' umgestellt: homogene "
                                  "Flaechen enthalten Druecke p <= 0, die auf einer Log-Achse nicht "
                                  "darstellbar waeren (config: diagram_settings.pvt.auto_pressure_axis=false "
                                  "schaltet das ab).")
    return requested, None


def pressure_ticks(pmin, pmax, pc, axis, max_ticks=9):
    """Tickpositionen z und Beschriftungen in ECHTEN Druckwerten [Pa] (mit SI-Praefix)."""
    def fmt(x):
        a = abs(x)
        for scale, unit in ((1e9, 'GPa'), (1e6, 'MPa'), (1e3, 'kPa')):
            if a >= scale:
                return f'{x / scale:g} {unit}'
        return f'{x:g} Pa'
    lo, hi = float(pmin), float(pmax)

    def candidates(mantissas):
        cands = {0.0}
        for e in range(2, 12):
            for m in mantissas:
                for val in (m * 10. ** e, -m * 10. ** e):
                    if lo <= val <= hi:
                        cands.add(val)
        return sorted(c for c in cands if lo <= c <= hi)

    vals = candidates((1.,))
    # Schmale Druckbereiche (z. B. nur Zweiphasenfläche): feinere 1-2-5-Teilung.
    if len([c for c in vals if c != 0.0]) < 3:
        vals = candidates((1., 2., 5.))
    if axis == 'log':
        vals = [c for c in vals if c > 0]
    while len(vals) > max_ticks:
        keep_zero = 0.0 in vals
        rest = [c for c in vals if c != 0.0][::2]
        vals = sorted(rest + ([0.0] if keep_zero else []))
    vals_a = np.array(vals, dtype=float)
    return pressure_to_z(vals_a, pc, axis), [fmt(v) for v in vals_a]


# ---------------------------------------------------------------------------
# P7: Export und Kennzeichnung
# ---------------------------------------------------------------------------

EXPORT_COLUMNS = ('surface', 'T_K', 'v_m3_per_mol', 'p_Pa', 'source', 'stability', 'display',
                  'filter', 'eval_status')


def export_rows(fields):
    """Zeilen (Tupel) aller Knoten; `source` ist 'HEOS' oder 'Ersatz' (Ersatzpunkte nie als HEOS)."""
    rows = []
    for fld in fields:
        n_T, n_s = fld.shape
        for i in range(n_T):
            for j in range(n_s):
                rows.append((fld.name, float(fld.T[i]), float(fld.v[i, j]), float(fld.p_raw[i, j]),
                             SRC_LABELS[int(fld.source[i, j])], STAB_LABELS[int(fld.stability[i, j])],
                             bool(fld.display[i, j]), FILTER_LABELS[int(fld.filter_code[i, j])],
                             EVAL_LABELS[int(fld.eval_status[i, j])]))
    return rows


def surface_as_field(name, points, T=None, source=SRC_HEOS):
    """Hilfsfeld fuer bestehende (n_T,n_v,3)-HEOS-Arrays (metastabile Flaechen) zum Export."""
    pts = np.asarray(points, dtype=float)
    n_T, n_s = pts.shape[:2]
    fld = _empty_field(name, pts[:, 0, 0] if T is None else T, n_s, source)
    fld.v[:], fld.p_raw[:] = pts[..., 1], pts[..., 2]
    finite = np.isfinite(pts).all(axis=2)
    fld.eval_status[finite] = EVAL_OK
    fld.display[:] = finite
    fld.filter_code[finite] = FILTER_NONE
    fld.filter_code[~finite] = FILTER_NOT_EVALUATED
    fld.row_status[:] = np.where(finite.all(axis=1), ROW_OK, ROW_NO_SPINODAL)
    fld.stability[finite] = STAB_METASTABLE
    return fld


def write_csv(path, fields):
    import csv
    rows = export_rows(fields)
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        w = csv.writer(handle)
        w.writerow(EXPORT_COLUMNS)
        w.writerows(rows)
    return len(rows)


def hover_label(source, T, v, p):
    """Hover-Text: Ersatzpunkte tragen immer das Ersatz-Label, Schnittpunkte das Interpolations-Label."""
    base = f'T={T:.6g} K, v={v:.6g} m³/mol, p={p:.6g} Pa'
    if source == SRC_SUBSTITUTE:
        return f'{SUBSTITUTE_HOVER_LABEL} | {base}'
    if source == SRC_CLIPPED:
        return f'{CLIPPED_HOVER_LABEL} | {base}'
    return f'HEOS | {base}'


def homogeneous_fields_for_export(daten):
    """Alle homogenen Knotenfelder eines Berechnungsergebnisses (HEOS-Ranflaechen + Ersatzflaeche)."""
    fields = []
    for key in ('metastable_liquid_surface', 'metastable_vapor_surface'):
        arr = daten.get(key)
        if arr is not None and np.asarray(arr).ndim == 3:
            fields.append(surface_as_field(key, arr))
    raw = daten.get('homogeneous_raw_field')
    if raw is not None:
        fields.append(raw)
    sub = daten.get('replacement_field')
    if sub is not None:
        fields.append(sub)
    return fields


def export_homogeneous_csv(daten, path):
    """CSV-Export; Spalte `source` = HEOS | Ersatz. Rueckgabe: Zeilen- und Herkunftsstatistik."""
    fields = homogeneous_fields_for_export(daten)
    n = write_csv(path, fields)
    by_source = {label: 0 for label in SRC_LABELS.values()}
    for f in fields:
        for k, label in SRC_LABELS.items():
            by_source[label] += int(np.count_nonzero(f.source == k))
    return {'rows': n, 'by_source': by_source, 'source_sum': int(sum(by_source.values()))}


HOMOGENEOUS_CHOICE = ('unstable_surface', 'replacement_surface')


def resolve_homogeneous_choice(checked, changed):
    """HEOS-Rohfortsetzung und Ersatzflaeche nie gleichzeitig: wer zuletzt eingeschaltet wurde, gewinnt."""
    out = dict(checked)
    if out.get(changed) and changed in HOMOGENEOUS_CHOICE:
        for key in HOMOGENEOUS_CHOICE:
            if key != changed:
                out[key] = False
    return out
