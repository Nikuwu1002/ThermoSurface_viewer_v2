"""
Datei: geometry_clip.py
Projekt: ThermoSurface
Kategorie: Visualisierung / Geometrie

Beschreibung:
Schneidet vorliegende 3D-Geometrien an einem achsenparallelen Quader.

Aufgaben:
- Clipt Dreiecksflächen und Polylinien und prüft Punktpositionen.
- Kennzeichnet geometrisch interpolierte Schnittpunkte getrennt von Rohdaten.

Abhängigkeiten:
- Projektmodule (direkt): keine
- Bibliotheken: numpy

Verwendet von:
- Projektmodule: diagram_viewport.py
"""

import numpy as np

EPS = 1e-12


def _box_arrays(box):
    box = np.asarray(box, dtype=float)
    if box.shape != (6,) or not np.isfinite(box).all() or np.any(box[0::2] > box[1::2]):
        raise ValueError("Clipping-Box muss (xmin, xmax, ymin, ymax, zmin, zmax) mit min <= max sein.")
    return box[0::2], box[1::2]


def points_in_box(points, box, rtol=1e-9):
    """Boolesche Maske: Punkte innerhalb der Box (mit relativer Toleranz an den Flaechen)."""
    lo, hi = _box_arrays(box)
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    tol = rtol * np.maximum(hi - lo, 1e-300)
    with np.errstate(invalid="ignore"):
        return np.all((pts >= lo - tol) & (pts <= hi + tol), axis=1)


# ---------------------------------------------------------------- Dreiecke

def clip_triangles(points, triangles, box, point_data=None):
    """Schneidet Dreiecke an der Box.

    points (n,3), triangles (m,3) Indizes, point_data {name: (n,)-Array} (wird linear interpoliert).
    Rueckgabe dict: points (k,3), triangles (t,3), new (k,) bool (Schnittpunkte), orig_cell (t,) Index des
    Ursprungsdreiecks, point_data {name: (k,)}. Kein Ergebnis -> leere Arrays.
    """
    lo, hi = _box_arrays(box)
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    triangles = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    names = list(point_data or {})
    values = {name: np.asarray(point_data[name], dtype=float) for name in names}
    P = points[triangles]                                                  # (m,3,3)
    A = {name: values[name][triangles] for name in names}
    I = triangles.copy()                                                   # Originalindex je Ecke, -1 = neuer Schnittpunkt
    O = np.arange(len(triangles), dtype=np.int64)

    finite = np.isfinite(P).all(axis=(1, 2))
    P, I, O = P[finite], I[finite], O[finite]
    A = {k: v[finite] for k, v in A.items()}

    for axis in range(3):
        for sign, bound in ((1.0, lo[axis]), (-1.0, hi[axis])):
            P, A, I, O = _clip_plane(P, A, I, O, axis, sign, bound, EPS * max(hi[axis] - lo[axis], 1e-300))
            if not len(P):
                break
        if not len(P):
            break

    if len(P):
        edge1, edge2 = P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]
        area2 = np.linalg.norm(np.cross(edge1, edge2), axis=1)
        keep = area2 > 0.0
        P, I, O = P[keep], I[keep], O[keep]
        A = {k: v[keep] for k, v in A.items()}

    if not len(P):
        return {"points": np.empty((0, 3)), "triangles": np.empty((0, 3), dtype=np.int64),
                "new": np.empty(0, dtype=bool), "orig_cell": np.empty(0, dtype=np.int64),
                "point_data": {name: np.empty(0) for name in names}}

    # Unveraenderte Originalknoten behalten ihren Index (schnell); nur die wenigen Schnittpunkte werden
    # (exakt gleiche Koordinaten) untereinander zusammengefuehrt.
    flat_I = I.reshape(-1)
    original = flat_I >= 0
    used = np.unique(flat_I[original])
    remap = np.full(len(points), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    new_points = P.reshape(-1, 3)[~original]
    if len(new_points):
        unique_new, first_new, inverse_new = np.unique(new_points, axis=0, return_index=True, return_inverse=True)
        inverse_new = inverse_new.reshape(-1)
    else:
        unique_new, first_new, inverse_new = np.empty((0, 3)), np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    flat_tri = np.empty(len(flat_I), dtype=np.int64)
    flat_tri[original] = remap[flat_I[original]]
    flat_tri[~original] = len(used) + inverse_new
    out_data = {name: np.concatenate((values[name][used], A[name].reshape(-1)[~original][first_new]))
                for name in names}
    return {"points": np.vstack((points[used], unique_new)), "triangles": flat_tri.reshape(-1, 3),
            "new": np.r_[np.zeros(len(used), bool), np.ones(len(unique_new), bool)], "orig_cell": O,
            "point_data": out_data}


def _clip_plane(P, A, I, O, axis, sign, bound, tol):
    """Behaelt den Halbraum sign*(x[axis]-bound) >= 0; liefert (P, A, I, O) mit neuen Dreiecken."""
    d = sign * (P[:, :, axis] - bound)
    inside = d >= -tol
    count = inside.sum(axis=1)

    keep_all = count == 3
    one = count == 1
    two = count == 2
    out_P, out_I, out_O = [P[keep_all]], [I[keep_all]], [O[keep_all]]
    out_A = {k: [v[keep_all]] for k, v in A.items()}

    def lerp(i_idx, j_idx, rows):
        """Schnittpunkt auf Kante (i drinnen, j draussen) fuer die Dreiecke `rows`."""
        di = d[rows, i_idx]
        dj = d[rows, j_idx]
        t = np.clip(di / (di - dj), 0.0, 1.0)
        pi, pj = P[rows, i_idx], P[rows, j_idx]
        point = pi + t[:, None] * (pj - pi)
        point[:, axis] = bound
        values = {k: v[rows, i_idx] + t * (v[rows, j_idx] - v[rows, i_idx]) for k, v in A.items()}
        return point, values

    if one.any():
        rows = np.flatnonzero(one)
        i = np.argmax(inside[rows], axis=1)                               # einziger innerer Knoten
        j, k = (i + 1) % 3, (i + 2) % 3
        pij, aij = lerp(i, j, rows)
        pik, aik = lerp(i, k, rows)
        none = np.full(len(rows), -1, dtype=np.int64)
        out_P.append(np.stack((P[rows, i], pij, pik), axis=1))
        out_I.append(np.stack((I[rows, i], none, none), axis=1))
        out_O.append(O[rows])
        for name in A:
            out_A[name].append(np.stack((A[name][rows, i], aij[name], aik[name]), axis=1))

    if two.any():
        rows = np.flatnonzero(two)
        k = np.argmax(~inside[rows], axis=1)                              # einziger aeusserer Knoten
        i, j = (k + 1) % 3, (k + 2) % 3
        pjk, ajk = lerp(j, k, rows)
        pki, aki = lerp(i, k, rows)
        none = np.full(len(rows), -1, dtype=np.int64)
        pi_, pj_ = P[rows, i], P[rows, j]
        ii_, ij_ = I[rows, i], I[rows, j]
        out_P.append(np.stack((pi_, pj_, pjk), axis=1))
        out_I.append(np.stack((ii_, ij_, none), axis=1))
        out_P.append(np.stack((pi_, pjk, pki), axis=1))
        out_I.append(np.stack((ii_, none, none), axis=1))
        out_O.extend([O[rows], O[rows]])
        for name in A:
            ai_, aj_ = A[name][rows, i], A[name][rows, j]
            out_A[name].append(np.stack((ai_, aj_, ajk[name]), axis=1))
            out_A[name].append(np.stack((ai_, ajk[name], aki[name]), axis=1))

    return (np.concatenate(out_P), {k: np.concatenate(v) for k, v in out_A.items()},
            np.concatenate(out_I), np.concatenate(out_O))


# ---------------------------------------------------------------- Kurven

def clip_polylines(points, lines, box):
    """Schneidet Polylinien (Liste von Indexarrays) segmentweise an der Box.

    Rueckgabe dict: points (k,3), lines (Liste von Indexarrays), orig_line (Index der Ursprungslinie je Ergebnislinie).
    Zusammenhaengende, unveraenderte Segmente bleiben eine Polylinie; geschnittene Enden liegen auf der Boxflaeche.
    """
    lo, hi = _box_arrays(box)
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    starts, ends, owner, order = [], [], [], []
    for number, line in enumerate(lines):
        idx = np.asarray(line, dtype=np.int64)
        if len(idx) < 2:
            continue
        starts.append(points[idx[:-1]])
        ends.append(points[idx[1:]])
        owner.append(np.full(len(idx) - 1, number, dtype=np.int64))
        order.append(np.arange(len(idx) - 1, dtype=np.int64))
    empty = {"points": np.empty((0, 3)), "lines": [], "orig_line": np.empty(0, dtype=np.int64)}
    if not starts:
        return empty
    a, b = np.concatenate(starts), np.concatenate(ends)
    owner, order = np.concatenate(owner), np.concatenate(order)
    d = b - a
    t0 = np.zeros(len(a))
    t1 = np.ones(len(a))
    alive = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
    tol = EPS * np.maximum(hi - lo, 1e-300)
    with np.errstate(divide="ignore", invalid="ignore"):
        for axis in range(3):
            for p, q in ((-d[:, axis], a[:, axis] - lo[axis]), (d[:, axis], hi[axis] - a[:, axis])):
                parallel = p == 0.0
                alive &= ~(parallel & (q < -tol[axis]))
                r = np.where(parallel, 0.0, q / np.where(parallel, 1.0, p))
                entering = (p < 0) & ~parallel
                leaving = (p > 0) & ~parallel
                t0 = np.where(entering, np.maximum(t0, r), t0)
                t1 = np.where(leaving, np.minimum(t1, r), t1)
    alive &= (t1 - t0) > EPS
    kept = np.flatnonzero(alive)
    if not len(kept):
        return empty

    pa = a[kept] + t0[kept, None] * d[kept]
    pb = a[kept] + t1[kept, None] * d[kept]
    clipped_a = t0[kept] > 0.0
    clipped_b = t1[kept] < 1.0
    pa[clipped_a] = np.clip(pa[clipped_a], lo, hi)
    pb[clipped_b] = np.clip(pb[clipped_b], lo, hi)
    owner_k, order_k = owner[kept], order[kept]
    joined = np.zeros(len(kept), dtype=bool)
    joined[1:] = ((owner_k[1:] == owner_k[:-1]) & (order_k[1:] == order_k[:-1] + 1)
                  & ~clipped_b[:-1] & ~clipped_a[1:])

    out_points, out_lines, out_owner = [], [], []
    current = None
    for n in range(len(kept)):
        if not joined[n]:
            if current is not None:
                out_lines.append(np.asarray(current, dtype=np.int64))
            current = [len(out_points)]
            out_points.append(pa[n])
            out_owner.append(owner_k[n])
        out_points.append(pb[n])
        current.append(len(out_points) - 1)
    out_lines.append(np.asarray(current, dtype=np.int64))
    return {"points": np.asarray(out_points), "lines": out_lines, "orig_line": np.asarray(out_owner, dtype=np.int64)}
