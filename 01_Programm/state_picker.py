
"""
Datei: state_picker.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Geometrische Lokalisierung ausgewählter Punkte
auf den dreidimensionalen Zustandsflächen.

Aufgaben:
- Erkennung angeklickter Flächen mittels VTK-CellPicker
- Bestimmung baryzentrischer Koordinaten
- Zuordnung zu ursprünglichen Raster- und Meshpunkten
- Unterscheidung zwischen verschiedenen Flächenarten
- Behandlung geometrisch beschnittener Oberflächen
- Bereitstellung der Trefferinformationen als MeshPick

Abhängigkeiten:
- state_store.py: Rasterreferenzen und Flächenklassifizierung
- NumPy: Geometrische Berechnungen
- VTK und SciPy: Zelllokalisierung und Randprüfung

Verwendet von:
- gui.py: Auswahl von Zustandspunkten
- state_resolver.py: Weiterverarbeitung der MeshPick-Daten
"""

import logging
import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from state_store import (CRITICAL, HOMOGENEOUS, SATURATION, SINGLE_PHASE, SUBSTITUTE, TWO_PHASE, NodeRef, format_ref)

LOGGER = logging.getLogger("thermosurface.picker")

CLIPPED_CELL_ARRAY = "ts_orig_cell"   # Zellarray geclippter Netze (viewport_clipping.ORIG_CELL)

EXACT_NODE_WEIGHT = 1.0 - 1e-6   # Gewicht ab dem ein Treffer als exakter Gitterknoten gilt
PICK_EDGE_EPSILON = 1e-6
PICK_FALLBACK_RADIUS_PX = 7

# Actor-Schluessel, die echte Zustandsflaechen darstellen (Kurven, kritischer Punkt, Ticks sind nicht pickbar).
PICKABLE_ACTORS = ("surface", "zweiphasengebiet", "metastable_liquid_surface", "metastable_vapor_surface",
                   "unstable_surface", "replacement_surface")

# source_class
CLASS_HEOS = "heos"
CLASS_TWO_PHASE = "two_phase"
CLASS_HOMOGENEOUS = "homogeneous"
CLASS_SUBSTITUTE = "synthetic_substitute"


@dataclass
class MeshPick:
    diagram_id: str
    actor_key: str
    cell_id: int
    source_point_ids: list                  # 3 x NodeRef | None
    weights: tuple                          # barycentrische Gewichte (Summe 1)
    display_xyz: tuple
    physical_hint: dict = field(default_factory=dict)
    source_class: str = CLASS_HEOS
    fluid: str = ""
    sub_index: Optional[int] = None
    message: str = ""
    clipped: bool = False                   # True: Treffer auf einem am Ansichtsbereich beschnittenen Netz

    @property
    def is_exact_node(self):
        return max(self.weights) >= EXACT_NODE_WEIGHT

    @property
    def exact_ref(self):
        if not self.is_exact_node:
            return None
        return self.source_point_ids[int(np.argmax(self.weights))]

    @property
    def fully_located(self):
        return all(ref is not None for ref in self.source_point_ids)

    def point_ids(self):
        return [format_ref(r) if r is not None else None for r in self.source_point_ids]


class PickError(ValueError):
    """Treffer nicht auswertbar (degeneriertes Dreieck, nicht pickbarer Actor, unbekannte Quellknoten)."""


def barycentric_weights(triangle_xyz, position):
    """Baryzentrische Koordinaten von `position` bezueglich des Dreiecks; auf die Ebene projiziert."""
    a, b, c = (np.asarray(p, dtype=float) for p in triangle_xyz)
    p = np.asarray(position, dtype=float)
    v0, v1, v2 = b - a, c - a, p - a
    d00, d01, d11 = v0 @ v0, v0 @ v1, v1 @ v1
    d20, d21 = v2 @ v0, v2 @ v1
    denom = d00 * d11 - d01 * d01
    if not np.isfinite(denom) or abs(denom) <= 1e-30 * max(d00 * d11, 1e-300):
        raise PickError("Degeneriertes Dreieck: keine eindeutigen Gewichte")
    w1 = (d11 * d20 - d01 * d21) / denom
    w2 = (d00 * d21 - d01 * d20) / denom
    w = np.array([1.0 - w1 - w2, w1, w2])
    if (w < -PICK_EDGE_EPSILON).any():
        raise PickError("Position liegt ausserhalb des Dreiecks")
    w = np.clip(w, 0.0, None)
    return tuple(float(x) for x in w / w.sum())


def classify(index, refs):
    """Quellklasse des Treffers aus den Flaechenarten der Eckknoten."""
    kinds = [index.surface(r.surface).kind for r in refs if r is not None]
    main = [k for k in kinds if k not in (SATURATION, CRITICAL)]
    kind = main[0] if main else (kinds[0] if kinds else SINGLE_PHASE)
    return {TWO_PHASE: CLASS_TWO_PHASE, HOMOGENEOUS: CLASS_HOMOGENEOUS, SUBSTITUTE: CLASS_SUBSTITUTE}.get(kind, CLASS_HEOS)


def pick_in_cell(index, actor_key, sub_index, cell_xyz, position, cell_id=-1):
    """Lokalisiert `position` im Dreieck `cell_xyz` (3x3 Darstellungskoordinaten) und ordnet Quellknoten zu."""
    if actor_key not in PICKABLE_ACTORS:
        raise PickError(f"Actor {actor_key!r} ist keine pickbare Zustandsflaeche")
    cell_xyz = np.asarray(cell_xyz, dtype=float)
    if cell_xyz.shape != (3, 3):
        raise PickError("Nur Dreieckszellen sind pickbar")
    weights = barycentric_weights(cell_xyz, position)
    located = index.locate(actor_key, sub_index, cell_xyz)
    refs = [ref for ref, _ in located]
    adapter = index.adapter
    hint = {
        "axes": adapter.axes,
        "native_interpolated": tuple(float(x) for x in adapter.from_display(np.asarray(position, dtype=float), index.daten)),
        "vertex_native": [index.native(r).tolist() if r is not None else None for r in refs],
        "vertex_distance": [float(d) for _, d in located],
        "vertex_kinds": [index.surface(r.surface).kind if r is not None else None for r in refs],
        "vertex_phase": [index.surface(r.surface).phase if r is not None else None for r in refs],
        "vertex_source": [index.surface(r.surface).node_source(r.index) if r is not None else None for r in refs],
        "vertex_aux_v": [index.surface(r.surface).node_aux("v", r.index) if r is not None else None for r in refs],
    }
    pick = MeshPick(index.diagram_id, actor_key, int(cell_id), refs, weights, tuple(float(x) for x in position),
                    physical_hint=hint, fluid=index.fluid, sub_index=sub_index)
    if all(r is not None for r in refs):
        pick.source_class = classify(index, refs)
    else:
        pick.message = "Mindestens ein Dreiecksknoten ist keinem Quellknoten zuordenbar"
    return pick


# ---------------------------------------------------------------- VTK-Anbindung

def _same(a, b):
    if a is b:
        return True
    try:
        return a.GetAddressAsString("") == b.GetAddressAsString("")
    except Exception:
        return False


def find_actor(actors, picked_actor, picked_mapper=None):
    """(actor_key, sub_index) des getroffenen Actors innerhalb des Actor-Dictionaries eines Viewports."""
    for key, objects in actors.items():
        if objects is None:
            continue
        listed = isinstance(objects, (list, tuple))
        for i, obj in enumerate(objects if listed else [objects]):
            if obj is None:
                continue
            if _same(obj, picked_actor) or (picked_mapper is not None and _same(obj.GetMapper(), picked_mapper)):
                return key, (i if listed else None)
    return None, None


def _pick_at_viewport(viewport, index, x, y):
    from vtkmodules.vtkRenderingCore import vtkCellPicker

    plotter = viewport.plotter
    picker = vtkCellPicker()
    picker.SetTolerance(0.0)
    if not picker.Pick(x, y, 0, plotter.renderer):
        raise PickError("Kein Treffer auf einer Flaeche")
    key, sub_index = find_actor(viewport.actors, picker.GetActor(), picker.GetMapper())
    if key is None or key not in PICKABLE_ACTORS:
        raise PickError("Treffer liegt nicht auf einer pickbaren Zustandsflaeche (Kurve/Punkt/Beschriftung)")
    data = picker.GetMapper().GetInput()
    cell_id = picker.GetCellId()
    cell = data.GetCell(cell_id)
    if cell is None or cell.GetNumberOfPoints() != 3:
        raise PickError("Getroffene Zelle ist kein Dreieck")
    original_cells = data.GetCellData().GetArray(CLIPPED_CELL_ARRAY)
    if original_cells is not None:
        # Beschnittenes Netz (TS-VP05): Schnittpunkte sind Renderer-Interpolationen und keine HEOS-Knoten. Die
        # Zustandsrekonstruktion geht daher vom unveraenderten Ursprungsdreieck aus (die Pickposition liegt darin).
        triangle = viewport.clip_registry.base_triangle(picker.GetActor(), int(original_cells.GetValue(cell_id)))
        if triangle is None:
            raise PickError("Beschnittene Zelle ist keinem Ursprungsdreieck zuordenbar")
        xyz = [tuple(float(c) for c in point) for point in triangle]
    else:
        xyz = [data.GetPoint(cell.GetPointId(k)) for k in range(3)]
    try:
        pick = pick_in_cell(index, key, sub_index, xyz, picker.GetPickPosition(), cell_id)
    except PickError as exc:
        if str(exc) == "Position liegt ausserhalb des Dreiecks":
            exc.hit = (key, sub_index, int(cell_id))
        raise
    pick.clipped = original_cells is not None
    return pick


def _nearby_screen_offsets():
    directions = ((1.0, 0.0), (0.7071, 0.7071), (0.0, 1.0), (-0.7071, 0.7071),
                  (-1.0, 0.0), (-0.7071, -0.7071), (0.0, -1.0), (0.7071, -0.7071))
    offsets = set()
    for radius in (2, 4, 6, PICK_FALLBACK_RADIUS_PX):
        for dx, dy in directions:
            offset = round(dx * radius), round(dy * radius)
            if math.hypot(*offset) <= PICK_FALLBACK_RADIUS_PX:
                offsets.add(offset)
    yield from sorted(offsets, key=lambda offset: (offset[0] ** 2 + offset[1] ** 2, offset))


def _outside_surface_hull(viewport, x, y):
    """Ein aeusserer Klick ist nur ausserhalb der konvexen Projektion aller sichtbaren Flaechen sicher."""
    from scipy.spatial import ConvexHull, QhullError
    from vtkmodules.util.numpy_support import vtk_to_numpy

    renderer = viewport.plotter.renderer
    camera = renderer.GetActiveCamera()
    transform = camera.GetCompositeProjectionTransformMatrix(renderer.GetTiledAspectRatio(), 0, 1)
    projection = np.array([[transform.GetElement(i, j) for j in range(4)] for i in range(4)])
    origin = np.asarray(renderer.GetOrigin(), dtype=float)
    size = np.asarray(renderer.GetSize(), dtype=float)
    projected = []
    for key in PICKABLE_ACTORS:
        objects = viewport.actors.get(key)
        if objects is None:
            continue
        for actor in objects if isinstance(objects, (tuple, list)) else (objects,):
            if actor is None or not actor.GetVisibility():
                continue
            data = actor.GetMapper().GetInput()
            points = data.GetPoints() if data is not None else None
            if points is None or not points.GetNumberOfPoints():
                continue
            xyz = vtk_to_numpy(points.GetData())
            matrix = actor.GetMatrix()
            model = np.array([[matrix.GetElement(i, j) for j in range(4)] for i in range(4)])
            homogeneous = np.column_stack((xyz, np.ones(len(xyz))))
            clip = homogeneous @ (projection @ model).T
            if not np.isfinite(clip).all() or np.any(clip[:, 3] <= 0):
                LOGGER.debug("point_pick fallback rejected: projection is not finite or crosses the camera")
                return False
            projected.append(origin + (clip[:, :2] / clip[:, 3, None] + 1) * size / 2)
    if not projected:
        return False
    try:
        hull = ConvexHull(np.concatenate(projected))
    except QhullError:
        LOGGER.debug("point_pick fallback rejected: degenerate projected surface")
        return False
    return bool(np.any(hull.equations[:, :2] @ np.array((x, y)) + hull.equations[:, 2] > 1e-6))


def _screen_depth_matches(viewport, x, y, xyz):
    renderer = viewport.plotter.renderer
    renderer.SetWorldPoint(*xyz, 1.0)
    renderer.WorldToDisplay()
    depth = renderer.GetDisplayPoint()[2]
    return abs(float(renderer.GetZ(int(x), int(y))) - depth) <= 1e-5


def _blocked_by_nonpickable_surface(viewport, x, y, xyz):
    from vtkmodules.vtkRenderingCore import vtkCellPicker

    renderer = viewport.plotter.renderer
    picker = vtkCellPicker()
    picker.SetTolerance(0.0)
    picker.PickFromListOn()
    candidates = []
    try:
        for actor in renderer.GetActors():
            if not actor.GetVisibility() or actor.GetPickable() or actor.GetProperty().GetOpacity() < 1.0:
                continue
            mapper = actor.GetMapper()
            data = mapper.GetInput() if mapper is not None else None
            if data is None or not data.GetNumberOfCells():
                continue
            if data.IsA("vtkPolyData"):
                has_surface = data.GetNumberOfPolys() > 0 or data.GetNumberOfStrips() > 0
            else:
                has_surface = any(data.GetCell(i).GetCellDimension() >= 2 for i in range(data.GetNumberOfCells()))
            if not has_surface:
                continue
            actor.PickableOn()
            candidates.append(actor)
            picker.AddPickList(actor)
        if not candidates or not picker.Pick(x, y, 0, renderer):
            return False
        point = np.asarray(xyz, dtype=float)
        blocker = np.asarray(picker.GetPickPosition(), dtype=float)
        direction = np.asarray(renderer.GetActiveCamera().GetDirectionOfProjection(), dtype=float)
        tolerance = 1e-10 * max(1.0, np.linalg.norm(point), np.linalg.norm(blocker))
        return float((blocker - point) @ direction) < -tolerance
    finally:
        for actor in candidates:
            actor.PickableOff()


def pick_from_viewport(viewport, index, x, y):
    """Pickt exakt und verwendet nur bei geometrisch sicherer Randnaehe benachbarte Strahlen."""
    try:
        pick = _pick_at_viewport(viewport, index, x, y)
        if _blocked_by_nonpickable_surface(viewport, x, y, pick.display_xyz):
            raise PickError("Treffer liegt hinter sichtbarer Geometrie")
        return pick
    except PickError as original_error:
        hit = getattr(original_error, "hit", None)
        reason = str(original_error)
        if reason == "Kein Treffer auf einer Flaeche":
            renderer = viewport.plotter.renderer
            width, height = renderer.GetRenderWindow().GetSize()
            if not (0 <= x < width and 0 <= y < height) or renderer.GetZ(int(x), int(y)) < 1.0:
                raise
            if not _outside_surface_hull(viewport, x, y):
                raise
        elif reason != "Position liegt ausserhalb des Dreiecks" or hit is None:
            raise

    for dx, dy in _nearby_screen_offsets():
        try:
            pick = _pick_at_viewport(viewport, index, x + dx, y + dy)
        except PickError:
            continue
        if hit is not None and (pick.actor_key, pick.sub_index, pick.cell_id) != hit:
            continue
        if hit is None and (not pick.fully_located or not _screen_depth_matches(viewport, x + dx, y + dy,
                                                                                pick.display_xyz)):
            continue
        if _blocked_by_nonpickable_surface(viewport, x + dx, y + dy, pick.display_xyz):
            continue
        pick.message = f"screen-near fallback ({dx:+d}, {dy:+d} px)"
        LOGGER.debug("point_pick fallback_hit diagram=%s offset=(%d,%d) radius=%d",
                     pick.diagram_id, dx, dy, PICK_FALLBACK_RADIUS_PX)
        return pick

    LOGGER.debug("point_pick miss diagram=%s radius=%d px", index.diagram_id, PICK_FALLBACK_RADIUS_PX)
    raise PickError(f"Keine pickbare Oberflaeche innerhalb von {PICK_FALLBACK_RADIUS_PX} px")
