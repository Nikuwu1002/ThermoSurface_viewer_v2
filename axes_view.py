"""
Datei: axes_view.py
Projekt: ThermoSurface
Kategorie: Visualisierung / Achsen

Beschreibung:
Zeichnet den physikalisch beschrifteten 3D-Achsenrahmen für die Diagramm-Viewports.

Aufgaben:
- Ermittelt Geometriegrenzen und berechnet Achsenticks und Titel.
- Passt Kamera und Standardansicht an die dargestellte Box an.

Abhängigkeiten:
- Projektmodule (direkt): diagram_axes.py
- Bibliotheken: numpy, pyvista

Verwendet von:
- Projektmodule: diagram_viewport.py, gui.py, visualization.py
"""
import logging

import numpy as np
import pyvista as pv

import diagram_axes as axes

LOGGER = logging.getLogger("thermosurface.axes")

FRAME = "ts_axes_frame"
TICKS = "ts_axes_ticks"
TITLES = "ts_axes_titles"
LEGACY_KEYS = ("pressure_ticks",)

_EDGES = ((0, 1), (2, 3), (4, 5), (6, 7), (0, 2), (1, 3), (4, 6), (5, 7), (0, 4), (1, 5), (2, 6), (3, 7))
DEFAULT_VIEW_VECTOR = (0.0, -1.0, 0.25)
DEFAULT_VIEW_UP = (0.0, 0.0, 1.0)


def visible_bounds(actors, include_hidden=False):
    """Gemeinsame Bounds der Geometrie-Actors (ohne Achsenbeschriftung); None, wenn keine endliche Geometrie existiert."""
    boxes = []
    for key, objects in (actors or {}).items():
        if key in LEGACY_KEYS or objects is None:
            continue
        for actor in (objects if isinstance(objects, (list, tuple)) else [objects]):
            if actor is None or not (include_hidden or actor.visibility):
                continue
            box = np.asarray(actor.GetBounds(), dtype=float)
            if box.shape == (6,) and np.isfinite(box).all() and np.all(box[0::2] <= box[1::2]):
                boxes.append(box)
    if not boxes:
        return None
    boxes = np.array(boxes)
    return tuple(np.column_stack((boxes[:, 0::2].min(axis=0), boxes[:, 1::2].max(axis=0))).ravel())


def auto_bounds(diagram_id, actors):
    """Automatische Sichtgrenzen: p-v-T folgt wie bisher den sichtbaren Ebenen, alle anderen der gesamten Geometrie."""
    return visible_bounds(actors, include_hidden=diagram_id != "pvt")


def _frame_mesh(box):
    xmin, xmax, ymin, ymax, zmin, zmax = box
    corners = np.array([(x, y, z) for z in (zmin, zmax) for y in (ymin, ymax) for x in (xmin, xmax)], dtype=float)
    lines = np.array([[2, a, b] for a, b in _EDGES], dtype=np.int64).ravel()
    return pv.PolyData(corners, lines=lines)


def _camera_vectors(plotter):
    try:
        camera = plotter.camera
        position = np.array(camera.position, dtype=float)
        focal = np.array(camera.focal_point, dtype=float)
        up = np.array(camera.up, dtype=float)
    except Exception:
        return None
    return position, focal, up


def choose_edges(box, camera=None):
    """x/y-Kanten liegen auf der unteren Boxseite; die vordere Ecke und z-Kante folgen der Kamera.

    Rueckgabe {achse: (fixe Koordinaten (3,), Auswaertsrichtung (3,))}; die Koordinate der Achse selbst wird variiert.
    """
    lo, hi = np.array(box[0::2], dtype=float), np.array(box[1::2], dtype=float)
    center = 0.5 * (lo + hi)
    size = np.maximum(hi - lo, 1e-12)
    if camera is None:
        position = center + np.array([1.0, -1.0, 1.0]) * size * 3.0
        right = np.array([1.0, 1.0, 0.0])
    else:
        position, focal, up = camera
        forward = focal - position
        right = np.cross(forward, up)
        if not np.linalg.norm(right) > 0:
            right = np.array([1.0, 0.0, 0.0])

    near = np.where(position > center, hi, lo)           # vordere Ecke je Achse
    far = np.where(position > center, lo, hi)
    sign_near = np.where(position > center, 1.0, -1.0)    # Auswaertsrichtung der vorderen Seite
    z_bottom = lo[2]
    z_out = -1.0

    edges = {
        0: (np.array([np.nan, near[1], z_bottom]), np.array([0.0, sign_near[1], 0.5 * z_out])),
        1: (np.array([near[0], np.nan, z_bottom]), np.array([sign_near[0], 0.0, 0.5 * z_out])),
    }
    candidates = [(near[0], far[1]), (far[0], near[1])]
    screen = [float(np.dot(np.array([x, y, center[2]]) - position, right)) for x, y in candidates]
    x_e, y_e = candidates[int(np.argmin(screen))]
    edges[2] = (np.array([x_e, y_e, np.nan]),
                np.array([1.0 if x_e == hi[0] else -1.0, 1.0 if y_e == hi[1] else -1.0, 0.0]))
    return edges


def _tick_points(diagram_id, daten, box, edges):
    """Tickpositionen (Darstellungskoordinaten) und Labeltexte je Achse."""
    physical = axes.display_bounds_to_physical(diagram_id, daten, box)
    size = np.maximum(np.array(box[1::2]) - np.array(box[0::2]), 1e-12)
    pad = 0.035 * size
    result = {}
    for index, spec in enumerate(axes.axes_for(diagram_id)):
        lo, hi = physical[index]
        kind = axes.axis_scale(diagram_id, index, daten)
        values = axes.physical_ticks(kind, lo, hi, min_abs=axes.tick_floor(diagram_id, index, daten))
        if not values:
            result[index] = (np.empty((0, 3)), [])
            continue
        coords = axes.to_display(diagram_id, daten, index, np.asarray(values) / spec.scale)
        ok = np.isfinite(coords) & (coords >= box[2 * index] - 1e-9 * size[index]) \
            & (coords <= box[2 * index + 1] + 1e-9 * size[index])
        coords = coords[ok]
        labels = [axes.format_tick(v) for v, keep in zip(values, ok) if keep]
        fixed, outward = edges[index]
        points = np.tile(fixed, (len(coords), 1))
        points[:, index] = coords
        points += outward * pad
        result[index] = (points, labels)
    return result


def apply_axes(plotter, diagram_id, daten, actors=None, box=None, camera=None):
    """Zeichnet Rahmen, physikalische Ticks und Titel. `box` (VTK-Koordinaten) oder automatisch aus den Actors.

    Rueckgabe: verwendete Box oder None (nichts zu zeichnen). Die Kamera wird nicht veraendert.
    """
    if box is None:
        box = auto_bounds(diagram_id, actors)
        if box is None:
            return None
    box = tuple(float(x) for x in box)
    try:
        plotter.remove_bounds_axes()
    except Exception:
        LOGGER.debug("Alte Achsen konnten nicht entfernt werden", exc_info=True)
    if actors is not None:
        for key in LEGACY_KEYS:
            if actors.get(key) is not None:
                plotter.remove_actor(key, render=False)
            actors[key] = None

    edges = choose_edges(box, camera if camera is not None else _camera_vectors(plotter))
    plotter.add_mesh(_frame_mesh(box), color="dimgray", line_width=1, name=FRAME, pickable=False,
                     reset_camera=False, render=False)

    ticks = _tick_points(diagram_id, daten, box, edges)
    all_points = [p for p, _ in ticks.values() if len(p)]
    all_labels = [t for _, labels in ticks.values() for t in labels]
    if all_points:
        plotter.add_point_labels(np.vstack(all_points), all_labels, show_points=False, shape=None, font_size=9,
                                 text_color="black", always_visible=True, name=TICKS, render=False,
                                 pickable=False, reset_camera=False)
    else:
        plotter.remove_actor(TICKS, render=False)

    size = np.maximum(np.array(box[1::2]) - np.array(box[0::2]), 1e-12)
    title_points, title_texts = [], []
    for index in range(3):
        fixed, outward = edges[index]
        point = np.array(fixed, dtype=float)
        point[index] = 0.5 * (box[2 * index] + box[2 * index + 1])
        point = point + outward * 0.11 * size
        title_points.append(point)
        title_texts.append(axes.axis_title(diagram_id, index, daten))
    plotter.add_point_labels(np.array(title_points), title_texts, show_points=False, shape=None, font_size=11,
                             bold=True, text_color="black", always_visible=True, name=TITLES, render=False,
                             pickable=False, reset_camera=False)
    return box


def fit_camera(plotter, box):
    """Kamera auf Mittelpunkt und Ausdehnung der Box einpassen; Blickrichtung, view_up und Projektion bleiben."""
    plotter.reset_camera(render=False, bounds=tuple(float(x) for x in box))


def set_default_view(plotter, box):
    """Setzt die Standardansicht von der -Y-Seite und passt die Kamera an die Geometriebox an."""
    plotter.view_vector(DEFAULT_VIEW_VECTOR, viewup=DEFAULT_VIEW_UP, render=False,
                        bounds=tuple(float(x) for x in box))
