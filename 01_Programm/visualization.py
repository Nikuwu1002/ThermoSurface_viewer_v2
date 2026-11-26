"""
Datei: visualization.py
Projekt: ThermoSurface
Kategorie: Visualisierung / 3D-Rendering

Beschreibung:
Rendert vorberechnete Zustandsflächen aller sieben Diagrammtypen mit PyVista.

Aufgaben:
- Erzeugt Flächenmeshes, Isolinien, Phasengrenzen und zusätzliche Actor-Layer.
- Nutzt einheitliche Koordinaten, Achsen und physikalische Hover-Metadaten.

Abhängigkeiten:
- Projektmodule (direkt): axes_view.py, diagram_coords.py, homogeneous_pvt.py,
  hover_data.py, zustandsflaeche.py
- Bibliotheken: numpy, pyvista

Verwendet von:
- Projektmodule: diagram_registry.py, diagram_viewport.py, gui.py
"""

import pyvista as pv
import numpy as np
import logging

from hover_data import attach_hover_metadata, format_mesh_hover_point
from homogeneous_pvt import (node_validity, grid_triangles, drop_degenerate,
                             mesh_topology_report)
from zustandsflaeche import (PRESSURE_AXES, pressure_to_z, z_to_pressure, pressure_ticks, assemble_shared_vertices,
                             surface_triangles, critical_fan, fan_applicable, hover_label, SRC_SUBSTITUTE)
LOGGER = logging.getLogger("thermosurface.visualization")

# Koordinatentransformationen liegen PyVista-frei in diagram_coords.py (gemeinsam mit der Zustandsrekonstruktion).
from diagram_coords import to_xyz, to_xyz_pht, to_xyz_pst, to_xyz_hsp, to_xyz_uvt  # noqa: E402,F401
from diagram_coords import to_xyz_gpt, to_xyz_fvt  # noqa: E402,F401
from axes_view import apply_axes, visible_bounds  # noqa: E402


def _pht_surface_mesh(surface, critical):
    """
    Erstellt ein PyVista-Mesh aus einem
    regulären p-h-T-Oberflächenarray.
    """

    surface = np.asarray(surface, dtype=float)

    if (surface.ndim != 3 or surface.shape[-1] != 3 or min(surface.shape[:2]) < 2):
        raise ValueError("p-h-T-Oberfläche benötigt Form (n_T, n_p, 3).")

    xyz = to_xyz_pht(surface, critical)

    if not np.isfinite(xyz).all():
        raise ValueError("Ungültige p-h-T-Koordinaten.")

    grid = pv.StructuredGrid(
        xyz[..., 0],
        xyz[..., 1],
        xyz[..., 2])
    attach_hover_metadata(grid, surface.reshape((-1, 3), order="F"), "pht")

    return grid.extract_surface().triangulate()
    

def create_surface(daten):

    surface = np.asarray(daten["surface"], dtype=float)

    # Erwartete Form: (Anzahl T, Anzahl v, 3)
    if (surface.ndim != 3 or surface.shape[2] != 3 or min(surface.shape[:2]) < 2):
        raise ValueError("Ungültige Arrayform für surface")

    # Alle Punkte in Darstellungskoordinaten umrechnen
    xyz = to_xyz(surface, daten["critical_point"], daten.get("pressure_axis", "log"))

    # Nur gültige Punkte berücksichtigen
    valid = (
        np.isfinite(surface).all(axis=2)
        & (surface[..., 0] > 0)
        & (surface[..., 1] > 0)
        & (surface[..., 2] > 0)
        & np.isfinite(xyz).all(axis=2))

    # Ungültige Positionen vorläufig durch 0 ersetzen
    safe_xyz = np.where(valid[..., None], xyz, 0.0)

    # Strukturiertes Gitter erzeugen
    grid = pv.StructuredGrid(
        safe_xyz[..., 0],
        safe_xyz[..., 1],
        safe_xyz[..., 2])

    # Nur Zellen mit vier gültigen Eckpunkten behalten
    valid_cells = (
        valid[:-1, :-1]
        & valid[1:, :-1]
        & valid[:-1, 1:]
        & valid[1:, 1:])

    phase_nodes = np.asarray(daten["phase_mask"], dtype=bool)

    if phase_nodes.shape != valid.shape:
        raise ValueError("Ungültige Dimension der Phasenmaske")

    phase_cells = (
        phase_nodes[:-1, :-1]
        & phase_nodes[1:, :-1]
        & phase_nodes[:-1, 1:]
        & phase_nodes[1:, 1:]
        & valid_cells)

    grid.cell_data["two_phase"] = (phase_cells.ravel(order="F").astype(np.uint8))
    ids = np.flatnonzero(valid_cells.ravel(order="F"))

    if ids.size == 0:
        raise ValueError("Keine gültigen Oberflächenzellen")

    return grid.extract_cells(ids)

def _draw_curves(plotter, curves, critical, color, min_z=None, pressure_axis="log"):

    actors = []

    for curve in curves:

        curve = np.asarray(curve, dtype=float)

        # Erwartete Form einer einzelnen Kurve: (N, 3)
        if curve.ndim != 2 or curve.shape[1] != 3:
            raise ValueError("Ungültige Linienform")

        xyz = to_xyz(curve, critical, pressure_axis)

        # Gültige Punkte erkennen
        valid = np.isfinite(xyz).all(axis=1)

        # Zusammenhängende gültige Abschnitte bestimmen
        edges = np.flatnonzero(np.diff(np.r_[False, valid, False].astype(np.int8)))

        for start, end in zip(edges[::2], edges[1::2]):
            if end - start < 2:
                continue
            poly = pv.lines_from_points(xyz[start:end])

            # Optional: Linie an unterer z-Grenze abschneiden
            if min_z is not None:
                poly = poly.clip( normal=(0, 0, 1), origin=(0, 0, float(min_z)), invert=False)

                if poly.n_cells == 0:
                    continue

            actor = plotter.add_mesh(poly, color=color, line_width=3)
            actors.append(actor)

    return actors

def create_isotherm_lines(plotter, daten):
    return _draw_curves(plotter, daten["isothermen"], daten["critical_point"], color="red", pressure_axis=daten.get("pressure_axis", "log"))

def create_isochore_lines(plotter, daten):
    return _draw_curves(plotter, daten["isochoren"], daten["critical_point"], color="green", pressure_axis=daten.get("pressure_axis", "log"))

def create_isobare_lines(plotter, daten):
    return _draw_curves(plotter, daten["isobaren"], daten["critical_point"], color="blue", pressure_axis=daten.get("pressure_axis", "log"))

def add_critical_point(plotter, daten):

    Tc, pc, vc = daten["critical_point"]

    # In die Standardreihenfolge (T, v, p) bringen
    point = np.array([[Tc, vc, pc]])

    xyz = to_xyz(point, daten["critical_point"], daten.get("pressure_axis", "log"))

    actor = plotter.add_mesh(pv.PolyData(xyz), color="black", point_size=12, render_points_as_spheres=True)

    return actor

def add_binodals(plotter, daten):

    critical = daten["critical_point"]

    liquid_actors = _draw_curves(plotter, [daten["binodal_liquid"]], critical, color="blue", pressure_axis=daten.get("pressure_axis", "log"))

    vapor_actors = _draw_curves(plotter, [daten["binodal_vapor"]], critical, color="red", pressure_axis=daten.get("pressure_axis", "log"))

    return liquid_actors + vapor_actors



def add_two_phase_surface(plotter, daten):
    """
    Erstellt die Gleichgewichts-Zweiphasenfläche
    als eigenständiges PyVista-Mesh.

    Eingabe:
        daten["two_phase_surface"]: (n_T, n_Q, 3)
        daten["critical_point"]: (Tc, pc, vc)

    Rückgabe:
        PyVista Actor
    """

    # 1. Thermodynamische Daten einlesen
    surface = np.asarray(daten["two_phase_surface"], dtype=float)

    # Fehlende Sättigungsdaten sind kein fataler Visualisierungsfehler.
    if surface.shape == (0, 0, 3):
        return None

    # 2. Arraystruktur überprüfen
    if (surface.ndim != 3 or surface.shape[2] != 3 or min(surface.shape[:2]) < 2):
        raise ValueError("Zweiphasenfläche muss die Form (n_T, n_Q, 3) haben.")

    if (not np.isfinite(surface).all() or np.any(surface <= 0)):
        raise ValueError("Zweiphasenfläche enthält ungültige Zustandswerte.")

    # 3. Physikalische Daten in 3D-Koordinaten umrechnen
    xyz = to_xyz(surface, daten["critical_point"], daten.get("pressure_axis", "log"))

    if not np.isfinite(xyz).all():
        raise ValueError("Ungültige Koordinaten nach Transformation.")

    # 4. Strukturiertes PyVista-Gitter erstellen
    grid = pv.StructuredGrid(
        xyz[..., 0],
        xyz[..., 1],
        xyz[..., 2])

    # 5. Gittertopologie überprüfen
    n_T, n_Q = surface.shape[:2]

    expected_points = n_T * n_Q
    expected_cells = (n_T - 1) * (n_Q - 1)

    if (grid.n_points != expected_points or grid.n_cells != expected_cells):
        raise ValueError("PyVista-Gitter hat unerwartete Dimensionen.")

    # 6. Mesh im Plotter darstellen
    # Vorhandenes strukturiertes Netz in Dreiecke umwandeln
    attach_hover_metadata(grid, surface.reshape((-1, 3), order="F"), "pvt")
    mesh = grid.extract_surface().triangulate()

    # Kritischen Punkt als einzelnen Knoten ergänzen
    Tc, pc, vc = daten["critical_point"]
    critical_tvp = np.array([[Tc, vc, pc]], dtype=float)

    # Letztes unterkritisches Q-Profil
    last_row = surface[-1, :, :]

    if last_row[0, 0] >= Tc:
        raise ValueError("Die Zweiphasenfläche darf keine " "vollständige kritische Q-Zeile enthalten.")

    # Durch eine obere Druckgrenze abgeschnittene Fläche endet am letzten
    # verbleibenden Q-Profil, nicht am kritischen Punkt.
    if not daten.get("two_phase_reaches_critical", True):
        return plotter.add_mesh(mesh, color="orange", opacity=0.85, show_edges=False, show_scalar_bar=False, lighting=False)

    # Letztes Q-Profil und kritischen Punkt verbinden
    cap_points = np.vstack((last_row, critical_tvp))

    cap_xyz = to_xyz(cap_points, daten["critical_point"], daten.get("pressure_axis", "log"))

    n_Q = len(last_row)

    # Dreieckszellen zum kritischen Punkt
    cap_faces = np.column_stack((
        np.full(n_Q - 1, 3, dtype=np.int64),
        np.arange(n_Q - 1),
        np.arange(1, n_Q),
        np.full(n_Q - 1, n_Q)
    )).ravel()

    cap_mesh = pv.PolyData(cap_xyz, cap_faces)
    attach_hover_metadata(cap_mesh, cap_points, "pvt")

    # Beide Netze verbinden
    mesh = mesh.merge(cap_mesh, merge_points=True)
    actor = plotter.add_mesh(mesh, color="orange", opacity=0.85, show_edges=False, show_scalar_bar=False, lighting=False)

    return actor

def _homogeneous_triangles(values, critical, pressure_axis):
    """Gemeinsame Knoten-/Dreieckslogik fuer Renderer und Diagnose."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 3 or values.shape[-1] != 3:
        raise ValueError("Homogene Flaeche benoetigt (n_T,n_v,3)")
    nr, nc = values.shape[:2]
    valid = node_validity(values, pressure_axis)
    xyz = to_xyz(values, critical, pressure_axis)
    valid &= np.isfinite(xyz).all(axis=2)
    flat = np.where(valid[..., None], xyz, 0.0).reshape((-1, 3))
    triangles, cells = grid_triangles(valid)
    triangles, cells, n_degenerate = drop_degenerate(triangles, cells, flat)
    return valid, flat, triangles, cells, n_degenerate


def homogeneous_surface_mesh(points, critical, pressure_axis="log", critical_cap=False):
    """Vernetzt zusammenhaengende 3- oder 4-Knoten-Zellen ohne NaN-Bruecken.

    4 gueltige Eckpunkte -> 2 Dreiecke; 3 Eckpunkte -> 1 Dreieck.
    Nie Knoten interpolieren oder ueber fehlende Temperaturzeilen verbinden.
    Die Trennspalte der instabilen Flaeche (NaN) wird nicht ueberbrueckt.
    critical_cap: Dreiecksfaecher von der letzten Zeile zum kritischen Punkt (AP20.6);
    nur anwenden, wenn die Annaeherung validiert wurde (daten['critical_closure']).
    """
    values = np.asarray(points, dtype=float)
    if values.ndim != 3 or values.shape[-1] != 3:
        raise ValueError("Homogene Flaeche benoetigt (n_T,n_v,3)")
    if min(values.shape[:2]) < 2:
        return pv.PolyData()
    valid, flat, triangles, _cells, _ = _homogeneous_triangles(values, critical, pressure_axis)
    if critical_cap:
        Tc, pc, vc = critical
        apex = valid.size
        flat = np.vstack((flat, to_xyz(np.array([[Tc, vc, pc]]), critical, pressure_axis)))
        last = np.arange(valid.shape[1]) + (valid.shape[0] - 1) * valid.shape[1]
        both = valid[-1, :-1] & valid[-1, 1:]
        fan = np.column_stack((last[:-1][both], np.full(int(both.sum()), apex), last[1:][both]))
        if len(fan) and np.isfinite(flat[apex]).all():
            triangles = np.vstack((triangles, fan)) if len(triangles) else fan
    if not len(triangles):
        return pv.PolyData()
    used = np.unique(triangles.ravel())
    remap = np.full(len(flat), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    faces = np.column_stack((np.full(len(triangles), 3, dtype=np.int64), remap[triangles]))
    mesh = pv.PolyData(flat[used], faces.ravel())
    if not mesh.n_cells or not np.isfinite(mesh.points).all():
        raise ValueError("Unzulaessige homogene Netzgeometrie")
    return mesh


def homogeneous_mesh_diagnostics(points, critical, pressure_axis="log"):
    """AP19.4/AP20.5: Topologie des Netzes ohne VTK (Komponenten, offene Kanten, Faltungen)."""
    values = np.asarray(points, dtype=float)
    if values.ndim != 3 or min(values.shape[:2]) < 2:
        return {"triangles": 0}
    valid, flat, triangles, cells, n_degenerate = _homogeneous_triangles(
        values, critical, pressure_axis)
    report = mesh_topology_report(triangles, cells, valid, xy=flat[:, :2])
    report["degenerate_removed"] = n_degenerate
    return report


HOMOGENEOUS_STYLE = {
    "metastable_liquid_surface": ("steelblue", 0.8),
    "metastable_vapor_surface": ("mediumseagreen", 0.8),
    "unstable_surface": ("salmon", 0.5),
    "replacement_surface": ("gold", 0.75),
}
_HOMOGENEOUS_LINKS = (
    ("metastable_liquid_surface", "unstable_surface"),
    ("metastable_liquid_surface", "replacement_surface"),
    ("unstable_surface", "metastable_vapor_surface"),
    ("replacement_surface", "metastable_vapor_surface"),
)


def homogeneous_shared_meshes(daten):
    """Alle homogenen Flaechen auf EINER gemeinsamen Vertextabelle (identische Randvertex-Indizes).

    Liefert (meshes, info). Metastabile/instabile HEOS-Flaechen: bestehende Zellregel (>=3 Ecken);
    Ersatzflaeche: feste Diagonale, Dreieck nur bei drei gezeigten Ecken. Jedes Mesh traegt die
    Punktdaten T, v, p_raw, source (0=HEOS, 1=Ersatz) fuer Hover/Export; Rohwerte bleiben unveraendert.
    """
    critical = daten["critical_point"]
    axis = daten.get("pressure_axis", "log")
    grids = {}
    for key in HOMOGENEOUS_STYLE:
        arr = daten.get(key)
        if arr is not None and np.asarray(arr).ndim == 3 and np.asarray(arr).shape[-1] == 3:
            grids[key] = np.asarray(arr, dtype=float)
    if not grids:
        return {}, {}
    shape_rows = {g.shape[0] for g in grids.values()}
    if len(shape_rows) != 1:
        raise ValueError("Homogene Flaechen haben unterschiedliche Temperaturzeilen")
    links = [(a, b) for a, b in _HOMOGENEOUS_LINKS if a in grids and b in grids]
    sv = assemble_shared_vertices(grids, links, critical)
    xyz = to_xyz(sv.vertices, critical, axis)
    finite = np.isfinite(xyz).all(axis=1)
    safe = xyz.copy()
    if finite.any():
        safe[~finite] = xyz[finite][0]
    source = np.zeros(len(sv.vertices), dtype=np.uint8)
    field = daten.get("replacement_field")
    if field is not None and "replacement_surface" in grids:
        interior = (field.source == SRC_SUBSTITUTE) & field.display
        source[sv.ids["replacement_surface"][interior]] = SRC_SUBSTITUTE
    meshes, info = {}, {}
    closure = bool(daten.get("critical_closure", {}).get("applied", False))
    for key, grid in grids.items():
        valid = node_validity(grid, axis)
        valid &= finite[sv.ids[key]]
        ids = sv.ids[key]
        if key == "replacement_surface":
            shown = valid
            if field is not None:
                shown = shown & field.display
            tri = surface_triangles(shown, ids)
            fan_ok, last = (fan_applicable(field.row_status) if field is not None else (False, None))
        else:
            local, _cells = grid_triangles(valid)
            tri = ids.ravel()[local] if len(local) else local
            fan_ok, last = (closure and bool(valid[-1].all()), grid.shape[0] - 1)
            shown = valid
        if fan_ok and sv.apex_id is not None and finite[sv.apex_id]:
            fan = critical_fan(shown[last], ids[last], sv.apex_id)
            if len(fan):
                tri = np.vstack((tri, fan)) if len(tri) else fan
        if len(tri):
            a, b, c = (safe[tri[:, k]] for k in range(3))
            area2 = np.linalg.norm(np.cross(b - a, c - a), axis=1)
            tri = tri[np.isfinite(area2) & (area2 > 0)]
        if not len(tri):
            continue
        faces = np.column_stack((np.full(len(tri), 3, dtype=np.int64), tri)).ravel()
        mesh = pv.PolyData(safe, faces)
        if hasattr(mesh, "point_data"):
            mesh.point_data["T"] = sv.vertices[:, 0]
            mesh.point_data["v"] = sv.vertices[:, 1]
            mesh.point_data["p_raw"] = sv.vertices[:, 2]
            mesh.point_data["source"] = source
        meshes[key] = mesh
        info[key] = {"triangles": int(len(tri)), "fan": bool(fan_ok)}
    info["shared"] = {"mismatches": sv.mismatches, "vertices": int(len(sv.vertices))}
    return meshes, info


def hover_text_for_point(mesh, point_id):
    """Hover-Text zu einem gepickten Punkt (Ersatzpunkte immer mit Ersatz-Label)."""
    text = format_mesh_hover_point(mesh, point_id)
    if text is not None:
        return text
    pd = mesh.point_data
    return hover_label(int(pd["source"][point_id]), float(pd["T"][point_id]),
                       float(pd["v"][point_id]), float(pd["p_raw"][point_id]))


def pick_hover_text(plotter, x, y):
    """Hover: Text zum Punkt unter dem Mauszeiger (nur Flaechen mit Punktdaten T, v, p_raw, source)."""
    from vtkmodules.vtkRenderingCore import vtkCellPicker
    picker = vtkCellPicker()
    picker.SetTolerance(0.005)
    if not picker.Pick(x, y, 0, plotter.renderer):
        return None
    mapper = picker.GetMapper()
    data = mapper.GetInput() if mapper is not None else None
    if data is None:
        return None
    cell = data.GetCell(picker.GetCellId())
    if cell is None:
        return None
    pos = np.array(picker.GetPickPosition())
    best, best_d = None, np.inf
    for k in range(cell.GetNumberOfPoints()):
        pid = cell.GetPointId(k)
        d = float(np.linalg.norm(np.array(data.GetPoint(pid)) - pos))
        if d < best_d:
            best, best_d = pid, d
    text = format_mesh_hover_point(data, best)
    if text is not None:
        return text
    pd = data.GetPointData()
    if any(pd.GetArray(name) is None for name in ("source", "T", "v", "p_raw")):
        return None
    return hover_label(int(pd.GetArray("source").GetValue(best)), pd.GetArray("T").GetValue(best),
                       pd.GetArray("v").GetValue(best), pd.GetArray("p_raw").GetValue(best))


def add_homogeneous_surface(plotter, daten, key, *, color, opacity, meshes=None):
    if meshes is None:
        meshes, _ = homogeneous_shared_meshes(daten)
    mesh = meshes.get(key)
    if mesh is None or not mesh.n_cells:
        return None
    return plotter.add_mesh(mesh, color=color, opacity=opacity,
                            show_edges=False, smooth_shading=False,
                            lighting=True, name=key)


def add_pressure_ticks(plotter, daten):
    """Veraltet (TS-VP06): die zentrale Achsenaktualisierung (axes_view.apply_axes) zeichnet alle Ticks."""
    return None


def _visible_bounds(actors):
    """Gemeinsame Bounds aller sichtbaren Geometrie-Actors (ohne Achsenbeschriftung)."""
    return visible_bounds(actors)


def update_pvt_axes(plotter, actors, daten):
    """Achsen und Ticks auf die aktuell sichtbaren Ebenen anpassen (zentrale Achsenaktualisierung, TS-VP06).

    Die Kamera bleibt unveraendert; es werden nur Achsenrahmen, Ticks und Titel neu gesetzt.
    """
    apply_axes(plotter, "pvt", daten, actors)


def render_pvt(plotter, daten):
    """p-v-T: randkonforme Einphasenflaeche und optionales Zweiphasennetz."""
    # AP2: die korrekte Flaeche ist immer der sichtbare surface-Actor;
    # das alte regulaere, an der Binodale nicht randkonforme Netz entfaellt.
    single_phase_grid = create_single_phase_surface(None, daten)
    if single_phase_grid.n_cells == 0:
        raise ValueError("p-v-T: Keine gueltigen einphasigen Oberflaechenzellen.")

    plotter.clear()
    actors = {}
    actors["surface"] = plotter.add_mesh(
        single_phase_grid, color="lightblue", opacity=0.65, show_edges=True
    )
    # Vorhandenen GUI-Actor-Vertrag beibehalten; kein doppeltes Mesh zeichnen.
    actors["surface_single_phase"] = None

    # Thermodynamisch getrennte homogene Flaechen auf gemeinsamer Vertextabelle.
    meshes, info = homogeneous_shared_meshes(daten)
    for key, (color, opacity) in HOMOGENEOUS_STYLE.items():
        actors[key] = add_homogeneous_surface(plotter, daten, key, color=color, opacity=opacity,
                                              meshes=meshes)
    # HEOS-Rohfortsetzung und Ersatzflaeche nie gleichzeitig deckungsgleich: Voreinstellung = Ersatz.
    if actors["unstable_surface"] is not None:
        actors["unstable_surface"].visibility = False
    daten["homogeneous_mesh_info"] = info
    plotter._thermosurface_meshes = meshes   # fuer Hover (siehe hover_text_for_point)

    actors["zweiphasengebiet"] = add_two_phase_surface(plotter, daten)
    if actors["zweiphasengebiet"] is not None:
        actors["zweiphasengebiet"].visibility = False

    actors["isothermen"] = create_isotherm_lines(plotter, daten)
    actors["isochoren"] = create_isochore_lines(plotter, daten)
    actors["isobaren"] = create_isobare_lines(plotter, daten)
    actors["binodale"] = add_binodals(plotter, daten)
    actors["critical_point"] = add_critical_point(plotter, daten)
    # Die Spinodale ist eine separate homogene EOS-Stabilitaetsgrenze.
    # Das Einphasenmesh definiert keine physikalische Untergrenze der Kurve.
    actors["spinodale"] = add_spinodals(plotter, daten, min_z=None)

    actors["pressure_ticks"] = None
    plotter.reset_camera()
    apply_axes(plotter, "pvt", daten, actors)
    plotter.render()
    LOGGER.debug("AP2: Randkonforme Einphasenflaeche: Punkte=%d, Zellen=%d",
                 single_phase_grid.n_points, single_phase_grid.n_cells)
    return actors

def add_spinodals(plotter, daten, min_z=None):

    critical = daten["critical_point"]

    # Identische numerische Randknoten wie die neuen metastabilen Flaechen.
    # Fallback auf historische Kurven nur ohne ausreichende neue Stützstellen.
    def choose(branch):
        candidate = np.asarray(daten.get("homogeneous_spinodal_" + branch, []))
        surface_key = ('metastable_liquid_surface' if branch == 'liquid'
                       else 'metastable_vapor_surface')
        info = daten.get('homogeneous_validation', {}).get('surfaces', {}).get(surface_key, {})
        # Hat die neue homogene Oberflaeche Zellen, MUSS ihr Rand verwendet
        # werden. Keine andere Legacy-Spinodale geometrisch ansetzen.
        if info.get('candidate_triangles', 0) > 0:
            return candidate
        present = (np.isfinite(candidate).all(axis=1) &
                   ((candidate[:, 2] > 0) if daten.get('pressure_axis', 'log') == 'log' else True)) if candidate.ndim == 2 and candidate.shape[1] == 3 else []
        if np.count_nonzero(present) >= 2:
            return candidate
        return daten['spinodale_' + branch]
    axis = daten.get("pressure_axis", "log")
    liquid_actors = _draw_curves(plotter, [choose("liquid")], critical, color="purple", min_z=min_z, pressure_axis=axis)
    vapor_actors = _draw_curves(plotter, [choose("vapor")], critical, color="darkorange", min_z=min_z, pressure_axis=axis)

    return liquid_actors + vapor_actors

def _mesh_from_rows(rows, critical, pressure_axis="log"):
    """
    Vernetzt Temperaturzeilen mit unterschiedlicher
    Anzahl an Volumenpunkten.

    Jede Zeile enthält Punkte (T, v, p),
    sortiert nach aufsteigendem v.
    """

    if len(rows) < 2:
        return pv.PolyData()

    rows = [np.asarray(r, dtype=float) for r in rows]

    for row in rows:
        if row.ndim != 2 or row.shape[1] != 3:
            raise ValueError("Ungültige Temperaturzeile.")
        if len(row) == 0:
            raise ValueError("Leere Temperaturzeile.")
        if not np.isfinite(row).all():
            raise ValueError("Nichtendliche Zustandskoordinaten.")
        if np.any(np.diff(row[:, 1]) <= 0):
            raise ValueError("Volumenwerte nicht aufsteigend.")

    points = np.concatenate(rows, axis=0)

    offsets = np.cumsum([0] + [len(r) for r in rows[:-1]])

    # Anzahl der benötigten Dreiecke vorab bestimmen
    n_triangles = sum(len(a) + len(b) - 2 for a, b in zip(rows[:-1], rows[1:]))

    faces = np.empty((n_triangles, 4), dtype=np.int64)

    k = 0

    for r in range(len(rows) - 1):

        a = rows[r]
        b = rows[r + 1]

        # Normierte Position innerhalb jeder Volumenzeile
        da = max(a[-1, 1] - a[0, 1], np.finfo(float).tiny)
        db = max(b[-1, 1] - b[0, 1], np.finfo(float).tiny)

        sa = (a[:, 1] - a[0, 1]) / da
        sb = (b[:, 1] - b[0, 1]) / db

        i = 0
        j = 0

        while i < len(a) - 1 or j < len(b) - 1:

            ai = offsets[r] + i
            bj = offsets[r + 1] + j

            if (j == len(b) - 1 or (i < len(a) - 1 and sa[i + 1] <= sb[j + 1])):
                faces[k] = (3, ai, bj, ai + 1)
                i += 1

            else:
                faces[k] = (3, ai, bj, bj + 1)
                j += 1

            k += 1

    xyz = to_xyz(points, critical, pressure_axis)

    if not np.isfinite(xyz).all():
        raise ValueError("Ungültige PyVista-Koordinaten.")

    mesh = pv.PolyData(xyz, faces.ravel())
    attach_hover_metadata(mesh, points, "pvt")
    # Punktdaten fuer den Hover (siehe pick_hover_text); Rohwerte T, v, p, source=HEOS
    mesh.point_data["T"] = points[:, 0]
    mesh.point_data["v"] = points[:, 1]
    mesh.point_data["p_raw"] = points[:, 2]
    mesh.point_data["source"] = np.zeros(len(points), dtype=np.uint8)
    return mesh


def create_single_phase_surface(grid, daten):
    """Randkonforme Fluessig-, Dampf- und Hochtemperatur-Teilnetze.

    Siede-/Taupunkte werden direkt aus denselben Binodalarrays wie die
    Zweiphasenflaeche uebernommen. ``grid`` bleibt aus API-Gruenden erhalten.
    """

    surface = np.asarray(daten["surface"], dtype=float)

    liquid = np.asarray(daten["binodal_liquid"], dtype=float)
    vapor = np.asarray(daten["binodal_vapor"], dtype=float)

    critical = daten["critical_point"]
    Tc, pc, vc = critical

    if len(liquid) < 3 or len(vapor) < 3:
        raise ValueError("AP2: Fuer die randkonforme Flaeche fehlen Binodalpunkte.")

    # Kritischen Punkt zunächst aus der
    # Sättigungszuordnung ausschließen
    sat_map = {
        float(l[0]): (l, g)
        for l, g in zip(liquid, vapor)
        if l[0] < Tc}

    liquid_rows = []
    vapor_rows = []
    supercritical_rows = []

    critical_tvp = np.array([Tc, vc, pc], dtype=float)

    def valid_row(points):
        """
        Ungültige Randpunkte entfernen.
        Innere Datenlücken nicht überbrücken.
        """
        points = np.asarray(points, dtype=float)

        valid = (
            np.isfinite(points).all(axis=1)
            & (points[:, 1] > 0)
            & (points[:, 2] > 0))

        indices = np.flatnonzero(valid)

        if len(indices) == 0:
            return np.empty((0, 3))

        if indices[-1] - indices[0] + 1 != len(indices):
            raise ValueError("Innere ungültige Punkte in einer " "Oberflächenzeile.")

        return points[indices[0]:indices[-1] + 1]

    for row in surface:

        T = float(row[0, 0])
        v = row[:, 1]

        # -----------------------------------
        # 1. Unterkritische Zustände
        # -----------------------------------
        if T < Tc:
            if T not in sat_map:
                raise ValueError(f"Kein Sättigungszustand für T={T} K.")

            sat_l, sat_g = sat_map[T]

            v_l = sat_l[1]
            v_g = sat_g[1]

            if not (v[0] < v_l < v_g < v[-1]):
                raise ValueError(f"Sättigungsvolumina bei T={T} K " "liegen außerhalb des " "Oberflächen-Volumenrasters.")

            # Bestehende HEOS-Punkte beibehalten
            left = valid_row(row[v < v_l])
            right = valid_row(row[v > v_g])

            # Exakte Sättigungspunkte als Randknoten
            left = np.vstack((left, sat_l))
            right = np.vstack((sat_g, right))

            liquid_rows.append(left)
            vapor_rows.append(right)

        # -----------------------------------
        # 2. Kritische Temperatur
        # -----------------------------------
        elif T == Tc:

            left = valid_row(row[v < vc])
            right = valid_row(row[v > vc])

            left = np.vstack((left, critical_tvp))
            right = np.vstack((critical_tvp, right))

            # Gesamte kritische Temperaturzeile
            full = np.vstack((left[:-1], critical_tvp, right[1:]))

            liquid_rows.append(left)
            vapor_rows.append(right)
            supercritical_rows.append(full)

        # -----------------------------------
        # 3. Überkritische Zustände
        # -----------------------------------
        else:
            supercritical_rows.append(valid_row(row))

    # -----------------------------------
    # 4. Die drei Teilnetze erstellen
    # -----------------------------------
    meshes = []

    for rows in (liquid_rows, vapor_rows, supercritical_rows):
        if len(rows) >= 2:
            mesh = _mesh_from_rows(rows, critical, daten.get("pressure_axis", "log"))

            if mesh.n_cells > 0:
                meshes.append(mesh)

    if not meshes:
        return pv.PolyData()

    # Doppelte identische Randknoten innerhalb
    # der einphasigen Teilnetze zusammenführen
    return pv.merge(meshes, merge_points=True)

def _pht_two_phase_mesh(daten):
    """
    Erzeugt die Zweiphasenfläche ohne
    entartete kritische Q-Zeile.
    """

    surface = np.asarray(daten["two_phase_surface"], dtype=float)

    critical = daten["critical_point"]

    # Letzte kritische Q-Zeile aus dem
    # regulären Oberflächennetz entfernen.
    subcritical = surface[:-1]

    mesh = _pht_surface_mesh(subcritical, critical)

    # Letztes unterkritisches Q-Profil
    last_row = subcritical[-1]

    Tc, pc, hc = critical

    critical_point = np.array([[Tc, hc, pc]], dtype=float)

    # Randprofil und kritischen Punkt verbinden
    cap_points = np.vstack((last_row, critical_point))

    cap_xyz = to_xyz_pht(cap_points, critical)

    n_Q = len(last_row)

    cap_faces = np.column_stack(
        (   np.full(n_Q - 1, 3, dtype=np.int64),
            np.arange(n_Q - 1),
            np.arange(1, n_Q),
            np.full(n_Q - 1, n_Q))).ravel()

    cap_mesh = pv.PolyData(cap_xyz, cap_faces)
    attach_hover_metadata(cap_mesh, cap_points, "pht")

    return mesh.merge(cap_mesh, merge_points=True)

def _pht_add_curve(plotter, points, critical, color):
    """Zeichnet eine zusammenhängende p-h-T-Kurve."""

    points = np.asarray(points, dtype=float)

    if (points.ndim != 2 or points.shape[1] != 3 or len(points) < 2):
        return None

    xyz = to_xyz_pht(points, critical)

    if not np.isfinite(xyz).all():
        raise ValueError("Ungültige Koordinaten der p-h-T-Kurve.")

    poly = pv.lines_from_points(xyz)

    return plotter.add_mesh(poly, color=color, line_width=3)


def _pht_add_critical_point(plotter, critical):
    """Zeichnet den kritischen Punkt in p-h-T."""

    Tc, pc, hc = critical

    point = np.array([[Tc, hc, pc]], dtype=float)

    xyz = to_xyz_pht(point, critical)

    return plotter.add_mesh(pv.PolyData(xyz), color="black", point_size=12, render_points_as_spheres=True)

def render_pht(plotter, daten):
    """
    Visualisiert berechnete p-h-T-Daten.

    Erwartet das Rückgabe-Dictionary aus
    coolprop_heos_pht.berechne_pht_daten().
    """

    critical = daten["critical_point"]

    # ---------------------------------------
    # 1. Geometrien vorbereiten
    # ---------------------------------------

    meshes = {}

    for name in ("surface_liquid", "surface_vapor", "surface_high_temperature"):
        meshes[name] = _pht_surface_mesh(daten[name], critical)

    phase_mesh = _pht_two_phase_mesh(daten)

    # ---------------------------------------
    # 2. Vorherige Darstellung entfernen
    # ---------------------------------------

    plotter.clear()

    # ---------------------------------------
    # 3. Actor-Dictionary erzeugen
    # ---------------------------------------

    actors = {
        "surface": [],
        "surface_single_phase": None,
        "zweiphasengebiet": None,
        "binodale": [],
        "critical_point": None,
        "isothermen": [],
        "isobaren": [],
        "isochoren": [],
        "spinodale": []
    }

    # ---------------------------------------
    # 4. Drei einphasige Teilflächen
    # ---------------------------------------

    colors = {
        "surface_liquid": "#47ADAC",
        "surface_vapor": "#F09C9C",
        "surface_high_temperature": "#EB0E0E"}

    for name, mesh in meshes.items():
        actor = plotter.add_mesh(mesh, color=colors[name], opacity=0.65, show_edges=True, lighting=True, smooth_shading=True)
        actors["surface"].append(actor)

    # ---------------------------------------
    # 5. Zweiphasengebiet
    # ---------------------------------------

    actors["zweiphasengebiet"] = plotter.add_mesh(phase_mesh, color="orange", opacity=0.85, show_edges=False, lighting=True, smooth_shading=True)

    # ---------------------------------------
    # 6. Siede- und Taulinie
    # ---------------------------------------

    liquid_actor = _pht_add_curve(plotter, daten["binodal_liquid"], critical, "blue")
    vapor_actor = _pht_add_curve(plotter, daten["binodal_vapor"], critical, "red")

    actors["binodale"] = [liquid_actor, vapor_actor]

    # ---------------------------------------
    # 7. Kritischer Punkt
    # ---------------------------------------

    actors["critical_point"] = _pht_add_critical_point(plotter, critical)

    # ---------------------------------------
    # 7b. Isolinien (T, p, v = konstant; NaN-getrennte Abschnitte)
    # ---------------------------------------

    for key, color in (("isothermen", "red"), ("isobaren", "blue"), ("isochoren", "green")):
        for curve in daten.get(key, ()):
            actors[key].extend(_new_curve(plotter, curve, critical, to_xyz_pht, color))

    # ---------------------------------------
    # 8. Achsen
    # ---------------------------------------

    plotter.reset_camera()
    apply_axes(plotter, "pht", daten, actors)

    # ---------------------------------------
    # 9. Darstellung aktualisieren
    # ---------------------------------------

    plotter.render()

    return actors

# ============================================================
# p-s-T: Entropie statt Enthalpie als zweite physikalische Achse
# ============================================================

def _pst_surface_mesh(surface, critical):
    surface = np.asarray(surface, dtype=float)

    if surface.ndim != 3 or surface.shape[-1] != 3 or min(surface.shape[:2]) < 2:
        raise ValueError("p-s-T-Oberflaeche erwartet (n_T, n_p, 3).")
    
    xyz = to_xyz_pst(surface, critical)

    if not np.isfinite(xyz).all():
        raise ValueError("p-s-T-Mesh enthaelt nichtendliche XYZ-Koordinaten.")
    
    grid = pv.StructuredGrid(xyz[..., 0], xyz[..., 1], xyz[..., 2])
    attach_hover_metadata(grid, surface.reshape((-1, 3), order="F"), "pst")
    return grid.extract_surface().triangulate()


def _pst_two_phase_mesh(daten):
    """Ohne degenerierte kritische Q-Zeile, mit Dreiecksfaecher zur Spitze."""
    surface = np.asarray(daten["two_phase_surface"], dtype=float)
    critical = daten["critical_point"]
    mesh = _pst_surface_mesh(surface[:-1], critical)
    last_row = surface[-2]
    Tc, pc, sc = critical
    pts = np.vstack((last_row, np.array([[Tc, sc, pc]])))
    n_Q = len(last_row)
    faces = np.column_stack((np.full(n_Q-1, 3, dtype=np.int64), np.arange(n_Q-1), np.arange(1, n_Q), np.full(n_Q-1, n_Q))).ravel()
    cap = pv.PolyData(to_xyz_pst(pts, critical), faces)
    attach_hover_metadata(cap, pts, "pst")

    return mesh.merge(cap, merge_points=True)


def _pst_add_curve(plotter, points, critical, color):
    """NaN-getrennte Linienabschnitte unterstuetzen (z.B. Isochoren-Clipping)."""
    pts = np.asarray(points, dtype=float)

    if pts.ndim != 2 or pts.shape[-1] != 3:
        raise ValueError("p-s-T-Linie braucht (n,3).")
    
    xyz = to_xyz_pst(pts, critical)
    valid = np.isfinite(xyz).all(axis=1)
    boundaries = np.flatnonzero(np.diff(np.r_[False, valid, False].astype(np.int8)))

    actors = []

    for start, end in zip(boundaries[::2], boundaries[1::2]):
        if end - start >= 2:
            actors.append(plotter.add_mesh(pv.lines_from_points(xyz[start:end]), color=color, line_width=3))

    return actors


def render_pst(plotter, daten):
    """Eigenstaendiger p-s-T-Renderer mit dem bestehenden Actor-Vertrag."""
    critical = daten["critical_point"]
    surfaces = [_pst_surface_mesh(daten[name], critical) for name in ("surface_liquid", "surface_vapor", "surface_high_temperature")]
    phase_mesh = _pst_two_phase_mesh(daten)
    plotter.clear()
    actors = {
        "surface":      [],
        "surface_single_phase": None,
        "zweiphasengebiet": None, 
        "binodale":     [], "critical_point": None,
        "isothermen":   [], 
        "isobaren":     [], 
        "isochoren":    [], 
        "spinodale":    [],}
    
    for mesh in surfaces:
        actors["surface"].append(plotter.add_mesh(mesh, color="lightblue", opacity=0.65, show_edges=True, lighting=True, smooth_shading=True))
        
    actors["zweiphasengebiet"] = plotter.add_mesh(phase_mesh, color="orange", opacity=0.85, show_edges=False, lighting=True, smooth_shading=True)

    for name, color in (("binodal_liquid", "blue"), ("binodal_vapor", "red")):
        actors["binodale"].extend(_pst_add_curve(plotter, daten[name], critical, color))

    Tc, pc, sc = critical
    critical_xyz = to_xyz_pst(np.array([[Tc, sc, pc]]), critical)

    actors["critical_point"] = plotter.add_mesh(pv.PolyData(critical_xyz), color="black", point_size=12, render_points_as_spheres=True)

    for key, color in (("isothermen", "red"), ("isobaren", "blue"), ("isochoren", "green")):
        for curve in daten[key]:
            actors[key].extend(_pst_add_curve(plotter, curve, critical, color))

    plotter.reset_camera()
    apply_axes(plotter, "pst", daten, actors)
    plotter.render()

    return actors


# ============================================================
# h-s-p und u-v-T: gemeinsame, diagrammspezifisch transformierte
# PyVista-Darstellung. Die drei Bestandsrenderer bleiben unangetastet.
# ============================================================

def _new_surface_mesh(surface, critical, transform, diagram_id=None):
    points = np.asarray(surface, dtype=np.float64)

    if (points.ndim != 3 or points.shape[-1] != 3 or min(points.shape[:2]) < 2):
        raise ValueError("h-s-p/u-v-T: Oberflaeche benoetigt (n1,n2,3).")
    
    xyz = transform(points, critical)

    if not np.isfinite(xyz).all():
        raise ValueError("h-s-p/u-v-T: Oberflaechen-Mesh nicht endlich.")
    
    grid = pv.StructuredGrid(xyz[..., 0], xyz[..., 1], xyz[..., 2])
    if diagram_id is not None:
        attach_hover_metadata(grid, points.reshape((-1, 3), order="F"), diagram_id)
    return grid.extract_surface().triangulate()


def _new_phase_mesh(daten, transform, diagram_id=None):
    """Kritische Q-Zeile weglassen und saubere Dreiecksspitze einfuegen."""
    two_phase = np.asarray(daten["two_phase_surface"], dtype=np.float64)
    critical = daten["critical_point"]

    if (two_phase.ndim != 3 or two_phase.shape[-1] != 3 or two_phase.shape[0] < 3 or two_phase.shape[1] < 2):
        raise ValueError("Zweiphasenflaeche braucht >= 3 T- und >= 2 Q-Zeilen.")
    
    mesh = _new_surface_mesh(two_phase[:-1], critical, transform, diagram_id)
    last = two_phase[-2]
    n_q = len(last)
    points = np.vstack((last, np.asarray(critical, dtype=float)[None, :]))
    faces = np.column_stack((np.full(n_q - 1, 3, dtype=np.int64), np.arange(n_q - 1), np.arange(1, n_q), np.full(n_q - 1, n_q))).ravel()
    cap = pv.PolyData(transform(points, critical), faces)
    if diagram_id is not None:
        attach_hover_metadata(cap, points, diagram_id)

    return mesh.merge(cap, merge_points=True)


def _new_curve(plotter, curve, critical, transform, color):
    points = np.asarray(curve, dtype=np.float64)

    if points.ndim != 2 or points.shape[-1] != 3:
        raise ValueError("h-s-p/u-v-T-Linie braucht (n,3).")
    
    xyz = transform(points, critical)
    finite = np.isfinite(xyz).all(axis=1)
    bounds = np.flatnonzero(np.diff(np.r_[False, finite, False].astype(np.int8)))

    actors = []

    for start, stop in zip(bounds[::2], bounds[1::2]):
        if stop - start >= 2:
            polyline = pv.lines_from_points(xyz[start:stop])
            actors.append(plotter.add_mesh(polyline, color=color, line_width=3))

    return actors


def _render_new(plotter, daten, transform, diagram_id):
    """Gemeinsames Actors-Dictionary fuer Checkboxen aus gui.py."""
    critical = daten["critical_point"]
    meshes = [_new_surface_mesh(daten[key], critical, transform, diagram_id)
              for key in ("surface_liquid", "surface_vapor", "surface_high_temperature")]
    phase_mesh = _new_phase_mesh(daten, transform, diagram_id)
    plotter.clear()

    actors = {
        "surface":      [], 
        "surface_single_phase": None,
        "zweiphasengebiet": None, 
        "binodale":     [],
        "critical_point": None, 
        "isothermen":   [],
        "isobaren":     [], 
        "isochoren":    [], 
        "spinodale":    [],}
    
    for mesh in meshes:
        actors["surface"].append(plotter.add_mesh(mesh, color="lightblue", opacity=0.65, show_edges=True, lighting=True, smooth_shading=True))

    actors["zweiphasengebiet"] = plotter.add_mesh(phase_mesh, color="orange", opacity=0.85, show_edges=False, lighting=True, smooth_shading=True)

    for key, color in (("binodal_liquid", "blue"), ("binodal_vapor", "red")):
        actors["binodale"].extend(_new_curve(plotter, daten[key], critical, transform, color))

    actors["critical_point"] = plotter.add_mesh(pv.PolyData(transform(np.asarray(critical, dtype=float)[None, :], critical)), color="black", point_size=12, render_points_as_spheres=True)

    for key, color in (("isothermen", "red"), ("isobaren", "blue"), ("isochoren", "green")):
        for curve in daten.get(key, ()):
            actors[key].extend(_new_curve(plotter, curve, critical, transform, color))

    plotter.reset_camera()
    apply_axes(plotter, diagram_id, daten, actors)
    plotter.render()
    return actors


def render_hsp(plotter, daten):
    """3D-Mollier-Diagramm h-s-p (aus coolprop_heos_hsp)."""
    return _render_new(plotter, daten, to_xyz_hsp, "hsp")


def render_uvt(plotter, daten):
    """3D-innere-Energie-Diagramm u-v-T (aus coolprop_heos_uvt)."""
    return _render_new(plotter, daten, to_xyz_uvt, "uvt")


# ============================================================
# g-p-T und f-v-T: freie Enthalpie / freie Energie.
# Gleiche Bausteine wie h-s-p/u-v-T (_new_surface_mesh, _new_curve, _new_phase_mesh); zusaetzlich die optionalen
# homogenen Ebenen (metastabil/instabil) und Spinodalen, die nie mit den Gleichgewichtsflaechen vermischt werden.
# ============================================================

_FREE_ENERGY_HOMOGENEOUS = ("metastable_liquid_surface", "metastable_vapor_surface", "unstable_surface")


def _masked_grid_mesh(values, critical, transform, diagram_id):
    """Netz aus (n,m,3)-Knoten; nur Zellen mit >= 3 darstellbaren Ecken, Luecken (NaN, p <= 0 auf Log-Achse) bleiben offen."""
    values = np.asarray(values, dtype=np.float64)

    if values.ndim != 3 or values.shape[-1] != 3 or min(values.shape[:2]) < 2:
        return None

    xyz = transform(values, critical)
    valid = np.isfinite(xyz).all(axis=2)
    flat = np.where(valid[..., None], xyz, 0.0).reshape((-1, 3))
    triangles, cells = grid_triangles(valid)
    triangles, cells, _ = drop_degenerate(triangles, cells, flat)

    if not len(triangles):
        return None

    used = np.unique(triangles.ravel())
    remap = np.full(len(flat), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    faces = np.column_stack((np.full(len(triangles), 3, dtype=np.int64), remap[triangles]))

    mesh = pv.PolyData(flat[used], faces.ravel())
    attach_hover_metadata(mesh, values.reshape((-1, 3))[used], diagram_id, "HEOS_homogeneous")
    return mesh


def _render_free_energy(plotter, daten, transform, diagram_id, *, two_phase, surface_colors, binodal_colors):
    """Actors-Dictionary fuer g-p-T/f-v-T (Schluessel wie die anderen Diagramme + homogene Ebenen + Spinodalen)."""
    critical = daten["critical_point"]
    meshes = [_new_surface_mesh(daten[key], critical, transform, diagram_id)
              for key in ("surface_liquid", "surface_vapor", "surface_high_temperature")]
    phase_mesh = _new_phase_mesh(daten, transform, diagram_id) if two_phase else None
    plotter.clear()

    actors = {"surface": [], "surface_single_phase": None, "zweiphasengebiet": None, "binodale": [],
              "critical_point": None, "isothermen": [], "isobaren": [], "isochoren": [], "spinodale": []}

    for mesh, color in zip(meshes, surface_colors):
        actors["surface"].append(plotter.add_mesh(mesh, color=color, opacity=0.65, show_edges=True, lighting=True,
                                                  smooth_shading=True))

    if phase_mesh is not None:
        actors["zweiphasengebiet"] = plotter.add_mesh(phase_mesh, color="orange", opacity=0.85, show_edges=False,
                                                      lighting=True, smooth_shading=True)

    for key, color in binodal_colors:
        actors["binodale"].extend(_new_curve(plotter, daten[key], critical, transform, color))

    actors["critical_point"] = plotter.add_mesh(pv.PolyData(transform(np.asarray(critical, dtype=float)[None, :], critical)),
                                                color="black", point_size=12, render_points_as_spheres=True)

    for key, color in (("isothermen", "red"), ("isobaren", "blue"), ("isochoren", "green")):
        for curve in daten.get(key, ()):
            actors[key].extend(_new_curve(plotter, curve, critical, transform, color))

    for key in _FREE_ENERGY_HOMOGENEOUS:
        if key not in daten:
            continue

        mesh = _masked_grid_mesh(daten[key], critical, transform, diagram_id)

        if mesh is not None:
            color, opacity = HOMOGENEOUS_STYLE[key]
            actors[key] = plotter.add_mesh(mesh, color=color, opacity=opacity, show_edges=True, lighting=True,
                                           smooth_shading=False)

    for key, color in (("spinodal_liquid", "purple"), ("spinodal_vapor", "darkorange")):
        if key in daten:
            actors["spinodale"].extend(_new_curve(plotter, daten[key], critical, transform, color))

    plotter.reset_camera()
    apply_axes(plotter, diagram_id, daten, actors)
    plotter.render()

    return actors


def render_gpt(plotter, daten):
    """3D-Gibbs-Diagramm g-p-T (aus coolprop_heos_gpt). Zweiphasengebiet = Sattigungslinie (keine Flaeche)."""
    axis = daten.get("pressure_axis", "log")
    raster = daten.get("raster") or {}
    window = (raster["p_min"], raster["p_max"]) if "p_min" in raster and "p_max" in raster else None

    def transform(points, critical):
        return to_xyz_gpt(points, critical, axis, window)

    return _render_free_energy(plotter, daten, transform, "gpt",
                               two_phase=False, surface_colors=("#47ADAC", "#F09C9C", "#EB0E0E"),
                               binodal_colors=(("binodal_liquid", "blue"),))


def render_fvt(plotter, daten):
    """3D-Helmholtz-Diagramm f-v-T (aus coolprop_heos_fvt)."""
    return _render_free_energy(plotter, daten, to_xyz_fvt, "fvt",
                               two_phase=True, surface_colors=("lightblue",) * 3,
                               binodal_colors=(("binodal_liquid", "blue"), ("binodal_vapor", "red")))
