
"""
Datei: state_store.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Zentrale Verwaltung von Rasterreferenzen, Koordinatenadaptern
und berechneten thermodynamischen Zuständen.

Aufgaben:
- Registrierung diagrammspezifischer Zustandsadapter
- Zuordnung von Rasterpunkten zu eindeutigen Referenzen
- Verwaltung von Flächen-, Phasen- und Quelldaten
- Aufbau durchsuchbarer Diagrammindizes
- Speicherung berechneter Zustände im RAM-Cache
- Optionale persistente Speicherung über HDF5

Abhängigkeiten:
- thermo_state.py: Thermodynamische Zustandsdatensätze
- diagram_coords.py: Koordinatentransformationen
- NumPy und SciPy: Datenverarbeitung und Punktsuche
- h5py: Optionale HDF5-Speicherung

Verwendet von:
- gui.py: Initialisierung und Verwaltung des Zustandsspeichers
- state_picker.py / state_resolver.py
- state_projection.py
"""

import hashlib
import logging
from collections import OrderedDict, namedtuple
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

import diagram_coords as dc
from thermo_state import SCHEMA_VERSION, UNITS, ThermoState

LOGGER = logging.getLogger("thermosurface.state_store")

NodeRef = namedtuple("NodeRef", "surface index")

# Flaechenarten
SINGLE_PHASE = "single_phase"
SATURATION = "saturation"          # Binodalpunkte (Q = 0 bzw. 1)
CRITICAL = "critical"
TWO_PHASE = "two_phase"
HOMOGENEOUS = "homogeneous"        # metastabile/instabile HEOS-Fortsetzung (p-v-T)
SUBSTITUTE = "substitute"          # Ersatzflaeche (Modell, keine EOS)

# Prioritaet bei gleichen Koordinaten: exakte Saettigungs-/kritische Knoten vor den Gitterflaechen.
_PRIORITY = {CRITICAL: 0, SATURATION: 1, SINGLE_PHASE: 2, TWO_PHASE: 2, HOMOGENEOUS: 2, SUBSTITUTE: 2}


def format_ref(ref):
    return f"{ref.surface}:" + ",".join(str(int(i)) for i in ref.index)


@dataclass
class SourceSurface:
    """Eine referenzierte Quellflaeche eines Diagramms (View auf das Datenarray)."""
    name: str                       # Schluessel im Ergebnis-Dictionary bzw. Pseudoname
    kind: str
    array: np.ndarray               # (n, m, 3) Gitter oder (n, 3) Kurve; native Spalten gemaess Adapter.axes
    actor_key: Optional[str] = None
    phase: Optional[str] = None     # liquid | vapor | supercritical | unstable (Hinweis, kein Gleichgewichtsurteil)
    q_values: Optional[np.ndarray] = None   # TWO_PHASE: Dampfgehalt je Spalte
    q_fixed: Optional[float] = None         # SATURATION: 0.0 (Siedelinie) / 1.0 (Taulinie)
    source_codes: Optional[np.ndarray] = None  # (n, m): 0 = HEOS, 1 = Ersatz (nur SUBSTITUTE)
    note: str = ""
    # Zusatzgroessen je Knoten (gleiche Form wie array ohne letzte Achse), z. B. {"v": ...} fuer die homogenen
    # Ebenen von g-p-T, deren Quellspalten (T, p, g) das Volumen nicht enthalten.
    aux: Optional[dict] = None

    @property
    def n_nodes(self):
        return int(np.prod(self.array.shape[:-1]))

    def node(self, index):
        return np.asarray(self.array[tuple(int(i) for i in index)], dtype=float)

    def node_aux(self, name, index):
        if not self.aux or name not in self.aux:
            return None
        value = float(self.aux[name][tuple(int(i) for i in index)])
        return value if np.isfinite(value) else None

    def node_source(self, index):
        if self.source_codes is None:
            return 0
        return int(self.source_codes[tuple(int(i) for i in index)])


# ============================================================ Adapter (Koordinatenadapter je Diagramm)

class StateAdapter:
    """Diagramm-spezifische Abbildung Daten <-> Darstellung <-> Quellflaechen."""

    diagram_id = ""
    axes = ("", "", "")                     # native Spaltenreihenfolge der Punktarrays
    two_phase_axis = "T"                    # Temperatur- bzw. Druckspalte des Zweiphasengitters ('T' | 'p')
    mass_axes = frozenset()                 # Achsen in J/kg bzw. J/(kg K) (massenbezogen)
    surface_names = ()                      # (daten-Key, Phase) der einphasigen Teilflaechen

    def pressure_axis(self, daten):
        return "log"

    def to_display(self, points, daten):
        raise NotImplementedError

    def from_display(self, xyz, daten):
        raise NotImplementedError

    def critical_native(self, daten):
        raise NotImplementedError

    def surfaces(self, daten):
        raise NotImplementedError

    # --- gemeinsame Hilfen
    def _saturation(self, daten):
        out = []
        for key, q in (("binodal_liquid", 0.0), ("binodal_vapor", 1.0)):
            arr = daten.get(key)
            if arr is not None and np.asarray(arr).ndim == 2 and len(arr):
                out.append(SourceSurface(key, SATURATION, np.asarray(arr, dtype=float), q_fixed=q))
        return out

    def _critical_surface(self, daten):
        return SourceSurface("critical_point", CRITICAL, np.asarray(self.critical_native(daten), dtype=float)[None, :])

    def _two_phase(self, daten):
        arr = daten.get("two_phase_surface")
        if arr is None or np.asarray(arr).ndim != 3 or np.asarray(arr).size == 0:
            return None
        return SourceSurface("two_phase_surface", TWO_PHASE, np.asarray(arr, dtype=float), actor_key="zweiphasengebiet",
                             q_values=np.asarray(daten["Q_values"], dtype=float))

    def _with_boundaries(self, surface, daten, actor_key):
        surface.actor_key = actor_key
        return [self._critical_surface(daten), *self._saturation(daten), surface]


class _PvtAdapter(StateAdapter):
    diagram_id = "pvt"
    axes = ("T", "v", "p")

    def pressure_axis(self, daten):
        return daten.get("pressure_axis", "log")

    def to_display(self, points, daten):
        return dc.to_xyz(points, daten["critical_point"], self.pressure_axis(daten))

    def from_display(self, xyz, daten):
        return dc.from_xyz(xyz, daten["critical_point"], self.pressure_axis(daten))

    def critical_native(self, daten):
        Tc, pc, vc = daten["critical_point"]
        return (Tc, vc, pc)

    def surfaces(self, daten):
        out = []
        main = SourceSurface("surface", SINGLE_PHASE, np.asarray(daten["surface"], dtype=float), actor_key="surface")
        out += [self._critical_surface(daten), *self._saturation(daten), main]
        for s in out[:-1]:
            s.actor_key = "surface"
        two = self._two_phase(daten)
        if two is not None:
            out.append(two)
        for key, phase in (("metastable_liquid_surface", "liquid"), ("metastable_vapor_surface", "vapor"),
                           ("unstable_surface", "unstable")):
            arr = daten.get(key)
            if arr is not None and np.asarray(arr).ndim == 3:
                out.append(SourceSurface(key, HOMOGENEOUS, np.asarray(arr, dtype=float), actor_key=key, phase=phase))
        arr = daten.get("replacement_surface")
        if arr is not None and np.asarray(arr).ndim == 3:
            field_ = daten.get("replacement_field")
            codes = np.asarray(field_.source) if field_ is not None else None
            out.append(SourceSurface("replacement_surface", SUBSTITUTE, np.asarray(arr, dtype=float),
                                     actor_key="replacement_surface", source_codes=codes,
                                     note="Ersatzflaeche (Modell), keine EOS-Auswertung"))
        return out


class _SpecificAdapter(StateAdapter):
    """pht, pst, hsp, uvt: drei einphasige Teilflaechen + Zweiphasenflaeche (gleiche Datenstruktur)."""
    surface_names = (("surface_liquid", "liquid"), ("surface_vapor", "vapor"), ("surface_high_temperature", "supercritical"))

    def surfaces(self, daten):
        out = [self._critical_surface(daten), *self._saturation(daten)]
        for s in out:
            s.actor_key = "surface"
        for key, phase in self.surface_names:
            arr = daten.get(key)
            if arr is not None and np.asarray(arr).ndim == 3:
                out.append(SourceSurface(key, SINGLE_PHASE, np.asarray(arr, dtype=float), actor_key="surface", phase=phase))
        two = self._two_phase(daten)
        if two is not None:
            out.append(two)
        return out


class _PhtAdapter(_SpecificAdapter):
    diagram_id = "pht"
    axes = ("T", "h", "p")
    mass_axes = frozenset({"h"})

    def to_display(self, points, daten):
        return dc.to_xyz_pht(points, daten["critical_point"])

    def from_display(self, xyz, daten):
        return dc.from_xyz_pht(xyz, daten["critical_point"])

    def critical_native(self, daten):
        Tc, pc, hc = daten["critical_point"]
        return (Tc, hc, pc)


class _PstAdapter(_SpecificAdapter):
    diagram_id = "pst"
    axes = ("T", "s", "p")
    mass_axes = frozenset({"s"})

    def to_display(self, points, daten):
        return dc.to_xyz_pst(points, daten["critical_point"])

    def from_display(self, xyz, daten):
        return dc.from_xyz_pst(xyz, daten["critical_point"])

    def critical_native(self, daten):
        Tc, pc, sc = daten["critical_point"]
        return (Tc, sc, pc)


class _HspAdapter(_SpecificAdapter):
    diagram_id = "hsp"
    axes = ("h", "s", "p")
    two_phase_axis = "p"
    mass_axes = frozenset({"h", "s"})

    def to_display(self, points, daten):
        return dc.to_xyz_hsp(points, daten["critical_point"])

    def from_display(self, xyz, daten):
        return dc.from_xyz_hsp(xyz, daten["critical_point"])

    def critical_native(self, daten):
        return tuple(daten["critical_point"])


class _UvtAdapter(_SpecificAdapter):
    diagram_id = "uvt"
    axes = ("T", "v", "u")
    mass_axes = frozenset({"u"})

    def to_display(self, points, daten):
        return dc.to_xyz_uvt(points, daten["critical_point"])

    def from_display(self, xyz, daten):
        return dc.from_xyz_uvt(xyz, daten["critical_point"])

    def critical_native(self, daten):
        return tuple(daten["critical_point"])


class _FreeEnergyMixin:
    """Homogene (metastabile/instabile) HEOS-Ebenen von g-p-T und f-v-T als eigene Quellflaechen."""
    homogeneous_surfaces = (("metastable_liquid_surface", "liquid"), ("metastable_vapor_surface", "vapor"),
                            ("unstable_surface", "unstable"))

    def _homogeneous(self, daten):
        out = []
        volumes = daten.get("homogeneous_v") or {}
        for key, phase in self.homogeneous_surfaces:
            arr = daten.get(key)
            if arr is not None and np.asarray(arr).ndim == 3:
                aux = {"v": np.asarray(volumes[key], dtype=float)} if key in volumes else None
                out.append(SourceSurface(key, HOMOGENEOUS, np.asarray(arr, dtype=float), actor_key=key, phase=phase,
                                         aux=aux))
        return out


class _GptAdapter(_FreeEnergyMixin, _SpecificAdapter):
    """g-p-T: Zweiphasengebiet = Saettigungslinie; dort ist Q aus (T,p,g) nicht bestimmt.

    Die Siede-/Taulinie ist deshalb keine Pick-Quelle: Randknoten der Fluessigkeits- und Dampfflaeche werden als
    einphasige Randzustaende (Phase vorgegeben) rekonstruiert.
    """
    diagram_id = "gpt"
    axes = ("T", "p", "g")
    mass_axes = frozenset({"g"})
    saturation_degenerate = True

    def pressure_axis(self, daten):
        return daten.get("pressure_axis", "log")

    @staticmethod
    def pressure_window(daten):
        raster = daten.get("raster") or {}
        return (raster["p_min"], raster["p_max"]) if "p_min" in raster and "p_max" in raster else None

    def to_display(self, points, daten):
        return dc.to_xyz_gpt(points, daten["critical_point"], self.pressure_axis(daten), self.pressure_window(daten))

    def from_display(self, xyz, daten):
        return dc.from_xyz_gpt(xyz, daten["critical_point"], self.pressure_axis(daten))

    def critical_native(self, daten):
        return tuple(daten["critical_point"])

    def surfaces(self, daten):
        out = [self._critical_surface(daten)]
        out[0].actor_key = "surface"
        for key, phase in self.surface_names:
            arr = daten.get(key)
            if arr is not None and np.asarray(arr).ndim == 3:
                out.append(SourceSurface(key, SINGLE_PHASE, np.asarray(arr, dtype=float), actor_key="surface", phase=phase))
        return out + self._homogeneous(daten)


class _FvtAdapter(_FreeEnergyMixin, _SpecificAdapter):
    diagram_id = "fvt"
    axes = ("T", "v", "f")
    mass_axes = frozenset({"f"})

    def to_display(self, points, daten):
        return dc.to_xyz_fvt(points, daten["critical_point"])

    def from_display(self, xyz, daten):
        return dc.from_xyz_fvt(xyz, daten["critical_point"])

    def critical_native(self, daten):
        return tuple(daten["critical_point"])

    def surfaces(self, daten):
        return super().surfaces(daten) + self._homogeneous(daten)


ADAPTERS = {}


def register_adapter(adapter):
    ADAPTERS[adapter.diagram_id] = adapter
    return adapter


for _adapter in (_PvtAdapter(), _PhtAdapter(), _PstAdapter(), _HspAdapter(), _UvtAdapter(), _GptAdapter(),
                 _FvtAdapter()):
    register_adapter(_adapter)


def adapter_for(diagram_id):
    try:
        return ADAPTERS[diagram_id]
    except KeyError:
        raise ValueError(f"Kein Zustandsadapter fuer Diagramm {diagram_id!r} registriert") from None


# ============================================================ Index je Diagrammergebnis

class DiagramIndex:
    """Referenzen auf alle berechneten Oberflaechenpunkte eines Diagrammergebnisses (lazy aufgebaut)."""

    MATCH_RTOL = 1e-5   # relative Toleranz (float32-faehig) der Knotenzuordnung im Darstellungsraum

    def __init__(self, fluid, diagram_id, daten):
        self.fluid = fluid
        self.diagram_id = diagram_id
        self.daten = daten
        self.adapter = adapter_for(diagram_id)
        self._surfaces = None
        self._trees = {}
        self._grid_hash = None

    # ---------------------------------------------------------------- Flaechen
    @property
    def surfaces(self):
        if self._surfaces is None:
            self._surfaces = {s.name: s for s in self.adapter.surfaces(self.daten)}
        return self._surfaces

    def surface(self, name):
        return self.surfaces[name]

    def candidates(self, actor_key, sub_index=None):
        """Quellflaechen eines gerenderten Actors in Zuordnungsreihenfolge (sub_index: Teilflaeche der Actor-Liste)."""
        out = []
        for s in self.surfaces.values():
            if s.actor_key != actor_key and s.kind != CRITICAL:
                continue
            if sub_index is not None and s.kind == SINGLE_PHASE and self.adapter.surface_names:
                names = [n for n, _ in self.adapter.surface_names]
                if sub_index < len(names) and s.name != names[sub_index]:
                    continue
            out.append(s)
        return sorted(out, key=lambda s: _PRIORITY[s.kind])

    @property
    def n_nodes(self):
        return sum(int(np.count_nonzero(np.isfinite(s.array).all(axis=-1))) for s in self.surfaces.values())

    def iter_refs(self):
        """Alle endlichen Oberflaechenpunkte (Gitterflaechen, Saettigungs- und kritischer Punkt) als NodeRef."""
        for s in self.surfaces.values():
            finite = np.isfinite(s.array).all(axis=-1)
            for idx in zip(*np.nonzero(finite)):
                yield NodeRef(s.name, tuple(int(i) for i in idx))

    def native(self, ref):
        """Unveraenderter Quelldatensatz des Knotens (native Spalten gemaess adapter.axes)."""
        return self.surface(ref.surface).node(ref.index)

    # ---------------------------------------------------------------- Zuordnung Darstellung -> Knoten
    def _tree(self, name):
        if name not in self._trees:
            from scipy.spatial import cKDTree
            s = self.surface(name)
            flat = s.array.reshape(-1, 3)
            xyz = self.adapter.to_display(flat, self.daten)
            finite = np.isfinite(xyz).all(axis=1)
            ids = np.flatnonzero(finite)
            self._trees[name] = (cKDTree(xyz[ids]) if len(ids) else None, ids, xyz)
        return self._trees[name]

    def locate(self, actor_key, sub_index, xyz_points):
        """Ordnet Darstellungskoordinaten den Quellknoten zu: Liste von (NodeRef | None, Abstand)."""
        xyz_points = np.atleast_2d(np.asarray(xyz_points, dtype=float))
        cands = self.candidates(actor_key, sub_index)
        out = []
        for point in xyz_points:
            best = (None, np.inf)
            for s in cands:
                tree, ids, xyz = self._tree(s.name)
                if tree is None:
                    continue
                dist, k = tree.query(point)
                if dist < best[1] - 1e-15:
                    flat = int(ids[k])
                    best = (NodeRef(s.name, tuple(int(i) for i in np.unravel_index(flat, s.array.shape[:-1]))), float(dist))
                    if dist == 0.0:
                        break
            scale = max(1.0, float(np.max(np.abs(point))))
            out.append(best if best[1] <= self.MATCH_RTOL * scale else (None, best[1]))
        return out

    # ---------------------------------------------------------------- Rasterhash
    @property
    def grid_hash(self):
        """Hash von Fluid, Einheiten-/Achsenvertrag und allen referenzierten Rasterarrays (lazy, nur bei Bedarf)."""
        if self._grid_hash is None:
            h = hashlib.sha1()
            h.update(f"{self.fluid}|{self.diagram_id}|{self.adapter.axes}|{self.adapter.pressure_axis(self.daten)}".encode())
            for name in sorted(self.surfaces):
                s = self.surfaces[name]
                h.update(name.encode())
                h.update(str(s.array.shape).encode())
                h.update(np.ascontiguousarray(s.array).tobytes())
            self._grid_hash = h.hexdigest()
        return self._grid_hash


# ============================================================ optionale HDF5-Datei

class StateDiskCache:
    """Eine einzige stabile HDF5-Datei; Gruppen je (Fluid, Rasterhash). Wird bei Folgelaeufen aktualisiert.

    Gueltigkeit: Schemaversion (Datei), CoolProp-Version und Einheiten (Gruppe). Abweichende Gruppen werden
    verworfen (nie falsch zugeordnet). h5py ist optional: fehlt es, bleibt der Cache inaktiv.
    """

    def __init__(self, path, coolprop_version):
        self.path = path
        self.coolprop_version = str(coolprop_version)
        self.available = False
        try:
            import h5py  # noqa: F401
            self.available = True
        except Exception:
            LOGGER.warning("h5py nicht verfuegbar: Festplatten-Zustandscache deaktiviert")

    def _open(self, mode):
        import h5py
        import os
        os.makedirs(os.path.dirname(os.path.abspath(self.path)) or ".", exist_ok=True)
        f = h5py.File(self.path, mode)
        if "schema_version" not in f.attrs:
            f.attrs["schema_version"] = SCHEMA_VERSION
        elif int(f.attrs["schema_version"]) != SCHEMA_VERSION:
            for name in list(f.keys()):
                del f[name]
            f.attrs["schema_version"] = SCHEMA_VERSION
        return f

    @staticmethod
    def _group_name(fluid, grid_hash):
        return f"{fluid}/{grid_hash}"

    def load(self, fluid, grid_hash):
        """Gueltige Zustaende der Gruppe als dict id -> ThermoState (leer bei fehlender/veralteter Gruppe)."""
        if not self.available:
            return {}
        import os
        if not os.path.exists(self.path):
            return {}
        with self._open("a") as f:
            name = self._group_name(fluid, grid_hash)
            if name not in f:
                return {}
            g = f[name]
            if g.attrs.get("coolprop_version") != self.coolprop_version or g.attrs.get("units") != repr(sorted(UNITS.items())):
                del f[name]
                return {}
            states = {}
            for text in g["states"].asstr()[:]:
                st = ThermoState.from_json(text)
                states[st.id] = st
            return states

    def save(self, fluid, grid_hash, states):
        """Schreibt/aktualisiert die Gruppe in derselben Datei."""
        if not self.available or not states:
            return False
        import h5py
        with self._open("a") as f:
            name = self._group_name(fluid, grid_hash)
            if name in f:
                del f[name]
            g = f.create_group(name)
            g.attrs["coolprop_version"] = self.coolprop_version
            g.attrs["units"] = repr(sorted(UNITS.items()))
            g.attrs["fluid"] = fluid
            g.attrs["grid_hash"] = grid_hash
            g.attrs["source"] = "HEOS (ThermoSurface)"
            g.create_dataset("states", data=[s.to_json() for s in states.values()], dtype=h5py.string_dtype())
        return True


# ============================================================ Speicher

@dataclass
class StateStoreConfig:
    max_ram_states: int = 5000          # RAM-Budget (Anzahl Zustandsdatensaetze)
    disk_enabled: bool = False
    disk_path: str = "cache/thermosurface_states.h5"
    disk_threshold: int = 1000          # ab dieser Anzahl je Gruppe wird in die Datei geschrieben


class StateStore:
    """Zentrale Verwaltung: Diagrammindizes, Zustandscache (RAM-LRU, optional Datei), Zaehler."""

    def __init__(self, config=None, coolprop_version="n/a"):
        self.config = config or StateStoreConfig()
        self._indices = {}
        self._states = OrderedDict()
        self._state_groups = {}          # state id -> (fluid, grid_hash)
        self._disk = None
        self._disk_loaded = set()
        self.coolprop_version = coolprop_version
        self.evaluations = 0             # tatsaechliche HEOS-Auswertungen (vom Resolver erhoeht)
        self.cache_hits = 0
        self.disk_reads = 0
        self.disk_writes = 0

    # ---------------------------------------------------------------- Indizes
    def register(self, fluid, diagram_id, daten, key=None):
        """Registriert ein Diagrammergebnis (nur Referenzen, keine Auswertung, kein Plattenzugriff)."""
        index = DiagramIndex(fluid, diagram_id, daten)
        self._indices[(fluid, diagram_id)] = (key, index)
        return index

    def index(self, fluid, diagram_id):
        try:
            return self._indices[(fluid, diagram_id)][1]
        except KeyError:
            raise KeyError(f"Diagramm {diagram_id} fuer {fluid} nicht registriert") from None

    def has(self, fluid, diagram_id):
        return (fluid, diagram_id) in self._indices

    def ensure(self, fluid, diagram_id, daten, key=None):
        """Index zum aktuellen Ergebnisobjekt (neu registriert, falls sich `daten` geaendert hat)."""
        entry = self._indices.get((fluid, diagram_id))
        if entry is not None and entry[1].daten is daten:
            return entry[1]
        return self.register(fluid, diagram_id, daten, key)

    def clear_indices(self):
        self._indices.clear()

    # ---------------------------------------------------------------- Zustandscache
    def _disk_cache(self):
        if self._disk is None and self.config.disk_enabled:
            self._disk = StateDiskCache(self.config.disk_path, self.coolprop_version)
        return self._disk

    def get(self, state_id, fluid=None, grid_hash=None):
        """Zustand aus RAM; bei aktivem Datei-Cache einmalig je Gruppe aus der Datei nachladen."""
        state = self._states.get(state_id)
        if state is not None:
            self._states.move_to_end(state_id)
            self.cache_hits += 1
            return state
        disk = self._disk_cache()
        if disk is not None and disk.available and fluid is not None and (fluid, grid_hash) not in self._disk_loaded:
            self._disk_loaded.add((fluid, grid_hash))
            loaded = disk.load(fluid, grid_hash)
            self.disk_reads += 1
            for sid, st in loaded.items():
                self._remember(sid, st, fluid, grid_hash)
            state = self._states.get(state_id)
            if state is not None:
                self.cache_hits += 1
        return state

    def _remember(self, state_id, state, fluid, grid_hash):
        self._states[state_id] = state
        self._states.move_to_end(state_id)
        self._state_groups[state_id] = (fluid, grid_hash)
        while len(self._states) > self.config.max_ram_states:
            old, _ = self._states.popitem(last=False)
            self._state_groups.pop(old, None)

    def put(self, state, grid_hash=None):
        self._remember(state.id, state, state.fluid, grid_hash)
        disk = self._disk_cache()
        if disk is not None and disk.available:
            group = [s for sid, s in self._states.items() if self._state_groups.get(sid) == (state.fluid, grid_hash)]
            if len(group) >= self.config.disk_threshold:
                self.flush()
        return state

    def flush(self):
        """Aktualisiert die einzige Cachedatei (nur bei aktivem Datei-Cache)."""
        disk = self._disk_cache()
        if disk is None or not disk.available:
            return 0
        groups = {}
        for sid, st in self._states.items():
            groups.setdefault(self._state_groups.get(sid, (st.fluid, None)), {})[sid] = st
        written = 0
        for (fluid, grid_hash), states in groups.items():
            if grid_hash is None:
                continue
            existing = disk.load(fluid, grid_hash)
            existing.update(states)
            if disk.save(fluid, grid_hash, existing):
                written += len(states)
                self.disk_writes += 1
        return written

    def __len__(self):
        return len(self._states)

    def stats(self):
        return {"states": len(self._states), "evaluations": self.evaluations, "cache_hits": self.cache_hits,
                "disk_reads": self.disk_reads, "disk_writes": self.disk_writes,
                "indices": len(self._indices)}
