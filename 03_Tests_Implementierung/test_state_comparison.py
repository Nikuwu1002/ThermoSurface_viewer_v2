"""TS-B01..B07: Zustandsmodell, Store, Picking, Rekonstruktion, Auswahl, Vergleichsgeometrie.

Rekonstruierte Thermodynamik wird quantitativ gegen direkte, unabhaengige CoolProp-Aufrufe geprueft
(Reinstoffe Water und CO2). Nur der abschliessende GUI-Test benoetigt ein OpenGL-Fenster.
"""
import os
import tempfile
import time
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np

import thermo_state as ts
from selection_manager import FluidMismatch, NotAddable, SelectionManager

try:
    os.environ.setdefault("QT_API", "pyside6")
    import CoolProp
    import CoolProp.CoolProp as CP
    from CoolProp import AbstractState
    import main as app_main
    import state_picker
    from diagram_registry import DIAGRAMS
    from state_picker import PickError, barycentric_weights, pick_in_cell
    from state_resolver import StateResolver
    from state_store import (ADAPTERS, SINGLE_PHASE, NodeRef, StateDiskCache, StateStore, StateStoreConfig,
                             adapter_for, format_ref)
    from comparison_view import build_geometry, axis_scale
    CP_SKIP = None
except Exception as exc:  # pragma: no cover
    CP_SKIP = f"CoolProp/main nicht verfuegbar: {exc}"

try:
    import h5py  # noqa: F401
    H5_SKIP = None
except Exception:
    H5_SKIP = "h5py nicht installiert"

_DATA = {}


def get_data(fluid, diagram_id):
    key = (fluid, diagram_id)
    if key not in _DATA:
        _, raster, visu = app_main.lade_konfiguration()
        if diagram_id == "pvt":
            raster, visu, _meta = app_main.generate_stoff_grid(fluid, raster, visu)
        _DATA[key] = app_main.berechne_auswahl(fluid, raster, visu, diagram_type=diagram_id)
    return _DATA[key]


def direct_state(fluid):
    return AbstractState("HEOS", fluid)


def make(fluid, diagram_id, **store_kwargs):
    store = StateStore(StateStoreConfig(**store_kwargs) if store_kwargs else None, coolprop_version=CoolProp.__version__)
    index = store.register(fluid, diagram_id, get_data(fluid, diagram_id))
    return store, index, StateResolver(store)


def cell(index, surface, nodes, weights):
    """MeshPick aus drei benachbarten Knoten (Darstellungskoordinaten per Adapter, keine VTK-Geometrie noetig)."""
    s = index.surface(surface)
    tri = index.adapter.to_display(np.array([s.node(n) for n in nodes]), index.daten)
    pos = np.asarray(weights) @ tri
    names = [n for n, _ in index.adapter.surface_names]
    sub_index = names.index(surface) if surface in names else None
    return pick_in_cell(index, s.actor_key, sub_index, tri, pos, 0)


@unittest.skipIf(CP_SKIP, CP_SKIP or "")
class B01StateModel(unittest.TestCase):
    def test_mass_to_molar_conversion_is_J_per_kg_times_kg_per_mol(self):
        M = direct_state("Water").molar_mass()
        self.assertAlmostEqual(ts.mass_to_molar(1.0e6, M), 1.0e6 * M, places=9)
        self.assertAlmostEqual(ts.molar_to_mass(ts.mass_to_molar(2.5e5, M), M), 2.5e5, places=9)
        self.assertIsNone(ts.mass_to_molar(None, M))
        # HEOS: hmass * M == hmolar
        st = direct_state("Water")
        st.update(CP.PT_INPUTS, 1e6, 400.0)
        self.assertAlmostEqual(ts.mass_to_molar(st.hmass(), st.molar_mass()) / st.hmolar(), 1.0, places=12)

    def test_unknown_fields_are_none_not_zero_and_serializable(self):
        state = ts.ThermoState(id="TS-x", fluid="Water", source_diagram="pvt", point_id="p", T=300.0, p=1e5)
        self.assertIsNone(state.h_molar)
        self.assertIn("h_molar", state.missing_fields())
        again = ts.ThermoState.from_json(state.to_json())
        self.assertEqual(again, state)
        self.assertEqual(again.units["g_molar"], "J/mol")

    def test_synthetic_state_cannot_carry_heos_energies(self):
        with self.assertRaises(ValueError):
            ts.ThermoState(id="a", fluid="Water", source_diagram="pvt", point_id="a", eos_source=ts.SOURCE_SYNTHETIC,
                           validity=ts.SYNTHETIC_SUBSTITUTE, g_molar=1.0)
        with self.assertRaises(ValueError):
            ts.ThermoState(id="a", fluid="Water", source_diagram="pvt", point_id="a", validity="bogus")

    def test_state_ids_are_deterministic_and_unique(self):
        a = ts.make_state_id("Water", "pvt", "node", "surface:1,2")
        self.assertEqual(a, ts.make_state_id("Water", "pvt", "node", "surface:1,2"))
        self.assertNotEqual(a, ts.make_state_id("Water", "pvt", "node", "surface:1,3"))
        self.assertNotEqual(a, ts.make_state_id("CO2", "pvt", "node", "surface:1,2"))


@unittest.skipIf(CP_SKIP, CP_SKIP or "")
class B02Store(unittest.TestCase):
    def test_every_surface_point_has_a_unique_reference_without_copying_arrays(self):
        for diagram_id in DIAGRAMS.ids():
            store, index, _ = make("Water", diagram_id)
            refs = list(index.iter_refs())
            self.assertEqual(len(refs), len(set(refs)), diagram_id)
            self.assertEqual(len(refs), index.n_nodes)
            for name, surface in index.surfaces.items():
                source = get_data("Water", diagram_id)
                if name in source and isinstance(source[name], np.ndarray):
                    self.assertTrue(np.shares_memory(surface.array, source[name]), f"{diagram_id}:{name} kopiert")
            ref = refs[len(refs) // 2]
            np.testing.assert_array_equal(index.native(ref), index.surface(ref.surface).array[ref.index])
            self.assertEqual(store.evaluations, 0)  # Referenzieren verursacht keine Zustandsauswertung

    def test_repeated_query_costs_at_most_one_heos_evaluation(self):
        store, index, resolver = make("Water", "pht")
        ref = NodeRef("surface_liquid", (10, 10))
        first = resolver.resolve_node(index, ref)
        evaluations = store.evaluations
        self.assertGreaterEqual(evaluations, 1)
        second = resolver.resolve_node(index, ref)
        self.assertEqual(store.evaluations, evaluations)
        self.assertIs(first, second)
        self.assertEqual(store.cache_hits, 1)

    def test_inactive_comparison_does_nothing_expensive(self):
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "states.h5")
            store, index, _ = make("Water", "uvt", disk_enabled=True, disk_path=path)
            store.flush()
            self.assertEqual(store.evaluations, 0)
            self.assertFalse(os.path.exists(path))
            self.assertIsNone(index._grid_hash)   # Rasterhash wird nur bei Bedarf berechnet

    def test_grid_hash_changes_with_grid_and_fluid(self):
        _, a, _ = make("Water", "pht")
        _, b, _ = make("CO2", "pht")
        self.assertNotEqual(a.grid_hash, b.grid_hash)
        modified = dict(get_data("Water", "pht"))
        modified["surface_liquid"] = get_data("Water", "pht")["surface_liquid"] * 1.0001
        store = StateStore()
        self.assertNotEqual(a.grid_hash, store.register("Water", "pht", modified).grid_hash)

    def test_lru_budget(self):
        store = StateStore(StateStoreConfig(max_ram_states=2))
        for k in range(4):
            store.put(ts.ThermoState(id=f"s{k}", fluid="Water", source_diagram="pvt", point_id="x"))
        self.assertEqual(len(store), 2)
        self.assertIsNone(store.get("s0"))
        self.assertIsNotNone(store.get("s3"))

    @unittest.skipIf(H5_SKIP, H5_SKIP or "")
    def test_single_hdf5_file_is_updated_and_validated(self):
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "cache", "thermosurface_states.h5")
            config = dict(disk_enabled=True, disk_path=path, disk_threshold=10 ** 9)
            store, index, resolver = make("Water", "pht", **config)
            node = NodeRef("surface_vapor", (8, 8))
            state = resolver.resolve_node(index, node)
            self.assertEqual(store.flush(), 1)
            self.assertTrue(os.path.exists(path))
            size = os.path.getsize(path)

            # Neustart: gleiche Datei, Zustand ohne HEOS-Auswertung wieder da
            store2, index2, resolver2 = make("Water", "pht", **config)
            again = resolver2.resolve_node(index2, node)
            self.assertEqual(store2.evaluations, 0)
            self.assertEqual(again.h_molar, state.h_molar)
            self.assertEqual(store2.disk_reads, 1)
            # weiterer Punkt: dieselbe Datei wird aktualisiert, keine zweite angelegt
            resolver2.resolve_node(index2, NodeRef("surface_vapor", (9, 9)))
            store2.flush()
            self.assertEqual(sorted(os.listdir(os.path.dirname(path))), ["thermosurface_states.h5"])
            self.assertGreaterEqual(os.path.getsize(path), size)

            # geaendertes Raster -> anderer Hash -> keine falsch zugeordneten Zustaende
            modified = dict(get_data("Water", "pht"))
            modified["surface_vapor"] = get_data("Water", "pht")["surface_vapor"] * 1.001
            store3 = StateStore(StateStoreConfig(**config), coolprop_version=CoolProp.__version__)
            index3 = store3.register("Water", "pht", modified)
            StateResolver(store3).resolve_node(index3, node)
            self.assertGreaterEqual(store3.evaluations, 1)

            # andere CoolProp-Version -> Gruppe verworfen
            store4 = StateStore(StateStoreConfig(**config), coolprop_version="0.0-test")
            index4 = store4.register("Water", "pht", get_data("Water", "pht"))
            StateResolver(store4).resolve_node(index4, node)
            self.assertGreaterEqual(store4.evaluations, 1)


@unittest.skipIf(CP_SKIP, CP_SKIP or "")
class B03Picking(unittest.TestCase):
    def test_only_roundoff_scale_barycentric_overshoot_is_clamped(self):
        triangle = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        weights = barycentric_weights(triangle, np.array([-5e-7, 0.4, 0.0]))
        self.assertAlmostEqual(sum(weights), 1.0)
        self.assertTrue(all(weight >= 0.0 for weight in weights))
        with self.assertRaises(PickError):
            barycentric_weights(triangle, np.array([-1e-5, 0.4, 0.0]))

    def test_screen_near_fallback_stops_at_first_nearby_surface_hit(self):
        index = SimpleNamespace(diagram_id="pvt")
        error = PickError("Position liegt ausserhalb des Dreiecks")
        error.hit = ("surface", None, 3)
        hit = SimpleNamespace(diagram_id="pvt", actor_key="surface", sub_index=None, cell_id=3, message="",
                              display_xyz=(0.0, 0.0, 0.0))
        with patch.object(state_picker, "_pick_at_viewport",
                          side_effect=[error, hit]) as pick_at, \
             patch.object(state_picker, "_blocked_by_nonpickable_surface", return_value=False):
            result = state_picker.pick_from_viewport(object(), index, 100, 200)
        self.assertIs(result, hit)
        self.assertEqual(pick_at.call_args_list[1].args[2:], (99, 199))
        self.assertIn("fallback", result.message)

    def test_fallback_does_not_jump_to_a_different_surface_or_cell(self):
        error = PickError("Position liegt ausserhalb des Dreiecks")
        error.hit = ("surface", None, 3)
        other_layer = SimpleNamespace(diagram_id="pvt", actor_key="unstable_surface", sub_index=None,
                                      cell_id=3, message="")
        other_cell = SimpleNamespace(diagram_id="pvt", actor_key="surface", sub_index=None,
                                    cell_id=4, message="")
        same_cell = SimpleNamespace(diagram_id="pvt", actor_key="surface", sub_index=None,
                                    cell_id=3, message="", display_xyz=(0.0, 0.0, 0.0))
        with patch.object(state_picker, "_pick_at_viewport",
                          side_effect=[error, other_layer, other_cell, same_cell]), \
             patch.object(state_picker, "_blocked_by_nonpickable_surface", return_value=False):
            result = state_picker.pick_from_viewport(object(), SimpleNamespace(diagram_id="pvt"), 10, 10)
        self.assertIs(result, same_cell)

    def test_screen_near_fallback_does_not_override_non_surface_hits(self):
        with patch.object(state_picker, "_pick_at_viewport",
                          side_effect=PickError("Treffer liegt nicht auf einer pickbaren Zustandsflaeche")) as pick_at:
            with self.assertRaises(PickError):
                state_picker.pick_from_viewport(object(), SimpleNamespace(diagram_id="pvt"), 10, 10)
        pick_at.assert_called_once()

    def test_exact_grid_node_returns_unchanged_source_data(self):
        _, index, _ = make("Water", "pht")
        nodes = [(10, 10), (10, 11), (11, 10)]
        pick = cell(index, "surface_liquid", nodes, (1.0, 0.0, 0.0))
        self.assertTrue(pick.is_exact_node)
        self.assertEqual(pick.exact_ref, NodeRef("surface_liquid", (10, 10)))
        np.testing.assert_array_equal(index.native(pick.exact_ref), index.surface("surface_liquid").array[10, 10])
        self.assertEqual(pick.source_class, "heos")

    def test_interior_hit_has_valid_weights_and_phase_information(self):
        _, index, _ = make("Water", "pst")
        pick = cell(index, "surface_vapor", [(5, 5), (5, 6), (6, 5)], (0.2, 0.5, 0.3))
        self.assertAlmostEqual(sum(pick.weights), 1.0, places=12)
        self.assertTrue(all(0 <= w <= 1 for w in pick.weights))
        self.assertTrue(pick.fully_located)
        self.assertEqual(pick.physical_hint["vertex_phase"], ["vapor"] * 3)
        self.assertFalse(pick.is_exact_node)

    def test_substitute_surface_is_identified_as_synthetic(self):
        _, index, _ = make("Water", "pvt")
        pick = cell(index, "replacement_surface", [(10, 5), (10, 6), (11, 5)], (0.4, 0.3, 0.3))
        self.assertEqual(pick.source_class, "synthetic_substitute")
        self.assertEqual(pick.actor_key, "replacement_surface")

    def test_curves_and_unknown_actors_are_not_pickable(self):
        _, index, _ = make("Water", "pht")
        with self.assertRaises(PickError):
            pick_in_cell(index, "binodale", None, np.zeros((3, 3)), np.zeros(3))
        with self.assertRaises(PickError):
            pick_in_cell(index, "surface", 0, np.zeros((3, 3)), np.zeros(3))  # degeneriertes Dreieck

    def test_unmatched_vertices_are_reported_not_invented(self):
        _, index, resolver = make("Water", "pht")
        tri = np.array([[0.1, 0.2, 0.3], [0.2, 0.2, 0.3], [0.1, 0.3, 0.3]])
        pick = pick_in_cell(index, "surface", 0, tri, tri.mean(axis=0))
        self.assertFalse(pick.fully_located)
        state = resolver.resolve(pick, index)
        self.assertEqual(state.validity, ts.NUMERICAL_FAILURE)
        self.assertIsNone(state.T)


@unittest.skipIf(CP_SKIP, CP_SKIP or "")
class B04Reconstruction(unittest.TestCase):
    """Rekonstruktion gegen direkte HEOS-Referenzen (Water, CO2)."""

    FLUIDS = ("Water", "CO2")
    TOL = 1e-8

    def assert_matches_heos(self, state, T, p, label=""):
        ref = direct_state(state.fluid)
        ref.update(CP.PT_INPUTS, p, T)
        M = ref.molar_mass()
        for name, value in (("h_molar", ref.hmolar()), ("s_molar", ref.smolar()), ("u_molar", ref.umolar()),
                            ("g_molar", ref.gibbsmolar()), ("f_molar", ref.helmholtzmolar()),
                            ("rho_molar", ref.rhomolar())):
            self.assertAlmostEqual(getattr(state, name) / value, 1.0, delta=1e-7, msg=f"{label} {name}")
        self.assertAlmostEqual(state.v_molar * state.rho_molar, 1.0, places=12)
        self.assertEqual(state.units["h_molar"], "J/mol")
        self.assertGreater(M, 0)

    def test_exact_nodes_all_diagrams_reproduce_source_data(self):
        """B-T01: Gitterknoten -> Quelldaten, Einheit (molar), Herkunft, Fluid."""
        for fluid in self.FLUIDS:
            for diagram_id, surface, node in (("pht", "surface_liquid", (12, 8)), ("pht", "surface_vapor", (12, 8)),
                                              ("pht", "surface_high_temperature", (5, 20)), ("pst", "surface_vapor", (10, 10)),
                                              ("hsp", "surface_liquid", (10, 10)), ("uvt", "surface_high_temperature", (5, 30)),
                                              ("pvt", "surface", (20, 300))):
                with self.subTest(fluid=fluid, diagram=diagram_id, surface=surface):
                    store, index, resolver = make(fluid, diagram_id)
                    ref = NodeRef(surface, node)
                    state = resolver.resolve_node(index, ref)
                    native = dict(zip(index.adapter.axes, index.native(ref)))
                    M = direct_state(fluid).molar_mass()
                    self.assertEqual(state.fluid, fluid)
                    self.assertEqual(state.source_diagram, diagram_id)
                    self.assertEqual(state.eos_source, ts.SOURCE_HEOS)
                    self.assertIn(state.validity, (ts.VALID, ts.PARTIALLY_VALID), state.notes)
                    if "T" in native:
                        self.assertAlmostEqual(state.T / native["T"], 1.0, delta=1e-6)
                    if "p" in native:
                        self.assertAlmostEqual(state.p / native["p"], 1.0, delta=1e-6)
                    for axis, field_name in (("h", "h_molar"), ("s", "s_molar"), ("u", "u_molar")):
                        if axis in native:
                            self.assertAlmostEqual(getattr(state, field_name) / ts.mass_to_molar(native[axis], M), 1.0,
                                                   delta=1e-6, msg=axis)
                    if "v" in native:
                        self.assertAlmostEqual(state.v_molar / native["v"], 1.0, delta=1e-6)
                    # die pvt-Gleichgewichtsflaeche enthaelt das Plateau im Zweiphasengebiet: dort ist Q definiert
                    if state.phase == "twophase":
                        self.assertTrue(0.0 <= state.Q <= 1.0)
                    else:
                        self.assertIsNone(state.Q)
                    self.assertEqual(state.grid_ref["surface"], surface)
                    self.assert_matches_heos(state, state.T, state.p, diagram_id)

    def test_single_phase_interior_points_all_regions(self):
        """B-T02: Punkte zwischen Stuetzstellen: Rekonstruktion mit Plausibilitaetspruefung und Status."""
        cases = {"pht": ("surface_liquid", "surface_vapor", "surface_high_temperature"),
                 "pst": ("surface_liquid", "surface_vapor", "surface_high_temperature"),
                 "hsp": ("surface_liquid", "surface_vapor", "surface_high_temperature"),
                 "uvt": ("surface_liquid", "surface_vapor", "surface_high_temperature")}
        for fluid in self.FLUIDS:
            for diagram_id, surfaces in cases.items():
                store, index, resolver = make(fluid, diagram_id)
                for surface in surfaces:
                    with self.subTest(fluid=fluid, diagram=diagram_id, surface=surface):
                        pick = cell(index, surface, [(6, 7), (6, 8), (7, 7)], (0.3, 0.3, 0.4))
                        state = resolver.resolve(pick, index)
                        self.assertIn(state.validity, (ts.VALID, ts.PARTIALLY_VALID), state.notes)
                        self.assertIsNotNone(state.h_molar)
                        self.assertIsNotNone(state.g_molar)
                        # Plausibilisierung: HEOS-Zustand liegt in der Naehe des geometrisch interpolierten Punkts
                        target = dict(zip(index.adapter.axes, pick.physical_hint["native_interpolated"]))
                        M = direct_state(fluid).molar_mass()
                        extent = np.ptp(np.array(pick.physical_hint["vertex_native"]), axis=0)
                        heos = {"T": state.T, "p": state.p, "v": state.v_molar, "h": state.h_molar / M,
                                "s": state.s_molar / M, "u": state.u_molar / M}
                        for k, axis in enumerate(index.adapter.axes):
                            tol = 0.25 * extent[k] + 1e-5 * abs(target[axis])
                            self.assertLessEqual(abs(heos[axis] - target[axis]), tol + 1e-12, axis)

    def test_pvt_single_phase_interior_matches_direct_heos(self):
        for fluid in self.FLUIDS:
            with self.subTest(fluid=fluid):
                store, index, resolver = make(fluid, "pvt")
                pick = cell(index, "surface", [(10, 50), (10, 51), (11, 50)], (0.3, 0.3, 0.4))
                state = resolver.resolve(pick, index)
                T, v, _ = pick.physical_hint["native_interpolated"]
                ref = direct_state(fluid)
                ref.update(CP.DmolarT_INPUTS, 1.0 / v, T)
                self.assertAlmostEqual(state.p / ref.p(), 1.0, delta=1e-9)
                self.assertAlmostEqual(state.h_molar / ref.hmolar(), 1.0, delta=1e-9)

    def test_saturated_liquid_vapor_and_mixture_with_given_Q(self):
        """B-T03: Q = 0, 0.5, 1: p, v, h, s stimmen mit HEOS ueberein; Mischungsregel erfuellt."""
        for fluid in self.FLUIDS:
            for diagram_id in DIAGRAMS.ids():
                if diagram_id == "gpt":
                    continue   # g-p-T: Zweiphasengebiet = Saettigungslinie, keine Flaeche (siehe test_gpt_fvt.py)
                with self.subTest(fluid=fluid, diagram=diagram_id):
                    store, index, resolver = make(fluid, diagram_id)
                    two = index.surface("two_phase_surface")
                    qv = two.q_values
                    j0, j1 = int(np.argmin(qv)), int(np.argmax(qv))
                    self.assertEqual(qv[j0], 0.0)
                    self.assertEqual(qv[j1], 1.0)
                    i = 5
                    liquid = resolver.resolve_node(index, NodeRef("two_phase_surface", (i, j0)))
                    vapor = resolver.resolve_node(index, NodeRef("two_phase_surface", (i, j1)))
                    self.assertEqual((liquid.Q, vapor.Q), (0.0, 1.0))
                    j = int(np.searchsorted(qv, 0.5)) - 1
                    t = (0.5 - qv[j]) / (qv[j + 1] - qv[j])
                    pick = cell(index, "two_phase_surface", [(i, j), (i, j + 1), (i + 1, j)], (1 - t, t, 0.0))
                    mix = resolver.resolve(pick, index)
                    self.assertEqual(mix.validity, ts.VALID, mix.notes)
                    self.assertAlmostEqual(mix.Q, 0.5, places=12)
                    self.assertEqual(mix.phase, "twophase")
                    ref = direct_state(fluid)
                    for q, st in ((0.0, liquid), (0.5, mix), (1.0, vapor)):
                        ref.update(CP.QT_INPUTS, q, st.T)
                        self.assertAlmostEqual(st.p / ref.p(), 1.0, delta=1e-9)
                        self.assertAlmostEqual(st.v_molar / (1.0 / ref.rhomolar()), 1.0, delta=1e-9)
                        self.assertAlmostEqual(st.h_molar / ref.hmolar(), 1.0, delta=1e-9)
                        self.assertAlmostEqual(st.s_molar / ref.smolar(), 1.0, delta=1e-9)
                    # Mischungsregel bei gleicher Temperatur (gleicher Zeilenindex, andere Q-Spalte)
                    ref.update(CP.QT_INPUTS, 0.0, mix.T)
                    hl, sl, vl = ref.hmolar(), ref.smolar(), 1.0 / ref.rhomolar()
                    ref.update(CP.QT_INPUTS, 1.0, mix.T)
                    hv, sv, vv = ref.hmolar(), ref.smolar(), 1.0 / ref.rhomolar()
                    self.assertAlmostEqual(mix.h_molar, 0.5 * hl + 0.5 * hv, delta=1e-6 * abs(mix.h_molar) + 1e-6)
                    self.assertAlmostEqual(mix.v_molar / (0.5 * vl + 0.5 * vv), 1.0, delta=1e-9)
                    self.assertAlmostEqual(mix.s_molar, 0.5 * sl + 0.5 * sv, delta=1e-6 * abs(mix.s_molar) + 1e-6)

    def test_saturation_nodes_in_single_phase_surfaces_have_Q_0_and_1(self):
        _, index, resolver = make("Water", "pht")
        self.assertEqual(resolver.resolve_node(index, NodeRef("binodal_liquid", (10,))).Q, 0.0)
        self.assertEqual(resolver.resolve_node(index, NodeRef("binodal_vapor", (10,))).Q, 1.0)

    def test_T_and_p_alone_at_saturation_do_not_define_Q(self):
        _, index, resolver = make("Water", "pht")
        s = index.surface("binodal_liquid")
        T, h, p = s.node((10,))
        ref = direct_state("Water")
        ref.update(CP.QT_INPUTS, 0.5, T)
        self.assertAlmostEqual(ref.p() / p, 1.0, delta=1e-6)   # gleiche (T, p) fuer jedes Q
        ref.update(CP.QT_INPUTS, 0.9, T)
        self.assertAlmostEqual(ref.p() / p, 1.0, delta=1e-6)
        # Der Resolver nutzt fuer Zweiphasenpunkte (T,Q) und verlangt Q; der Saettigungsknoten hat festes Q
        self.assertEqual(index.surface("binodal_liquid").q_fixed, 0.0)

    def test_critical_fan_is_not_unique_and_critical_node_is_valid(self):
        for diagram_id in ("pht", "hsp"):
            with self.subTest(diagram=diagram_id):
                _, index, resolver = make("Water", diagram_id)
                two = index.surface("two_phase_surface")
                last = len(two.array) - 2
                crit = index.surface("critical_point")
                tri = index.adapter.to_display(np.array([two.node((last, 3)), two.node((last, 4)), crit.node((0,))]),
                                               index.daten)
                pick = pick_in_cell(index, "zweiphasengebiet", None, tri, np.array([0.3, 0.3, 0.4]) @ tri, 0)
                state = resolver.resolve(pick, index)
                self.assertEqual(state.validity, ts.NOT_UNIQUE)
                self.assertIsNone(state.Q)
                self.assertIsNone(state.h_molar)
                node = resolver.resolve_node(index, NodeRef("critical_point", (0,)))
                self.assertEqual(node.stability, "critical")
                self.assertIn(node.validity, (ts.VALID, ts.PARTIALLY_VALID))

    def test_metastable_and_unstable_homogeneous_states_are_marked(self):
        for fluid in self.FLUIDS:
            with self.subTest(fluid=fluid):
                _, index, resolver = make(fluid, "pvt")
                found = {}
                for name in ("metastable_liquid_surface", "metastable_vapor_surface", "unstable_surface"):
                    arr = index.surface(name).array
                    finite = np.isfinite(arr).all(axis=-1)
                    i, j = [(a, b) for a, b in zip(*np.nonzero(finite)) if 0 < b < arr.shape[1] - 2 and a + 1 < arr.shape[0]
                            and finite[a + 1, b] and finite[a, b + 1]][3]
                    pick = cell(index, name, [(i, j), (i, j + 1), (i + 1, j)], (0.4, 0.3, 0.3))
                    state = resolver.resolve(pick, index)
                    found[name] = state
                    self.assertEqual(state.eos_source, ts.SOURCE_HEOS_HOMOGENEOUS)
                    self.assertIn(state.stability, ("metastable", "unstable"))
                    self.assertEqual(state.native_input["flash"], "DmolarT")
                    self.assertIn("keine Gleichgewichtsinterpretation", state.notes[0])
                    self.assertIsNotNone(state.g_molar)
                self.assertEqual(found["metastable_liquid_surface"].native_input["phase_imposed"], "liquid")
                self.assertEqual(found["metastable_vapor_surface"].stability, "metastable")

    def test_substitute_surface_gets_no_invented_heos_values(self):
        """B-T04: Ersatzflaeche -> synthetisch markiert, keine erfundenen HEOS-Eigenschaften."""
        for fluid in self.FLUIDS:
            with self.subTest(fluid=fluid):
                _, index, resolver = make(fluid, "pvt")
                pick = cell(index, "replacement_surface", [(10, 5), (10, 6), (11, 5)], (0.4, 0.3, 0.3))
                state = resolver.resolve(pick, index)
                self.assertEqual(state.validity, ts.SYNTHETIC_SUBSTITUTE)
                self.assertEqual(state.eos_source, ts.SOURCE_SYNTHETIC)
                self.assertTrue(state.is_synthetic)
                for name in ("h_molar", "s_molar", "u_molar", "f_molar", "g_molar", "rho_molar", "Q"):
                    self.assertIsNone(getattr(state, name), name)
                self.assertIsNotNone(state.T)
                self.assertIsNotNone(state.v_molar)
                self.assertEqual(state.missing_fields().count("g_molar"), 1)
                # Der HEOS-Randknoten der Ersatzflaeche bleibt ein regulaerer, unveraenderter HEOS-Wert
                edge = resolver.resolve_node(index, NodeRef("replacement_surface", (10, 0)))
                self.assertEqual(edge.eos_source, ts.SOURCE_HEOS_HOMOGENEOUS)
                self.assertIsNotNone(edge.g_molar)

    def test_negative_pressure_is_kept_but_flagged(self):
        _, index, resolver = make("Water", "pvt")
        arr = index.surface("metastable_liquid_surface").array
        negative = np.argwhere(np.isfinite(arr[..., 2]) & (arr[..., 2] < 0))
        self.assertGreater(len(negative), 0, "Testdaten enthalten keine negativen Druecke")
        i, j = next((a, b) for a, b in negative if a + 1 < arr.shape[0] and b + 1 < arr.shape[1]
                    and np.isfinite(arr[a + 1, b]).all() and np.isfinite(arr[a, b + 1]).all())
        state = resolver.resolve_node(index, NodeRef("metastable_liquid_surface", (int(i), int(j))))
        self.assertLess(state.p, 0)
        self.assertEqual(state.validity, ts.PARTIALLY_VALID)
        self.assertTrue(any("negativ" in n for n in state.notes))

    def test_critical_nearness_is_flagged(self):
        _, index, resolver = make("Water", "uvt")
        node = resolver.resolve_node(index, NodeRef("critical_point", (0,)))
        self.assertEqual(node.phase, "critical_point")
        self.assertEqual(node.stability, "critical")
        self.assertTrue(any("kritischer Punkt" in n for n in node.notes))
        # Zustand sehr nahe am kritischen Punkt: kontrollierte Warnung statt Scheinpraezision
        state = resolver._state("Water")
        state.update(CP.DmolarT_INPUTS, state.rhomolar_critical() * 1.001, state.T_critical() * 1.0001)
        notes = []
        validity = resolver._warnings(state, resolver._read(state), ts.VALID, notes)
        self.assertEqual(validity, ts.PARTIALLY_VALID)
        self.assertTrue(any("kritische Naehe" in n for n in notes))

    def test_not_computable_point_gets_status_and_no_values(self):
        _, index, resolver = make("Water", "pht")

        class Broken:
            def update(self, *a, **k):
                raise RuntimeError("HEOS kaputt")

            def __getattr__(self, name):
                return lambda *a, **k: 1.0

        resolver._states["Water"] = Broken()
        pick = cell(index, "surface_liquid", [(6, 7), (6, 8), (7, 7)], (0.3, 0.3, 0.4))
        state = resolver.resolve(pick, index)
        self.assertEqual(state.validity, ts.NUMERICAL_FAILURE)
        self.assertIsNone(state.T)
        self.assertIsNone(state.h_molar)
        self.assertTrue(state.notes)

    def test_corrupted_dataset_is_reported_not_silently_fixed(self):
        data = dict(get_data("Water", "pht"))
        broken = np.array(data["surface_liquid"], copy=True)
        broken[7, 7] = np.nan
        data["surface_liquid"] = broken
        store = StateStore()
        index = store.register("Water", "pht", data)
        resolver = StateResolver(store)
        self.assertFalse(any(r.surface == "surface_liquid" and r.index == (7, 7) for r in index.iter_refs()))
        state = resolver.resolve_node(index, NodeRef("surface_liquid", (7, 7)))
        self.assertIn(state.validity, (ts.OUT_OF_RANGE, ts.NUMERICAL_FAILURE))
        self.assertIsNone(state.h_molar)

    def test_out_of_range_temperature_is_rejected(self):
        _, index, resolver = make("Water", "uvt")
        pick = cell(index, "surface_liquid", [(6, 7), (6, 8), (7, 7)], (0.3, 0.3, 0.4))
        pick.physical_hint["native_interpolated"] = (5.0, pick.physical_hint["native_interpolated"][1],
                                                     pick.physical_hint["native_interpolated"][2])
        state = resolver.resolve(pick, index)
        self.assertEqual(state.validity, ts.OUT_OF_RANGE)

    def test_every_diagram_has_an_adapter(self):
        self.assertEqual(set(ADAPTERS), set(DIAGRAMS.ids()))


class B05Selection(unittest.TestCase):
    def state(self, fluid="Water", diagram="pvt", validity=ts.VALID, **kw):
        return ts.ThermoState(id=f"id-{diagram}-{time.perf_counter_ns()}", fluid=fluid, source_diagram=diagram,
                              point_id="p", validity=validity, T=300.0, p=1e5, v_molar=1e-3, h_molar=1.0, **kw)

    def test_numbering_origin_color_add_delete_reset(self):
        m = SelectionManager()
        events = []
        m.subscribe(events.append)
        p1 = m.add(self.state(diagram="pvt"), (0, 0, 0), {"T [K]": 300.0})
        p2 = m.add(self.state(diagram="hsp"))
        p3 = m.add(self.state(diagram="pvt"))
        self.assertEqual([p.selection_index for p in m.points], [1, 2, 3])
        self.assertNotEqual(p1.color, p2.color)
        self.assertEqual(p1.color, p3.color)
        self.assertIn("aus pvt", p1.description())
        m.remove(2)
        p4 = m.add(self.state(diagram="pht"))
        self.assertEqual([p.selection_index for p in m.points], [1, 3, 4])     # Nummer 2 wird nie wiederverwendet
        self.assertEqual([p.selection_order for p in m.points], sorted(p.selection_order for p in m.points))
        self.assertEqual(len({p.selection_order for p in m.points}), 3)
        self.assertEqual(len(m.points_of("pvt")), 2)
        m.reset()
        self.assertEqual(len(m), 0)
        self.assertIsNone(m.fluid)
        self.assertEqual(m.add(self.state()).selection_index, 5)
        self.assertEqual(events[:3], ["added"] * 3)
        self.assertIn("removed", events)
        self.assertIn("reset", events)

    def test_unsubscribe_stops_lifecycle_notifications(self):
        manager = SelectionManager()
        events = []
        callback = events.append
        manager.subscribe(callback)
        manager.unsubscribe(callback)
        manager.add(self.state())
        self.assertEqual(events, [])

    def test_reorder_changes_only_connection_order(self):
        m = SelectionManager()
        for _ in range(3):
            m.add(self.state())
        self.assertTrue(m.move(3, -1))
        self.assertEqual([p.selection_index for p in m.points], [1, 3, 2])
        self.assertFalse(m.move(1, -1))
        self.assertEqual({p.selection_index for p in m.points}, {1, 2, 3})

    def test_only_same_fluid_and_only_usable_states(self):
        m = SelectionManager()
        m.add(self.state("Water"))
        with self.assertRaises(FluidMismatch):
            m.add(self.state("CO2"))
        with self.assertRaises(NotAddable):
            m.add(self.state(validity=ts.NUMERICAL_FAILURE))
        with self.assertRaises(NotAddable):
            m.add(self.state(validity=ts.NOT_UNIQUE))
        m.add(ts.ThermoState(id="syn", fluid="Water", source_diagram="pvt", point_id="x", T=1.0, p=1.0, v_molar=1.0,
                             eos_source=ts.SOURCE_SYNTHETIC, validity=ts.SYNTHETIC_SUBSTITUTE))
        self.assertEqual(len(m), 2)
        self.assertTrue(m.clear_for_fluid_change("CO2"))
        self.assertEqual(len(m), 0)

    def test_selection_is_not_persisted(self):
        import selection_manager
        source = open(selection_manager.__file__, encoding="utf-8").read()
        for forbidden in ("open(", "json.dump", "h5py", "pickle", "sqlite"):
            self.assertNotIn(forbidden, source)


class B06Geometry(unittest.TestCase):
    def point(self, m, diagram, **values):
        state = ts.ThermoState(id=f"s{len(m)}", fluid="Water", source_diagram=diagram, point_id="x", T=values.get("T", 300.0),
                               p=values.get("p", 1e5), v_molar=values.get("v", 1e-3), h_molar=values.get("h", 1000.0),
                               s_molar=values.get("s", 50.0), u_molar=values.get("u", 900.0),
                               f_molar=values.get("f", -5000.0), g_molar=values.get("g", -4900.0))
        return m.add(state, None, {"T [K]": state.T})

    def setUp(self):
        self.m = SelectionManager()
        for k, diagram in enumerate(("pvt", "hsp", "pht")):
            self.point(self.m, diagram, T=300 + 50 * k, p=1e5 * 10 ** k, v=1e-3 / 10 ** k, h=1000.0 * (k + 1), g=-4000.0 - k)

    def test_three_points_two_diagram_types_switch_axes(self):
        for axes in (("p", "v", "T"), ("h", "s", "p"), ("g", "T", "u")):
            geometry = build_geometry(self.m.points, axes)
            self.assertEqual(len(geometry.valid_rows), 3)
            self.assertEqual(geometry.segments, [(0, 1), (1, 2)])
            self.assertEqual([r.selection_index for r in geometry.rows], [1, 2, 3])
            self.assertEqual({r.color for r in geometry.rows}.__len__(), 3)
            self.assertEqual(geometry.warnings, [])

    def test_si_values_scales_and_units(self):
        geometry = build_geometry(self.m.points, ("p", "h", "T"))
        self.assertEqual(geometry.scales[0].mode, "log10")           # p > 0: logarithmische Druckskala
        self.assertIn("log10(p / Pa)", geometry.scales[0].title)
        self.assertAlmostEqual(geometry.rows[1].xyz[0], np.log10(1e6), places=12)
        self.assertAlmostEqual(geometry.rows[1].xyz[1], 2.0, places=12)   # h in kJ/mol: 2000 J/mol
        self.assertIn("kJ/mol", geometry.scales[1].title)
        self.assertAlmostEqual(geometry.rows[2].xyz[2], 400.0, places=12)  # T unveraendert in K

    def test_nonpositive_pressure_switches_to_linear_scale(self):
        scale = axis_scale("p", [1e5, -2e6, 3e6])
        self.assertEqual(scale.mode, "linear")
        self.assertIn("p ≤ 0", scale.title)
        self.assertAlmostEqual(scale.apply(-2e6), -2.0)

    def test_duplicate_axes_rejected(self):
        with self.assertRaises(ValueError):
            build_geometry(self.m.points, ("p", "p", "T"))
        with self.assertRaises(ValueError):
            build_geometry(self.m.points, ("p", "v", "x"))

    def test_invalid_point_breaks_the_line_and_warns(self):
        m = SelectionManager()
        self.point(m, "pvt")
        m.add(ts.ThermoState(id="syn", fluid="Water", source_diagram="pvt", point_id="x", T=350.0, p=2e5, v_molar=2e-3,
                             eos_source=ts.SOURCE_SYNTHETIC, validity=ts.SYNTHETIC_SUBSTITUTE))
        self.point(m, "pht")
        geometry = build_geometry(m.points, ("T", "h", "g"))
        self.assertEqual(len(geometry.valid_rows), 2)
        self.assertEqual(geometry.segments, [])                      # keine Verbindung ueber den ungueltigen Punkt
        self.assertEqual(len(geometry.warnings), 1)
        self.assertIn("nicht berechnet", geometry.warnings[0])
        full = build_geometry(m.points, ("T", "p", "v"))             # Ersatzpunkt hat T, p, v: Verbindung geschlossen
        self.assertEqual(full.segments, [(0, 1), (1, 2)])

    def test_axis_change_does_not_touch_states(self):
        before = [p.state.to_json() for p in self.m.points]
        build_geometry(self.m.points, ("p", "v", "T"))
        build_geometry(self.m.points, ("h", "s", "p"))
        self.assertEqual(before, [p.state.to_json() for p in self.m.points])

    def test_no_interpolated_surface_is_built(self):
        import comparison_view
        source = open(comparison_view.__file__, encoding="utf-8").read()
        for forbidden in ("delaunay", "StructuredGrid", "interpolate", "triangulate", "add_mesh(surf"):
            self.assertNotIn(forbidden, source)


if __name__ == "__main__":
    unittest.main()
