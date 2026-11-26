
"""
Datei: state_resolver.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Rekonstruktion vollständiger thermodynamischer Zustände
aus geometrisch ausgewählten Punkten mittels CoolProp HEOS.

Aufgaben:
- Auswertung unabhängiger thermodynamischer Zustandsgrößen
- Durchführung diagrammspezifischer HEOS-Flashberechnungen
- Behandlung einphasiger und zweiphasiger Zustände
- Unterscheidung homogener und synthetischer Zustände
- Plausibilitätsprüfung der rekonstruierten Ergebnisse
- Nutzung des Zustandscaches zur Vermeidung erneuter Berechnungen

Abhängigkeiten:
- state_picker.py: Geometrische Trefferinformationen
- state_store.py: Rasterreferenzen und Zustandsspeicher
- thermo_state.py: Einheitliches Zustandsdatenmodell
- coolprop_heos_common.py: Gemeinsame HEOS-Funktionen
- CoolProp: Thermodynamische Zustandsrekonstruktion

Verwendet von:
- gui.py: Rekonstruktion angeklickter Zustandspunkte
"""

import logging

import numpy as np

import thermo_state as ts
from state_picker import (CLASS_HEOS, CLASS_HOMOGENEOUS, CLASS_SUBSTITUTE, CLASS_TWO_PHASE, EXACT_NODE_WEIGHT,
                          MeshPick, PickError)
from state_store import (CRITICAL, HOMOGENEOUS, SATURATION, SINGLE_PHASE, SUBSTITUTE, TWO_PHASE, DiagramIndex,
                         NodeRef, format_ref)

LOGGER = logging.getLogger("thermosurface.resolver")

# CoolProp-Phasenindex -> Text (iphase_*)
PHASE_NAMES = {0: "liquid", 1: "supercritical", 2: "supercritical_gas", 3: "supercritical_liquid",
               4: "critical_point", 5: "gas", 6: "twophase", 7: "unknown", 8: "not_imposed"}

# Unabhaengige Groessen je Diagramm fuer einphasige Quellflaechen: (Flash, Eingabeachsen)
SINGLE_PHASE_FLASH = {
    "pvt": ("DmolarT", ("v", "T")),
    "pht": ("HmassP", ("h", "p")),
    "pst": ("PSmass", ("s", "p")),
    "hsp": ("HmassP", ("h", "p")),
    "uvt": ("DmolarT", ("v", "T")),
    "gpt": ("PT", ("p", "T")),
    "fvt": ("DmolarT", ("v", "T")),
}

PLAUSIBILITY_EXTENT_FRACTION = 0.25    # erlaubte Abweichung relativ zur Ausdehnung des Dreiecks (nur Interpolation)
NODE_RTOL = 1e-6                       # exakter Knoten: relative Abweichung gegen die Quelldaten
CRITICAL_NEARNESS = (2e-3, 5e-2)       # |T-Tc|/Tc und |rho-rhoc|/rhoc fuer die Kritische-Naehe-Warnung


class StateResolver:
    """Rekonstruiert ThermoState-Datensaetze; erhoeht `store.evaluations` je HEOS-update."""

    def __init__(self, store, state_factory=None, cp_module=None):
        self.store = store
        self._cp = cp_module
        self._factory = state_factory
        self._states = {}
        self._molar_mass = {}

    # ---------------------------------------------------------------- HEOS-Zugriff
    @property
    def cp(self):
        if self._cp is None:
            import CoolProp.CoolProp as cp
            self._cp = cp
        return self._cp

    def _state(self, fluid):
        if fluid not in self._states:
            if self._factory is not None:
                self._states[fluid] = self._factory(fluid)
            else:
                from coolprop_heos_common import _new_state
                self._states[fluid] = _new_state(fluid)
        return self._states[fluid]

    def _coolprop_version(self):
        try:
            import CoolProp
            return CoolProp.__version__
        except Exception:
            return "n/a"

    def _update(self, state, flash, a, b, phase=None):
        """Genau ein HEOS-update (mit optionaler Phasenvorgabe); zaehlt als Auswertung."""
        cp = self.cp
        if flash == "PT":
            # (p, T): ohne Phase scheitert HEOS auf der Saettigungslinie (+-1e-4 %); der Aufrufer faellt dann
            # bewusst auf die Phase der Quellflaeche zurueck. pt_update_robust behandelt zusaetzlich den kritischen Rand.
            from coolprop_heos_common import pt_update_robust
            self.store.evaluations += 1
            pt_update_robust(state, b, a, phase)
            return
        inputs = {"DmolarT": cp.DmolarT_INPUTS, "HmassP": cp.HmassP_INPUTS, "PSmass": cp.PSmass_INPUTS,
                  "QT": cp.QT_INPUTS, "PQ": cp.PQ_INPUTS}[flash]
        self.store.evaluations += 1
        if phase is not None:
            state.specify_phase(phase)
        try:
            # CoolProp-Reihenfolge: DmolarT(rho,T), HmassP(h,p), PSmass(p,s), QT(Q,T), PQ(p,Q)
            if flash == "DmolarT":
                state.update(inputs, a, b)
            elif flash == "HmassP":
                state.update(inputs, a, b)
            elif flash == "PSmass":
                state.update(inputs, b, a)       # (s, p) -> (p, s)
            elif flash == "QT":
                state.update(inputs, a, b)
            elif flash == "PQ":
                state.update(inputs, b, a)       # (Q, p) -> (p, Q)
        finally:
            if phase is not None:
                state.unspecify_phase()

    def _read(self, state):
        """Alle molaren Groessen des aktuellen HEOS-Zustands."""
        rho = float(state.rhomolar())
        out = {
            "T": float(state.T()), "p": float(state.p()), "rho": rho, "v": 1.0 / rho,
            "h": float(state.hmolar()), "s": float(state.smolar()), "u": float(state.umolar()),
            "g": float(state.gibbsmolar()), "f": float(state.helmholtzmolar()),
            "molar_mass": float(state.molar_mass()),
            "phase_index": int(state.phase()), "Q": float(state.Q()),
        }
        if not all(np.isfinite(out[k]) for k in ("T", "p", "rho", "h", "s", "u", "g", "f")):
            raise ValueError("HEOS lieferte nicht endliche Werte")
        return out

    # ---------------------------------------------------------------- Zustandsaufbau
    def _build(self, pick_or_ref, fluid, diagram_id, props, *, validity=ts.VALID, eos_source=ts.SOURCE_HEOS,
               stability="stable", native_input=None, notes=None, grid_ref=None, provenance=None,
               state_id=None, point_id=None, phase_name=None, q_defined=False):
        phase_index = props["phase_index"]
        phase = phase_name or PHASE_NAMES.get(phase_index, "unknown")
        if q_defined or phase_index == 6:
            q = props["Q"] if 0.0 <= props["Q"] <= 1.0 else None
        else:
            q = None
        prov = {"backend": "HEOS", "coolprop_version": self._coolprop_version(),
                "reference_state": ts.REFERENCE_STATE_NOTE,
                "molar_mass_kg_per_mol": props["molar_mass"], **(provenance or {})}
        return ts.ThermoState(
            id=state_id, fluid=fluid, source_diagram=diagram_id, point_id=point_id or state_id, grid_ref=grid_ref,
            T=props["T"], p=props["p"], v_molar=props["v"], h_molar=props["h"], s_molar=props["s"],
            u_molar=props["u"], f_molar=props["f"], g_molar=props["g"], rho_molar=props["rho"], Q=q, phase=phase,
            stability=stability, eos_source=eos_source, validity=validity, native_input=native_input or {},
            provenance=prov, notes=list(notes or []))

    def _failure(self, fluid, diagram_id, state_id, point_id, status, message, grid_ref=None, provenance=None,
                 native_input=None):
        return ts.ThermoState(
            id=state_id, fluid=fluid, source_diagram=diagram_id, point_id=point_id, grid_ref=grid_ref,
            phase="nicht berechnet", validity=status, native_input=native_input or {},
            provenance={"backend": "HEOS", "coolprop_version": self._coolprop_version(), **(provenance or {})},
            notes=[message])

    # ---------------------------------------------------------------- Plausibilitaet / Warnungen
    def _plausibility(self, index, props, targets, extent, exact, M):
        """Vergleicht HEOS-Werte (in Quelleinheiten) mit Zielgroessen; Liste von Hinweisen bei Abweichung."""
        adapter = index.adapter
        critical = np.asarray(adapter.critical_native(index.daten), dtype=float)
        notes = []
        heos = {"T": props["T"], "v": props["v"], "p": props["p"],
                "h": props["h"] / M, "s": props["s"] / M, "u": props["u"] / M,
                "g": props["g"] / M, "f": props["f"] / M}
        for k, axis in enumerate(adapter.axes):
            target = targets[k]
            if target is None or not np.isfinite(target):
                continue
            scale = max(abs(target), abs(critical[k]), 1e-12)
            tol = NODE_RTOL * scale if exact else PLAUSIBILITY_EXTENT_FRACTION * abs(extent[k]) + NODE_RTOL * scale
            dev = abs(heos[axis] - target)
            if dev > tol:
                notes.append(f"Plausibilitaet: {axis} HEOS={heos[axis]:.6g} gegen "
                             f"{'Quelldaten' if exact else 'interpoliert'} {target:.6g} (Abweichung {dev:.3g} > {tol:.3g})")
        return notes

    def _warnings(self, state, props, validity, notes, hint_phase=None, allow_negative=False):
        """Phasenrand-, Kritische-Naehe- und Negativdruck-Hinweise."""
        if props["p"] <= 0 and not allow_negative:
            notes.append("p <= 0: nicht im stabilen Bereich; Wert aus HEOS ohne Gleichgewichtsinterpretation")
            validity = ts.PARTIALLY_VALID
        if props["p"] <= 0 and allow_negative:
            notes.append("negativer/verschwindender Druck der homogenen HEOS-Fortsetzung (metastabil/zugfest)")
            validity = ts.PARTIALLY_VALID
        Tc, rhoc = float(state.T_critical()), float(state.rhomolar_critical())
        if (abs(props["T"] - Tc) / Tc < CRITICAL_NEARNESS[0]) and (abs(props["rho"] - rhoc) / rhoc < CRITICAL_NEARNESS[1]):
            notes.append("kritische Naehe: Phasen- und Ableitungsangaben mit erhoehter Unsicherheit")
            validity = ts.PARTIALLY_VALID
        if hint_phase is not None:
            actual = PHASE_NAMES.get(props["phase_index"], "unknown")
            ok = {"liquid": ("liquid", "supercritical_liquid", "supercritical"),
                  "vapor": ("gas", "supercritical_gas", "supercritical"),
                  "supercritical": ("supercritical", "supercritical_gas", "supercritical_liquid", "gas", "liquid")}[hint_phase]
            if actual == "twophase":
                notes.append("Phasenrandnaehe: der geometrisch interpolierte Punkt liegt im Zweiphasengebiet "
                             f"(Q={props['Q']:.4g}); Zustand ueber die unabhaengigen Groessen eindeutig bestimmt")
                validity = ts.PARTIALLY_VALID
            elif actual not in ok:
                notes.append(f"Phasenrandnaehe: Quellflaeche {hint_phase}, HEOS-Phase {actual}")
                validity = ts.PARTIALLY_VALID
        return validity

    # ---------------------------------------------------------------- Oeffentliche API
    def resolve(self, pick: MeshPick, index: DiagramIndex):
        """Zustand zu einem Treffer (mit Cache ueber den Zustands-Store)."""
        ghash = index.grid_hash if self.store.config.disk_enabled else None
        if pick.is_exact_node and pick.exact_ref is not None:
            sid = ts.make_state_id(pick.fluid, pick.diagram_id, "node", format_ref(pick.exact_ref))
        else:
            sid = ts.make_state_id(pick.fluid, pick.diagram_id, pick.actor_key, pick.sub_index,
                                   [format_ref(r) if r else None for r in pick.source_point_ids],
                                   [round(w, 9) for w in pick.weights])
        cached = self.store.get(sid, pick.fluid, ghash)
        if cached is not None:
            return cached
        state = self._resolve_uncached(pick, index, sid)
        return self.store.put(state, ghash)

    def resolve_node(self, index: DiagramIndex, ref: NodeRef):
        """Zustand zu einem exakten Rasterknoten (unveraenderte Quelldaten)."""
        sid = ts.make_state_id(index.fluid, index.diagram_id, "node", format_ref(ref))
        ghash = index.grid_hash if self.store.config.disk_enabled else None
        cached = self.store.get(sid, index.fluid, ghash)
        if cached is not None:
            return cached
        state = self._resolve_node(index, ref, sid, pick=None)
        return self.store.put(state, ghash)

    # ---------------------------------------------------------------- interne Faelle
    def _resolve_uncached(self, pick, index, sid):
        fluid, diagram_id = pick.fluid, pick.diagram_id
        provenance = {"actor": pick.actor_key, "cell_id": pick.cell_id, "source_class": pick.source_class,
                      "source_points": pick.point_ids(), "weights": list(pick.weights),
                      "viewport_clipped": bool(getattr(pick, "clipped", False))}
        if not pick.fully_located:
            return self._failure(fluid, diagram_id, sid, "pick", ts.NUMERICAL_FAILURE,
                                 pick.message or "Quellknoten nicht zuordenbar", provenance=provenance)
        if pick.is_exact_node and pick.exact_ref is not None:
            return self._resolve_node(index, pick.exact_ref, sid, pick)
        try:
            if pick.source_class == CLASS_SUBSTITUTE:
                return self._synthetic(index, pick, sid, provenance)
            if pick.source_class == CLASS_TWO_PHASE:
                return self._two_phase(index, pick, sid, provenance)
            if pick.source_class == CLASS_HOMOGENEOUS:
                return self._homogeneous(index, pick, sid, provenance)
            return self._single_phase(index, pick, sid, provenance)
        except (ValueError, RuntimeError, OverflowError) as exc:
            return self._failure(fluid, diagram_id, sid, "pick", ts.NUMERICAL_FAILURE,
                                 f"HEOS-Auswertung fehlgeschlagen: {exc}", provenance=provenance)

    def _extent(self, pick):
        nodes = np.array([v for v in pick.physical_hint["vertex_native"]], dtype=float)
        return nodes.max(axis=0) - nodes.min(axis=0)

    def _targets_native(self, pick):
        return list(pick.physical_hint["native_interpolated"])

    @staticmethod
    def _aux_volume(pick):
        values = pick.physical_hint.get("vertex_aux_v")
        if not values or any(x is None for x in values):
            return None
        return float(np.dot(pick.weights, values))

    # --- einphasig (stabil), Flash aus unabhaengigen Quellgroessen
    def _single_phase(self, index, pick, sid, provenance, targets=None, exact=False, node_ref=None):
        fluid, diagram_id = pick.fluid, pick.diagram_id
        adapter = index.adapter
        state = self._state(fluid)
        M = float(state.molar_mass())
        flash, axes = SINGLE_PHASE_FLASH[diagram_id]
        targets = targets if targets is not None else self._targets_native(pick)
        by_axis = dict(zip(adapter.axes, targets))
        a_name, b_name = axes
        a = by_axis[a_name] if a_name != "v" else 1.0 / by_axis["v"]
        b = by_axis[b_name]
        if not (np.isfinite(a) and np.isfinite(b)) or (a_name == "v" and by_axis["v"] <= 0):
            return self._failure(fluid, diagram_id, sid, "pick", ts.OUT_OF_RANGE, "Unabhaengige Groessen nicht endlich/positiv",
                                 provenance=provenance)
        if b_name == "T" and not (state.Tmin() <= b <= state.Tmax()):
            return self._failure(fluid, diagram_id, sid, "pick", ts.OUT_OF_RANGE,
                                 f"T={b:.6g} K ausserhalb des HEOS-Bereichs [{state.Tmin():.4g}, {state.Tmax():.4g}] K",
                                 provenance=provenance)
        notes = []
        hint = self._phase_hint(pick)
        at_saturation = False
        if flash == "PT" and a <= 0:
            return self._failure(fluid, diagram_id, sid, "pick", ts.OUT_OF_RANGE,
                                 f"p={a:.6g} Pa <= 0: kein stabiler Gleichgewichtszustand", provenance=provenance)
        try:
            self._update(state, flash, a, b)
        except (ValueError, RuntimeError):
            if flash != "PT" or hint not in ("liquid", "vapor"):
                raise
            # Auf der Saettigungslinie (g-p-T) ist (T,p) allein kein eindeutiger Gleichgewichtszustand:
            # Zustand der Seite der Quellflaeche, Q ist aus (T,p,g) nicht bestimmt.
            cp = self.cp
            self._update(state, flash, a, b, phase=cp.iphase_liquid if hint == "liquid" else cp.iphase_gas)
            at_saturation = True
            notes.append(f"auf der Saettigungslinie: gesaettigter Zustand der {'Fluessigkeits' if hint == 'liquid' else 'Dampf'}"
                         "seite; Q ist aus (T, p, g) nicht bestimmbar (g_l = g_g)")
        props = self._read(state)
        validity = self._warnings(state, props, ts.VALID, notes, hint_phase=hint)
        extent = self._extent(pick) if not exact else [0.0] * 3
        notes += self._plausibility(index, props, targets, extent, exact, M)
        if any(n.startswith("Plausibilitaet") for n in notes):
            validity = ts.PARTIALLY_VALID
        stability = "saturated" if props["phase_index"] == 6 or at_saturation else "stable"
        return self._build(pick, fluid, diagram_id, props, validity=validity, stability=stability,
                           native_input={"flash": flash, a_name: float(a), b_name: float(b)}, notes=notes,
                           grid_ref=self._grid_ref(pick, node_ref), provenance=provenance, state_id=sid,
                           point_id=format_ref(node_ref) if node_ref else "pick")

    def _phase_hint(self, pick):
        phases = {p for p in pick.physical_hint.get("vertex_phase", []) if p}
        return next(iter(phases)) if len(phases) == 1 else None

    def _grid_ref(self, pick, node_ref=None):
        if node_ref is not None:
            return {"surface": node_ref.surface, "index": list(node_ref.index)}
        if pick is None:
            return None
        return {"actor": pick.actor_key, "cell_id": pick.cell_id, "nodes": pick.point_ids(),
                "weights": list(pick.weights)}

    # --- Zweiphasenflaeche
    def _two_phase(self, index, pick, sid, provenance):
        fluid, diagram_id = pick.fluid, pick.diagram_id
        adapter = index.adapter
        kinds = pick.physical_hint["vertex_kinds"]
        if CRITICAL in kinds:
            return self._failure(fluid, diagram_id, sid, "pick", ts.NOT_UNIQUE,
                                 "Faecher zum kritischen Punkt: Dampfgehalt Q ist dort nicht eindeutig; "
                                 "bitte einen Punkt im regulaeren Q-Netz oder den kritischen Knoten waehlen",
                                 provenance=provenance)
        qs = []
        for ref in pick.source_point_ids:
            surf = index.surface(ref.surface)
            qs.append(float(surf.q_values[ref.index[1]]))
        Q = float(np.dot(pick.weights, qs))
        param_axis = adapter.two_phase_axis
        param = pick.physical_hint["native_interpolated"][adapter.axes.index(param_axis)]
        return self._two_phase_state(index, pick, sid, provenance, param_axis, param, Q, exact=False)

    def _two_phase_state(self, index, pick, sid, provenance, param_axis, param, Q, exact, node_ref=None, targets=None):
        fluid, diagram_id = pick.fluid, pick.diagram_id
        state = self._state(fluid)
        M = float(state.molar_mass())
        if not (0.0 <= Q <= 1.0):
            return self._failure(fluid, diagram_id, sid, "pick", ts.OUT_OF_RANGE, f"Q={Q:.4g} ausserhalb [0,1]",
                                 provenance=provenance)
        if param_axis == "T":
            if not (state.Ttriple() < param < state.T_critical()):
                return self._failure(fluid, diagram_id, sid, "pick", ts.OUT_OF_RANGE,
                                     f"T={param:.6g} K nicht im Sattigungsbereich", provenance=provenance)
            self._update(state, "QT", Q, param)
        else:
            if not (0 < param < state.p_critical()):
                return self._failure(fluid, diagram_id, sid, "pick", ts.OUT_OF_RANGE,
                                     f"p={param:.6g} Pa nicht im Sattigungsbereich", provenance=provenance)
            self._update(state, "PQ", Q, param)
        props = self._read(state)
        props["Q"] = Q if props["Q"] < 0 else props["Q"]
        notes = []
        validity = self._warnings(state, props, ts.VALID, notes)
        targets = targets if targets is not None else list(pick.physical_hint["native_interpolated"])
        extent = [0.0] * 3 if exact else self._extent(pick)
        notes += self._plausibility(index, props, targets, extent, exact, M)
        if any(n.startswith("Plausibilitaet") for n in notes):
            validity = ts.PARTIALLY_VALID
        stability = "saturated" if Q in (0.0, 1.0) else "stable"
        if props["phase_index"] != 6 and Q not in (0.0, 1.0):
            notes.append("HEOS meldet keine Zweiphasen-Phase")
        return self._build(pick, fluid, diagram_id, props, validity=validity, stability=stability, q_defined=True,
                           phase_name="twophase" if 0.0 < Q < 1.0 else ("saturated_liquid" if Q == 0.0 else "saturated_vapor"),
                           native_input={"flash": "QT" if param_axis == "T" else "PQ", param_axis: float(param), "Q": float(Q)},
                           notes=notes, grid_ref=self._grid_ref(pick, node_ref), provenance=provenance, state_id=sid,
                           point_id=format_ref(node_ref) if node_ref else "pick")

    # --- homogene HEOS-Fortsetzung (metastabil/instabil)
    def _homogeneous(self, index, pick, sid, provenance, exact=False, node_ref=None, phase_name=None, targets=None):
        fluid, diagram_id = pick.fluid, pick.diagram_id
        state = self._state(fluid)
        M = float(state.molar_mass())
        cp = self.cp
        targets = targets if targets is not None else list(pick.physical_hint["native_interpolated"])
        by_axis = dict(zip(index.adapter.axes, targets))
        T = by_axis["T"]
        v = by_axis.get("v")
        if v is None:
            # g-p-T: die Quellspalten (T, p, g) enthalten das Volumen nicht; es steht je Knoten als Zusatzgroesse
            # der homogenen Ebene (Parametrisierung (T, v)) und wird wie die Gewichte lokal interpoliert.
            v = self._aux_volume(pick)
            if v is None:
                return self._failure(fluid, diagram_id, sid, "pick", ts.NUMERICAL_FAILURE,
                                     "Volumen der homogenen Quellknoten nicht verfuegbar", provenance=provenance)
        if not (np.isfinite(T) and np.isfinite(v) and T > 0 and v > 0):
            return self._failure(fluid, diagram_id, sid, "pick", ts.OUT_OF_RANGE, "Ungueltige Temperatur/Volumen",
                                 provenance=provenance)
        phase_name = phase_name or next((p for p in pick.physical_hint["vertex_phase"] if p), "vapor")
        phase = cp.iphase_liquid if phase_name == "liquid" else cp.iphase_gas
        self._update(state, "DmolarT", 1.0 / v, T, phase=phase)
        props = self._read(state)
        state.specify_phase(phase)
        try:
            slope = float(state.first_partial_deriv(cp.iP, cp.iDmolar, cp.iT))
        finally:
            state.unspecify_phase()
        stability = "metastable" if slope > 0 else "unstable"
        notes = [f"homogene HEOS-Fortsetzung (Phase {'fluessig' if phase_name == 'liquid' else 'gasfoermig'} vorgegeben), "
                 f"(dp/drho)_T={slope:.4g}: {stability}; keine Gleichgewichtsinterpretation"]
        validity = self._warnings(state, props, ts.VALID, notes, allow_negative=True)
        extent = [0.0] * 3 if exact else self._extent(pick)
        plaus = self._plausibility(index, props, targets, extent, exact, M)
        notes += plaus
        if plaus:
            validity = ts.PARTIALLY_VALID
        phase_label = "liquid_metastable_continuation" if phase_name == "liquid" else "gas_homogeneous_continuation"
        return self._build(pick, fluid, diagram_id, props, validity=validity, eos_source=ts.SOURCE_HEOS_HOMOGENEOUS,
                           stability=stability, phase_name=phase_label,
                           native_input={"flash": "DmolarT", "T": float(T), "v": float(v), "phase_imposed": phase_name},
                           notes=notes, grid_ref=self._grid_ref(pick, node_ref), provenance=provenance, state_id=sid,
                           point_id=format_ref(node_ref) if node_ref else "pick")

    # --- Ersatzflaeche
    def _synthetic(self, index, pick, sid, provenance, node_ref=None, native=None):
        values = native if native is not None else pick.physical_hint["native_interpolated"]
        by_axis = dict(zip(index.adapter.axes, values))
        return ts.ThermoState(
            id=sid, fluid=pick.fluid, source_diagram=pick.diagram_id,
            point_id=format_ref(node_ref) if node_ref else "pick", grid_ref=self._grid_ref(pick, node_ref),
            T=float(by_axis["T"]), p=float(by_axis["p"]), v_molar=float(by_axis["v"]),
            phase="synthetic_substitute", stability="nicht bestimmt", eos_source=ts.SOURCE_SYNTHETIC,
            validity=ts.SYNTHETIC_SUBSTITUTE,
            native_input={"T": float(by_axis["T"]), "v": float(by_axis["v"]), "p_model": float(by_axis["p"])},
            provenance={"backend": "none (Modell)", "coolprop_version": self._coolprop_version(), **provenance},
            notes=["Ersatzflaeche (Hermite-Modell zwischen den HEOS-Spinodalpunkten): p ist ein Modellwert, "
                   "keine EOS-Auswertung; h, s, u, f, g sind nicht berechnet"])

    # ---------------------------------------------------------------- exakter Knoten
    def _resolve_node(self, index, ref, sid, pick):
        surface = index.surface(ref.surface)
        native = index.native(ref)
        fluid, diagram_id = index.fluid, index.diagram_id
        adapter = index.adapter
        provenance = {"source_class": CLASS_HEOS, "node": format_ref(ref), "kind": surface.kind}
        if pick is not None:
            provenance.update({"actor": pick.actor_key, "cell_id": pick.cell_id, "weights": list(pick.weights)})
        fake = pick or MeshPick(diagram_id, surface.actor_key or "", -1, [ref] * 3, (1.0, 0.0, 0.0),
                                (0.0, 0.0, 0.0), fluid=fluid,
                                physical_hint={"vertex_phase": [surface.phase] * 3, "vertex_kinds": [surface.kind] * 3,
                                               "vertex_native": [native.tolist()] * 3,
                                               "vertex_aux_v": [surface.node_aux("v", ref.index)] * 3,
                                               "native_interpolated": tuple(native)})
        try:
            if surface.kind == CRITICAL:
                state = self._state(fluid)
                self._update(state, "DmolarT", float(state.rhomolar_critical()), float(state.T_critical()))
                props = self._read(state)
                props["phase_index"] = 4
                return self._build(fake, fluid, diagram_id, props, validity=ts.VALID, stability="critical",
                                   native_input={"flash": "DmolarT", "T": float(state.T_critical()),
                                                 "rho": float(state.rhomolar_critical())},
                                   notes=["kritischer Punkt (Rand der Zustandsflaechen)"], grid_ref=self._grid_ref(None, ref),
                                   provenance=provenance, state_id=sid, point_id=format_ref(ref), phase_name="critical_point")
            if surface.kind == SATURATION:
                axis = "T" if "T" in adapter.axes else "p"
                param = float(native[adapter.axes.index(axis)])
                return self._two_phase_state(index, fake, sid, provenance, axis, param, float(surface.q_fixed), True,
                                             node_ref=ref, targets=list(native))
            if surface.kind == TWO_PHASE:
                axis = adapter.two_phase_axis
                param = float(native[adapter.axes.index(axis)])
                return self._two_phase_state(index, fake, sid, provenance, axis, param,
                                             float(surface.q_values[ref.index[1]]), True, node_ref=ref, targets=list(native))
            if surface.kind == HOMOGENEOUS:
                phase_name = "liquid" if surface.phase == "liquid" else "vapor"
                provenance["source_class"] = CLASS_HOMOGENEOUS
                return self._homogeneous(index, fake, sid, provenance, exact=True, node_ref=ref, phase_name=phase_name,
                                         targets=list(native))
            if surface.kind == SUBSTITUTE:
                if surface.node_source(ref.index) == 1:
                    provenance["source_class"] = CLASS_SUBSTITUTE
                    return self._synthetic(index, fake, sid, provenance, node_ref=ref, native=native)
                # HEOS-Spinodalpunkt am Rand der Ersatzflaeche: unveraenderter HEOS-Wert
                provenance["source_class"] = CLASS_HOMOGENEOUS
                phase_name = "liquid" if ref.index[1] == 0 else "vapor"
                return self._homogeneous(index, fake, sid, provenance, exact=True, node_ref=ref, phase_name=phase_name,
                                         targets=list(native))
            return self._single_phase(index, fake, sid, provenance, targets=list(native), exact=True, node_ref=ref)
        except (ValueError, RuntimeError, OverflowError) as exc:
            return self._failure(fluid, diagram_id, sid, format_ref(ref), ts.NUMERICAL_FAILURE,
                                 f"HEOS-Auswertung fehlgeschlagen: {exc}", grid_ref=self._grid_ref(None, ref),
                                 provenance=provenance)
