
"""
Datei: thermo_state.py
Projekt: ThermoSurface
Autor: Niklas King

Beschreibung:
Einheitliches Datenmodell zur Speicherung und Weitergabe
thermodynamischer Zustände unabhängig vom Diagrammtyp.

Aufgaben:
- Definition der ThermoState-Datenstruktur
- Speicherung thermodynamischer Zustandsgrößen in SI-Einheiten
- Verwaltung von Phasen-, Stabilitäts- und Gültigkeitsinformationen
- Kennzeichnung der Herkunft berechneter Zustände
- Umrechnung zwischen molaren und massenspezifischen Größen
- Erzeugung eindeutiger Zustandskennungen
- Serialisierung und Deserialisierung von Zustandsdaten

Abhängigkeiten:
- Python-Standardbibliothek: dataclasses, hashlib und json
- Keine weiteren projektspezifischen Module

Verwendet von:
- state_store.py / state_resolver.py
- selection_manager.py / comparison_view.py
- state_projection.py / gui.py
"""

import hashlib
import json
from dataclasses import dataclass, field, asdict, replace
from typing import Optional

# ---------------------------------------------------------------- Statuscodes (Datenvertrag Kap. 6)
VALID = "valid"
PARTIALLY_VALID = "partially_valid"
NOT_UNIQUE = "not_unique"
UNSUPPORTED_EOS = "unsupported_eos"
SYNTHETIC_SUBSTITUTE = "synthetic_substitute"
OUT_OF_RANGE = "out_of_range"
NUMERICAL_FAILURE = "numerical_failure"
STATUS_CODES = (VALID, PARTIALLY_VALID, NOT_UNIQUE, UNSUPPORTED_EOS, SYNTHETIC_SUBSTITUTE, OUT_OF_RANGE,
                NUMERICAL_FAILURE)
# Status, bei denen die HEOS-Eigenschaften vollstaendig und belastbar sind (ggf. mit Hinweis).
USABLE_STATUS = (VALID, PARTIALLY_VALID)

# ---------------------------------------------------------------- EOS-Herkunft
SOURCE_HEOS = "HEOS"                                  # stabiler Gleichgewichtszustand (HEOS)
SOURCE_HEOS_HOMOGENEOUS = "HEOS_homogeneous"          # homogene (metastabile/instabile) HEOS-Fortsetzung
SOURCE_SYNTHETIC = "synthetic_substitute"             # Modell-/Ersatzflaeche, keine EOS-Auswertung
EOS_SOURCES = (SOURCE_HEOS, SOURCE_HEOS_HOMOGENEOUS, SOURCE_SYNTHETIC)

SCHEMA_VERSION = 1
REFERENCE_STATE_NOTE = "CoolProp-Standardreferenzzustand des Fluids (unveraendert)"

# Einheiten der Quelldaten je Achse in den Diagrammen (Originalkoordinaten der Auswahl)
SOURCE_AXIS_UNITS = {"T": "K", "p": "Pa", "v": "m³/mol", "h": "J/kg", "s": "J/(kg K)", "u": "J/kg",
                     "g": "J/kg", "f": "J/kg"}

# Zustandsgroessen, die als Achsen des Vergleichsdiagramms verfuegbar sind: Feld -> (Label, Einheit)
QUANTITIES = {
    "T": ("T", "K"),
    "p": ("p", "Pa"),
    "v": ("v", "m³/mol"),
    "h": ("h", "J/mol"),
    "s": ("s", "J/(mol K)"),
    "u": ("u", "J/mol"),
    "f": ("f", "J/mol"),
    "g": ("g", "J/mol"),
}
QUANTITY_FIELD = {"T": "T", "p": "p", "v": "v_molar", "h": "h_molar", "s": "s_molar",
                  "u": "u_molar", "f": "f_molar", "g": "g_molar"}

UNITS = {"T": "K", "p": "Pa", "v_molar": "m3/mol", "h_molar": "J/mol", "s_molar": "J/(mol K)",
         "u_molar": "J/mol", "f_molar": "J/mol", "g_molar": "J/mol", "rho_molar": "mol/m3", "Q": "-"}

VALUE_FIELDS = ("T", "p", "v_molar", "h_molar", "s_molar", "u_molar", "f_molar", "g_molar", "rho_molar", "Q")


def mass_to_molar(value_per_kg, molar_mass_kg_per_mol):
    """Einmalige zentrale Umrechnung J/kg -> J/mol bzw. J/(kg K) -> J/(mol K): x[J/kg] * M[kg/mol]."""
    if value_per_kg is None:
        return None
    return float(value_per_kg) * float(molar_mass_kg_per_mol)


def molar_to_mass(value_per_mol, molar_mass_kg_per_mol):
    """Umkehrung (nur fuer Plausibilitaetsvergleiche mit massenbezogenen Quelldaten)."""
    if value_per_mol is None:
        return None
    return float(value_per_mol) / float(molar_mass_kg_per_mol)


def make_state_id(fluid, source_diagram, *parts):
    """Deterministische, eindeutige ID: Stoff + Diagramm + Herkunftsmerkmale (Gitterindizes/Gewichte/Eingaben)."""
    text = json.dumps([fluid, source_diagram, *parts], default=repr, sort_keys=True)
    return f"TS-{source_diagram}-{hashlib.sha1(text.encode('utf-8')).hexdigest()[:10]}"


@dataclass
class ThermoState:
    id: str
    fluid: str
    source_diagram: str
    point_id: str                         # stabile Rasterreferenz bzw. 'pick'
    grid_ref: Optional[dict] = None       # Raster-/Zellenreferenz (surface, indices, weights, cell)
    T: Optional[float] = None
    p: Optional[float] = None
    v_molar: Optional[float] = None
    h_molar: Optional[float] = None
    s_molar: Optional[float] = None
    u_molar: Optional[float] = None
    f_molar: Optional[float] = None
    g_molar: Optional[float] = None
    rho_molar: Optional[float] = None
    Q: Optional[float] = None
    phase: str = "nicht berechnet"
    stability: str = "nicht bestimmt"     # stable | metastable | unstable | saturated | critical | nicht bestimmt
    eos_source: str = SOURCE_HEOS
    validity: str = VALID
    native_input: dict = field(default_factory=dict)   # Eingaben des Flashs (Name -> Wert, SI bzw. Quelleinheit)
    units: dict = field(default_factory=lambda: dict(UNITS))
    provenance: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)
    schema_version: int = SCHEMA_VERSION

    def __post_init__(self):
        if self.validity not in STATUS_CODES:
            raise ValueError(f"Unbekannter Statuscode: {self.validity}")
        if self.eos_source not in EOS_SOURCES:
            raise ValueError(f"Unbekannte EOS-Herkunft: {self.eos_source}")
        if self.eos_source == SOURCE_SYNTHETIC:
            # Ersatzwerte duerfen nie als regulaere HEOS-Energien auftreten.
            for name in ("h_molar", "s_molar", "u_molar", "f_molar", "g_molar", "rho_molar", "Q"):
                if getattr(self, name) is not None:
                    raise ValueError(f"Synthetischer Ersatzzustand darf {name} nicht tragen")

    # ------------------------------------------------------------ Zugriff
    @property
    def is_usable(self):
        return self.validity in USABLE_STATUS

    @property
    def is_synthetic(self):
        return self.eos_source == SOURCE_SYNTHETIC

    def value(self, quantity):
        """Wert einer Achsengroesse ('T','p','v','h','s','u','f','g') oder None (= nicht berechnet)."""
        return getattr(self, QUANTITY_FIELD[quantity])

    def missing_fields(self):
        return [name for name in VALUE_FIELDS if getattr(self, name) is None]

    def with_changes(self, **changes):
        return replace(self, **changes)

    # ------------------------------------------------------------ Serialisierung
    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, data):
        data = dict(data)
        data.pop("schema_version", None)
        return cls(**data)

    def to_json(self):
        return json.dumps(self.to_dict(), sort_keys=True, default=repr)

    @classmethod
    def from_json(cls, text):
        return cls.from_dict(json.loads(text))
