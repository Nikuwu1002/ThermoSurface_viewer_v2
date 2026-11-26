"""
Datei: diagram_cache.py
Projekt: ThermoSurface
Kategorie: Hilfsmodul / Cache

Beschreibung:
Speichert berechnete Diagrammdatensätze in einem speicherbegrenzten RAM-Cache.

Aufgaben:
- Erzeugt stabile Cache-Schlüssel aus Stoff, Modell und Rechenparametern.
- Verdrängt ältere Einträge nach LRU-Prinzip und meldet Speichergrenzen.

Abhängigkeiten:
- Projektmodule (direkt): keine
- Bibliotheken: CoolProp, numpy

Verwendet von:
- Projektmodule: calculation_manager.py, gui.py
"""

import hashlib
import json
import logging
from collections import OrderedDict
from dataclasses import dataclass
from threading import RLock

import numpy as np

LOGGER = logging.getLogger("thermosurface.cache")

MODEL = "heos"
SCHEMA_VERSION = 1


def solver_version():
    try:
        import CoolProp
        return f"CoolProp {CoolProp.__version__}|schema {SCHEMA_VERSION}"
    except Exception:
        return f"CoolProp n/a|schema {SCHEMA_VERSION}"


def _canonical(value):
    if isinstance(value, dict):
        return {str(k): _canonical(v) for k, v in sorted(value.items(), key=lambda item: str(item[0]))}
    if isinstance(value, (list, tuple)):
        return [_canonical(v) for v in value]
    if isinstance(value, np.ndarray):
        return {"__ndarray__": hashlib.sha1(np.ascontiguousarray(value).tobytes()).hexdigest(),
                "shape": list(value.shape), "dtype": str(value.dtype)}
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, float):
        return repr(value)
    return value


def params_hash(params):
    text = json.dumps(_canonical(params), sort_keys=True, default=repr)
    return hashlib.sha1(text.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class DiagramKey:
    fluid: str
    diagram_id: str
    params_hash: str
    model: str = MODEL
    solver_version: str = ""

    @classmethod
    def build(cls, fluid, diagram_id, params):
        return cls(str(fluid), str(diagram_id), params_hash(params), MODEL, solver_version())


def estimate_nbytes(value, _seen=None):
    """Grobe Speicherabschaetzung (numpy-Arrays zaehlen exakt, sonst rekursiv/Objekt-__dict__)."""
    _seen = set() if _seen is None else _seen
    if id(value) in _seen:
        return 0
    _seen.add(id(value))
    if isinstance(value, np.ndarray):
        return int(value.nbytes)
    if isinstance(value, dict):
        return sum(estimate_nbytes(v, _seen) for v in value.values())
    if isinstance(value, (list, tuple, set)):
        return sum(estimate_nbytes(v, _seen) for v in value)
    if hasattr(value, "__dict__"):
        return estimate_nbytes(vars(value), _seen)
    return 64


class DiagramResultCache:
    """LRU-Cache. Bei Ueberschreitung des Budgets werden aelteste Eintraege verdraengt und gemeldet
    (nie stillschweigend); der zuletzt gespeicherte Eintrag bleibt immer erhalten."""

    def __init__(self, max_bytes=2 * 1024 ** 3, warn_bytes=None):
        self.max_bytes = int(max_bytes)
        self.warn_bytes = int(warn_bytes) if warn_bytes is not None else int(0.75 * self.max_bytes)
        self._items = OrderedDict()
        self._sizes = {}
        self._lock = RLock()
        self.hits = 0
        self.misses = 0
        self.evictions = []
        self.warnings = []

    def __len__(self):
        return len(self._items)

    def __contains__(self, key):
        return key in self._items

    @property
    def total_bytes(self):
        return sum(self._sizes.values())

    def get(self, key):
        with self._lock:
            if key in self._items:
                self._items.move_to_end(key)
                self.hits += 1
                return self._items[key]
            self.misses += 1
            return None

    def peek(self, key):
        return self._items.get(key)

    def put(self, key, daten):
        """Speichert `daten`; gibt Liste der dabei verdraengten Schluessel zurueck."""
        with self._lock:
            self._items[key] = daten
            self._items.move_to_end(key)
            self._sizes[key] = estimate_nbytes(daten)
            evicted = []
            while self.total_bytes > self.max_bytes and len(self._items) > 1:
                old, _ = next(iter(self._items.items()))
                if old == key:
                    break
                del self._items[old]
                self._sizes.pop(old, None)
                evicted.append(old)
            if evicted:
                self.evictions.extend(evicted)
                LOGGER.warning("Cache-Budget %.0f MB ueberschritten: %d Ergebnis(se) verdraengt (werden bei Bedarf neu berechnet)",
                               self.max_bytes / 1e6, len(evicted))
            if self.total_bytes > self.warn_bytes:
                message = (f"Cache belegt {self.total_bytes / 1e6:.0f} MB "
                           f"(Warnschwelle {self.warn_bytes / 1e6:.0f} MB, Budget {self.max_bytes / 1e6:.0f} MB)")
                self.warnings.append(message)
                LOGGER.warning(message)
            return evicted

    def clear(self):
        with self._lock:
            self._items.clear()
            self._sizes.clear()

    def stats(self):
        return {"entries": len(self._items), "bytes": self.total_bytes, "hits": self.hits,
                "misses": self.misses, "evictions": len(self.evictions)}
