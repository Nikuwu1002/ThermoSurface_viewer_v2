from pathlib import Path
import yaml

CONFIG_PATH = Path(__file__).resolve().parent / "config.yaml"

with CONFIG_PATH.open("r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

MAX_POINTS = config["limits"]["max_points"]
# Optionale Cache-Grenzen (Default, falls eine aeltere config.yaml sie nicht enthaelt).
CACHE_MAX_BYTES = int(float(config["limits"].get("cache_max_mb", 2048)) * 1024 ** 2)
CACHE_WARN_BYTES = int(float(config["limits"].get("cache_warn_mb", 1536)) * 1024 ** 2)

# Zustandscache des Vergleichsmoduls (optional in config.yaml; Datei-Cache standardmaessig aus).
_state_cache = config.get("state_cache") or {}
STATE_CACHE = {
    "max_ram_states": int(_state_cache.get("max_ram_states", 5000)),
    "disk_enabled": bool(_state_cache.get("disk_enabled", False)),
    "disk_path": str(CONFIG_PATH.parent / _state_cache.get("disk_path", "cache/thermosurface_states.h5")),
    "disk_threshold": int(_state_cache.get("disk_threshold", 1000)),
}

if (type(MAX_POINTS) is not int or MAX_POINTS < 2):
    raise ValueError("limits.max_points muss eine ganze Zahl >= 2 sein.")