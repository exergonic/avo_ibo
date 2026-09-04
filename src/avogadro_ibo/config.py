"""Persistent configuration for avo_ibo (method, basis, etc.)."""

import json
from pathlib import Path
from . import CALCS_DIR

CONFIG_PATH = CALCS_DIR / "config.json"

_DEFAULT_CONFIG = {
    "method": "wb97x-d",
    "basis": "def2-TZVP",
    "iboview_style": True,
    # Empty means the plugin default (PLUGIN_DIR/calcs).  The settings
    # file itself always stays there, so moving the output home never
    # orphans this setting.
    "calcs_dir": "",
}

METHODS = ["HF", "B3LYP", "PBE", "PBE0", "wB97X-D", "MN15-L", "M06-2X"]
BASIS_SETS = [
    "cc-pVDZ",
    "aug-cc-pVDZ",
    "cc-pVTZ",
    "def2-SVP",
    "def2-SVPD",
    "def2-TZVP",
]


def load_config():
    if not CONFIG_PATH.exists():
        return dict(_DEFAULT_CONFIG)
    try:
        return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return dict(_DEFAULT_CONFIG)


def save_config(config):
    CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    CONFIG_PATH.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")


def resolve_output_dir(explicit=None):
    """Root folder for per-molecule calc dirs.

    Precedence: per-call ``explicit`` > persisted ``calcs_dir`` > plugin
    default ``CALCS_DIR``.  Empty string means default.  Resolved
    per-call (not a mutated global) so every consumer — adapter,
    links, CLI — sees the same home.  The folder is created, parents
    included; failures raise RuntimeError with a human-readable message
    (the adapter turns it into plugin JSON ``{"error": ...}``).
    """
    raw = explicit if explicit else load_config().get("calcs_dir", "")
    root = CALCS_DIR if not raw else Path(raw).expanduser()
    try:
        root.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        raise RuntimeError(
            f"Could not create calculations folder '{root}': {e}. "
            f"Check the 'Run calculations in' path and permissions."
        )
    return root


def get_config_options():
    config = load_config()
    method_default = config.get("method", "HF")
    basis_default = config.get("basis", "cc-pVDZ")
    return {
        "method": {
            "type": "stringList",
            "label": "SCF Method",
            "values": METHODS,
            "default": METHODS.index(method_default)
            if method_default in METHODS
            else METHODS.index("HF"),
            "order": 1.0,
        },
        "basis": {
            "type": "stringList",
            "label": "Basis Set",
            "values": BASIS_SETS,
            "default": BASIS_SETS.index(basis_default)
            if basis_default in BASIS_SETS
            else BASIS_SETS.index("cc-pVDZ"),
            "order": 2.0,
        },
        "calcs_dir": {
            "type": "string",
            "label": "Run calculations in",
            "default": config.get("calcs_dir") or str(CALCS_DIR),
            "order": 2.5,
        },
        "iboview_style": {
            "type": "boolean",
            "label": "IboView-like isosurface",
            "default": config.get("iboview_style", True),
            "order": 3.0,
        },
        "memory_note": {
            "type": "text",
            "label": "Note",
            "default": (
                "\nRecommended presets:\n"
                "  General use (recommended):\t\twB97X-D / def2-TZVP\n"
                "  Small / quick preview:\t\tHF / cc-pVDZ\n"
                "  Charged / anions:\t\t\twB97x-D / aug-cc-pVDZ\n"
                "  Transition metals:\t\t\tMN15-L / def2-SVP\n"
                "\n"
                "Memory requirements increase with system size and basis set.\n"
                "Larger systems (30+ atoms) with triple-zeta or\n"
                "diffuse basis may require significant memory.\n"
                "Switch to a smaller basis if the calculation fails."
            ),
            "order": 99.0,
        },
    }


def update_config(avo_input):
    options = avo_input.get("options", {})
    config = load_config()
    changed = False
    for key in ("method", "basis", "iboview_style"):
        if key in options:
            config[key] = options[key]
            changed = True
    if "calcs_dir" in options:
        raw = (options["calcs_dir"] or "").strip()
        if raw == "":
            if config.pop("calcs_dir", None) is not None:
                changed = True
        else:
            probe = Path(raw).expanduser()
            try:
                probe.mkdir(parents=True, exist_ok=True)
            except OSError as e:
                return {"cjson": avo_input.get("cjson", {}),
                        "error": f"Could not create folder '{raw}': {e}."}
            # Store absolute: a relative path would later resolve against
            # Avogadro's unpredictable working directory.
            saved = str(probe.resolve())
            if config.get("calcs_dir") != saved:
                config["calcs_dir"] = saved
                changed = True
    if changed:
        save_config(config)
    return {"cjson": avo_input.get("cjson", {})}
