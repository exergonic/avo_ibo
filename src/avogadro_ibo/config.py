"""Persistent configuration for avo_ibo (method, basis, etc.)."""

import json
from . import CALCS_DIR

CONFIG_PATH = CALCS_DIR / "config.json"

_DEFAULT_CONFIG = {
    "method": "wb97x-d",
    "basis": "def2-TZVP",
    "iboview_style": False,
}

METHODS = ["hf", "b3lyp", "pbe", "pbe0", "wb97x-d", "mn15-l", "m06-2x"]
BASIS_SETS = [
    "cc-pVDZ",
    "cc-pVTZ",
    "aug-cc-pVDZ",
    "aug-cc-pVTZ",
    "def2-SVP",
    "def2-TZVP",
    "def2-SVPD",
    "def2-TZVPD",
    "6-31G(d,p)",
    "6-311G(d,p)",
    "6-31+G(d,p)",
    "6-311+G(d,p)",
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


def get_config_options():
    config = load_config()
    method_default = config.get("method", "hf")
    basis_default = config.get("basis", "cc-pVDZ")
    return {
        "method": {
            "type": "stringList",
            "label": "SCF Method",
            "values": METHODS,
            "default": METHODS.index(method_default)
            if method_default in METHODS
            else METHODS.index("hf"),
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
        "charge": {
            "type": "integer",
            "label": "Charge",
            "minimum": -10,
            "default": config.get("charge", 0),
            "order": 3.0,
        },
        "mult": {
            "type": "integer",
            "label": "Spin Multiplicity",
            "minimum": 1,
            "maximum": 1,
            "default": config.get("mult", 1),
            "order": 4.0,
        },
        "iboview_style": {
            "type": "boolean",
            "label": "IboView-like isosurface (truncate repolarization tails)",
            "default": config.get("iboview_style", False),
            "order": 5.0,
        },
        "memory_note": {
            "type": "text",
            "label": "Note",
            "default": (
                "\nRecommended presets:\n"
                "  General use (recommended):\t\twB97x-D / def2-TZVP\n"
                "  Small / quick preview:\t\tHF / cc-pVDZ\n"
                "  Charged / anions:\t\t\twB97x-D / aug-cc-pVDZ\n"
                "  Transition metals:\t\t\tMN15-L / def2-TZVP\n"
                "\n"
                "Only singlet spin multiplicities are currently supported.\n"
                "Calculations with higher spin multiplicities are not yet supported.\n"
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
    for key in ("method", "basis", "charge", "mult", "iboview_style"):
        if key in options:
            config[key] = options[key]
            changed = True
    if changed:
        save_config(config)
    return {"cjson": avo_input.get("cjson", {})}
