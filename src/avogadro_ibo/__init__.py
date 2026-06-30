import argparse
import json
import logging
import shutil
import sys
import traceback
from pathlib import Path

# Make sure stdout stream is always Unicode, as Avogadro expects
sys.stdout.reconfigure(encoding="utf-8")

PLUGIN_DIR = Path(__file__).resolve().parent.parent.parent
CALCS_DIR = PLUGIN_DIR / "calcs"
TEMP_DIR = CALCS_DIR / "last"

logger = logging.getLogger(__name__)


def _prepare_calc_dir():
    """Ensure TEMP_DIR exists and clear its contents."""
    TEMP_DIR.mkdir(parents=True, exist_ok=True)
    for child in TEMP_DIR.iterdir():
        if child.is_file():
            child.unlink(missing_ok=True)
        elif child.is_dir():
            shutil.rmtree(child, ignore_errors=True)


def main():
    logging.basicConfig(
        level=logging.DEBUG,
        stream=sys.stderr,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    logger.debug("CLI input: " + " ".join(sys.argv))

    parser = argparse.ArgumentParser(
        description="Avogadro Intrinsic Bond Orbital plugin"
    )
    parser.add_argument("feature", nargs="?", default="ibo", help="Feature to run")
    parser.add_argument("--lang", default="en_US", help="Language")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    args = parser.parse_args()
    logger.debug(f"Parsed args: {args}")

    raw = sys.stdin.read()
    logger.debug(f"Read {len(raw)} bytes from stdin")
    data = json.loads(raw)
    logger.debug(f"Input keys: {', '.join(data.keys())}")

    cjson = data.get("cjson", {})
    options = data.get("options", {})
    charge = data.get("charge", 0)
    spin = data.get("spin", 1)

    try:
        if args.feature == "ibo":
            _prepare_calc_dir()
            from .calcs import compute_ibo

            result = compute_ibo(cjson, options, charge, spin, debug=args.debug)
        elif args.feature == "open":
            from .links import open_calcs_dir

            result = open_calcs_dir(cjson)
        elif args.feature == "test-molden":
            from .test_plugins import test_molden

            result = test_molden(cjson, options, charge, spin)
        elif args.feature == "test-energy":
            from .test_plugins import test_energy

            result = test_energy(cjson, options, charge, spin)
        else:
            result = {"error": f"Unknown feature: {args.feature}"}
    except Exception as e:
        limit = None if args.debug else 3
        result = {"error": "".join(traceback.format_exception(e, limit=limit))}
        logger.exception("Unhandled exception")

    print(json.dumps(result))


if __name__ == "__main__":
    main()
