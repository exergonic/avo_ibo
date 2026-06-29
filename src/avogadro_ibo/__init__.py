import argparse
import json
import logging
import sys
import traceback

# Make sure stdout stream is always Unicode, as Avogadro expects
sys.stdout.reconfigure(encoding="utf-8")

logger = logging.getLogger(__name__)


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
    logger.debug(f"Input keys: {list(data.keys())}")

    cjson = data.get("cjson", {})
    options = data.get("options", {})
    charge = data.get("charge", 0)
    spin = data.get("spin", 1)

    try:
        if args.feature == "ibo":
            from .calcs import compute_ibo

            result = compute_ibo(cjson, options, charge, spin, debug=args.debug)
        else:
            result = {"error": f"Unknown feature: {args.feature}"}
    except Exception as e:
        limit = None if args.debug else 3
        result = {"error": "".join(traceback.format_exception(e, limit=limit))}
        logger.exception("Unhandled exception")

    json.dump(result, sys.stdout, indent=2)
    sys.stdout.flush()
    logger.debug(f"Output keys: {list(result.keys())}")
    with open("_last_output.json", "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)


if __name__ == "__main__":
    main()
