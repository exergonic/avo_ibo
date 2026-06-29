import argparse
import json
import sys

def main():
    parser = argparse.ArgumentParser(description="Avogadro Intrinsic Bond Orbital plugin")
    parser.add_argument("feature", nargs="?", default="ibo", help="Feature to run")
    parser.add_argument("--lang", default="en_US", help="Language")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    args = parser.parse_args()
    raw = sys.stdin.read()
    data = json.loads(raw)
    cjson = data.get("cjson", {})
    options = data.get("options", {})
    charge = data.get("charge", 0)
    spin = data.get("spin", 1)
    if args.feature == "ibo":
        from .calcs import compute_ibo
        result = compute_ibo(cjson, options, charge, spin, debug=args.debug)
    else:
        result = {"error": f"Unknown feature: {args.feature}"}
    json.dump(result, sys.stdout, indent=2)
    sys.stdout.flush()

if __name__ == "__main__":
    main()