"""Standalone CLI: python -m avogadro_ibo <file.xyz> [--method ...] [--basis ...] [--charge ...] [--spin ...]"""

import argparse
import sys
from pathlib import Path
from .calcs import compute_ibo

_ELEMENT_NUMBERS = {symbol: Z for Z, symbol in enumerate(
    ["X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
     "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
     "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
     "Ga", "Ge", "As", "Se", "Br", "Kr",
     "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
     "In", "Sn", "Sb", "Te", "I"], start=0)}


def _parse_xyz(path):
    lines = Path(path).read_text().strip().splitlines()
    natoms = int(lines[0].strip())
    coords, numbers = [], []
    for line in lines[2:2+natoms]:
        parts = line.strip().split()
        numbers.append(_ELEMENT_NUMBERS.get(parts[0].capitalize(), 0))
        coords.extend(map(float, parts[1:4]))
    return coords, numbers


def main():
    parser = argparse.ArgumentParser(
        description="Compute IBOs from an XYZ file"
    )
    parser.add_argument("xyz_file", help="Input XYZ file")
    parser.add_argument("--method", default=None,
                        help="SCF method (hf, b3lyp, pbe, ...)")
    parser.add_argument("--basis", default=None,
                        help="Basis set (cc-pVDZ, def2-TZVP, ...)")
    parser.add_argument("--charge", type=int, default=None,
                        help="Total charge (default: config or 0)")
    parser.add_argument("--spin", "--mult", type=int, default=None, dest="spin",
                        help="Spin multiplicity (default: config or 1)")
    parser.add_argument("--iboview-style", action="store_true",
                        help="Truncate repolarization tails (IboView-like isosurface)")
    args = parser.parse_args()

    try:
        import psi4
    except ImportError:
        print("Error: Psi4 is required but not installed.", file=sys.stderr)
        print("Install it via pixi, or via conda (conda install psi4).", file=sys.stderr)
        sys.exit(1)

    from .config import load_config
    _cfg = load_config()
    charge = args.charge if args.charge is not None else _cfg.get("charge", 0)
    spin = args.spin if args.spin is not None else _cfg.get("spin", 1)

    coords, numbers = _parse_xyz(args.xyz_file)
    mol_name = Path(args.xyz_file).stem
    cjson = {
        "name": mol_name,
        "atoms": {"coords": {"3d": coords}, "elements": {"number": numbers}},
        "properties": {"totalCharge": charge, "totalSpinMultiplicity": spin},
    }
    options = {}
    if args.method:
        options["method"] = args.method
    if args.basis:
        options["basis"] = args.basis
    if args.iboview_style:
        options["iboview_style"] = True
    result = compute_ibo(cjson, options, charge=charge, spin=spin)
    print(result.get("message", "done"))


if __name__ == "__main__":
    main()
