"""Standalone CLI: python -m avogadro_ibo <file.xyz>"""

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
    if len(sys.argv) != 2:
        print("Usage: python -m avogadro_ibo <file.xyz>", file=sys.stderr)
        sys.exit(1)
    try:
        import psi4
    except ImportError:
        print("Error: Psi4 is required but not installed.", file=sys.stderr)
        print("Install it via pixi, or via conda (conda install psi4).", file=sys.stderr)
        sys.exit(1)
    coords, numbers = _parse_xyz(sys.argv[1])
    cjson = {
        "atoms": {"coords": {"3d": coords}, "elements": {"number": numbers}},
        "properties": {"totalCharge": 0, "totalSpinMultiplicity": 1},
    }
    result = compute_ibo(cjson, {}, charge=0, spin=1)
    print(result.get("message", "done"))

if __name__ == "__main__":
    main()
