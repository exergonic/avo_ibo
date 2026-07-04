"""
Thin Psi4 glue layer — delegates all IAO/IBO logic to ``pyibo``.

BSD 3-Clause License
Copyright (c) 2025-2026, Billy Wayne McCann
SPDX-License-Identifier: BSD-3-Clause
"""

import numpy as np


# -- helpers ----------------------------------------------------------------

def _option(options, key, default):
    v = options.get(key, default)
    if isinstance(v, str):
        return v.strip()
    return v


def _molden_header(wfn):
    """Return everything before [MO] from Psi4 Molden output."""
    import tempfile
    from pathlib import Path
    tmp = Path(tempfile.mktemp(suffix=".molden"))
    psi4_molden = __import__("psi4").molden
    psi4_molden(wfn, str(tmp))
    text = tmp.read_text(encoding="utf-8")
    tmp.unlink(missing_ok=True)
    idx = text.find("[MO]")
    if idx == -1:
        raise RuntimeError("Psi4 Molden output has no [MO] section")
    return text[:idx]


def _ao_perm_and_scale(wfn):
    """
    Return (perm, scale) to convert Psi4-internal AO ordering to
    Molden-standard ordering, including normalization for Cartesian
    d and f functions.
    """
    D_PERM = [0, 3, 5, 1, 2, 4]
    D_NORM = [1.0, 1.0, 1.0, 1.0 / np.sqrt(3.0),
              1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)]
    F_PERM = [0, 6, 9, 3, 1, 2, 5, 8, 7, 4]
    F_NORM = [
        1.0, 1.0, 1.0,
        1.0 / np.sqrt(5.0), 1.0 / np.sqrt(5.0), 1.0 / np.sqrt(5.0),
        1.0 / np.sqrt(5.0), 1.0 / np.sqrt(5.0), 1.0 / np.sqrt(5.0),
        1.0 / np.sqrt(15.0),
    ]

    n_AO = wfn.basisset().nbf()
    perm = np.arange(n_AO)
    scale = np.ones(n_AO)
    bas = wfn.basisset()
    ao = 0
    for sh in range(bas.nshell()):
        am = bas.shell(sh).am
        nf = bas.shell(sh).nfunction
        if am == 2:
            for i in range(6):
                perm[ao + i] = ao + D_PERM[i]
                scale[ao + i] = D_NORM[i]
        elif am == 3:
            for i in range(10):
                perm[ao + i] = ao + F_PERM[i]
                scale[ao + i] = F_NORM[i]
        ao += nf
    return perm, scale


# -- top-level entry point --------------------------------------------------

def compute_ibo(cjson, options, charge, spin, debug=False):
    import logging
    import psi4
    import pyibo
    from . import TEMP_DIR

    _psi_logger = logging.getLogger("psi4")
    _psi_logger.propagate = False
    _psi_logger.setLevel(logging.DEBUG if debug else logging.WARNING)
    _psi_handler = logging.FileHandler(str(TEMP_DIR / "psi4.log"), mode="w")
    _psi_handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    )
    _psi_logger.addHandler(_psi_handler)
    for _name in ["psi4.core", "psi4.driver"]:
        logging.getLogger(_name).setLevel(logging.WARNING)

    psi4.set_output_file(str(TEMP_DIR / "psi4.log"), append=True)

    # -- Parse input geometry and options
    atoms = cjson["atoms"]
    coords_raw = atoms["coords"]
    if isinstance(coords_raw, dict):
        coords = coords_raw["3d"]
    else:
        coords = coords_raw
    elem = atoms["elements"]["number"]
    charge_val = int(cjson.get("properties", {}).get("totalCharge", charge))
    spin_val = int(cjson.get("properties", {}).get("totalSpinMultiplicity", spin))
    ref = "uhf" if spin_val != 1 else "rhf"

    geom_lines = "\n".join(
        f"  {elem[i]:3d}  {coords[3*i]:12.8f}  {coords[3*i+1]:12.8f} "
        f"{coords[3*i+2]:12.8f}"
        for i in range(len(elem))
    )
    charge_tag = f"charge {charge_val}" if charge_val != 0 else ""
    mol_spec = f"{charge_tag}\n{geom_lines}\nno_com\nno_reorient\nsymmetry c1"
    mol = psi4.geometry(mol_spec)

    basis = _option(options, "basis", "cc-pVDZ")
    method = _option(options, "method", "hf")
    psi4.set_options({
        "basis": basis,
        "scf_type": "df",
        "reference": ref,
        "e_convergence": 1e-8,
        "d_convergence": 1e-8,
        "puream": 0,
    })
    energy, wfn = psi4.energy(method, return_wfn=True)

    # -- Extract numpy arrays from Psi4
    Ca = wfn.Ca().np
    nocc = wfn.doccpi()[0] + wfn.soccpi()[0]
    mints = psi4.core.MintsHelper(wfn.basisset())
    S_full = mints.ao_overlap().np
    bas_min = psi4.core.BasisSet.build(mol, "BASIS", "STO-3G", puream=0)
    S_min = mints.ao_overlap(bas_min, bas_min).np
    S12 = mints.ao_overlap(wfn.basisset(), bas_min).np
    F_AO = wfn.Fa().np

    # -- Build shell descriptions for the minimal basis
    shells = []
    for sh in range(bas_min.nshell()):
        shell = bas_min.shell(sh)
        shells.append((shell.ncenter, shell.am, shell.nfunction))

    # -- Run pyibo pipeline
    result = pyibo.compute_ibos(
        S_full, S12, S_min, Ca, nocc, F_AO, elem, shells,
        method=method, basis=basis, ref=ref,
    )

    # -- Write IAO-basis Molden
    header_text = _molden_header(wfn)
    C_AO_all = result["C_IAO"] @ result["C_IAO_all"]
    perm, scale = _ao_perm_and_scale(wfn)
    C_AO_molden = C_AO_all[perm, :] * scale[:, np.newaxis]
    pyibo.write_iao_molden(
        TEMP_DIR / "ibo.molden", header_text,
        C_AO_molden, result["occ_all"], result["energies_all"],
        result["n_orb"], result["n_ao"],
    )
    molden_text = (TEMP_DIR / "ibo.molden").read_text(encoding="utf-8")

    # -- Canonical Molden (reference)
    canon_path = TEMP_DIR / "canonical.molden"
    psi4.molden(wfn, str(canon_path))

    # -- Assemble analysis message
    msg = result["analysis_text"] + result["total_wiberg_text"]
    cjson["atoms"]["partialCharges"] = [round(c, 4) for c in result["net_charges"]]
    analysis_path = TEMP_DIR / "ibos.txt"
    analysis_path.write_text(msg, encoding="utf-8")

    return {
        "readProperties": True,
        "moleculeFormat": "molden",
        "molden": molden_text,
        "cjson": cjson,
        "message":
            "IBO analysis saved to calcs/last/ibos.txt\n"
            "Canonical MOs: canonical.molden",
    }
