import tempfile
from pathlib import Path


def _option(options, key, default):
    v = options.get(key, default)
    if isinstance(v, str):
        return v.strip()
    return v


def compute_ibo(cjson, options, charge, spin, debug=False):
    import psi4
    import numpy as np

    log_path = Path(tempfile.mktemp(suffix=".log"))
    psi4.set_output_file(str(log_path), append=False)
    molden_path = None
    try:
        atoms = cjson["atoms"]
        coords = atoms["coords"]
        elem = atoms["elements"]["number"]
        charge_val = int(cjson.get("properties", {}).get("totalCharge", charge))
        spin_val = int(cjson.get("properties", {}).get("totalSpinMultiplicity", spin))
        ref = "uhf" if spin_val != 1 else "rhf"

        geom_lines = "\n".join(
            f"  {elem[i]:3d}  {coords[3*i]:12.8f}  {coords[3*i+1]:12.8f}  {coords[3*i+2]:12.8f}"
            for i in range(len(elem))
        )
        charge_tag = f"charge {charge_val}" if charge_val != 0 else ""
        mol_spec = f"{charge_tag}\n{geom_lines}\nno_com\nno_reorient\nsymmetry c1"
        mol = psi4.geometry(mol_spec)

        basis = _option(options, "basis", "def2-SVP")
        method = _option(options, "method", "hf")
        psi4.set_options({
            "basis": basis, "scf_type": "df", "reference": ref,
            "e_convergence": 1e-8, "d_convergence": 1e-8,
        })
        energy, wfn = psi4.energy(method, return_wfn=True)

        Ca = wfn.Ca()
        nocc = wfn.doccpi()[0] + wfn.soccpi()[0]
        mints = psi4.core.MintsHelper(wfn.basisset())
        bas_min = psi4.core.BasisSet.build(mol, "BASIS", "STO-3G", puream=0)
        S_full = mints.ao_overlap().np
        S_mixed = mints.ao_overlap(wfn.basisset(), bas_min).np
        S_min = mints.ao_overlap(bas_min, bas_min).np
        C_occ = Ca.np[:, :nocc]
        S_min_inv = np.linalg.inv(S_min)
        X = S_mixed @ S_min_inv @ S_mixed.T @ C_occ
        XtSX = X.T @ S_full @ X
        evals, evecs = np.linalg.eigh(XtSX)
        X_ortho = X @ (evecs @ np.diag(evals ** -0.5) @ evecs.T)
        new_Ca = Ca.np.copy()
        new_Ca[:, :nocc] = X_ortho
        wfn.Ca().copy(psi4.core.Matrix.from_array(new_Ca))

        molden_path = Path(tempfile.mktemp(suffix=".molden"))
        psi4.molden(wfn, str(molden_path))
        molden_text = molden_path.read_text(encoding="utf-8")
    finally:
        if molden_path:
            molden_path.unlink(missing_ok=True)
        try:
            log_path.unlink(missing_ok=True)
        except PermissionError:
            pass

    return {
        "moldenData": molden_text,
        "message": f"Computed IAO orbitals ({method}/{basis})",
    }
