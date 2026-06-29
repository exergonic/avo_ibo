import json
import subprocess
import tempfile
from pathlib import Path


PSI4_SCRIPT_TEMPLATE = '''import psi4
import json
import numpy as np

mol = psi4.geometry("""
{geometry}
no_com
no_reorient
symmetry c1
""")
psi4.set_options({{"basis": "{basis}", "scf_type": "df", "reference": "{reference}", "e_convergence": 1e-8, "d_convergence": 1e-8}})

energy, wfn = psi4.energy("{method}", return_wfn=True)
Ca = wfn.Ca()
nmo = Ca.shape[0]
nocc = wfn.doccpi()[0] + wfn.soccpi()[0]

mints = psi4.core.MintsHelper(wfn.basisset())

bas_min = psi4.core.BasisSet.build(mol, "BASIS", "STO-3G", puream=0)
bas_full = wfn.basisset()
S_full = mints.ao_overlap().np
S_mixed = mints.ao_overlap(bas_full, bas_min).np
S_min = mints.ao_overlap(bas_min, bas_min).np
C_occ = Ca.np[:, :nocc]
S_min_inv = np.linalg.inv(S_min)
X = S_mixed @ S_min_inv @ S_mixed.T @ C_occ
XtSX = X.T @ S_full @ X
evals, evecs = np.linalg.eigh(XtSX)
X_ortho = X @ (evecs @ np.diag(evals ** -0.5) @ evecs.T)
new_Ca = Ca.np.copy()
new_Ca[:, :nocc] = X_ortho
Ca = psi4.core.Matrix.from_array(new_Ca)
wfn.Ca().copy(Ca)

psi4.molden(wfn, "{molden_path}")
print("MOLDEN_OK")'''


def _make_p4_geometry(cjson):
    atoms = cjson["atoms"]
    coords = atoms["coords"]
    elem = atoms["elements"]["number"]
    nat = len(elem)
    lines = []
    for i in range(nat):
        x, y, z = coords[3 * i : 3 * i + 3]
        lines.append(f"  {elem[i]:3d}  {x:12.8f}  {y:12.8f}  {z:12.8f}")
    return "\n".join(lines) + "\n"


def _option(options, key, default):
    v = options.get(key, default)
    if isinstance(v, str):
        return v.strip()
    return v


def compute_ibo(cjson, options, charge, spin, debug=False):
    geometry = _make_p4_geometry(cjson)
    basis = _option(options, "basis", "def2-SVP")
    method = _option(options, "method", "hf")
    charge_val = int(cjson.get("properties", {}).get("totalCharge", charge))
    spin_val = int(cjson.get("properties", {}).get("totalSpinMultiplicity", spin))

    ref = "uhf" if spin_val != 1 else "rhf"
    charge_tag = f"charge {charge_val}" if charge_val != 0 else ""
    geom_with_charge = f"{charge_tag}\n{geometry}" if charge_tag else geometry

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".py", delete=False, encoding="utf-8"
    ) as f_script:
        script = PSI4_SCRIPT_TEMPLATE.format(
            geometry=geom_with_charge,
            basis=basis,
            method=method,
            reference=ref,
            molden_path="{molden_path}",
        )
        f_script.write(script)
        script_path = f_script.name

    molden_path_obj = Path(script_path).with_suffix(".molden")
    molden_path_str = molden_path_obj.as_posix()
    script_final = script.replace("{molden_path}", molden_path_str)

    with open(script_path, "w", encoding="utf-8") as f:
        f.write(script_final)

    try:
        proc = subprocess.run(
            ["psi4", script_path], capture_output=True, text=True, timeout=300
        )
        stdout = proc.stdout
        stderr = proc.stderr
        rc = proc.returncode
    except subprocess.TimeoutExpired:
        return {"error": "Psi4 calculation timed out after 300 seconds"}
    except FileNotFoundError:
        return {"error": "Psi4 executable not found"}
    finally:
        Path(script_path).unlink(missing_ok=True)

    if debug:
        import sys as _sys

        print(f"RC={rc}", file=_sys.stderr)
        print(f"STDOUT={stdout[:500]}", file=_sys.stderr)
        print(f"STDERR={stderr[:500]}", file=_sys.stderr)

    if rc != 0:
        return {"error": f"Psi4 failed (rc={rc}): {stdout[:500]}"}
    if "MOLDEN_OK" not in stdout:
        return {"error": f"Psi4 did not produce Molden output: {stdout[:500]}"}
    if not molden_path_obj.exists():
        return {"error": f"Molden file not found at {molden_path_str}"}

    molden_text = molden_path_obj.read_text(encoding="utf-8")
    molden_path_obj.unlink(missing_ok=True)

    return {
        "moldenData": molden_text,
        "message": f"Computed IAO orbitals ({method}/{basis})",
    }
