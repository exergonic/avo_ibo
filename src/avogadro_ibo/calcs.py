"""
IAO/IBO construction and Pipek-Mezey localization for Avogadro 2.

References:
  G. Knizia, JCTC 2013, 9, 4834-4843.  DOI: 10.1021/ct400687b
  ("Intrinsic Atomic Orbitals: An Unbiased Bridge between Quantum
   Theory and Chemical Concepts.")

Paper equation numbers and appendix references refer to the above.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Helper: extract per-function atom index and angular momentum from a BasisSet
# ---------------------------------------------------------------------------

def _get_basis_maps(basis):
    """
    Return arrays mapping each basis function in *basis* to its atom center
    (0-indexed) and angular momentum (0=s, 1=p, 2=d, ...).
    """
    atom_of = []
    am_of = []
    for sh in range(basis.nshell()):
        shell = basis.shell(sh)
        atom = shell.ncenter          # 0-indexed atom index
        am = shell.am                 # angular momentum quantum number
        nfunc = shell.nfunction       # total number of functions in shell
        for _ in range(nfunc):
            atom_of.append(atom)
            am_of.append(am)
    return np.array(atom_of, dtype=np.int32), np.array(am_of, dtype=np.int32)


# ---------------------------------------------------------------------------
# IAO construction   (Appendix C of Knizia JCTC 2013)
# ---------------------------------------------------------------------------

def _build_iao_basis(S, S12, S_min, C_occ):
    """
    Construct the Intrinsic Atomic Orbital (IAO) basis following the
    IAO/2014 algorithm (implemented in IboView's MakeIaoBasisNew).

    Parameters
    ----------
    S     : (n_AO, n_AO)  full AO overlap matrix
    S12   : (n_AO, n_min) overlap between full AO and minimal basis
    S_min : (n_min, n_min) minimal-basis overlap matrix
    C_occ : (n_AO, n_occ)  occupied MO coefficients

    Returns
    -------
    C_IAO     : (n_AO, n_min)  IAO coefficients, orthonormal w.r.t. S
    C_IAO_occ : (n_min, n_occ) occupied MO coefficients in the IAO basis
    """
    from scipy.linalg import cho_factor, cho_solve

    n_AO, n_occ = C_occ.shape
    n_min = S12.shape[1]

    # (1) Projector from minimal basis to AO basis: P12 = S^{-1} @ S12
    L_S, low_S = cho_factor(S)
    P12 = cho_solve((L_S, low_S), S12)            # (n_AO, n_min)

    # (2) Occupied MOs expressed in the minimal basis
    C_occ_min = S12.T @ C_occ                     # (n_min, n_occ)

    # (3) Solve S_min @ C_tilde = C_occ_min
    L_min, low_min = cho_factor(S_min)
    C_tilde = cho_solve((L_min, low_min), C_occ_min)   # (n_min, n_occ)

    # (4) Metric in the occupied space
    S_tilde = C_occ_min.T @ C_tilde                    # (n_occ, n_occ)

    # (5) Solve S_tilde @ C_tilde_2bar^T = C_tilde^T
    L_tilde, low_tilde = cho_factor(S_tilde)
    C_tilde_2bar_T = cho_solve((L_tilde, low_tilde), C_tilde.T)  # (n_occ, n_min)
    C_tilde_2bar = C_tilde_2bar_T.T                          # (n_min, n_occ)

    # (6) Residual part of occupied MOs beyond the minimal projection
    T4 = C_occ - P12 @ C_tilde_2bar                  # (n_AO, n_occ)

    # (7) Construct IAO coefficients
    C_IAO = P12 + T4 @ C_occ_min.T                   # (n_AO, n_min)

    # (8) Symmetric (Loewdin) orthogonalisation of IAOs
    #     Find M^{-1/2} where M = C_IAO^T @ S @ C_IAO
    metric = C_IAO.T @ S @ C_IAO                     # (n_min, n_min)
    evals, evecs = np.linalg.eigh(metric)
    C_IAO = C_IAO @ (evecs @ np.diag(evals ** -0.5) @ evecs.T)

    # Express the occupied MOs in the orthonormal IAO basis.
    # Since IAOs span the occupied space (by construction),
    # C_IAO @ C_IAO_occ = C_occ should hold exactly.
    C_IAO_occ = C_IAO.T @ S @ C_occ                # (n_min, n_occ)

    return C_IAO, C_IAO_occ


# ---------------------------------------------------------------------------
# Pipek-Mezey localisation in the IAO basis   (eq 4 and Appendix D)
# ---------------------------------------------------------------------------

def _localize_ibos(C_occ, atom_of, max_iter=2048, conv=1e-12):
    """
    Localise the occupied orbitals in the IAO basis by maximising

        L = sum_A sum_i  [n_A(i)]^4

    where n_A(i) = sum_{mu in A} C(mu,i)^2 is the electron population of
    orbital i on atom A (in the orthonormal IAO basis).

    This is the IBO functional of Knizia JCTC 2013 with exponent p=4.
    The procedure follows the standard Pipek-Mezey Jacobi sweep
    (Appendix D of the paper), but in the IAO basis.

    To avoid local minima in the p=4 functional when the IAO basis
    differs from the SCF basis (e.g. STO-3G IAOs but def2-SVP SCF),
    the p=2 functional (convex, no local minima) is run first as a
    warm-start, then p=4 refines.

    Parameters
    ----------
    C_occ   : (n_IAO, n_occ)  occupied coefficients in IAO basis
              (modified in place)
    atom_of : (n_IAO,)         atom index for each IAO basis function
    max_iter: int              maximum sweeps per functional
    conv    : float            gradient-norm convergence threshold

    Returns
    -------
    n_sweeps : total sweeps performed (p=2 warmup + p=4 refine)
    """
    n_IAO, n_occ = C_occ.shape
    n_atoms = int(np.max(atom_of)) + 1

    # --- Random perturbation to break symmetry (IboView uses 18 deg) ---
    rng = np.random.default_rng()
    U, _ = np.linalg.qr(rng.normal(0, 1, (n_occ, n_occ)))
    C_occ[:] = C_occ @ U

    total_sweeps = 0

    for exponent, inv_p in ((2, 0.5), (4, 0.25)):
        # Run p=exponent localisation (p=2 warmup, then p=4 refine)
        for _ in range(max_iter):
            grad_norm = 0.0

            for i in range(1, n_occ):
                for j in range(i):
                    ci = C_occ[:, i]
                    cj = C_occ[:, j]

                    Qii = np.zeros(n_atoms, dtype=np.float64)
                    Qjj = np.zeros(n_atoms, dtype=np.float64)
                    Qij = np.zeros(n_atoms, dtype=np.float64)
                    np.add.at(Qii, atom_of, ci * ci)
                    np.add.at(Qjj, atom_of, cj * cj)
                    np.add.at(Qij, atom_of, ci * cj)

                    if exponent == 2:
                        # Pipek-Mezey p=2 (convex)
                        Aij = 0.0
                        Bij = 0.0
                        for A in range(n_atoms):
                            qii, qjj, qij = Qii[A], Qjj[A], Qij[A]
                            Aij += -qii * qii - qjj * qjj + 2.0 * qij * qij
                            Bij += 4.0 * qij * (qii - qjj)
                        if abs(Aij) <= conv:
                            continue
                        phi = 0.5 * np.arctan2(Bij, -Aij)
                        grad_term = 2.0
                    else:
                        # Pipek-Mezey p=4 (eq 4 / Appendix D)
                        Aij = 0.0
                        Bij = 0.0
                        for A in range(n_atoms):
                            qii, qjj, qij = Qii[A], Qjj[A], Qij[A]
                            qii_2 = qii * qii; qjj_2 = qjj * qjj; qij_2 = qij * qij
                            Aij += (-qii_2 * qii_2 - qjj_2 * qjj_2
                                    + 6.0 * (qii_2 + qjj_2) * qij_2
                                    + qii_2 * qii * qjj + qii * qjj_2 * qjj)
                            Bij += 4.0 * qij * (qii_2 * qii - qjj_2 * qjj)
                        if abs(Aij) <= conv:
                            continue
                        phi = 0.25 * np.arctan2(Bij, -Aij)
                        grad_term = 4.0

                    cs = np.cos(phi)
                    sn = np.sin(phi)

                    for mu in range(n_IAO):
                        old_i = C_occ[mu, i]
                        old_j = C_occ[mu, j]
                        C_occ[mu, i] = cs * old_i + sn * old_j
                        C_occ[mu, j] = cs * old_j - sn * old_i

                    grad_norm += (grad_term * phi * Bij) ** 2

            grad_norm = np.sqrt(grad_norm) / n_occ
            total_sweeps += 1

            if grad_norm < conv:
                break

    return total_sweeps


# ---------------------------------------------------------------------------
# Resolve on-atom degeneracies that PM cannot separate
# ---------------------------------------------------------------------------

def _resolve_on_atom_mixing(C_occ, atom_of, F_IAO, dom_threshold=0.99):
    """
    Diagonalise F_IAO within each group of occupied orbitals that share
    the same dominant atom and have DOM > *dom_threshold*.

    The PM functional uses only atomic populations n_A(i), so two
    orbitals on the same atom (e.g. O 2s and O lone pair) are
    degenerate in the functional — any rotation within the subspace
    gives the same L value.  This routine breaks that degeneracy by
    the aufbau principle: the eigenvectors of F_IAO within the
    subspace give the lowest-energy (most s-like) to highest-energy
    (most p-like) orbitals.

    Parameters are modified in-place.
    """
    n_IAO, n_occ = C_occ.shape
    n_atoms = int(np.max(atom_of)) + 1

    sq = C_occ ** 2
    pop = np.zeros((n_occ, n_atoms), dtype=np.float64)
    for i in range(n_occ):
        np.add.at(pop[i], atom_of, sq[:, i])

    # Identify same-atom, high-DOM groups
    groups = {}
    for i in range(n_occ):
        order = np.argsort(-pop[i])
        top_A = order[0]
        dom_val = pop[i, top_A] ** 2 + pop[i, order[1]] ** 2
        if dom_val > dom_threshold:
            groups.setdefault(top_A, []).append(i)

    for atom, indices in groups.items():
        n_g = len(indices)
        if n_g < 2:
            continue
        C_block = C_occ[:, indices]                     # (n_IAO, n_g)
        Fb = C_block.T @ (F_IAO @ C_block)              # (n_g, n_g)
        evals, evecs = np.linalg.eigh(Fb)
        C_occ[:, indices] = C_block @ evecs


# ---------------------------------------------------------------------------
# IBO analysis table
# ---------------------------------------------------------------------------

# Periodic-table lookup for element symbols
_ELEM_SYMBOLS = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
]


def _analyze_ibos(C_IAO_all, occ_all, energies_all, nocc,
                  atom_of, am_of, elem, method, basis, ref):
    """
    Build a formatted IBO analysis table covering all IAO-basis orbitals.

    For each orbital (occupied IBO or valence-virtual IAO), compute:
      - per-atom populations from IAO coefficients
      - DOM (largest two n_A fractions summed)
      - s/p/d character on the dominant atom
    """
    n_IAO, n_orb = C_IAO_all.shape
    n_atoms = len(elem)

    lines = []
    orbid_labels = [""] * n_orb
    lines.append(f"IBO Analysis  ({method}/{basis}, {ref.upper()})")
    lines.append("")
    header = (
        f"  {'#':>3}  {'Occ':>7}  {'Energy':>10}  {'DOM':>5}  "
        f"{'Type':>12}  {('Composition'):<38}  Hybrid"
    )
    lines.append(header)
    lines.append("-" * len(header))

    for orb in range(n_orb):
        oc = occ_all[orb]
        sq = C_IAO_all[:, orb] ** 2
        pop = np.zeros(n_atoms, dtype=np.float64)
        np.add.at(pop, atom_of, sq)

        # Dominant atom and its population
        order = np.argsort(-pop)
        top_A = order[0]
        top_B = order[1]
        dom = pop[top_A] ** 2 + pop[top_B] ** 2

        # s/p/d breakdown on the dominant atom
        s_char = _s_char(C_IAO_all[:, orb], am_of, top_A, atom_of)
        p_char = _p_char(C_IAO_all[:, orb], am_of, top_A, atom_of)
        d_char = _d_char(C_IAO_all[:, orb], am_of, top_A, atom_of)

        # Determine orbital type
        if oc > 1.5:
            # Occupied — classify like before
            if pop[top_A] > 0.99 and s_char > 0.75:
                orbid = f"{_elem_symbol(elem[top_A])}(Core)"
            elif pop[top_A] > 0.90:
                orbid = f"{_elem_symbol(elem[top_A])}(LP)"
            elif pop[top_A] + pop[top_B] > 0.75:
                def _p_frac(c, am_of, atom, atom_of):
                    _s = _s_char(c, am_of, atom, atom_of)
                    _p = _p_char(c, am_of, atom, atom_of)
                    _d = _d_char(c, am_of, atom, atom_of)
                    _t = _s + _p + _d
                    return _p / _t if _t > 0 else 0.0
                pfrac_A = _p_frac(C_IAO_all[:, orb], am_of, top_A, atom_of)
                pfrac_B = _p_frac(C_IAO_all[:, orb], am_of, top_B, atom_of)
                bond_type = "pi" if (pfrac_A > 0.85 and pfrac_B > 0.85) else "sigma"
                symA = _elem_symbol(elem[top_A])
                symB = _elem_symbol(elem[top_B])
                if top_A < top_B:
                    orbid = f"{symA}-{symB} {bond_type}"
                else:
                    orbid = f"{symB}-{symA} {bond_type}"
            elif pop[top_A] > 0.70:
                symA = _elem_symbol(elem[top_A])
                if s_char > 0.5:
                    orbid = f"{symA}(LP-s)"
                else:
                    orbid = f"{symA}(LP)"
            else:
                orbid = "Deloc"
        else:
            # Virtual — label as anti-bonding where appropriate
            if pop[top_A] + pop[top_B] > 0.60:
                symA = _elem_symbol(elem[top_A])
                symB = _elem_symbol(elem[top_B])
                if top_A < top_B:
                    orbid = f"{symA}-{symB} anti*"
                else:
                    orbid = f"{symB}-{symA} anti*"
            elif pop[top_A] > 0.50:
                orbid = f"{_elem_symbol(elem[top_A])}(virt)"
            else:
                orbid = "Virt"

        orbid_labels[orb] = orbid

        comp_parts = []
        for A in order[:4]:
            if pop[A] > 0.02:
                sym = _elem_symbol(elem[A])
                pct = pop[A] * 100.0
                comp_parts.append(f"{sym}({pct:.1f}%)")
        comp = " + ".join(comp_parts)

        s_pct = s_char * 100.0
        p_pct = p_char * 100.0
        d_pct = d_char * 100.0
        hybrid = ""
        if s_pct > 1 or p_pct > 1 or d_pct > 1:
            hybrid = f"{s_pct:.0f}% s + {p_pct:.0f}% p + {d_pct:.0f}% d"

        lines.append(
            f"  {orb+1:>3d}  {oc:>7.3f}  {energies_all[orb]:>10.6f}  "
            f"{dom:>5.3f}  {orbid:>12}  {comp:<38}  {hybrid}"
        )

    lines.append("")
    lines.append(f"Total electrons: {int(2 * nocc) if ref == 'rhf' else nocc}")
    return "\n".join(lines), orbid_labels


def _elem_symbol(Z):
    Z = int(round(Z))
    if Z < len(_ELEM_SYMBOLS):
        return _ELEM_SYMBOLS[Z]
    return f"E{Z}"


def _s_char(c, am_of, atom, atom_of):
    """Fraction of s-orbital (am=0) contribution on the given atom."""
    idx = np.where((atom_of == atom) & (am_of == 0))[0]
    return float(np.sum(c[idx] ** 2)) if len(idx) else 0.0


def _p_char(c, am_of, atom, atom_of):
    """Fraction of p-orbital (am=1) contribution on the given atom."""
    idx = np.where((atom_of == atom) & (am_of == 1))[0]
    return float(np.sum(c[idx] ** 2)) if len(idx) else 0.0


def _d_char(c, am_of, atom, atom_of):
    """Fraction of d-orbital (am=2) contribution on the given atom."""
    idx = np.where((atom_of == atom) & (am_of == 2))[0]
    return float(np.sum(c[idx] ** 2)) if len(idx) else 0.0


# ---------------------------------------------------------------------------
# Canonical MO delocalization analysis (atom-density matching in AO basis)
# ---------------------------------------------------------------------------

# -- (canonical MO delocalization analysis removed per user request in 2026-06-30;
# -- the canonical.molden file is written below for reference in Avogadro)--


# ---------------------------------------------------------------------------
# Molden writer using IAO-basis orbitals
# ---------------------------------------------------------------------------

def _write_iao_molden(path, wfn, C_AO, occ, energies, n_orb):
    """
    Write a Molden file whose [MO] section contains IAO-basis orbitals.

    The [Atoms] and [GTO] header sections are copied from Psi4's own Molden
    output; only the [MO] block is replaced with the IAO-basis orbitals.
    """
    import tempfile
    tmp = path.with_suffix(".molden.tmp")
    psi4_molden = __import__("psi4").molden
    psi4_molden(wfn, str(tmp))
    text = tmp.read_text(encoding="utf-8")
    tmp.unlink()

    # Keep everything before the [MO] section
    mo_tag = "[MO]"
    idx = text.find(mo_tag)
    if idx == -1:
        raise RuntimeError("Psi4 Molden output has no [MO] section")
    header = text[:idx]

    n_AO = C_AO.shape[0]
    lines = [header + "\n[MO]\n"]

    for i in range(n_orb):
        ei = energies[i]
        oi = occ[i]
        lines.append(f" Sym= A\n Ene= {ei:15.10f}\n Spin= Alpha\n"
                     f" Occup= {oi:14.10f}\n")
        for j in range(n_AO):
            lines.append(f"  {j + 1:>4d}  {C_AO[j, i]:16.10f}\n")

    path.write_text("".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def _option(options, key, default):
    v = options.get(key, default)
    if isinstance(v, str):
        return v.strip()
    return v


def compute_ibo(cjson, options, charge, spin, debug=False):
    import logging
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

    import psi4
    import numpy as np

    psi4.set_output_file(str(TEMP_DIR / "psi4.log"), append=True)

    # -- Parse input geometry and options -----------------------------------
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
    # NOTE: puream=0 gives Cartesian basis functions, which is what the
    # paper assumes.  Changing this would affect the IAO construction.
    energy, wfn = psi4.energy(method, return_wfn=True)

    # -- Extract occupied coefficients and overlap matrices ----------------
    Ca = wfn.Ca()
    nocc = wfn.doccpi()[0] + wfn.soccpi()[0]
    mints = psi4.core.MintsHelper(wfn.basisset())

    # Full-AO overlap, full-minimal cross overlap, minimal overlap
    S_full = mints.ao_overlap().np
    bas_min = psi4.core.BasisSet.build(mol, "BASIS", "STO-3G", puream=0)
    S_min = mints.ao_overlap(bas_min, bas_min).np
    S12 = mints.ao_overlap(wfn.basisset(), bas_min).np

    C_occ = Ca.np[:, :nocc].copy()           # (n_AO, n_occ)

    # -- Build IAO basis (Appendix C) --------------------------------------
    C_IAO, C_IAO_occ = _build_iao_basis(S_full, S12, S_min, C_occ)

    # -- Atom / angular-momentum map for the minimal basis -----------------
    # Each IAO inherits the atom and AM from the minimal-basis function
    # it was built from.
    atom_of, am_of = _get_basis_maps(bas_min)  # both (n_min,)

    # -- Pipek-Mezey localisation in IAO basis (eq 4 / Appendix D) --------
    _localize_ibos(C_IAO_occ, atom_of,
                   max_iter=2048, conv=1e-12)

    # -- Compute orbital energies from Fock matrix -------------------------
    F_AO = wfn.Fa().np                              # (n_AO, n_AO)
    F_IAO = C_IAO.T @ F_AO @ C_IAO                 # (n_min, n_min)

    # -- Resolve on-atom degeneracies that PM cannot separate --------------
    # The PM functional uses atomic populations n_A(i), so orbitals on the
    # same atom with DOM ≈ 1 are degenerate (O 2s and O lone pair mix
    # arbitrarily).  Diagonalise F_IAO within each such subspace to restore
    # energy ordering (s-rich lowest, p-rich highest).
    _resolve_on_atom_mixing(C_IAO_occ, atom_of, F_IAO)

    occ_energies = np.array(
        [C_IAO_occ[:, i].dot(F_IAO @ C_IAO_occ[:, i]) for i in range(nocc)]
    )

    # -- Valence-virtual IAOs via SVD  (IboView MakeValenceVirtuals) -------
    C_all = Ca.np
    C_vir = C_all[:, nocc:]                         # (n_AO, n_vir)
    SIbVir = C_IAO.T @ S_full @ C_vir               # (n_min, n_vir)
    U_svd, Sigma, _ = np.linalg.svd(SIbVir, full_matrices=False)
    n_val_vir = int(np.sum(Sigma > 1e-8))
    U_val = U_svd[:, :n_val_vir]                    # (n_min, n_val_vir)

    # -- Localize the virtual block too (IboView localizes ALL case blocks) ---
    if n_val_vir > 1:
        _localize_ibos(U_val, atom_of,
                       max_iter=2048, conv=1e-12)

    vir_energies = np.array(
        [U_val[:, i].dot(F_IAO @ U_val[:, i]) for i in range(n_val_vir)]
    )

    # -- Combined IAO-basis orbital set, sorted by energy ------------------
    C_IAO_all = np.hstack([C_IAO_occ, U_val])       # (n_min, n_orb)
    occ_all = np.array([2.0] * nocc + [0.0] * n_val_vir)
    energies_all = np.concatenate([occ_energies, vir_energies])

    order = np.argsort(energies_all)
    C_IAO_all = C_IAO_all[:, order]
    occ_all = occ_all[order]
    energies_all = energies_all[order]

    # -- Write Molden with IAO-basis orbitals ------------------------------
    C_AO_all = C_IAO @ C_IAO_all                    # (n_AO, n_orb)
    molden_path = TEMP_DIR / "ibo.molden"
    _write_iao_molden(molden_path, wfn, C_AO_all, occ_all, energies_all,
                      C_AO_all.shape[1])
    molden_text = molden_path.read_text(encoding="utf-8")

    # -- Canonical Molden (for reference in Avogadro's MO surface dialog) ---
    canon_path = TEMP_DIR / "canonical.molden"
    _psi_molden = __import__("psi4").molden
    _psi_molden(wfn, str(canon_path))

    # -- IBO analysis table -------------------------------------------------
    msg, labels = _analyze_ibos(C_IAO_all, occ_all, energies_all, nocc,
                                atom_of, am_of, elem, method, basis, ref)
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
