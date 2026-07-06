"""
IAO/IBO construction and Pipek-Mezey localization for Avogadro 2.

BSD 3-Clause License
Copyright (c) 2025-2026, Billy Wayne McCann
SPDX-License-Identifier: BSD-3-Clause

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
    (0-indexed), angular momentum (0=s, 1=p, 2=d, ...), principal quantum
    number n, and orbital subtype label.

    The principal quantum number n is inferred from shell ordering per atom
    (1s→2s→2p→3s→3p→3d→4s→4p), which matches the STO-3G shell layout.

    The subtype label identifies:
      - p-functions: "px", "py", "pz"  (Psi4 Cartesian order: x, y, z)
      - d-functions: "dxx","dxy","dxz","dyy","dyz","dzz" (Psi4 order)
      - s-functions: "" (empty)
    """
    atom_of = []
    am_of = []
    n_of = []
    dtype_of = []

    # Track the next n to assign for each (atom, am) pair
    #   s(am=0): start at n=1, increment by 1 per s shell
    #   p(am=1): start at n=2, increment by 1 per p shell
    #   d(am=2): start at n=3, increment by 1 per d shell
    # This matches STO-3G's aufbau ordering of shells.
    _P_AM_START = {0: 1, 1: 2, 2: 3}
    _P_SUBTYPE = {
        1: ['px', 'py', 'pz'],              # Psi4 Cartesian p order: x, y, z
        2: ['dxx', 'dxy', 'dxz', 'dyy', 'dyz', 'dzz'],  # Psi4 Cartesian d order
    }

    next_n_per_atom = {}

    for sh in range(basis.nshell()):
        shell = basis.shell(sh)
        atom = shell.ncenter
        am = shell.am
        nfunc = shell.nfunction

        # Determine principal quantum number for this shell
        key = (atom, am)
        next_n = next_n_per_atom.get(key, _P_AM_START.get(am, 1))
        shell_n = next_n
        next_n_per_atom[key] = next_n + 1

        subtypes = _P_SUBTYPE.get(am, [''] * nfunc)
        for f_idx in range(nfunc):
            atom_of.append(atom)
            am_of.append(am)
            n_of.append(shell_n)
            dtype_of.append(subtypes[f_idx] if f_idx < len(subtypes) else '')

    return (np.array(atom_of, dtype=np.int32),
            np.array(am_of, dtype=np.int32),
            np.array(n_of, dtype=np.int32),
            dtype_of)


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

def _localize_ibos(C_occ, atom_of, max_iter=2048, conv=1e-12,
                   exponents=(2, 4), cayley_deg=0.0, seed=42):
    """
    Localise the occupied orbitals in the IAO basis by maximising

        L = Σ_A Σ_i  [n_A(i)]^p

    where n_A(i) = Σ_{μ ∈ A} C(μ,i)² is the electron population of
    orbital i on atom A (in the orthonormal IAO basis) and p is the PM
    exponent.

    The procedure follows the standard Pipek-Mezey Jacobi sweep
    (Appendix D of the paper), but in the IAO basis.

    Parameters
    ----------
    C_occ   : (n_IAO, n_occ)  coefficients in IAO basis (modified in place)
    atom_of : (n_IAO,)         atom index for each IAO basis function
    max_iter: int              maximum sweeps per functional
    conv    : float            gradient-norm convergence threshold
    exponents: tuple of PM exponents to apply sequentially.
              (2,)  matches IboView GUI default (p=2 only).
              (2, 4)  adds p=4 refinement (default; sharper convergence
                      for bond-direction p-vector alignment).
    cayley_deg: float          Cayley random rotation angle in degrees
              (IboView: 18°).  Set to 0 (default) to skip — the fixed
              sequential sweep order starting from canonical MOs gives
              the best energy degeneracy for symmetric molecules.
    seed    : int              RNG seed (only used if cayley_deg > 0).

    Returns
    -------
    n_sweeps : total sweeps performed
    """
    n_IAO, n_occ = C_occ.shape
    n_atoms = int(np.max(atom_of)) + 1

    # Cayley random rotation (IboView's RotateVectorsRandomly)
    if cayley_deg > 0:
        rng = np.random.default_rng(seed)
        sigma = cayley_deg * np.pi / 180.0
        A = rng.normal(0, sigma, (n_occ, n_occ))
        A = (A - A.T) / 2  # anti-symmetric
        U = np.linalg.solve(np.eye(n_occ) - 0.5 * A,
                            np.eye(n_occ) + 0.5 * A)
        C_occ[:] = C_occ @ U

    total_sweeps = 0

    for exponent in exponents:
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
                    elif exponent == 4:
                        # Pipek-Mezey p=4 (eq 4).  The published Appendix D
                        # 2x2 update formulas contain a production error
                        # (confirmed by Knizia at https://sites.psu.edu/knizia/software/).
                        # These formulas match the corrected reference
                        # implementation (ibo-ref).
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
                    else:
                        raise ValueError(f"Unsupported PM exponent: {exponent}")

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
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I",
]


def _d_spherical_weights(c, atom_idx, atom_of, am_of):
    """
    Compute weights for each spherical d-type from Cartesian d coefficients
    on a given atom.  The 6 Cartesian d-functions in Psi4 (puream=0) are
    ordered: xx, xy, xz, yy, yz, zz.  We project the coefficient vector
    onto the five spherical harmonic directions.

    Returns dict of {name: weight} where weight = squared projection.
    """
    idx = np.where((atom_of == atom_idx) & (am_of == 2))[0]
    if len(idx) < 6:
        return {}
    c = np.asarray(c, dtype=np.float64)
    # Psi4 Cartesian d order: xx, xy, xz, yy, yz, zz
    c_xx, c_xy, c_xz, c_yy, c_yz, c_zz = c[idx[:6]]
    return {
        'dxy':    c_xy ** 2,
        'dxz':    c_xz ** 2,
        'dyz':    c_yz ** 2,
        'dz2':    (-c_xx - c_yy + 2 * c_zz) ** 2,
        'dx2y2':  (c_xx - c_yy) ** 2,
    }


def _hybrid_str(c, am_of, atom_of, func_n, func_dtype, top_atom):
    """
    Build a specific hybrid label for the dominant atom, e.g.
        "57% 4s + 43% 3dz²"
        "100% 1s"
        "100% 4pz"
        "83% 3s + 17% 3pz"
        "46% 4s + 54% 3d"
    """
    c = np.asarray(c, dtype=np.float64)
    pA = float(np.sum(c[np.where(atom_of == top_atom)] ** 2))
    if pA < 1e-12:
        return ""

    parts = []
    for am_label, am_val in [('s', 0), ('p', 1), ('d', 2)]:
        idx_am = np.where((atom_of == top_atom) & (am_of == am_val))[0]
        if len(idx_am) == 0:
            continue
        total_am = float(np.sum(c[idx_am] ** 2))
        pct = total_am / pA * 100.0
        if pct < 1.0:
            continue

        # Find dominant n within this l-subspace
        n_counts = {}
        for fi in idx_am:
            n_key = func_n[fi]
            n_counts[n_key] = n_counts.get(n_key, 0.0) + c[fi] ** 2
        dominant_n = max(n_counts, key=n_counts.get)

        # Determine dominant subtype
        subtype = ""
        if am_val == 1:   # p-orbitals: px, py, pz
            st_counts = {}
            for fi in idx_am:
                st = func_dtype[fi]
                if st:
                    st_counts[st] = st_counts.get(st, 0.0) + c[fi] ** 2
            if st_counts:
                top_st = max(st_counts, key=st_counts.get)
                if st_counts[top_st] > 0.5 * total_am:
                    subtype = top_st  # e.g. "pz"
        elif am_val == 2:   # d-orbitals: dxy, dxz, dyz, dz2, dx2y2
            d_weights = _d_spherical_weights(c, top_atom, atom_of, am_of)
            if d_weights:
                top_st = max(d_weights, key=d_weights.get)
                if d_weights[top_st] > 0.5 * max(d_weights.values()):
                    subtype = top_st  # e.g. "dz2"

        label = str(dominant_n) + (subtype if subtype else am_label)
        parts.append(f"{pct:.0f}% {label}")

    return " + ".join(parts)


def _wiberg_per_ibo(pop, occ, A, B):
    """
    Per-IBO contribution to the Wiberg bond order between atoms A and B.

    In the orthonormal IAO basis, the density contribution from a single IBO
    with coefficient vector c_k is D^{(k)} = occ_k · c_k c_k^T.  The Wiberg
    index between A and B from this IBO is:

        W_AB^{(k)} = Σ_{i∈A} Σ_{j∈B} (D^{(k)}_ij)²

    which simplifies (by the independence of i and j sums) to:

        W_AB^{(k)} = occ_k² · P_A · P_B

    where P_X = Σ_{i∈X} c_{k,i}² is the Mulliken population on atom X.

    For RHF occupied (occ=2): W_AB = 4 · P_A · P_B, ranging from 0 (pure
    ionic, no shared density) to 1 (pure covalent 2c-2e bond with 50/50
    sharing).
    """
    return float(occ ** 2 * pop[A] * pop[B])


def _ionic_pct(pop, A, B):
    """
    Percent ionic character between atoms A and B from per-atom populations.

        Ionic% = |P_A - P_B| / (P_A + P_B) × 100

    Ranges from 0% (pure covalent, equal sharing) to 100% (pure ionic,
    all density on one atom).
    """
    num = abs(pop[A] - pop[B])
    den = pop[A] + pop[B]
    return num / den * 100.0 if den > 1e-12 else 0.0


def _analyze_ibos(C_IAO_all, occ_all, energies_all, nocc,
                  atom_of, am_of, func_n, func_dtype,
                  elem, method, basis, ref):
    """
    Build a formatted IBO analysis table covering all IAO-basis orbitals.

    For each orbital (occupied IBO or valence-virtual IAO), compute:
      - per-atom populations from IAO coefficients
      - DOM (largest two n_A fractions summed)
      - per-IBO Wiberg bond order (W_AB) and percent ionic character
      - specific nl/subtype hybrid label on the dominant atom
    """
    n_IAO, n_orb = C_IAO_all.shape
    n_atoms = len(elem)

    lines = []
    orbid_labels = [""] * n_orb
    atom_pop = np.zeros(n_atoms, dtype=np.float64)  # accumulated per-atom electron counts
    lines.append(f"IBO Analysis  ({method}/{basis}, {ref.upper()})")
    lines.append("")
    header = (
        f"  {'#':>3}  {'Occ':>7}  {'Energy':>10}  {'DOM':>5}  "
        f"{'Type':>12}  {('Composition'):<38}  {'Hybrid':<25}  "
        f"{'W_AB':>7}  {'Ion%':>5}"
    )
    lines.append(header)
    lines.append("-" * len(header))

    for orb in range(n_orb):
        oc = occ_all[orb]
        sq = C_IAO_all[:, orb] ** 2
        pop = np.zeros(n_atoms, dtype=np.float64)
        np.add.at(pop, atom_of, sq)
        if oc > 1.5:
            atom_pop += pop * oc       # accumulate electron count per atom

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
            elif pop[top_A] + pop[top_B] > 0.75 and pop[top_B] > 0.02:
                def _p_frac(c, am_of, atom, atom_of):
                    _s = _s_char(c, am_of, atom, atom_of)
                    _p = _p_char(c, am_of, atom, atom_of)
                    _d = _d_char(c, am_of, atom, atom_of)
                    _t = _s + _p + _d
                    return _p / _t if _t > 0 else 0.0
                pfrac_A = _p_frac(C_IAO_all[:, orb], am_of, top_A, atom_of)
                pfrac_B = _p_frac(C_IAO_all[:, orb], am_of, top_B, atom_of)
                bond_type = "pi" if (pfrac_A > 0.85 and pfrac_B > 0.85) else "sigma"
                a, b = sorted([top_A, top_B])
                symA = _elem_symbol(elem[a])
                symB = _elem_symbol(elem[b])
                orbid = f"{symA}-{symB} {bond_type}"
            elif pop[top_A] > 0.70:
                symA = _elem_symbol(elem[top_A])
                if s_char > 0.5:
                    orbid = f"{symA}(LP-s)"
                else:
                    orbid = f"{symA}(LP)"
            else:
                # Check for 2-electron 3-center bond
                if len(order) >= 3 and pop[order[2]] > 0.10:
                    atoms = sorted([order[0], order[1], order[2]],
                                  key=lambda i: (_elem_symbol(elem[i]), i))
                    syms = '-'.join(_elem_symbol(elem[a]) for a in atoms)
                    orbid = f"{syms} 2e3c"
                else:
                    orbid = "Deloc"
        else:
            # Virtual — label as anti-bonding counterpart
            if len(order) >= 3 and pop[order[2]] > 0.08:
                # 3-center antibonding counterpart
                atoms = sorted([order[0], order[1], order[2]],
                              key=lambda i: (_elem_symbol(elem[i]), i))
                syms = '-'.join(_elem_symbol(elem[a]) for a in atoms)
                orbid = f"{syms} 2e3c anti*"
            elif pop[top_A] + pop[top_B] > 0.60 and pop[top_B] > 0.02:
                a, b = sorted([top_A, top_B])
                symA = _elem_symbol(elem[a])
                symB = _elem_symbol(elem[b])
                orbid = f"{symA}-{symB} anti*"
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

        # For hybrid label, prefer a non-H atom when it carries meaningful density
        hybrid_atom = top_A
        if elem[top_A] == 1:
            for A in order[1:]:
                if elem[A] != 1 and pop[A] > 0.02:
                    hybrid_atom = A
                    break
        hybrid = _hybrid_str(C_IAO_all[:, orb], am_of, atom_of,
                              func_n, func_dtype, hybrid_atom)
        if hybrid_atom != top_A:
            hybrid = f"{_elem_symbol(elem[hybrid_atom])}: {hybrid}"

        # Per-IBO Wiberg bond order and ionic character (between top_A, top_B)
        w_ab = _wiberg_per_ibo(pop, oc, top_A, top_B)
        if oc > 1.5 and w_ab > 0.001:
            ion_str = f"{_ionic_pct(pop, top_A, top_B):.1f}"
        else:
            ion_str = "---"

        lines.append(
            f"  {orb+1:>3d}  {oc:>7.3f}  {energies_all[orb]:>10.6f}  "
            f"{dom:>5.3f}  {orbid:>12}  {comp:<38}  {hybrid:<25}  "
            f"{w_ab:>7.3f}  {ion_str:>5}"
        )

    lines.append("")
    lines.append(f"Total electrons: {int(2 * nocc) if ref == 'rhf' else nocc}")

    charge_section = _format_charge_decomposition(atom_pop, elem)
    lines.append(charge_section)

    net_charges = [
        float(int(round(elem[A])) - atom_pop[A]) for A in range(n_atoms)
    ]
    return "\n".join(lines), orbid_labels, net_charges


def _format_total_wiberg(C_IAO_occ, atom_of, elem):
    """
    Compute total Wiberg bond order matrix in the orthonormal IAO basis and
    format a compact table of atom pairs with significant bond order.

    The density matrix in the IAO basis (RHF) is:

        D = 2 · C_IAO_occ @ C_IAO_occ.T

    where C_IAO_occ has shape (n_min, n_occ).  The Wiberg index between
    atoms A and B is:

        W_AB = Σ_{i∈A} Σ_{j∈B} D_ij²

    Only pairs with W_AB > 0.01 are shown.
    """
    # RHF: D = 2 · C_occ @ C_occ.T  (occupied density in orthonormal IAO basis)
    D = 2.0 * C_IAO_occ @ C_IAO_occ.T       # (n_min, n_min)
    D_sq = D ** 2

    n_atoms = len(elem)
    W = np.zeros((n_atoms, n_atoms), dtype=np.float64)
    for i in range(D_sq.shape[0]):
        ai = atom_of[i]
        for j in range(D_sq.shape[1]):
            aj = atom_of[j]
            W[ai, aj] += D_sq[i, j]

    pairs = []
    for A in range(n_atoms):
        for B in range(A + 1, n_atoms):
            w = W[A, B]
            if w > 0.01:
                symA = _elem_symbol(elem[A])
                symB = _elem_symbol(elem[B])
                pairs.append((symA, symB, w))

    if not pairs:
        return ""

    lines = ["", "", "--- Total Wiberg Bond Orders ---"]
    for symA, symB, w in sorted(pairs, key=lambda x: -x[2]):
        lines.append(f"  {symA}-{symB}    {w:>7.3f}")
    return "\n".join(lines)


def _format_charge_decomposition(atom_pop, elem):
    """
    Format a charge decomposition table from accumulated IAO populations.

    For each atom A:
        Q_A = Σ_k occ_k · P_A^{(k)}   (total electrons on atom A)
        Net charge = Z_A - Q_A

    ``atom_pop[A]`` is Q_A as a float.  ``elem`` gives atomic numbers.
    """
    lines = ["", "--- Charge Decomposition ---"]
    header = f"  {'Atom':>4}  {'Z':>3}  {'Pop':>8}  {'Net Charge':>10}"
    lines.append(header)
    lines.append("-" * len(header))
    total_pop = 0.0
    total_z = 0
    for A in range(len(elem)):
        Z = int(round(elem[A]))
        pop = atom_pop[A]
        net = Z - pop
        sym = _elem_symbol(Z)
        lines.append(f"  {sym:>4}  {Z:>3d}  {pop:>8.3f}  {net:>+10.3f}")
        total_pop += pop
        total_z += Z
    lines.append("-" * len(header))
    lines.append(
        f"Total:  {total_z:>3d}  {total_pop:>8.3f}  "
        f"{total_z - total_pop:>+10.3f}"
    )
    return "\n".join(lines)


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


# -- (canonical MO deloc analysis removed 2026-06-30; canonical.molden below)--


# ---------------------------------------------------------------------------
# Molden writer using IAO-basis orbitals
# ---------------------------------------------------------------------------

def _write_iao_molden(path, wfn, C_AO, occ, energies, n_orb):
    """
    Write a Molden file whose [MO] section contains IAO-basis orbitals.

    The [Atoms] and [GTO] header sections are copied from Psi4's own Molden
    output; only the [MO] block is replaced with the IAO-basis orbitals.

    The [MO] section is padded with zero-energy dummy orbitals up to n_AO
    total entries so Avogadro's MO slot count matches the [GTO] basis set
    size, preventing uninitialised-slot noise.
    """
    import tempfile
    import numpy as np
    tmp = path.with_suffix(".molden.tmp")
    psi4_molden = __import__("psi4").molden
    psi4_molden(wfn, str(tmp))
    text = tmp.read_text(encoding="utf-8")
    tmp.unlink()

    # Build index permutation and re-scaling vector to convert Psi4 internal AO
    # order to Molden standard.  Psi4's CCA convention uses unnormalized Cartesian
    # Gaussians (off-diagonal d/f have self-overlap < 1); the Molden/Gaussian
    # convention includes the angular normalization factor.  For each shell we
    # apply both reordering and renormalisation.
    #
    # Psi4 Cartesian d:  xx, xy, xz, yy, yz, zz
    # Molden standard d: xx, yy, zz, xy, xz, yz
    # Off-diagonal d (xy, xz, yz) need 1/√3 scaling.
    D_PERM = [0, 3, 5, 1, 2, 4]
    D_NORM = [1.0, 1.0, 1.0, 1.0 / np.sqrt(3.0),
              1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)]
    #
    # Psi4 Cartesian f:  xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
    # Molden standard f: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
    # Off-diagonal f (all but xxx/yyy/zzz) need 1/√(15) or 1/√(3) scaling
    # depending on the triple-exponent pattern.
    # (F-support included for forward-compatibility; cc-pVDZ does not have f.)
    F_PERM = [0, 6, 9, 3, 1, 2, 5, 8, 7, 4]
    F_NORM = [
        1.0, 1.0, 1.0,                              # xxx, yyy, zzz (diag)
        1.0 / np.sqrt(5.0),                          # xyy
        1.0 / np.sqrt(5.0),                          # xxy
        1.0 / np.sqrt(5.0),                          # xxz
        1.0 / np.sqrt(5.0),                          # xzz
        1.0 / np.sqrt(5.0),                          # yzz
        1.0 / np.sqrt(5.0),                          # yyz
        1.0 / np.sqrt(15.0),                         # xyz
    ]

    n_AO = C_AO.shape[0]
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

    # Keep everything before the [MO] section
    mo_tag = "[MO]"
    idx = text.find(mo_tag)
    if idx == -1:
        raise RuntimeError("Psi4 Molden output has no [MO] section")
    header = text[:idx]

    lines = [header + "\n[MO]\n"]

    for i in range(n_orb):
        ei = energies[i]
        oi = occ[i]
        lines.append(f" Sym= A\n Ene= {ei:15.10f}\n Spin= Alpha\n"
                     f" Occup= {oi:14.10f}\n")
        coeffs = C_AO[perm, i] * scale
        for j in range(n_AO):
            lines.append(f"  {j + 1:>4d}  {coeffs[j]:16.10f}\n")

    # Pad with dummy orbitals so Avogadro's MO slot count matches [GTO]
    for i in range(n_orb, n_AO):
        lines.append(f" Sym= A\n Ene= {0.0:15.10f}\n Spin= Alpha\n"
                     f" Occup= {0.0:14.10f}\n")
        for j in range(n_AO):
            lines.append(f"  {j + 1:>4d}  {0.0:16.10f}\n")

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
    from .config import load_config as _load_config

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

    _cfg = _load_config()
    basis = _option(options, "basis", _cfg.get("basis", "cc-pVDZ"))
    method = _option(options, "method", _cfg.get("method", "hf"))
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
    atom_of, am_of, func_n, func_dtype = _get_basis_maps(bas_min)  # (n_min,) each

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
    C_vir = Ca.np[:, nocc:]                         # (n_AO, n_vir)
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
    msg, labels, net_charges = _analyze_ibos(
        C_IAO_all, occ_all, energies_all, nocc,
        atom_of, am_of, func_n, func_dtype,
        elem, method, basis, ref,
    )
    cjson["atoms"]["partialCharges"] = [round(c, 4) for c in net_charges]
    # -- Total Wiberg bond order section -----------------------------------
    wiberg_section = _format_total_wiberg(
        C_IAO_all[:, :nocc], atom_of, elem
    )
    msg += wiberg_section
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
