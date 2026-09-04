"""
IAO/IBO construction and Pipek-Mezey localization for Avogadro 2.

BSD 3-Clause License
Copyright (c) 2025-2026, Billy Wayne McCann
SPDX-License-Identifier: BSD-3-Clause

References:
  G. Knizia, JCTC 2013, 9, 4834-4843.  DOI: 10.1021/ct400687b
  ("Intrinsic Atomic Orbitals: An Unbiased Bridge between Quantum
   Theory and Chemical Concepts.")
  W. D. Derricotte and F. A. Evangelista, JCTC 2017, 13, 5984-5999.
  DOI: 10.1021/acs.jctc.7b00493 ("Localized Intrinsic Valence Virtual
  Orbitals ...").  Valence-virtual construction (their eqs 3-5).

Paper equation numbers and appendix references refer to Knizia 2013
unless marked D&E2017.
"""

from dataclasses import dataclass

import numpy as np
import warnings

# 1 Hartree in electron-volts (CODATA 2018) and in kcal/mol
# (2625.4996394799 kJ/mol ÷ 4.184).  Only used for the human-facing
# HOMO-LUMO gap line; all internal energies stay in Ha.
HA_TO_EV = 27.211386245988
HA_TO_KCAL = 627.5094740631

# Floor for the smallest kept VVO singular value.  Anything below the
# floor means the "valence" label is suspect (pathological SCF basis,
# near-linear dependence).  Set 20x below the observed suite minima
# (all 1.000 to three decimals, cc-pVDZ through aug-cc-pVDZ) and 5x
# above junk scale (~0.01): it never fires on sane input.
VVO_MIN_SIGMA = 0.05


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
        1: ["px", "py", "pz"],  # Psi4 Cartesian p order: x, y, z
        2: ["dxx", "dxy", "dxz", "dyy", "dyz", "dzz"],  # Psi4 Cartesian d order
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

        subtypes = _P_SUBTYPE.get(am, [""] * nfunc)
        for f_idx in range(nfunc):
            atom_of.append(atom)
            am_of.append(am)
            n_of.append(shell_n)
            dtype_of.append(subtypes[f_idx] if f_idx < len(subtypes) else "")

    return (
        np.array(atom_of, dtype=np.int32),
        np.array(am_of, dtype=np.int32),
        np.array(n_of, dtype=np.int32),
        dtype_of,
    )


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
    try:
        L_S, low_S = cho_factor(S)
    except np.linalg.LinAlgError:
        S_work = S + np.eye(S.shape[0], dtype=np.float64) * 1e-12
        L_S, low_S = cho_factor(S_work)
    P12 = cho_solve((L_S, low_S), S12)  # (n_AO, n_min)

    # (2) Occupied MOs expressed in the minimal basis
    C_occ_min = S12.T @ C_occ  # (n_min, n_occ)

    # (3) Solve S_min @ C_tilde = C_occ_min
    try:
        L_min, low_min = cho_factor(S_min)
    except np.linalg.LinAlgError:
        S_min_work = S_min + np.eye(S_min.shape[0], dtype=np.float64) * 1e-12
        L_min, low_min = cho_factor(S_min_work)
    C_tilde = cho_solve((L_min, low_min), C_occ_min)  # (n_min, n_occ)

    # (4) Metric in the occupied space
    S_tilde = C_occ_min.T @ C_tilde  # (n_occ, n_occ)

    # (5) Solve S_tilde @ C_tilde_2bar^T = C_tilde^T
    try:
        L_tilde, low_tilde = cho_factor(S_tilde)
    except np.linalg.LinAlgError:
        S_tilde = S_tilde + np.eye(S_tilde.shape[0], dtype=np.float64) * 1e-12
        L_tilde, low_tilde = cho_factor(S_tilde)
    C_tilde_2bar_T = cho_solve((L_tilde, low_tilde), C_tilde.T)  # (n_occ, n_min)
    C_tilde_2bar = C_tilde_2bar_T.T  # (n_min, n_occ)

    # (6) Residual part of occupied MOs beyond the minimal projection
    T4 = C_occ - P12 @ C_tilde_2bar  # (n_AO, n_occ)

    # (7) Construct IAO coefficients
    C_IAO = P12 + T4 @ C_occ_min.T  # (n_AO, n_min)

    # (8) Symmetric (Loewdin) orthogonalisation of IAOs
    #     Find M^{-1/2} where M = C_IAO^T @ S @ C_IAO
    metric = C_IAO.T @ S @ C_IAO  # (n_min, n_min)
    evals, evecs = np.linalg.eigh(metric)
    evals = np.maximum(
        evals, 1e-14
    )  # guard against near-zero from near-linear-dependence
    C_IAO = C_IAO @ (evecs @ np.diag(evals**-0.5) @ evecs.T)

    # Express the occupied MOs in the orthonormal IAO basis.
    # Since IAOs span the occupied space (by construction),
    # C_IAO @ C_IAO_occ = C_occ should hold exactly.
    C_IAO_occ = C_IAO.T @ S @ C_occ  # (n_min, n_occ)

    return C_IAO, C_IAO_occ


# ---------------------------------------------------------------------------
# Pipek-Mezey localisation in the IAO basis   (eq 4 and Appendix D)
# ---------------------------------------------------------------------------


def _localize_ibos(
    C_occ, atom_of, max_iter=2048, conv=1e-12, exponents=(2, 4), cayley_deg=0.0, seed=42
):
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
        U = np.linalg.solve(np.eye(n_occ) - 0.5 * A, np.eye(n_occ) + 0.5 * A)
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
                        # Pipek-Mezey p=2 (Appendix D)
                        # A_ij = Σ_A [-2(q_ii² + q_jj²) + 4·q_ii·q_jj + 4·q_ij²]
                        # B_ij = Σ_A 4·q_ij·(q_ii - q_jj)
                        # φ = 0.25·atan2(B, -A)  [from tan(4φ) = B/-A]
                        Aij = np.sum(
                            -2.0 * Qii * Qii
                            - 2.0 * Qjj * Qjj
                            + 4.0 * Qii * Qjj
                            + 4.0 * Qij * Qij
                        )
                        Bij = np.sum(4.0 * Qij * (Qii - Qjj))
                        if abs(Aij) <= conv:
                            continue
                        phi = 0.25 * np.arctan2(Bij, -Aij)
                        grad_term = 2.0
                    elif exponent == 4:
                        # Pipek-Mezey p=4 (eq 4).  The published Appendix D
                        # 2x2 update formulas contain a production error
                        # (confirmed by Knizia at https://sites.psu.edu/knizia/software/).
                        # These formulas match the corrected reference
                        # implementation (ibo-ref).
                        qii_2 = Qii * Qii
                        qjj_2 = Qjj * Qjj
                        qij_2 = Qij * Qij
                        Aij = np.sum(
                            -qii_2 * qii_2
                            - qjj_2 * qjj_2
                            + 6.0 * (qii_2 + qjj_2) * qij_2
                            + qii_2 * Qii * Qjj
                            + Qii * qjj_2 * Qjj
                        )
                        Bij = np.sum(4.0 * Qij * (qii_2 * Qii - qjj_2 * Qjj))
                        if abs(Aij) <= conv:
                            continue
                        phi = 0.25 * np.arctan2(Bij, -Aij)
                        grad_term = 4.0
                    else:
                        raise ValueError(f"Unsupported PM exponent: {exponent}")

                    cs = np.cos(phi)
                    sn = np.sin(phi)

                    old_i = C_occ[:, i].copy()
                    old_j = C_occ[:, j].copy()
                    C_occ[:, i] = cs * old_i + sn * old_j
                    C_occ[:, j] = cs * old_j - sn * old_i

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

    sq = C_occ**2
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
        C_block = C_occ[:, indices]  # (n_IAO, n_g)
        Fb = C_block.T @ (F_IAO @ C_block)  # (n_g, n_g)
        evals, evecs = np.linalg.eigh(Fb)
        C_occ[:, indices] = C_block @ evecs


# ---------------------------------------------------------------------------
# Resolve bond-flat PM degeneracies (sigma/pi vs banana bonds)
# ---------------------------------------------------------------------------


def _resolve_flat_degeneracies(C_occ, atom_of, F_IAO, flat_tol=1e-6,
                               fock_tol=1e-8, pm_exponent=4):
    """
    Rotate PM-functionally-degenerate orbital pairs to their Fock-diagonal
    basis.

    The PM functional measures per-atom populations only.  When two orbitals
    share identical population vectors n_A(i) = n_A(j) on every atom — the
    {sigma, pi} plane of a symmetric two-centre bond, its antibond
    counterpart, or two orbitals on one atom — every rotation within the
    pair leaves L unchanged.  Jacobi sweeps neither prefer nor repair such
    mixtures, so on symmetry-broken geometries the converged picture
    (sigma+pi vs two banana bonds) is decided by the SCF-seeded trajectory,
    not by the functional.

    Detection mirrors the PM sweep itself: a pair is *flat* when rotating
    it by 45 degrees changes L by less than *flat_tol* relative to |L|.
    Among flat pairs, only *Fock-coupled* ones (|F_ij| > *fock_tol*) are
    touched; they are rotated by the minimal Jacobi angle that zeroes F_ij,
    so the aufbau (energy) ordering emerges without perturbing anything
    else.  Flat pairs that are already Fock-diagonal yield phi = 0 exactly
    and are left byte-identical.  Non-flat pairs (any real population
    asymmetry) sit in a steep PM bowl, never satisfy the tolerance, and are
    never modified — the same principle as _resolve_on_atom_mixing, extended
    from one atom to two.

    Parameters are modified in place.

    Returns
    -------
    int : number of pairs rotated.
    """
    n_IAO, n_occ = C_occ.shape
    if n_occ < 2:
        return 0
    n_atoms = int(np.max(atom_of)) + 1
    p = pm_exponent
    if p not in (2, 4):
        raise ValueError(f"Unsupported PM exponent: {p}")

    FC = F_IAO @ C_occ  # (n_IAO, n_occ), for Fock couplings

    def _pair_pops(a, b):
        pa = np.zeros(n_atoms)
        pb = np.zeros(n_atoms)
        np.add.at(pa, atom_of, a * a)
        np.add.at(pb, atom_of, b * b)
        return pa, pb

    n_rotated = 0
    for i in range(1, n_occ):
        ci = C_occ[:, i]
        Fi = FC[:, i]
        tii = float(ci.dot(Fi))
        for j in range(i):
            cj = C_occ[:, j]

            # PM functional value of the pair now and at a 45-degree
            # rotation.  Other orbitals are unaffected by the pair
            # rotation, so pair-only L decides flatness.
            pi, pj = _pair_pops(ci, cj)
            L0 = float(np.sum(pi**p) + np.sum(pj**p))
            if L0 <= 0.0:
                continue
            r2 = 1.0 / np.sqrt(2.0)
            pr, ps = _pair_pops(r2 * (ci + cj), r2 * (cj - ci))
            L45 = float(np.sum(pr**p) + np.sum(ps**p))

            # Converged PM makes every pair a local maximum; any change in
            # L beyond tolerance means the functional actively distinguishes
            # the pair (covers both directions in case PM exited at
            # max_iter before full convergence).
            if abs(L0 - L45) / abs(L0) > flat_tol:
                continue

            # Flat pair: break the tie only if the Fock matrix actually
            # couples the two states.  Purely relative tolerance: SCF dust
            # (~1e-8 relative) must not trigger a microscopic rotation,
            # while genuinely coupled flat pairs (ratio ~0.4 in the failing
            # ethene case) fire robustly.  Already-diagonal pairs stay
            # byte-identical.
            tij = float(ci.dot(FC[:, j]))
            tjj = float(cj.dot(FC[:, j]))
            if abs(tij) <= fock_tol * max(abs(tii), abs(tjj)):
                continue

            # Minimal 2x2 Jacobi rotation zeroing F_ij:
            #   <i'|F|j'> = cs (t_jj - t_ii) + (c^2 - s^2) t_ij = 0
            phi = 0.5 * np.arctan2(2.0 * tij, tii - tjj)
            c_, s_ = np.cos(phi), np.sin(phi)
            # Copy before writing: ci/cj and Fi are views into C_occ/FC;
            # the second assignment must read the ORIGINAL columns.
            old_i = ci.copy()
            old_j = cj.copy()
            old_fi = Fi.copy()
            old_fj = FC[:, j].copy()
            C_occ[:, i] = c_ * old_i + s_ * old_j
            C_occ[:, j] = -s_ * old_i + c_ * old_j
            FC[:, i] = c_ * old_fi + s_ * old_fj
            FC[:, j] = -s_ * old_fi + c_ * old_fj
            n_rotated += 1

    return n_rotated


# ---------------------------------------------------------------------------
# IBO analysis table
# ---------------------------------------------------------------------------

# Periodic-table lookup for element symbols
_ELEM_SYMBOLS = [
    "X",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
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
        "dxy": c_xy**2,
        "dxz": c_xz**2,
        "dyz": c_yz**2,
        "dz2": (-c_xx - c_yy + 2 * c_zz) ** 2,
        "dx2y2": (c_xx - c_yy) ** 2,
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
    for am_label, am_val in [("s", 0), ("p", 1), ("d", 2)]:
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
        dominant_n = max(n_counts, key=lambda k: n_counts[k])

        # Determine dominant subtype
        subtype = ""
        if am_val == 1:  # p-orbitals: px, py, pz
            st_counts = {}
            for fi in idx_am:
                st = func_dtype[fi]
                if st:
                    st_counts[st] = st_counts.get(st, 0.0) + c[fi] ** 2
            if st_counts:
                top_st = max(st_counts, key=lambda k: st_counts[k])
                if st_counts[top_st] > 0.5 * total_am:
                    subtype = top_st  # e.g. "pz"
        elif am_val == 2:  # d-orbitals: dxy, dxz, dyz, dz2, dx2y2
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
    return float(occ**2 * pop[A] * pop[B])


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


def _classify_orbital(oc, pop, order, top_A, top_B, s_char, p_char, d_char,
                      elem, am_of, atom_of, func_n, c):
    """Return a classification label string for one IAO-basis orbital.

    Classifies occupied orbitals as Core, LP, σ/π bond, 2e3c, or Deloc,
    and virtual orbitals as the corresponding antibond (*) type.

    Thresholds (DOM-based, matching IboView defaults):
      Core:  DOM > 0.99 + s-character > 0.75 on n=1 (1s only; 2s/3s are valence)
      LP:    DOM > 0.90 on one atom
      σ/π:   DOM_shared > 0.75, both atoms carry density (>0.02);
             π if p-fraction > 0.85 on both atoms
      LP-s:  DOM > 0.70, s-character > 0.5 (transitional)
      2e3c:  3rd atom carries >10% density, 4th atom <3%
      Deloc: everything else (multi-atom delocalisation)
      Virtual antibond thresholds follow the same logic with slightly
      looser cut-offs (0.60 vs 0.75 for shared density).
    """
    if oc > 1.5:
        # Determine principal quantum number of the dominant s-contribution
        # on the top atom.  Only n=1 (1s) is a true core; 2s/3s are valence.
        _s_idx = np.where((atom_of == top_A) & (am_of == 0))[0]
        _s_n = int(func_n[_s_idx[np.argmax(c[_s_idx] ** 2)]]) if len(_s_idx) else 0
        if pop[top_A] > 0.99 and s_char > 0.75 and _s_n == 1:
            return f"{_elem_symbol(elem[top_A])}(Core)"
        elif pop[top_A] > 0.90:
            return f"{_elem_symbol(elem[top_A])}(LP)"
        elif pop[top_A] + pop[top_B] > 0.75 and pop[top_B] > 0.02:
            pfrac_A = _p_frac(c, am_of, top_A, atom_of)
            pfrac_B = _p_frac(c, am_of, top_B, atom_of)
            bond_type = "π" if (pfrac_A > 0.85 and pfrac_B > 0.85) else "σ"
            a, b = sorted([top_A, top_B])
            symA = _elem_symbol(elem[a])
            symB = _elem_symbol(elem[b])
            return f"{symA}-{symB} {bond_type}"
        elif pop[top_A] > 0.70:
            symA = _elem_symbol(elem[top_A])
            if s_char > 0.5:
                return f"{symA}(LP-s)"
            else:
                return f"{symA}(LP)"
        else:
            if len(order) >= 3 and pop[order[2]] > 0.10 and (len(order) < 4 or pop[order[3]] <= 0.03):
                atoms = sorted(
                    [order[0], order[1], order[2]],
                    key=lambda i: (_elem_symbol(elem[i]), i),
                )
                syms = "-".join(_elem_symbol(elem[a]) for a in atoms)
                return f"{syms} 2e3c"
            else:
                return "Deloc"
    else:
        if pop[top_A] + pop[top_B] > 0.75 and pop[top_B] > 0.02:
            pfrac_A = _p_frac(c, am_of, top_A, atom_of)
            pfrac_B = _p_frac(c, am_of, top_B, atom_of)
            bond_type = "π" if (pfrac_A > 0.85 and pfrac_B > 0.85) else "σ"
            a, b = sorted([top_A, top_B])
            symA = _elem_symbol(elem[a])
            symB = _elem_symbol(elem[b])
            return f"{symA}-{symB} {bond_type}*"
        elif len(order) >= 3 and pop[order[2]] > 0.08:
            pfrac_top = _p_frac(c, am_of, top_A, atom_of)
            symA = _elem_symbol(elem[top_A])
            return f"{symA} π*" if pfrac_top > 0.85 else f"{symA} anti*"
        elif pop[top_A] + pop[top_B] > 0.60 and pop[top_B] > 0.02:
            a, b = sorted([top_A, top_B])
            symA = _elem_symbol(elem[a])
            symB = _elem_symbol(elem[b])
            return f"{symA}-{symB} anti*"
        elif pop[top_A] > 0.50:
            return f"{_elem_symbol(elem[top_A])}(virt)"
        else:
            return "Virt"


def _analyze_ibos(
    C_IAO_all,
    occ_all,
    energies_all,
    nocc,
    atom_of,
    am_of,
    func_n,
    func_dtype,
    elem,
    method,
    basis,
    ref,
    mol_name="",
):
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
    atom_pop = np.zeros(
        n_atoms, dtype=np.float64
    )  # accumulated per-atom electron counts
    lines.append(f"IBO Analysis: {mol_name}  ({method}/{basis}, {ref.upper()})")
    lines.append("")

    # Pre-pass to size the Composition column to fit the widest entry
    comp_width = len("Composition")
    for orb in range(n_orb):
        sq = C_IAO_all[:, orb] ** 2
        pop = np.zeros(n_atoms, dtype=np.float64)
        np.add.at(pop, atom_of, sq)
        order = np.argsort(-pop)
        comp_parts = []
        for A in order[:4]:
            if pop[A] > 0.005:
                sym = f"{_elem_symbol(elem[A])}{A + 1}"
                pct = pop[A] * 100.0
                comp_parts.append(f"{sym}({pct:.1f}%)")
        comp = " + ".join(comp_parts)
        comp_width = max(comp_width, len(comp))
    comp_width += 2

    header = (
        f"  {'#':>3}  {'Occ':>7}  {'Energy':>10}  "
        f"{'Type':>16}  {{:<{comp_width}}}  {'Hybrid':<22}  "
        f"{'Ion%':>5}  {'H/L':>9}"
    ).format("Composition")
    lines.append(header)
    lines.append("-" * len(header))

    # Identify degenerate manifolds: groups of consecutive orbitals
    # with energy differences < 1e-4 Ha.  Within symmetric molecules,
    # the PM functional leaves symmetry-equivalent bonds with small
    # residual energy splittings (~1e-5 Ha for benzene C-H σ).
    DEG_THRESH = 2e-4  # ~0.13 kcal/mol; catches all PM convergence noise
    deg_ranges = []
    is_degen = np.zeros(n_orb, dtype=bool)
    group_start = None
    for i in range(n_orb):
        if i > 0 and abs(energies_all[i] - energies_all[i - 1]) < DEG_THRESH:
            if group_start is None:
                group_start = i - 1
            is_degen[i] = True
        else:
            if group_start is not None and i - group_start > 1:
                deg_ranges.append((group_start, i))
                for j in range(group_start, i):
                    is_degen[j] = True
            group_start = None
    if group_start is not None and n_orb - group_start > 1:
        deg_ranges.append((group_start, n_orb))
        for j in range(group_start, n_orb):
            is_degen[j] = True

    for orb in range(n_orb):
        oc = occ_all[orb]
        sq = C_IAO_all[:, orb] ** 2
        pop = np.zeros(n_atoms, dtype=np.float64)
        np.add.at(pop, atom_of, sq)
        if oc > 1.5:
            # RHF-only: occ=2.0 per occupied orbital.  UHF would need
            # separate alpha/beta occupancy arrays — unreachable due
            # to the closed-shell guard in compute_ibo.
            atom_pop += pop * oc  # accumulate electron count per atom

        # Dominant atom and its population
        order = np.argsort(-pop)
        top_A = order[0]
        top_B = order[1]

        # s/p/d breakdown on the dominant atom (single pass)
        s_char, p_char, d_char = _spd_frac(C_IAO_all[:, orb], am_of, top_A, atom_of)

        # Determine orbital type
        orbid = _classify_orbital(oc, pop, order, top_A, top_B, s_char, p_char,
                                  d_char, elem, am_of, atom_of, func_n,
                                  C_IAO_all[:, orb])

        orbid_labels[orb] = orbid

        comp_parts = []
        for A in order[:4]:
            if pop[A] > 0.005:
                sym = f"{_elem_symbol(elem[A])}{A + 1}"
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
        hybrid = _hybrid_str(
            C_IAO_all[:, orb], am_of, atom_of, func_n, func_dtype, hybrid_atom
        )
        if hybrid_atom != top_A:
            hybrid = f"{_elem_symbol(elem[hybrid_atom])}: {hybrid}"

        # Per-IBO Wiberg bond order and ionic character (between top_A, top_B)
        w_ab = _wiberg_per_ibo(pop, oc, top_A, top_B)
        if oc > 1.5 and w_ab > 0.001:
            ion_str = f"{_ionic_pct(pop, top_A, top_B):.1f}"
        else:
            ion_str = "---"

        hl = ""
        if orb == nocc - 1:
            hl = "<- HOMO"
        elif orb == nocc:
            hl = "<- LUMO"
        degen_tag = " †" if is_degen[orb] else ""
        lines.append(
            f"  {orb + 1:>3d}  {oc:>7.3f}  {energies_all[orb]:>10.6f}  "
            f"{orbid:>{16 - len(degen_tag)}}{degen_tag}  "
            f"{comp:<{comp_width}}  {hybrid:<22}  "
            f"{ion_str:>5}  {hl:>9}"
        )

    # Footnote for degenerate manifolds
    if deg_ranges:
        deg_groups = []
        for start, end in deg_ranges:
            deg_groups.append(f"{start + 1}-{end}")
        lines.append(
            f"  † Orbitals {', '.join(deg_groups)} form degenerate manifolds (ΔE < {DEG_THRESH:.0e} Ha).\n"
            f"  Small energy differences within each manifold are PM convergence noise and do not indicate true energy splittings."
        )

    lines.append("")
    # RHF-only: total = 2 × nocc.  The UHF branch (else nocc) is
    # unreachable — closed-shell guard prevents open-shell in compute_ibo.
    lines.append(f"Total electrons: {int(2 * nocc) if ref == 'rhf' else nocc}")

    # Frontier-orbital summary.  HOMO/LUMO are found by occupancy, not by
    # rank position: energies_all is ascending but not guaranteed to be a
    # strict occupied-then-virtual prefix after sorting (see IBOResult).
    occ_idx = np.where(occ_all > 1.5)[0]
    vir_idx = np.where(occ_all < 0.5)[0]
    if len(occ_idx) and len(vir_idx):
        homo_i = int(occ_idx[np.argmax(energies_all[occ_idx])])
        lumo_i = int(vir_idx[np.argmin(energies_all[vir_idx])])
        gap_ha = energies_all[lumo_i] - energies_all[homo_i]
        gap_ev = gap_ha * HA_TO_EV
        gap_kcal = gap_ha * HA_TO_KCAL
        lines.append("")
        lines.append("--- Frontier Orbital Energies ---")
        lines.append(
            f"  HOMO ({orbid_labels[homo_i]}, orb {homo_i + 1}): {energies_all[homo_i]:>10.6f} Ha"
        )
        lines.append(
            f"  LUMO ({orbid_labels[lumo_i]}, orb {lumo_i + 1}): {energies_all[lumo_i]:>10.6f} Ha"
        )
        lines.append(
            f"  HOMO-LUMO gap: {gap_ha:>10.6f} Ha = {gap_ev:>7.3f} eV = {gap_kcal:>8.1f} kcal/mol"
        )

    charge_section = _format_charge_decomposition(atom_pop, elem)
    lines.append(charge_section)

    net_charges = [float(int(round(elem[A])) - atom_pop[A]) for A in range(n_atoms)]
    return "\n".join(lines), orbid_labels, net_charges


# Orbital-pair interference terms with |term| at or above this value get
# their own lines in the "Significant orbital-pair interference" detail
# section.  At 0.01 the section fires only on genuine delocalization
# chemistry (3c-2e bridges, conjugated π) and stays silent on ordinary
# single/double bonds.
PAIR_DETAIL_THRESH = 0.01


def _format_wiberg(C_IAO_occ, atom_of, am_of, elem, labels=None):
    """Single Wiberg table: exact density-matrix totals decomposed σ/π.

    The density Wiberg (as originally reported) is

        W_AB = Σ_{i∈A,j∈B} D²_ij,   D = 2·C_occ·C_occᵀ  (RHF)

    Expanding the square over orbitals gives an exact sum over orbital
    pairs,

        W_AB = 4 Σ_{k,l} G^A_kl G^B_kl,   G^A_kl = Σ_{i∈A} c_ki c_li

    with diagonal terms (k=l) the per-IBO shares occ²·P_A·P_B and
    off-diagonal terms the inter-orbital interference.  Orbitals are
    classed σ or π purely by p-fraction (p > 0.85 on both dominant
    atoms → π, no population gate), and the interference is folded
    into its class: (σ,σ) → σ, (π,π) → π, (σ,π) → split 50/50.
    Total = σ + π exactly; the folded interference is echoed in a
    parenthesised column for transparency.

    A follow-on detail section lists individual orbital-pair
    interference terms with |term| >= PAIR_DETAIL_THRESH, grouped by
    bond in table order, so the reader can see which orbitals drive a
    bond's interference.  ``labels`` supplies the per-orbital table
    labels (column k ↔ analysis-table row k+1); when omitted, bare orb
    numbers are shown.  Bonds whose every pair term falls below the
    threshold get no detail lines, and the section is omitted entirely
    when empty.
    """
    n_occ = C_IAO_occ.shape[1]
    n_atoms = len(elem)

    # Classify each occupied orbital σ/π (p-fraction rule only).
    is_pi = np.zeros(n_occ, dtype=bool)
    for k in range(n_occ):
        c = C_IAO_occ[:, k]
        sq = c**2
        pop = np.zeros(n_atoms, dtype=np.float64)
        np.add.at(pop, atom_of, sq)
        order = np.argsort(-pop)
        A, B = int(order[0]), int(order[1])
        pa = _p_frac(c, am_of, A, atom_of)
        pb = _p_frac(c, am_of, B, atom_of)
        is_pi[k] = pa > 0.85 and pb > 0.85

    # Per-atom orbital-pair overlap blocks G^A_kl.
    G = np.zeros((n_atoms, n_occ, n_occ))
    for a in range(n_atoms):
        Ca = C_IAO_occ[atom_of == a, :]
        G[a] = Ca.T @ Ca

    sigma = np.zeros((n_atoms, n_atoms), dtype=np.float64)
    pi = np.zeros((n_atoms, n_atoms), dtype=np.float64)
    # Folded interference, tracked per class: int_sigma = σσ + ½σπ,
    # int_pi = ππ + ½σπ.  These are what the σ/π columns actually
    # contain beyond their diagonal shares (4·P_A·P_B).
    int_sigma = np.zeros((n_atoms, n_atoms), dtype=np.float64)
    int_pi = np.zeros((n_atoms, n_atoms), dtype=np.float64)
    # Significant off-diagonal terms, kept for the detail section:
    # (A, B) with A < B -> [(k, l, value)] with |value| >= threshold.
    pair_detail = {}

    for k in range(n_occ):
        # Diagonal (per-orbital share): 4·P_A·P_B
        contrib = 4.0 * np.einsum("a,b->ab", G[:, k, k], G[:, k, k])
        if is_pi[k]:
            pi += contrib
        else:
            sigma += contrib
        for l in range(k + 1, n_occ):
            # Off-diagonal (interference): 8·G^A_kl·G^B_kl
            contrib = 8.0 * np.einsum("a,b->ab", G[:, k, l], G[:, k, l])
            if is_pi[k] and is_pi[l]:
                pi += contrib
                int_pi += contrib
            elif is_pi[k] or is_pi[l]:
                # σ-π cross term: split evenly (contribution is symmetric)
                half = 0.5 * contrib
                sigma += half
                pi += half
                int_sigma += half
                int_pi += half
            else:
                sigma += contrib
                int_sigma += contrib
            # Record significant pair terms for the detail section.
            big = np.abs(contrib) >= PAIR_DETAIL_THRESH
            if big.any():
                for A in range(n_atoms):
                    for B in range(A + 1, n_atoms):
                        if big[A, B]:
                            pair_detail.setdefault((A, B), []).append(
                                (k, l, float(contrib[A, B]))
                            )

    rows = []
    total = sigma + pi
    for A in range(n_atoms):
        for B in range(A + 1, n_atoms):
            if total[A, B] > 0.01:
                symA = _elem_symbol(elem[A])
                symB = _elem_symbol(elem[B])
                rows.append(
                    (symA, A, symB, B, total[A, B], sigma[A, B], pi[A, B],
                     int_sigma[A, B], int_pi[A, B])
                )

    if not rows:
        return ""

    lines = [
        "",
        "",
        "--- Wiberg Bond Orders (σ/π, density) ---",
        "  W_AB = Σ_{i∈A,j∈B} D²_ij (density Wiberg); σ + π = total exactly.",
        "  σ, π columns include their class's folded interference; the",
        "  parenthetical reports the same interference as (σ-part, π-part).",
        f"  {'Bond':<10}{'Total':>8}{'σ':>8}{'π':>8}{'  (interference)':>28}",
    ]
    for symA, a, symB, b, t, s, p, is_, ip in sorted(rows, key=lambda x: -x[4]):
        # Kill floating-point -0.000 noise in the σ/π columns (the
        # interference parts keep their genuine sign).
        s_disp = 0.0 if abs(s) < 5e-4 else s
        p_disp = 0.0 if abs(p) < 5e-4 else p
        if abs(is_) + abs(ip) >= 5e-4:
            # Both parts below the noise floor -> omit the parenthetical
            # entirely; the row is pure diagonal shares.
            is_disp = 0.0 if abs(is_) < 5e-4 else is_
            ip_disp = 0.0 if abs(ip) < 5e-4 else ip
            lines.append(
                f"  {symA}{a+1}-{symB}{b+1:<7}{t:>8.3f}{s_disp:>8.3f}{p_disp:>8.3f}"
                f"  ({is_disp+ip_disp:+.3f}: σ{is_disp:+.3f}, π{ip_disp:+.3f})"
            )
        else:
            lines.append(
                f"  {symA}{a+1}-{symB}{b+1:<7}{t:>8.3f}{s_disp:>8.3f}{p_disp:>8.3f}"
            )

    # Detail section: significant orbital-pair interference terms, grouped
    # by bond in table order.  Omitted entirely when nothing clears the bar
    # (ordinary single/double bonds), so quiet molecules gain no lines.
    def _orb_tag(k):
        if labels is not None and k < len(labels) and labels[k]:
            return f"orb{k+1}({labels[k]})"
        return f"orb{k+1}"

    detail = []
    for symA, a, symB, b, t, s, p, is_, ip in sorted(rows, key=lambda x: -x[4]):
        terms = pair_detail.get((a, b))
        if not terms:
            continue
        for k, l, v in sorted(terms, key=lambda t_: -abs(t_[2])):
            detail.append(
                f"  {symA}{a+1}-{symB}{b+1}: {_orb_tag(k)} × {_orb_tag(l)}: {v:+.4f}"
            )
    if detail:
        lines.append("")
        lines.append(
            "--- Significant orbital-pair interference "
            f"(|term| ≥ {PAIR_DETAIL_THRESH:g}) ---"
        )
        lines.append("  Positive pair terms add to the bond order; negative pair")
        lines.append("  terms subtract from it. Signs are relative to the listed")
        lines.append("  bond — the same pair may contribute oppositely elsewhere.")
        lines.append("  Orbitals are numbered as in the analysis table above.")
        lines.extend(detail)
    return "\n".join(lines)


def _format_charge_decomposition(atom_pop, elem):
    """
    Format a charge decomposition table from accumulated IAO populations.

    For each atom A:
        Q_A = Σ_k occ_k · P_A^{(k)}   (total electrons on atom A)
        Net charge = Z_A - Q_A

    ``atom_pop[A]`` is Q_A as a float.  ``elem`` gives atomic numbers.
    RHF-only: assumes occ=2.0 per occupied orbital.  UHF would need
    separate alpha/beta populations — unreachable due to the
    closed-shell guard in compute_ibo.
    """
    lines = ["", "--- Charge Decomposition ---"]
    header = f"  {'Atom':>5}  {'Z':>3}  {'Pop':>8}  {'Net Charge':>10}"
    lines.append(header)
    lines.append("-" * len(header))
    total_pop = 0.0
    total_z = 0
    for A in range(len(elem)):
        Z = int(round(elem[A]))
        pop = atom_pop[A]
        net = Z - pop
        sym = _elem_symbol(Z)
        lines.append(f"  {sym}{A+1:<3}  {Z:>3d}  {pop:>8.3f}  {net:>+10.3f}")
        total_pop += pop
        total_z += Z
    lines.append("-" * len(header))
    lines.append(
        f"Total:  {total_z:>3d}  {total_pop:>8.3f}  {total_z - total_pop:>+10.3f}"
    )
    return "\n".join(lines)


def _elem_symbol(Z):
    Z = int(round(Z))
    if Z < len(_ELEM_SYMBOLS):
        return _ELEM_SYMBOLS[Z]
    return f"E{Z}"


def _spd_frac(c, am_of, atom, atom_of):
    """Return (s_char, p_char, d_char) for the given atom in one pass."""
    idx = np.where(atom_of == atom)[0]
    if len(idx) == 0:
        return 0.0, 0.0, 0.0
    am_atom = am_of[idx]
    c_atom = c[idx] ** 2
    s = float(np.sum(c_atom[am_atom == 0])) if np.any(am_atom == 0) else 0.0
    p = float(np.sum(c_atom[am_atom == 1])) if np.any(am_atom == 1) else 0.0
    d = float(np.sum(c_atom[am_atom == 2])) if np.any(am_atom == 2) else 0.0
    return s, p, d


def _p_frac(c, am_of, atom, atom_of):
    """p/s/d ratio on the given atom; 0 if no density."""
    s, p, d = _spd_frac(c, am_of, atom, atom_of)
    total = s + p + d
    return p / total if total > 0 else 0.0


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
    # NOTE: The [GTO] header is copied from Psi4's own molden() output, so
    # the primitive coefficients come from Psi4 (original_coef convention).
    # If any future code path writes [GTO] directly using shell.coef(p),
    # note that shell.coef(p) includes primitive normalization — the Molden
    # reader would re-apply normalization, producing incorrect basis
    # function values.  Use shell.original_coef(p) for direct GTO output.
    import psi4

    tmp = path.with_suffix(".molden.tmp")
    psi4.molden(wfn, str(tmp))
    text = tmp.read_text(encoding="utf-8")
    tmp.unlink()

    # Build index permutation and re-scaling vector to convert Psi4 internal AO
    # order to Molden standard.  Psi4's CCA convention uses unnormalized Cartesian
    # Gaussians (off-diagonal d/f have self-overlap < 1); the Molden/Gaussian
    # convention includes the angular normalization factor.
    #
    # Psi4 Cartesian d:  xx, xy, xz, yy, yz, zz
    # Molden standard d: xx, yy, zz, xy, xz, yz
    # Off-diagonal d (xy, xz, yz) need 1/√3 scaling.
    D_PERM = [0, 3, 5, 1, 2, 4]
    #
    # Psi4 Cartesian f:  xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
    # Molden standard f: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
    # (F-support included for forward-compatibility; cc-pVDZ does not have f.)
    F_PERM = [0, 6, 9, 3, 1, 2, 5, 8, 7, 4]

    n_AO = C_AO.shape[0]
    perm = np.arange(n_AO)
    scale = np.ones(n_AO)
    bas = wfn.basisset()
    ao = 0
    for sh in range(bas.nshell()):
        am = bas.shell(sh).am
        nf = bas.shell(sh).nfunction
        if am == 2:
            # Psi4 d order: xx, xy, xz, yy, yz, zz
            # Diagonal (xx/yy/zz): scale 1.0; off-diagonal: 1/3
            d_norm_in = [
                1.0,
                1.0 / np.sqrt(3),
                1.0 / np.sqrt(3),
                1.0,
                1.0 / np.sqrt(3),
                1.0,
            ]
            for i in range(6):
                perm[ao + i] = ao + D_PERM[i]
                scale[ao + i] = d_norm_in[i]
        elif am == 3:
            # Psi4 f order: xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
            f_norm_in = [
                1.0,
                1.0 / np.sqrt(5),
                1.0 / np.sqrt(5),
                1.0 / np.sqrt(5),
                1.0 / np.sqrt(15),
                1.0 / np.sqrt(5),
                1.0,
                1.0 / np.sqrt(5),
                1.0 / np.sqrt(5),
                1.0,
            ]
            for i in range(10):
                perm[ao + i] = ao + F_PERM[i]
                scale[ao + i] = f_norm_in[i]
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
        lines.append(f" Sym= A\n Ene= {ei:15.10f}\n Spin= Alpha\n Occup= {oi:14.10f}\n")
        coeffs = (C_AO[:, i] * scale)[perm]
        for j in range(n_AO):
            lines.append(f"  {j + 1:>4d}  {coeffs[j]:16.10f}\n")

    # Pad with dummy orbitals so Avogadro's MO slot count matches [GTO]
    for i in range(n_orb, n_AO):
        lines.append(
            f" Sym= A\n Ene= {0.0:15.10f}\n Spin= Alpha\n Occup= {0.0:14.10f}\n"
        )
        for j in range(n_AO):
            lines.append(f"  {j + 1:>4d}  {0.0:16.10f}\n")

    path.write_text("".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Molecule name helpers
# ---------------------------------------------------------------------------


def _mol_formula(numbers):
    """Molecular formula from atomic number list, preserving first-occurrence order."""
    from collections import Counter

    counts = Counter(numbers)
    seen = set()
    parts = []
    for Z in numbers:
        if Z not in seen:
            seen.add(Z)
            c = counts[Z]
            parts.append(f"{_ELEM_SYMBOLS[Z]}{c if c > 1 else ''}")
    return "".join(parts)


def _sanitize_name(name):
    """Sanitize a string for use as a filesystem directory name."""
    import re

    s = re.sub(r"[^a-zA-Z0-9]", "_", str(name))
    s = re.sub(r"_+", "_", s).strip("_")
    return s[:50] or "molecule"


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------


def _option(options, key, default):
    v = options.get(key, default)
    if isinstance(v, str):
        return v.strip()
    return v


def _write_input_xyz(path, coords, elem, mol_name):
    """Write an XYZ file from CJSON-style coordinate/element arrays."""
    n_atoms = len(elem)
    lines = [f"{n_atoms}\n", f"{mol_name}\n"]
    for i in range(n_atoms):
        sym = _elem_symbol(elem[i])
        lines.append(
            f"{sym:<3s}  {coords[3*i]:12.8f}  {coords[3*i+1]:12.8f}  {coords[3*i+2]:12.8f}\n"
        )
    path.write_text("".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Typed API (issue #3): pure core + IBOResult
# ---------------------------------------------------------------------------


@dataclass
class IBOResult:
    """Typed result of the IAO/IBO pipeline (see ``compute_ibo_data``).

    Arrays are plain numpy arrays; no files are written and no persistent
    configuration is consulted by the producing function.

    Attributes
    ----------
    C_IAO : (n_AO, n_min) array
        Orthonormal IAO coefficients in the AO basis (Knizia App. C).
    C_IAO_all : (n_min, n_orb) array
        Localized occupied + valence-virtual orbitals in the IAO basis,
        sorted ascending by energy.
    C_AO_all : (n_AO, n_orb) array
        Full-AO projection, ``C_IAO @ C_IAO_all``.
    occupations : (n_orb,) array
        2.0 for occupied orbitals, 0.0 for valence virtuals.
    energies : (n_orb,) array
        Fock-diagonal orbital energies (Hartree), ascending.
    atom_of, am_of, func_n : (n_min,) arrays
        Per-IAO atom index, angular momentum, principal quantum number.
    func_dtype : list[str]
        Per-IAO subtype label ("px"/"py"/"pz"/"dxx"/..., "" for s).
    elements : list[int]
        Atomic numbers, input order.
    coords : list[float]
        Flattened xyz coordinates, input order and units.
    n_occ : int
        Number of occupied orbitals (occupations == 2.0 identifies them;
        energies are ascending, so they are not a prefix after sorting).
    partial_charges : (n_atoms,) array
        IBO charge-decomposition net charges.
    labels : list[str]
        Per-orbital classification labels from the analysis pass.
    analysis_text : str
        Full human-readable analysis table (incl. Wiberg section) as
        written to ``ibos.txt`` by the renderer.
    method, basis, reference : str
        SCF level used.
    mol_name : str
        Display name (cjson name or formula fallback), truncated to 50.
    mol_spec : str
        Psi4 geometry specification string (charge/spin + xyz +
        no_com/no_reorient), usable to rebuild the molecule.
    scf_energy : float
        Total SCF energy (Hartree).
    wfn : object
        Opaque Psi4 wavefunction handle for renderers that need the
        working-basis [GTO] section or canonical MOs.  Not part of the
        stable typed surface.
    """

    C_IAO: np.ndarray
    C_IAO_all: np.ndarray
    C_AO_all: np.ndarray
    occupations: np.ndarray
    energies: np.ndarray
    atom_of: np.ndarray
    am_of: np.ndarray
    func_n: np.ndarray
    func_dtype: list
    elements: list
    coords: list
    n_occ: int
    partial_charges: np.ndarray
    labels: list
    analysis_text: str
    method: str
    basis: str
    reference: str
    mol_name: str
    mol_spec: str
    scf_energy: float
    wfn: object


def compute_ibo_data(cjson, options, charge=0, spin=1, psi4_output=None):
    """Pure typed core of the IBO pipeline.

    Runs the SCF (Psi4, in-process), builds the IAO basis, localizes
    occupied and valence-virtual blocks (Pipek-Mezey p=2 -> p=4),
    applies both degeneracy resolutions, and performs the composition
    analysis.  Writes no project files, touches no persistent
    configuration, and installs no logging handlers -- those are
    renderer/adapter concerns layered on top by :func:`compute_ibo`.

    Parameters mirror :func:`compute_ibo`: ``cjson`` supplies geometry
    (and optionally charge/spin via ``properties``); ``options`` may carry
    ``method``, ``basis``, ``iboview_style``.  Unspecified options fall
    back to library defaults (hf/cc-pVDZ) rather than the user config
    file -- callers wanting config persistence should merge it in
    beforehand.

    Psi4 routes its primary output through process-global state; without
    a destination it fails on Windows ("PSIOManager cannot get a mirror
    file handle").  Pass ``psi4_output`` to control that destination;
    by default a private temporary file is used so library callers need
    not care.

    Returns
    -------
    IBOResult

    Raises
    ------
    ValueError
        On open-shell input (the pipeline is RHF-only).
    RuntimeError
        If the Psi4 SCF fails.
    """
    atoms = cjson["atoms"]
    coords_raw = atoms["coords"]
    coords = coords_raw["3d"] if isinstance(coords_raw, dict) else coords_raw
    elem = atoms["elements"]["number"]

    elem_raw = cjson.get("atoms", {}).get("elements", {}).get("number", [])
    mol_name = cjson.get("name", "") or _mol_formula(elem_raw) or "molecule"
    mol_name = mol_name[:50]

    charge_val = int(cjson.get("properties", {}).get("totalCharge", charge))
    spin_val = int(cjson.get("properties", {}).get("totalSpinMultiplicity", spin))
    if spin_val != 1:
        raise ValueError(
            f"Open-shell systems are not supported (spin multiplicity "
            f"= {spin_val}). The IAO pipeline is RHF-only — all "
            f"occupied MOs are treated as doubly occupied, beta spin "
            f"is ignored, and charge/spin decomposition would be "
            f"incorrect."
        )
    ref = "rhf"

    geom_lines = "\n".join(
        f"  {elem[i]:3d}  {coords[3 * i]:12.8f}  {coords[3 * i + 1]:12.8f} "
        f"{coords[3 * i + 2]:12.8f}"
        for i in range(len(elem))
    )
    mol_spec = (
        f"{charge_val} {spin_val}\n{geom_lines}\nno_com\nno_reorient"
    )

    import psi4

    # Psi4's primary-output routing is process-global; without a
    # destination the PSIO manager fails on Windows.  Route to a private
    # temp file by default (the adapter passes its calc-dir log instead).
    if psi4_output is None:
        import tempfile

        _tmp_out = tempfile.NamedTemporaryFile(
            prefix="avo_ibo_psi4_", suffix=".log", delete=False
        )
        _tmp_out.close()
        psi4_output = _tmp_out.name
    psi4.set_output_file(str(psi4_output), append=True)

    # Register the molecule as Psi4's active geometry and force C1 so the
    # SCF, the minimal-basis build, and the IBO pipeline all see identical
    # AO orderings (reset_point_group also reorders shells otherwise).
    mol = psi4.geometry(mol_spec)
    mol.reset_point_group("c1")

    basis = _option(options, "basis", "cc-pVDZ")
    method = _option(options, "method", "hf")
    psi4.set_options(
        {
            "basis": basis,
            "scf_type": "df",
            "reference": ref,
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
            "puream": 0,
        }
    )
    # NOTE: puream=0 gives Cartesian basis functions, which is what the
    # paper assumes.  Changing this would affect the IAO construction.
    try:
        scf_energy, wfn = psi4.energy(method, return_wfn=True)
    except Exception as e:
        raise RuntimeError(f"Psi4 SCF failed for {method}/{basis}.") from e

    # -- Extract occupied coefficients and overlap matrices ----------------
    Ca = wfn.Ca()
    nocc = wfn.doccpi()[0] + wfn.soccpi()[0]
    mints = psi4.core.MintsHelper(wfn.basisset())

    S_full = mints.ao_overlap().np
    bas_min = psi4.core.BasisSet.build(mol, "BASIS", "STO-3G", puream=0)
    S_min = mints.ao_overlap(bas_min, bas_min).np
    S12 = mints.ao_overlap(wfn.basisset(), bas_min).np

    C_occ = Ca.np[:, :nocc].copy()  # (n_AO, n_occ)

    # -- Build IAO basis (Appendix C) --------------------------------------
    C_IAO, C_IAO_occ = _build_iao_basis(S_full, S12, S_min, C_occ)

    atom_of, am_of, func_n, func_dtype = _get_basis_maps(bas_min)

    # -- Pipek-Mezey localisation in IAO basis (eq 4 / Appendix D) --------
    _localize_ibos(C_IAO_occ, atom_of, max_iter=2048, conv=1e-12)

    # -- Compute orbital energies from Fock matrix -------------------------
    F_AO = wfn.Fa().np  # (n_AO, n_AO)
    F_IAO = C_IAO.T @ F_AO @ C_IAO  # (n_min, n_min)

    # -- Resolve on-atom degeneracies that PM cannot separate --------------
    # PM cannot separate orbitals on the same atom with DOM ~ 1 (e.g. O 2s
    # vs lone pair); Fock-diagonalise within each such subspace.
    _resolve_on_atom_mixing(C_IAO_occ, atom_of, F_IAO)

    # -- Resolve bond-flat PM degeneracies (sigma/pi vs banana bonds) ------
    # PM cannot distinguish orbitals sharing identical per-atom populations
    # — the {sigma, pi} plane of a symmetric bond.  See NOTES.md.
    _resolve_flat_degeneracies(C_IAO_occ, atom_of, F_IAO)

    occ_energies = np.array(
        [C_IAO_occ[:, i].dot(F_IAO @ C_IAO_occ[:, i]) for i in range(nocc)]
    )

    # -- Valence virtuals via SVD (D&E2017, eqs 3-5) ----------------------
    # S^IbVir_aρ = <φ_a|ψ_ρ>: canonical virtuals against the IAOs; the SVD
    # brings the two spaces into maximum coincidence, and the VVOs are the
    # first N_VVO = n_min − n_occ columns of U.  The count is structural:
    # a σ-threshold (as in IboView's MakeValenceVirtuals) can in principle
    # admit near-null squatters where diffuse manifolds overlap the IAO
    # space; the count rules that out by construction.  (No threshold
    # over-keep has been observed in any tested regime — cc-pVDZ through
    # aug-cc-pVDZ — so this is paper-parity plus guarantee, not a repair.)
    C_vir = Ca.np[:, nocc:]  # (n_AO, n_vir)
    SIbVir = C_IAO.T @ S_full @ C_vir  # (n_min, n_vir)
    U_svd, Sigma, _ = np.linalg.svd(SIbVir, full_matrices=False)
    n_val_vir = C_IAO.shape[1] - nocc
    if n_val_vir > U_svd.shape[1]:
        # SCF basis smaller than the minimal basis: not enough virtuals
        # to span the valence complement; keep what exists.
        warnings.warn(
            f"VVO count {n_val_vir} exceeds available "
            f"{U_svd.shape[1]} virtuals; keeping all"
        )
        n_val_vir = U_svd.shape[1]
    elif Sigma[n_val_vir - 1] < VVO_MIN_SIGMA:
        warnings.warn(
            f"smallest kept VVO singular value "
            f"{Sigma[n_val_vir - 1]:.3f} < {VVO_MIN_SIGMA}"
        )
    U_val = U_svd[:, :n_val_vir]  # (n_min, n_val_vir)

    # -- Localize the virtual block too (IboView localizes ALL case blocks) ---
    if n_val_vir > 1:
        _localize_ibos(U_val, atom_of, max_iter=2048, conv=1e-12)

    # -- Resolve bond-flat degeneracies in the virtual block -----------------
    # Same resolver as the occupied block (it is block-agnostic); safe here
    # because the fixed-count VVO construction above admits no junk
    # columns.  Distorted geometries yield σ*+π* instead of two
    # σ*-mixtures; equilibrium geometries are already Fock-diagonal and
    # stay byte-identical.
    _resolve_flat_degeneracies(U_val, atom_of, F_IAO)

    vir_energies = np.array(
        [U_val[:, i].dot(F_IAO @ U_val[:, i]) for i in range(n_val_vir)]
    )

    # -- Combined IAO-basis orbital set, sorted by energy ------------------
    C_IAO_all = np.hstack([C_IAO_occ, U_val])  # (n_min, n_orb)
    occ_all = np.array([2.0] * nocc + [0.0] * n_val_vir)
    energies_all = np.concatenate([occ_energies, vir_energies])

    order = np.argsort(energies_all)
    C_IAO_all = C_IAO_all[:, order]
    occ_all = occ_all[order]
    energies_all = energies_all[order]
    C_AO_all = C_IAO @ C_IAO_all  # (n_AO, n_orb)

    # -- Composition analysis ----------------------------------------------
    msg, labels, net_charges = _analyze_ibos(
        C_IAO_all,
        occ_all,
        energies_all,
        nocc,
        atom_of,
        am_of,
        func_n,
        func_dtype,
        elem,
        method,
        basis,
        ref,
        mol_name,
    )
    msg += _format_wiberg(C_IAO_all[:, :nocc], atom_of, am_of, elem, labels)

    return IBOResult(
        C_IAO=C_IAO,
        C_IAO_all=C_IAO_all,
        C_AO_all=C_AO_all,
        occupations=occ_all,
        energies=energies_all,
        atom_of=atom_of,
        am_of=am_of,
        func_n=func_n,
        func_dtype=func_dtype,
        elements=list(elem),
        coords=list(coords),
        n_occ=int(nocc),
        partial_charges=np.asarray(net_charges),
        labels=labels,
        analysis_text=msg,
        method=method,
        basis=basis,
        reference=ref,
        mol_name=mol_name,
        mol_spec=mol_spec,
        scf_energy=float(scf_energy),
        wfn=wfn,
    )


def compute_ibo(cjson, options, charge, spin, debug=False):
    """Avogadro adapter: typed core + renderers, preserving the exact
    plugin JSON contract (files under {output-dir}/{name}_NNN/, molden strings,
    message).  See :func:`compute_ibo_data` for the typed surface."""
    import logging
    from .config import load_config as _load_config
    from .config import resolve_output_dir as _resolve_output_dir

    # Determine molecule name and create output directory
    atoms_data = cjson.get("atoms", {})
    elem_raw = atoms_data.get("elements", {}).get("number", [])
    mol_name = cjson.get("name", "") or _mol_formula(elem_raw) or "molecule"
    mol_name = mol_name[:50]
    safe_name = _sanitize_name(mol_name)
    out_root = _resolve_output_dir(options.get("calcs_dir"))
    counter = 1
    while (out_root / f"{safe_name}_{counter:03d}").exists():
        counter += 1
    calc_dir = out_root / f"{safe_name}_{counter:03d}"
    calc_dir.mkdir(parents=True, exist_ok=True)

    _psi_logger = logging.getLogger("psi4")
    _psi_logger.propagate = False
    _psi_logger.setLevel(logging.DEBUG if debug else logging.WARNING)
    _psi_handler = logging.FileHandler(str(calc_dir / "psi4.log"), mode="w")
    _psi_handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    )
    _psi_logger.addHandler(_psi_handler)
    for _name in ["psi4.core", "psi4.driver"]:
        logging.getLogger(_name).setLevel(logging.WARNING)

    import psi4

    psi4.set_output_file(str(calc_dir / "psi4.log"), append=True)

    # -- Merge persistent config into options (adapter concern) ------------
    _cfg = _load_config()
    opts = dict(options)
    for key in ("basis", "method", "iboview_style"):
        if key not in opts and key in _cfg:
            opts[key] = _cfg[key]

    # Write input.xyz BEFORE the SCF so it survives SCF failures.
    atoms = cjson["atoms"]
    coords_raw = atoms["coords"]
    coords = coords_raw["3d"] if isinstance(coords_raw, dict) else coords_raw
    elem = atoms["elements"]["number"]
    _write_input_xyz(calc_dir / "input.xyz", coords, elem, mol_name)

    try:
        res = compute_ibo_data(cjson, opts, charge, spin,
                               psi4_output=calc_dir / "psi4.log")
    except (ValueError, RuntimeError):
        raise
    except Exception as e:
        raise RuntimeError(
            f"{e} Check {calc_dir.name}/psi4.log for details."
        ) from e

    # -- Write Molden with IAO-basis orbitals ------------------------------
    iboview_style = _option(opts, "iboview_style", True)
    molden_path = calc_dir / "ibo.molden"
    if iboview_style:
        sto_mol = psi4.geometry(res.mol_spec)
        psi4.set_options(
            {"basis": "STO-3G", "scf_type": "df", "reference": "rhf", "puream": 0}
        )
        try:
            _, wfn_sto = psi4.energy("hf", return_wfn=True)
        except Exception as e:
            raise RuntimeError(
                "HF/STO-3G (IboView-style rendering) failed. "
                f"Check {calc_dir.name}/psi4.log for details."
            ) from e
        _write_iao_molden(
            molden_path, wfn_sto, res.C_IAO_all, res.occupations,
            res.energies, res.C_IAO_all.shape[1],
        )
        psi4.set_options(
            {
                "basis": res.basis,
                "scf_type": "df",
                "reference": res.reference,
                "e_convergence": 1e-8,
                "d_convergence": 1e-8,
                "puream": 0,
            }
        )
    else:
        _write_iao_molden(
            molden_path, res.wfn, res.C_AO_all, res.occupations,
            res.energies, res.C_IAO_all.shape[1],
        )
    molden_text = molden_path.read_text(encoding="utf-8")

    # -- Canonical Molden (for reference in Avogadro's MO surface dialog) ---
    canon_path = calc_dir / "canonical.molden"
    psi4.molden(res.wfn, str(canon_path))

    # -- Analysis table ------------------------------------------------------
    analysis_path = calc_dir / "ibos.txt"
    analysis_path.write_text(res.analysis_text, encoding="utf-8")

    cjson["atoms"]["partialCharges"] = [
        round(c, 4) for c in res.partial_charges
    ]

    return {
        "readProperties": True,
        "moleculeFormat": "molden",
        "molden": molden_text,
        "cjson": cjson,
        "calcDir": str(calc_dir),
        "message": f"IBO analysis saved to {calc_dir.name}/ibos.txt\n"
        f"Canonical MOs: {calc_dir.name}/canonical.molden",
    }
