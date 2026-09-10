"""Pipek-Mezey localization + degeneracy resolvers (moved verbatim from calcs.py)."""


import numpy as np


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
