"""IAO/2014 construction."""


import numpy as np


def get_basis_maps(basis):
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


def build_iao_basis(S, S12, S_min, C_occ):
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
