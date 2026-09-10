"""IAO-basis Molden writer."""


import numpy as np


def write_iao_molden(path, wfn, C_AO, occ, energies, n_orb):
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
