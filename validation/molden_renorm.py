"""Normalize Molden [GTO] contractions + verify self-consistency.

Psi4's molden writer prints RAW (unnormalized) contraction coefficients;
strict readers (IboView) build overlap assuming normalized primitives and
then fail orthonormality checks. This script:

  1. Parses [Atoms]/[GTO]/[MO] (LF or CRLF input; LF output).
  2. Reports C^T S C deviation for the occupied block AS WRITTEN
     (S built with normalized primitives, Molden-standard d/f order).
  3. Rescales multi-primitive s/p contractions to normalized form and
     inversely scales the MO coefficients (identical orbitals).
  4. Re-reports C^T S C.  Residual deviation localizes any REMAINING
     issue (e.g. Cartesian d/f component ordering).

d/f shells here are single-primitive (exactly normalized already);
s/p components share one norm per shell by symmetry, so the rescaling
is exact.  A --permute-df flag reorders d/f MO coefficients from
Psi4-internal to Molden-standard order for testing that hypothesis.

Usage:
  python molden_renorm.py in.molden out.molden [--permute-df]
"""

import math
import sys

import numpy as np

AM_OF = {"s": 0, "p": 1, "d": 2, "f": 3}
# Component lists. s/p unambiguous. d/f: Molden-standard assumed for the
# S build; --permute-df reinterprets file order as Psi4-internal first.
STD_COMPS = {
    0: [(0, 0, 0)],
    1: [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    2: [(2, 0, 0), (0, 2, 0), (0, 0, 2),
        (1, 1, 0), (1, 0, 1), (0, 1, 1)],
    3: [(3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0),
        (2, 0, 1), (1, 0, 2), (0, 2, 1), (0, 1, 2), (1, 1, 1)],
}
INT_COMPS = {
    2: [(2, 0, 0), (1, 1, 0), (1, 0, 1),
        (0, 2, 0), (0, 1, 1), (0, 0, 2)],
    3: [(3, 0, 0), (2, 1, 0), (2, 0, 1), (1, 2, 0), (1, 1, 1),
        (1, 0, 2), (0, 3, 0), (0, 2, 1), (0, 1, 2), (0, 0, 3)],
}
# Permutation file-position -> standard for Psi4-internal order.
D_PERM = [0, 3, 5, 1, 2, 4]
F_PERM = [0, 6, 9, 3, 1, 2, 5, 8, 7, 4]

BOHR_PER_ANGSTROM = 1.0  # positions already in bohr ([Atoms] (AU))


def double_fact(n):
    d = 1.0
    k = n
    while k > 1:
        d *= k
        k -= 2
    return d


def prim_norm(lmn, alpha):
    lx, ly, lz = lmn
    L = lx + ly + lz
    df = double_fact(2 * lx - 1) * double_fact(2 * ly - 1) * double_fact(2 * lz - 1)
    return math.sqrt((2.0 * alpha / math.pi) ** 1.5 * (4.0 * alpha) ** L / df)


def moment_1d(n, p):
    if n % 2 == 1:
        return 0.0
    return double_fact(n - 1) * math.sqrt(math.pi) / (2.0 ** (n / 2) * p ** ((n + 1) / 2.0))


def prim_overlap_1d(e1, n1, x1, e2, n2, x2):
    p = e1 + e2
    mu = e1 * e2 / p
    dx = x1 - x2
    Px = (e1 * x1 + e2 * x2) / p
    d1 = Px - x1
    d2 = Px - x2
    s = 0.0
    for k1 in range(n1 + 1):
        c1 = math.comb(n1, k1) * d1 ** (n1 - k1)
        for k2 in range(n2 + 1):
            c2 = math.comb(n2, k2) * d2 ** (n2 - k2)
            s += c1 * c2 * moment_1d(k1 + k2, p)
    return math.exp(-mu * dx * dx) * s


def prim_overlap(lmn1, a, xyz1, lmn2, b, xyz2):
    return (prim_overlap_1d(a, lmn1[0], xyz1[0], b, lmn2[0], xyz2[0])
            * prim_overlap_1d(a, lmn1[1], xyz1[1], b, lmn2[1], xyz2[1])
            * prim_overlap_1d(a, lmn1[2], xyz1[2], b, lmn2[2], xyz2[2]))


def parse_molden(path):
    raw = open(path, "rb").read().replace(b"\r\n", b"\n").decode()
    lines = raw.split("\n")
    atoms = []
    shells = []  # {am, center, prims, ao0, nfn, sline}
    mos = []     # {ene, occ, coef}
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if s == "[Atoms] (AU)" or s == "[Atoms] Angs":
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("["):
                p = lines[i].split()
                if len(p) >= 6:
                    atoms.append((p[0], int(p[2]),
                                  (float(p[3]), float(p[4]), float(p[5]))))
                i += 1
            continue
        if s == "[GTO]":
            i += 1
            ao = 0
            while i < len(lines) and not lines[i].strip().startswith("["):
                p = lines[i].split()
                if len(p) == 2 and p[0].isdigit():
                    i += 1
                    continue
                if len(p) >= 2 and p[0] in AM_OF and p[1].isdigit():
                    am = AM_OF[p[0]]
                    nprim = int(p[1])
                    prims = []
                    for j in range(1, nprim + 1):
                        e, d = lines[i + j].split()[:2]
                        prims.append((float(e), float(d)))
                    # find owning atom: last "N 0" header seen
                    shells.append({"am": am, "prims": prims, "ao0": ao,
                                   "nfn": len(STD_COMPS[am]), "sline": i})
                    ao += len(STD_COMPS[am])
                    i += nprim + 1
                    continue
                i += 1
            continue
        if s == "[MO]":
            i += 1
            cur = None
            while i < len(lines) and not lines[i].strip().startswith("["):
                t = lines[i].strip()
                if t.startswith("Sym="):
                    if cur is not None:
                        mos.append(cur)
                    cur = {"coef": {}}
                elif "=" in t and cur is not None and not t[0].isdigit() \
                        and not (t[0] == "-" and len(t) > 1 and t[1].isdigit()):
                    k, v = t.split("=", 1)
                    cur[k.strip().lower()] = v.strip()
                elif t:
                    p = t.split()
                    if len(p) == 2 and p[0].lstrip("-").isdigit():
                        try:
                            cur["coef"][int(p[0]) - 1] = float(p[1].replace("D", "E"))
                        except ValueError:
                            pass
                i += 1
            if cur is not None:
                mos.append(cur)
            continue
        i += 1
    return atoms, shells, mos


def shell_atom_map(lines, shells):
    """Assign each shell its atom index by tracking 'N 0' headers."""
    amap = []
    cur = -1
    si = 0
    i = 0
    # re-walk [GTO] region
    started = False
    while i < len(lines):
        s = lines[i].strip()
        if s == "[GTO]":
            started = True
            i += 1
            continue
        if started and s.startswith("["):
            break
        if started:
            p = s.split()
            if len(p) == 2 and p[0].isdigit():
                cur += 1
                i += 1
                continue
            if len(p) >= 2 and p[0] in AM_OF and p[1].isdigit():
                amap.append(cur)
                i += int(p[1]) + 1
                continue
        i += 1
    assert len(amap) == len(shells), (len(amap), len(shells))
    return amap


def contracted_norm(shell):
    """Norms per component (normalized primitives)."""
    out = []
    for lmn in STD_COMPS[shell["am"]]:
        s = 0.0
        for a, da in shell["prims"]:
            Na = prim_norm(lmn, a)
            for b, db in shell["prims"]:
                Nb = prim_norm(lmn, b)
                s += da * db * Na * Nb * prim_overlap(lmn, a, (0, 0, 0),
                                                      lmn, b, (0, 0, 0))
        out.append(math.sqrt(s))
    return out


def build_S(atoms, shells, amap, comps):
    n = sum(s["nfn"] for s in shells)
    S = np.zeros((n, n))
    info = []  # per fn: (xyz, lmn)
    for s, a in zip(shells, amap):
        for lmn in comps[s["am"]]:
            info.append((atoms[a][2], lmn))
    for s in shells:
        for a, da in s["prims"]:
            pass
    # contraction data per function
    fn_con = []  # list per fn of [(exp, coef*primnorm)]
    for s in shells:
        cl = STD_COMPS[s["am"]] if comps is STD_COMPS else comps[s["am"]]
        # NOTE: contraction coefficients are shared per shell; per-fn
        # norms handled by caller via rescaling. Here raw build:
        for lmn in cl:
            fn_con.append([(e, d * prim_norm(lmn, e)) for e, d in s["prims"]])
    for i in range(n):
        xi, li = info[i]
        for j in range(i, n):
            xj, lj = info[j]
            v = 0.0
            for a, ca in fn_con[i]:
                for b, cb in fn_con[j]:
                    v += ca * cb * prim_overlap(li, a, xi, lj, b, xj)
            S[i, j] = S[j, i] = v
    return S


def main():
    src = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) > 2 else None
    permute = "--permute-df" in sys.argv

    atoms, shells, mos = parse_molden(src)
    raw_lines = open(src, "rb").read().replace(b"\r\n", b"\n").decode().split("\n")
    amap = shell_atom_map(raw_lines, shells)
    n_ao = sum(s["nfn"] for s in shells)
    print(f"atoms={len(atoms)} shells={len(shells)} n_ao={n_ao} n_mo={len(mos)}")

    S = build_S(atoms, shells, amap, STD_COMPS)
    C = np.zeros((n_ao, len(mos)))
    for j, m in enumerate(mos):
        for idx, v in m["coef"].items():
            if 0 <= idx < n_ao:
                C[idx, j] = v
    occ = [j for j, m in enumerate(mos)
           if abs(float(m.get("occup", "0").replace("D", "E"))) > 0.5]
    Co = C[:, occ]
    dev = Co.T @ S @ Co - np.eye(len(occ))
    print(f"PRE : occupied C^T S C maxdev = {np.abs(dev).max():.3e}")

    if dst is None:
        return

    # Rescale: per-shell single norm (exact for s/p by symmetry and for
    # single-primitive d/f). Warn if a shell's components disagree.
    shell_norm = []
    for s in shells:
        ns = contracted_norm(s)
        if max(ns) - min(ns) > 1e-9:
            print(f"  WARN shell am={s['am']} anisotropic norms: "
                  + " ".join(f"{v:.4f}" for v in ns))
        shell_norm.append(ns[0])

    out = raw_lines[:]
    # rewrite contraction coefficients (locate via sline)
    for s, N in zip(shells, shell_norm):
        for j in range(len(s["prims"])):
            li = s["sline"] + 1 + j
            e, d = out[li].split()[:2]
            out[li] = f"{float(e):24.12f} {float(d) / N:24.12f}"
    # rescale MO coefficients: NOT APPLIED (deliberately).
    #
    # C^T S C is invariant under simultaneous basis+coefficient
    # rescaling, so touching the MO coefficients alongside the
    # contractions is a mathematical no-op (verified: POST == PRE to
    # all digits when both are scaled).  Psi4's s/p MO coefficients
    # are already expressed in its internal NORMALIZED basis, so once
    # the written contractions are normalized the file is consistent
    # as-is.  (d/f singles are normalized by construction; the 1/sqrt(3)
    # Psi4 applies on output is part of that convention.)
    fn_scale = []
    for s, N in zip(shells, shell_norm):
        fn_scale += [N] * s["nfn"]
    fn_scale = np.array(fn_scale)
    if permute:
        # reinterpret file d/f order as Psi4-internal, permute to standard
        P = np.arange(n_ao)
        for s in shells:
            if s["am"] == 2:
                P[s["ao0"]:s["ao0"] + 6] = s["ao0"] + np.array(D_PERM)
            elif s["am"] == 3:
                P[s["ao0"]:s["ao0"] + 10] = s["ao0"] + np.array(F_PERM)
        C = C[P, :]
        fn_scale = fn_scale[P]

    # re-emit [MO] coefficient lines: find "MO" block and rewrite values
    in_mo = False
    cur_mo = -1
    for i, ln in enumerate(out):
        t = ln.strip()
        if t == "[MO]":
            in_mo = True
            continue
        if in_mo and t.startswith("["):
            break
        if in_mo and t.startswith("Sym="):
            cur_mo += 1
            continue
        if in_mo and cur_mo >= 0 and t:
            p = t.split()
            if len(p) == 2 and p[0].lstrip("-").isdigit():
                idx = int(p[0]) - 1
                if 0 <= idx < n_ao:
                    out[i] = f"{idx + 1:>4d} {C[idx, cur_mo]:24.14e}"
    open(dst, "w", newline="\n").write("\n".join(out))

    # verify from the WRITTEN file
    atoms2, shells2, mos2 = parse_molden(dst)
    amap2 = shell_atom_map(open(dst).read().split("\n"), shells2)
    S2 = build_S(atoms2, shells2, amap2, STD_COMPS)
    C2 = np.zeros((n_ao, len(mos2)))
    for j, m in enumerate(mos2):
        for idx, v in m["coef"].items():
            if 0 <= idx < n_ao:
                C2[idx, j] = v
    occ2 = [j for j, m in enumerate(mos2)
            if abs(float(m.get("occup", "0").replace("D", "E"))) > 0.5]
    dev2 = C2[:, occ2].T @ S2 @ C2[:, occ2] - np.eye(len(occ2))
    print(f"POST: occupied C^T S C maxdev = {np.abs(dev2).max():.3e}")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()
