"""Unit tests for orbital type classification (σ / π / δ).

Synthetic two-center Fe–C IAO coefficient vectors: six Cartesian d
functions on Fe (Psi4 order xx, xy, xz, yy, yz, zz) plus one s on C.
No SCF — these assert the Type-column rule in isolation.
"""

import numpy as np

from avogadro_ibo.analysis import _classify_orbital, _spd_frac


# Fe: 6 Cartesian d; C: 1 s.
_ATOM_OF = np.array([0, 0, 0, 0, 0, 0, 1])
_AM_OF = np.array([2, 2, 2, 2, 2, 2, 0])
_FUNC_N = np.array([3, 3, 3, 3, 3, 3, 2])
_ELEM = np.array([26, 6])  # Fe, C
_XX, _XY, _XZ, _YY, _YZ, _ZZ, _C_S = range(7)


def _two_center(d_indices, fe_pop=0.84, c_pop=0.16):
    """Unit-normalized Fe(d) + C(s) coefficient vector and populations."""
    c = np.zeros(7, dtype=np.float64)
    w = np.sqrt(fe_pop / len(d_indices))
    for i in d_indices:
        c[i] = w
    c[_C_S] = np.sqrt(c_pop)
    pop = np.array([fe_pop, c_pop])
    order = np.argsort(-pop)
    s_char, p_char, d_char = _spd_frac(c, _AM_OF, order[0], _ATOM_OF)
    return dict(
        oc=2.0,
        pop=pop,
        order=order,
        top_A=int(order[0]),
        top_B=int(order[1]),
        s_char=s_char,
        p_char=p_char,
        d_char=d_char,
        elem=_ELEM,
        am_of=_AM_OF,
        atom_of=_ATOM_OF,
        func_n=_FUNC_N,
        c=c,
    )


def _label(**kwargs):
    return _classify_orbital(**kwargs)


def test_dxy_metal_ligand_is_delta():
    """e₂′ (dxy) two-center overlap is δ, not the σ default."""
    assert _label(**_two_center([_XY])) == "Fe-C δ"


def test_dx2y2_metal_ligand_is_delta():
    """e₂′ (dx²−y²) two-center overlap is δ."""
    # (xx − yy); equal opposite signs → dx2y2 spherical weight wins.
    kw = _two_center([_XX, _YY])
    kw["c"][_XX] = np.sqrt(0.42)
    kw["c"][_YY] = -np.sqrt(0.42)
    assert _label(**kw) == "Fe-C δ"


def test_dz2_metal_ligand_stays_sigma():
    """a₁′ (dz²) along z is σ — the exclusion that keeps ligand-field σ intact."""
    assert _label(**_two_center([_ZZ])) == "Fe-C σ"


def test_dxz_metal_ligand_is_pi():
    """e₁′ (dxz) two-center overlap is π even though the metal side is pure d."""
    assert _label(**_two_center([_XZ])) == "Fe-C π"


def test_pure_p_two_center_is_pi():
    """Organic π (high p-fraction both atoms) is unchanged."""
    # Two atoms, three p each: atom 0 px+py+pz, atom 1 pz.
    atom_of = np.array([0, 0, 0, 1, 1, 1])
    am_of = np.array([1, 1, 1, 1, 1, 1])
    func_n = np.array([2, 2, 2, 2, 2, 2])
    c = np.array([0.0, 0.0, np.sqrt(0.5), 0.0, 0.0, np.sqrt(0.5)])
    pop = np.array([0.5, 0.5])
    order = np.array([0, 1])
    s_char, p_char, d_char = _spd_frac(c, am_of, 0, atom_of)
    assert (
        _classify_orbital(
            2.0, pop, order, 0, 1, s_char, p_char, d_char,
            np.array([6, 6]), am_of, atom_of, func_n, c,
        )
        == "C-C π"
    )


def test_metal_lp_stays_lp():
    """High-DOM metal (ferrocene dz² nonbonding) never reaches the two-center branch."""
    kw = _two_center([_ZZ], fe_pop=0.986, c_pop=0.014)
    assert _label(**kw) == "Fe(LP)"


def _fe_c_system():
    """One-orbital Fe(dxy)+C(s) system for the Wiberg folding test."""
    from avogadro_ibo.analysis import format_wiberg

    C = np.zeros((7, 1), dtype=np.float64)
    C[1, 0] = np.sqrt(0.84)  # Fe dxy
    C[6, 0] = np.sqrt(0.16)  # C s
    atom_of = np.array([0, 0, 0, 0, 0, 0, 1])
    am_of = np.array([2, 2, 2, 2, 2, 2, 0])
    return format_wiberg(C, atom_of, am_of, np.array([26, 6]))


def test_wiberg_has_delta_column():
    """Wiberg table carries σ + π + δ exactly; δ-labeled shares fold as δ."""
    out = _fe_c_system()
    header = next(ln for ln in out.splitlines() if "Bond" in ln and "Total" in ln)
    assert "δ" in header
    row = next(ln for ln in out.splitlines() if ln.strip().startswith("Fe1-C"))
    nums = [float(t) for t in row.split()[1:5]]
    total, s, p, d = nums
    assert d > 0.0
    assert abs(total - (s + p + d)) < 1e-9
