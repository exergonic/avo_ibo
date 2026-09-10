"""Pipeline API: IBOResult, compute_ibo_data, compute_ibo adapter (moved verbatim from calcs.py)."""


from dataclasses import dataclass

import numpy as np
import warnings

from .constants import VVO_MIN_SIGMA, _ELEM_SYMBOLS
from .iao import _build_iao_basis, _get_basis_maps
from .localize import _localize_ibos, _resolve_flat_degeneracies, _resolve_on_atom_mixing
from .analysis import _analyze_ibos, _elem_symbol, _format_wiberg
from .molden import _write_iao_molden

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
# Small helpers: options lookup and input-XYZ writer
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
