# avogadro-ibo

Avogadro 2 plugin (and standalone CLI) for computing and visualizing Intrinsic
Bond Orbitals (IBOs) using [Psi4](https://psicode.org).

IBOs (Knizia, *J. Chem. Theory Comput.* **2013**, *9*, 4834) provide a chemically
intuitive picture of bonding by localizing occupied molecular orbitals onto
atoms — revealing σ-bonds, π-bonds, lone pairs, and core orbitals directly
from first-principles HF/DFT.  See [`MATHEMATICS.md`](MATHEMATICS.md) for the
full mathematical derivation.

## How it works

1. **Avogadro** sends the current molecule geometry (CJSON) to the plugin.
2. **Psi4** runs a HF/DFT calculation (default: HF/cc-pVDZ, Cartesian
   functions, `puream=0`).  The SCF basis can be overridden via options.
3. **Intrinsic Atomic Orbitals (IAOs)** are constructed using the full
   IAO/2014 algorithm (Knizia Appendix C, matching IboView's
   `MakeIaoBasisNew`).  A minimal basis (STO-3G) defines the split between
   occupied atomic orbitals and virtual residual.
4. **Occupied IBOs** are obtained via Pipek-Mezey localization (p=4, eq 4)
   in the IAO basis (Jacobi sweeps, conv 1e-12, max 2048 iterations).
   A warm-start strategy (p=2 first, then p=4 refine) avoids local minima.
5. **On-atom degeneracy resolution**: The PM functional uses only atomic
   populations n_A(i), so orbitals on the same atom with DOM ≈ 1 (e.g.
   O 2s and O lone pair) are degenerate — any rotation within that
   subspace yields the same L value.  After PM, a post-processing step
   diagonalises the Fock matrix within each same-atom, high-DOM subspace,
   restoring the aufbau ordering (s-rich lowest, p-rich highest).
6. **Valence-virtual IAOs** are constructed via SVD of `C_IAO^T @ S @ C_vir`
   (IboView `MakeValenceVirtuals`), keeping singular values > 1e-8, then
   localized with the same PM p=4 procedure.  This gives clean anti-bond
   σ\* and π\* orbitals for visualisation.
7. Orbital energies are computed as ε_i = C_i^T @ F_AO @ C_i (Fock-matrix
   diagonal elements in the AO basis), matching IboView's
   `MakeOrbitalEnergies_General`.
8. All n_min IAO-basis orbitals (occupied + valence-virtual) are written
   to a **Molden file** with custom [MO] sections, and returned to Avogadro.
9. An **analysis table** (`calcs/last/ibos.txt`) lists each orbital's
   occupancy, energy, DOM, assigned type (σ/π/LP/Core/anti\*), atomic
   composition, and s/p/d hybridization.
10. Double-click any IBO in Avogadro's **Data Sets** panel to view its
    isosurface.

### Algorithm pipeline (inside `compute_ibo`)

```
CJSON geom → Psi4 HF/cc-pVDZ → C_occ (MO coefficients)
    ↓
C_IAO, C_IAO_occ = _build_iao_basis(S_full, S12, S_min, C_occ)
    ↓
_localize_ibos(C_IAO_occ)       ← PM p=2→p=4 warm-start
    ↓
F_IAO = C_IAO^T @ F_AO @ C_IAO
    ↓
_resolve_on_atom_mixing(C_IAO_occ)  ← Fock diag per same-atom group
    ↓
occ_energies[i] = C_i^T @ F_IAO @ C_i
    ↓
Valence-virtual: SVD(SIbVir) → U_val → _localize_ibos(U_val)
    ↓
Combine, sort by energy → Molden output → analysis table
```

## Installation

All paths require [pixi](https://pixi.sh).  From the repo root:

```powershell
# Step 1: install dependencies (psi4, numpy) via pixi
pixi install

# Step 2: install the package as editable (creates entry points)
pixi run pip install -e .
```

### Usage as standalone CLI

```powershell
pixi run python -m avogadro_ibo molecule.xyz
```

Reads an XYZ file, runs the full IBO pipeline, and writes results to
`calcs/last/` (ibo.molden, canonical.molden, ibos.txt, psi4.log).
Neutral singlet only — for charged or open-shell systems use
Avogadro or IboView.

### Usage as Avogadro 2 plugin

The plugin will appear under **Extensions → Intrinsic Bond Orbitals →
Compute IBOs** in Avogadro 2, provided that Avogadro can discover it.

Avogadro 2 scans `%LOCALAPPDATA%\OpenChemistry\Avogadro\plugins\` for
Python script plugins.  Create a symlink:

```powershell
New-Item -ItemType SymbolicLink -Path "$env:LOCALAPPDATA\OpenChemistry\Avogadro\plugins\avo_ibo" `
  -Target "C:\path\to\avo_ibo"
```

1. Open a molecule in Avogadro 2.
2. **Extensions → Intrinsic Bond Orbitals → Compute IBOs**.
3. Wait for the calculation to complete (Psi4 runs in-process).
4. Double-click any orbital in **Data Sets** to toggle its isosurface.
5. The analysis table is available in `calcs/last/ibos.txt`.
6. **Extensions → Intrinsic Bond Orbitals → Go to files …** opens the
   `calcs/` directory in Explorer.

### Tests

```powershell
pixi run test
```

Runs a CLI integration suite against `tests/files/` (water, methane, ethene, ammonia). Validates: correct number of IAOs, Molden entry count matches n_AO, occupancy partitioning (occ/vir), on-atom degeneracy resolution (water), C-H sigma degeneracy (methane), and lone pair detection (ammonia). Tests call the CLI as a black box — no production code modified for testability.

### Options

The `options` dict passed from Avogadro (or via stdin) can contain:

| Key      | Default   | Description |
|----------|-----------|-------------|
| `basis`  | `cc-pVDZ` | SCF basis set name (Psi4 format) |
| `method` | `hf`      | Quantum chemistry method (hf, b3lyp, etc.) |
