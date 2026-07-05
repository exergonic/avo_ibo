# avo_ibo

`avo_ibo` is a plugin for Avogadro 2 (and standalone CLI) that computes and
visualizes Intrinsic Bond Orbitals (IBOs) using
[Psi4](https://psicode.org) for SCF and wavefunction generation, then
performs post-SCF IAO/2014 orbital construction (Knizia, *J. Chem. Theory
Comput.* **2013**, *9*, 4834) and Pipek-Mezey localization directly in the
plugin.

## Capabilities

* Occupied and valence-virtual IBOs — Pipek-Mezey localization (p=2 warmup + p=4 refinement, conv 1e-12)
* On-atom degeneracy resolution via post-PM Fock diagonalization
* Valence-virtual IBOs via SVD projection of canonical virtual MOs onto IAO space
* IAO-basis Molden export with Fock-diagonal energies for Avogadro rendering
* Analysis table (`calcs/last/ibos.txt`) — occupancy, energy, DOM, bond type, atomic composition, s/p/d hybridization
* Standalone CLI (`python -m avogadro_ibo molecule.xyz`) and Avogadro in-app mode


## Quick Start

### Avogadro Plugin (easiest)
Download this repo from GitHub (Code → Download ZIP).  Unzip it, then in
Avogadro click **Extensions → Manage Plugins... → Install from Directory...**
and choose the extracted `avo_ibo` folder.

**Extensions → Intrinsic Bond Orbitals → Compute IBOs**.
Orbitals appear in the **Molecular Orbitals** panel.

### Development Setup
Requires [pixi](https://pixi.sh).

```powershell
git clone https://github.com/exergonic/avo_ibo.git
cd avo_ibo
pixi install
pixi run test
```

The lock file uses v6 format (safe to ignore — the version of pixi
bundled with Avogadro v0.66.0 reads v6 natively).

### Standalone CLI
```powershell
pixi run python -m avogadro_ibo molecule.xyz
```

Writes to `calcs/last/` (ibo.molden, canonical.molden, ibos.txt, psi4.log).

### pip (no pixi)

```powershell
pip install git+https://github.com/exergonic/avo_ibo.git
```

Psi4 must be installed separately via conda.

## Data location

Calculations are run in `calcs/last/` in the current working directory:

* `ibo.molden` — IAO-basis orbitals for visualization
* `ibos.txt` — analysis table with per-orbital data
* `canonical.molden` — canonical MOs for reference
* `psi4.log` — Psi4 SCF output

## License

BSD 3-Clause. See [LICENSE](LICENSE).
