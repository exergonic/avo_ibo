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
* Valence-virtual IBOs via Single Value Decomposition projection of canonical virtual MOs onto IAO space
* IAO-basis Molden export with Fock-diagonal energies for Avogadro rendering
* Analysis table (`ibos.txt`) — occupancy, energy, bond type, atomic composition, s/p/d hybridization, partial charges, bond orders
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

```shell
git clone https://github.com/exergonic/avo_ibo.git
cd avo_ibo
pixi install
pixi run test
```

The lock file uses v6 format (safe to ignore — the version of pixi
bundled with Avogadro v0.66.0 reads v6 natively).

### Standalone CLI
```shell
pixi run python -m avogadro_ibo molecule.xyz
```

Writes to `calcs/` (ibo.molden, canonical.molden, ibos.txt, psi4.log).

### pip (no pixi)

```shell
pip install git+https://github.com/exergonic/avo_ibo.git
```

Psi4 must be installed separately via conda.

## Data location

Calculations are run in `calcs/` in the current working directory:

* `input.xyz` - the input molecule used for calculations
* `ibos.txt` — analysis table with per-orbital data
* `ibo.molden` — IAO-basis orbitals for visualization
* `canonical.molden` — canonical MOs for reference visualization
* `psi4.log` — Psi4 SCF output

## Limitations and Considerations

* **Closed-shell only.**  
  The IAO/IBO pipeline treats all occupied
  orbitals as doubly occupied (RHF-style).  Open-shell systems
  (radicals, triplet states, broken-symmetry calculations) are not
  supported.  The SCF will still run, but the orbital construction,
  analysis, and Molden output will be invalid.

* **Symmetric molecules.**  
  Pipek-Mezey localization uses fixed sequential Jacobi sweeps.  For
  highly symmetric molecules, symmetry-equivalent orbitals may show
  small (sub-milliHartree) energy splittings (a known consequence of
  the orthogonality constraint — see Knizia JCTC 2013 and
  `mathematics.md`).

* **Analysis Table vs. Isosurface Tails**  
  The IBO analysis table (`ibos.txt`) reports orbital compositions in the
  IAO basis, where populations are clean and bond assignments are crisp.
  The Molden isosurfaces are rendered in the full SCF basis via the
  projection `C_AO = C_IAO @ C_IAO_all`, which correctly includes the
  IAO repolarization components.

  These two representations are slightly inconsistent by construction:
  small density tails visible on non-dominant atoms in the isosurface
  are physically real repolarization contributions, not rendering
  artifacts or bugs.  The analysis table intentionally omits these for
  clarity of chemical interpretation.  This discrepancy is
  mathematically unavoidable and is present in all IAO-based
  implementations.

  These small tails represent the repolarization of each intrinsic atomic orbital in 
  response to the molecular environment. It is the same physics that makes 
  bonds polar and atoms non-spherical in molecules. The analysis table reports
  populations in the compressed IAO basis for chemical clarity; the
  isosurface renders the full physical wavefunction including these
  repolarization contributions.


## License

BSD 3-Clause. See [LICENSE](LICENSE).
