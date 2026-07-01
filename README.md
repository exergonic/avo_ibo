# avo_ibo

[![Status: Stable](https://img.shields.io/badge/status-stable-green.svg)](https://github.com/exergonic/avo_ibo)
[![Environment: Pixi](https://img.shields.io/badge/env-pixi-blue.svg)](https://pixi.sh)
[![Python: >=3.11](https://img.shields.io/badge/python-%E2%89%A53.11-green.svg)](https://python.org)
[![License: BSD-3](https://img.shields.io/badge/license-BSD--3-yellow.svg)](LICENSE)
[![Psi4](https://img.shields.io/badge/dependency-Psi4-purple.svg)](https://psicode.org)

Avogadro 2 plugin (and standalone CLI) for computing and visualizing Intrinsic
Bond Orbitals (IBOs).  It uses [Psi4](https://psicode.org) for SCF and wavefunction
generation, then performs post-SCF orbital construction inside the plugin as a
complete implementation of the IAO/2014 algorithm (Knizia, *J. Chem. Theory
Comput.* **2013**, *9*, 4834).  This provides a single workflow for computation,
localization, and interactive orbital inspection in Avogadro.

## 🎯 Why use avo_ibo?

- **Integrated workflow** — compute, localize, and inspect orbitals without switching tools or managing intermediate files.
- **Reproducible pipeline** — one in-process computation with structured logging and output saved to a calculation directory.
- **Complete IAO implementation** — full post-Psi4 orbital construction and localization owned by the plugin, not delegated to external libraries.
- **Multi-purpose outputs** — Molden files for isosurface visualization, analysis tables for assignment and interpretation, and debug logs for troubleshooting.
- **Teaching and research** — display chemically intuitive orbitals that match textbook bonding models (σ, π, lone pairs) for both education and publication.

## ✨ Features

- 🧬 **Complete IAO/2014 implementation in-plugin** — post-Psi4 IAO construction, Pipek-Mezey localization, and virtual IBO generation are implemented directly in avo_ibo, not delegated to Psi4 internals
- ⚡ **Occupied and valence-virtual IBOs** — both localized blocks with Pipek-Mezey (p=4, conv 1e-12, max 2048 iter) for direct chemical interpretation
- 🔬 **On-atom degeneracy resolution** — post-PM Fock diagonalisation recovers physically sensible energy ordering for same-atom orbitals (e.g., O 2s LP vs O 2p LP)
- 🎯 **Valence-virtual via SVD** — SVD-projected virtual IBOs with automatic localization for anti-bond σ\* and π\* isosurfaces
- 📄 **Molden export** — IAO-basis orbitals with Fock-diagonal energies, padded for robust Avogadro rendering
- 📊 **Analysis table** (`calcs/last/ibos.txt`) — occupancy, energy, DOM, bond type, atomic composition, and s/p/d hybridization for assignment and reporting
- 🖥️ **Headless and integrated** — standalone CLI (`python -m avogadro_ibo molecule.xyz`) plus Avogadro plugin mode
- 🚀 **Single in-process pipeline** — no subprocess or file-based IPC; JSON communication with Avogadro

## 📦 What You Get

Each computation produces main outputs in `calcs/last/`:

- **`ibo.molden`** — IAO-basis orbitals ready for visualization in Avogadro or other Molden viewers.
- **`ibos.txt`** — analysis table with occupancy, energy (Ha), degree of monoatomicity (DOM), bond type (σ/π/LP/core), atomic composition, and s/p/d hybridization for each orbital.
- **`canonical.molden`** — canonical MOs for reference (Psi4 output).
- **`psi4.log`** — full Psi4 SCF output for debugging and reproducibility.

## 🔄 Workflow in 3 Steps

1. **Run electronic structure** — Psi4 computes the SCF wavefunction (HF or DFT).
2. **Build and localize IBOs in-plugin** — avo_ibo performs IAO basis construction, Pipek-Mezey localization, on-atom degeneracy resolution, and SVD-based valence-virtual generation.
3. **Inspect and analyze** — Orbitals appear in Avogadro's Molecular Orbitals panel, with analysis outputs written for downstream use.

## 🚀 Quick Start

Requires [pixi](https://pixi.sh). 

If you are using the Avogadro plugin, the pixi bundled with Avogadro is already on your PATH.

```powershell
git clone https://github.com/exergonic/avo_ibo.git
cd avo_ibo
pixi install
pixi run python -m pip install -e .
```
If pixi complains of `no module named pip`, just install it with:
```shell
 pixi run python -m ensurepip
```
Pixi may warn that the lock file uses an older format (v6):
>  WARN the lock file is up-to-date but uses an older format (v6), run `pixi lock` to upgrade to v7 for improved reproducibility

This is safe to ignore. Avogadro's bundled pixi (v0.66.0) reads v6 lock files natively. If using a standalone pixi that upgraded the lock, run `pixi lock --no-update` to keep v6.

### ✅ Tests

```powershell
pixi run test
```

7 CLI integration tests against water, methane, ethene, and ammonia — validates
IAO counts, Molden structure, occupancy partitioning, and per-molecule orbital
patterns.

### 💻 Standalone CLI

```powershell
pixi run python -m avogadro_ibo molecule.xyz
```

Reads an XYZ file, runs the full IBO pipeline, writes to `calcs/last/`
(ibo.molden, canonical.molden, ibos.txt, psi4.log).
Neutral singlet only — for charged or open-shell systems use Avogadro or IboView.

### 🧪 Avogadro 2 plugin

Create a symlink so Avogadro discovers the plugin:

```powershell
New-Item -ItemType SymbolicLink -Path "$env:LOCALAPPDATA\OpenChemistry\Avogadro\plugins\avo_ibo" `
  -Target "C:\path\to\avo_ibo"
```

Then in Avogadro: **Extensions → Intrinsic Bond Orbitals → Compute IBOs**.
The computed IBOs appear in the **Molecular Orbitals** panel where you can
double-click any orbital to toggle its isosurface.  All results are also
saved to `calcs/last/` for reference.

### 🐍 pip (no pixi)

```powershell
pip install git+https://github.com/exergonic/avo_ibo.git
```

> **Note:** Psi4 is not on PyPI — install it via `conda install psi4`.
> This pulls in `numpy` (also required by `avogadro_ibo`) as a transitive dependency.

That's it. Run a calculation with:

```powershell
python -m avogadro_ibo molecule.xyz
```
Output files are stored in in the `.\calcs\last\` directory.

## 📋 Scope and Limitations

**Supported systems:**
- Neutral singlet closed-shell molecules only.
- Tested for small to medium-sized systems (up to ~50 atoms) with standard quantum chemistry basis sets.

**Method and basis support:**
- Hartree-Fock and DFT (B3LYP, PBE, ωB97X, etc.). [not yet implemented]
- Common basis sets: cc-pVDZ, cc-pVTZ, def2-SVP, def2-TZVP, 6-31G\*, etc. [not yet implemented]
- Requires STO-3G availability for minimal-basis IAO construction.

**Known limitations:**
- Open-shell and charged systems are not supported in the current implementation.
- The IAO basis is exact and lossless for the occupied space, which implies zero off-diagonal Fock coupling between occupied and virtual blocks — **no NBO-style E(2) donor-acceptor analysis is possible** in this representation. This is a feature of the IAO design, not a bug. See [mathematics.md §9](mathematics.md#9-why-occ-vir-delocalization-analysis-is-impossible-in-the-iao-basis) for the proof.

## 📚 Further Reading

- **[`mathematics.md`](mathematics.md)** — Full mathematical derivation of the IAO/2014 algorithm, Pipek-Mezey localization, and on-atom degeneracy resolution, with proofs and references to the original literature.
- **[`tutorial.md`](tutorial.md)** — Architecture, algorithm walkthrough, and practical development pitfalls for building Avogadro 2 Pixi plugins.
- **[Troubleshooting](TROUBLESHOOTING.md)** — Common installation issues, runtime errors, and diagnostic steps.
- **[Why no donor-acceptor analysis?](mathematics.md#9-why-occ-vir-delocalization-analysis-is-impossible-in-the-iao-basis)** — IAO exactness makes the occupied-virtual Fock coupling structurally zero by the spectral theorem, preventing E(2)-style delocalization analysis. This is mathematically rigorous, not a numerical artifact.

## 🔗 See Also

- **[IboView](http://www.iboview.org/)** — Gerhard Knizia's reference IBO implementation (standalone tool with its own workflow).
- **[NBO 7.0](https://nbo7.chem.wisc.edu/)** — Natural Bond Orbital analysis for alternative localization schemes and donor-acceptor energetics.
- **[Avogadro 2](https://avogadro.cc)** — Molecular editor and visualization framework.
- **[Psi4](https://psicode.org)** — Open-source quantum chemistry engine.

## ⚖️ License

BSD 3-Clause.  See [LICENSE](LICENSE).
