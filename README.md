# avo_ibo

[![Status: Stable](https://img.shields.io/badge/status-stable-green.svg)](https://github.com/exergonic/avo_ibo)
[![Environment: Pixi](https://img.shields.io/badge/env-pixi-blue.svg)](https://pixi.sh)
[![Python: >=3.11](https://img.shields.io/badge/python-%E2%89%A53.11-green.svg)](https://python.org)
[![License: BSD-3](https://img.shields.io/badge/license-BSD--3-yellow.svg)](LICENSE)
[![Psi4](https://img.shields.io/badge/dependency-Psi4-purple.svg)](https://psicode.org)

Avogadro 2 plugin (and standalone CLI) for computing and visualizing Intrinsic
Bond Orbitals (IBOs) via [Psi4](https://psicode.org).  IBOs (Knizia, *J. Chem.
Theory Comput.* **2013**, *9*, 4834) provide a chemically intuitive picture of
bonding — σ-bonds, π-bonds, lone pairs, and core orbitals — from first-principles
HF/DFT.

## ✨ Features

- 🧬 **Full IAO/2014 pipeline** — depolarized + repolarized IAO construction (Knizia Appendix C), matching IboView's `MakeIaoBasisNew`
- ⚡ **Pipek-Mezey localization** (p=4, conv 1e-12, max 2048 iter) for both occupied and valence-virtual blocks
- 🔬 **On-atom degeneracy resolution** — post-PM Fock diagonalisation separates same-atom orbitals (e.g. O 2s LP vs O 2p LP) by energy
- 🎯 **Valence-virtual IAOs** via SVD (`C_IAO^T @ S @ C_vir`), localizing anti-bond σ\* and π\* orbitals for isosurface viewing
- 📄 **Molden output** — IAO-basis orbitals with Fock-diagonal energies, zero-energy padded to `n_AO` slots for robust Avogadro rendering
- 📊 **Analysis table** (`calcs/last/ibos.txt`) — occupancy, energy, DOM, bond type, atomic composition, and s/p/d hybridization for every IBO
- 🖥️ **Standalone CLI** — `python -m avogadro_ibo molecule.xyz` for headless IBO computation
- 🚀 **Psi4 in-process** — single Python process, no subprocess or file-based IPC

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

The Avogadro plugin symlink still works the same way; Avogadro calls the
`avogadro-ibo` entry point under the hood.

### 💻 Standalone CLI

```powershell
pixi run python -m avogadro_ibo molecule.xyz
```

Reads an XYZ file, runs the full IBO pipeline, writes to `calcs/last/`
(ibo.molden, canonical.molden, ibos.txt, psi4.log).
Neutral singlet only — for charged or open-shell systems use Avogadro or IboView.

### ✅ Tests

```powershell
pixi run test
```

7 CLI integration tests against water, methane, ethene, and ammonia — validates
IAO counts, Molden structure, occupancy partitioning, and per-molecule orbital
patterns.

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


## 📚 Further Reading

- [`mathematics.md`](mathematics.md) — full derivation of the IAO/IBO pipeline
- [`tutorial.md`](tutorial.md) — architecture, algorithm walkthrough, and gotchas
- [Why no donor-acceptor analysis?](https://github.com/exergonic/avo_ibo/blob/master/mathematics.md#9-why-occ-vir-delocalization-analysis-is-impossible-in-the-iao-basis) — IAO exactness makes F_ov = 0 by spectral theorem, killing NBO-style E(2) analysis

## ⚖️ License

BSD 3-Clause.  See [LICENSE](LICENSE).
