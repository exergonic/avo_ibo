# avogadro-ibo

Avogadro 2 plugin for computing and visualizing Intrinsic Bond Orbitals (IBOs) using [pySCF](https://pyscf.org).

IBOs (Knizia, *J. Chem. Theory Comput.* **2013**, *9*, 4834) provide a chemically intuitive picture of bonding by localizing occupied molecular orbitals onto atoms — revealing σ-bonds, π-bonds, lone pairs, and core orbitals directly from first-principles DFT.

## How it works

1. Avogadro sends the current molecule geometry to the plugin
2. pySCF runs a DFT calculation (PBE/def2-SVP)
3. Intrinsic Atomic Orbitals (IAOs) are constructed
4. Occupied MOs are localized to IBOs using Knizia's algorithm
5. IBO coefficients are written to a Molden file and returned to Avogadro
6. Double-click any IBO in the Data Sets panel to view its isosurface

## Installation

```bash
cd avo_ibo
pixi run pip install -e .
```

Then the plugin will appear under **Extensions → Intrinsic Bond Orbitals → Compute IBOs** in Avogadro 2.
