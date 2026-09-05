# Allene — orthogonal π systems and axial chirality

```
Input: allene.xyz  (D2d)
Run:   pixi run python -m avogadro_ibo examples/allene.xyz --method wB97X-D --basis def2-TZVP
```

The two π bonds resolve by symmetry into perpendicular planes:

```
   10    2.000   -0.370755           C-C π †  C3(53.9%) + C2(45.0%)               100% 2pz                  9.1           
   11    2.000   -0.370754           C-C π †  C1(53.9%) + C2(45.0%)               100% 2py                  9.1    <- HOMO
   12    0.000    0.196556          C-C π* †  C2(52.8%) + C1(44.4%) + C3(2.1%)    100% 2py                  ---    <- LUMO
   13    0.000    0.196557          C-C π* †  C2(52.8%) + C3(44.4%) + C1(2.1%)    100% 2pz                  ---           
```

100% 2pz vs 100% 2py, degenerate to 1e-6 Ha, π* mirroring π — the two
ends of the molecule are electronically isolated, which is the quantum
basis of allene's axial chirality. C–C orders are 1.989 each: two
independent, unconjugated double bonds. The central carbon sits nearly
neutral (+0.026) at the intersection of the two π systems while the
terminals carry −0.189 from C–H donation.

Notably, allene's detail section fires genuine σπ mixing lines —
`C2-C3: orb6(C-H σ) × orb10(C-C π): -0.0247` alongside `+0.0123`
partners — the bent-bond early warning working as designed in a
molecule where σ and π manifolds legitimately interpenetrate off-axis.
Compare ethene, whose σπ channel reads exactly zero: planarity protects
it there; D2d exposes it here.
