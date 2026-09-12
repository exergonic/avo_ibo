# Diborane — three-centre two-electron bridges

```
Input: diborane.xyz  (wB97X-D/6-31G(d,p) geometry)
Run:   pixi run python -m avogadro_ibo examples/diborane.xyz --method wB97X-D --basis def2-TZVP
```

The two bridging hydrogens refuse the two-centre picture. That is
the point of the `2e3c` label: the classifier was not told "this is
diborane." It sees a third atom above 10% and a fourth below 3%,
and writes `B-B-H 2e3c`. Hydrogen carries the largest share, the
two borons split the rest nearly equally, and the pair is nearly
degenerate — against four ordinary terminal B–H σ bonds (0.984
each):

```
    #      Occ      Energy              Type  Composition              Hybrid                   Ion%        H/L
    3    2.000   -0.610304      B-B-H 2e3c †  H5(45.1%) + B1(27.3%) + B6(27.3%)               B: 22% 2s + 78% 2pz      24.6           
    4    2.000   -0.610275      B-B-H 2e3c †  H3(45.1%) + B1(27.3%) + B6(27.3%)               B: 22% 2s + 78% 2pz      24.7           
```

A 2e3c bond is not a full two-centre bond to anyone. Each bridge
leg is about half an order (0.482); the B–B contact (0.634) is the
two bridges reinforcing each other across the borons. The detail
section attributes every one of those parentheticals to a *single*
orbital pair, both directions:

```
  B1-B6: orb3(B-B-H 2e3c) × orb4(B-B-H 2e3c): +0.0134
  B1-H5: orb3(B-B-H 2e3c) × orb4(B-B-H 2e3c): -0.0134
  ...
  H3-H5: orb3(B-B-H 2e3c) × orb4(B-B-H 2e3c): +0.0134
```

Negative on a shared B–H leg means the two bridges are
orthogonalizing against each other — they compete for the same
atoms. Positive on B–B (and on the H···H contact between the
bridges) means the same pair cooperates across the contact. One
orbital pair, three signs, each chemically legible. When a later
page (alkenes, the water dimer) shows the same compete/cooperate
pattern, it is this reading: the aggregate column only says "minus
here, plus there"; the detail says which orbitals, and why.

![Diborane 3c-2e bridge](img/diborane_2e3c_bond.png)
