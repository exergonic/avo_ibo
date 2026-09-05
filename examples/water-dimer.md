# Water dimer — a hydrogen bond, quantified

```
Input: examples/water-dimer.xyz  (O···O ≈ 2.91 Å; ORCA-optimized geometry)
Run:   pixi run python -m avogadro_ibo examples/water-dimer.xyz --method wB97X-D --basis def2-TZVP
```

Donor O1–H2 points at acceptor O4 (H2···O4 ≈ 1.94 Å). The table reads
the textbook hydrogen-bond picture back, with numbers on every arrow:

```
    #      Occ      Energy              Type  Composition                         Hybrid                   Ion%        H/L
    6    2.000   -0.764132             O-H σ  O1(64.2%) + H2(35.8%)               25% 2s + 75% 2px         28.4
    9    2.000   -0.552705             O(LP)  O4(97.0%) + H2(2.2%) + O1(0.8%)     25% 2s + 75% 2px         95.6
   14    0.000    0.494422            O-H σ*  H2(62.0%) + O1(35.0%) + O4(3.0%)    O: 24% 2s + 76% 2px       ---
```

Column guide for newcomers: `Occ` is occupancy, `Ion%` percent ionic
character from the population split, `H/L` HOMO/LUMO markers — all
defined in the [reading guide](how-to-read-ibos.md).

Orbital 9 is the smoking gun: the acceptor lone pair is not on O4
alone — 2.2% of it lives on the *donor's hydrogen*, with a tail
reaching the donor oxygen (0.8%). The back-channel shows too: the
donor σ* (orb 14) carries 3.0% acceptor character. 2.2% of a
lone pair across the contact is covalency you can point at:

![Orbital 9 — acceptor lone pair reaching across the H-bond](img/water-dimer_H-bond.png)
*Orbital 9: the acceptor (O4) lone pair. The major lobe spans the H···O
contact toward the donor; the minor lobe (opposite phase) sits on the
donor oxygen — 97.0% / 2.2% / 0.8%, as tabulated above.*

The donor O–H pays for it. O1–H2 drops to Wiberg 0.867 against
0.932–0.942 for the three free O–H bonds, and polarizes up (64.2/35.8,
ionicity 28.4 vs ~24–26) — the red-shift signature. Meanwhile the
H···O contact itself carries a real bond order, H2–O4 = 0.075, all σ:

```
  Bond         Total       σ       π              (interference)
  O1-H2         0.867   0.867   0.000  (-0.052: σ-0.052, π+0.000)
  H2-O4         0.075   0.075   0.000  (-0.012: σ-0.012, π+0.000)
  O1-O4         0.042   0.042   0.000  (+0.010: σ+0.010, π+0.000)
```

And the interference detail names the mechanism — one orbital pair,
both directions:

```
  O1-H2: orb6(O-H σ) × orb9(O(LP)): -0.0510
  H2-O4: orb6(O-H σ) × orb9(O(LP)): -0.0116
```

LP→σ* donation weakening the donor while building the bridge. The same
pair signs a third bond, too: O1–O4's parenthetical (+0.010) folds from
orb6×orb9 at +0.0095 — just under the detail section's |term| ≥ 0.01
print threshold, so the shipped table shows the number without naming
the pair. One orbital pair, three bonds, three signs (−0.051, −0.012,
+0.010): competition on both legs, cooperation across the O···O
contact. Charges
agree: donor O1 goes more negative (−0.540 vs monomer −0.494),
acceptor O4 less (−0.459) — and the bonded hydrogen H2 (+0.239) is
*less* positive than the free hydrogens (+0.259/+0.260), because
donation into σ* puts density back onto H character. Pure
electrostatics would polarize it the other way.

---
*Geometry optimized in ORCA 6.1.1 (wB97X-D3/def2-TZVP), confirmed
minimum by frequency calculation (no imaginary modes). IBO analysis
with Psi4 1.11 (wB97X-D/def2-TZVP, RHF), avo_ibo 0.4.0; Pipek–Mezey
localization (p = 2 → 4) per Knizia 2013.*
