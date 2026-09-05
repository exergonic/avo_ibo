# How to read an `ibos.txt` file

Every calculation in this directory produces an `ibos.txt` analysis table
with the same six sections. This guide walks through each one using water
(`water.xyz`, wB97X-D/def2-TZVP) as the running example. Orbital pictures
in the text refer to the matching `ibo.molden` file, viewable in Avogadro's
Molecular Orbitals panel.

Run anything yourself with:

```
pixi run python -m avogadro_ibo examples/water.xyz --method wB97X-D --basis def2-TZVP
```

Results land in `calcs/water_NNN/` (numbered per run).

## 1. The orbital table

```
    #      Occ      Energy              Type  Composition              Hybrid                   Ion%        H/L
---------------------------------------------------------------------------------------------------------------
    1    2.000  -19.235944           O(Core)  O1(100.0%)               100% 1s                   ---           
    2    2.000   -0.794692           O-H σ †  O1(62.3%) + H2(37.7%)    23% 2s + 77% 2py         24.7           
    3    2.000   -0.794692           O-H σ †  O1(62.3%) + H3(37.7%)    23% 2s + 77% 2py         24.7           
    4    2.000   -0.637900             O(LP)  O1(100.0%)               55% 2s + 45% 2pz          ---           
    5    2.000   -0.401123             O(LP)  O1(100.0%)               100% 2px                  ---    <- HOMO
    6    0.000    0.436958          O-H σ* †  H2(62.3%) + O1(37.7%)    O: 22% 2s + 78% 2py       ---    <- LUMO
    7    0.000    0.436958          O-H σ* †  H3(62.3%) + O1(37.7%)    O: 22% 2s + 78% 2py       ---           
```

- **`#`** — orbital index. Occupied orbitals come first in the listing
  (sorted by energy within each block); rows 6–7 here are the
  valence-virtual (antibonding) orbitals, `Occ = 0.000`.
- **`Occ`** — 2.000 (doubly occupied) or 0.000 (virtual). The pipeline
  is closed-shell only, so nothing in between ever appears.
- **`Energy`** — the orbital's Fock expectation value in Hartree, not a
  canonical eigenvalue. Ordering is aufbau: core deep, lone pairs and
  bonds in the middle, antibonds positive.
- **`Type`** — the classifier's verdict: `Core`, `(LP)`, `σ`/`π` bonds
  (e.g. `O-H σ`), `anti*`/`π*` virtuals, `2e3c` for three-centre bonds,
  `Deloc` for anything spread over more than two centres. A `†` marks
  membership in a degenerate manifold (see §2).
- **`Composition`** — per-atom populations as percentages: orbital 2 is
  62.3% on O1, 37.7% on H2. Atoms below 0.5% are omitted, so a lone
  pair reads `O1(100.0%)` even though orthogonalization tails exist.
- **`Hybrid`** — s/p/d breakdown on the dominant atom: orbital 4 is
  55% 2s + 45% 2p (an sp-hybridised lone pair), orbital 5 is 100% 2px
  (a pure p lone pair — and the HOMO). For bonds involving hydrogen,
  the heavy atom's hybrid is shown, prefixed when it isn't the top
  atom (`O: 22% 2s + 78% 2py` on row 6).
- **`Ion%`** — percent ionic character from the population asymmetry:
  orbital 2 is 24.7% ionic (62.3/37.7 split). Shown only for two-centre
  bonds (`---` elsewhere).
- **`H/L`** — HOMO/LUMO markers. Selected by occupancy, not position.

Read orbital 4 vs 5 as the classic payoff: two lone pairs on the same
oxygen, split by 0.24 Ha into s-rich (−0.638) and pure-p (−0.401). The
PM functional alone cannot separate them (same atom, same populations);
the post-PM Fock diagonalisation does.

## 2. Degenerate manifolds (†)

```
  † Orbitals 2-3, 6-7 form degenerate manifolds (ΔE < 2e-04 Ha).
  Small energy differences within each manifold are PM convergence noise and do not indicate true energy splittings.
```

Symmetry-equivalent orbitals (the two O–H bonds) converge to energies
that differ by ~1e-5 Ha of optimizer noise. The footnote tells you which
groupings to treat as exactly degenerate. Only trust splittings larger
than the stated threshold as chemistry.

## 3. Frontier Orbital Energies

```
--- Frontier Orbital Energies ---
  HOMO (O(LP), orb 5):  -0.401123 Ha
  LUMO (O-H σ*, orb 6):   0.436958 Ha
  HOMO-LUMO gap:   0.838081 Ha =  22.805 eV =    525.9 kcal/mol
```

The headline numbers: which orbitals bracket the gap, with their table
labels and row numbers for cross-reference, and the gap in three units.
(HF gaps run large — virtuals are unrelaxed. Compare gaps across
molecules and methods, never against experiment directly.)

## 4. Charge Decomposition

```
--- Charge Decomposition ---
   Atom    Z       Pop  Net Charge
----------------------------------
  O1      8     8.494      -0.494
  H2      1     0.753      +0.247
  H3      1     0.753      +0.247
----------------------------------
Total:   10    10.000      +0.000
```

Net atomic charges from summed IBO populations. The `Total` row is the
sanity check: populations must sum to the electron count (10 here) and
the net to the molecular charge (0). If it doesn't, something failed
upstream — do not trust the rest of the file.

## 5. Wiberg Bond Orders (σ/π, density)

```
--- Wiberg Bond Orders (σ/π, density) ---
  W_AB = Σ_{i∈A,j∈B} D²_ij (density Wiberg); σ + π = total exactly.
  σ, π columns include their class's folded interference; the
  parenthetical reports the same interference as (σ-part, π-part).
  Bond         Total       σ       π              (interference)
  O1-H3         0.939   0.939   0.000
  O1-H2         0.939   0.939   0.000
```

Each bond's total order decomposed into σ and π shares that sum exactly
to the total. The parenthetical shows the folded-in inter-orbital
interference split by class — water's rows have none (pure diagonal
shares), so no parenthetical appears. See ethene (`C1-C2: 2.028, σ
1.028, π 1.000 (+0.013: σ+0.013, π+0.000)`) for a row where it does, and
the diborane entry for what it means chemically.

## 6. Significant orbital-pair interference (when present)

Bonds whose orbital-pair interference reaches |term| ≥ 0.01 get a
follow-on section naming the responsible pairs:

```
  B1-H5: orb3(B-B-H 2e3c) × orb4(B-B-H 2e3c): -0.0145
```

Positive pair terms add to the bond order; negative pair terms subtract
from it. Signs are relative to the listed bond — the same pair may
contribute oppositely elsewhere (diborane's bridge×bridge pair is
−0.0145 on every B–H leg but +0.0145 across the B–B contact). Orbital
numbers match the table rows above. Ordinary molecules (water, ethene,
benzene) print no such section — silence means no bond reaches
|term| ≥ 0.01. When the section does print, a closing footnote counts
near-miss terms in [0.005, 0.01) (count plus largest term), so a
rounded parenthetical like (+0.010) always names its pair — see the
water-dimer entry, where O1–O4's +0.0095 earns exactly that footnote.

### How the thresholds are being calibrated

The 0.01 print floor keeps the section to genuine delocalization
chemistry; the 0.005 near-miss band (half the floor) catches terms
that rounding would otherwise orphan. Three gallery cases calibrate
them: alkene C–H σ × C=C π donation prints cleanly at −0.0150/−0.0128
([alkenes.md](alkenes.md)); the water-dimer H-bond contact sits at
+0.0095, inside the footnote band; tert-butyl's finer hyperconjugative
pairs reach only −0.0048, below the band, and are quoted in that
entry's prose instead. The pattern so far: strong donation into a
genuine π* acceptor clears the floor with room to spare, while weak
donation into distant or symmetrized acceptors lands in or under the
gray zone. If future examples contradict that pattern, the floor moves
— the band edge is evidence-based, not sacred.

The full derivation behind all six sections is in
`../mathematics/mathematics.md` (§9 for the analysis table and Wiberg
decomposition).
