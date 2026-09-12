# Formaldehyde — the lone pair that reaches into the C–H bonds

Formaldehyde's HOMO is drawn in every textbook as an oxygen in-plane
lone pair. The IBO table agrees — and then quantifies the fine print
the textbooks omit. That lone pair leaks onto both C–H bonds. The
pair term, −0.0625, is the largest hyperconjugative fingerprint in
this gallery.

```
Input: formaldehyde.xyz  (ORCA 6.1.1 wB97X-D3/def2-TZVP opt + freq, no imaginary modes)
Run:   pixi run python -m avogadro_ibo examples/formaldehyde.xyz --method wB97X-D --basis def2-TZVP
```

```
    #      Occ      Energy              Type  Composition              Hybrid                   Ion%        H/L
    8    2.000   -0.411669             O(LP)  O2(93.5%) + C1(3.1%) + H3(1.7%) + H4(1.7%)    100% 2py                 93.6    <- HOMO
```

Unmistakably the in-plane LP: almost all on oxygen, pure 2py, with
small tails on carbon and each hydrogen.

```
  Bond         Total       σ       π       δ              (interference)
  C1-O2         2.038   0.975   1.063   0.000  (+0.008: σ+0.004, π+0.004, δ+0.000)
  C1-H3         0.925   0.955  -0.030   0.000  (-0.072: σ-0.040, π-0.033, δ+0.000)
  C1-H4         0.925   0.955  -0.030   0.000  (-0.072: σ-0.040, π-0.033, δ+0.000)
  O2-H3         0.062   0.000   0.062   0.000  (-0.004: σ-0.002, π-0.002, δ+0.000)
  O2-H4         0.062   0.000   0.062   0.000  (-0.004: σ-0.002, π-0.002, δ+0.000)
```

The C–H bonds pay for the visit. Each drops to 0.925 — the largest
depletion anywhere in the gallery — and picks up a trace of π
character a plain single bond has no business having. The detail
section names the source directly: `C-H σ × O(LP): −0.0625` on each
leg, all subtractive — the lone pair donating straight into the C–H
σ* framework. The through-space O···H contacts (0.062, pure π) are
the same donation read from the other end. Donation reads as
depletion, once more: the donors give, the donors weaken.

The C=O bond, meanwhile, over-delivers — π character above 1.0, more
than a textbook full double bond, fed by the same oxygen that's
donating sideways into the C–H bonds. The charges undersell it:
C +0.207, O −0.280 is a modest split for the canonical polar double
bond. The orbitals tell a more active story than the charges do —
the density stays shared even while the lone pair is visibly
reaching into three different places at once.

![Formaldehyde HOMO — in-plane lone pair](img/formaldehyde_homo.png)
*Orbital 8 (HOMO): the oxygen in-plane lone pair, 93.5% O with 3.1%
on C and 1.7% on each H — the C–H tails rendered.*

---
*Geometry optimized in ORCA 6.1.1 (wB97X-D3/def2-TZVP), confirmed
minimum by frequency calculation (no imaginary modes). IBO analysis
with Psi4 1.11 (wB97X-D/def2-TZVP, RHF), avo_ibo 0.4.0; Pipek–Mezey
localization (p = 2 → 4) per Knizia 2013.*
