# Alkenes: hyperconjugation across the substitution ladder

Textbook organic chemistry ranks alkene stability ethene < propene <
isobutene and credits hyperconjugation: methyl C–H σ bonds donating
into the C=C π* acceptor. These three molecules render that ranking as
IBO numbers at uniform wB97X-D/def2-TZVP — the same donation physics
as the tert-butyl cation's 6+3 split ([carbocations.md](carbocations.md)),
now with a π acceptor instead of an empty p. Reproduce with the
commands shown. (Ethene's section stays silent throughout — no bond
reaches the interference floor — while propene and isobutene fire
detail lines plus footnotes; [how-to-read-ibos.md](how-to-read-ibos.md#6-significant-orbital-pair-interference-when-present)
explains the display convention.)

## Ethene, C₂H₄ — the zero-donor control

```
Input: ethene.xyz  (ORCA 6.1.1 wB97X-D3/def2-TZVP opt + freq, no imaginary modes)
Run:   pixi run python -m avogadro_ibo examples/ethene.xyz --method wB97X-D --basis def2-TZVP
```

The baseline, and deliberately boring:

```
  C1-C2         2.034   1.034   1.000  (+0.017: σ+0.017, π+0.000)
  C1-H3         0.978   0.978   0.000  (-0.009: σ-0.009, π+0.000)   (×4)
```

A full double bond — σ 1.034, π exactly 1.000 — with four equivalent
C–H bonds and no interference detail section at all. Nothing donates
because there is nothing to donate. Every deviation below is measured
against this page.

## Propene, C₃H₆ — one methyl donor

```
Input: propene.xyz  (ORCA 6.1.1 wB97X-D3/def2-TZVP opt + freq, no imaginary modes)
Run:   pixi run python -m avogadro_ibo examples/propene.xyz --method wB97X-D --basis def2-TZVP
```

```
  C2-C3         1.984   1.013   0.971  (-0.017: σ-0.002, π-0.015)
  C1-C2         1.040   1.028   0.012  (+0.016: σ+0.007, π+0.009)
  C1-H4         0.973   0.976  -0.003  (-0.007: σ-0.004, π-0.003)
  C1-H6         0.973   0.976  -0.003  (-0.007: σ-0.004, π-0.003)
  C1-H5         0.983   0.983   0.000  (-0.007: σ-0.007, π+0.000)
```

Donation reads on both ends: the π order drops 1.000 → 0.971 while the
methyl–vinyl single rises to 1.040 carrying genuine π character
(0.012) — a single bond with a hundredth of a double bond in it. Two
of the three methyl C–H bonds show σ/π mixing (−0.003 π) from overlap
with the acceptor; the third, held away from it, stays pure σ at
0.983. And the through-space contacts print:

```
  C3-H4         0.022   0.012   0.010  (+0.010: σ+0.005, π+0.005)
  C3-H6         0.022   0.012   0.010  (+0.010: σ+0.005, π+0.005)
```

the hyperconjugative reach made quantitative.

The detail section fires at full threshold — no probe needed:

```
  C2-C3: orb10(C-H σ) × orb12(C-C π): -0.0150
  C2-C3: orb9(C-H σ) × orb12(C-C π): -0.0150
  C3-H4: orb10(C-H σ) × orb12(C-C π): +0.0100
  C3-H6: orb9(C-H σ) × orb12(C-C π): +0.0100
  (11 terms in [0.005, 0.01) omitted; largest: C1-C2: orb10(C-H σ) × orb12(C-C π) = +0.0087)
```

The familiar compete/cooperate signature: the donor pairs erode the
double bond (−0.0150) while building the contacts (+0.0100) — the same
sign pattern as diborane's three-sign story, now written by a methyl
group and a π bond.

## Isobutene, C₄H₈ — two methyl donors

```
Input: isobutene.xyz  (2-methylpropene; ORCA 6.1.1 wB97X-D3/def2-TZVP opt + freq, no imaginary modes)
Run:   pixi run python -m avogadro_ibo examples/isobutene.xyz --method wB97X-D --basis def2-TZVP
```

```
  C1-C4         1.936   0.997   0.939  (-0.041: σ-0.016, π-0.026)
  C1-C2         1.028   1.016   0.012  (+0.008: σ+0.000, π+0.008)
  C1-C3         1.028   1.016   0.012  (+0.008: σ+0.000, π+0.008)
  C2-H5         0.973   0.976  -0.003  (-0.008: σ-0.005, π-0.003)
  C2-H6         0.973   0.976  -0.003  (-0.008: σ-0.005, π-0.003)
  C3-H8         0.973   0.976  -0.003  (-0.008: σ-0.005, π-0.003)
  C3-H9         0.973   0.976  -0.003  (-0.008: σ-0.005, π-0.003)
  C2-H7         0.984   0.984   0.000  (-0.006: σ-0.006, π+0.000)
  C3-H10        0.984   0.984   0.000  (-0.006: σ-0.006, π+0.000)
```

The trend continues monotonically: π 1.000 → 0.971 → 0.939, with two
methyl–vinyl singles each carrying the same 0.012 of π character
propene showed once. The donor C–H bonds split as before — four mixed
(0.973, −0.003 π) against two pure-σ spectators (0.984) — and the
detail section prints four aligned-donor pairs:

```
  C1-C4: orb13(C-H σ) × orb16(C-C π): -0.0128
  C1-C4: orb12(C-H σ) × orb16(C-C π): -0.0128
  C1-C4: orb11(C-H σ) × orb16(C-C π): -0.0128
  C1-C4: orb10(C-H σ) × orb16(C-C π): -0.0128
  (20 terms in [0.005, 0.01) omitted; largest: C4-H6: orb13(C-H σ) × orb16(C-C π) = +0.0097)
```

backed by four 0.021 (+0.010) vinyl···H contacts:

```
  C4-H6         0.021   0.010   0.011  (+0.010: σ+0.005, π+0.005)
  C4-H8         0.021   0.010   0.011  (+0.010: σ+0.005, π+0.005)
  C4-H9         0.021   0.010   0.011  (+0.010: σ+0.005, π+0.005)
  C4-H5         0.021   0.010   0.011  (+0.010: σ+0.005, π+0.005)
```

Twice the donors, twice the fingerprints — at ~15% smaller per-donor
magnitude (−0.0128 vs propene's −0.0150), consistent with two methyls
sharing one π* acceptor. Total π depletion scales evenly (0.029 →
0.061); per-donor terms do not, and the page no longer claims they do.

## The comparison

| Property | Ethene | Propene | Isobutene |
|----------|--------|---------|-----------|
| C=C total (σ / π) | 2.034 (1.034 / 1.000) | 1.984 (1.013 / 0.971) | 1.936 (0.997 / 0.939) |
| Methyl–vinyl single (π part) | — | 1.040 (0.012) | 1.028 ×2 (0.012) |
| Donor C–H range | 0.978 (pure σ) | 0.973–0.983 (2 mixed) | 0.973–0.984 (4 mixed) |
| Max detail term | — (silent) | −0.0150 | −0.0128 ×4 |
| Near-miss footnote | — | 11 terms | 20 terms |

The stability ordering falls out of the π column directly: each methyl
costs the double bond ~0.03 of π order and parks ~0.012 of π character
on its single bond.

## References

- R. S. Mulliken, C. A. Rieke, W. G. Brown, "Hyperconjugation,"
  *J. Am. Chem. Soc.* **1941**, *63*, 41–56,
  DOI:[10.1021/ja01846a008](https://doi.org/10.1021/ja01846a008) —
  the original proposal that alkyl groups donate into unsaturated
  systems.

---
*Geometries optimized in ORCA 6.1.1 (wB97X-D3/def2-TZVP), confirmed
minima by frequency calculation (no imaginary modes); PubChem 3D
starting coordinates (CIDs 6325, 8252, 8255), DFT-relaxed before
analysis. IBO analysis with Psi4 1.11 (wB97X-D/def2-TZVP, RHF),
avo_ibo 0.4.0; Pipek–Mezey localization (p = 2 → 4) per Knizia 2013.*
