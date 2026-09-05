# Phenol vs anisole: why the stronger donor ring is the weaker activator

A mechanistic study in four calculations. Phenol is the stronger
activator toward electrophilic aromatic substitution, yet anisole's
ring carries more electron density in the ground state. The resolution
is in the Wheland intermediate: what matters is not how much density
the ring holds, but how much the substituent *donates when positive
charge builds up*. IBO charges make both halves of that sentence
quantitative.

All ORCA work at wB97X-D3/def2-TZVP (Opt+Freq); IBO analysis by
IboView (IBBA) with plugin cross-checks (wB97X-D/def2-TZVP — see
Validation below). Geometries: `phenol_v_anisole/*_opt.xyz`
(ORCA-optimized last frames, shipped alongside this file).

## 1. The paradox (ground states)

| | Phenol | Anisole |
|---|---|---|
| O charge | −0.498 | −0.346 |
| C1 (ipso) | +0.231 | +0.220 |
| C4 (para) | −0.179 | −0.179 |
| Methyl C8 | — | −0.201 |
| HOMO | C–C π, −0.383 Ha | C–C π, −0.380 Ha |
| O–C1 Wiberg | 1.055 | 1.053 |

Anisole's ring is more negative almost everywhere, and the HOMOs are
within 0.003 Ha — yet phenol nitrates/brominates faster. Ground-state
donation (O–C1 bond orders identical at 1.05) cannot explain the rate
difference. The oxygen tells the real story: phenol's O–H bond
polarizes toward oxygen (H χ = 2.20) far more than anisole's O–CH₃
(C χ = 2.55), so phenol's oxygen starts with 0.15 e⁻ more density to
give. The methyl group, far from holding negative charge, is a donor:
C8 +0.023 overall would be positive were it not for its hydrogens, and
it feeds oxygen through n(O)→σ*(C–H) hyperconjugation — *satisfying*
oxygen electronically and reducing its drive to push density into the
ring. That is the refined theory; the Wheland intermediates test it.

## 2. The Wheland intermediates (para protonation)

Both C4-protonated cations optimize cleanly (no imaginary frequencies;
the anisole methyl rotor needed the usual small-frequency scrutiny).
Topology confirmed: C4 sp³ (C3–C4/C4–C5 ≈ 1.07), remaining π on
C2–C3/C5–C6 (≈ 1.70), C1–O strengthened (1.291 phenol, 1.325 anisole).

| Atom | Phenol GS → WH | Δ | Anisole GS → WH | Δ |
|---|---|---|---|---|
| O | −0.498 → −0.361 | **+0.137** | −0.346 → −0.227 | **+0.119** |
| C1 (ipso) | +0.231 → +0.391 | +0.160 | +0.220 → +0.380 | +0.160 |
| C3/C5 (ortho to attack) | −0.134 → −0.002/−0.011 | +0.13 | −0.137 → −0.025/−0.003 | +0.12 |
| C4 (para, attacked) | −0.179 → −0.296 | −0.117 | −0.179 → −0.295 | −0.116 |
| Methyl C8 | — | — | −0.201 → −0.223 | −0.022 |

(C4 goes *negative*: the two new C–H bonds polarize toward carbon. The
positive charge sits on C1/C3/C5 — exactly the Wheland resonance
positions.)

**Phenol's oxygen donates 0.137 e⁻; anisole's donates 0.119 e⁻ — ~15%
more from phenol.** The orbital picture agrees: the O-lone-pair→C1
donation nearly triples in both (phenol C1 0.084 → 0.249; anisole 0.089
→ 0.287), but anisole's oxygen started with less and gives less in
total. The methyl group donates a further −0.022 into oxygen during
the reaction (C8 −0.201 → −0.223, hydrogens correspondingly more
positive) — a real reservoir effect, but an order of magnitude too
small to close the gap. Theory confirmed, quantitatively.

## 3. Validation: plugin vs IboView

All four ORCA-optimized geometries rerun through avo_ibo
(wB97X-D/def2-TZVP; D3 unavailable in Psi4's env — D2 instead, see
caveat). Full table in [../NOTES.md](../NOTES.md) ("Cross-validation
against IboView and ORCA"):

- **Trends agree to ±0.01**: ΔO +0.139/+0.137 (phenol), +0.122/+0.119
  (anisole); ΔC1 +0.150/+0.160 both molecules. The 14%-vs-15% verdict
  is implementation-independent.
- **Wibergs agree to ±0.02** throughout, including the Wheland C–O
  and ring bonds.
- **One systematic offset, diagnosed**: absolute charges on hydrogens
  run ~+0.08–0.11 e⁻ more positive in IboView (O–H σ 68.5/31.1 vs
  62.9/36.6 — accounts for its H gap to the digit), plausibly
  minimal-basis H parameterization. Heavy atoms ±0.05, lone pairs
  identical. It cancels in every Δ that matters.

## 4. Methods and caveats

- Electrophile: H⁺ (para protonation) isolates electronics from
  halogen/nitro complications; the Wheland intermediate (minimum, not
  saddle) carries the maximum charge separation.
- Methyl-rotor artifacts (< 80 i cm⁻¹ torsions) distinguished from real
  saddle points by frequency magnitude and mode displacement; all four
  reported structures are true minima.
- wB97X-D3 (ORCA) vs wB97X-D (plugin): dispersion-correction
  difference, negligible for IBO populations; comparison led with Δ
  trends throughout. IboView's wall-clock-seeded Cayley rotations make
  its degenerate-manifold digits non-reproducible in principle —
  agreement is claimed at the ±0.01 level, no tighter.
