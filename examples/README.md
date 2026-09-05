# avo_ibo Examples

Pre-computed IBO analyses demonstrating the chemical insights that IAO/IBO
decomposition provides. New to the output format? Start with
[how-to-read-ibos.md](how-to-read-ibos.md), which walks through every
section of `ibos.txt` using water as the running example.

Reproduce anything here with:

```
pixi run python -m avogadro_ibo examples/<molecule>.xyz --method wB97X-D --basis def2-TZVP [--charge N --spin 1]
```

Results land in `calcs/<molecule>_NNN/`. All outputs below are
wB97X-D/def2-TZVP unless noted.

## Methods and provenance

Small rigid molecules at MMFF94-optimised geometries are adequate inputs
— the IBO picture is robust to the small geometric differences between
MMFF94 and DFT minima. Suspect cases are flagged in their entries:
`ozone.xyz` was recovered from a prior calculation with a short O–O
distance (re-optimise before production use), and the cyclopropenyl
anion pair deliberately contrasts a constrained planar geometry against
the relaxed nonplanar one. Each entry lists its input file and charge
state; orbitals visualised from the matching `ibo.molden`.

## Spotlights

| Molecule | Chemical feature | File |
|----------|-----------------|------|
| Methyl → norbornyl | Bare, bridged, classical, and the historical 3c–2e system; hyperconjugation 6+3 | [carbocations.md](carbocations.md) |
| Ethene → isobutene | Hyperconjugation ladder: π 1.000 → 0.939, methyl–vinyl π character | [alkenes.md](alkenes.md) |
| Cyclopropenyl ± | Aromatic thirds vs symmetry-broken anion; planar/nonplanar pair | [aromaticity-cyclopropenyl.md](aromaticity-cyclopropenyl.md) |
| Cyclopropane | Bent bonds, sp4.6 strain hybrids | [cyclopropane.md](cyclopropane.md) |
| Cyclooctatetraene | Tub hyperconjugation, Karplus-type dihedral dependence | [cot.md](cot.md) |
| Ozone | Charge separation, terminal O···O half-bond | [ozone.md](ozone.md) |
| Malonaldehyde enol | H-bond fingerprint, competing π channels | [malonaldehyde.md](malonaldehyde.md) |
| Allene | Orthogonal π systems, axial chirality | [allene.md](allene.md) |
| Diborane | 3c–2e bridges, cooperation/competition | [diborane.md](diborane.md) |
| Methylamine | Amine lone pair | [methylamine.md](methylamine.md) |
| Water dimer | H-bond covalency: delocalized LP, 0.075 H···O order | [water-dimer.md](water-dimer.md) |

Deferred for later spotlights (data on disk, no writeup yet):
cyclobutadiene, formaldehyde, ferrocene (wB97X-D geometry, slightly
desymmetrised from D5h — the angle is symmetry-broken IBOs from
near-symmetric input; existing table at MN15-L/def2-SVP needs a
wB97X-D/def2-TZVP rerun first).
Separate mechanistic study (own writeup, now alongside the validation
record): [phenol-vs-anisole.md](../validation/phenol-vs-anisole.md)
— why phenol out-activates
anisole in EAS, via Wheland-intermediate IBO charges, with a
plugin-vs-IboView cross-validation appendix.

## First look: seven small molecules

Condensed results. See the linked spotlights and the reading guide for
the full treatment.

### Sulfur Trioxide (SO₃) — hypervalency and ionicity gradients

```
Input: SO3.xyz  (4 atoms, D3h)
```

S–O σ bonds (17.5% ionic) vs S–O π bonds (63.4% ionic) vs O 2s lone
pairs (99.4% ionic): the σ framework is relatively covalent while the π
system is highly polar, density overwhelmingly on oxygen. Charges
S +2.047 / O −0.682 each — the resonance picture (two S=O plus one
coordinate S→O) as numbers. The σ/π Wiberg split reads S–O 1.310 (σ
0.825 + π 0.485) with heavy folded interference (−0.508, mostly π):
three-centre π sharing spread over three equivalent bonds. O–O 0.288
(σ 0.051 + π 0.237) — oxygens communicating through the delocalised π
system, more than double benzene's meta coupling (0.116).

### Water (H₂O) — on-atom degeneracy resolution

Covered in full in [how-to-read-ibos.md](how-to-read-ibos.md): O 2s vs
O 2p lone pairs split by Fock diagonalisation (−0.638 vs −0.401 Ha),
O–H 0.939 all-σ, gap 22.8 eV.

### Methane (CH₄) — sp³ hybridisation

```
    2    2.000   -0.560401           C-H σ †  C1(52.9%) + H5(47.1%)    24% 2s + 76% 2p           5.8
```

24% s + 76% p ≈ sp³, all four bonds degenerate, 5.8% ionic toward H.
Charges C −0.230 / H +0.058; C–H Wiberg 0.997; gap 29.6 eV.

### Ethene (C₂H₄) — σ/π separation

```
    3    2.000   -0.729415             C-C σ  C1(50.0%) + C2(50.0%)    37% 2s + 63% 2pz          0.0
```

C–C σ (37% 2s, sp²-like) at −0.73 Ha; C–C π (100% 2p) is the HOMO; four
degenerate C–H σ (32% 2s, 6.6% ionic). Wiberg C=C 2.033 (σ 1.033 + π
1.000); gap 14.0 eV. See [carbocations.md](carbocations.md) for what
happens when this double bond gets protonated.

### Ammonia (NH₃) — lone pair character

N LP (orb 5, HOMO): 100% on N, 17% 2s + 83% 2pz — predominantly p-type,
consistent with the pyramidal geometry. N–H σ: 28% 2s, 17.7% ionic
(ΔEN ≈ 0.8). Charges N −0.530 / H +0.177; N–H Wiberg 0.969; gap
24.1 eV.

### Benzene (C₆H₆) — delocalised π system

Three occupied π orbitals (100% pz, "Deloc" type, 50/22.2/22.2/5.6
composition — annealing-mode superpositions of ethylenic π bonds),
HOMO at −0.374 Ha. C–C Wiberg 1.444 (σ 1.000 + π 0.444): the aromatic
third of the π system per bond, landing in the π column with no gate
tricks. Meta C–C 0.116 — almost pure π (0.111) through-bond coupling.
Gap 15.2 eV.

### ZnCl₂ — 3d¹⁰ transition metal

Zn d¹⁰ shell intact as pure-atom `Zn(LP)` orbitals (the lowest mixes 8%
4s into 3dz² by axial symmetry). Zn–Cl σ: Cl 68.6% / Zn 31.3%, 37.4%
ionic; Cl 3p LPs 91–99% ionic. Charges Zn +0.386 / Cl −0.193. Wiberg
Zn–Cl 1.172 (σ 0.852 + π 0.319) — above 1 from Cl→Zn σ-donation,
consistent with formal Zn²⁺. Gap 13.2 eV.
