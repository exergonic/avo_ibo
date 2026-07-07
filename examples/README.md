# avo_ibo Examples

This directory contains pre-computed IBO analyses demonstrating the chemical
insights that IAO/IBO decomposition provides.  Each `.xyz` input is provided
so you can reproduce the results with:

    pixi run python -m avogadro_ibo examples/molecule.xyz


## Table of Contents

| Molecule | Chemical Feature | File |
|----------|-----------------|------|
| Sulfur trioxide | Hypervalency, π ionicity gradient | `SO3.xyz` |
| Water | On-atom degeneracy resolution | `water.xyz` |
| Methane | sp³ hybridization | `methane.xyz` |
| Ethene | σ/π separation | `ethene.xyz` |
| Benzene | Delocalized π system | `benzene.xyz` |
| Ammonia | Lone pair character | `ammonia.xyz` |
| ZnCl₂ | 3d¹⁰ transition metal | `zncl2.xyz` |


## Sulfur Trioxide (SO₃) — Hypervalency and Ionicity Gradients

```
Input: molecules/SO3.xyz  (4 atoms, D3h symmetry)
Run:   pixi run python -m avogadro_ibo examples/SO3.xyz --method hf --basis cc-pVDZ
```

```
     #      Occ      Energy              Type  Composition                         Hybrid                 Ion%     H/L
-------------------------------------------------------------------------------------------------------------------------
     1    2.000  -92.297548           S(Core)  S1(100.0%)                         100% 1s                 ---    
     2    2.000  -20.644180         O(Core) †  O3(100.0%)                         100% 1s                 ---    
     3    2.000  -20.644180         O(Core) †  O4(100.0%)                         100% 1s                 ---    
     4    2.000  -20.644172         O(Core) †  O2(100.0%)                         100% 1s                 ---    
     5    2.000   -9.239443             S(LP)  S1(100.0%)                         100% 2s                 ---    
     6    2.000   -6.960065             S(LP)  S1(100.0%)                         100% 2pz                ---    
     7    2.000   -6.957016           S(LP) †  S1(100.0%)                         100% 2py                ---    
     8    2.000   -6.957016           S(LP) †  S1(100.0%)                         100% 2px                ---    
     9    2.000   -1.231823           S-O σ †  O2(57.0%) + S1(42.4%)              21% 2s + 79% 2px       14.7   
    10    2.000   -1.231642           S-O σ †  O4(57.0%) + S1(42.5%)              21% 2s + 79% 2py       14.6   
    11    2.000   -1.231642           S-O σ †  O3(57.0%) + S1(42.5%)              21% 2s + 79% 2py       14.6   
    12    2.000   -1.044461           O(LP) †  O3(99.7%)                          79% 2s + 21% 2py       99.6   
    13    2.000   -1.044461           O(LP) †  O4(99.7%)                          79% 2s + 21% 2py       99.6   
    14    2.000   -1.044386           O(LP) †  O2(99.7%)                          79% 2s + 21% 2px       99.6   
    15    2.000   -0.624919           S-O π †  O2(81.9%) + S1(16.3%)              100% 2pz               66.8   
    16    2.000   -0.624864           S-O π †  O3(81.9%) + S1(16.3%)              100% 2pz               66.8   
    17    2.000   -0.624864           S-O π †  O4(81.9%) + S1(16.3%)              100% 2pz               66.8   
    18    2.000   -0.609773           O(LP) †  O2(92.4%) + S1(5.2%)               100% 2py               89.4   
    19    2.000   -0.609734           O(LP) †  O3(92.4%) + S1(5.2%)               100% 2px               89.4   
    20    2.000   -0.609734           O(LP) †  O4(92.4%) + S1(5.2%)               100% 2px               89.4  <- HOMO
    21    0.000    0.099285              S π*  S1(51.1%) + O2(16.3%) + O4(16.3%) + O3(16.3%)  100% 3pz  ---  <- LUMO
    22    0.000    0.486388          S-O σ* †  S1(52.2%) + O4(39.2%) + O2(4.3%) + O3(4.3%)  27% 3s + 73% 3py ---    
    23    0.000    0.486388          S-O σ* †  S1(52.2%) + O3(39.2%) + O2(4.3%) + O4(4.3%)  27% 3s + 73% 3py ---    
    24    0.000    0.486572          S-O σ* †  S1(52.2%) + O2(39.2%) + O3(4.3%) + O4(4.3%)  27% 3s + 73% 3px ---    

Total electrons: 40

--- Charge Decomposition ---
  Atom    Z       Pop  Net Charge
---------------------------------
     S   16    13.849      +2.151
     O    8     8.717      -0.717
     O    8     8.717      -0.717
     O    8     8.717      -0.717
---------------------------------
Total:   40    40.000      +0.000

--- Total Wiberg Bond Orders ---
  S-O      1.313
  S-O      1.313
  S-O      1.313
  O-O      0.241
  O-O      0.241
  O-O      0.241
```

### Analysis

**The ionicity gradient.**  The S-O σ bonds (orbitals 9-11) are 14.7% ionic, the
S-O π bonds (15-17) are 66.8% ionic, and the oxygen 2p lone pairs (18-20) are
89.4% ionic.  The σ framework is relatively covalent, but the π system is
highly polar (electron density overwhelmingly on oxygen).  The 2s lone pairs
(12-14) are essentially oxide-like at 99.6% ionic.

The formal charge is S=+2.151, O=−0.717 each — consistent with the resonance
structures students draw (two S=O double bonds + one S→O coordinate bond).

**Wiberg O-O bond order of 0.241.**  The three oxygens communicate through the
delocalised π system.  For comparison, benzene's meta C-C through-space bonding
is 0.115 — SO₃'s O-O is more than double that, reflecting how strongly
polarised the S-O π system is toward oxygen, creating significant π-orbital
interaction between the oxygen centres.

**S 2p core orbitals (6-8).**  Orbital 6 (2pz) is non-degenerate while 7-8
(2px, 2py) are a degenerate pair.  This is chemically meaningful: 2pz is
perpendicular to the molecular plane and interacts differently with the Σ
bonding framework than the in-plane 2px/2py.


## Water (H₂O) — On-Atom Degeneracy Resolution

```
     #      Occ      Energy              Type  Composition         Hybrid              Ion%     H/L
------------------------------------------------------------------------------------------------------
     1    2.000  -20.534460           O(Core)  O1(100.0%)         100% 1s              ---    
     2    2.000   -0.916099           O-H σ †  O1(61.7%)+H2(38.3%) 23% 2s+77% 2py      23.3   
     3    2.000   -0.916099           O-H σ †  O1(61.7%)+H3(38.3%) 23% 2s+77% 2py      23.3   
     4    2.000   -0.788503             O(LP)  O1(100.0%)         55% 2s+45% 2pz        ---    
     5    2.000   -0.493602             O(LP)  O1(100.0%)         100% 2px              ---  <- HOMO
     6    0.000    0.587632          O-H σ* †  H3(61.6%)+O1(38.3%) O:22%2s+78%2py      ---  <- LUMO
     7    0.000    0.587632          O-H σ* †  H2(61.6%)+O1(38.3%) O:22%2s+78%2py      ---    

Net charges: O=−0.466, H=+0.233 each
```

The PM functional alone cannot separate orbitals on the same atom — O 2s and
O 2p lone pairs have DOM=1 on the same oxygen and would mix arbitrarily.
On-atom Fock diagonalisation (orbital 4 at −0.79 Ha, 55% 2s) resolves them:
the s-rich LP is 0.3 Ha lower than the pure 2p LP (orbital 5, −0.49 Ha),
matching the aufbau principle.


## Methane (CH₄) — sp³ Hybridisation

```
     #      Occ      Energy          Type  Composition       Hybrid          Ion%     H/L
-------------------------------------------------------------------------------------------------------
     1    2.000  -11.193312           C(Core)  C1(100.0%)   100% 1s          ---    
     2    2.000   -0.648272           C-H σ †  C1(51.8%)+H5(48.2%) 25%2s+75%2p  3.7   
     3    2.000   -0.648272           C-H σ †  C1(51.8%)+H3(48.2%) 25%2s+75%2p  3.7   
     4    2.000   -0.648272           C-H σ †  C1(51.8%)+H2(48.2%) 25%2s+75%2p  3.7   
     5    2.000   -0.648272           C-H σ †  C1(51.8%)+H4(48.2%) 25%2s+75%2p  3.7  <- HOMO

Net charges: C=−0.147, H=+0.037 each
Wiberg C-H: 0.999
```

The 25% 2s + 75% 2p hybridisation is exactly sp³ (s:p = 1:3).  All four C-H
bonds are degenerate in energy (−0.648 Ha), confirming PM localisation finds
the symmetric solution for a tetrahedral molecule.  The 3.7% ionic character
indicates slightly polar C→H bonding consistent with the electronegativity
difference.


## Ethene (C₂H₄) — σ/π Separation

```
     #      Occ      Energy          Type  Composition     Hybrid           Ion%     H/L
---------------------------------------------------------------------------------------------------
     1    2.000  -11.207594         C(Core) †  C1(100.0%)  100% 1s          ---    
     2    2.000  -11.207594         C(Core) †  C2(100.0%)  100% 1s          ---    
     3    2.000   -0.866714             C-C σ  C1(50.0%)+C2(50.0%) 38%2s+62%2pz  0.0   
     4    2.000   -0.684189           C-H σ †  C1(52.3%)+H4(47.3%) 32%2s+68%2py  5.0   
     5    2.000   -0.684189           C-H σ †  C2(52.3%)+H6(47.3%) 32%2s+68%2py  5.0   
     6    2.000   -0.684189           C-H σ †  C1(52.3%)+H3(47.3%) 32%2s+68%2py  5.0   
     7    2.000   -0.684189           C-H σ †  C2(52.3%)+H5(47.3%) 32%2s+68%2py  5.0   
     8    2.000   -0.374804             C-C π  C2(50.0%)+C1(50.0%) 100%2px      0.0  <- HOMO

Wiberg C-C: 2.028 (σ + π contributions)
```

The C-C σ bond (orbital 3) at −0.87 Ha uses 38% 2s + 62% 2p (sp²-like). The
C-C Σ bond (orbital 8) at −0.37 Ha is 100% 2p — pure p-orbital side-on
overlap.  The total Wiberg C-C bond order of 2.028 correctly sumarises the
σ+π double bond.

Four C-H σ bonds are all degenerate at −0.684 Ha with 32% 2s + 68% 2p (sp²-like
hybridisation), each carrying 5.0% ionic character toward hydrogen.


## Ammonia (NH₃) — Lone Pair Character

```
     #      Occ      Energy          Type  Composition     Hybrid           Ion%     H/L
---------------------------------------------------------------------------------------------------
     1    2.000  -15.495619           N(Core)  N1(100.0%)  100% 1s          ---    
     2    2.000   -0.796391           N-H σ †  N1(57.8%)+H2(42.2%) 28%2s+72%2py  15.7   
     3    2.000   -0.796273           N-H σ †  N1(57.8%)+H4(42.2%) 28%2s+72%2px  15.7   
     4    2.000   -0.796273           N-H σ †  N1(57.8%)+H3(42.2%) 28%2s+72%2px  15.7   
     5    2.000   -0.444227             N(LP)  N1(100.0%)  18%2s+82%2pz    ---  <- HOMO

Net charges: N=−0.471, H=+0.157 each
Wiberg N-H: 0.975
```

The N lone pair (orbital 5, HOMO) is 18% 2s + 82% 2p — predominantly p-type.
This is consistent with the molecular geometry (pyramidal, ∠HNH ≈ 107°) and
the expectation that the lone pair occupies an orbital with more p character
in a non-planar molecule.  The N-H bonds show 15.7% ionic character,
consistent with the Pauling electronegativity difference ΔEN≈0.8.


## Benzene (C₆H₆) — Delocalised π System

```
     #      Occ      Energy          Type  Composition                         Hybrid     Ion%     H/L
------------------------------------------------------------------------------------------------------------------
    13    2.000   -0.693554           C-H σ †  C6(52.3%)+H12(47.2%)            32%2s+68%2px  5.1   
   ...
    19    2.000   -0.388091           Deloc †  C2(50.0%)+C3(22.2%)+C1(22.2%)+C5(5.6%) 100%2pz  38.5   
    20    2.000   -0.388091           Deloc †  C6(50.0%)+C5(22.2%)+C1(22.2%)+C3(5.6%) 100%2pz  38.5   
    21    2.000   -0.388085           Deloc †  C4(50.0%)+C3(22.2%)+C5(22.2%)+C1(5.6%) 100%2pz  38.5  <- HOMO
    22    0.000    0.318236            C π* †  C1(50.0%)+C6(22.2%)+C2(22.2%)+C4(5.6%) 100%2pz  ---  <- LUMO
    23    0.000    0.318242            C π* †  C5(50.0%)+C6(22.2%)+C4(22.2%)+C2(5.6%) 100%2pz  ---    
    24    0.000    0.318242            C π* †  C3(50.0%)+C2(22.2%)+C4(22.2%)+C6(5.6%) 100%2pz  ---    

Wiberg C-C (ring): 1.444  (between adjacent carbons)
Wiberg C-C (meta): 0.115  (through-bond coupling)
```

The three occupied π orbitals are 100% pz and delocalised over the ring.  The
"50% + 22.2% + 22.2% + 5.6%" composition tells you they're anealing-mode
superpositions of localised ethylene-like π bonds.  The 38.5% ionic character
on these "Deloc" labels is a technical artefact — the PM functional cannot
further localise a fully delocalised system, and the 50-22-22 composition
represents a natural stationary point of the PM functional on the π subspace.

The adjacent C-C Wiberg order of 1.444 (σ + one-third of the π system per
bond) is consistent with a bond order of ∼1.5.  The meta C-C Wiberg (0.115)
is through-bond coupling — non-neighbour carbons communicating via the
σ framework.


## ZnCl₂ — 3d¹⁰ Transition Metal

```
Input: molecules/zncl2.xyz  (3 atoms, D∞h symmetry)
Run:   pixi run python -m avogadro_ibo examples/zncl2.xyz --method hf --basis cc-pVDZ
```

```
     #      Occ      Energy          Type  Composition          Hybrid                Ion%     H/L
----------------------------------------------------------------------------------------------------------------
    26    2.000   -0.846207            Zn(LP)  Zn2(100.0%)     5% 4s + 95% 3dz²       ---    
    27    2.000   -0.706144         Zn-Cl σ †  Cl3(71.0%)+Zn2(29.0%) 17%3s+83%3pz    41.9   
    28    2.000   -0.706144         Cl-Zn σ †  Cl1(71.0%)+Zn2(29.0%) 17%3s+83%3pz    41.9   
    29    2.000   -0.467437          Cl(LP) †  Cl3(96.2%)+Zn2(3.8%) 100%3px           92.5   
    30    2.000   -0.467437          Cl(LP) †  Cl3(96.2%)+Zn2(3.8%) 100%3py           92.5   
    31    2.000   -0.467437          Cl(LP) †  Cl1(96.2%)+Zn2(3.8%) 100%3py           92.5   
    32    2.000   -0.467437          Cl(LP) †  Cl1(96.2%)+Zn2(3.8%) 100%3px           92.5  <- HOMO

Wiberg Zn-Cl: 1.104
Net charges: Zn=+0.529, Cl=−0.264 each
```

Zinc has a full d¹⁰ shell — all five 3d orbitals (22-26) appear as pure-atom
"Zn(LP)" with 100% d-character.  The lowest d-orbital (26) is 5% 4s + 95% 3dz²
— the dz² mixes with the 4s due to its Σ-symmetry along the molecular axis.
The Zn-Cl σ bonds (27-28) are 41.9% ionic, polarised toward chlorine.  The
chlorine 3p lone pairs (29-32) are 92.5% ionic — essentially Cl⁻ with a tiny
ligand-to-metal σ-donation tail.

The Wiberg order of 1.104 (slightly above 1) reflects the σ-donation from Cl
to Zn, consistent with the formal oxidation state of Zn²⁺ in a linear complex.
