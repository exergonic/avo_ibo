# Ferrocene — ten equivalent bonds, and where the partition stops

Ferrocene is the gallery's transition-metal exam: ten Fe–C contacts
related by D5h symmetry, d-electron donation both ways, and a
localization functional with no knowledge of point groups. It passes —
with one honest failure that extends the tool's known limitation into
new territory.

```
Input: ferrocene.xyz  (ORCA 6.1.1 wB97X-D3/def2-TZVP opt + freq from exact D5h, no imaginary modes)
Run:   pixi run python -m avogadro_ibo examples/ferrocene.xyz --method wB97X-D --basis def2-TZVP
```

Two legs, following Knizia 2013 Appendix E. The unoptimized
construction from experimental distances is disclosed, not hidden.
Leg 1 is exact D5h from Fe–C 2.064 / C–C 1.440 / C–H 1.08 Å.
Leg 2 is ORCA-relaxed from that start (Fe–C 2.050, C–C 1.4196).
D5h stays intact to 1e-4 with no symmetry constraints — the
functional respects the group unassisted, so no constrained rerun
was needed.

**The density totals respect D5h perfectly, in both legs:**

```
  Bond         Total       σ       π              (interference)
  Fe1-C18        0.534   0.533   0.001  (-0.012: σ-0.008, π-0.004)
  Fe1-C6         0.534   0.239   0.295  (-0.041: σ-0.024, π-0.018)
  Fe1-C20        0.534   0.459   0.075  (-0.047: σ-0.068, π+0.022)
  Fe1-C2         0.534   0.533   0.001  (-0.012: σ-0.008, π-0.004)
```

Ten contacts, identical to three decimals — half-bonds. The charges
match: Fe −0.238, every carbon −0.055, every hydrogen +0.080. The
donation pairs read cleanly too. Six "Deloc" orbitals at ~22% Fe
bridge ring π and metal d. The nonbonding dz² (orb 46, 98.6%) sits
apart exactly where ligand-field theory puts it. The HOMO/LUMO is an
Fe–C σ/σ* pair (−0.292/+0.163 Ha).

![Ferrocene dz² — nonbonding a1g orbital](img/ferrocene_dz2.png)
*Orbital 46: the nonbonding dz², 98.6% Fe — the a1g orbital
ligand-field theory predicts, untouched by either ring.*

**But the σ/π partition does not respect D5h — and that is the
finding.** Same-total bonds decompose incompatibly: the table's first
row is almost pure σ, the second is mixed. Worse, the assignment
*permutes between legs*. Leg 1 puts the π-heavy character on Fe1–C20;
Leg 2 moves it to Fe1–C6. Totals are density-protected and
reproducible. The partition is localization-trajectory-dependent
wherever symmetry-equivalent bonds share overlapping Fe-d/C-p
subspaces. This is the PM symmetry limitation (benzene C–H split,
SO₃ d-asymmetry — the gallery's standing caveat) in a new place: not
energies or coefficients, but the decomposition itself.

The writeup rule that follows: quote totals and compositions for
ferrocene, never per-bond σ/π splits. Switching functionals would not
help. Boys localization cannot separate σ from π by construction — it
returns banana bonds for multiples — so it would abolish the
distinction rather than symmetrize it. A documented limitation with a
characterized failure mode beats a hidden one.

---
*Leg-1 geometry: exact D5h from experimental parameters (Fe–C 2.064,
C–C 1.440, C–H 1.08 Å), unoptimized per Knizia 2013 Appendix E.
Leg-2 geometry optimized in ORCA 6.1.1 (wB97X-D3/def2-TZVP),
confirmed minimum by frequency calculation (no imaginary modes).
IBO analysis with Psi4 1.11 (wB97X-D/def2-TZVP, RHF), avo_ibo 0.4.0;
Pipek–Mezey localization (p = 2 → 4) per Knizia 2013.*
