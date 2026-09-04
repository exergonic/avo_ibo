# Validation

We checked this plugin against two independent IBO implementations:
Knizia's [IboView](http://www.iboview.org) and ORCA 6.1.1 (`%loc` /
`IAOIBO`). Five molecules, one level of theory (wB97X-D/def2-TZVP),
same geometry in every leg. The full record lives in `NOTES.md`;
this page states what a user needs to know.

## Charges agree to 0.006

IAO partial charges, ours vs ORCA with the same minimal basis (STO-3G):

| Molecule | Atom | Ours | ORCA | Difference |
|---|---|---|---|---|
| Water | O / H | −0.494 / +0.247 | −0.4896 / +0.2448 | 0.004 |
| Ethene | C / H | −0.131 / +0.066 | −0.1286 / +0.0643 | 0.002 |
| Ammonia | N / H | −0.530 / +0.177 | −0.5236 / +0.1746 | 0.006 |
| Benzene | C / H | −0.065 / +0.065 | −0.0636 / +0.0636 | 0.001 |
| Carbon monoxide | C / O | +0.069 / −0.069 | +0.0679 / −0.0679 | 0.001 |

Worst difference across all molecules: **0.006 electrons** (nitrogen).
That residual is SCF-grid and convergence trivia between two different
programs, not method error.

## The minimal basis is a convention, not physics

IboView builds its IAOs over the MINAO basis; we build ours over
STO-3G (native to Psi4); ORCA defaults to yet another choice (atomic
SCF orbitals) and offers all of them as options. ORCA's own manual
calls its default charges only "very similar to the original IAO
charges." So cross-program comparisons must fix the convention first:

- Same convention (STO-3G vs STO-3G): digit agreement, table above.
- Different convention (MINAO vs STO-3G): systematic offsets, always
  in the same direction (MINAO more polarized) — 0.16 on O–H, ~0.2 on
  N–H, 0.12 on C–O. Bond assignments, trends, and differences between
  molecules are unaffected.

We do not chase MINAO. STO-3G is documented, robust, and
characterized; agreement-with-one-program is not correctness.

## Bond orders agree in pattern; our σ/π split stands alone

Our density-Wiberg totals track ORCA's Mayer bond orders (C–C 2.03 vs
1.96, benzene C–C 1.44 vs 1.41) — same multiplicities, digits differ
because the two quantities are built in different bases. No other
program partitions σ from π; IboView's own orbitals confirm our
partition independently (CO σ-like 65/35 vs our 62/38, π-like 72/28
vs our 71/29).

## What we do not claim

- Localized orbital *energies* are not comparable across programs
  (definition- and basin-dependent). One exception: symmetry-protected
  orbitals must match, and water's HOMO does (−0.401 vs −0.414).
- Feeding our Molden files to IboView works for small hydrides but
  degrades with heavy f-basis content (Psi4 writes raw contraction
  coefficients; see `validation/molden_renorm.py` and `NOTES.md`).

Bottom line: same convention gives the same numbers to 0.006.
The implementation is sound; the remaining differences are named
conventions with measured sizes.

## "But who is right?"

We all are. An IAO charge has no single true value waiting to be
found — it is a well-defined number *within a stated convention*
(minimal basis, partitioning scheme). Fix the convention and every
program here gives the same answer. So the question dissolves: pick a
convention, state it, and compare within it.
