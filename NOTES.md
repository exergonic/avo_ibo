# NOTES

Engineering notebook for avo_ibo: open items, design decisions, and the
reasoning behind them. Claims here are checked against the code; when the
code changes, update the note or delete it.

Status docs live in [AGENTS.md](AGENTS.md); derivations in
[mathematics/mathematics.md](mathematics/mathematics.md); a hands-on
plugin-development guide in [tutorial.md](tutorial.md).

## Open items

- **Typed API** ([#3](https://github.com/exergonic/avo_ibo/issues/3)):
  `compute_ibo()` returns display strings and side-effect files today; a
  typed return would make the core usable as a library.
- **pixi-pack distribution** ([#7](https://github.com/exergonic/avo_ibo/issues/7)):
  self-extracting environment archive so users need neither pixi nor conda.
  Not yet shipped — the v6-lock manual-editing dance (AGENTS.md Gotcha 28)
  is the friction point this would remove.
- **Bond-flat degeneracy resolution** (below): implemented for the occupied
  block (`_resolve_flat_degeneracies`); virtual block intentionally deferred
  until the SVD valence-virtual columns get junk-hygiene (see NOTE in
  `compute_ibo`).

## Next step: resolve bond-flat PM degeneracies

### Observation (2026-08-24, ethene_024 vs ethene_025)

The same molecule localized to σ+π on an idealized D₂h geometry but to two
banana bonds on a slightly distorted one. The vectorization of the PM loops
(commit 02737df) was suspected and exonerated: run against a pre-vectorization
checkout it reproduces both outcomes bit-for-bit. Method/basis is likewise
irrelevant (factorial-tested). Geometry alone flips the result.

### Mechanism

The PM functional sees per-atom populations only. For a symmetric two-center
bond, every rotation inside the {σ, π} plane keeps each orbital's split at
50/50, so the whole plane is stationary. Measured on the distorted ethene
(hf/cc-pVDZ IAOs): L varies by ~6×10⁻⁸ relative across the full 45° rotation,
with a *microscopic* maximum at the banana end — the optimizer legitimately
prefers bananas there by 10⁻⁷-level margins decided by SCF noise steering the
Jacobi trajectory. Idealized geometries escape because their canonical MOs are
symmetry-pure, so B_ij = 0 and the sweeps never rotate at all. Which picture
you get is trajectory luck; IboView lands on σ/π for the same input for the
same non-reason. The same flatness holds for σ*/π* (LUMO assignment).

### Proposed fix (implemented 2026-08-24, occupied block)

Mirror `_resolve_on_atom_mixing`, which already breaks the *same-atom* flat
subspaces (O 2s vs lone pair) by Fock-diagonalizing them. The bond case is
the two-centre analogue — `_resolve_flat_degeneracies`:

1. After PM convergence, for each occupied pair (i, j), compute the PM
   functional value of the pair now and at a 45-degree rotation.
2. If |ΔL|/|L| exceeds ~10⁻⁶ the functional distinguishes the pair; skip.
   Otherwise treat it as functionally degenerate and, only if the Fock
   coupling is real (|F_ij| > 10⁻⁶·max(|F_ii|,|F_jj|), purely relative so
   SCF dust never triggers), rotate by the minimal Jacobi angle zeroing
   F_ij. Aufbau ordering falls out: real σ/π pairs are strongly
   Fock-split (~0.25 Ha here) while banana pairs are nearly degenerate
   (~5×10⁻³ Ha). Already-diagonal pairs rotate by exactly φ = 0 and stay
   byte-identical.

Implementation notes learned the hard way:

- Column aliasing: `C_occ[:, i]` is a view; both rotation writes must read
  copies of the original columns. This bit even after studying the
  vectorization commit that fixed the same class of bug in the PM sweep.
- The virtual block must NOT get this pass yet: the SVD valence-virtual
  block retains near-null-singular-value residual columns (σ ≈ 0.01) that
  are functionally degenerate with real antibonds; rotating across them
  mixes π* with delocalized junk and scrambles the energy sort.

Verification (all on hf/cc-pVDZ unless noted):

- Inertness on idealized suite geometries: byte-identical `ibo.molden` and
  `ibos.txt` for water/methane/ethene/ammonia/benzene/zncl2, pre- vs
  post-change code (12/12 sha256 matches).
- Full pytest suite: 13/13 passed post-change.
- Real molecules at their own wB97X-D3/def2-TZVP equilibria (PubChem
  structures re-optimized per protocol): norbornene and cyclohexene give
  old-code == new-code output exactly, with clean σ frameworks, pure-p π
  HOMO and π* LUMO — the pass detects nothing to fix at true minima of
  these Cₛ/C₂ alkenes. Cyclooctatetraene (D₂d tub) was still converging at
  write time; check `avo_real` experiment notes if it later shows a fire.
- The bad ethene geometry from this investigation now yields σ + π + pure
  π* instead of bananas, with textbook charges — but no fixture was added:
  a test fixture encodes a requirement, not an incident (decision below).

The trigger is narrow: flat directions exist only when two orbitals share
identical population vectors, i.e. same-atom pairs (already handled) and
same-two-atom ½–½ pairs (this fix). Pairs with any population asymmetry sit
in steep bowls and can never fire. Expected suite behavior is numerical
inertness: benzene/water/methane/ammonia pairs are already Fock-diagonal;
the pass detects them and rotates nothing. It acts exactly where the suite
doesn't test — symmetry-broken geometries of symmetric molecules, i.e. real
sketcher input.

Risks: the 10⁻⁶ relative tolerances are the knobs (flat directions vary by
~10⁻⁸ relative, curved pairs by ~10⁻²); diffuse delocalized systems
(ferrocene-class) may have near-flat pairs that get aufbau-ordered within
their manifold — still valid PM-stationary IBOs, but a behavior change worth
eyeballing.

## Standing decisions

| Decision | Choice | Why |
|---|---|---|
| PM schedule | p=2 warmup → p=4 refine, fixed sweep order, no Cayley rotation | Cayley (IboView's RotateVectorsRandomly) measurably worsens symmetry-equivalent splits here (AGENTS.md Gotcha 20); fixed order gives best degeneracy |
| PM convergence | grad norm < 1e-12, max 2048 sweeps | Tighter than typical; keeps near-degenerate manifolds honest |
| Basis conventions | Cartesian (`puream=0`) everywhere, SCF basis user-selectable | Psi4 spherical-harmonics Molden output silently misrenders in Avogadro 2 (Gotcha 16); Cartesian matches the paper's assumption |
| Pipeline symmetry | C1 for the IBO pipeline; native point group only for `canonical.molden` | Individual IBOs need not respect irreps; canonical MOs render prettier when they do |
| Formula ordering | First occurrence, not Hill | Preserves input atom order (SO₃ reads better than O3S) |
| Closed-shell only | Reject spin ≥ 2 at entry | RHF-style double occupancy is baked into IAO/PM/analysis; UHF is a rewrite, not a flag (README "Limitations") |
| Donor–acceptor analysis | None; structurally impossible | Occupied block diagonalizes F^IAO, so its occ–vir block vanishes identically — mathematics.md §9 |
| Test fixtures encode requirements, not incidents | No fixture from the 2026-08 ethene investigation | Running an analysis on a bad geometry is not a requirement; revisit if a real molecule (e.g. a non-planar conjugated system at its own optimum) needs this code path covered |
| Structure provenance | PubChem 3D structures are untrusted until re-optimized by us (wB97X-D3/def2-TZVP via ORCA) | Unvetted PubChem geometries caused the 2026-08-24 ethene confusion; see also user protocol |

## Verified against code (2026-08-24)

- 13 pytest CLI integration tests (`pixi run test`); counts parametrized over
  water/methane/ethene/ammonia/benzene/zncl2.
- Element tables extend through iodine (Z=53) in both `calcs.py`
  (`_ELEM_SYMBOLS`) and `__main__.py` (`_ELEMENT_NUMBERS`).
- Charge decomposition and total Wiberg bond orders are implemented
  (`_analyze_ibos`, `_format_total_wiberg`) — they shipped in v0.4 despite
  being written up as "future" work in older planning docs.
- Defaults live in `calcs/config.json` via `config.py`:
  wB97X-D/def2-TZVP, `iboview_style=True`; the Avogadro Options dialog
  edits them persistently.
