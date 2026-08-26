# NOTES

Engineering notebook for avo_ibo: open items, design decisions, and the
reasoning behind them. Claims here are checked against the code; when the
code changes, update the note or delete it.

Status docs live in [AGENTS.md](AGENTS.md); derivations in
[mathematics/mathematics.md](mathematics/mathematics.md); a hands-on
plugin-development guide in [tutorial.md](tutorial.md).

## Open items

- **pixi-pack distribution**: deliberately deferred (2026-08-26).  The
  ~1 GB self-extracting environment archive isn't justified at current
  adoption; revisit only if usage grows.  The working recipe stays in
  [tutorial.md §17](tutorial.md).
- **Virtual-block junk-column hygiene** (see IboView audit below): either
  adopt their overlap-weighted energy scheme or exclude weak SVD columns
  from ordering-sensitive paths; unlocks a σ*/π* tie-break later.
- **Bond-flat degeneracy resolution**: implemented for the occupied
  block (`_resolve_flat_degeneracies`); virtual block intentionally deferred
  until the junk-column hygiene above.

## Ideas parked for the future (2026-08-26, not active work items)

Ranked roughly by value/cost; none is planned.  Feature-complete verdict
vs. Knizia's single-geometry palette: yes (see open-items discussion of
2026-08-26); the trajectory/electron-flow class is a second plugin, not an
extension.

- **σ/π Wiberg decomposition**: IMPLEMENTED 2026-08-26
  (`_format_wiberg_by_type`, section "Wiberg Bond Orders by Type" in
  ibos.txt).  Per-orbital W_AB = occ²·P_A·P_B summed by σ/π class; the
  header states it differs from the density-matrix total (cross-orbital
  D² terms), so both are shown.
- **Charge decomposition by orbital class**: split of Q_A into core/LP/
  bonding contributions per atom.  Low cost (all machinery exists).
- **d-electron count / oxidation-state readout for metals**: sum of
  d-character IBO populations on a metal centre; low cost, useful for
  the ZnCl₂-style tests.  Requires care on d/σ mixing (SO₃ d-polarization
  class of ambiguity).
- **HOMO–LUMO gap summary block**: energies already computed; a headline
  ΔE in the table footer.  Trivial cost.
- **Aromaticity index (multi-centre bond order)** e.g. 6-centre index for
  benzene: research-flavoured, medium cost; probably not worth it for a
  student audience.
- **Trajectory / electron-flow analysis** (Knizia–Klein 2015 Angew
  "Electron Flow in Reaction Mechanisms"; 2018 cPCET-vs-HAT): requires
  consuming pre-computed IRC paths (ORCA/Gaussian output or xyz series) —
  NOT running IRCs, which stays out of scope.  A second plugin class:
  per-frame SCF exists, but the Avogadro one-call-per-molecule protocol
  needs a file-reading command and per-frame rendering.  Not planned
  unless adoption or teaching demand materialises.

## Typed API (issue #3, landed 2026-08-25)

`compute_ibo_data(cjson, options, charge=0, spin=1) -> IBOResult` is the
pure core: no project files, no config reads, no logger wiring.  Options
fall back to library defaults (hf/cc-pVDZ); callers wanting persistence
merge `config.load_config()` beforehand.  `compute_ibo()` remains the
Avogadro adapter (calc dirs, psi4 logging, molden rendering, JSON
contract) and is verified byte-for-byte identical pre/post refactor on
all six suite molecules.

Psi4 global-state lessons now encoded in the core:

- Primary-output routing is mandatory on Windows — without
  `psi4.set_output_file(...)` the PSIO manager dies ("cannot get a mirror
  file handle").  Core defaults to a private temp file; adapter passes
  its calc-dir log.
- Molecule registration (`psi4.geometry` + `reset_point_group("c1")`)
  must happen before SCF inside the core — it defines AO ordering for
  everything downstream.
- Psi4 aux files (`timer.dat`) land in cwd; library callers in read-only
  directories will fail regardless of output routing.

## Bond-flat PM degeneracies: resolved for the occupied block

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
the two-centre analogue — `_resolve_flat_degeneracies` (formal derivation:
[mathematics/mathematics.md](mathematics/mathematics.md) §6):

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
| PM schedule | p=2 warmup → p=4 refine, fixed sweep order, no Cayley rotation | Determinism is a feature: IboView's `LocalizeVectors` applies an *unconditionally*, wall-clock-seeded 18° Cayley perturbation to every block it localizes (CtOrbLoc.cpp `if (1)` guard, seed = `g_LastSeedOffset++ + 1e6·Second()`), so its output on functionally-degenerate manifolds is non-reproducible in principle. Fixed order from canonical MOs gives best degeneracy here (Gotcha 20). |
| PM convergence | grad norm < 1e-12, max 2048 sweeps | Matches modern IboView exactly (`LocOpt.ThrLoc = 1e-12`, hard override at IvIao.cpp:471 of the older `min(1e-6, 1e-2·ThrGrad)`; see cgk's 2019 comment on core/lone-pair mixing under the loose value) |
| Basis conventions | Cartesian (`puream=0`) everywhere, SCF basis user-selectable | Psi4 spherical-harmonics Molden output silently misrenders in Avogadro 2 (Gotcha 16); Cartesian matches the paper's assumption |
| Pipeline symmetry | C1 for the IBO pipeline; native point group only for `canonical.molden` | Individual IBOs need not respect irreps; canonical MOs render prettier when they do |
| Formula ordering | First occurrence, not Hill | Preserves input atom order (SO₃ reads better than O3S) |
| Closed-shell only | Reject spin ≥ 2 at entry | RHF-style double occupancy is baked into IAO/PM/analysis; UHF is a rewrite, not a flag (README "Limitations") |
| Donor–acceptor analysis | None; structurally impossible | Occupied block diagonalizes F^IAO, so its occ–vir block vanishes identically — mathematics.md §10 |
| Test fixtures encode requirements, not incidents | COT tub fixture added (2026-08-24); synthetic-perturbation test rejected | The bad ethene geometry itself is not a requirement; COT is — a real non-planar conjugated π system at its own equilibrium, guarding σ/π character and tie-breaker inertness. If trigger-path coverage is ever needed, generate the perturbation programmatically at test runtime from ethene.xyz (synthetic stimulus, clearly labeled), never as a memorialized "molecule". |
| Structure provenance | PubChem 3D structures are untrusted until re-optimized by us (wB97X-D3/def2-TZVP via ORCA) | Unvetted PubChem geometries caused the 2026-08-24 ethene confusion; see also user protocol |

## IboView source audit (2026-08-24, ibo-view.20211019-RevA)

Findings from reading `IvIao.cpp` and `MicroScf/CtOrbLoc.cpp` directly,
relevant to parity claims and future work:

1. **Valence-virtual SVD cutoff is identical**: keep σ > 1e-8
   (`AllocAndComputeSvd` + loop at IvIao.cpp:111).  Our port was faithful;
   the near-null residual columns exist in IboView too.
2. **Why they don't scramble on junk columns**: virtual-block energies are
   *not* Fock diagonals.  `MakeOrbitalEnergies_General` computes
   ε_new = Σₖ |⟨refₖ|new⟩|² · ε_canonₖ — overlap-weighted averages of
   canonical MO energies over the full reference set.  Comment admits the
   scheme needs cleanup.  Consequence for us: our Fock-diagonal energies
   plus cross-block energy sorting is where junk columns become dangerous;
   any future fix should either adopt an IboView-style energy estimate or
   exclude weak columns from ordering-sensitive paths.
3. **They localize all four case blocks** through the same
   `LocalizeVectors` (RHF: occupied + valence-virtual), so our pre-existing
   virtual localization matches their structure.
4. **Random perturbation is unconditional and time-seeded**: the `if (1)`
   guard at CtOrbLoc.cpp applies an 18° Cayley rotation to every localized
   block; `FRandomNumberGenerator` seeds from wall clock (`CxRandom.h`).
   Bit-parity with IBOView on degenerate manifolds is therefore impossible
   in principle — and their own output varies run to run there.
5. **The Hessian gate handles one-center cases only**: pairs are skipped
   when |A_ij| ≤ ThrLoc; the comment cites lone-pair/same-center mixing.
   No two-center bond-flat handling exists in either implementation — our
   `_resolve_flat_degeneracies` fills that gap rather than duplicating one.

## Verified against code (updated 2026-08-26)

- 15 pytest CLI integration tests (`pixi run test`); counts parametrized over
  water/methane/ethene/ammonia/benzene/zncl2/cyclooctatetraene.
- Element tables extend through iodine (Z=53) in both `calcs.py`
  (`_ELEM_SYMBOLS`) and `__main__.py` (`_ELEMENT_NUMBERS`).
- Charge decomposition and total Wiberg bond orders are implemented
  (`_analyze_ibos`, `_format_total_wiberg`) — they shipped in v0.4 despite
  being written up as "future" work in older planning docs.
- Defaults live in `calcs/config.json` via `config.py`:
  wB97X-D/def2-TZVP, `iboview_style=True`; the Avogadro Options dialog
  edits them persistently.
