# NOTES

Engineering notebook for avo_ibo: open items, design decisions, and the
reasoning behind them. Claims here are checked against the code; when the
code changes, update the note or delete it.

Status docs live in [AGENTS.md](AGENTS.md); derivations in
[mathematics/mathematics.md](mathematics/mathematics.md); a hands-on
plugin-development guide in [tutorial.md](tutorial.md).

## Cross-validation against IboView and ORCA (2026-09-03 to 09-04; record, not gate)

Three programs computed the same chemistry: this plugin (Psi4 SCF plus
IAO/PM pipeline), Knizia's IboView, and ORCA 6.1.1 (`%loc` / `IAOIBO`).
Five molecules at wB97X-D/def2-TZVP, identical geometries in every leg.
Verdict: our implementation reproduces both programs wherever the
comparison is well-posed. The user-facing summary is
[validation/Validation.md](validation/Validation.md); inputs and result
logs live under [validation/](validation/) (`iboview/`,
`molden_bridge/`, `orca/`).

Two standing rules govern everything below. Snapshots are dated
records, never suite gates — IboView's wall-clock-seeded Cayley
rotations make cross-implementation golden tests flaky by
construction. And comparisons lead with Δ trends, where functional
systematics cancel; the EAS mechanism study itself stays out of this
repo.

### Pilot: phenol and anisole Wheland geometries

Four ORCA-optimized geometries from the Wheland study (phenol,
para-protonated phenol, anisole, para-protonated anisole;
`orca_calcs/phenol_v_anisole`, last `.trj.xyz` frames) went through the
plugin at wB97X-D/def2-TZVP and were compared against existing IboView
IBBA logs (wB97X-D3/def2-TZVP). One hard lesson first: an early run on
unoptimized inputs produced a 0.15 e⁻ phantom offset and a 7 kcal/mol
energy gap before the geometry mismatch was caught. Same geometry is a
precondition for any cross-implementation claim; inputs preserved at
`validation/phenol_v_anisole/*_opt.xyz` (untracked).

| | Ph GS (plug/IBV) | Ph WH | An GS | An WH |
|---|---|---|---|---|
| O | −0.342/−0.498 | −0.203/−0.361 | −0.265/−0.346 | −0.144/−0.227 |
| C1 (ipso) | +0.171/+0.231 | +0.321/+0.391 | +0.165/+0.220 | +0.312/+0.380 |
| C4 (para) | −0.090/−0.179 | −0.092/−0.296 | −0.089/−0.179 | −0.091/−0.295 |
| C8 (methyl) | — | — | +0.020/−0.201 | +0.010/−0.223 |
| ΔO (WH−GS) | +0.139/+0.137 | | +0.122/+0.119 | |
| ΔC1 | +0.150/+0.160 | | +0.147/+0.160 | |
| ΔC8 | | | −0.010/−0.022 | |

Wheland Wibergs agree to ±0.02 throughout: C1–O 1.312/1.291,
C2–C3 1.702/1.698, C3–C4 1.068/1.070 (phenol); C1–O 1.347/1.325,
O–C8 0.912/0.895 (anisole).

Three findings. First, trends agree while absolutes carry a systematic
X–H offset: every significant absolute gap sits on a hydrogen
(~+0.08–0.11 e⁻ more positive in IboView per X–H bond). The phenolic
O–H σ reads 68.5/31.1 vs our 62.9/36.6 (the 0.11 e⁻ accounts for the H
gap exactly); methyl C8's 0.22 gap is three hydrogens × ~0.08; the C4 Δ
"discrepancy" (−0.002 vs −0.117) is two new C–H bonds × ~0.09 — same
phenomenon, not new physics. The heavy-atom framework agrees to ~±0.05
and lone-pair compositions are identical (O 2pz LP 1.872/0.078 vs
1.876/0.084). Second, the EAS mechanistic conclusion is robust across
implementations: phenol O donation exceeds anisole's by 14% here (0.139
vs 0.122) against 15% in IboView (0.137 vs 0.119) — same verdict, same
magnitude. Third, a functional caveat: plugin side ran wB97X-D, IboView
side wB97X-D3. Psi4 shells D3 out to an `s-dftd3` binary that isn't
installed, and installing it would dirty the locked distribution
environment for fourth-decimal charges. Revisit only on reviewer
demand.

### Feeding our Molden files to IboView: two file-format traps

Before any IboView leg could run, our `canonical.molden` tripped two
traps in IboView's reader. Both were ours to fix; both are now handled
by `validation/molden_renorm.py`:

(a) CRLF line endings. `canonical.molden` is 2560/2560 CRLF, and
IboView's `is_whitespace_cxp1()` accepts only space/tab, so the CR
survives `str_trim` and line 1 mismatches (`IvOrbitalFile.cpp:1222`,
`CxParse1.cpp:180`). Convert CRLF to LF.

(b) Unnormalized contractions. Psi4's molden writer prints RAW
contraction coefficients (O def2-TZVP s-shell self-overlaps 0.1441 /
0.3881 as written; verified against S built from normalized
primitives). IboView builds overlap assuming normalized primitives —
renormalizing only ORCA-sourced files, detected via an `orca_2mkl`
`[Title]` tag we don't emit — so its S is inconsistent with our
coefficients, and its own sanity check reports rmsd 0.19 with 6.27 e⁻
instead of 10. The fix normalizes multi-primitive s/p contractions
(exact by symmetry; d/f singles already normalized) and re-emits LF.
Occupied C^T S C improves 0.71 → 4.9e-3 (150×). Permuting d/f to or
from Psi4-internal order makes it WORSE (2.6e-2), so Psi4 already
writes Molden-standard order; the residual 5e-3 lives in d/f tails, and
IboView warns-but-continues, so comparisons allow ±0.02.

One thing NOT to do: scale the MO coefficients alongside the
contractions. C^T S C is invariant under simultaneous rescaling
(verified POST == PRE to all digits) — Psi4's s/p coefficients are
already in its normalized basis; only the written contractions are raw.

### Same-density water test: the gap is the minimal basis

With the bridge fixed, our wB97X-D water MOs went through both
pipelines: O −0.494 (ours) vs −0.658 (IboView); O–H σ 62.9/36.6 vs
67.6/32.4; Wiberg 0.939 vs 0.904; lone pairs agree. Every gap sits on
H. The mechanism is confirmed, not suspected: IboView ships
`bases/minao.libmol` and builds MINAO IAOs, while we build over Psi4
STO-3G. MINAO's H 1s is a 5-primitive function (exponents 33.87→0.10)
against STO-3G's 3-primitive — a much better atomic function, hence
systematically different X–H polarization. The residual file slop
(5e-3) contributes ~0.01, not 0.16. Two caveats: IboView's
orbital-by-orbital mapping differs (its unconditional 18° Cayley kick
on imperfect input found another PM basin — the banana-bond phenomenon
observed cross-program), but total charges are rotation-invariant and
unaffected; and its virtuals are still garbage from d/f tails, so
compare occupied orbitals only.

### Decision: stay on STO-3G (2026-09-04)

ORCA's manual (§9.2.5, `orca_loc`) confirms the minimal basis is a soft
convention, not physics. ORCA's own default is SCF_SV — IAOs from
converged atomic SCF MOs, "instead of the MINI or STO-3G basis sets as
in the original method" — with `IAOBasis` offering STO_3G, MINI,
ANO_SZ, ANO_RCC_MB, and MINAO_AUTO_PP. ORCA's own default charges are
admitted to be only "very similar to the original IAO charges." If the
reference implementation's successors treat the choice as a parameter,
chasing MINAO would buy agreement with one program's default at the
cost of churning every number in this project — golden values,
examples, docs — without changing any chemistry (assignments, trends,
Δs all unaffected). ORCA's other improvement, IAOBOYS (Foster–Boys
instead of PM), changes orbital shapes, not charges: IAO-Mulliken
charges are localization-invariant, so it is irrelevant to the gap and
inapplicable to our PM-based σ/π edifice. The conformance claim is
therefore: same-convention agreement ±0.02; cross-convention, trends
and Δs.

### Third leg: ORCA reproduces us to 0.006

The decisive experiment: ORCA `LocMet IAOIBO` with `IAOBasis STO_3G` —
same minimal basis, same PM functional — must reproduce our charges.
Recipe (`validation/orca/*.inp`):

```
%loc
  LocMet IAOIBO
  IAOBasis STO_3G
  T_CORE -99.9
  PrintLevel 3
end
```

Gotcha: `%loc` excludes the core by default. The first run's
populations summed to 8.0 e⁻ and charges to +2.0 (the O 1s silently
dropped from the window); `T_CORE -99.9` fixes it, and the checksum is
post-localization charges summing to 0.000.

| Molecule | Atom | Ours | ORCA | Δ |
|---|---|---|---|---|
| Water | O / H | −0.494 / +0.247 | −0.4896 / +0.2448 | 0.004 |
| Ethene | C / H | −0.131 / +0.066 | −0.1286 / +0.0643 | 0.002 |
| Ammonia | N / H | −0.530 / +0.177 | −0.5236 / +0.1746 | 0.006 |
| Benzene | C / H | −0.065 / +0.065 | −0.0636 / +0.0636 | 0.001 |
| Carbon monoxide | C / O | +0.069 / −0.069 | +0.0679 / −0.0679 | 0.001 |

Worst Δ: 0.006 on nitrogen — grid noise between two programs, not
method error. Convention isolated; implementation exonerated.

Two companion lessons. Mayer bond orders (default SCF population
output, pairs above 0.1, no keyword needed) track our density-Wiberg
pattern but not digits — Mayer 1.9628 vs 2.033 on C–C, 0.9481 vs 0.969
on N–H, 1.4115–24 vs 1.444 on benzene C–C — because Mayer is built in
the full non-orthogonal AO basis while ours comes from the density in
IAO space; expect ~0.05 offsets on multiples, and ORCA offers no σ/π
partition. And partitioning matters ~50× more than program: ORCA's own
full-basis Mulliken gives ethene C −0.2546 against its IAO −0.1286 on
the identical density.

### IboView legs: ammonia passes, benzene exceeds the bridge

Ammonia (user-loaded bridge file): N −0.723/H +0.255 vs ours
N −0.530/H +0.177 (Δ = 0.19, direction as predicted: N–H σ 62.9/37.1
vs ours 58.8/41.2); Wiberg 0.9387 vs 0.969 (Δ = 0.03, same as water's
0.035). The log's electron leak (0.043, rmsd 1.64e-2, worst of the
hydrides) bounds the true IboView-N to [−0.77, −0.72] — N–H is
genuinely more convention-sensitive than O–H, not just leak.

Benzene voids — and the void is the finding. IboView counts 40.22 of
42 electrons (4.2% leak), biasing every atom ≈ +0.15 and flipping
carbon's sign. The leak progression (water 0.025, ammonia 0.043,
benzene 1.78: six times the f-channels, forty times the leak) fingers
the bridge's d/f-tail residual, not the pipeline — and qualitatively
the run is perfect (clean σ orbitals, delocalized four-center π
orbitals, C–H 58/42 in the MINAO direction, para C···C 0.109 vs our
0.116). Lesson: the bridge is valid for first-row hydrides (leak <
0.5%) but out of its depth at benzene-scale f-content. Benzene's
validation rests on the ORCA leg (Δ = 0.001), which never touches our
Molden writer.

### Carbon monoxide: the problem child passes

Ours (Psi4 density): C +0.069/O −0.069, σ 0.942, π 1.657, total 2.599.
ORCA (own density): C +0.0679 (Δ = 0.001 — tightest yet); Mayer 2.4549
(Mayer/Wiberg divergence grows with bond order; the pattern still reads
triple). The IboView leg ran on ORCA's OWN `.loc.molden`, which reads
cleanly with no warnings and an exact 14.000 total — confirming the
earlier residuals were our writer's contractions, and giving a pure
same-density convention comparison: MINAO C +0.1916 vs STO-3G +0.068,
gap 0.124 in the familiar direction. Best of all, IboView's own
orbitals confirm our σ/π assignment compositionally (σ-like 64.9/35.1
vs ours 62/38; π-like 72.4/27.6 vs ours 70.7/29.3). Caveats, both
contained: Mayer diverges more on triple bonds (never promised digits),
and IboView-on-molden core energies are unreliable (−16.8/−12.9 vs true
−19.3/−10.4, the same anomaly as water's −16.56) while
valence/charges/Wibergs are unaffected. (IAO charges needn't reproduce
CO's reversed dipole sign; the s-rich C(LP) HOMO is the C-end
basicity.)

### A trap for future validators: `orca_2mkl` writes canonical energies

`orca_2mkl` on a `.loc` file carries the localized vectors but leaves
the CANONICAL eigenvalues in the Molden `[MO]` block — verified
bit-for-bit against ORCA's own ORBITAL ENERGIES table. Ours are
localized expectations (`calcs.py:1536`: C_i·F_IAO·C_i). Comparing the
two is apples-to-oranges; one apparent σ-degeneracy breaking dissolved
into canonical 3a1 vs 1b1. LMO energies are therefore NOT validation
targets (basin-, seed-, and definition-dependent; ORCA doesn't print
LMO F_ii). One rigorous exception: symmetry-protected orbitals must
match, and water's HOMO does (ORCA canonical 1b1 −0.4136 vs our p-LP
−0.4011, Δ = 0.013 ≈ functional/grid). Charges and Wibergs stay the
validators — rotation-invariant, basin-independent. (Related: ORCA
`%loc` offers `Random 0` for fixed-seed testing. We already have that
property — deterministic Jacobi sweeps, no random kick — unlike
IboView's wall-clock seed.)

## Open items

- **pixi-pack distribution**: deliberately deferred (2026-08-26).  The
  ~1 GB self-extracting environment archive isn't justified at current
  adoption; revisit only if usage grows.  The working recipe stays in
  [tutorial.md §17](tutorial.md).
- **Virtual-block junk-column hygiene** (see IboView audit below): either
  adopt their overlap-weighted energy scheme or exclude weak SVD columns
  from ordering-sensitive paths; unlocks a σ*/π* tie-break later.
  THIRD OPTION (2026-09-04, from Derricotte & Evangelista JCTC 2017,
  DOI 10.1021/acs.jctc.7b00493): our SVD block already IS their VVO
  construction (overlap S between canonical virtuals and IAOs, SVD,
  rotate virtuals by U — their eqs 3–5). The difference is selection:
  they take EXACTLY N_VVO = N_IAO − N_occ columns by construction and
  never threshold, so junk columns cannot occur; we keep all σ > 1e-8,
  which admits σ≈0.01 squatters. Recommended direction: fixed-count
  VVO (1-line change + min-kept-σ warning), IboView-energy-weighting as
  fallback. Scope: σ*/π* tie-break in the virtual block only (pictures
  + Type labels at distorted geometries); occupied output must stay
  byte-identical, but virtual `ibos.txt` lines, LUMO identity, and
  example gap numbers must all be re-verified (junk may currently own
  some LUMOs). Verification mirrors the occupied tie-break: suite
  18/18, runtime-generated (never fixture) distorted-geometry
  banana→σ*/π* repair check, examples gap re-grep.
- **Bond-flat degeneracy resolution**: implemented for the occupied
  block (`_resolve_flat_degeneracies`); virtual block intentionally deferred
  until the junk-column hygiene above.

## Examples gallery (2026-09-03; tracked in-repo since 2026-09-04)

`examples/README.md` is now an index over a how-to-read guide
(`how-to-read-ibos.md`, section-by-section through the current format
on water) plus per-molecule spotlights, each with run command, verified
numbers, and orbital images in `examples/img/`:

- First-look set refreshed to wB97X-D/def2-TZVP: SO₃, water, methane,
  ethene, ammonia, benzene, ZnCl₂; methylamine writeup added.
- New: carbocations (bridged ethylium vs classical t-butyl),
  cyclopropenyl cation + planar/nonplanar anion pair (antiaromatic →
  nonaromatic, HOMO +0.029 → −0.090 Ha), cyclopropane bent bonds,
  COT tub hyperconjugation, ozone, malonaldehyde H-bond, allene
  orthogonal-π, diborane bridges.
- All 19 molecules regenerated with the current pipeline at the
  project standard; every prose number re-verified (two level-carry
  mistakes caught in drafting: always re-grep, never carry).
- Geometry provenance: small-molecule MMFF94 inputs accepted as
  adequate; fixed three mislabelled comment lines found along the way
  (ethylium "ethane", methyl cation "methane", t-butyl
  "2-methylpropan-1-ium"; anion file labelled "cyclopropenium").
  The old `methyl_amine_001/analysis.txt` describes *allene*
  (central C +0.026 — matches to the digit) and was reassigned
  accordingly. Cyclopropane/ozone/allene geometries recovered from
  prior-calculation molden files (bohr → Å).
- Suspects flagged in entries, not silently used: ozone O–O 1.11 Å
  vs 1.27 exp (multireference-sensitive — re-optimise before
  production use). Ferrocene is NOT suspect-geometry (wB97X-D
  optimised, per user correction 2026-09-03) — it is slightly
  desymmetrised from D5h, and the writeup angle is precisely that:
  near-symmetric input giving symmetry-broken IBO tables to readers
  who assume perfect symmetry. Existing IBO table on disk at
  MN15-L/def2-SVP; needs a wB97X-D/def2-TZVP rerun before any
  writeup.
- Deferred spotlights (data on disk): anisole, phenol, cyclobutadiene,
  formaldehyde.

## Ideas parked for the future (2026-08-26, not active work items)

Ranked roughly by value/cost; none is planned.  Feature-complete verdict
vs. Knizia's single-geometry palette: yes (see open-items discussion of
2026-08-26); the trajectory/electron-flow class is a second plugin, not an
extension.

- **σ/π Wiberg decomposition**: IMPLEMENTED 2026-08-26
  (`_format_wiberg`, section "Wiberg Bond Orders (σ/π, density)" in
  ibos.txt).  One consolidated table replaces both old Wiberg sections
  (density-matrix total + per-IBO by-type).  Total is the original
  density Wiberg W_AB = Σ D²_ij, decomposed EXACTLY via the orbital-pair
  expansion W_AB = 4 Σ_kl G^A_kl G^B_kl (mathematics.md §9.4): diagonal
  terms are the per-IBO shares, off-diagonal terms are interference,
  folded into σ or π by class (mixed σ/π pairs split 50/50), so
  σ + π = total identically.  A parenthesised (interference) column
  echoes the folded interference, partitioned as (σ-part, π-part):
  the σσ terms folded into σ, the ππ terms folded into π, σπ cross
  terms split 50/50.  Rows whose interference is < 5e-4 show no
  parenthetical at all (pure diagonal shares).  Orbitals classed by
  p-fraction alone (no two-centre gate): benzene's aromatic π now
  reads π = 0.444/bond, diborane's 2e3c bridges classify σ as they
  should (H is s-only).  Verified: totals bit-identical to the old
  density table; σ+π=total to ~1e-12 on all suite molecules.
- **Detail section (2026-08-26)**: "Significant orbital-pair
  interference (|term| ≥ 0.01)" lists individual (k,l) interference
  terms grouped by bond in table order, e.g. diborane's bridge×bridge
  pair at -0.0145 on every B-H leg and +0.0145 across B-B.  Threshold
  PAIR_DETAIL_THRESH = 0.01 in calcs.py: fires only on delocalization
  chemistry, silent on ordinary molecules (ethene's largest pair terms
  are ±0.0053 — correctly excluded).  Truth-in-labelling rules agreed
  2026-08-26: positive pair terms *add to* the bond order, negative
  pair terms *subtract from* it (compositional claims only — deleting
  an orbital re-solves the SCF, so no counterfactual "would weaken by
  X" language); negative ≠ antibonding (both partners may be bonding
  orbitals — the sign is the phase relationship of overlap patterns on
  that bond); pair terms are defined in the IBO frame (remixing within
  a degenerate manifold redistributes them while the bond total is
  fixed).
- **Key finding (2026-08-26): σπ interference is identically zero on
  shipped output.**  The bond-flat and on-atom Fock resolvers run
  BEFORE the Wiberg analysis and canonicalize every flat pair into its
  σ/π eigenframe, so by the time `_format_wiberg` sees the orbitals,
  different-symmetry orbitals are already decoupled: the σπ cross term
  contributes nothing.  Deliberately disabling the resolver to
  reproduce the old banana-bond state shows WHY this matters: the two
  banana orbitals are ~50/50 σ/π mixtures (p-frac ≈ 0.81 on each
  centre — below the 0.85 π threshold), so they classify as σ and the
  π column would read 0.000 while σ inflates to ~2.0 — a silently
  wrong table.  The resolver is therefore load-bearing for the
  σ/π table, not just for the orbital pictures.  Do not "simplify" it
  away.
- **Charge decomposition by orbital class**: split of Q_A into core/LP/
  bonding contributions per atom.  Low cost (all machinery exists).
- **d-electron count / oxidation-state readout for metals**: sum of
  d-character IBO populations on a metal centre; low cost, useful for
  the ZnCl₂-style tests.  Requires care on d/σ mixing (SO₃ d-polarization
  class of ambiguity).
- **HOMO–LUMO gap summary block**: IMPLEMENTED 2026-08-26 (section
  "Frontier Orbital Energies" in ibos.txt).  HOMO/LUMO selected by
  occupancy, not rank (energies_all ascending is not guaranteed to be a
  strict occ-then-vir prefix); gap in Ha, eV and kcal/mol (CODATA 2018
  HA_TO_EV / HA_TO_KCAL).  Ethene test: π-HOMO / π*-LUMO,
  gap = LUMO − HOMO.
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
