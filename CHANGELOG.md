# Changelog

All notable changes to avo_ibo are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); version numbers
track `pyproject.toml` (`avogadro-ibo`) and `CITATION.cff`.

## [Unreleased]

### Added
- Bond-flat PM degeneracy resolution (`_resolve_flat_degeneracies`,
  occupied block): distorted geometries of symmetric molecules now yield
  σ+π deterministically instead of trajectory-dependent banana bonds.
  Formal derivation in `mathematics/mathematics.md` §6; engineering
  record in `NOTES.md`.
- Cyclooctatetraene (D₂d tub) test fixture with π-system tests, guarding
  σ/π character and tie-breaker inertness at a real non-planar
  equilibrium. IboView source audit recorded in `NOTES.md`.
- Typed API (closes #3): `compute_ibo_data()` pure core returning the
  `IBOResult` dataclass; `compute_ibo()` kept as a thin Avogadro adapter
  with byte-for-byte identical outputs (24/24 sha256 across six
  molecules).
- Consolidated Wiberg table ("Wiberg Bond Orders (σ/π, density)"): the
  original density Wiberg decomposed exactly via the orbital-pair
  expansion (new `mathematics.md` §9.4), so σ + π = total identically;
  folded interference echoed parenthesised as (σ-part, π-part).
  Replaces both previous Wiberg sections.
- "Significant orbital-pair interference" detail section (|term| ≥ 0.01,
  `PAIR_DETAIL_THRESH`): individual (k,l) interference terms grouped by
  bond with table-cross-referencing orbital labels; omitted entirely on
  quiet molecules.
- "Frontier Orbital Energies" section in `ibos.txt`: HOMO/LUMO selected
  by occupancy (not rank), gap reported in Ha, eV, and kcal/mol (CODATA
  2018 constants `HA_TO_EV`, `HA_TO_KCAL`).

### Changed
- `mathematics.md`: new §6 (bond-flat derivation); §§6–9 renumbered to
  §§7–10; all `calcs.py:<line>` anchors replaced with stable
  function-name anchors.
- pixi-pack self-extracting distribution deferred pending adoption
  (closed #7): the ~1 GB artifact isn't justified at current usage; the
  working recipe stays in `tutorial.md` §17.
- Vectorized Pipek–Mezey localization sweeps (same results, faster).
- Charge/multiplicity removed from persistent config options; dead
  config loading removed from the CLI entry point.

### Fixed
- `max(d, key=d.get)` → `max(d, key=lambda k: d[k])` at two sites:
  identical semantics, silences Pylance/ty overload complaints
  (`dict.get`'s `value | None` stub return).
- Cross-document section references repaired after the mathematics
  renumbering; `[MO]`-section counting in tests hardened via regex.

## [0.4.0] - 2026-07-22

First plugin-distribution release: repo stripped to shipping essentials
(internal docs untracked from the distribution tree), distribution
renamed to `avogadro-ibo`, v6 `pixi.lock` hand-maintained for Avogadro's
bundled pixi v0.66.0. Charge decomposition and total Wiberg bond orders
shipped in the analysis table; IboView-style rendering made default;
`input.xyz` saved to each calc directory. Registered in the Avogadro
Plugin Index (closed #9).

## [0.3.0] - 2026-07-07

IboView-style rendering via a dummy STO-3G SCF; Cartesian d/f Molden
ordering and normalization corrected to Psi4's internal convention;
p=2 Pipek–Mezey angle formulas fixed; degenerate-manifold footnote
(`DEG_THRESH = 2e-4 Ha`); orbital classification extracted into
`_classify_orbital`. Closed #2, #6.

## [0.2.0] - 2026-07-06

Analysis-table improvements (heavy-atom hybrid fallback, 2e3c labels,
canonical symmetry labels); math comments standardised to Unicode
notation; install documentation reworked. Closed #1.

## [0.1.0] - 2026-07-04

Initial release: IAO/2014 construction (Knizia JCTC 2013) via in-process
Psi4, Pipek–Mezey localization (p=2 warmup → p=4 refine), Molden export
for Avogadro isosurfaces, `ibos.txt` analysis table (occupancy, energy,
DOM, bond type, composition, hybridization), Avogadro menu commands
(Compute IBOs, Options, Go to files) plus standalone CLI. Closed #4, #5.

[Unreleased]: https://github.com/exergonic/avo_ibo/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/exergonic/avo_ibo/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/exergonic/avo_ibo/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/exergonic/avo_ibo/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/exergonic/avo_ibo/releases/tag/v0.1.0
