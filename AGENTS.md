# AGENTS.md — Session Context

## Goal
Build a Pixi-based Avogadro 2 Python script plugin that computes Intrinsic
Atomic Orbitals (IAOs / Knizia IBOs) via Psi4 and returns Molden data for
interactive isosurface viewing with a full IAO-basis analysis table.

## Audience

This plugin is for:
* **Students** who wish to get a visual on a system, though a chemically
  correct one;
* **Researchers** who want to perform a preliminary analysis, but one that
  is technically sound and even publishable so long as all of the limitations
  and caveats are well understood.

Serious professionals can always run their own SCF calculations and then use
[IboView](http://www.iboview.org/).

## Status

### What Works
- **Compute IBOs**: IAO/2014 construction + Pipek-Mezey localization (p=2→p=4,
  Jacobi sweeps, conv 1e-12, max 2048 iter, no symmetry breaking).
  Produces localized bond orbitals (σ, π, LP, core) for both occupied and
  **virtual** blocks. Energies from Fock-matrix diagonal elements.
- **Valence-virtual IAOs**: SVD projection of canonical virtual MOs onto IAO
  space (IboView `MakeValenceVirtuals`), keeping components with σ > 1e-8,
  then localized with same PM p=2→p=4 procedure.
- **On-atom degeneracy resolution**: After PM, `_resolve_on_atom_mixing()`
  diagonalises F_IAO within each same-atom, DOM≈1.0 subspace.  PM cannot
  separate orbitals on the same atom (it only measures atomic populations),
  so O 2s + O lone pair mix arbitrarily.  The Fock diagonalisation restores
  the aufbau ordering (s-rich lowest, p-rich highest).
- **IAO-basis Molden writer**: `_write_iao_molden()` copies Psi4's [Atoms]/[GTO]
  headers but replaces [MO] with IAO-basis orbitals, using Fock-diagonal energies.
  Produces `n_min` orbitals (minimal-basis size).
- **Analysis table** (`ibos.txt` in `calcs/last/`): per-orbital occupancy,
  energy, DOM, bond type (σ/π/LP/Core/anti*), composition (top atoms + %),
  s/p/d hybridization on dominant atom.
- **Go to files …** menu command opens `calcs/` in Explorer.
- Per-calculation directory (`calcs/last/`, avo_xtb-style): cleared before
  each `"ibo"` run; contains `psi4.log`, `ibo.molden`, `canonical.molden`,
  `ibos.txt`.
- **No donor/acceptor delocalisation analysis**: The occupied block
  diagonalises $\mathbf{F}^{\rm IAO}$ (spectral theorem ⇒ $\mathbf{F}^{\rm IAO}_{ov}=0$),
  making both overlap-based and Fock-based occ-vir analysis structurally
  impossible.  See `mathematics.md §9` for proof.
- **Standalone CLI**: `python -m avogadro_ibo <file.xyz>` computes IBOs from an XYZ
  file without Avogadro, writing results to `calcs/last/`.
- Signal discipline: all debug output before final `print(json.dumps(...))`.
- **3d transition metals**: Works with no basis change — cc-pVDZ + STO-3G
  cover Sc–Zn. Verified with ZnCl₂ (linear, D∞h, closed-shell d¹⁰).
- **Element symbol table**: `_ELEM_SYMBOLS` in `calcs.py` and
  `_ELEMENT_NUMBERS` in `__main__.py` both extend through I (Z=53).
  All common organic elements (Br, I, etc.) display correctly in the
  analysis table.
- **Test suite**: 10 pytest CLI integration tests in `tests/test_cli.py`
  (`pixi run test`). Validates counts, Molden structure, occupancy, on-atom
  resolution, 3d metal support, and per-molecule orbital patterns.

### Tested Molecules (hf/cc-pVDZ)
- **Methane** (9 IAOs): C core + 4 C-H σ (sp³-like) + 4 C-H σ*
- **Water** (7 IAOs): O core + O 2s (s-rich LP) + O 2p (pure LP) + 2 O-H σ + 2 O-H σ*
  — on-atom resolution separates O 2s (−0.79 Ha) from O 2p LP (−0.49 Ha) cleanly.
- **Ethene**: 2 C cores + 4 C-H σ + C-C σ + C-C π.
- **Ammonia**: N core + 3 N-H σ + N LP.
- **Benzene**: 6 C cores + 6 C-C σ + 6 C-H σ + 3 π (delocalized).  The PM
  functional (p=4) cannot fully resolve the 6 C-H σ p-orbital alignment.
  Small off-axis density (pointing toward adjacent carbons) is a known
  consequence of the orthogonality constraint between σ bonds sharing a
  carbon's p-subspace — both Knizia's ibo-ref.py and IboView produce
  similar solutions.  C-H energies split by ~2e⁻⁵ Ha.
- **Formaldehyde**: 2 O/C cores + 2 C-H σ + C-O σ + C-O π + 2 O LPs.
- **ZnCl₂** (2026-07-01): Zn 3d¹⁰ closed-shell, linear D∞h, hf/cc-pVDZ.
  37 IAOs (Zn: 16, Cl: 9 each).  Verifies 3d transition metal support.

### What's Next (User's Intent)
- Additional analysis features (bond order, charge decomposition, IRC tracking)
- Multiple method/basis support

## Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Build backend | `uv_build>=0.10.2,<0.11.0` | Avogadro-bundled pixi v0.66.0 doesn't support hatchling |
| Lock file version | v6 | Bundled pixi v0.66.0 can't read v7 |
| Package install | `pixi install` + `pip install -e .` (manual, not `[tool.pixi.pypi-dependencies]`) | pypi-dependencies requires lock v7 |
| Plugin discovery | Symlink in `%LOCALAPPDATA%\OpenChemistry\Avogadro\plugins\` | Avogadro scans this directory |
| Psi4 integration | In-process (`import psi4` in compute function) | Simpler, faster, better error handling than subprocess |
| SCF basis (default) | cc-pVDZ, `puream=0` (Cartesian) | Better wavefunction than def2-SVP for IBO isosurfaces |
| Minimal basis for IAO | STO-3G, `puream=0` (Cartesian) | MINAO unavailable in Psi4; STO-3G adequate |
| IAO construction | IAO/2014 algorithm (IboView's `MakeIaoBasisNew`) | Full depolarized + repolarized IAO build, not simplified |
| Localization functional | Pipek-Mezey p=2 warmup + p=4 refinement | p=2 avoids local minima; p=4 sharpens bond-direction p-vector alignment |
| Localization convergence | Gradient norm < 1e-12, max 2048 sweeps | Matches IboView default |
| Symmetry breaking | None (canonical MOs used directly) | Fixed sweep order + p=4 refinement converges nearest symmetric local max; Cayley perturbation (18°) tested but INCREASED split — see Gotcha 20 |
| On-atom degeneracy fix | Post-PM Fock diagonalization (`_resolve_on_atom_mixing`) | PM cannot separate same-atom orbitals; F_IAO eigenvalues restore aufbau |
| Virtual IAOs via SVD | `C_IAO^T @ S @ C_vir` → SVD → keep σ > 1e-8 | Matches IboView `MakeValenceVirtuals` |
| Virtual localization | Apply same PM p=2→p=4 to virtual block | IboView localizes ALL iCase blocks |
| Orbital energies | ε_i = C_i^T @ F_AO @ C_i (diagonal of F_IAO) | Matches IboView `MakeOrbitalEnergies_General` |
| Molden output | Custom [MO] section with IAO-basis orbitals | Psi4 writes canonical MOs; we replace with n_min IAO orbitals |
| Analysis table destination | Written to `calcs/last/ibos.txt`, message points to file | Avoids blocking popup on every run |
| Stderr discipline | Absolutely no `{` or `[` on stderr before JSON output | `stripLeadingNonJson` finds first `{`/`[` in merged channel |
| `_prepare_calc_dir()` call | Only inside `"ibo"` branch | "Go to files" shouldn't wipe previous results |
| Test isolation | Tests call CLI as black box, check `calcs/last/` | No test-specific production code paths; CLI, plugin, and tests are orthogonal |

## Gotchas Hit

1. **Build backend**: hatchling → Avogadro's pixi can't install. Fix: use `uv_build`.
2. **Lock file v7**: User pixi creates v7 locks, bundled pixi v0.66.0 reads only v6.
3. **Missing entry point**: `pip install -e .` needed to create `.exe` shim.
4. **`logger.debug()` after `print(json.dumps(...))`**: Trailing stderr after JSON
   breaks `QJsonDocument::fromJson`. Fix: all debug before final `print()`.
5. **`logger.debug(f"Input keys: {list(data.keys())}")`**: Produces `[...]` on stderr
   → `[` found by `stripLeadingNonJson` before JSON `{`. Fix: use `', '.join()`.
6. **Psi4 INFO logs on stderr**: Contain `{` chars. Fix: redirect to file,
   `propagate=False`.
7. **CJSON coords as dict**: Sometimes `{"3d": [...]}` instead of flat array.
8. **UTF-8 on Windows**: Avogadro expects UTF-8 stdout.
9. **Wrong plugin directory**: Avogadro scans `%LOCALAPPDATA%`, not `%APPDATA%`.
10. **Virtual block: PM localization needed**: SVD-projected virtuals are
    delocalized without localization. Fix: call `_localize_ibos(U_val, ...)`.
11. **`n_occ` typo**: Was used instead of `nocc`. Fix: use `nocc` (defined
    as `wfn.doccpi()[0] + wfn.soccpi()[0]`).
12. **PM cannot separate same-atom orbitals**: PM functional uses atomic
    populations n_A(i) only; two orbitals with DOM≈1 on the same atom
    (e.g. O 2s and O lone pair) are degenerate in the PM functional —
    any rotation within the subspace gives identical L.  Fix: post-PM
    Fock diagonalisation within each same-atom DOM≈1 subspace.
13. **Nondeterministic on-atom mixing**: Without the Fock fix, water/def2-SVP
    produces different O 2s/O 2p mixtures each run (random QR seed).  The
    Fock diagonalisation makes the output deterministic.
14. **Canonical MO delocalisation analysis removed (2026-06-30)**:  The
    occupied block $\mathbf{C}^{\rm IAO,occ}$ diagonalises
    $\mathbf{F}^{\rm IAO}$ (it is a unitary rotation of the canonical
    occupied MOs).  By the spectral theorem, the off-diagonal occ-vir
    block of $\mathbf{F}^{\rm IAO}$ is identically zero.  This kills
    **both** overlap-based and Fock-based donor/acceptor analysis
    (including NBO-style E2).  It is not a numerical artifact — it is a
    mathematical identity proven in `mathematics.md §9`.  The plugin now
    focuses on the core IBO/Molden pipeline; `canonical.molden` is kept
    as a reference for Avogadro's MO surface dialog.  The attempted
    approaches (atom-set density, anti* pairing, within-subspace top-2)
    were each too noisy or chemically misleading for a focused plugin.
15. **IAO vs NBO: exactness forbids E(2)**: NBO gets non-zero E(2)
    because Lewis NBOs are NOT exact eigenvectors of F (occupancies
    1.90–1.99, not exactly 2.0) — the leakage creates genuine off-diag
    Fock elements.  IAO/IBO was designed for exact, lossless occupied
    representation, so F_ov = 0 by theorem.  No bug, no missing trick —
    a feature of the IAO design goal.  See `mathematics.md §9.5`.
16. **Psi4 Molden output + spherical harmonics silently wrong in
    Avogadro 2**: Psi4's `psi4.molden()` with default `puream` (spherical
    5D/7F) produces files where the [GTO] shell layout and [MO]
    coefficient ordering misalign under Avogadro 2's parser.  Orbitals
    appear to load but are chemically incorrect.  Using Cartesian
    (`puream=0`) fixes it.  Ethene B3LYP/cc-pVDZ with `puream=1` was
    the reproduction case.  See `tutorial.md §13`.
17. **[MO] padding to match [GTO] slot count (2026-06-30)**: Avogadro
    allocates one MO slot per basis function in [GTO].  If [MO] has fewer
    entries, the extra slots show uninitialised noise.  Fix: pad [MO] with
    zero-energy dummy orbitals up to n_AO total entries.  (The real fix
    would be in Avogadro's Molden reader — count [MO] entries, not [GTO]
    functions — but the padding approach is robust for any Molden reader.)
    Filed as [OpenChemistry/avogadrolibs #2890](https://github.com/OpenChemistry/avogadrolibs/issues/2890).
18. **`$$` display math broken by GitHub Markdown preprocessor (2026-06-30)**: GitHub uses MathJax v3 but runs a Markdown preprocessor
    that [mangles LaTeX inside `$$...$$` blocks](https://stackoverflow.com/a/77726873) (per MathJax core dev Davide Cervone).
    The `\Vert\nabla L\Vert = \frac{1}{n_{\mathrm{occ}}} \sqrt{\sum_{i<j} \bigl(p \,\phi_{ij}\, B_{ij}^{(p)}\bigr)^2}`
    equation at mathematics.md §4.5 triggered a "Missing close brace" error on GitHub because the preprocessor corrupted
    braces around `\mathrm{occ}`, `\sum_{i<j}`, `\phi_{ij}`, or `B_{ij}^{(p)}`.  Initial fix: use ` ```math ```` code block
    syntax instead of `$$...$$` — the code block bypasses the Markdown preprocessor entirely.
    **Update 2026-06-30**: Even with ` ```math ```` code block, the equation STILL failed on GitHub — the bug is
    in MathJax v3 itself (or GitHub's MathJax configuration), not the Markdown preprocessor.
    **Ultimate fix**: use `\lVert`/`\rVert` instead of `\Vert`, and `\lt` instead of `<` inside `\sum_{i\lt j}`,
    still inside a ` ```math ```` code block.  Both changes together resolved it; either alone may suffice but
    both are applied for robustness.  See commit `850b117`.
19. **Geometric p-orbital alignment fix attempted and failed (2026-07-01)**:  
    Rotating each σ bond's p-vector on the dominant atom to point along the
    bond direction should, in theory, make symmetry-equivalent bonds degenerate.
    In practice, F_IAO diagonal elements for px/py/pz differ per atom in the
    molecular frame because the Fock matrix is not isotropic — rotating to
    different bond directions on different atoms samples different diagonal
    elements, creating splits of ~0.9 Ha.  This approach cannot work without
    also accounting for the atomic-position-dependent Fock structure.
    The aufbau diagonalisation (`_resolve_on_atom_mixing`) remains the only
    correct fix for same-atom degeneracy (DOM≈1 subspaces).  The 2e⁻⁵ Ha
    C-H σ split in benzene is a known, chemically negligible limitation of
    the PM functional (cf. Knizia's ibo-ref.py and IboView — neither
    addresses it).
20. **Cayley random rotation (IboView) did not fix benzene asymmetry
    (2026-07-01)**:  Adding a pre-PM Cayley random rotation (18°) to match
    IboView's `RotateVectorsRandomly` actually INCREASED the C-H σ energy
    split from 2e⁻⁵ Ha to 1.2e⁻⁴ Ha.  IboView uses the same fixed sequential
    sweep order as our `_localize_ibos()`; the symmetric results the user
    sees in IboView may stem from a different Psi4 version, basis set, or
    SCF convergence — not from a different PM algorithm.  The smallest split
    (2e⁻⁵ Ha) is achieved with the original p=2 warmup + p=4 refinement, no
    symmetry breaking, and no Cayley perturbation.  IboView library defaults
    to p=4 for the same reason: p=4's stronger penalisation of mixed atomic
    populations subtly sharpens p-vector alignment along bond directions.
21. **Cartesian d-function ordering + normalization in Molden writer
    (2026-07-01)**: Psi4's internal d-function ordering (xx, xy, xz, yy, yz, zz)
    differs from Molden standard (xx, yy, zz, xy, xz, yz).  Additionally, Psi4's
    CCA convention omits the angular normalization factor for Cartesian Gaussians
    with multiple non-zero exponents, giving off-diagonal d-components (xy, xz, yz)
    a self-overlap of 1/3 vs 1 for diagonal components.  Our `_write_iao_molden()`
    copied the [GTO] header from Psi4 (implicitly Molden-standard) but wrote [MO]
    coefficients in Psi4 internal ordering with CCA normalization, causing Avogadro
    to misread d-contributions.  Fix: apply `D_PERM = [0, 3, 5, 1, 2, 4]` (ordering)
    and `D_NORM = [1, 1, 1, 1/√3, 1/√3, 1/√3]` (normalization) to each d-shell
    when writing coefficients.  Identified by colleague via C-H σ d-function
    fingerprint.  Energies (from F_IAO diagonalisation) unaffected — purely a
    Molden writer bug.  See `calcs.py:_write_iao_molden()`.
22. **Cartesian f-function reordering is forward-compat only (2026-07-01)**:
    `F_PERM = [0, 6, 9, 3, 1, 2, 5, 8, 7, 4]` and `F_NORM` (1/√5 for
    two-off-diagonal f, 1/√15 for xyz) were added to `_write_iao_molden()` at
    the same time as the d-function fix.  Not tested — cc-pVDZ has no f
    functions — but provided so that cc-pVTZ and higher bases work without a
    silent corruption bug.
23. **`_ELEMENT_NUMBERS` in `__main__.py` truncated (2026-07-01)**: The CLI
    parser (`python -m avogadro_ibo <file.xyz>) used a separate element-number
    lookup that only went to Ar (Z=18).  Any element beyond Ar (Zn, Br, I, etc.)
    silently mapped to Z=0 (ghost atom), producing wrong coordinates for
    molecules containing those elements.  Avogadro's cjson path (stdin/stdout
    protocol) was unaffected — it sends atomic numbers directly.  Fix: extend
    `_ELEMENT_NUMBERS` in `__main__.py` to match `_ELEM_SYMBOLS` in `calcs.py`,
    covering all elements through I (Z=53).  See `__main__.py:7`.
24. **D₂h symmetry splitting in diborane (2026-07-06)**: The 2e3c bridge bonds
    and terminal B-H σ bonds in D₂h diborane show small energy splittings
    (4×10⁻⁵ Ha for bridges, 7×10⁻⁵ Ha for B-H σ).  Same class of limitation
    as the benzene C-H σ split (Gotcha 19/20): the PM functional with fixed
    sequential sweep order converges to the nearest local maximum, which for
    symmetric molecules with overlapping bond subspaces is not guaranteed to
    be the perfectly symmetric solution.  The orthogonality constraint between
    σ bonds sharing a boron's p-subspace drives the splitting.  No fix known —
    the aufbau diagonalisation (`_resolve_on_atom_mixing`) addresses only
    same-atom DOM≈1 subspaces, not this orthogonal-bond coupling.

## Relevant Files

- `C:\Users\mccan\Documents\Code\avo_ibo\pyproject.toml` — 2 menu-commands (`ibo`, `open`)
- `C:\Users\mccan\Documents\Code\avo_ibo\pixi.toml` — Pixi environment
- `C:\Users\mccan\Documents\Code\avo_ibo\src\avogadro_ibo\__init__.py` — CLI entry point,
  stdin/stdout JSON I/O, dispatching for `ibo`/`open`, `_prepare_calc_dir()`
- `C:\Users\mccan\Documents\Code\avo_ibo\src\avogadro_ibo\calcs.py` — All IAO logic:
  - `_get_basis_maps()` — per-function atom/AM mapping from Psi4 BasisSet
  - `_build_iao_basis()` — IAO/2014 construction (Appendix C)
  - `_localize_ibos()` — PM Jacobi sweep with p=2→p=4 (eq 4, Appendix D)
  - `_resolve_on_atom_mixing()` — post-PM Fock diag within same-atom DOM≈1 subspaces
  - `_analyze_ibos()` — IBO table: occupancy, energy, DOM, type, composition, hybrid
  - `_write_iao_molden()` — Molden writer using IAO-basis orbitals (Fock energies)
  - `compute_ibo()` — top-level: Psi4 run, IAO build, occ+vir PM localize, on-atom
    Fock resolve, Fock energies, Molden write, analysis table
- `C:\Users\mccan\Documents\Code\avo_ibo\src\avogadro_ibo\links.py` — `open_calcs_dir()` (Explorer)
- `C:\Users\mccan\Documents\Code\avo_ibo\src\avogadro_ibo\__main__.py` — standalone CLI entry (`python -m avogadro_ibo <file.xyz>`)
- `C:\Users\mccan\Documents\Code\avo_ibo\mathematics.md` — full mathematical derivation with 2013 paper references
- `C:\Users\mccan\Documents\Code\avo_ibo\tutorial.md` — Avogadro 2 plugin architecture and gotchas
- `C:\Users\mccan\Documents\Code\avo_ibo\calcs\last\` — Runtime artifacts
  (psi4.log, ibo.molden, canonical.molden, ibos.txt)
- `C:\Users\mccan\OneDrive\Documents\Code\avo_ibo\tests\test_cli.py` — 10 CLI integration tests
- `C:\Users\mccan\OneDrive\Documents\Code\avo_ibo\tests\files\` — XYZ reference inputs
- `C:\Users\mccan\AppData\Local\OpenChemistry\Avogadro\plugins\avo_ibo` → symlink

### External Paths
- **Avogadro 2**: `C:\Program Files\Avogadro2\bin\avogadro2.exe`
- **Bundled pixi**: `C:\Program Files\Avogadro2\bin\pixi.exe` (v0.66.0)
- **User pixi**: `C:\Users\mccan\.pixi\bin\pixi.exe` (v0.70.2)
- **IboView source**: `C:\Users\mccan\OneDrive\Downloads\ibo-view.20211019-RevA\`
- **Reference plugin**: `%LOCALAPPDATA%\OpenChemistry\Avogadro\plugins\avo_xtb\`

## Avogadro Async Pipeline (Key Insight)

```
QProcess (MergedChannels) → finished signal
  → PythonScript::processFinished()
    → PythonScript::finished signal
      → InterfaceScript::commandFinished()
        → InterfaceScript::finished signal
          → Command::processFinished()
            → InterfaceScript::processCommand()
              → reads QProcess::readAll()
              → stripLeadingNonJson()
              → QJsonDocument::fromJson()
              → asyncResponse()
```

`stripLeadingNonJson` finds first `{` or `[` in merged stdout+stderr buffer.
Any stderr containing these characters before the JSON on stdout will break
parsing. This is the single most critical constraint.
