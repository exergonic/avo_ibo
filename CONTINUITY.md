# Continuity File — Cold-Start Context

Read this first when resuming work on `avo_ibo`.  It bridges the
high-level notes in `AGENTS.md` with the implementation details
needed to make changes and run tests without rediscovery.

Related files:
- `AGENTS.md` — goal, status, key decisions, gotchas
- `MATHEMATICS.md` — full mathematical derivation with paper references
- `tutorial.md` — Avogadro 2 plugin architecture and gotchas (generic)

---

## 1. Project Card

| Key | Value |
|-----|-------|
| **Root** | `<repo_root>` (was `C:\Users\mccan\Documents\Code\avogadroapp`, now `avo_ibo/`) |
| **Package name** | `avogadro_ibo` |
| **Entry point** | `avogadro_ibo:main` → shim `avogadro-ibo.exe` |
| **Build** | `uv_build>=0.10.2,<0.11.0` (not hatchling) |
| **Env** | Pixi via `[tool.pixi]` in `pyproject.toml` (no separate `pixi.toml`) |
| **Lock file** | v6 (Avogadro bundled pixi v0.66.0 cannot read v7) |
| **Install** | `pixi run pip install -e .` (creates `.exe` shim) |
| **Plugin dir** | `%LOCALAPPDATA%\OpenChemistry\Avogadro\plugins\avo_ibo` → symlink to root |
| **Psi4** | In-process `import psi4` (not subprocess) |
| **Default SCF** | HF/cc-pVDZ, Cartesian (`puream=0`) |
| **Minimal basis** | STO-3G, Cartesian (`puream=0`) |

---

## 2. File Map

### `pyproject.toml` (43 lines)
Node config.  Two menu-commands: `ibo` (Compute IBOs) and `open` (Go to files …).
Modern Avogadro format with `[tool.pixi.workspace]`, `[tool.pixi.dependencies]`,
and `[tool.avogadro]` sections.  No `[tool.pixi.pypi-dependencies]` (would
require lock v7).

### `src/avogadro_ibo/__init__.py` (78 lines)
**CLI entry point.**  Reads JSON from stdin, dispatches by `args.feature`.

| Constant | Value |
|----------|-------|
| `PLUGIN_DIR` | `Path(__file__).resolve().parent.parent.parent` = repo root |
| `CALCS_DIR` | `PLUGIN_DIR / "calcs"` |
| `TEMP_DIR` | `CALCS_DIR / "last"` |

| Function | Lines | Called by | Signature |
|----------|-------|-----------|-----------|
| `_prepare_calc_dir()` | 19–26 | `main()` when `feature="ibo"` | `() -> None` |
| `main()` | 29–78 | `[project.scripts]` entry | `() -> None` |

**Critical discipline**: `print(json.dumps(result))` is the **last** stdout
action.  All debug goes to stderr *before* it.  Avogadro merges channels
and `stripLeadingNonJson()` finds the first `{` or `[` — any stderr
containing those chars before JSON breaks parsing.

### `src/avogadro_ibo/calcs.py` (595 lines)
**All IAO/IBO logic.**  One public function: `compute_ibo()`.

| Function | Lines | Purpose |
|----------|-------|---------|
| `_get_basis_maps(basis)` | 19–34 | Per-function atom index + AM from a Psi4 BasisSet |
| `_build_iao_basis(S, S12, S_min, C_occ)` | 41–99 | IAO/2014 construction (Appendix C) |
| `_localize_ibos(C_occ, atom_of, ...)` | 106–208 | PM Jacobi sweep p=2→p=4 warm-start (Eq 4, App D) |
| `_resolve_on_atom_mixing(C_occ, atom_of, F_IAO, ...)` | 215–254 | Post-PM Fock diag within same-atom DOM≈1 subspaces |
| `_analyze_ibos(C_IAO_all, occ_all, ...)` | 268–375 | Build formatted ibos.txt table |
| `_elem_symbol(Z)` | 378–382 | Z → symbol lookup (Z ≤ 18) |
| `_s_char(c, am_of, atom, atom_of)` | 385–389 | s-fraction on given atom |
| `_p_char(c, am_of, atom, atom_of)` | 391–395 | p-fraction on given atom |
| `_d_char(c, am_of, atom, atom_of)` | 397–400 | d-fraction on given atom |
| `_write_iao_molden(path, wfn, C_AO, occ, energies, n_orb)` | 407–439 | Replace [MO] with IAO-basis orbitals |
| `_option(options, key, default)` | 446–450 | Safe dict lookup for options |
| `compute_ibo(cjson, options, charge, spin, debug=False)` | 453–595 | **Top-level**: Psi4 → IAO build → PM → Fock resolve → SVD vir → Molden → analysis |

**Imports**: `numpy`, `scipy.linalg`, `psi4`, `logging`, `tempfile`,
`pathlib.Path`.

### `src/avogadro_ibo/links.py` (23 lines)
| Function | Lines | Purpose |
|----------|-------|---------|
| `open_calcs_dir(cjson)` | 10–23 | Opens `CALCS_DIR` in platform file manager |

---

## 3. Data Flow

```
CJSON geom
    │
    ▼
Psi4 HF/cc-pVDZ, puream=0
    │  Ca  (n_AO × n_AO)
    │  Fa  (n_AO × n_AO)
    │  S   (n_AO × n_AO)  ← mints.ao_overlap()
    │  S12 (n_AO × n_min) ← mints.ao_overlap(SCF, STO-3G)
    │  S_min (n_min × n_min) ← mints.ao_overlap(STO-3G, STO-3G)
    │  nocc = doccpi[0] + soccpi[0]
    ▼
C_occ = Ca[:, :nocc]       (n_AO × n_occ)
    │
    ▼
C_IAO, C_IAO_occ = _build_iao_basis(S, S12, S_min, C_occ)
    │  C_IAO     (n_AO × n_min)  ← orthonormal IAOs in AO basis
    │  C_IAO_occ (n_min × n_occ) ← occ MOs in IAO basis
    ▼
_localize_ibos(C_IAO_occ, atom_of)     ← modifies in-place
    │  random QR break symmetry
    │  p=2 warm-start, p=4 refine
    ▼
F_AO = Fa().np               (n_AO × n_AO)
F_IAO = C_IAO^T @ F_AO @ C_IAO  (n_min × n_min)
    │
    ▼
_resolve_on_atom_mixing(C_IAO_occ, atom_of, F_IAO)   ← modifies in-place
    │  groups same-atom, DOM>0.99 orbitals
    │  diagonalises F_block, rotates block
    ▼
occ_energies[i] = C_IAO_occ[:,i]^T @ F_IAO @ C_IAO_occ[:,i]
    │
    ▼
C_vir = Ca[:, nocc:]         (n_AO × n_vir)
SIbVir = C_IAO^T @ S @ C_vir  (n_min × n_vir)
U_svd, Sigma, _ = svd(SIbVir)
U_val = U_svd[:, Sigma > 1e-8]  (n_min × n_val_vir)
if n_val_vir > 1: _localize_ibos(U_val, atom_of)
vir_energies[i] = U_val[:,i]^T @ F_IAO @ U_val[:,i]
    │
    ▼
C_IAO_all = hstack([C_IAO_occ, U_val])   (n_min × n_orb)
energies_all = concat([occ, vir])
order = argsort(energies_all)
C_IAO_all = C_IAO_all[:, order]          ← sorted by energy
    │
    ▼
C_AO_all = C_IAO @ C_IAO_all   (n_AO × n_orb)   ← back to AO basis
_write_iao_molden(path, wfn, C_AO_all, occ_all, energies_all, n_orb)
_analyze_ibos(C_IAO_all, ...) → ibos.txt
    │
    ▼
JSON response to Avogadro: { molden, cjson, message, readProperties, moleculeFormat }
```

### Array shapes for water (STO-3G IAO, cc-pVDZ SCF)

| Object | Shape | Notes |
|--------|-------|-------|
| AO basis | n_AO = 24 | cc-pVDZ: O(9s5p1d) + H(4s1p) × 2 |
| Minimal basis | n_min = 7 | STO-3G: O(1s,2s,2p×3) + H(1s) × 2 |
| n_occ | 5 | O: 1s²2s²2p⁴ → 5 doubly-occupied MOs |
| n_vir | 19 | 24 − 5 |
| n_val_vir | ~2 | SVD σ > 1e-8 |
| n_orb out | 7 | 5 occ + 2 val-vir |

---

## 4. Test Suite

### Prerequisites

```powershell
pixi run pip install -e .
# Creates .pixi\envs\default\Scripts\avogadro-ibo.exe
```

### Test all molecules

```powershell
cd C:\Users\mccan\Documents\Code\avo_ibo

pixi run python -c "
import json, sys
sys.path.insert(0, 'src')

molecules = {
    'water': {
        'coords': [0.0, 0.0, 0.117, 0.0, 0.757, -0.469, 0.0, -0.757, -0.469],
        'numbers': [8, 1, 1],
    },
    'methane': {
        'coords': [0.0,0.0,0.0, 0.0,0.0,1.089, 0.0,1.089,0.0, 1.089,0.0,0.0, -0.629,-0.629,-0.629],
        'numbers': [6, 1, 1, 1, 1],
    },
    'ethene': {
        'coords': [0.0,0.0,0.662, 0.0,0.0,-0.662, 0.0,0.924,1.238, 0.0,-0.924,1.238, 0.0,0.924,-1.238, 0.0,-0.924,-1.238],
        'numbers': [6, 6, 1, 1, 1, 1],
    },
    'ammonia': {
        'coords': [0.0,0.0,0.116, 0.0,0.941,-0.272, 0.815,-0.471,-0.272, -0.815,-0.471,-0.272],
        'numbers': [7, 1, 1, 1],
    },
    'formaldehyde': {
        'coords': [0.0,0.0,-0.537, 0.0,0.0,0.693, 0.0,0.937,-1.137, 0.0,-0.937,-1.137],
        'numbers': [6, 8, 1, 1],
    },
}

from avogadro_ibo.calcs import compute_ibo

for name, mol in molecules.items():
    cjson = {'atoms': {'coords': {'3d': mol['coords']}, 'elements': {'number': mol['numbers']}},
             'properties': {'totalCharge': 0, 'totalSpinMultiplicity': 1}}
    result = compute_ibo(cjson, {}, charge=0, spin=1)
    molden = result.get('molden', '')
    print(f'{name}: {len(molden)}b molden, msg={result.get(\"message\",\"\")}')" 2>&1 | Select-String -NotMatch 'WARN|^20'
```

### Expected IBO patterns

| Molecule | n_min | Occupied IBOs | Virtual IBOs |
|----------|-------|---------------|--------------|
| Water | 7 | O(Core), O 2s (s-rich LP), O 2p (pure LP), 2× O-H σ | 2× O-H σ\* |
| Methane | 9 | C(Core), 4× C-H σ (sp³) | 4× C-H σ\* |
| Ethene | 14 | 2× C(Core), 4× C-H σ, C-C σ, C-C π | C-C π\*, 4× C-H σ\*, C-C σ\* |
| Ammonia | 8 | N(Core), 3× N-H σ, N(LP) sp³ | 3× N-H σ\* |
| Formaldehyde | 12 | O(Core), C(Core), C-O σ, 2× C-H σ, C-O π, O(LP s-rich), O(LP pure p) | C-O π\*, 2× C-H σ\*, C-O σ\* |

### Verification checklist

After any change:
1. Test water (smallest, fastest) — check O 2s vs O 2p LP are separate
2. Test methane — check 4 degenerate C-H σ bonds
3. Verify analysis table (`calcs/last/ibos.txt`) DOM values and hybridisation
4. Verify Molden output has correct `[MO]` section with `Ene=` lines
5. Run in Avogadro: **Extensions → Intrinsic Bond Orbitals → Compute IBOs**

---

## 5. Development Survival

### Edit → Test → Verify loop

```powershell
pixi run pip install -e .    # after any code change
pixi run python -c "..."     # inline test (see §4)
```

No need to reinstall pixi environment unless dependencies change.

### Debug tips

- **Nondeterministic output?** The PM random QR seed makes IBOs vary across
  runs (same-atom mixing).  The Fock diagonalisation fix
  (`_resolve_on_atom_mixing`) removes this *for existing molecules*, but
  verify on any new molecule type.
- **Psi4 crash?** Check `calcs/last/psi4.log` for the full output.
- **Avogadro not showing IBOs?** The Molden output must have `Sym= A` for every
  orbital and correct `Ene=` / `Occup=` lines.  Check via CLI test first.
- **Stderr / JSON parsing failure?** Look for `{` or `[` in stderr before
  `print(json.dumps(...))`.  Use the psi4 log redirection pattern.

### Pitfalls

1. **Lock file v6**: If you run `pixi install` with a new pixi (≥0.70), it
   upgrades to v7 and Avogadro's bundled pixi (v0.66.0) can't read it.
   Use the bundled pixi to generate the lock: `& "C:\Program
   Files\Avogadro2\bin\pixi.exe" install`
2. **`puream` mismatch**: The IAO construction assumes Cartesian functions
   (`puream=0`).  Changing this will break the orbital counting in
   `_get_basis_maps` and the IAO build.
3. **Non-linear dependencies**: `scipy.linalg.cho_factor` and `.cho_solve`
   are used in `_build_iao_basis`.  If you add new linear algebra, keep
   it in the `scipy.linalg` family.
4. **Psi4 32-bit on Windows**: The Psi4 conda package may be 32-bit,
   limiting memory to ~2 GB.  Large basis sets (triple-zeta or larger)
   may cause out-of-memory.
5. **Avogadro error popup**: If Avogadro shows an error popup rather than
   orbitals, the JSON response likely has an error field.  The debug_log/
   pattern can capture input/output.

### Useful CLI one-liners

```powershell
# Run the plugin as Avogadro would
pixi run --as-is avogadro-ibo ibo < debug_log\input.json

# Quick energy + Molden check
pixi run python -c "from avogadro_ibo.calcs import compute_ibo; ..."

# Clear calcs/last/ and rerun
Remove-Item -Recurse -Force calcs\last\ && pixi run python ...
```

---

## 6. Critical Context (Non-Obvious Constraints)

1. **MergedChannels**: Avogadro merges stdout+stderr.  `stripLeadingNonJson`
   finds the first `{` or `[`.  **Nothing on stderr may contain these chars
   before the JSON output.**  This means:
   - Suppress Psi4 Python logging to stderr (`propagate=False`, redirect to
     file).
   - No `logger.debug(f"...{dict}...")` or `f"...{list}..."`.
   - All debug before `print(json.dumps(...))`.

2. **PM cannot separate same-atom orbitals**: The PM functional
   $L = \sum_A \sum_i n_A(i)^4$ uses only atomic populations.  Orbitals
   with DOM ≈ 1.0 on the same atom (e.g. O 2s and O lone pair) are
   degenerate — any rotation within the subspace gives identical $L$.
   **Post-PM Fock diagonalisation is mandatory** for correct s/p separation.

3. **IAO basis ≠ SCF basis**: IAOs are built from the STO-3G minimal basis
   regardless of the SCF basis.  The IAO construction (Appendix C)
   orthonormalises them via Löwdin orthogonalisation.  The q-factor
   acceptance threshold affects basis size.

4. **Valence-virtual via SVD**: The number of valence-virtual IAOs
   ($n_{\rm val\,vir}$) depends on the condition number of
   $\mathbf{C}_{\rm IAO}^T \mathbf{S} \mathbf{C}_{\rm vir}$.
   Truncation at $\sigma > 10^{-8}$ typically gives $n_{\rm min} - n_{\rm occ}$
   virtuals for simple molecules, but may differ for large systems.

5. **The random QR seed is global**: `_localize_ibos` uses
   `np.random.default_rng()` without a seed parameter.  Results are
   nondeterministic across runs *for the on-atom mixing* (fixed by
   Fock diagonalisation), but the warm-start p=2 path is convex and
   deterministic.

---

## 7. Open Items

- **Additional analysis**: bond order, charge decomposition, IRC tracking
- **Multiple method/basis**: currently only HF/cc-pVDZ is tested
- **Psi4 in-process memory limit**: consider subprocess mode for large systems
- **On-atom mixing edge cases**: verify with transition metals (d-orbitals)
- **DOM threshold tuning**: `dom_threshold=0.99` in `_resolve_on_atom_mixing`
  may need adjustment for highly polar bonds
