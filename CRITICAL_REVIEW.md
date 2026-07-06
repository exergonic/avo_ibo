Now I have a thorough understanding of all source files. Here is my critical review.

---

## Critical Review — `avo_ibo` (commit `f9a101b`)

### Findings by Priority

---

### HIGH — Correctness Issues

**1. Config key mismatch: `spin` vs `mult` breaks multiplicity persistence**

`config.py:101` stores user-selected multiplicity as `config["mult"]`, but `__init__.py:67` reads `_cfg.get("spin", spin)` and `get_config_options:74` reads `config.get("spin", 1)`. The keys don't match — multiplicity saved by the options dialog is **never applied**. (Charge uses consistent `"charge"` keys and works correctly.)

Fix: normalize all three sites to use `"spin"` or `"mult"`.

**2. Open-shell is architecturally unsupported, but no guard exists**

`compute_ibo` sets `ref = "uhf"` for spin > 1 (correct for Psi4), but the IAO pipeline is pure RHF:
- `occ_all = np.array([2.0] * nocc + ...)` — assigns occupancies 2.0 to *all* occupied, even singly-occupied alpha orbitals (`calcs.py:990`)
- `nocc = wfn.doccpi()[0] + wfn.soccpi()[0]` — counts alpha MOs only, ignores beta (`calcs.py:934`)
- `C_occ = Ca.np[:, :nocc]` — α MOs only, β block never read (`calcs.py:943`)
- `C_occ @ C_occ^T` density formula assumes closed-shell factor of 2, would be wrong for UHF (`calcs.py:635`)

A user who sets spin=3 will get a UHF SCF (correct), then the IAO pipeline will silently produce chemically wrong results. There is no `assert spin == 1` or error to prevent this.

**3. Primitive `shell.coef(p)` includes normalization — used for GTO coefficients?**

The current code uses Psi4's `molden()` output for the [GTO] section (via the temp file), so the GTO coefficients are written by Psi4 itself — this is safe. However, if any future code path writes GTO coefficients directly using `shell.coef()` (as the earlier STO-3G-based code did), the coefficients would include primitive normalisation that the Molden reader would re-apply. This happened in the reverted STO-3G approach. The current code path is fine, but the risk remains if any direct GTO-writing path is reintroduced.

---

### MEDIUM — Mathematical & Numerical

**4. IAO orthogonalization has no eigenvalue guard**

`calcs.py:134-135`:
```python
evals, evecs = np.linalg.eigh(metric)
C_IAO = C_IAO @ (evecs @ np.diag(evals ** -0.5) @ evecs.T)
```
If any eigenvalue of the IAO metric drops below ~1e-15 due to near-linear-dependence, `evals ** -0.5` produces NaN or infinity. No `clip(min=...)` guard. For well-separated atoms this won't occur, but for systems with ghost atoms, custom minimal bases, or STO-3G on heavy atoms, it could.

Fix: `np.maximum(evals, 1e-14) ** -0.5`.

**5. PM convergence: duplicate j-loop skips i < j pairs redundantly**

`calcs.py:202-203`: `for i in range(1, n_occ)` then `for j in range(i)`. This gives pairs (i=1,j=0), (i=2,j=0), (i=2,j=1), etc. This is correct and matches Appendix D's `for i=2..N, for j=1..i-1`.

However, the convergence criterion `abs(Aij) <= conv` at lines 222/241 skips the rotation but doesn't reset the pair's contribution to `grad_norm`. If the final sweep has all Aij ≈ 0 but one pair with Aij ≈ conv+ε, that single rotation contributes `(p*ϕ*Bij)²` and the algorithm exits. This is correct behavior — one final rotation then convergence. No issue.

**6. No `cho_factor` fallback for near-singular overlap matrices**

`calcs.py:107,114,121` — three Cholesky decompositions. If any fails (ill-conditioned overlap, ghost basis functions), the entire computation crashes with a LapackError. Psi4 guarantees positive-definite overlap for normal calculations, but edge cases (diffuse bases at long range, transition metal STO-3G) could trigger this.

**7. SVD truncation at 1e-8 is hardcoded, not parameterized**

`calcs.py:976`: `n_val_vir = int(np.sum(Sigma > 1e-8))`. For all tested molecules, Sigma has a clean gap: 4 values ≈ 1.0, rest ≈ 1e-15. This works perfectly, but the threshold is an empirical constant not justified in comments.

---

### MEDIUM — Code Architecture & Convention

**8. Monolithic `_analyze_ibos` is 170 lines, mixing four concerns**

`calcs.py:446-615` combines: orbital classification (π/σ/LP/core/anti*), DOM computation, type string canonicalization, hybrid-string construction, ionic character, HOMO/LUMO marking, charge decomposition, and total Wiberg formatting. Splitting into `_classify_ibos`, `_format_ibo_table`, `_format_charge_decomposition`, `_format_total_wiberg` would improve testability and readability.

**9. `__import__` antipattern used in two sites**

`calcs.py:751` and `1007`: `psi4_molden = __import__("psi4").molden`. This is a code smell — `__import__` is a CPython implementation detail. The standard form is `from psi4 import molden` (but that would break lazy loading). At minimum, use `import psi4; molden = psi4.molden` or move the import to the top of the function.

**10. `tempfile` imported but `tempfile` module unused**

`calcs.py:748`: `import tempfile` — imported but never referenced (neither `tempfile.mkstemp` nor `tempfile.NamedTemporaryFile`). The temp file is created manually at line 750 via path manipulation.

**11. Functions `_s_char`, `_p_char`, `_d_char`, `_p_frac` duplicate logic**

`calcs.py:703-727`: These four tiny functions each repeat `np.where((atom_of == atom) & (am_of == ...))`. They're called once per orbital per atom in the analysis loop. Computing all three in one pass would be faster and eliminate code duplication.

---

### LOW — Edge Cases & Polish

**12. Total electron count wrong for UHF**  
`calcs.py:607`: `int(2 * nocc) if ref == 'rhf' else nocc` — for UHF, the actual total is `2 * doccpi[0] + soccpi[0]`. Mitigated by issue #2 (open-shell shouldn't be used), but the formula is present in the code.

**13. No error handling for Psi4 calculation failure**  
`psi4.energy()` at line 930 has no try/except. If Psi4 runs out of memory, fails SCF convergence, or the basis set isn't found, the plugin dies with a raw traceback.

**14. `_ELEM_SYMBOLS` and `_ELEMENT_NUMBERS` are duplicated**  
`calcs.py:319-326` and `__main__.py:8-14` contain identical element tables through Z=53. A single-source-of-truth module constant would prevent divergence (historically, `_ELEMENT_NUMBERS` was truncated at Ar — Gotcha 23).

**15. DOM value is counterintuitive — users may misinterpret**  
Knizia's DOM = Σ n_A(i)² gives d = 0.500 for a perfect 50/50 bond and d = 0.857 for 92.4% on one atom (O LP). First-time users expect the "degree of mono/di-atomicity" to read as a percentage. A footnote in the ibos.txt header explaining the definition would help.

---

### Summary Table

| # | Priority | Category | Issue |
|---|----------|----------|-------|
| 1 | **HIGH** | Config bug | spin/mult key mismatch — multiplicity never applied |
| 2 | **HIGH** | Open-shell | UHF pipeline produces silently wrong IBOs |
| 3 | **HIGH** | Latent risk | GTO coefficient normalisation if future code writes GTO directly |
| 4 | MEDIUM | Numerical | No eigenvalue guard in IAO Löwdin orthogonalisation |
| 5 | MEDIUM | Robustness | No Cholesky fallback for ill-conditioned overlap |
| 6 | MEDIUM | Architecture | Monolithic `_analyze_ibos` (170 lines, 4 concerns) |
| 7 | MEDIUM | Code quality | `__import__` antipattern in two sites |
| 8 | MEDIUM | Code quality | `import tempfile` imported but unused |
| 9 | MEDIUM | Efficiency | `_s_char/_p_char/_d_char` triple `np.where` per orbital |
| 10 | LOW | Edge case | Wrong UHF electron count formula |
| 11 | LOW | Robustness | No try/except around `psi4.energy()` |
| 12 | LOW | Maintainability | Duplicate element tables in calcs.py and __main__.py |
| 13 | LOW | UX | DOM defintion counterintuitive for new users |

---

On the positive side: the IAO construction (Appendix C), PM Jacobi sweeps (Appendix D with corrected formulas), on-atom Fock diagonalization, and D_PERM/D_NORM Molden reordering are all **mathematically correct** and match the Knizia 2013/2014 papers and IboView reference implementation. The SVD-based valence-virtual IAO construction matches `MakeValenceVirtuals`. The spectral-theorem proof that F_IAO has no occ-vir coupling (**mathematics.md §9**) is correct and a noteworthy contribution.

The project is fundamentally sound but has one active bug (#1) and one silent correctness risk (#2) that should be addressed before wider distribution.