# Troubleshooting

## Installation Issues

**Problem:** Creating the Avogadro plugin symlink fails with "Access Denied."

**Solution:**
1. Run PowerShell as Administrator.
2. Then run:
   ```powershell
   New-Item -ItemType SymbolicLink -Path "$env:LOCALAPPDATA\OpenChemistry\Avogadro\plugins\avo_ibo" `
     -Target "C:\path\to\avo_ibo"
   ```

### `pixi` not found

**Problem:** Running `pixi install` fails with "command not found."

**Solution:**
- Install pixi from [pixi.sh](https://pixi.sh)
- If using the Avogadro plugin, the pixi bundled with Avogadro should already be on your PATH. Verify by running `pixi --version`.
- If not, add Avogadro's pixi to your PATH: typically `C:\Program Files\Avogadro2\bin\pixi.exe`.

### Psi4 import fails: `ModuleNotFoundError: No module named 'psi4'`

**Problem:** The plugin runs but crashes with a Psi4 import error.

**Solution:**
1. Ensure Psi4 is installed in the pixi environment:
   ```powershell
   pixi run python -c "import psi4; print(psi4.__version__)"
   ```
   If this fails, install Psi4:
   ```powershell
   pixi add psi4
   ```
2. If using `pip install` without pixi, install Psi4 via conda:
   ```powershell
   conda install psi4 -c psi4
   ```
3. Restart your Python environment or Avogadro after installing Psi4.

### Lock file version warning

**Problem:** Pixi warns: "WARN the lock file is up-to-date but uses an older format (v6)."

**Solution:** This is safe to ignore if you are using Avogadro's bundled pixi (v0.66.0), which only reads v6 lock files.

If you have a newer standalone pixi and want to upgrade the lock file:
```powershell
pixi lock
```

To keep v6 format (for compatibility with Avogadro's bundled pixi):
```powershell
pixi lock --no-update
```


## Runtime Issues

### Plugin does not appear in Avogadro Extensions menu

**Problem:** After installing the plugin, the "Intrinsic Bond Orbitals" menu item is missing.

**Solution:**
1. Verify the symlink exists:
   ```powershell
   ls "$env:LOCALAPPDATA\OpenChemistry\Avogadro\plugins\"
   ```
   It should list `avo_ibo`.
2. Restart Avogadro completely (close and reopen).
3. Check that the symlink target path is correct:
   ```powershell
   (Get-Item "$env:LOCALAPPDATA\OpenChemistry\Avogadro\plugins\avo_ibo").Target
   ```

### Computation fails: "molecule contains no atoms"

**Problem:** Running a computation in Avogadro fails immediately with this error.

**Solution:**
1. Load a molecule into Avogadro before triggering the plugin.
2. Verify the molecule has atoms: use Avogadro's View menu to inspect the structure.
3. Save the structure as XYZ or similar format and use the standalone CLI to debug:
   ```powershell
   pixi run python -m avogadro_ibo molecule.xyz
   ```

### Computation times out or crashes after long runtime

**Problem:** For large molecules, the computation hangs or crashes.

**Causes and solutions:**
1. **Molecule too large:** The Pipek-Mezey localization (max 2048 Jacobi sweeps) scales as O(n³) for n occupied orbitals. Systems larger than ~100 atoms may be slow.
   - Use a smaller basis set (e.g., 6-31G\* instead of cc-pVTZ).
   - Use a faster method (e.g., HF instead of DFT).
   - Run on a faster CPU.

2. **Memory exhaustion:** Psi4 may run out of memory for very large systems.
   - Increase system RAM or reduce system size.
   - Check `calcs/last/psi4.log` for Psi4 memory errors.

3. **Numerical issues in Pipek-Mezey convergence:** If the gradient norm does not converge below 1e-12, the algorithm may fail.
   - Try a smaller basis set or simpler method.
   - Check `calcs/last/psi4.log` for warnings about eigenvalues or convergence.

### "Neutral singlet only" limitation

**Problem:** Your molecule is charged or open-shell, and the plugin rejects it.

**Solution:** avo_ibo currently supports neutral singlets only. For other electronic states:
- Use IboView standalone (supports open-shell systems).
- Use NBO 7.0 for donor-acceptor analysis (does not require exactness).
- File a feature request on [GitHub](https://github.com/exergonic/avo_ibo/issues).

## Avogadro Integration Issues

### Orbitals do not appear in Molecular Orbitals panel

**Problem:** Plugin runs without error, but no orbitals appear in Avogadro.

**Solution:**
1. Check that the Molden file was created:
   ```powershell
   ls calcs\last\
   ```
   You should see `ibo.molden`.
2. Verify the Molden file is not empty:
   ```powershell
   wc -l calcs\last\ibo.molden
   ```
3. Try loading the file manually in Avogadro:
   - File → Open → select `ibo.molden`
   - Go to View → Molecular Orbitals to inspect.
4. Check the debug log (Avogadro console or stderr) for parsing errors.

### Orbital geometry is incorrect or distorted

**Problem:** Orbitals render in Avogadro but show unexpected shapes or nodal structure.

**Solution:**
1. Verify the Molden file uses the correct basis set:
   ```powershell
   head -20 calcs\last\ibo.molden
   ```
   Should show `[GTO]` section with proper basis function ordering.
2. Check that `puream` (spherical harmonics) and `Cartesian` basis settings match:
   - Psi4 default: `puream=False` (Cartesian).
   - Avogadro expects Cartesian orbitals. If Psi4 used `puream=True`, recompute with `puream=False`.
3. Verify occupancy and energy values in the `ibos.txt` analysis table make chemical sense.

## CLI Issues

### Standalone CLI: "XYZ file not found"

**Problem:** Running `pixi run python -m avogadro_ibo molecule.xyz` fails.

**Solution:**
1. Use the absolute or correct relative path:
   ```powershell
   pixi run python -m avogadro_ibo "C:\path\to\molecule.xyz"
   ```
2. Verify the file exists:
   ```powershell
   Test-Path "C:\path\to\molecule.xyz"
   ```
3. Ensure the file is valid XYZ format (first line = atom count, second line = comment).

### CLI runs but produces no output

**Problem:** Command completes silently with no results in `calcs/last/`.

**Solution:**
1. Check for error output (redirected to stderr):
   ```powershell
   pixi run python -m avogadro_ibo molecule.xyz 2>&1 | Tee-Object output.log
   ```
2. Check `calcs/last/psi4.log` for Psi4 convergence failures.
3. Ensure molecule is neutral and singlet (no charges, unpaired electrons).

## Analysis and Output Issues

### Analysis table (`ibos.txt`) is empty or incomplete

**Problem:** The `ibos.txt` file is missing or does not contain expected data.

**Solution:**
1. Verify the plugin computation completed without errors (check Avogadro message panel).
2. Confirm `calcs/last/` was created:
   ```powershell
   ls calcs\last\
   ```
3. Check `calcs/last/psi4.log` for Psi4 errors.
4. Ensure the occupation and energy parsing logic matched your Psi4 version.

### Bond type assignment looks wrong

**Problem:** σ bonds are labeled as lone pairs, or vice versa.

**Solution:**
1. Check the degree of monoatomicity (DOM) in `ibos.txt`: orbitals with DOM ≈ 1.0 on a single atom are likely lone pairs.
2. Inspect the atomic composition column: orbitals strongly localized on two atoms are σ/π bonds.
3. Examine the isosurface shape in Avogadro to confirm chemical intuition.
4. If assignment is consistently wrong, file an issue with the molecule structure and output files.

## Reporting Issues

If you encounter a bug that is not covered here:

1. **Collect diagnostics:**
   - The input XYZ file.
   - The full contents of `calcs/last/psi4.log`.
   - The `ibos.txt` and Molden file.
   - The exact error message from Avogadro or the CLI.

2. **File an issue on GitHub:**
   - [avo_ibo Issues](https://github.com/exergonic/avo_ibo/issues)
   - Include the diagnostic files and a minimal example to reproduce the issue.

3. **Check existing documentation:**
   - [`mathematics.md`](mathematics.md) — full algorithmic details.
   - [`tutorial.md`](tutorial.md) — development and architecture.

