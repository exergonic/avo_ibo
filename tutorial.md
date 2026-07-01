# Avogadro 2 Pixi Python Plugin Tutorial

Everything we learned building a working Psi4 IAO/IBO plugin for Avogadro 2
on Windows. This guide covers the gotchas that aren't in any single doc.

See [`mathematics.md`](mathematics.md) for the full mathematical derivation
of the IAO/IBO pipeline, and [`README.md`](README.md) for usage.

---

## Quick Start

Bootstrap a new plugin from scratch:

```powershell
# 1. Create project
mkdir my_plugin && cd my_plugin
mkdir -p src\my_plugin

# 2. Write pyproject.toml, pixi.toml, src/my_plugin/__init__.py
#    (use the templates in sections 2, 3, 5 below)

# 3. Install with the bundled pixi (produces lock file v6)
& "C:\Program Files\Avogadro2\bin\pixi.exe" install

# 4. Install the editable entry point
pixi run pip install -e .

# 5. Create the symlink so Avogadro finds the plugin
#    (run PowerShell as Admin)
New-Item -ItemType SymbolicLink -Path "$env:LOCALAPPDATA\OpenChemistry\Avogadro\plugins\my_plugin" -Target "C:\path\to\my_plugin"

# 6. Test from CLI
pixi run --as-is avogadro-my-plugin test < debug_log\input.json

# 7. Launch Avogadro and click the menu command
Start-Process "C:\Program Files\Avogadro2\bin\avogadro2.exe"
```

The avogadro2.exe path is `C:\Program Files\Avogadro2\bin\avogadro2.exe`.
The bundled pixi path is `C:\Program Files\Avogadro2\bin\pixi.exe`.

---

## Algorithm Walkthrough

### How it works

1. **Avogadro** sends the current molecule geometry (CJSON) to the plugin.
2. **Psi4** runs a HF/DFT calculation (default: HF/cc-pVDZ, Cartesian
   functions, `puream=0`).  The SCF basis can be overridden via options.
3. **Intrinsic Atomic Orbitals (IAOs)** are constructed using the full
   IAO/2014 algorithm (Knizia Appendix C, matching IboView's
   `MakeIaoBasisNew`).  A minimal basis (STO-3G) defines the split between
   occupied atomic orbitals and virtual residual.
4. **Occupied IBOs** are obtained via Pipek-Mezey localization (p=4, eq 4)
   in the IAO basis (Jacobi sweeps, conv 1e-12, max 2048 iterations).
   A warm-start strategy (p=2 first, then p=4 refine) avoids local minima.
5. **On-atom degeneracy resolution**: The PM functional uses only atomic
   populations n_A(i), so orbitals on the same atom with DOM ≈ 1
   (Degree of Mono-centricity; the sum of squares of the top two
   atomic populations; 1.00 = one-atom orbital, see mathematics.md §8.1) (e.g.
   O 2s and O lone pair) are degenerate — any rotation within that
   subspace yields the same L value.  After PM, a post-processing step
   diagonalises the Fock matrix within each same-atom, high-DOM subspace,
   restoring the aufbau ordering (s-rich lowest, p-rich highest).
6. **Valence-virtual IAOs** are constructed via SVD of `C_IAO^T @ S @ C_vir`
   (IboView `MakeValenceVirtuals`), keeping singular values > 1e-8, then
   localized with the same PM p=4 procedure.  This gives clean anti-bond
   σ\* and π\* orbitals for visualisation.
7. Orbital energies are computed as ε_i = C_i^T @ F_AO @ C_i (Fock-matrix
   diagonal elements in the AO basis), matching IboView's
   `MakeOrbitalEnergies_General`.
8. All n_min IAO-basis orbitals (occupied + valence-virtual) are written
   to a **Molden file** with custom [MO] sections, and returned to Avogadro.
9. An **analysis table** (`calcs/last/ibos.txt`) lists each orbital's
   occupancy, energy, DOM, assigned type (σ/π/LP/Core/anti\*), atomic
   composition, and s/p/d hybridization.
10. Double-click any IBO in Avogadro's **Molecular Orbitals** panel to view its
    isosurface.

### Algorithm pipeline (inside `compute_ibo`)

```
CJSON geom → Psi4 HF/cc-pVDZ → C_occ (MO coefficients)
    ↓
C_IAO, C_IAO_occ = _build_iao_basis(S_full, S12, S_min, C_occ)
    ↓
_localize_ibos(C_IAO_occ)       ← PM p=2→p=4 warm-start
    ↓
F_IAO = C_IAO^T @ F_AO @ C_IAO
    ↓
_resolve_on_atom_mixing(C_IAO_occ)  ← Fock diag per same-atom group
    ↓
occ_energies[i] = C_i^T @ F_IAO @ C_i
    ↓
Valence-virtual: SVD(SIbVir) → U_val → _localize_ibos(U_val)
    ↓
Combine, sort by energy → Molden output → analysis table
```

### Options

The `options` dict passed from Avogadro (or via stdin) can contain:

| Key      | Default   | Description |
|----------|-----------|-------------|
| `basis`  | `cc-pVDZ` | SCF basis set name (Psi4 format) |
| `method` | `hf`      | Quantum chemistry method (hf, b3lyp, etc.) |

---

## 1. Project Structure

```
my_plugin/
├── pyproject.toml          # build config + avogadro metadata
├── pixi.toml               # pixi environment (optional if using bundled pixi)
├── pixi.lock               # lock file (v6!)
├── README.md               # project overview
├── tutorial.md             # this tutorial
├── mathematics.md          # mathematical derivation
├── AGENTS.md               # session context for LLMs
├── calcs/
│   └── last/               # per-run artifacts
└── src/
    └── my_plugin/
        ├── __init__.py     # main(): CLI entry point + dispatch
        ├── calcs.py        # computation functions
        └── links.py        # auxiliary actions (open directory, etc.)
```

---

## 2. `pyproject.toml`

```toml
[project]
name = "avo_ibo"
version = "0.1.0"
requires-python = ">=3.12"
dependencies = [
    "numpy",
    "psi4",
]

[build-system]
requires = ["uv_build>=0.10.2,<0.11.0"]
build-backend = "uv_build"

[project.scripts]
avogadro-ibo = "avogadro_ibo:main"

[tool.avogadro]
identifier-prefix = "IBO"

[[tool.avogadro.menu-commands]]
identifier = "ibo"
label = "Compute IBOs"
command = "avogadro-ibo ibo"

[[tool.avogadro.menu-commands]]
identifier = "open"
label = "Go to files …"
command = "avogadro-ibo open"
```

### Critical details

- **Build backend**: Use `uv_build` (not hatchling). The bundled pixi v0.66.0
  doesn't know hatchling.
- **Identifier-prefix**: Short prefix shown in Avogadro's plugin menus.
- **Each menu-command** gets a unique `identifier` and a `command` string.
  The first word of `command` must match an entry in `[project.scripts]`.
- **Arguments after the command** (e.g., `ibo`, `open`) are passed as
  `sys.argv[1]` — dispatch on them in `main()`.
- Use only the two menu commands you actually need.  Remove boilerplate
  test commands from the template.

---

## 3. `pixi.toml` — Environment Setup

```toml
[project]
name = "avo_ibo"
version = "0.1.0"
channels = ["conda-forge"]
platforms = ["win-64"]

[tasks]

[dependencies]
python = ">=3.12"
pip = "*"
psi4 = "*"
numpy = "*"
```

### Lock file v6 — Critical Gotcha

The **Avogadro-bundled pixi (v0.66.0)** only reads **lock file v6**.
Newer pixi versions create v7 lock files. If someone runs your plugin
through Avogadro, pixi v0.66.0 will fail with an unrecognised version error.

**Solutions** (pick one):

1. **Use bundled pixi exclusively**: run `& "C:\Program Files\Avogadro2\bin\pixi.exe" install`
   instead of your user pixi. The bundled pixi writes v6 locks.
2. **Or** generate lock with bundled pixi once, then use user pixi for development.
3. **Or** manually edit the lock file header if needed.

### Do NOT use `[tool.pixi.pypi-dependencies]`

The `pypi-dependencies` feature requires lock file v7+. Use conda-forge
dependencies in `[dependencies]` and manually `pip install -e .` for the
entry point.

### pip install -e . for Entry Point

After `pixi install`, run:

```powershell
pixi run python -m pip install -e .
```

This creates the `.exe` shim under `.pixi\envs\default\Scripts\avogadro-ibo.exe`.
Without it, `[project.scripts]` won't produce a runnable command.

---

## 4. Plugin Registration (Symlink)

Avogadro scans `%LOCALAPPDATA%\OpenChemistry\Avogadro\plugins\` for
packages. Create a symlink to your project:

```powershell
# From an admin PowerShell
New-Item -ItemType SymbolicLink -Path "$env:LOCALAPPDATA\OpenChemistry\Avogadro\plugins\avo_ibo" -Target "C:\path\to\avo_ibo"
```

### Symlink vs Junction

- **Symlink**: Preferred. Avogadro follows symlinks. Requires admin or
  Developer Mode.
- **Junction**: Works but some tools handle them differently.
- **Copy**: Also works but you must re-copy after every change.

---

## 5. Entry Point (`__init__.py`)

```python
import argparse
import json
import logging
import sys
import traceback
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

logger = logging.getLogger(__name__)

def main():
    logging.basicConfig(
        level=logging.DEBUG,
        stream=sys.stderr,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    parser = argparse.ArgumentParser()
    parser.add_argument("feature", nargs="?", default="ibo")
    parser.add_argument("--lang", default="en_US")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    raw = sys.stdin.read()
    data = json.loads(raw)

    cjson = data.get("cjson", {})
    options = data.get("options", {})
    charge = data.get("charge", 0)
    spin = data.get("spin", 1)

    try:
        if args.feature == "ibo":
            from .calcs import compute_ibo
            result = compute_ibo(cjson, options, charge, spin)
        elif args.feature == "open":
            from .links import open_calcs_dir
            result = open_calcs_dir(cjson)
        else:
            result = {"error": f"Unknown feature: {args.feature}"}
    except Exception as e:
        result = {"error": "".join(traceback.format_exception(e))}
        logger.exception("Unhandled exception")

    print(json.dumps(result))

if __name__ == "__main__":
    main()
```

### Key details

- **stdin/stdout**: Avogadro writes JSON to stdin, reads JSON from stdout.
- **`sys.stdout.reconfigure(encoding="utf-8")`**: Required on Windows;
  Avogadro expects UTF-8.
- **`print(json.dumps(result))`** is the *last* thing that touches stdout
  — no debug prints after it.
- **Error handling**: Catch all exceptions, return as `{"error": ...}` JSON.
- **Dispatch by `args.feature`**: Each menu-command's identifier maps to
  a feature string.  Keep the dispatch table simple — one `elif` per
  feature.

---

## 6. Per-Calculation Directory Pattern

For a clean file system, keep per-run artifacts in a dedicated directory
(e.g. `calcs/last/`).  Clear it before each `"ibo"` run so the user always
finds the latest result:

```python
from pathlib import Path
import shutil

CALCS_DIR = Path(__file__).parent.parent.parent / "calcs"
TEMP_DIR = CALCS_DIR / "last"

def _prepare_calc_dir():
    if TEMP_DIR.exists():
        shutil.rmtree(TEMP_DIR)
    TEMP_DIR.mkdir(parents=True, exist_ok=True)
```

Call `_prepare_calc_dir()` **only** inside the `"ibo"` branch — the
"open" feature should never wipe previous results.

---

## 7. The `{` / `[` Stderr Problem (Three Gotchas)

**This is the single most important thing to understand.**

Avogadro uses `MergedChannels` — stdout and stderr are merged into one
buffer. After the process finishes, `stripLeadingNonJson()` is called:

```cpp
QByteArray stripLeadingNonJson(QByteArray data) {
    int start = data.indexOf('{');
    if (start == -1)
        start = data.indexOf('[');
    if (start != -1)
        data = data.mid(start);
    return data;
}
```

It finds the **first** `[` or `{` in the merged buffer and tries to parse
JSON starting there. If **anything** on stderr contains `[` or `{`
*before* the JSON `print()` to stdout, parsing fails.

### Gotcha 1: Trailing stderr after JSON print

```python
# BAD: stderr has content AFTER the JSON is printed
print(json.dumps(result))      # stdout: JSON
logger.debug("Output keys ...")  # stderr: appended after JSON
```

Even though `stripLeadingNonJson` finds the correct `{`, the trailing
stderr data causes `QJsonDocument::fromJson` to fail.

**Fix:** Move all debug output *before* the `print()`.

### Gotcha 2: `[` in debug strings

```python
# BAD: produces ['charge', 'cjson', ...] on stderr
logger.debug(f"Input keys: {list(data.keys())}")
# stderr: "Input keys: ['charge', 'cjson']"
```

The `[` at position 13 is found by `stripLeadingNonJson` before the
JSON `{`. The function tries to parse `['charge', 'cjson']` as JSON.

**Fix:**
```python
# GOOD: no square brackets
logger.debug(f"Input keys: {', '.join(data.keys())}")
```

### Gotcha 3: Psi4 logging with `{` on stderr

Psi4 writes INFO-level Python logging messages to stderr like:

```
PLANNING Atomic: keywords={'D_CONVERGENCE': 1e-08, ...}
Compute energy(): method=hf, basis=cc-pVDZ, ...
```

The `{` in `keywords={'D_CONVERGENCE': ...}` is found by
`stripLeadingNonJson` first.

**Fix:** Redirect Psi4 logging away from stderr into a file:

```python
import logging

_psi_logger = logging.getLogger("psi4")
_psi_logger.propagate = False
_psi_logger.setLevel(logging.WARNING)
_psi_handler = logging.FileHandler(str(log_path), mode="w")
_psi_logger.addHandler(_psi_handler)

# Also set sub-loggers
for name in ["psi4.core", "psi4.driver"]:
    logging.getLogger(name).setLevel(logging.WARNING)
```

---

## 8. The General Rule

> **Nothing on stderr may contain `{` or `[` before the JSON output
> on stdout.**

This means:

- No `logger.debug(f"...{dict}...")` calls (dict repr has `{}`)
- No `logger.debug(f"...{list}...")` calls (list repr has `[]`)
- No `traceback.format_exc()` on stderr before `print()` (tracebacks
  don't generally contain `{`/`[`, but be careful)
- Suppress or redirect any library that logs to stderr with `{`/`[`
  (numpy, scipy, psi4, etc.)

---

## 9. Debug Log Pattern

Save `input.json`, `output.json`, and any computation artifacts to a
`debug_log/` directory. This is invaluable for troubleshooting because
Avogadro swallows stderr and you can't see what happened:

```python
_DEBUG_DIR = Path(__file__).parent.parent.parent / "debug_log"

_DEBUG_DIR.mkdir(parents=True, exist_ok=True)
(_DEBUG_DIR / "input.json").write_text(raw, encoding="utf-8")

# ... compute ...

(_DEBUG_DIR / "output.json").write_text(
    json.dumps(result, indent=2), encoding="utf-8"
)
print(json.dumps(result))
```

---

## 10. CJSON Coordinate Format

The CJSON coordinates can come in two forms:

```python
# Form 1: flat array (most common)
coords = cjson["atoms"]["coords"]  # [x1, y1, z1, x2, y2, z2, ...]

# Form 2: dict with "3d" key
coords_raw = cjson["atoms"]["coords"]
if isinstance(coords_raw, dict):
    coords = coords_raw["3d"]
else:
    coords = coords_raw
```

Always handle both.

---

## 11. "Go to files …" Feature

Auxiliary menu commands (like opening the output directory) are simpler:
they don't need Psi4 or heavy computation.

```python
# links.py
from pathlib import Path
from . import CALCS_DIR

def open_calcs_dir(cjson):
    import subprocess
    subprocess.Popen(["explorer", str(CALCS_DIR)])
    return {
        "readProperties": True,
        "cjson": cjson,
        "message": f"Opened {CALCS_DIR} in Explorer",
    }
```

This is registered as a second menu-command in `pyproject.toml` with
`identifier = "open"`.

---

## 12. Manual Pipeline Testing

Test the full pipeline from command line to verify:

```powershell
# Run with test input
pixi run --as-is avogadro-ibo ibo < debug_log\input.json

# Check exit code (0 = success)
$LASTEXITCODE

# Check output.json for valid JSON
Get-Content debug_log\output.json | python -m json.tool
```

This simulates what Avogadro does (minus the MergedChannels), and
saves `input.json`/`output.json` for inspection.

---

## 13. Psi4-Specific Notes

### In-process Psi4 (recommended)

Import Psi4 in-process inside the compute function (not at module level)
to avoid import-time side effects:

```python
def compute_ibo(...):
    import psi4
    import numpy as np
```

### Redirect Psi4 output

```python
psi4.set_output_file(str(TEMP_DIR / "psi4.log"), append=True)
```

This captures Psi4's C++ output (scf iterations, energy, etc.) to a file.
Combine with the Python logging redirect above for full capture.

### `puream=0` (Cartesian) is required for Avogadro

Psi4's Molden writer (`psi4.molden()`) with the default spherical harmonic
basis (`puream=1`, or unspecified) produces files whose [GTO] shell layout
and [MO] coefficient ordering are silently misaligned under Avogadro 2's
Molden parser.  Orbitals appear to load — they have isosurfaces and energy
values — but are chemically wrong.  The mismatch is invisible unless you
already know what shape to expect (e.g. checking that ethene's π orbital
actually looks like a π orbital).

**Always use Cartesian functions (`puream=0`)** when generating Molden
files for Avogadro.  Our plugin does this by default for both the SCF and
minimal IAO basis sets.  Ethene B3LYP/cc-pVDZ is the canonical
reproduction: run `psi4.molden()` with `puream=1` and load the result
into Avogadro — the occupied π orbital appears distorted or missing its
characteristic nodal plane.

### Psi4 Molden output (IAO basis)

Use a custom Molden writer that replaces the [MO] section with
IAO-basis orbitals:

```python
def _write_iao_molden(path, wfn, C_AO, occ, energies, n_orb):
    import tempfile
    tmp = path.with_suffix(".molden.tmp")
    psi4_molden = __import__("psi4").molden
    psi4_molden(wfn, str(tmp))
    text = tmp.read_text(encoding="utf-8")
    tmp.unlink()

    mo_tag = "[MO]"
    idx = text.find(mo_tag)
    header = text[:idx]

    n_AO = C_AO.shape[0]
    lines = [header + "\n[MO]\n"]
    for i in range(n_orb):
        ei = energies[i]
        oi = occ[i]
        lines.append(f" Sym= A\n Ene= {ei:15.10f}\n Spin= Alpha\n"
                     f" Occup= {oi:14.10f}\n")
        for j in range(n_AO):
            lines.append(f"  {j + 1:>4d}  {C_AO[j, i]:16.10f}\n")

    # Pad to n_AO orbitals so Avogadro's slot count matches [GTO]
    for i in range(n_orb, n_AO):
        lines.append(f" Sym= A\n Ene= {0.0:15.10f}\n Spin= Alpha\n"
                     f" Occup= {0.0:14.10f}\n")
        for j in range(n_AO):
            lines.append(f"  {j + 1:>4d}  {0.0:16.10f}\n")

    path.write_text("".join(lines), encoding="utf-8")
```

This padding works around an upstream bug: Avogadro allocates one MO
display slot per [GTO] basis function rather than per [MO] entry,
so files with fewer orbitals than basis functions show uninitialised
noise in the extra slots.  Filed as
[OpenChemistry/avogadrolibs #2890](https://github.com/OpenChemistry/avogadrolibs/issues/2890).

Return in the JSON response:

```python
{
    "readProperties": True,
    "moleculeFormat": "molden",
    "molden": molden_text,
    "cjson": cjson,
    "message": "IBO analysis saved to calcs/last/ibos.txt",
}
```

The `molden` field provides the orbitals. `moleculeFormat: "molden"` tells
Avogadro to read orbital data from the molden string. `readProperties: True`
enables property display (energy). `cjson` provides the 3D structure.

The analysis table is written separately to `calcs/last/ibos.txt` (not
embedded in the popup message) to avoid blocking the user.

---

## 14. Test Plugins (No Heavy Dependencies)

Create separate test functions that don't require Psi4. This lets you
verify the Avogadro pipeline works before debugging computation logic:

```python
def test_molden(cjson, options, charge, spin):
    """Hardcoded STO-3G methane molden data — no psi4 needed."""
    molden_text = [...]  # embed or generate molden data
    return {
        "readProperties": True,
        "moleculeFormat": "molden",
        "molden": molden_text,
        "cjson": cjson,
        "message": "Test Molden (STO-3G methane)",
    }


def test_energy(cjson, options, charge, spin):
    """Return input CJSON with a fake energy — tests readProperties."""
    cjson["properties"]["totalEnergy"] = -39.5
    return {
        "readProperties": True,
        "cjson": cjson,
        "message": "Test Energy (-39.5 hartree)",
    }
```

---

## 15. Giving Context to Another LLM (or Future You)

To give an LLM full context for continuing work on this plugin, share
**these files**:

1. **`tutorial.md`** — explains the architecture and gotchas.
2. **`mathematics.md`** — full mathematical derivation with paper references.
3. **`AGENTS.md`** — session context: current status, pending items,
   key decisions, gotchas hit.

Keep `AGENTS.md` updated as you go. A good format:

```markdown
# AGENTS.md — Session Context

## Goal
Briefly: what are we building?

## Status
- What works
- What's broken
- What's next

## Key Decisions
- Build backend: uv_build
- Lock file: v6 (for Avogadro's bundled pixi v0.66.0)
- ... etc

## Gotchas Hit
- #1: ...
- #2: ...

## Relevant Files
- mathematics.md — mathematical derivation
- README.md — usage
- tutorial.md — architecture and gotchas
```

When starting a new session, present `AGENTS.md` + relevant snippets
from `tutorial.md` and `mathematics.md` as context. The LLM will be up
to speed in one read.

---

## 16. GitHub MathJax Rendering

`mathematics.md` uses MathJax v3 to render LaTeX display math on GitHub.
GitHub's pipeline has two independent failure modes for display math:

### 16.1 Markdown preprocessor brace mangling

GitHub runs a Markdown preprocessor *before* MathJax that [mangles LaTeX
inside `$$...$$` blocks](https://stackoverflow.com/a/77726873) (confirmed
by MathJax core dev Davide Cervone).  The preprocessor tries to protect
special characters (`_`, `{`, `}`, `\`) but its algorithm is buggy and
can corrupt braces, causing MathJax's "Missing close brace" error.

**Fix**: use ` ```math ```` code block syntax instead of `$$...$$`.
The code block bypasses the Markdown preprocessor entirely.

### 16.2 MathJax `\Vert` + `<` parsing error

Even inside a ` ```math ```` code block, the combination of `\Vert` and bare
`<` (in `\sum_{i<j}`) triggers a "Missing close brace" error in MathJax v3
itself (or GitHub's specific MathJax configuration).  The equation

```math
\Vert\nabla L\Vert = \frac{1}{n_{\mathrm{occ}}} \sqrt{\sum_{i<j} \bigl(p \,\phi_{ij}\, B_{ij}^{(p)}\bigr)^2}
```

produces the error even with ` ```math ```` isolation.

**Fix**: use `\lVert`/`\rVert` instead of `\Vert`, and `\lt` instead of `<`:

```math
\lVert\nabla L\rVert = \frac{1}{n_{\mathrm{occ}}} \sqrt{\sum_{i\lt j} \bigl(p \,\phi_{ij}\, B_{ij}^{(p)}\bigr)^2}
```

Both changes together are required; either alone may suffice but both are
applied for robustness.

---

## 17. Distributing the Environment with pixi-pack

[`pixi-pack`](https://github.com/quantco/pixi-pack) bundles a pixi environment
into a self-extracting archive.  Users on a clean machine run a single script
and get the full computation environment — no pixi, no conda, no Python install
required.

### Building the pack

Install pixi-pack (once):

```powershell
pixi global install pixi-pack pixi-unpack
```

From the project directory, create the self-extracting executable:

```powershell
pixi-pack --create-executable
```

This reads the lockfile, downloads all conda packages, and produces
`environment.ps1` (~500 MB).  The pack contains Psi4, numpy, scipy, pytest,
ruff — everything declared in `[tool.pixi.dependencies]`.

> **Note:** The `avo_ibo` package itself is installed via `pip install -e .`
> and is **not** included in the pack (it is not a conda package).  It is
> tiny (~100 KB) and added separately.

### What the user sees

The self-extracting script produces three items:

| File / Dir | Size | Purpose |
|------------|------|---------|
| `environment.ps1` | 500 MB | Self-extracting archive (run once) |
| `env/` | 1.9 GB | The extracted conda environment |
| `activate.bat` | — | Activates the environment (sets PATH, CONDA_PREFIX) |

### End-to-end workflow on a clean machine

```powershell
# 1. Extract the Psi4 environment (~2 min)
.\environment.ps1

# 2. Activate it
.\activate.bat

# 3. Install the thin plugin layer (once avo_ibo is on PyPI)
pip install avogadro-ibo

# 4. Run a calculation
python -m avogadro_ibo molecule.xyz
```

The environment at `env/` contains a complete Python 3.14 + Psi4 installation.
After activation, all standard Python and conda commands work:

- `python -c "import psi4; print(psi4.__version__)"`
- `pip list` shows all installed packages
- `pytest tests/` runs the test suite (if included)

### Rebuilding for updates

When the lockfile changes (new dependencies, updated versions):

```powershell
pixi install --locked     # recreate env from lockfile
pixi-pack --create-executable
```

The new `environment.ps1` is a drop-in replacement — users just re-run it.

---

## Summary of Gotchas

| # | Gotcha | Symptom | Fix |
|---|--------|---------|-----|
| 1 | `[tool.pixi.pypi-dependencies]` | Lock file v7 not readable by Avogadro's pixi v0.66.0 | Remove it; use conda deps + `pip install -e .` |
| 2 | Entry point missing | "command not found" | `pixi run pip install -e .` |
| 3 | `logger.debug` after `print` | Trailing non-JSON after JSON | Move all debug output before `print()` |
| 4 | `list(data.keys())` in debug output | `[` found first by `stripLeadingNonJson` | Use `', '.join(data.keys())` |
| 5 | Psi4 logging to stderr | `{` found first by `stripLeadingNonJson` | Redirect psi4 logger to file, `propagate=False` |
| 6 | CJSON coords as dict | IndexError on `coords[3*i]` | Handle both flat array and `{"3d": [...]}` forms |
| 7 | UTF-8 on Windows | UnicodeDecodeError in Avogadro | `sys.stdout.reconfigure(encoding="utf-8")` |
| 8 | hatchling build backend | pixi can't install | Use `uv_build` |
| 9 | Plugin not showing up in Avogadro | Wrong plugin directory | Use `%LOCALAPPDATA%\OpenChemistry\Avogadro\plugins\` |
| 10 | Psi4 Molden with spherical harmonics | Orbitals load but are chemically wrong in Avogadro | Use `puream=0` (Cartesian) |
| 11 | [MO] entries < [GTO] basis count | Extra MO slots show junk noise in Avogadro | Pad [MO] with zero-energy dummies up to n_AO ([#2890](https://github.com/OpenChemistry/avogadrolibs/issues/2890)) |
| 12 | GitHub MathJax `$$` + `\Vert` + `<` | "Missing close brace" error on `mathematics.md` display math | ` ```math ```` code block + `\lVert`/`\rVert` + `\lt` (§16) |
| 13 | pip-installed packages missing from pixi-pack | `avo_ibo` not available after unpacking | `pixi-pack` reads the lockfile, not site-packages; install plugin separately with `pip install avogadro-ibo` (§17) |
