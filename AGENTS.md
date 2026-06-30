# AGENTS.md — Session Context

## Goal
Build a Pixi-based Avogadro 2 Python script plugin that computes Intrinsic
Atomic Orbitals (IAOs / Knizia IBOs) using Psi4 and returns them as Molden
data for interactive isosurface viewing.

## Status

### What Works
- **test-molden**: Menu command displays hardcoded STO-3G methane orbitals
  (Molecular Orbitals menu appears, isosurfaces renderable)
- **test-energy**: Menu command returns input CJSON with `totalEnergy: -39.5`,
  `readProperties: true` enables energy display in Avogadro
- **Compute IBOs (real Psi4)**: Menu command runs Psi4 HF → IAO projection →
  Molden output end-to-end; orbitals appear in Avogadro's MO menu
- All three menu commands show message popups on completion
- Manual CLI pipeline (`pixi run --as-is avogadro-ibo ibo`) works with exit 0

### What's Broken / None
- No known bugs. The plugin is functional for HF/DFT IAO computation.

### What's Next (User's Intent)
- Expand the plugin with additional features (user mentioned plans to expand)

## Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Build backend | `uv_build>=0.10.2,<0.11.0` | Avogadro-bundled pixi v0.66.0 doesn't support hatchling |
| Lock file version | v6 | Bundled pixi v0.66.0 can't read v7 |
| Package install | `pip install -e .` (manual, not `[tool.pixi.pypi-dependencies]`) | pypi-dependencies requires lock v7 |
| Plugin discovery | Symlink in `%LOCALAPPDATA%\OpenChemistry\Avogadro\plugins\` | Avogadro scans this directory |
| Psi4 integration | In-process (`import psi4` in compute function) | Simpler, faster, better error handling than subprocess |
| Basis for IAO projection | STO-3G | MINAO unavailable in Psi4 |
| Output format | `{"readProperties": true, "moleculeFormat": "molden", "molden": ..., "cjson": ..., "message": ...}` | Matches avo_xtb reference plugin |
| Stderr discipline | Absolutely no `{` or `[` on stderr before JSON output | `stripLeadingNonJson` finds first `{`/`[` in merged stdout+stderr channel |

## Gotchas Hit

1. **Build backend**: hatchling → Avogadro's pixi can't install. Fix: use `uv_build`.
2. **Lock file v7**: User pixi creates v7 locks, bundled pixi v0.66.0 reads only v6. Fix: use bundled pixi to install or manually maintain v6 lock.
3. **Missing entry point**: `pip install -e .` needed to create `.exe` shim under `.pixi/envs/default/Scripts/`. Fix: always run after `pixi install`.
4. **`logger.debug()` after `print(json.dumps(...))`**: Trailing stderr data after JSON causes `QJsonDocument::fromJson` to fail. Fix: all debug output before final `print()`.
5. **`logger.debug(f"Input keys: {list(data.keys())}")`**: Produces `['charge', 'cjson', ...]` on stderr — `[` found by `stripLeadingNonJson` before JSON `{`. Fix: use `', '.join(data.keys())`.
6. **Psi4 INFO logs on stderr**: `keywords={'D_CONVERGENCE': 1e-08, ...}` contains `{` found first by `stripLeadingNonJson`. Fix: redirect psi4 Python logger to file, `propagate=False`.
7. **CJSON coords as dict**: Sometimes `coords` is `{"3d": [x, y, z, ...]}` instead of flat array. Fix: check `isinstance(coords_raw, dict)` and extract `["3d"]`.
8. **UTF-8 on Windows**: Avogadro expects UTF-8 stdout. Fix: `sys.stdout.reconfigure(encoding="utf-8")`.
9. **Wrong plugin directory**: Avogadro scans `%LOCALAPPDATA%\OpenChemistry\Avogadro\plugins\`, not `%APPDATA%`.

## Relevant Files

- `C:\Users\mccan\Documents\Code\avo_ibo\pyproject.toml` — Project config, build system, avogadro metadata, 3 menu-commands
- `C:\Users\mccan\Documents\Code\avo_ibo\pixi.toml` — Pixi environment (`pixi.toml` not `pixi.toml`)
- `C:\Users\mccan\Documents\Code\avo_ibo\src\avogadro_ibo\__init__.py` — CLI entry point, stdin/stdout JSON I/O, command dispatch
- `C:\Users\mccan\Documents\Code\avo_ibo\src\avogadro_ibo\calcs.py` — In-process Psi4 IAO computation + Molden generation
- `C:\Users\mccan\Documents\Code\avo_ibo\src\avogadro_ibo\test_plugins.py` — test-molden and test-energy (no Psi4 dependency)
- `C:\Users\mccan\Documents\Code\avo_ibo\AVOGADRO2_PLUGIN_TUTORIAL.md` — Full tutorial with architecture, gotchas, and bootstrap steps
- `C:\Users\mccan\Documents\Code\avo_ibo\debug_log\` — Runtime artifacts (input.json, output.json, psi4.log, result.molden)

### External Paths

- **Avogadro 2**: `C:\Program Files\Avogadro2\bin\avogadro2.exe`
- **Bundled pixi**: `C:\Program Files\Avogadro2\bin\pixi.exe` (v0.66.0)
- **User pixi**: `C:\Users\mccan\.pixi\bin\pixi.exe` (v0.70.2)
- **Plugin symlink**: `%LOCALAPPDATA%\OpenChemistry\Avogadro\plugins\avo_ibo` → project dir
- **Reference plugin**: `%LOCALAPPDATA%\OpenChemistry\Avogadro\plugins\avo_xtb\`
- **Avogadro source (2.0.0 tag)**: Used to read `pythonscript.cpp`, `interfacescript.cpp`, `packagemanager.cpp`, `command.cpp` for async pipeline understanding

## Avogadro's Async Pipeline (Key Insight)

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
Any stderr output containing these characters before the JSON on stdout
will break parsing. This is the single most critical constraint.
