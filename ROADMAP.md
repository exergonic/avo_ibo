# Roadmap

## v0.2 — standalone `ibo` library
Extract the core IAO/IBO math into a pip-installable package.
`import ibo; ibo.compute_ibos(wfn)` — no CJSON required.
See [#1](https://github.com/exergonic/avo_ibo/issues/1).

## v0.3 — convenience API
`read_xyz()`, `compute_ibo_from_xyz()`, typed return values.
See [#2](https://github.com/exergonic/avo_ibo/issues/2),
[#3](https://github.com/exergonic/avo_ibo/issues/3).

## v0.4 — pixi-pack distribution
Use [`pixi-pack`](https://github.com/quantco/pixi-pack) to create a
self-extracting archive containing the full conda environment (Psi4 + numpy +
scipy + avo_ibo). Users run a single `.ps1`/`.sh` script — no pixi, no conda
needed.  See [#4](https://github.com/exergonic/avo_ibo/issues/4).

## v0.5 — enhanced analysis
Bond orders, charge decomposition.
See [#5](https://github.com/exergonic/avo_ibo/issues/5),
[#6](https://github.com/exergonic/avo_ibo/issues/6).

## Future
IRC tracking, charged/open-shell systems.
See [#7](https://github.com/exergonic/avo_ibo/issues/7),
[#8](https://github.com/exergonic/avo_ibo/issues/8).
