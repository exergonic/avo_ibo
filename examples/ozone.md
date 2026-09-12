# Ozone — charge separation and through-bond coupling

```
Input: ozone.xyz  (ORCA 6.1.1 wB97X-D3/def2-TZVP opt + freq, no imaginary modes)
Run:   pixi run python -m avogadro_ibo examples/ozone.xyz --method wB97X-D --basis def2-TZVP
```

The table tells ozone's reactivity story at a glance. The central
oxygen is positive (+0.376), the terminals negative (−0.188 each) —
the famous charge separation behind its electrophilicity, now as
numbers.

The O–O bonds are mixed σ and π (total 1.420). The surprise is the
terminal O···O contact: two atoms not bonded in any Lewis structure,
yet carrying over half a bond (**0.550**), essentially pure π. The
detail section names the leak — `orb9 × orb10: +0.2174` — the
3-centre π bond reaching across the terminal pair. For scale,
benzene's meta C–C is 0.116 and SO₃'s O–O is 0.288. Ozone's bent
π pathway nearly doubles SO₃'s threefold-symmetric one.

The LUMO is a bound, low-lying O π* (−0.056 Ha), hence the
extraordinary reactivity. The terminal-oxygen lone pairs are the
nucleophilic sites for 1,3-dipolar cycloaddition. Electrophile and
nucleophile in one table.

![Ozone LUMO](img/ozone_lumo.png)

---
*Geometry optimized in ORCA 6.1.1 (wB97X-D3/def2-TZVP), confirmed
minimum by frequency calculation (no imaginary modes): O–O 1.2357 Å
vs 1.272 experimental, angle 118.0° vs 116.8°. IBO analysis with
Psi4 1.11 (wB97X-D/def2-TZVP, RHF), avo_ibo 0.4.0; Pipek–Mezey
localization (p = 2 → 4) per Knizia 2013.*
