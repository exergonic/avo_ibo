# Ozone — charge separation and through-bond coupling

```
Input: ozone.xyz  (ORCA 6.1.1 wB97X-D3/def2-TZVP opt + freq, no imaginary modes)
Run:   pixi run python -m avogadro_ibo examples/ozone.xyz --method wB97X-D --basis def2-TZVP
```

The table tells ozone's whole reactivity story at a
glance. Charges: central O +0.376, terminals −0.188 each — the famous
charge separation behind its electrophilicity, straight out of the
resonance structures but now as numbers. The O–O bonds read 1.420 (σ
0.954 + π 0.466), and the terminal O···O contact — two atoms not bonded
in any Lewis structure — carries **0.550**, over half a bond: σ 0.002
+ π 0.548, essentially pure π, with the detail section showing
`orb9 × orb10: +0.2174`,
the 3-centre π bond leaking across the terminal pair. For scale:
benzene's meta C–C is 0.116, SO₃'s O–O 0.288 — ozone's direct bent
π pathway nearly doubles SO₃'s threefold-symmetric one.

The frontier picture completes it: LUMO is the symmetric O π* at
−0.056 Ha — bound and low-lying, hence the extraordinary
reactivity — while the terminal-oxygen lone pairs (−0.488 Ha, 97.5%
2p)
are the nucleophilic sites for 1,3-dipolar cycloaddition. Electrophile
and nucleophile, HOMO and LUMO, in one table.

![Ozone LUMO](img/ozone_lumo.png)

---
*Geometry optimized in ORCA 6.1.1 (wB97X-D3/def2-TZVP), confirmed
minimum by frequency calculation (no imaginary modes): O–O 1.2357 Å
vs 1.272 experimental, angle 118.0° vs 116.8°. IBO analysis with
Psi4 1.11 (wB97X-D/def2-TZVP, RHF), avo_ibo 0.4.0; Pipek–Mezey
localization (p = 2 → 4) per Knizia 2013.*
