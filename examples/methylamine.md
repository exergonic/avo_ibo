# Methylamine — amine lone pair and C–N bonding

```
Input: methylamine.xyz
Run:   pixi run python -m avogadro_ibo examples/methylamine.xyz --method wB97X-D --basis def2-TZVP
```

The amine lone pair (orb 9, HOMO) is 98.7% on nitrogen, 24% 2s + 76%
2p — predominantly p-type, as expected for a pyramidal amine. Water
and formaldehyde taught that a lone pair with tails is already
donating. The small tails onto carbon (0.6%) and one hydrogen
(0.5%) are that donation into the methyl framework. Compare
ammonia's LP (18% s) in the first-look set: methylation pushes
slightly more s character into the lone pair. Charges: N −0.424,
C −0.046, amine H's positive. The C–N bond (Wiberg 1.026,
essentially single) is polarised toward nitrogen without drama.

The three C–H bonds split into a degenerate pair (orbs 6–7) plus a
single orb 8 distinguished by its 2py composition. The C–H bond
antiperiplanar to the lone pair is the one that can donate into
the LP's acceptor direction — the same alignment test as propene's
methyl. Small, but resolved: when three otherwise identical bonds
split, look at the dihedral to the donor.

![Methylamine HOMO](img/methylamine_HOMO.png)
![C-H hyperconjugation](img/hyperconjugation.png)
