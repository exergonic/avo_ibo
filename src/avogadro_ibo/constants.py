"""Shared numeric constants for the IBO pipeline."""


# 1 Hartree in electron-volts (CODATA 2018) and in kcal/mol
# (2625.4996394799 kJ/mol ÷ 4.184).  Only used for the human-facing
# HOMO-LUMO gap line; all internal energies stay in Ha.
HA_TO_EV = 27.211386245988
HA_TO_KCAL = 627.5094740631

# Floor for the smallest kept VVO singular value.  Anything below the
# floor means the "valence" label is suspect (pathological SCF basis,
# near-linear dependence).  Set 20x below the observed suite minima
# (all 1.000 to three decimals, cc-pVDZ through aug-cc-pVDZ) and 5x
# above junk scale (~0.01): it never fires on sane input.
VVO_MIN_SIGMA = 0.05

# Periodic-table lookup for element symbols
_ELEM_SYMBOLS = [
    "X",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
]

# Orbital-pair interference terms with |term| at or above this value get
# their own lines in the "Significant orbital-pair interference" detail
# section.  At 0.01 the section fires only on genuine delocalization
# chemistry (3c-2e bridges, conjugated π) and stays silent on ordinary
# single/double bonds.
PAIR_DETAIL_THRESH = 0.01
# Near-miss band [PAIR_NEAR_THRESH, PAIR_DETAIL_THRESH): counted in a
# closing footnote (count plus largest term) whenever the section
# prints — so a parenthetical like (+0.010) always names its pair
# (cf. dimer O1-O4 at +0.0095). The footnote never creates a section
# on its own: quiet molecules stay fully silent (see ethene test).
PAIR_NEAR_THRESH = 0.005
