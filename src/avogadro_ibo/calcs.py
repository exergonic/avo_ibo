"""IAO/IBO construction and Pipek-Mezey localization for Avogadro 2.

BSD 3-Clause License
Copyright (c) 2025-2026, Billy Wayne McCann
SPDX-License-Identifier: BSD-3-Clause

Backward-compatible re-export shim (implementation lives in api.py and
sibling modules).  Kept so Avogadro entry points and existing imports
(`from .calcs import compute_ibo`, `from .calcs import _ELEM_SYMBOLS`)
keep working unchanged.

References:
  G. Knizia, JCTC 2013, 9, 4834-4843.  DOI: 10.1021/ct400687b
  ("Intrinsic Atomic Orbitals: An Unbiased Bridge between Quantum
   Theory and Chemical Concepts.")
  W. D. Derricotte and F. A. Evangelista, JCTC 2017, 13, 5984-5999.
  DOI: 10.1021/acs.jctc.7b00493 ("Localized Intrinsic Valence Virtual
  Orbitals ...").  Valence-virtual construction (their eqs 3-5).

Paper equation numbers and appendix references refer to Knizia 2013
unless marked D&E2017.
"""
from .api import (
    IBOResult,
    compute_ibo,
    compute_ibo_data,
    _mol_formula,
    _option,
    _sanitize_name,
    _write_input_xyz,
)
from .constants import _ELEM_SYMBOLS

__all__ = [
    "IBOResult",
    "compute_ibo",
    "compute_ibo_data",
    "_ELEM_SYMBOLS",
    "_mol_formula",
    "_option",
    "_sanitize_name",
    "_write_input_xyz",
]
