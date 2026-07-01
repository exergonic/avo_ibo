"""CLI integration tests for avogadro_ibo.

Patterned after avo_xtb: uses tests/files/ for reference inputs.
Runs the standalone CLI (python -m avogadro_ibo <file.xyz>) and
validates outputs in the project's calcs/last/ directory.
"""

import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_DIR = Path(__file__).resolve().parent.parent
FILES_DIR = PROJECT_DIR / "tests" / "files"
CALCS_LAST = PROJECT_DIR / "calcs" / "last"

# (xyz, label, n_AO, n_IAO, n_occ, n_vir)
# Basis counts: Cartesian cc-pVDZ (H=5, C/N/O=15); STO-3G min (H=1, C/N/O=5)
MOLECULES = [
    ("water.xyz", "water", 25, 7, 5, 2),
    ("methane.xyz", "methane", 35, 9, 5, 4),
    ("ethene.xyz", "ethene", 50, 14, 8, 6),
    ("ammonia.xyz", "ammonia", 30, 8, 5, 3),
    ("benzene.xyz", "benzene", 120, 36, 21, 15),
]


def parse_ibos(path):
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    data = []
    for line in lines:
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("---") or s.startswith("IBO"):
            continue
        cols = line.split()
        if cols and cols[0].isdigit():
            data.append({
                "idx": int(cols[0]),
                "occ": float(cols[1]),
                "ene": float(cols[2]),
                "dom": float(cols[3]),
            })
    return data


@pytest.mark.parametrize("xyz,label,n_ao,n_iao,n_occ,n_vir", MOLECULES)
def test_cli_counts(xyz, label, n_ao, n_iao, n_occ, n_vir):
    xyz_path = FILES_DIR / xyz
    assert xyz_path.exists(), f"Missing test file: {xyz_path}"

    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0, f"CLI failed:\nstdout:{result.stdout}\nstderr:{result.stderr}"

    assert (CALCS_LAST / "ibo.molden").exists(), "ibo.molden not found"
    assert (CALCS_LAST / "canonical.molden").exists(), "canonical.molden not found"
    assert (CALCS_LAST / "ibos.txt").exists(), "ibos.txt not found"

    molden_text = (CALCS_LAST / "ibo.molden").read_text(encoding="utf-8")
    mo_count = molden_text.count("[MO]")
    assert mo_count == 1, "Expected exactly one [MO] section"

    ene_lines = [l for l in molden_text.splitlines() if l.startswith(" Ene=")]
    assert len(ene_lines) == n_ao, (
        f"Expected {n_ao} MO entries (n_AO), got {len(ene_lines)}"
    )

    ibos = parse_ibos(CALCS_LAST / "ibos.txt")
    assert len(ibos) == n_iao, (
        f"Expected {n_iao} IAO orbitals, got {len(ibos)}"
    )

    occ_ibos = [o for o in ibos if abs(o["occ"] - 2.0) < 1e-4]
    vir_ibos = [o for o in ibos if abs(o["occ"]) < 1e-4]
    assert len(occ_ibos) == n_occ, (
        f"Expected {n_occ} occupied IBOs, got {len(occ_ibos)}"
    )
    assert len(vir_ibos) == n_vir, (
        f"Expected {n_vir} virtual IBOs, got {len(vir_ibos)}"
    )


def test_water_on_atom_resolution():
    """On-atom Fock diagonalisation separates O 2s LP from O 2p LP by energy."""
    xyz_path = FILES_DIR / "water.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0

    ibos = parse_ibos(CALCS_LAST / "ibos.txt")
    lps = [o for o in ibos if abs(o["occ"] - 2.0) < 1e-4 and o["dom"] > 0.99 and o["ene"] > -10]
    assert len(lps) == 2, f"Expected 2 O LPs with DOM≈1, got {len(lps)}"

    lps_sorted = sorted(lps, key=lambda o: o["ene"])
    s_lp, p_lp = lps_sorted[0], lps_sorted[1]
    assert s_lp["ene"] < p_lp["ene"], "s-rich LP should be lower in energy than p-pure LP"

    text = (CALCS_LAST / "ibos.txt").read_text(encoding="utf-8")
    assert "55% s + 45% p" in text or "100% s" in text
    assert "100% p" in text


def test_methane_pattern():
    """Methane: 1 core + 4 identical C-H sigma (degenerate energies)."""
    xyz_path = FILES_DIR / "methane.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0

    ibos = parse_ibos(CALCS_LAST / "ibos.txt")
    occ = [o for o in ibos if abs(o["occ"] - 2.0) < 1e-4]
    assert len(occ) == 5

    core = next(o for o in occ if o["dom"] > 0.99)
    assert "C(Core)" in (CALCS_LAST / "ibos.txt").read_text(encoding="utf-8")

    ch_sigmas = [o for o in occ if o["dom"] < 0.99]
    assert len(ch_sigmas) == 4
    energies = {round(o["ene"], 3) for o in ch_sigmas}
    assert len(energies) == 1, f"C-H sigma energies should be degenerate, got {energies}"


def test_ammonia_lp():
    """Ammonia: 1 core + 3 N-H sigma + 1 N LP (DOM≈1)."""
    xyz_path = FILES_DIR / "ammonia.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0

    ibos = parse_ibos(CALCS_LAST / "ibos.txt")
    occ = [o for o in ibos if abs(o["occ"] - 2.0) < 1e-4]
    assert len(occ) == 5

    lps = [o for o in occ if o["dom"] > 0.99 and o["ene"] > -10]
    assert len(lps) == 1, f"Expected 1 N LP, got {len(lps)}"


def test_benzene_symmetry():
    """
    Benzene: all 6 C-H sigma bonds have essentially degenerate energies.

    The PM functional cannot resolve intra-atom p-orbital alignment for
    symmetry-equivalent bonds (a known limitation — see Knizia JCTC 2013,
    ibo-ref.py, and mathematics.md).  The split is < 1e-4 Ha.
    """
    xyz_path = FILES_DIR / "benzene.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=300,
    )
    assert result.returncode == 0

    ibos = parse_ibos(CALCS_LAST / "ibos.txt")
    occ = [o for o in ibos if abs(o["occ"] - 2.0) < 1e-4]
    assert len(occ) == 21, f"Expected 21 occupied IBOs, got {len(occ)}"

    text = (CALCS_LAST / "ibos.txt").read_text(encoding="utf-8")
    ch_energies = []
    for line in text.splitlines():
        s = line.strip()
        if not s or not s[0].isdigit():
            continue
        if "C-H sigma" not in line:
            continue
        cols = line.split()
        ch_energies.append(float(cols[2]))

    assert len(ch_energies) == 6, f"Expected 6 C-H sigma, got {len(ch_energies)}"

    emax = max(ch_energies)
    emin = min(ch_energies)
    assert (emax - emin) < 1e-4, (
        f"C-H sigma energies split by {emax - emin:.6f} Ha (limit < 1e-4)"
    )
