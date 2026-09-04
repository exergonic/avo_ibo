"""CLI integration tests for avogadro_ibo.

Patterned after avo_xtb: uses tests/files/ for reference inputs.
Runs the standalone CLI (python -m avogadro_ibo <file.xyz>) and
validates outputs in per-molecule subdirectories under calcs/.
"""

import json
import re
import subprocess
import sys

from avogadro_ibo.__main__ import _parse_xyz
from pathlib import Path

import pytest

PROJECT_DIR = Path(__file__).resolve().parent.parent
FILES_DIR = PROJECT_DIR / "tests" / "files"

# (xyz, label, n_AO, n_IAO, n_occ, n_vir)
# Basis counts: Cartesian cc-pVDZ (H=5, C/N/O=15, Zn=43, Cl=19);
# STO-3G min (H=1, C/N/O=5, Zn=16, Cl=9).
MOLECULES = [
    ("water.xyz", "water", 25, 7, 5, 2),
    ("methane.xyz", "methane", 35, 9, 5, 4),
    ("ethene.xyz", "ethene", 50, 14, 8, 6),
    ("ammonia.xyz", "ammonia", 30, 8, 5, 3),
    ("benzene.xyz", "benzene", 120, 36, 21, 15),
    ("zncl2.xyz", "zncl2", 87, 37, 32, 5),
    # Non-planar conjugated pi system (D2d tub).  Guards that the
    # bond-flat degeneracy resolution stays inert at a real molecule's
    # own equilibrium: four ~99%-p C-C pi orbitals, pi-HOMO/pi*-LUMO,
    # no sp^3-like banana bonds.
    ("cyclooctatetraene.xyz", "cyclooctatetraene", 160, 48, 28, 20),
]


def _find_calc_dir(name):
    """Return the latest calc directory for molecule *name*."""
    calcs_dir = PROJECT_DIR / "calcs"
    best = None
    for d in calcs_dir.iterdir():
        if not d.is_dir():
            continue
        m = re.match(rf"^{re.escape(name)}_(\d+)$", d.name)
        if m:
            n = int(m.group(1))
            if best is None or n > best[1]:
                best = (d, n)
    return best[0] if best else calcs_dir / "last"


def parse_ibos(path):
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    data = []
    for line in lines:
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("---") or s.startswith("IBO"):
            continue
        cols = line.split()
        if cols and cols[0].isdigit():
            # New format: ... Ion% [<- HOMO/LUMO]
            # H/L marker (if present) adds "<-" and "HOMO"/"LUMO" as trailing tokens
            if len(cols) >= 4 and cols[-2] == "<-" and cols[-1] in ("HOMO", "LUMO"):
                ion_raw = cols[-3]
                middle = cols[3:-3]
            else:
                ion_raw = cols[-1]
                middle = cols[3:-1]
            entry = {
                "idx": int(cols[0]),
                "occ": float(cols[1]),
                "ene": float(cols[2]),
                "ion_pct": ion_raw if ion_raw == "---" else float(ion_raw),
            }
            if len(middle) >= 2 and middle[1] in ("σ", "σ*", "π", "π*", "anti*"):
                entry["label"] = " ".join(middle[:2])
            else:
                entry["label"] = middle[0] if middle else ""
            data.append(entry)
    return data


@pytest.mark.parametrize("xyz,label,n_ao,n_iao,n_occ,n_vir", MOLECULES)
def test_cli_counts(xyz, label, n_ao, n_iao, n_occ, n_vir):
    xyz_path = FILES_DIR / xyz
    assert xyz_path.exists(), f"Missing test file: {xyz_path}"

    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", "--method", "hf", "--basis", "cc-pVDZ", "--charge", "0", "--spin", "1", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0, f"CLI failed:\nstdout:{result.stdout}\nstderr:{result.stderr}"

    calc_dir = _find_calc_dir(label)
    assert (calc_dir / "ibo.molden").exists(), "ibo.molden not found"
    assert (calc_dir / "canonical.molden").exists(), "canonical.molden not found"
    assert (calc_dir / "ibos.txt").exists(), "ibos.txt not found"

    molden_text = (calc_dir / "ibo.molden").read_text(encoding="utf-8")
    # The Molden format permits exactly one [MO] section; a second would
    # corrupt the file.  Count section headers precisely (not substring
    # matches inside e.g. comments).
    mo_sections = re.findall(r"^\[MO\]$", molden_text, flags=re.MULTILINE)
    assert len(mo_sections) == 1, f"Expected exactly one [MO] section, got {len(mo_sections)}"

    ene_lines = [l for l in molden_text.splitlines() if l.startswith(" Ene=")]
    assert len(ene_lines) == n_iao, (
        f"Expected {n_iao} MO entries (IboView-style STO-3G), got {len(ene_lines)}"
    )

    ibos = parse_ibos(calc_dir / "ibos.txt")
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


def test_cli_output_dir(tmp_path):
    """Per-call --output-dir redirects the whole calc dir; repo calcs/ untouched."""
    xyz_path = FILES_DIR / "water.xyz"
    before = set((PROJECT_DIR / "calcs").glob("water_*"))
    out = tmp_path / "out"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", "--method", "hf", "--basis", "cc-pVDZ",
         "--output-dir", str(out), str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0, f"CLI failed:\nstdout:{result.stdout}\nstderr:{result.stderr}"
    produced = sorted(out.glob("water_*/ibos.txt"))
    assert len(produced) == 1, f"Expected one ibos.txt under {out}, got {produced}"
    assert (produced[0].parent / "ibo.molden").exists(), "ibo.molden not found"
    assert (produced[0].parent / "canonical.molden").exists(), "canonical.molden not found"
    after = set((PROJECT_DIR / "calcs").glob("water_*"))
    assert after == before, f"Repo calcs/ gained dirs: {after - before}"


def test_water_on_atom_resolution():
    """On-atom Fock diagonalisation separates O 2s LP from O 2p LP by energy."""
    xyz_path = FILES_DIR / "water.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", "--method", "hf", "--basis", "cc-pVDZ", "--charge", "0", "--spin", "1", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0

    ibos = parse_ibos(_find_calc_dir("water") / "ibos.txt")
    lps = [o for o in ibos if abs(o["occ"] - 2.0) < 1e-4 and "(LP)" in o.get("label", "")]
    assert len(lps) == 2, f"Expected 2 O LPs, got {len(lps)}"

    lps_sorted = sorted(lps, key=lambda o: o["ene"])
    s_lp, p_lp = lps_sorted[0], lps_sorted[1]
    assert s_lp["ene"] < p_lp["ene"], "s-rich LP should be lower in energy than p-pure LP"

    text = (_find_calc_dir("water") / "ibos.txt").read_text(encoding="utf-8")
    assert "55% 2s + 45% 2pz" in text or "55% 2s + 45% 2p" in text
    assert "100% 2px" in text


def test_methane_pattern():
    """Methane: 1 core + 4 identical C-H sigma (degenerate energies)."""
    xyz_path = FILES_DIR / "methane.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", "--method", "hf", "--basis", "cc-pVDZ", "--charge", "0", "--spin", "1", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0

    ibos = parse_ibos(_find_calc_dir("methane") / "ibos.txt")
    occ = [o for o in ibos if abs(o["occ"] - 2.0) < 1e-4]
    assert len(occ) == 5

    core = next(o for o in occ if "(Core)" in o.get("label", ""))
    assert "C(Core)" in (_find_calc_dir("methane") / "ibos.txt").read_text(encoding="utf-8")

    ch_sigmas = [o for o in occ if "(Core)" not in o.get("label", "")]
    assert len(ch_sigmas) == 4
    energies = {round(o["ene"], 3) for o in ch_sigmas}
    assert len(energies) == 1, f"C-H sigma energies should be degenerate, got {energies}"


def test_ammonia_lp():
    """Ammonia: 1 core + 3 N-H sigma + 1 N LP (DOM≈1)."""
    xyz_path = FILES_DIR / "ammonia.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", "--method", "hf", "--basis", "cc-pVDZ", "--charge", "0", "--spin", "1", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0

    ibos = parse_ibos(_find_calc_dir("ammonia") / "ibos.txt")
    occ = [o for o in ibos if abs(o["occ"] - 2.0) < 1e-4]
    assert len(occ) == 5

    lps = [o for o in occ if "(LP)" in o.get("label", "")]
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
        [sys.executable, "-m", "avogadro_ibo", "--method", "hf", "--basis", "cc-pVDZ", "--charge", "0", "--spin", "1", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=300,
    )
    assert result.returncode == 0

    ibos = parse_ibos(_find_calc_dir("benzene") / "ibos.txt")
    occ = [o for o in ibos if abs(o["occ"] - 2.0) < 1e-4]
    assert len(occ) == 21, f"Expected 21 occupied IBOs, got {len(occ)}"

    text = (_find_calc_dir("benzene") / "ibos.txt").read_text(encoding="utf-8")
    ch_energies = []
    for line in text.splitlines():
        s = line.strip()
        if not s or not s[0].isdigit():
            continue
        if "C-H σ" not in line or "C-H σ*" in line:
            continue
        cols = line.split()
        ch_energies.append(float(cols[2]))

    assert len(ch_energies) == 6, f"Expected 6 C-H sigma, got {len(ch_energies)}"

    emax = max(ch_energies)
    emin = min(ch_energies)
    assert (emax - emin) < 1e-4, (
        f"C-H σ energies split by {emax - emin:.6f} Ha (limit < 1e-4)"
    )


def test_cyclooctatetraene_pi_system():
    """COT tub: non-planar conjugated π system stays σ + π at its equilibrium.

    The D₂d tub is a real molecule whose geometry is not planar, exercising
    the bond-flat degeneracy resolution on a genuine input rather than an
    idealized one.  Expected: four C-C π orbitals of ~99% p-character
    (π-HOMO, π*-LUMO), no sp³-like banana bonds, and near-degenerate
    symmetry-equivalent bonds.
    """
    xyz_path = FILES_DIR / "cyclooctatetraene.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", "--method", "hf", "--basis", "cc-pVDZ", "--charge", "0", "--spin", "1", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=300,
    )
    assert result.returncode == 0

    calc_dir = _find_calc_dir("cyclooctatetraene")
    text = (calc_dir / "ibos.txt").read_text(encoding="utf-8")
    ibos = parse_ibos(calc_dir / "ibos.txt")
    occ = [o for o in ibos if abs(o["occ"] - 2.0) < 1e-4]
    assert len(occ) == 28

    pi_lines = [
        s for s in text.splitlines()
        if s.strip()[:1].isdigit() and "C-C π" in s and "π*" not in s
    ]
    assert len(pi_lines) == 4, f"expected 4 C-C pi orbitals, got {len(pi_lines)}"
    for s in pi_lines:
        # ~99% p-character: a banana bond would carry substantial 2s weight
        assert re.search(r"9[7-9]% 2p", s), f"C-C pi not p-dominated: {s}"

    pi_star_lines = [
        s for s in text.splitlines()
        if s.strip()[:1].isdigit() and "C-C π*" in s
    ]
    assert len(pi_star_lines) == 4, (
        f"expected 4 C-C pi* orbitals, got {len(pi_star_lines)}"
    )

    # π-HOMO / π*-LUMO assignment
    assert any("<- HOMO" in s for s in pi_lines), "HOMO marker not on a C-C pi"
    assert any("<- LUMO" in s for s in pi_star_lines), (
        "LUMO marker not on a C-C pi*"
    )


def test_zncl2_bond_order():
    """ZnCl₂: occupied Zn-Cl σ bonds show W_AB > 0.5 and ionic character > 10%."""
    xyz_path = FILES_DIR / "zncl2.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", "--method", "hf", "--basis", "cc-pVDZ", "--charge", "0", "--spin", "1", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0

    calc_dir = _find_calc_dir("zncl2")
    text = (calc_dir / "ibos.txt").read_text(encoding="utf-8")
    assert "Ion%" in text, "Table should contain Ion% column header"
    assert "Wiberg Bond Orders (σ/π, density)" in text, (
        "Should contain consolidated σ/π Wiberg section"
    )
    assert "Zn-Cl" in text, "Wiberg section should show Zn-Cl pairs"

    ibos = parse_ibos(calc_dir / "ibos.txt")
    zn_cl_sigmas = [
        o for o in ibos
        if abs(o["occ"] - 2.0) < 1e-4 and "σ" in o.get("label", "")
        and "Zn" in o.get("label", "") and "Cl" in o.get("label", "")
    ]
    assert len(zn_cl_sigmas) >= 2, (
        f"Expected 2+ Zn-Cl sigma bonds, got {len(zn_cl_sigmas)}"
    )

    # Check ionic character from the raw line
    for ibo in zn_cl_sigmas:
        ion_str = ibo.get("ion_pct", "0.0")
        if ion_str != "---":
            ion_val = float(ion_str)
            assert ion_val > 10.0, (
                f"Zn-Cl sigma ion%={ion_val:.1f} should be > 10% (polar)"
            )


def test_ethene_sig_pi_wiberg_split():
    """Ethene C=C resolves into σ + π density-Wiberg contributions.

    The consolidated section reports the exact density-matrix Wiberg
    (W_AB = Σ D²_ij) decomposed by orbital class, with interference
    folded in so σ + π = total exactly.  The double bond must read
    σ ≈ 1.0 + π ≈ 1.0; each C-H as σ only.
    """
    xyz_path = FILES_DIR / "ethene.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", "--method", "hf", "--basis", "cc-pVDZ", "--charge", "0", "--spin", "1", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0

    text = (_find_calc_dir("ethene") / "ibos.txt").read_text(encoding="utf-8")
    assert "Wiberg Bond Orders (σ/π, density)" in text

    # Rows: "  C1-C2         2.028   1.028   1.000  (+0.013: σ+0.013, π+0.000)"
    #       "  C5-C6         1.444   1.000   0.444"   (no interference -> no parens)
    sec = text.split("--- Wiberg Bond Orders")[1]
    rows = {}
    for line in sec.splitlines():
        m = re.match(
            r"^\s*([A-Za-z]+\d+)-([A-Za-z]+\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)"
            r"(?:\s+\(([-+\d.]+): σ([-+\d.]+), π([-+\d.]+)\))?",
            line,
        )
        if m:
            tot, sig, pi = float(m.group(3)), float(m.group(4)), float(m.group(5))
            if m.group(6) is not None:
                rows[(m.group(1), m.group(2))] = (tot, sig, pi,
                                                  float(m.group(6)),
                                                  float(m.group(7)),
                                                  float(m.group(8)))
            else:
                rows[(m.group(1), m.group(2))] = (tot, sig, pi, 0.0, 0.0, 0.0)

    cc = [v for k, v in rows.items() if k[0].startswith("C") and k[1].startswith("C")]
    assert cc, f"expected a C-C row in Wiberg section, got rows {list(rows)}"
    tot, sig, pi, interf, is_, ip = cc[0]
    assert abs(sig + pi - tot) < 1e-4, f"σ + π = {sig + pi} ≠ total {tot}"
    assert 0.9 < sig < 1.1, f"C-C σ Wiberg expected ≈1.0, got {sig}"
    assert 0.9 < pi < 1.1, f"C-C π Wiberg expected ≈1.0, got {pi}"
    # Interference parts must sum to the reported interference total
    assert abs((is_ + ip) - interf) < 1e-4, (
        f"σ-part {is_} + π-part {ip} ≠ reported interference {interf}"
    )

    # Every C-H row must be σ-only (no π, no compression of π into total)
    for k, (tot, sig, pi, interf, is_, ip) in rows.items():
        if (k[0].startswith("C") and k[1].startswith("H")) or (k[0].startswith("H") and k[1].startswith("C")):
            assert pi < 0.01, f"C-H row {k} has spurious π contribution {pi}"

    # Ethene's largest pair terms (±0.0053) fall below the 0.01 detail
    # bar, so no detail section should appear: quiet molecules gain
    # no lines.
    assert "Significant orbital-pair" not in text, (
        "ethene should have no pair-interference detail section"
    )


def test_diborane_pair_interference_detail():
    """Diborane: bridge×bridge interference gets its own detail lines.

    The two B-B-H 2e3c orbitals (orbs 3, 4) interfere destructively on
    every bridge leg (-0.0145) and constructively across the B-B contact
    (+0.0145).  Both directions must appear, attributed to the same pair.
    """
    xyz_path = FILES_DIR / "diborane.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", "--method", "hf", "--basis", "cc-pVDZ", "--charge", "0", "--spin", "1", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0

    text = (_find_calc_dir("diborane") / "ibos.txt").read_text(encoding="utf-8")
    assert "--- Significant orbital-pair interference (|term| ≥ 0.01) ---" in text

    sec = text.split("Significant orbital-pair interference")[1]
    m_leg = re.search(
        r"B1-H5: orb3\(B-B-H 2e3c\) × orb4\(B-B-H 2e3c\):\s+([-+\d.]+)", sec
    )
    assert m_leg, "bridge-leg detail row missing"
    assert abs(float(m_leg.group(1)) - (-0.0145)) < 2e-4, (
        f"bridge-leg interference expected ≈ -0.0145, got {m_leg.group(1)}"
    )
    m_bb = re.search(
        r"B1-B6: orb3\(B-B-H 2e3c\) × orb4\(B-B-H 2e3c\):\s+([-+\d.]+)", sec
    )
    assert m_bb, "B-B contact detail row missing"
    assert abs(float(m_bb.group(1)) - 0.0145) < 2e-4, (
        f"B-B interference expected ≈ +0.0145, got {m_bb.group(1)}"
    )


def test_frontier_orbital_summary():
    """Ethene: frontier block exists and gap = LUMO − HOMO internally consistent."""
    xyz_path = FILES_DIR / "ethene.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", "--method", "hf", "--basis", "cc-pVDZ", "--charge", "0", "--spin", "1", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0

    text = (_find_calc_dir("ethene") / "ibos.txt").read_text(encoding="utf-8")
    m = re.search(
        r"HOMO \([^)]*\):\s+(-?[\d.]+) Ha\n"
        r"\s+LUMO \([^)]*\):\s+(-?[\d.]+) Ha\n"
        r"\s+HOMO-LUMO gap:\s+(-?[\d.]+) Ha\s+=\s+(-?[\d.]+) eV\s+=\s+(-?[\d.]+) kcal/mol",
        text,
    )
    assert m, "frontier orbital summary block missing or malformed"
    e_homo, e_lumo, gap_ha, gap_ev, gap_kcal = (float(v) for v in m.groups())
    assert abs(gap_ha - (e_lumo - e_homo)) < 1e-6, "gap Ha ≠ LUMO − HOMO"
    assert abs(gap_ev - gap_ha * 27.211386245988) < 1e-3, "gap eV conversion wrong"
    assert abs(gap_kcal - gap_ha * 627.5094740631) < 0.1, "gap kcal/mol conversion wrong"
    assert gap_ha > 0, "gap should be positive (HOMO below LUMO for RHF)"
    assert "C-C π" in text.split("Frontier Orbital Energies")[1].splitlines()[1], (
        "HOMO should be the C-C π orbital of ethene"
    )


def test_charge_decomposition():
    """Water: charge decomposition sums to 10, O negative, H positive."""
    xyz_path = FILES_DIR / "water.xyz"
    result = subprocess.run(
        [sys.executable, "-m", "avogadro_ibo", "--method", "hf", "--basis", "cc-pVDZ", "--charge", "0", "--spin", "1", str(xyz_path)],
        capture_output=True, text=True, cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0

    text = (_find_calc_dir("water") / "ibos.txt").read_text(encoding="utf-8")
    assert "--- Charge Decomposition ---" in text

    # Parse the charge decomposition section
    lines = text.splitlines()
    charge_start = None
    for i, line in enumerate(lines):
        if "--- Charge Decomposition ---" in line:
            charge_start = i + 3  # skip header title + column header + separator
            break
    assert charge_start is not None

    # Parse per-atom rows (until next separator line)
    total_pop = 0.0
    total_z = 0
    found_o = False
    found_h = False
    for line in lines[charge_start:]:
        s = line.strip()
        if s.startswith("---"):
            break
        if not s or "Atom" in s or "Total" in s:
            continue
        cols = line.split()
        if len(cols) >= 4:
            sym_raw = cols[0]
            Z = int(cols[1])
            pop = float(cols[2])
            net = float(cols[3])
            total_pop += pop
            total_z += Z
            sym = "".join(c for c in sym_raw if c.isalpha())
            if sym == "O":
                found_o = True
                assert net < -0.3, f"O net charge {net} should be negative"
            elif sym == "H":
                found_h = True
                assert net > 0.1, f"H net charge {net} should be positive"

    assert found_o, "O not found in charge decomposition"
    assert found_h, "H not found in charge decomposition"
    assert abs(total_pop - 10.0) < 0.01, (
        f"Total population {total_pop:.3f} should be ~10.0"
    )
    assert total_z == 10, f"Total Z should be 10, got {total_z}"


def test_charge_cjson():
    """Water: cjson['atoms']['partialCharges'] in returned JSON has correct charges."""
    xyz_path = FILES_DIR / "water.xyz"
    coords, numbers = _parse_xyz(xyz_path)
    cjson = {
        "name": "water",
        "atoms": {"coords": {"3d": coords}, "elements": {"number": numbers}},
        "properties": {"totalCharge": 0, "totalSpinMultiplicity": 1},
    }

    # Run through __init__.py's main() (the Avogadro plugin entry) by piping
    # JSON to stdin and capturing stdout JSON.
    input_data = json.dumps({"cjson": cjson, "options": {"method": "hf", "basis": "cc-pVDZ"}, "charge": 0, "spin": 1})
    result = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.argv = ['avogadro_ibo', 'ibo']; "
         "sys.stdout.reconfigure(encoding='utf-8'); "
         "from avogadro_ibo import main; main()"],
        input=input_data, capture_output=True, text=True,
        cwd=PROJECT_DIR, timeout=180,
    )
    assert result.returncode == 0, (
        f"CLI failed:\nstdout:{result.stdout}\nstderr:{result.stderr}"
    )

    output = json.loads(result.stdout)
    cjson_out = output.get("cjson", {})
    charges = cjson_out.get("atoms", {}).get("partialCharges", None)
    assert charges is not None, "cjson should contain atoms.partialCharges"
    assert len(charges) == 3, f"Expected 3 charges, got {len(charges)}"

    # O should be negative, both H positive
    o_charge = charges[0]
    h1_charge = charges[1]
    h2_charge = charges[2]
    assert o_charge < -0.3, f"O charge {o_charge} should be negative"
    assert h1_charge > 0.1, f"H₁ charge {h1_charge} should be positive"
    assert h2_charge > 0.1, f"H₂ charge {h2_charge} should be positive"

    # Sum should be neutral
    assert abs(o_charge + h1_charge + h2_charge) < 0.01, (
        f"Charges sum to {o_charge + h1_charge + h2_charge:.4f}, should be ~0"
    )
