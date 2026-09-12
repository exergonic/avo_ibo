"""IBO analysis table + Wiberg formatting."""


import numpy as np

from .constants import (
    HA_TO_EV,
    HA_TO_KCAL,
    PAIR_DETAIL_THRESH,
    PAIR_NEAR_THRESH,
    _ELEM_SYMBOLS,
)


def _d_spherical_weights(c, atom_idx, atom_of, am_of):
    """
    Compute weights for each spherical d-type from Cartesian d coefficients
    on a given atom.  The 6 Cartesian d-functions in Psi4 (puream=0) are
    ordered: xx, xy, xz, yy, yz, zz.  We project the coefficient vector
    onto the five spherical harmonic directions.

    Returns dict of {name: weight} where weight = squared projection.
    """
    idx = np.where((atom_of == atom_idx) & (am_of == 2))[0]
    if len(idx) < 6:
        return {}
    c = np.asarray(c, dtype=np.float64)
    # Psi4 Cartesian d order: xx, xy, xz, yy, yz, zz
    c_xx, c_xy, c_xz, c_yy, c_yz, c_zz = c[idx[:6]]
    return {
        "dxy": c_xy**2,
        "dxz": c_xz**2,
        "dyz": c_yz**2,
        "dz2": (-c_xx - c_yy + 2 * c_zz) ** 2,
        "dx2y2": (c_xx - c_yy) ** 2,
    }


def _hybrid_str(c, am_of, atom_of, func_n, func_dtype, top_atom):
    """
    Build a specific hybrid label for the dominant atom, e.g.
        "57% 4s + 43% 3dz²"
        "100% 1s"
        "100% 4pz"
        "83% 3s + 17% 3pz"
        "46% 4s + 54% 3d"
    """
    c = np.asarray(c, dtype=np.float64)
    pA = float(np.sum(c[np.where(atom_of == top_atom)] ** 2))
    if pA < 1e-12:
        return ""

    parts = []
    for am_label, am_val in [("s", 0), ("p", 1), ("d", 2)]:
        idx_am = np.where((atom_of == top_atom) & (am_of == am_val))[0]
        if len(idx_am) == 0:
            continue
        total_am = float(np.sum(c[idx_am] ** 2))
        pct = total_am / pA * 100.0
        if pct < 1.0:
            continue

        # Find dominant n within this l-subspace
        n_counts = {}
        for fi in idx_am:
            n_key = func_n[fi]
            n_counts[n_key] = n_counts.get(n_key, 0.0) + c[fi] ** 2
        dominant_n = max(n_counts, key=lambda k: n_counts[k])

        # Determine dominant subtype
        subtype = ""
        if am_val == 1:  # p-orbitals: px, py, pz
            st_counts = {}
            for fi in idx_am:
                st = func_dtype[fi]
                if st:
                    st_counts[st] = st_counts.get(st, 0.0) + c[fi] ** 2
            if st_counts:
                top_st = max(st_counts, key=lambda k: st_counts[k])
                if st_counts[top_st] > 0.5 * total_am:
                    subtype = top_st  # e.g. "pz"
        elif am_val == 2:  # d-orbitals: dxy, dxz, dyz, dz2, dx2y2
            d_weights = _d_spherical_weights(c, top_atom, atom_of, am_of)
            if d_weights:
                top_st = max(d_weights, key=lambda k: d_weights[k])
                if d_weights[top_st] > 0.5 * max(d_weights.values()):
                    subtype = top_st  # e.g. "dz2"

        label = str(dominant_n) + (subtype if subtype else am_label)
        parts.append(f"{pct:.0f}% {label}")

    return " + ".join(parts)


def _wiberg_per_ibo(pop, occ, A, B):
    """
    Per-IBO contribution to the Wiberg bond order between atoms A and B.

    In the orthonormal IAO basis, the density contribution from a single IBO
    with coefficient vector c_k is D^{(k)} = occ_k · c_k c_k^T.  The Wiberg
    index between A and B from this IBO is:

        W_AB^{(k)} = Σ_{i∈A} Σ_{j∈B} (D^{(k)}_ij)²

    which simplifies (by the independence of i and j sums) to:

        W_AB^{(k)} = occ_k² · P_A · P_B

    where P_X = Σ_{i∈X} c_{k,i}² is the Mulliken population on atom X.

    For RHF occupied (occ=2): W_AB = 4 · P_A · P_B, ranging from 0 (pure
    ionic, no shared density) to 1 (pure covalent 2c-2e bond with 50/50
    sharing).
    """
    return float(occ**2 * pop[A] * pop[B])


def _ionic_pct(pop, A, B):
    """
    Percent ionic character between atoms A and B from per-atom populations.

        Ionic% = |P_A - P_B| / (P_A + P_B) × 100

    Ranges from 0% (pure covalent, equal sharing) to 100% (pure ionic,
    all density on one atom).
    """
    num = abs(pop[A] - pop[B])
    den = pop[A] + pop[B]
    return num / den * 100.0 if den > 1e-12 else 0.0


def _classify_orbital(oc, pop, order, top_A, top_B, s_char, p_char, d_char,
                      elem, am_of, atom_of, func_n, c):
    """Return a classification label string for one IAO-basis orbital.

    Classifies occupied orbitals as Core, LP, σ/π bond, 2e3c, or Deloc,
    and virtual orbitals as the corresponding antibond (*) type.

    Thresholds (DOM-based, matching IboView defaults):
      Core:  DOM > 0.99 + s-character > 0.75 on n=1 (1s only; 2s/3s are valence)
      LP:    DOM > 0.90 on one atom
      σ/π:   DOM_shared > 0.75, both atoms carry density (>0.02);
             π if p-fraction > 0.85 on both atoms
      LP-s:  DOM > 0.70, s-character > 0.5 (transitional)
      2e3c:  3rd atom carries >10% density, 4th atom <3%
      Deloc: everything else (multi-atom delocalisation)
      Virtual antibond thresholds follow the same logic with slightly
      looser cut-offs (0.60 vs 0.75 for shared density).
    """
    if oc > 1.5:
        # Determine principal quantum number of the dominant s-contribution
        # on the top atom.  Only n=1 (1s) is a true core; 2s/3s are valence.
        _s_idx = np.where((atom_of == top_A) & (am_of == 0))[0]
        _s_n = int(func_n[_s_idx[np.argmax(c[_s_idx] ** 2)]]) if len(_s_idx) else 0
        if pop[top_A] > 0.99 and s_char > 0.75 and _s_n == 1:
            return f"{elem_symbol(elem[top_A])}(Core)"
        elif pop[top_A] > 0.90:
            return f"{elem_symbol(elem[top_A])}(LP)"
        elif pop[top_A] + pop[top_B] > 0.75 and pop[top_B] > 0.02:
            pfrac_A = _p_frac(c, am_of, top_A, atom_of)
            pfrac_B = _p_frac(c, am_of, top_B, atom_of)
            bond_type = "π" if (pfrac_A > 0.85 and pfrac_B > 0.85) else "σ"
            a, b = sorted([top_A, top_B])
            symA = elem_symbol(elem[a])
            symB = elem_symbol(elem[b])
            return f"{symA}-{symB} {bond_type}"
        elif pop[top_A] > 0.70:
            symA = elem_symbol(elem[top_A])
            if s_char > 0.5:
                return f"{symA}(LP-s)"
            else:
                return f"{symA}(LP)"
        else:
            if len(order) >= 3 and pop[order[2]] > 0.10 and (len(order) < 4 or pop[order[3]] <= 0.03):
                atoms = sorted(
                    [order[0], order[1], order[2]],
                    key=lambda i: (elem_symbol(elem[i]), i),
                )
                syms = "-".join(elem_symbol(elem[a]) for a in atoms)
                return f"{syms} 2e3c"
            else:
                return "Deloc"
    else:
        if pop[top_A] + pop[top_B] > 0.75 and pop[top_B] > 0.02:
            pfrac_A = _p_frac(c, am_of, top_A, atom_of)
            pfrac_B = _p_frac(c, am_of, top_B, atom_of)
            bond_type = "π" if (pfrac_A > 0.85 and pfrac_B > 0.85) else "σ"
            a, b = sorted([top_A, top_B])
            symA = elem_symbol(elem[a])
            symB = elem_symbol(elem[b])
            return f"{symA}-{symB} {bond_type}*"
        elif len(order) >= 3 and pop[order[2]] > 0.08:
            pfrac_top = _p_frac(c, am_of, top_A, atom_of)
            symA = elem_symbol(elem[top_A])
            return f"{symA} π*" if pfrac_top > 0.85 else f"{symA} anti*"
        elif pop[top_A] + pop[top_B] > 0.60 and pop[top_B] > 0.02:
            a, b = sorted([top_A, top_B])
            symA = elem_symbol(elem[a])
            symB = elem_symbol(elem[b])
            return f"{symA}-{symB} anti*"
        elif pop[top_A] > 0.50:
            return f"{elem_symbol(elem[top_A])}(virt)"
        else:
            return "Virt"


def analyze_ibos(
    C_IAO_all,
    occ_all,
    energies_all,
    nocc,
    atom_of,
    am_of,
    func_n,
    func_dtype,
    elem,
    method,
    basis,
    ref,
    mol_name="",
):
    """
    Build a formatted IBO analysis table covering all IAO-basis orbitals.

    For each orbital (occupied IBO or valence-virtual IAO), compute:
      - per-atom populations from IAO coefficients
      - DOM (largest two n_A fractions summed)
      - per-IBO Wiberg bond order (W_AB) and percent ionic character
      - specific nl/subtype hybrid label on the dominant atom
    """
    n_IAO, n_orb = C_IAO_all.shape
    n_atoms = len(elem)

    lines = []
    orbid_labels = [""] * n_orb
    atom_pop = np.zeros(
        n_atoms, dtype=np.float64
    )  # accumulated per-atom electron counts
    lines.append(f"IBO Analysis: {mol_name}  ({method}/{basis}, {ref.upper()})")
    lines.append("")

    # Pre-pass to size the Composition column to fit the widest entry
    comp_width = len("Composition")
    for orb in range(n_orb):
        sq = C_IAO_all[:, orb] ** 2
        pop = np.zeros(n_atoms, dtype=np.float64)
        np.add.at(pop, atom_of, sq)
        order = np.argsort(-pop)
        comp_parts = []
        for A in order[:4]:
            if pop[A] > 0.005:
                sym = f"{elem_symbol(elem[A])}{A + 1}"
                pct = pop[A] * 100.0
                comp_parts.append(f"{sym}({pct:.1f}%)")
        comp = " + ".join(comp_parts)
        comp_width = max(comp_width, len(comp))
    comp_width += 2

    header = (
        f"  {'#':>3}  {'Occ':>7}  {'Energy':>10}  "
        f"{'Type':>16}  {{:<{comp_width}}}  {'Hybrid':<22}  "
        f"{'Ion%':>5}  {'H/L':>9}"
    ).format("Composition")
    lines.append(header)
    lines.append("-" * len(header))

    # Identify degenerate manifolds: groups of consecutive orbitals
    # with energy differences < 1e-4 Ha.  Within symmetric molecules,
    # the PM functional leaves symmetry-equivalent bonds with small
    # residual energy splittings (~1e-5 Ha for benzene C-H σ).
    DEG_THRESH = 2e-4  # ~0.13 kcal/mol; catches all PM convergence noise
    deg_ranges = []
    is_degen = np.zeros(n_orb, dtype=bool)
    group_start = None
    for i in range(n_orb):
        if i > 0 and abs(energies_all[i] - energies_all[i - 1]) < DEG_THRESH:
            if group_start is None:
                group_start = i - 1
            is_degen[i] = True
        else:
            if group_start is not None and i - group_start > 1:
                deg_ranges.append((group_start, i))
                for j in range(group_start, i):
                    is_degen[j] = True
            group_start = None
    if group_start is not None and n_orb - group_start > 1:
        deg_ranges.append((group_start, n_orb))
        for j in range(group_start, n_orb):
            is_degen[j] = True

    # HOMO/LUMO by occupancy, not energy-sort rank.  energies_all is
    # ascending but not a guaranteed occupied-then-virtual prefix
    # (see IBOResult.n_occ).  Rank markers (orb == nocc-1 / nocc) would
    # lie if a valence virtual drops below an occupied IBO energy.
    occ_idx = np.where(occ_all > 1.5)[0]
    vir_idx = np.where(occ_all < 0.5)[0]
    homo_i = None
    lumo_i = None
    if len(occ_idx):
        homo_i = int(occ_idx[np.argmax(energies_all[occ_idx])])
    if len(vir_idx):
        lumo_i = int(vir_idx[np.argmin(energies_all[vir_idx])])

    for orb in range(n_orb):
        oc = occ_all[orb]
        sq = C_IAO_all[:, orb] ** 2
        pop = np.zeros(n_atoms, dtype=np.float64)
        np.add.at(pop, atom_of, sq)
        if oc > 1.5:
            # RHF-only: occ=2.0 per occupied orbital.  UHF would need
            # separate alpha/beta occupancy arrays — unreachable due
            # to the closed-shell guard in compute_ibo.
            atom_pop += pop * oc  # accumulate electron count per atom

        # Dominant atom and its population
        order = np.argsort(-pop)
        top_A = order[0]
        top_B = order[1]

        # s/p/d breakdown on the dominant atom (single pass)
        s_char, p_char, d_char = _spd_frac(C_IAO_all[:, orb], am_of, top_A, atom_of)

        # Determine orbital type
        orbid = _classify_orbital(oc, pop, order, top_A, top_B, s_char, p_char,
                                  d_char, elem, am_of, atom_of, func_n,
                                  C_IAO_all[:, orb])

        orbid_labels[orb] = orbid

        comp_parts = []
        for A in order[:4]:
            if pop[A] > 0.005:
                sym = f"{elem_symbol(elem[A])}{A + 1}"
                pct = pop[A] * 100.0
                comp_parts.append(f"{sym}({pct:.1f}%)")
        comp = " + ".join(comp_parts)

        # For hybrid label, prefer a non-H atom when it carries meaningful density
        hybrid_atom = top_A
        if elem[top_A] == 1:
            for A in order[1:]:
                if elem[A] != 1 and pop[A] > 0.02:
                    hybrid_atom = A
                    break
        hybrid = _hybrid_str(
            C_IAO_all[:, orb], am_of, atom_of, func_n, func_dtype, hybrid_atom
        )
        if hybrid_atom != top_A:
            hybrid = f"{elem_symbol(elem[hybrid_atom])}: {hybrid}"

        # Per-IBO Wiberg bond order and ionic character (between top_A, top_B)
        w_ab = _wiberg_per_ibo(pop, oc, top_A, top_B)
        if oc > 1.5 and w_ab > 0.001:
            ion_str = f"{_ionic_pct(pop, top_A, top_B):.1f}"
        else:
            ion_str = "---"

        hl = ""
        if homo_i is not None and orb == homo_i:
            hl = "<- HOMO"
        elif lumo_i is not None and orb == lumo_i:
            hl = "<- LUMO"
        degen_tag = " †" if is_degen[orb] else ""
        lines.append(
            f"  {orb + 1:>3d}  {oc:>7.3f}  {energies_all[orb]:>10.6f}  "
            f"{orbid:>{16 - len(degen_tag)}}{degen_tag}  "
            f"{comp:<{comp_width}}  {hybrid:<22}  "
            f"{ion_str:>5}  {hl:>9}"
        )

    # Footnote for degenerate manifolds
    if deg_ranges:
        deg_groups = []
        for start, end in deg_ranges:
            deg_groups.append(f"{start + 1}-{end}")
        lines.append(
            f"  † Orbitals {', '.join(deg_groups)} form degenerate manifolds (ΔE < {DEG_THRESH:.0e} Ha).\n"
            f"  Small energy differences within each manifold are PM convergence noise and do not indicate true energy splittings."
        )

    lines.append("")
    # RHF-only: total = 2 × nocc.  The UHF branch (else nocc) is
    # unreachable — closed-shell guard prevents open-shell in compute_ibo.
    lines.append(f"Total electrons: {int(2 * nocc) if ref == 'rhf' else nocc}")

    # Frontier block uses the same occupancy-based indices as the table H/L
    # markers above, so the arrows and the summary cannot disagree.
    if homo_i is not None and lumo_i is not None:
        gap_ha = energies_all[lumo_i] - energies_all[homo_i]
        gap_ev = gap_ha * HA_TO_EV
        gap_kcal = gap_ha * HA_TO_KCAL
        lines.append("")
        lines.append("--- Frontier Orbital Energies ---")
        lines.append(
            f"  HOMO ({orbid_labels[homo_i]}, orb {homo_i + 1}): {energies_all[homo_i]:>10.6f} Ha"
        )
        lines.append(
            f"  LUMO ({orbid_labels[lumo_i]}, orb {lumo_i + 1}): {energies_all[lumo_i]:>10.6f} Ha"
        )
        lines.append(
            f"  HOMO-LUMO gap: {gap_ha:>10.6f} Ha = {gap_ev:>7.3f} eV = {gap_kcal:>8.1f} kcal/mol"
        )

    charge_section = _format_charge_decomposition(atom_pop, elem)
    lines.append(charge_section)

    net_charges = [float(int(round(elem[A])) - atom_pop[A]) for A in range(n_atoms)]
    return "\n".join(lines), orbid_labels, net_charges




def format_wiberg(C_IAO_occ, atom_of, am_of, elem, labels=None):
    """Single Wiberg table: exact density-matrix totals decomposed σ/π.

    The density Wiberg (as originally reported) is

        W_AB = Σ_{i∈A,j∈B} D²_ij,   D = 2·C_occ·C_occᵀ  (RHF)

    Expanding the square over orbitals gives an exact sum over orbital
    pairs,

        W_AB = 4 Σ_{k,l} G^A_kl G^B_kl,   G^A_kl = Σ_{i∈A} c_ki c_li

    with diagonal terms (k=l) the per-IBO shares occ²·P_A·P_B and
    off-diagonal terms the inter-orbital interference.  Orbitals are
    classed σ or π purely by p-fraction (p > 0.85 on both dominant
    atoms → π, no population gate), and the interference is folded
    into its class: (σ,σ) → σ, (π,π) → π, (σ,π) → split 50/50.
    Total = σ + π exactly; the folded interference is echoed in a
    parenthesised column for transparency.

    A follow-on detail section lists individual orbital-pair
    interference terms with |term| >= PAIR_DETAIL_THRESH, grouped by
    bond in table order, so the reader can see which orbitals drive a
    bond's interference.  ``labels`` supplies the per-orbital table
    labels (column k ↔ analysis-table row k+1); when omitted, bare orb
    numbers are shown.  Bonds whose every pair term falls below the
    threshold get no detail lines, and the section is omitted entirely
    when empty.
    """
    n_occ = C_IAO_occ.shape[1]
    n_atoms = len(elem)

    # Classify each occupied orbital σ/π (p-fraction rule only).
    is_pi = np.zeros(n_occ, dtype=bool)
    for k in range(n_occ):
        c = C_IAO_occ[:, k]
        sq = c**2
        pop = np.zeros(n_atoms, dtype=np.float64)
        np.add.at(pop, atom_of, sq)
        order = np.argsort(-pop)
        A, B = int(order[0]), int(order[1])
        pa = _p_frac(c, am_of, A, atom_of)
        pb = _p_frac(c, am_of, B, atom_of)
        is_pi[k] = pa > 0.85 and pb > 0.85

    # Per-atom orbital-pair overlap blocks G^A_kl.
    G = np.zeros((n_atoms, n_occ, n_occ))
    for a in range(n_atoms):
        Ca = C_IAO_occ[atom_of == a, :]
        G[a] = Ca.T @ Ca

    sigma = np.zeros((n_atoms, n_atoms), dtype=np.float64)
    pi = np.zeros((n_atoms, n_atoms), dtype=np.float64)
    # Folded interference, tracked per class: int_sigma = σσ + ½σπ,
    # int_pi = ππ + ½σπ.  These are what the σ/π columns actually
    # contain beyond their diagonal shares (4·P_A·P_B).
    int_sigma = np.zeros((n_atoms, n_atoms), dtype=np.float64)
    int_pi = np.zeros((n_atoms, n_atoms), dtype=np.float64)
    # Significant off-diagonal terms, kept for the detail section:
    # (A, B) with A < B -> [(k, l, value)] with |value| >= threshold.
    pair_detail = {}
    pair_near = {}

    for k in range(n_occ):
        # Diagonal (per-orbital share): 4·P_A·P_B
        contrib = 4.0 * np.einsum("a,b->ab", G[:, k, k], G[:, k, k])
        if is_pi[k]:
            pi += contrib
        else:
            sigma += contrib
        for l in range(k + 1, n_occ):
            # Off-diagonal (interference): 8·G^A_kl·G^B_kl
            contrib = 8.0 * np.einsum("a,b->ab", G[:, k, l], G[:, k, l])
            if is_pi[k] and is_pi[l]:
                pi += contrib
                int_pi += contrib
            elif is_pi[k] or is_pi[l]:
                # σ-π cross term: split evenly (contribution is symmetric)
                half = 0.5 * contrib
                sigma += half
                pi += half
                int_sigma += half
                int_pi += half
            else:
                sigma += contrib
                int_sigma += contrib
            # Record significant pair terms for the detail section, plus
            # near-miss terms for the closing footnote.
            big = np.abs(contrib) >= PAIR_DETAIL_THRESH
            if big.any():
                for A in range(n_atoms):
                    for B in range(A + 1, n_atoms):
                        if big[A, B]:
                            pair_detail.setdefault((A, B), []).append(
                                (k, l, float(contrib[A, B]))
                            )
            near = (np.abs(contrib) >= PAIR_NEAR_THRESH) & ~big
            if near.any():
                for A in range(n_atoms):
                    for B in range(A + 1, n_atoms):
                        if near[A, B]:
                            pair_near.setdefault((A, B), []).append(
                                (k, l, float(contrib[A, B]))
                            )

    rows = []
    total = sigma + pi
    for A in range(n_atoms):
        for B in range(A + 1, n_atoms):
            if total[A, B] > 0.01:
                symA = elem_symbol(elem[A])
                symB = elem_symbol(elem[B])
                rows.append(
                    (symA, A, symB, B, total[A, B], sigma[A, B], pi[A, B],
                     int_sigma[A, B], int_pi[A, B])
                )

    if not rows:
        return ""

    lines = [
        "",
        "",
        "--- Wiberg Bond Orders (σ/π, density) ---",
        "  W_AB = Σ_{i∈A,j∈B} D²_ij (density Wiberg); σ + π = total exactly.",
        "  σ, π columns include their class's folded interference; the",
        "  parenthetical reports the same interference as (σ-part, π-part).",
        f"  {'Bond':<10}{'Total':>8}{'σ':>8}{'π':>8}{'  (interference)':>28}",
    ]
    for symA, a, symB, b, t, s, p, is_, ip in sorted(rows, key=lambda x: -x[4]):
        # Kill floating-point -0.000 noise in the σ/π columns (the
        # interference parts keep their genuine sign).
        s_disp = 0.0 if abs(s) < 5e-4 else s
        p_disp = 0.0 if abs(p) < 5e-4 else p
        if abs(is_) + abs(ip) >= 5e-4:
            # Both parts below the noise floor -> omit the parenthetical
            # entirely; the row is pure diagonal shares.
            is_disp = 0.0 if abs(is_) < 5e-4 else is_
            ip_disp = 0.0 if abs(ip) < 5e-4 else ip
            lines.append(
                f"  {symA}{a+1}-{symB}{b+1:<7}{t:>8.3f}{s_disp:>8.3f}{p_disp:>8.3f}"
                f"  ({is_disp+ip_disp:+.3f}: σ{is_disp:+.3f}, π{ip_disp:+.3f})"
            )
        else:
            lines.append(
                f"  {symA}{a+1}-{symB}{b+1:<7}{t:>8.3f}{s_disp:>8.3f}{p_disp:>8.3f}"
            )

    # Detail section: significant orbital-pair interference terms, grouped
    # by bond in table order.  Omitted entirely when nothing clears the bar
    # (ordinary single/double bonds), so quiet molecules gain no lines.
    def _orb_tag(k):
        if labels is not None and k < len(labels) and labels[k]:
            return f"orb{k+1}({labels[k]})"
        return f"orb{k+1}"

    detail = []
    for symA, a, symB, b, t, s, p, is_, ip in sorted(rows, key=lambda x: -x[4]):
        terms = pair_detail.get((a, b))
        if not terms:
            continue
        for k, l, v in sorted(terms, key=lambda t_: -abs(t_[2])):
            detail.append(
                f"  {symA}{a+1}-{symB}{b+1}: {_orb_tag(k)} × {_orb_tag(l)}: {v:+.4f}"
            )
    if detail:
        lines.append("")
        lines.append(
            "--- Significant orbital-pair interference "
            f"(|term| ≥ {PAIR_DETAIL_THRESH:g}) ---"
        )
        lines.append("  Positive pair terms add to the bond order; negative pair")
        lines.append("  terms subtract from it. Signs are relative to the listed")
        lines.append("  bond — the same pair may contribute oppositely elsewhere.")
        lines.append("  Orbitals are numbered as in the analysis table above.")
        lines.extend(detail)
        near_all = [
            (symA, a, symB, b, k, l, v)
            for symA, a, symB, b, t, s, p, is_, ip
            in sorted(rows, key=lambda x: -x[4])
            for k, l, v in pair_near.get((a, b), [])
        ]
        if near_all:
            n_near = len(near_all)
            sA, a, sB, b, k, l, v = max(near_all, key=lambda e: abs(e[6]))
            noun = "term" if n_near == 1 else "terms"
            lines.append(
                f"  ({n_near} {noun} in [{PAIR_NEAR_THRESH:g}, "
                f"{PAIR_DETAIL_THRESH:g}) omitted; largest: "
                f"{sA}{a+1}-{sB}{b+1}: {_orb_tag(k)} × {_orb_tag(l)} = {v:+.4f})"
            )
    return "\n".join(lines)


def _format_charge_decomposition(atom_pop, elem):
    """
    Format a charge decomposition table from accumulated IAO populations.

    For each atom A:
        Q_A = Σ_k occ_k · P_A^{(k)}   (total electrons on atom A)
        Net charge = Z_A - Q_A

    ``atom_pop[A]`` is Q_A as a float.  ``elem`` gives atomic numbers.
    RHF-only: assumes occ=2.0 per occupied orbital.  UHF would need
    separate alpha/beta populations — unreachable due to the
    closed-shell guard in compute_ibo.
    """
    lines = ["", "--- Charge Decomposition ---"]
    header = f"  {'Atom':>5}  {'Z':>3}  {'Pop':>8}  {'Net Charge':>10}"
    lines.append(header)
    lines.append("-" * len(header))
    total_pop = 0.0
    total_z = 0
    for A in range(len(elem)):
        Z = int(round(elem[A]))
        pop = atom_pop[A]
        net = Z - pop
        sym = elem_symbol(Z)
        lines.append(f"  {sym}{A+1:<3}  {Z:>3d}  {pop:>8.3f}  {net:>+10.3f}")
        total_pop += pop
        total_z += Z
    lines.append("-" * len(header))
    lines.append(
        f"Total:  {total_z:>3d}  {total_pop:>8.3f}  {total_z - total_pop:>+10.3f}"
    )
    return "\n".join(lines)


def elem_symbol(Z):
    Z = int(round(Z))
    if Z < len(_ELEM_SYMBOLS):
        return _ELEM_SYMBOLS[Z]
    return f"E{Z}"


def _spd_frac(c, am_of, atom, atom_of):
    """Return (s_char, p_char, d_char) for the given atom in one pass."""
    idx = np.where(atom_of == atom)[0]
    if len(idx) == 0:
        return 0.0, 0.0, 0.0
    am_atom = am_of[idx]
    c_atom = c[idx] ** 2
    s = float(np.sum(c_atom[am_atom == 0])) if np.any(am_atom == 0) else 0.0
    p = float(np.sum(c_atom[am_atom == 1])) if np.any(am_atom == 1) else 0.0
    d = float(np.sum(c_atom[am_atom == 2])) if np.any(am_atom == 2) else 0.0
    return s, p, d


def _p_frac(c, am_of, atom, atom_of):
    """p/s/d ratio on the given atom; 0 if no density."""
    s, p, d = _spd_frac(c, am_of, atom, atom_of)
    total = s + p + d
    return p / total if total > 0 else 0.0
