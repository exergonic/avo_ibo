#set page(margin: 1in)
#set text(font: "New Computer Modern", size: 11pt)

 Mathematics of the avogadro-ibo Plugin

This document walks through every computational step in
#link("src/avogadro_ibo/calcs.py")[`calcs.py`], connecting each operation to the
equations in Knizia, *J. Chem. Theory Comput.* *2013*, *9*, 4834–4843
(#link("https://doi.org/10.1021/ct400687b")[DOI: 10.1021/ct400687b]).

#line()

= 1 Notation

Let the *AO basis* (the SCF basis, e.g. cc-pVDZ) have $n_("AO")$
functions ${|chi_mu⟩}$, and the *minimal basis* (STO-3G,
$n_("min")$) have functions ${|tilde(chi)_alpha⟩}$.

#table(columns: (auto, auto, auto),
  [*Symbol*], [*Shape*], [*Meaning*],
  [$bold(S)$], [$n_("AO") × n_("AO")$], [AO overlap, $S_(mu nu) = ⟨chi_(mu)|chi_(nu)⟩$],
  [$bold(S)^(12)$], [$n_("AO") × n_("min")$], [Cross-overlap, $S^(12)_(mu alpha) = ⟨chi_(mu)|tilde(chi)_(alpha)⟩$],
  [$bold(S)^("min")$], [$n_("min") × n_("min")$], [Minimal-basis overlap, $S^("min")_(alpha b eta) = ⟨tilde(chi)_(alpha)|tilde(chi)_(b eta)⟩$],
  [$bold(C)^("occ")$], [$n_("AO") × n_("occ")$], [Occupied MO coefficients (columns = orbitals)],
  [$bold(C)^("IAO")$], [$n_("AO") × n_("min")$], [IAO coefficients in the AO basis],
  [$bold(C)^("IAO,occ")$], [$n_("min") × n_("occ")$], [Occupied MO coefficients in the IAO basis],
  [$bold(F)^("AO")$], [$n_("AO") × n_("AO")$], [Fock matrix in the AO basis],
  [$bold(F)^("IAO")$], [$n_("min") × n_("min")$], [Fock matrix in the IAO basis],
  [$n_A(i)$], [$RR$], [Mulliken population of orbital $i$ on atom $A$],
  [$delta_i$], [$RR$], [*D*iatomic *O*ccupation of orbital $i$],
)

Indices run over:
- $mu,nu$ — AO functions ($1 .. n_("AO")$)
- $alpha,beta$ — minimal-basis functions ($1 .. n_("min")$)
- $i,j$ — occupied orbitals ($1 .. n_("occ")$)
- $A,B$ — atoms ($1 .. n_("atom")$)

#line()

= 2 The SCF wavefunction

Psi4 solves the Hartree–Fock (or DFT) equations for the molecular geometry
from Avogadro:

$
bold(F)^("AO") bold(C) = bold(S) bold(C) bold(epsilon)
$

The default method is HF with the cc-pVDZ basis and Cartesian
($"puream"=0$) functions.  The occupied block
$bold(C)^("occ") in RR^(n_"AO") × n_("occ")$ is
extracted from the converged solution (`calcs.py:509`).

#line()

= 3 Intrinsic Atomic Orbitals (IAOs)

IAOs form an orthonormal, chemically balanced basis that exactly spans the
occupied space while keeping the atomic character of the minimal basis.
The construction follows *Appendix C* of the paper (IboView's
`MakeIaoBasisNew`).

== 3.1 Project the minimal basis into the full AO space

Each minimal-basis function $|tilde(chi)_alpha⟩$ is expanded in
the AO basis via the overlap:

$
bold(P)^(12) = bold(S)^(-1) bold(S)^(12),
quad quad
|tilde(chi)_alpha^("(AO)")⟩ = sum_(mu=1)^(n_"AO")
(bold(P)^(12))_(mu alpha)  |chi_mu⟩ .
$

`calcs.py:64–65`

== 3.2 Express occupied MOs in the minimal basis

The projection of each occupied orbital onto the minimal-basis functions
gives the minimal-basis representation:

$
(bold(C)^("occ,min"))_(alpha i)
= ⟨tilde(chi)_alpha|psi_i⟩
= sum_(mu) S^(12)_(mu alpha) C^("occ")_(mu i)
= (bold(S)^(12 T) bold(C)^("occ"))_(alpha i).
$

`calcs.py:68`

== 3.3 Depolarisation — solve $bold(S)^("min") tilde(bold(C)) = bold(C)^("occ,min")$

We find the coefficients $tilde(bold(C))$ that express the occupied
orbitals in the minimal basis:

$
bold(S)^("min") tilde(bold(C)) = bold(C)^("occ,min")
quad=>quad
tilde(bold(C)) = (bold(S)^("min"))^(-1) bold(C)^("occ,min").
$

$(bold(S)^("min"))^(-1)$ is obtained via Cholesky factorisation.
`calcs.py:70–72`

== 3.4 Occupied-space metric

The overlap of the projected occupied orbitals within the minimal basis is:

$
tilde(bold(S)) = (bold(C)^("occ,min"))^(T) tilde(bold(C))
= (bold(C)^("occ,min"))^(T) (bold(S)^("min"))^(-1)
  bold(C)^("occ,min").
$

`calcs.py:75`

== 3.5 Inverse metric — $tilde(bold(C))^((2))$

We solve for the coefficients $tilde(bold(C))^((2))$:

$
tilde(bold(S)) tilde(bold(C))^((2) T) = tilde(bold(C))^(T)
quad=>quad
tilde(bold(C))^((2)) = tilde(bold(C)) tilde(bold(S))^(-1).
$

`calcs.py:78–80`

== 3.6 Tight residual — $bold(T)^((4))$

The component of each occupied orbital that lies *orthogonal* to the
minimal-basis span is the tight residual:

$
bold(T)^((4)) = bold(C)^("occ")
   - bold(P)^(12) tilde(bold(C))^((2)).
$

`calcs.py:83`

The columns of $bold(T)^((4))$ are the occupied-space vectors that
the minimal basis alone cannot describe.

== 3.7 Construct the IAO coefficients

Each IAO is the corresponding projected minimal-basis function plus the
tight residual weighted by the overlap of that minimal-basis function with
the occupied orbitals:

$
bold(C)^("IAO") = bold(P)^(12)
   + bold(T)^((4)) (bold(C)^("occ,min"))^(T).
$

In matrix form: $(bold(C)^("IAO"))_(mu alpha)$
is the coefficient of AO function $|chi_mu⟩$ in IAO $|phi_alpha⟩$.

`calcs.py:86`

== 3.8 Symmetric (Löwdin) orthogonalisation

The IAOs are not yet orthonormal.  Define the metric of the IAO basis:

$
bold(M) = (bold(C)^("IAO"))^(T) bold(S) bold(C)^("IAO").
$

Diagonalise $bold(M) = bold(U) bold(m) bold(U)^(T)$ and apply
$bold(M)^(-1/2) = bold(U) bold(m)^(-1/2) bold(U)^(T)$:

$
bold(C)^("IAO") arrow
   bold(C)^("IAO")   bold(M)^(-1/2).
$

Now $⟨phi_(alpha)|phi_(b eta)⟩ = (bold(C)^("IAO"))^(T) bold(S) bold(C)^("IAO") = bold(I)$.

`calcs.py:90–92`

== 3.9 Express occupied orbitals in the IAO basis

The occupied MO coefficients in the orthonormal IAO basis are:

$
bold(C)^("IAO,occ") = (bold(C)^("IAO"))^(T)
                           bold(S) bold(C)^("occ").
$

Because IAOs span the occupied space (by construction), this is a
unitary rotation and $bold(C)^("IAO") bold(C)^("IAO,occ") = bold(C)^("occ")$ is exact.

`calcs.py:97`

#line()

= 4 Pipek–Mezey localisation in the IAO basis

The IBOs are obtained by maximising the generalised
Pipek–Mezey functional (*Eq. 4* of the paper):

$
L = sum_(A=1)^(n_"atom") sum_(i=1)^(n_"occ")
    [n_A(i)]^(4),
quad quad
n_A(i) = sum_(alpha in A) |C^("IAO,occ")_(alpha i)|^(2).
$

Here $n_A(i)$ is the Mulliken population of orbital $i$ on atom $A$,
computed directly from the squared IAO coefficients (the IAO basis is
orthonormal).

== 4.1 Initial guess: canonical MOs (no symmetry breaking)

The occupied block $bold(C)^("IAO,occ")$ is the identity in
IAO space — it already diagonalises $bold(F)^("IAO")$ and
provides a unique, deterministic starting point.  *No random
perturbation is applied before localisation.*

IboView applies an 18° Cayley rotation (`RotateVectorsRandomly`) to
break accidental near-degeneracies, but in our testing this *increased*
the benzene C-H σ split from $2×10^(-5)$ Ha to $1.2×10^(-4)$ Ha
(Gotcha 20 in AGENTS.md).  The symmetric converged solution is the
nearest local maximum for the fixed sequential sweep order
($i=1..n_"occ", j&lt;i$) — no perturbation improves it.

All PAO-like methods that use the PM functional on an orthogonal
minimal basis converge to the same nearest local maximum without
symmetry breaking, as long as the sweep order is fixed and
reproducible.  The line-reference at the bottom of this section
is therefore a stub; the actual localisation code begins at §4.2.

`calcs.py:147–148`

== 4.2 p=2 warm-start

The $p=2$ functional $L_2 = sum_A sum_i [n_A(i)]^2$ is convex and has
no local minima.  Running it first provides a good starting point for
$p=4$.  The Jacobi rotation angle for a pair $(i,j)$ at $p=2$ is:

$
A_(i j)^((2)) = sum_A
   [-n_i(A)^2 - n_j(A)^2 + 2  n_(i j)(A)^2],
quad quad
B_(i j)^((2)) = 4 sum_A n_(i j)(A) 
   [n_i(A) - n_j(A)],
$
$
phi = frac(1, 2) arctan (
   frac(B_(i j)^((2)), -A_(i j)^((2))) ),
$

where $n_(i j)(A) = sum_(alpha in A) C^("IAO,occ")_(alpha i) C^("IAO,occ")_(alpha j)$.

`calcs.py:175–186`

== 4.3 p=4 refine ($L$ from Eq. 4)

The $p=4$ functional replaces $n_A(i)^2$ with $n_A(i)^4$:

$
A_(i j)^((4)) = sum_A
   [-n_i(A)^4 - n_j(A)^4
         + 6 (n_i(A)^2 + n_j(A)^2)  n_(i j)(A)^2
         + n_i(A)^3 n_j(A) + n_i(A) n_j(A)^3],
$
$
B_(i j)^((4)) = 4 sum_A n_(i j)(A) 
   [n_i(A)^3 - n_j(A)^3],
$
$
phi = frac(1, 4) arctan (
   frac(B_(i j)^((4)), -A_(i j)^((4))) ).
$

The angle is $frac(1, 4)$ (not $frac(1, 2)$) because the leading
power in the trigonometric expansion is $4theta$.

`calcs.py:187–205`

== 4.4 Jacobi sweep algorithm

For each sweep:
1. Loop over all unique orbital pairs $(i, j)$ in *fixed sequential
   order*: for $i = 1 .. n_"occ"$, for $j = 0 .. i-1$.
2. For each pair, compute the $2×2$ rotation angle $phi$ that
   maximises $L$ using the current exponent $p$ (Eq. 5 for $p=2$,
   Eq. 6 for $p=4$).
3. Apply the rotation to columns $i$ and $j$ of
   $bold(C)^("IAO,occ")$.
4. Track the gradient norm $||nabla L||$; stop when it falls below
   $10^(-12)$ (IboView default).

The fixed sequential sweep order is essential for reproducibility.
Because each sweep revisits pairs in the same order, the algorithm
converges to the nearest local maximum of $L$ — which for symmetric
molecules (e.g. benzene) is the symmetric solution.  Randomising the
pair order (or applying an initial perturbation, §4.1) does not
improve the energy degeneracy of symmetry-equivalent bonds; it
merely changes which nearby local maximum is selected, typically
producing a worse (less symmetric) result.

`calcs.py:159–226`

== 4.5 Convergence criterion

After each sweep the gradient norm is:

$
||nabla L|| = frac(1, n_("occ"))sqrt(sum_(i< j) (p  phi_(i j)  B_(i j)^((p)))^2)
$

where $p$ is the current exponent (2 or 4).  Sweeps continue until
$||nabla L|| < 10^(-12)$ or 2048 sweeps are reached.

`calcs.py:218–224`

#line()

= 5 Fock resolution of on-atom degeneracy

The PM functional depends only on atomic populations $n_A(i)$, so any
rotation *within* a subspace of orbitals that all have $n_A(i) approx 1$
on the same atom leaves $L$ unchanged.  This is the *on-atom degeneracy*
problem — for example, the O 1s, O 2s, and O lone-pair orbitals of water
mix arbitrarily after PM localisation.

== 5.1 Detection

For each occupied orbital compute the DOM (*Eq. 2*):

$
delta_i = n_(A)(i)^2 + n_(B)(i)^2,
quad quad
A = "atom with largest " n_A(i),  B = "second largest".
$

If $delta_i > 0.99$, the orbital is essentially mono-atomic.
Orbitals sharing the same dominant atom $A$ and meeting this threshold
form a group.

`calcs.py:256–263`

== 5.2 Fock diagonalisation

For a group of $g$ orbitals ${i_k}_(k=1)^(g)$ on atom $A$, extract
the coefficient block $bold(C)_("block") in RR^(n_"min") × g$ and compute the projected Fock matrix:

$
bold(F)_("block")
= bold(C)_("block")^(T)
  bold(F)^("IAO")
  bold(C)_("block"),
quad quad
bold(F)^("IAO")
= (bold(C)^("IAO"))^(T)
  bold(F)^("AO")
  bold(C)^("IAO").
$

Diagonalise $bold(F)_("block") = bold(V) bold(lambda) bold(V)^(T)$ and rotate the block:

$
bold(C)_("block") arrow
   bold(C)_("block") bold(V).
$

The eigenvalues $bold(lambda)$ become the new orbital energies
(aufbau ordering), and the eigenvectors give the cleanly separated
orbitals (e.g. s-rich lowest, p-rich highest).

`calcs.py:247–254`

#line()

= 6 Valence-virtual IAOs via SVD

The IAO basis ($n_("min")$ functions) accounts for all occupied
orbitals but typically only a fraction of the virtual space.  The
valence-virtual IAOs span the part of the canonical virtual space
that has non-negligible overlap with the IAO basis.

== 6.1 SVD of the cross overlap

$
bold(S)^("IbVir") = (bold(C)^("IAO"))^(T)
   bold(S) bold(C)^("vir")
    in  RR^(n_"min") × n_("vir"),
$

where $bold(C)^("vir")$ are the canonical virtual MO coefficients.
Perform the singular value decomposition:

$
bold(S)^("IbVir") = bold(U) bold(Sigma) bold(V)^(T).
$

`calcs.py:550–552`

== 6.2 Truncation

Keep only singular values $sigma_k > 10^(-8)$:

$
n_("val vir") = |{k : sigma_k > 10^(-8)}|.
$

The valence-virtual IAO coefficients are:

$
bold(U)^("val") = bold(U)_(:,:n_"val vir").
$

`calcs.py:553–554`

== 6.3 Localise the virtual block

The SVD-projected virtuals are delocalised because each right-singular
vector mixes many canonical virtuals.  The same PM $p=2-> p=4$
procedure (Section 4) is applied to $bold(U)^("val")$.

`calcs.py:557–559`

#line()

= 7 Orbital energies

Orbital energies are not canonical eigenvalues after localisation.  They
are computed as the diagonal of the Fock matrix in the basis of the
(now rotated) orbitals:

$
epsilon_i
= bold(c)_i^(T) bold(F)^("AO") bold(c)_i,
quad quad
bold(c)_i = bold(C)^("IAO") bold(C)^("IAO,all")_(:,i),
$

where $bold(c)_i$ are the orbital coefficients in the original AO
basis, and $bold(C)^("IAO,all")$ is the concatenation of occupied
and valence-virtual blocks.  Equivalently, in the IAO basis:

$
epsilon_i
= (bold(C)^("IAO,all")_(:,i))^(T)
  bold(F)^("IAO")
  bold(C)^("IAO,all")_(:,i),
quad quad
bold(F)^("IAO") = (bold(C)^("IAO"))^(T)
                      bold(F)^("AO")
                      bold(C)^("IAO").
$

`calcs.py:534–546, 561–563`

#line()

= 8 Analysis table

The analysis written to `ibos.txt` computes per-orbital descriptors.

== 8.1 DOM (Eq. 2)

$
delta_i = n_A(i)^2 + n_B(i)^2,
quad quad
A = arg max _(A') n_(A')(i),  B = arg max _(A' != A) n_(A')(i).
$

A DOM of 1.00 means the orbital is entirely on one atom (core or lone
pair); 0.50 means a perfectly shared diatomic bond.

`calcs.py:298–301`

== 8.2 s/p/d character on the dominant atom

Fraction of the orbital density on the dominant atom $A$ that comes from
functions of each angular momentum $l$ (0=s, 1=p, 2=d):

$
s = sum_(alpha in A,  e l l_alpha = 0)
    (C^("IAO,all")_(alpha, i))^2,
quad quad
p = sum_(alpha in A,  e l l_alpha = 1)
    (C^("IAO,all")_(alpha, i))^2,
quad quad
d = sum_(alpha in A,  e l l_alpha = 2)
    (C^("IAO,all")_(alpha, i))^2.
$

The hybridisation label is "$s\%$ s + $p\%$ p + $d\%$ d".

`calcs.py:304–306, 385–400`

== 8.3 Bond-type classification

#table(columns: (auto, auto),
  [*Condition*], [*Label*],
  [$n_A(i) > 0.99$ and $s > 0.75$], [`A(Core)`],
  [$n_A(i) > 0.90$], [`A(LP)`],
  [$n_A(i) + n_B(i) > 0.75$, $p_A > 0.85$ and $p_B > 0.85$], [`A-B pi`],
  [$n_A(i) + n_B(i) > 0.75$, otherwise], [`A-B sigma`],
  [Virtual and $n_A(i)+n_B(i) > 0.60$], [`A-B anti*`],
)

`calcs.py:309–351`

#line()

= 9 Why occ-vir delocalization analysis is impossible in the IAO basis

A donor/acceptor delocalization analysis (like NBO second-order perturbation
theory) requires non-zero off-diagonal Fock matrix elements between occupied
and virtual orbitals.  In the IAO basis, these are *exactly zero* by
construction — not a numerical approximation, but a mathematical identity.

== 9.1 Occupied MOs are eigenvectors of $bold(F)^("IAO")$

Let $bold(c)_k$ be the $k$-th occupied canonical MO (column of
$bold(C)^("occ")$).  The Roothaan–Hall equation gives:

$
bold(F)^("AO") bold(c)_k = epsilon_k bold(S)
bold(c)_k. quad quad ("HF")
$

Because IAOs exactly span the occupied space (Section 3),
$bold(C)^("IAO") bold(C)^("IAO,occ") = bold(C)^("occ")$
is an exact identity.  Substitute $bold(c)_k = bold(C)^("IAO") bold(C)^("IAO,occ")_(:,k)$ into (HF):

$
bold(F)^("AO") bold(C)^("IAO")
bold(C)^("IAO,occ")_(:,k)
= epsilon_k bold(S) bold(C)^("IAO")
  bold(C)^("IAO,occ")_(:,k).
$

Left-multiply by $(bold(C)^("IAO"))^(T)$:

$
(bold(C)^("IAO"))^(T) bold(F)^("AO")
bold(C)^("IAO") bold(C)^("IAO,occ")_(:,k)
= epsilon_k (bold(C)^("IAO"))^(T) bold(S)
  bold(C)^("IAO") bold(C)^("IAO,occ")_(:,k).
$

The left factor is $bold(F)^("IAO") bold(C)^("IAO,occ")_(:,k)$.  The right factor uses the IAO orthonormality $(bold(C)^("IAO"))^(T) bold(S) bold(C)^("IAO") = bold(I)$ (Section 3.8), giving:

$
bold(F)^("IAO") bold(C)^("IAO,occ")_(:,k)
= epsilon_k bold(C)^("IAO,occ")_(:,k).
$

Each column of $bold(C)^("IAO,occ")$ is an exact eigenvector of
$bold(F)^("IAO")$ with eigenvalue $epsilon_k$.  This holds for
the *raw* canonical occupied MOs before any PM localisation.

== 9.2 Spectral theorem $=>$ $bold(F)^("IAO")_(o v) = bold(0)$

Since $bold(F)^("IAO")$ is a real symmetric matrix, eigenvectors
corresponding to different eigenvalues are orthogonal.  The occupied
eigenvectors (columns of $bold(C)^("IAO,occ")$) span a subspace
$cal(O) subset RR^(n_"min")$, and the virtual eigenvectors
span the orthogonal complement $cal(V) = cal(O)^(perp)$.

The valence-virtual IAO coefficients $bold(U)^("val")$ (Section 6)
live in $cal(V)$ because the SVD of $bold(S)^("IbVir")$ extracts
the part of the canonical virtual space that projects into the IAO subspace
but is orthogonal to the occupied space.  Therefore:

$
(bold(C)^("IAO,occ"))^(T)
bold(F)^("IAO")
bold(U)^("val")
= bold(0),
quad quad
(bold(U)^("val"))^(T)
bold(F)^("IAO")
bold(C)^("IAO,occ")
= bold(0).
$

The off-diagonal occupied-virtual block of $bold(F)^("IAO")$ is
identically zero.  This is *not* a numerical truncation — it follows
from the spectral theorem applied to the symmetric matrix
$bold(F)^("IAO")$.

== 9.3 Consequences for delocalization analysis

#table(columns: (auto, auto, auto, auto),
  [*Analysis type*], [*Mathematical requirement*], [*IAO-basis value*], [*Why*],
  [Overlap-based donor/acceptor], [$⟨psi_i^("occ")|psi_j^("vir")⟩$], [$0$], [$cal(O) perp cal(V)$ in orthonormal IAO basis],
  [Fock-based (NBO E2)], [$F_(i j)^2 / (epsilon_j - epsilon_i)$], [$0$], [$F_(o v) = 0$ (Section 9.2)],
  [Orbital mixing coefficient], [$F_(i j) / (epsilon_j - epsilon_i)$], [$0$], [Same reason],
)

The only way to obtain non-zero $F_(o v)$ would be to abandon the IAO
projection for virtuals and use the raw canonical virtual MOs.  However,
canonical virtuals contain diffuse and Rydberg character that pollutes
the chemically meaningful valence picture — exactly what the IAO
construction is designed to filter out.

== 9.4 Numerical verification

For water/cc-pVDZ ($n_("min")=7$, $n_("occ")=5$, $n_("val vir")=2$):

- $||bold(F)^("IAO")_(o v)||_F = 5.2 × 10^(-8)$ (machine precision)
- $||(bold(C)^("occ"))^(T) bold(C)^("vir")||_F = 0.0$ (exact, since
  canonical MOs are orthonormal)

Both are zero to within numerical noise, confirming the mathematical
identity.

== 9.5 Why NBO can report non-zero E(2) but IAO cannot

NBO's occupied Lewis-structure orbitals are *not* an exact, lossless
rotation of the canonical occupied MOs.  By design, each NBO carries
a small occupancy leakage into the virtual/Rydberg space — typical
Lewis NBO occupancies are 1.90–1.99, not exactly 2.0.  That leakage
is the hyperconjugative "donation" density, and it produces genuine
off-diagonal Fock elements between the occupied and virtual NBO sets.

IAO/IBO was deliberately designed to do the opposite — give an exact,
basis-independent, *lossless* representation of the occupied space
(Knizia's "unbiased bridge", 2013).  That exactness is a feature for
orbital shapes and isosurface visualisation, but it structurally
forbids any Fock-coupling-based donor-acceptor energy.  Real NBO-style
E(2) numbers require a different, deliberately non-exact localisation
scheme (Psi4 does not ship one natively).

#line()

= References

1. G. Knizia, "Intrinsic Atomic Orbitals: An Unbiased Bridge between
   Quantum Theory and Chemical Concepts", *J. Chem. Theory Comput.*
   *2013*, *9*, 4834–4843.  DOI: #link("https://doi.org/10.1021/ct400687b")[10.1021/ct400687b].
2. G. Knizia, J. E. M. N. Klein, "Electron Flow in Reaction
   Mechanisms—Revealed from First Principles", *Angew. Chem. Int. Ed.*
   *2015*, *54*, 5518–5522.  DOI: #link("https://doi.org/10.1002/anie.201410637")[10.1002/anie.201410637].
3. IboView source code, Rev. A 2021-10-19 (the reference implementation
   for the IAO/2014 algorithm).  Downloaded from
   #link("https://www.iboview.org")[https://www.iboview.org].

