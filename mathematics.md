# Mathematics of the avogadro-ibo Plugin

This document walks through every computational step in
[`calcs.py`](src/avogadro_ibo/calcs.py), connecting each operation to the
equations in Knizia, *J. Chem. Theory Comput.* **2013**, *9*, 4834–4843
([DOI: 10.1021/ct400687b](https://doi.org/10.1021/ct400687b)).

---

## 1 Notation

Let the **AO basis** (the SCF basis, e.g. cc-pVDZ) have $n_{\mathrm{AO}}$
functions $\{\lvert\chi_\mu\rangle\}$, and the **minimal basis** (STO-3G,
$n_{\mathrm{min}}$) have functions $\{\lvert\tilde\chi_\alpha\rangle\}$.

| Symbol | Shape | Meaning |
|--------|-------|---------|
| $\mathbf{S}$ | $n_{\mathrm{AO}} \times n_{\mathrm{AO}}$ | AO overlap, $S_{\mu\nu} = \langle\chi_\mu\vert\chi_\nu\rangle$ |
| $\mathbf{S}^{12}$ | $n_{\mathrm{AO}} \times n_{\mathrm{min}}$ | Cross-overlap, $S^{12}_{\mu\alpha} = \langle\chi_\mu\vert\tilde\chi_\alpha\rangle$ |
| $\mathbf{S}^{\mathrm{min}}$ | $n_{\mathrm{min}} \times n_{\mathrm{min}}$ | Minimal-basis overlap, $S^{\mathrm{min}}_{\alpha\beta} = \langle\tilde\chi_\alpha\vert\tilde\chi_\beta\rangle$ |
| $\mathbf{C}^{\mathrm{occ}}$ | $n_{\mathrm{AO}} \times n_{\mathrm{occ}}$ | Occupied MO coefficients (columns = orbitals) |
| $\mathbf{C}^{\mathrm{IAO}}$ | $n_{\mathrm{AO}} \times n_{\mathrm{min}}$ | IAO coefficients in the AO basis |
| $\mathbf{C}^{\mathrm{IAO},occ}$ | $n_{\mathrm{min}} \times n_{\mathrm{occ}}$ | Occupied MO coefficients in the IAO basis |
| $\mathbf{F}^{\mathrm{AO}}$ | $n_{\mathrm{AO}} \times n_{\mathrm{AO}}$ | Fock matrix in the AO basis |
| $\mathbf{F}^{\mathrm{IAO}}$ | $n_{\mathrm{min}} \times n_{\mathrm{min}}$ | Fock matrix in the IAO basis |
| $n_A(i)$ | $\mathbb{R}$ | Mulliken population of orbital $i$ on atom $A$ |
| $\delta_i$ | $\mathbb{R}$ | **D**iatomic **O**ccupation of orbital $i$ |

Indices run over:
- $\mu,\nu$ — AO functions ($1 \ldots n_{\mathrm{AO}}$)
- $\alpha,\beta$ — minimal-basis functions ($1 \ldots n_{\mathrm{min}}$)
- $i,j$ — occupied orbitals ($1 \ldots n_{\mathrm{occ}}$)
- $A,B$ — atoms ($1 \ldots n_{\mathrm{atom}}$)

---

## 2 The SCF wavefunction

Psi4 solves the Hartree–Fock (or DFT) equations for the molecular geometry
from Avogadro:

$$
\mathbf{F}^{\mathrm{AO}} \mathbf{C} = \mathbf{S} \mathbf{C} \mathbf{\varepsilon}
$$

The default method is HF with the cc-pVDZ basis and Cartesian
($\mathtt{puream}=0$) functions.  The occupied block
$\mathbf{C}^{\mathrm{occ}} \in \mathbb{R}^{n_{\mathrm{AO}} \times n_{\mathrm{occ}}}$ is
extracted from the converged solution (`calcs.py:509`).

---

## 3 Intrinsic Atomic Orbitals (IAOs)

IAOs form an orthonormal, chemically balanced basis that exactly spans the
occupied space while keeping the atomic character of the minimal basis.
The construction follows **Appendix C** of the paper (IboView's
`MakeIaoBasisNew`).

### 3.1 Project the minimal basis into the full AO space

Each minimal-basis function $\lvert\tilde\chi_\alpha\rangle$ is expanded in
the AO basis via the overlap:

$$
\mathbf{P}^{12} = \mathbf{S}^{-1} \mathbf{S}^{12},
\qquad
\vert\tilde\chi_\alpha^{\mathrm{(AO)}}\rangle = \sum_{\mu=1}^{n_{\mathrm{AO}}}
(\mathbf{P}^{12})_{\mu\alpha} \,\lvert\chi_\mu\rangle .
$$

`calcs.py:64–65`

### 3.2 Express occupied MOs in the minimal basis

The projection of each occupied orbital onto the minimal-basis functions
gives the minimal-basis representation:

$$
(\mathbf{C}^{\mathrm{occ},min})_{\alpha i}
= \langle\tilde\chi_\alpha\vert\psi_i\rangle
= \sum_{\mu} S^{12}_{\mu\alpha} C^{\mathrm{occ}}_{\mu i}
= (\mathbf{S}^{12\,T} \mathbf{C}^{\mathrm{occ}})_{\alpha i}.
$$

`calcs.py:68`

### 3.3 Depolarisation — solve $\mathbf{S}^{\mathrm{min}} \tilde{\mathbf{C}} = \mathbf{C}^{\mathrm{occ},min}$

We find the coefficients $\tilde{\mathbf{C}}$ that express the occupied
orbitals in the minimal basis:

$$
\mathbf{S}^{\mathrm{min}} \tilde{\mathbf{C}} = \mathbf{C}^{\mathrm{occ},min}
\quad\Longrightarrow\quad
\tilde{\mathbf{C}} = (\mathbf{S}^{\mathrm{min}})^{-1} \mathbf{C}^{\mathrm{occ},min}.
$$

$(\mathbf{S}^{\mathrm{min}})^{-1}$ is obtained via Cholesky factorisation.
`calcs.py:70–72`

### 3.4 Occupied-space metric

The overlap of the projected occupied orbitals within the minimal basis is:

$$
\tilde{\mathbf{S}} = (\mathbf{C}^{\mathrm{occ},min})^{T} \tilde{\mathbf{C}}
= (\mathbf{C}^{\mathrm{occ},min})^{T} (\mathbf{S}^{\mathrm{min}})^{-1}
  \mathbf{C}^{\mathrm{occ},min}.
$$

`calcs.py:75`

### 3.5 Inverse metric — $\tilde{\mathbf{C}}^{(2)}$

We solve for the coefficients $\tilde{\mathbf{C}}^{(2)}$:

$$
\tilde{\mathbf{S}} \tilde{\mathbf{C}}^{(2)\,T} = \tilde{\mathbf{C}}^{T}
\quad\Longrightarrow\quad
\tilde{\mathbf{C}}^{(2)} = \tilde{\mathbf{C}} \tilde{\mathbf{S}}^{-1}.
$$

`calcs.py:78–80`

### 3.6 Tight residual — $\mathbf{T}^{(4)}$

The component of each occupied orbital that lies *orthogonal* to the
minimal-basis span is the tight residual:

$$
\mathbf{T}^{(4)} = \mathbf{C}^{\mathrm{occ}}
   - \mathbf{P}^{12} \tilde{\mathbf{C}}^{(2)}.
$$

`calcs.py:83`

The columns of $\mathbf{T}^{(4)}$ are the occupied-space vectors that
the minimal basis alone cannot describe.

### 3.7 Construct the IAO coefficients

Each IAO is the corresponding projected minimal-basis function plus the
tight residual weighted by the overlap of that minimal-basis function with
the occupied orbitals:

$$
\mathbf{C}^{\mathrm{IAO}} = \mathbf{P}^{12}
   + \mathbf{T}^{(4)} (\mathbf{C}^{\mathrm{occ},min})^{T}.
$$

In matrix form: $(\mathbf{C}^{\mathrm{IAO}})_{\mu\alpha}$
is the coefficient of AO function $\lvert\chi_\mu\rangle$ in IAO $\lvert\phi_\alpha\rangle$.

`calcs.py:86`

### 3.8 Symmetric (Löwdin) orthogonalisation

The IAOs are not yet orthonormal.  Define the metric of the IAO basis:

$$
\mathbf{M} = (\mathbf{C}^{\mathrm{IAO}})^{T} \mathbf{S} \mathbf{C}^{\mathrm{IAO}}.
$$

Diagonalise $\mathbf{M} = \mathbf{U} \mathbf{m} \mathbf{U}^{T}$ and apply
$\mathbf{M}^{-1/2} = \mathbf{U} \mathbf{m}^{-1/2} \mathbf{U}^{T}$:

$$
\mathbf{C}^{\mathrm{IAO}} \leftarrow
   \mathbf{C}^{\mathrm{IAO}} \, \mathbf{M}^{-1/2}.
$$

Now $\langle\phi_\alpha\vert\phi_\beta\rangle
= (\mathbf{C}^{\mathrm{IAO}})^{T} \mathbf{S} \mathbf{C}^{\mathrm{IAO}} = \mathbf{I}$.

`calcs.py:90–92`

### 3.9 Express occupied orbitals in the IAO basis

The occupied MO coefficients in the orthonormal IAO basis are:

$$
\mathbf{C}^{\mathrm{IAO},occ} = (\mathbf{C}^{\mathrm{IAO}})^{T}
                           \mathbf{S} \mathbf{C}^{\mathrm{occ}}.
$$

Because IAOs span the occupied space (by construction), this is a
unitary rotation and $\mathbf{C}^{\mathrm{IAO}} \mathbf{C}^{\mathrm{IAO},occ}
= \mathbf{C}^{\mathrm{occ}}$ is exact.

`calcs.py:97`

---

## 4 Pipek–Mezey localisation in the IAO basis

The IBOs are obtained by maximising the generalised
Pipek–Mezey functional (**Eq. 4** of the paper):

$$
L = \sum_{A=1}^{n_{\mathrm{atom}}} \sum_{i=1}^{n_{\mathrm{occ}}}
    [n_A(i)]^{4},
\qquad
n_A(i) = \sum_{\alpha\in A} |C^{\mathrm{IAO},occ}_{\alpha i}|^{2}.
$$

Here $n_A(i)$ is the Mulliken population of orbital $i$ on atom $A$,
computed directly from the squared IAO coefficients (the IAO basis is
orthonormal).

### 4.1 Random symmetry breaking

A random orthogonal matrix $\mathbf{U}$ (QR decomposition of a matrix of
standard normal variates) is applied to the occupied block *before*
localisation:

$$
\mathbf{C}^{\mathrm{IAO},occ} \leftarrow
   \mathbf{C}^{\mathrm{IAO},occ} \mathbf{U}.
$$

This breaks any accidental near-degeneracies that could trap the
Jacobi sweep in a shallow local minimum.

`calcs.py:140–142`

### 4.2 p=2 warm-start

The $p=2$ functional $L_2 = \sum_A \sum_i [n_A(i)]^2$ is convex and has
no local minima.  Running it first provides a good starting point for
$p=4$.  The Jacobi rotation angle for a pair $(i,j)$ at $p=2$ is:

$$
A_{ij}^{(2)} = \sum_A
   \bigl[-n_i(A)^2 - n_j(A)^2 + 2\, n_{ij}(A)^2\bigr],
\qquad
B_{ij}^{(2)} = 4 \sum_A n_{ij}(A)\,
   \bigl[n_i(A) - n_j(A)\bigr],
$$
$$
\phi = \frac{1}{2} \arctan\!\left(
   \frac{B_{ij}^{(2)}}{-A_{ij}^{(2)}} \right),
$$

where $n_{ij}(A) = \sum_{\alpha\in A}
C^{\mathrm{IAO},occ}_{\alpha i} C^{\mathrm{IAO},occ}_{\alpha j}$.

`calcs.py:163–174`

### 4.3 p=4 refine ($L$ from Eq. 4)

The $p=4$ functional replaces $n_A(i)^2$ with $n_A(i)^4$:

$$
A_{ij}^{(4)} = \sum_A
   \bigl[-n_i(A)^4 - n_j(A)^4
         + 6\,(n_i(A)^2 + n_j(A)^2)\, n_{ij}(A)^2
         + n_i(A)^3 n_j(A) + n_i(A) n_j(A)^3\bigr],
$$
$$
B_{ij}^{(4)} = 4 \sum_A n_{ij}(A)\,
   \bigl[n_i(A)^3 - n_j(A)^3\bigr],
$$
$$
\phi = \frac{1}{4} \arctan\!\left(
   \frac{B_{ij}^{(4)}}{-A_{ij}^{(4)}} \right).
$$

The angle is $\frac{1}{4}$ (not $\frac{1}{2}$) because the leading
power in the trigonometric expansion is $4\theta$.

`calcs.py:176–189`

### 4.4 Jacobi sweep algorithm

For each sweep:
1. Loop over all unique orbital pairs $(i, j)$.
2. For each pair, compute the $2\times2$ rotation angle $\phi$ that
   maximises $L$.
3. Apply the rotation to columns $i$ and $j$ of
   $\mathbf{C}^{\mathrm{IAO},occ}$.
4. Track the gradient norm $\|\nabla L\|$; stop when it falls below
   $10^{-12}$ (IboView default).

`calcs.py:148–207`

### 4.5 Convergence criterion

After each sweep the gradient norm is:

$$
\|\nabla L\| = \frac{1}{n_{\mathrm{occ}}}
   \sqrt{\sum_{i<j} \bigl(p \,\phi_{ij}\, B_{ij}^{(p)}\bigr)^2}
$$

where $p$ is the current exponent (2 or 4).  Sweeps continue until
$\|\nabla L\| < 10^{-12}$ or 2048 sweeps are reached.

`calcs.py:202–206`

---

## 5 Fock resolution of on-atom degeneracy

The PM functional depends only on atomic populations $n_A(i)$, so any
rotation **within** a subspace of orbitals that all have $n_A(i) \approx 1$
on the same atom leaves $L$ unchanged.  This is the **on-atom degeneracy**
problem — for example, the O 1s, O 2s, and O lone-pair orbitals of water
mix arbitrarily after PM localisation.

### 5.1 Detection

For each occupied orbital compute the DOM (**Eq. 2**):

$$
\delta_i = n_{A}(i)^2 + n_{B}(i)^2,
\qquad
A = \text{atom with largest } n_A(i),\; B = \text{second largest}.
$$

If $\delta_i > 0.99$, the orbital is essentially mono-atomic.
Orbitals sharing the same dominant atom $A$ and meeting this threshold
form a group.

`calcs.py:239–245`

### 5.2 Fock diagonalisation

For a group of $g$ orbitals $\{i_k\}_{k=1}^{g}$ on atom $A$, extract
the coefficient block $\mathbf{C}_{\mathrm{block}} \in \mathbb{R}^{n_{\mathrm{min}}
\times g}$ and compute the projected Fock matrix:

$$
\mathbf{F}_{\mathrm{block}}
= \mathbf{C}_{\mathrm{block}}^{T}
  \mathbf{F}^{\mathrm{IAO}}
  \mathbf{C}_{\mathrm{block}},
\qquad
\mathbf{F}^{\mathrm{IAO}}
= (\mathbf{C}^{\mathrm{IAO}})^{T}
  \mathbf{F}^{\mathrm{AO}}
  \mathbf{C}^{\mathrm{IAO}}.
$$

Diagonalise $\mathbf{F}_{\mathrm{block}} = \mathbf{V} \mathbf{\lambda}
\mathbf{V}^{T}$ and rotate the block:

$$
\mathbf{C}_{\mathrm{block}} \leftarrow
   \mathbf{C}_{\mathrm{block}} \mathbf{V}.
$$

The eigenvalues $\mathbf{\lambda}$ become the new orbital energies
(aufbau ordering), and the eigenvectors give the cleanly separated
orbitals (e.g. s-rich lowest, p-rich highest).

`calcs.py:247–254`

---

## 6 Valence-virtual IAOs via SVD

The IAO basis ($n_{\mathrm{min}}$ functions) accounts for all occupied
orbitals but typically only a fraction of the virtual space.  The
valence-virtual IAOs span the part of the canonical virtual space
that has non-negligible overlap with the IAO basis.

### 6.1 SVD of the cross overlap

$$
\mathbf{S}^{\mathrm{IbVir}} = (\mathbf{C}^{\mathrm{IAO}})^{T}
   \mathbf{S} \mathbf{C}^{\mathrm{vir}}
   \;\in\; \mathbb{R}^{n_{\mathrm{min}} \times n_{\mathrm{vir}}},
$$

where $\mathbf{C}^{\mathrm{vir}}$ are the canonical virtual MO coefficients.
Perform the singular value decomposition:

$$
\mathbf{S}^{\mathrm{IbVir}} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^{T}.
$$

`calcs.py:550–552`

### 6.2 Truncation

Keep only singular values $\sigma_k > 10^{-8}$:

$$
n_{\mathrm{val\,vir}} = \bigl|\{k : \sigma_k > 10^{-8}\}\bigr|.
$$

The valence-virtual IAO coefficients are:

$$
\mathbf{U}^{\mathrm{val}} = \mathbf{U}_{:,:n_{\mathrm{val\,vir}}}.
$$

`calcs.py:553–554`

### 6.3 Localise the virtual block

The SVD-projected virtuals are delocalised because each right-singular
vector mixes many canonical virtuals.  The same PM $p=2\to p=4$
procedure (Section 4) is applied to $\mathbf{U}^{\mathrm{val}}$.

`calcs.py:557–559`

---

## 7 Orbital energies

Orbital energies are not canonical eigenvalues after localisation.  They
are computed as the diagonal of the Fock matrix in the basis of the
(now rotated) orbitals:

$$
\varepsilon_i
= \mathbf{c}_i^{T} \mathbf{F}^{\mathrm{AO}} \mathbf{c}_i,
\qquad
\mathbf{c}_i = \mathbf{C}^{\mathrm{IAO}} \mathbf{C}^{\mathrm{IAO},all}_{:,i},
$$

where $\mathbf{c}_i$ are the orbital coefficients in the original AO
basis, and $\mathbf{C}^{\mathrm{IAO},all}$ is the concatenation of occupied
and valence-virtual blocks.  Equivalently, in the IAO basis:

$$
\varepsilon_i
= (\mathbf{C}^{\mathrm{IAO},all}_{:,i})^{T}
  \mathbf{F}^{\mathrm{IAO}}
  \mathbf{C}^{\mathrm{IAO},all}_{:,i},
\qquad
\mathbf{F}^{\mathrm{IAO}} = (\mathbf{C}^{\mathrm{IAO}})^{T}
                      \mathbf{F}^{\mathrm{AO}}
                      \mathbf{C}^{\mathrm{IAO}}.
$$

`calcs.py:534–546, 561–563`

---

## 8 Analysis table

The analysis written to `ibos.txt` computes per-orbital descriptors.

### 8.1 DOM (Eq. 2)

$$
\delta_i = n_A(i)^2 + n_B(i)^2,
\qquad
A = \arg\max_{A'} n_{A'}(i),\; B = \arg\max_{A' \neq A} n_{A'}(i).
$$

A DOM of 1.00 means the orbital is entirely on one atom (core or lone
pair); 0.50 means a perfectly shared diatomic bond.

`calcs.py:298–301`

### 8.2 s/p/d character on the dominant atom

Fraction of the orbital density on the dominant atom $A$ that comes from
functions of each angular momentum $l$ (0=s, 1=p, 2=d):

$$
s = \sum_{\alpha \in A,\; \ell_\alpha = 0}
    \bigl(C^{\mathrm{IAO},all}_{\alpha, i}\bigr)^2,
\qquad
p = \sum_{\alpha \in A,\; \ell_\alpha = 1}
    \bigl(C^{\mathrm{IAO},all}_{\alpha, i}\bigr)^2,
\qquad
d = \sum_{\alpha \in A,\; \ell_\alpha = 2}
    \bigl(C^{\mathrm{IAO},all}_{\alpha, i}\bigr)^2.
$$

The hybridisation label is "${s}\%$ s + ${p}\%$ p + ${d}\%$ d".

`calcs.py:304–306, 385–400`

### 8.3 Bond-type classification

| Condition | Label |
|-----------|-------|
| $n_A(i) > 0.99$ and $s > 0.75$ | `A(Core)` |
| $n_A(i) > 0.90$ | `A(LP)` |
| $n_A(i) + n_B(i) > 0.75$, $p_A > 0.85$ and $p_B > 0.85$ | `A-B pi` |
| $n_A(i) + n_B(i) > 0.75$, otherwise | `A-B sigma` |
| Virtual and $n_A(i)+n_B(i) > 0.60$ | `A-B anti*` |

`calcs.py:309–351`

---

## 9 Why occ-vir delocalization analysis is impossible in the IAO basis

A donor/acceptor delocalization analysis (like NBO second-order perturbation
theory) requires non-zero off-diagonal Fock matrix elements between occupied
and virtual orbitals.  In the IAO basis, these are **exactly zero** by
construction — not a numerical approximation, but a mathematical identity.

### 9.1 Occupied MOs are eigenvectors of $\mathbf{F}^{\mathrm{IAO}}$

Let $\mathbf{c}_k$ be the $k$-th occupied canonical MO (column of
$\mathbf{C}^{\mathrm{occ}}$).  The Roothaan–Hall equation gives:

$$
\mathbf{F}^{\mathrm{AO}} \mathbf{c}_k = \varepsilon_k \mathbf{S}
\mathbf{c}_k. \tag{HF}
$$

Because IAOs exactly span the occupied space (Section 3),
$\mathbf{C}^{\mathrm{IAO}} \mathbf{C}^{\mathrm{IAO},occ} = \mathbf{C}^{\mathrm{occ}}$
is an exact identity.  Substitute $\mathbf{c}_k = \mathbf{C}^{\mathrm{IAO}}
\mathbf{C}^{\mathrm{IAO},occ}_{:,k}$ into (HF):

$$
\mathbf{F}^{\mathrm{AO}} \mathbf{C}^{\mathrm{IAO}}
\mathbf{C}^{\mathrm{IAO},occ}_{:,k}
= \varepsilon_k \mathbf{S} \mathbf{C}^{\mathrm{IAO}}
  \mathbf{C}^{\mathrm{IAO},occ}_{:,k}.
$$

Left-multiply by $(\mathbf{C}^{\mathrm{IAO}})^{T}$:

$$
(\mathbf{C}^{\mathrm{IAO}})^{T} \mathbf{F}^{\mathrm{AO}}
\mathbf{C}^{\mathrm{IAO}} \mathbf{C}^{\mathrm{IAO},occ}_{:,k}
= \varepsilon_k (\mathbf{C}^{\mathrm{IAO}})^{T} \mathbf{S}
  \mathbf{C}^{\mathrm{IAO}} \mathbf{C}^{\mathrm{IAO},occ}_{:,k}.
$$

The left factor is $\mathbf{F}^{\mathrm{IAO}}
\mathbf{C}^{\mathrm{IAO},occ}_{:,k}$.  The right factor uses the IAO
orthonormality $(\mathbf{C}^{\mathrm{IAO}})^{T} \mathbf{S}
\mathbf{C}^{\mathrm{IAO}} = \mathbf{I}$ (Section 3.8), giving:

$$
\mathbf{F}^{\mathrm{IAO}} \mathbf{C}^{\mathrm{IAO},occ}_{:,k}
= \varepsilon_k \mathbf{C}^{\mathrm{IAO},occ}_{:,k}.
$$

Each column of $\mathbf{C}^{\mathrm{IAO},occ}$ is an exact eigenvector of
$\mathbf{F}^{\mathrm{IAO}}$ with eigenvalue $\varepsilon_k$.  This holds for
the **raw** canonical occupied MOs before any PM localisation.

### 9.2 Spectral theorem $\Rightarrow$ $\mathbf{F}^{\mathrm{IAO}}_{ov} = \mathbf{0}$

Since $\mathbf{F}^{\mathrm{IAO}}$ is a real symmetric matrix, eigenvectors
corresponding to different eigenvalues are orthogonal.  The occupied
eigenvectors (columns of $\mathbf{C}^{\mathrm{IAO},occ}$) span a subspace
$\mathcal{O} \subset \mathbb{R}^{n_{\mathrm{min}}}$, and the virtual eigenvectors
span the orthogonal complement $\mathcal{V} = \mathcal{O}^{\perp}$.

The valence-virtual IAO coefficients $\mathbf{U}^{\mathrm{val}}$ (Section 6)
live in $\mathcal{V}$ because the SVD of $\mathbf{S}^{\mathrm{IbVir}}$ extracts
the part of the canonical virtual space that projects into the IAO subspace
but is orthogonal to the occupied space.  Therefore:

$$
(\mathbf{C}^{\mathrm{IAO},occ})^{T}
\mathbf{F}^{\mathrm{IAO}}
\mathbf{U}^{\mathrm{val}}
= \mathbf{0},
\qquad
(\mathbf{U}^{\mathrm{val}})^{T}
\mathbf{F}^{\mathrm{IAO}}
\mathbf{C}^{\mathrm{IAO},occ}
= \mathbf{0}.
$$

The off-diagonal occupied-virtual block of $\mathbf{F}^{\mathrm{IAO}}$ is
identically zero.  This is **not** a numerical truncation — it follows
from the spectral theorem applied to the symmetric matrix
$\mathbf{F}^{\mathrm{IAO}}$.

### 9.3 Consequences for delocalization analysis

| Analysis type | Mathematical requirement | IAO-basis value | Why |
|---|---|---|---|
| Overlap-based donor/acceptor | $\langle\psi_i^{\mathrm{occ}}\vert\psi_j^{\mathrm{vir}}\rangle$ | $0$ | $\mathcal{O} \perp \mathcal{V}$ in orthonormal IAO basis |
| Fock-based (NBO E2) | $F_{ij}^2 / (\varepsilon_j - \varepsilon_i)$ | $0$ | $F_{ov} = 0$ (Section 9.2) |
| Orbital mixing coefficient | $F_{ij} / (\varepsilon_j - \varepsilon_i)$ | $0$ | Same reason |

The only way to obtain non-zero $F_{ov}$ would be to abandon the IAO
projection for virtuals and use the raw canonical virtual MOs.  However,
canonical virtuals contain diffuse and Rydberg character that pollutes
the chemically meaningful valence picture — exactly what the IAO
construction is designed to filter out.

### 9.4 Numerical verification

For water/cc-pVDZ ($n_{\mathrm{min}}=7$, $n_{\mathrm{occ}}=5$, $n_{\mathrm{val\,vir}}=2$):

- $\|\mathbf{F}^{\mathrm{IAO}}_{ov}\|_F = 5.2 \times 10^{-8}$ (machine precision)
- $\|(\mathbf{C}^{\mathrm{occ}})^{T} \mathbf{C}^{\mathrm{vir}}\|_F = 0.0$ (exact, since
  canonical MOs are orthonormal)

Both are zero to within numerical noise, confirming the mathematical
identity.

### 9.5 Why NBO can report non-zero E(2) but IAO cannot

NBO's occupied Lewis-structure orbitals are **not** an exact, lossless
rotation of the canonical occupied MOs.  By design, each NBO carries
a small occupancy leakage into the virtual/Rydberg space — typical
Lewis NBO occupancies are 1.90–1.99, not exactly 2.0.  That leakage
is the hyperconjugative "donation" density, and it produces genuine
off-diagonal Fock elements between the occupied and virtual NBO sets.

IAO/IBO was deliberately designed to do the opposite — give an exact,
basis-independent, **lossless** representation of the occupied space
(Knizia's "unbiased bridge", 2013).  That exactness is a feature for
orbital shapes and isosurface visualisation, but it structurally
forbids any Fock-coupling-based donor-acceptor energy.  Real NBO-style
E(2) numbers require a different, deliberately non-exact localisation
scheme (Psi4 does not ship one natively).

---

## References

1. G. Knizia, "Intrinsic Atomic Orbitals: An Unbiased Bridge between
   Quantum Theory and Chemical Concepts", *J. Chem. Theory Comput.*
   **2013**, *9*, 4834–4843.  DOI: [10.1021/ct400687b](https://doi.org/10.1021/ct400687b).
2. G. Knizia, J. E. M. N. Klein, "Electron Flow in Reaction
   Mechanisms—Revealed from First Principles", *Angew. Chem. Int. Ed.*
   **2015**, *54*, 5518–5522.  DOI: [10.1002/anie.201410637](https://doi.org/10.1002/anie.201410637).
3. IboView source code, Rev. A 2021-10-19 (the reference implementation
   for the IAO/2014 algorithm).  Downloaded from
   [https://www.iboview.org](https://www.iboview.org).
