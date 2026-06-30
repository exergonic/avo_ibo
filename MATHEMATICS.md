# Mathematics of the avogadro-ibo Plugin

This document walks through every computational step in
[`calcs.py`](src/avogadro_ibo/calcs.py), connecting each operation to the
equations in Knizia, *J. Chem. Theory Comput.* **2013**, *9*, 4834–4843
([DOI: 10.1021/ct400687b](https://doi.org/10.1021/ct400687b)).

---

## 1 Notation

Let the **AO basis** (the SCF basis, e.g. cc-pVDZ) have $n_{\rm AO}$
functions $\{\lvert\chi_\mu\rangle\}$, and the **minimal basis** (STO-3G,
$n_{\rm min}$) have functions $\{\lvert\tilde\chi_\alpha\rangle\}$.

| Symbol | Shape | Meaning |
|--------|-------|---------|
| $\mathbf{S}$ | $n_{\rm AO} \times n_{\rm AO}$ | AO overlap, $S_{\mu\nu} = \langle\chi_\mu\vert\chi_\nu\rangle$ |
| $\mathbf{S}^{12}$ | $n_{\rm AO} \times n_{\rm min}$ | Cross-overlap, $S^{12}_{\mu\alpha} = \langle\chi_\mu\vert\tilde\chi_\alpha\rangle$ |
| $\mathbf{S}^{\rm min}$ | $n_{\rm min} \times n_{\rm min}$ | Minimal-basis overlap, $S^{\rm min}_{\alpha\beta} = \langle\tilde\chi_\alpha\vert\tilde\chi_\beta\rangle$ |
| $\mathbf{C}^{\rm occ}$ | $n_{\rm AO} \times n_{\rm occ}$ | Occupied MO coefficients (columns = orbitals) |
| $\mathbf{C}^{\rm IAO}$ | $n_{\rm AO} \times n_{\rm min}$ | IAO coefficients in the AO basis |
| $\mathbf{C}^{\rm IAO,occ}$ | $n_{\rm min} \times n_{\rm occ}$ | Occupied MO coefficients in the IAO basis |
| $\mathbf{F}^{\rm AO}$ | $n_{\rm AO} \times n_{\rm AO}$ | Fock matrix in the AO basis |
| $\mathbf{F}^{\rm IAO}$ | $n_{\rm min} \times n_{\rm min}$ | Fock matrix in the IAO basis |
| $n_A(i)$ | $\mathbb{R}$ | Mulliken population of orbital $i$ on atom $A$ |
| $\delta_i$ | $\mathbb{R}$ | **D**iatomic **O**ccupation of orbital $i$ |

Indices run over:
- $\mu,\nu$ — AO functions ($1 \ldots n_{\rm AO}$)
- $\alpha,\beta$ — minimal-basis functions ($1 \ldots n_{\rm min}$)
- $i,j$ — occupied orbitals ($1 \ldots n_{\rm occ}$)
- $A,B$ — atoms ($1 \ldots n_{\rm atom}$)

---

## 2 The SCF wavefunction

Psi4 solves the Hartree–Fock (or DFT) equations for the molecular geometry
from Avogadro:

$$
\mathbf{F}^{\rm AO} \mathbf{C} = \mathbf{S} \mathbf{C} \mathbf{\varepsilon}
$$

The default method is HF with the cc-pVDZ basis and Cartesian
($\mathtt{puream}=0$) functions.  The occupied block
$\mathbf{C}^{\rm occ} \in \mathbb{R}^{n_{\rm AO} \times n_{\rm occ}}$ is
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
\vert\tilde\chi_\alpha^{\rm (AO)}\rangle = \sum_{\mu=1}^{n_{\rm AO}}
(\mathbf{P}^{12})_{\mu\alpha} \,\lvert\chi_\mu\rangle .
$$

`calcs.py:64–65`

### 3.2 Express occupied MOs in the minimal basis

The projection of each occupied orbital onto the minimal-basis functions
gives the minimal-basis representation:

$$
(\mathbf{C}^{\rm occ,min})_{\alpha i}
= \langle\tilde\chi_\alpha\vert\psi_i\rangle
= \sum_{\mu} S^{12}_{\mu\alpha} C^{\rm occ}_{\mu i}
= (\mathbf{S}^{12\,T} \mathbf{C}^{\rm occ})_{\alpha i}.
$$

`calcs.py:68`

### 3.3 Depolarisation — solve $\mathbf{S}^{\rm min} \tilde{\mathbf{C}} = \mathbf{C}^{\rm occ,min}$

We find the coefficients $\tilde{\mathbf{C}}$ that express the occupied
orbitals in the minimal basis:

$$
\mathbf{S}^{\rm min} \tilde{\mathbf{C}} = \mathbf{C}^{\rm occ,min}
\quad\Longrightarrow\quad
\tilde{\mathbf{C}} = (\mathbf{S}^{\rm min})^{-1} \mathbf{C}^{\rm occ,min}.
$$

$(\mathbf{S}^{\rm min})^{-1}$ is obtained via Cholesky factorisation.
`calcs.py:70–72`

### 3.4 Occupied-space metric

The overlap of the projected occupied orbitals within the minimal basis is:

$$
\tilde{\mathbf{S}} = (\mathbf{C}^{\rm occ,min})^{T} \tilde{\mathbf{C}}
= (\mathbf{C}^{\rm occ,min})^{T} (\mathbf{S}^{\rm min})^{-1}
  \mathbf{C}^{\rm occ,min}.
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
\mathbf{T}^{(4)} = \mathbf{C}^{\rm occ}
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
\mathbf{C}^{\rm IAO} = \mathbf{P}^{12}
   + \mathbf{T}^{(4)} (\mathbf{C}^{\rm occ,min})^{T}.
$$

In matrix form: $(\mathbf{C}^{\rm IAO})_{\mu\alpha}$
is the coefficient of AO function $\lvert\chi_\mu\rangle$ in IAO $\lvert\phi_\alpha\rangle$.

`calcs.py:86`

### 3.8 Symmetric (Löwdin) orthogonalisation

The IAOs are not yet orthonormal.  Define the metric of the IAO basis:

$$
\mathbf{M} = (\mathbf{C}^{\rm IAO})^{T} \mathbf{S} \mathbf{C}^{\rm IAO}.
$$

Diagonalise $\mathbf{M} = \mathbf{U} \mathbf{m} \mathbf{U}^{T}$ and apply
$\mathbf{M}^{-1/2} = \mathbf{U} \mathbf{m}^{-1/2} \mathbf{U}^{T}$:

$$
\mathbf{C}^{\rm IAO} \leftarrow
   \mathbf{C}^{\rm IAO} \, \mathbf{M}^{-1/2}.
$$

Now $\langle\phi_\alpha\vert\phi_\beta\rangle
= (\mathbf{C}^{\rm IAO})^{T} \mathbf{S} \mathbf{C}^{\rm IAO} = \mathbf{I}$.

`calcs.py:90–92`

### 3.9 Express occupied orbitals in the IAO basis

The occupied MO coefficients in the orthonormal IAO basis are:

$$
\mathbf{C}^{\rm IAO,occ} = (\mathbf{C}^{\rm IAO})^{T}
                           \mathbf{S} \mathbf{C}^{\rm occ}.
$$

Because IAOs span the occupied space (by construction), this is a
unitary rotation and $\mathbf{C}^{\rm IAO} \mathbf{C}^{\rm IAO,occ}
= \mathbf{C}^{\rm occ}$ is exact.

`calcs.py:97`

---

## 4 Pipek–Mezey localisation in the IAO basis

The IBOs are obtained by maximising the generalised
Pipek–Mezey functional (**Eq. 4** of the paper):

$$
L = \sum_{A=1}^{n_{\rm atom}} \sum_{i=1}^{n_{\rm occ}}
    [n_A(i)]^{4},
\qquad
n_A(i) = \sum_{\alpha\in A} |C^{\rm IAO,occ}_{\alpha i}|^{2}.
$$

Here $n_A(i)$ is the Mulliken population of orbital $i$ on atom $A$,
computed directly from the squared IAO coefficients (the IAO basis is
orthonormal).

### 4.1 Random symmetry breaking

A random orthogonal matrix $\mathbf{U}$ (QR decomposition of a matrix of
standard normal variates) is applied to the occupied block *before*
localisation:

$$
\mathbf{C}^{\rm IAO,occ} \leftarrow
   \mathbf{C}^{\rm IAO,occ} \mathbf{U}.
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
C^{\rm IAO,occ}_{\alpha i} C^{\rm IAO,occ}_{\alpha j}$.

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
   $\mathbf{C}^{\rm IAO,occ}$.
4. Track the gradient norm $\|\nabla L\|$; stop when it falls below
   $10^{-12}$ (IboView default).

`calcs.py:148–207`

### 4.5 Convergence criterion

After each sweep the gradient norm is:

$$
\|\nabla L\| = \frac{1}{n_{\rm occ}}
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
the coefficient block $\mathbf{C}_{\rm block} \in \mathbb{R}^{n_{\rm min}
\times g}$ and compute the projected Fock matrix:

$$
\mathbf{F}_{\rm block}
= \mathbf{C}_{\rm block}^{T}
  \mathbf{F}^{\rm IAO}
  \mathbf{C}_{\rm block},
\qquad
\mathbf{F}^{\rm IAO}
= (\mathbf{C}^{\rm IAO})^{T}
  \mathbf{F}^{\rm AO}
  \mathbf{C}^{\rm IAO}.
$$

Diagonalise $\mathbf{F}_{\rm block} = \mathbf{V} \mathbf{\lambda}
\mathbf{V}^{T}$ and rotate the block:

$$
\mathbf{C}_{\rm block} \leftarrow
   \mathbf{C}_{\rm block} \mathbf{V}.
$$

The eigenvalues $\mathbf{\lambda}$ become the new orbital energies
(aufbau ordering), and the eigenvectors give the cleanly separated
orbitals (e.g. s-rich lowest, p-rich highest).

`calcs.py:247–254`

---

## 6 Valence-virtual IAOs via SVD

The IAO basis ($n_{\rm min}$ functions) accounts for all occupied
orbitals but typically only a fraction of the virtual space.  The
valence-virtual IAOs span the part of the canonical virtual space
that has non-negligible overlap with the IAO basis.

### 6.1 SVD of the cross overlap

$$
\mathbf{S}^{\rm IbVir} = (\mathbf{C}^{\rm IAO})^{T}
   \mathbf{S} \mathbf{C}^{\rm vir}
   \;\in\; \mathbb{R}^{n_{\rm min} \times n_{\rm vir}},
$$

where $\mathbf{C}^{\rm vir}$ are the canonical virtual MO coefficients.
Perform the singular value decomposition:

$$
\mathbf{S}^{\rm IbVir} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^{T}.
$$

`calcs.py:550–552`

### 6.2 Truncation

Keep only singular values $\sigma_k > 10^{-8}$:

$$
n_{\rm val\,vir} = \bigl|\{k : \sigma_k > 10^{-8}\}\bigr|.
$$

The valence-virtual IAO coefficients are:

$$
\mathbf{U}^{\rm val} = \mathbf{U}_{:,:n_{\rm val\,vir}}.
$$

`calcs.py:553–554`

### 6.3 Localise the virtual block

The SVD-projected virtuals are delocalised because each right-singular
vector mixes many canonical virtuals.  The same PM $p=2\to p=4$
procedure (Section 4) is applied to $\mathbf{U}^{\rm val}$.

`calcs.py:557–559`

---

## 7 Orbital energies

Orbital energies are not canonical eigenvalues after localisation.  They
are computed as the diagonal of the Fock matrix in the basis of the
(now rotated) orbitals:

$$
\varepsilon_i
= \mathbf{c}_i^{T} \mathbf{F}^{\rm AO} \mathbf{c}_i,
\qquad
\mathbf{c}_i = \mathbf{C}^{\rm IAO} \mathbf{C}^{\rm IAO,all}_{:,i},
$$

where $\mathbf{c}_i$ are the orbital coefficients in the original AO
basis, and $\mathbf{C}^{\rm IAO,all}$ is the concatenation of occupied
and valence-virtual blocks.  Equivalently, in the IAO basis:

$$
\varepsilon_i
= (\mathbf{C}^{\rm IAO,all}_{:,i})^{T}
  \mathbf{F}^{\rm IAO}
  \mathbf{C}^{\rm IAO,all}_{:,i},
\qquad
\mathbf{F}^{\rm IAO} = (\mathbf{C}^{\rm IAO})^{T}
                      \mathbf{F}^{\rm AO}
                      \mathbf{C}^{\rm IAO}.
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
    \bigl(C^{\rm IAO,all}_{\alpha, i}\bigr)^2,
\qquad
p = \sum_{\alpha \in A,\; \ell_\alpha = 1}
    \bigl(C^{\rm IAO,all}_{\alpha, i}\bigr)^2,
\qquad
d = \sum_{\alpha \in A,\; \ell_\alpha = 2}
    \bigl(C^{\rm IAO,all}_{\alpha, i}\bigr)^2.
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
