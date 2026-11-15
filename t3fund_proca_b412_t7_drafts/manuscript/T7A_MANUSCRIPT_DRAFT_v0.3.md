# Tree-Level Proca Integration and Independent IBP Classification of Higher-Derivative Operators in a Two-Scalar EFT

**Dimitar Kretski**
Center for Hydro- and Aerodynamics, Bulgarian Academy of Sciences, Varna, Bulgaria · ORCID 0000-0001-5108-2243

*Manuscript draft v0.3 (2026-11-15). Technical Markdown version for internal checking against the repository; not yet formatted for Zenodo or a journal.*

*Change from v0.2: the editorial note `AUTHOR TO CONFIRM` in §1.1 is resolved by the author's decision of 2026-11-15; the paragraph now states the scope of the README statement explicitly. No result, formula or limitation is changed.*

*Changes from v0.1 (including two editorial corrections after the review of 2026-11-15, 10:48, with the version number kept): verified bibliography and project-document references (§1, References); terminology of the cross-check clarified, and the basis of the equality $\ker E|_{V_S} = W_S$ stated explicitly (§6); explicit statement that only the rank-one subfamily of B4 is reached (§3.5, §7.2); relation to the project's existing Proca verdicts (§1.1); reproducibility table re-verified against the repository blobs (§8).*

> **Status of the results.** Everything reported here is a **verified draft result** of the T7a work package of the project `Kretski/B4-BEC-Lambda-Bridge`, taken exclusively from the commits 61fefff, a72c323, e7ab8d4, 98780c9, 4074199 and 2fc4d16. T7a is **not registered**. The results are not a registered outcome, they do not evaluate any acceptance criterion of a registered protocol, and they are not evidence that the EFT or its UV completion is consistent (perturbative, unitary, or valid in any specified regime). Those questions belong to the separate work package B4.12, which is open and is not used anywhere in this manuscript.

---

## Abstract

We study a 1+1-dimensional two-scalar effective field theory (EFT) whose four-derivative contact sector, $\tfrac12[p(f_tf_x)^2 - 2\kappa(f_tf_x)(g_tg_x) + s(g_tg_x)^2]$, is postulated in the B4 model. We ask whether this operator structure is reproduced exactly when a massive Proca field, coupled to the light scalars through $A_1 J$ with $J = g_\phi f_t f_x + g_\psi g_t g_x$, is integrated out at tree level. Starting from the full action, with the Gauss constraint retained, we derive the exact tree-level kernel $K(\omega,k) = (1+k^2/M^2)/(M^2+k^2-\omega^2)$. The leading order reproduces the B4 structure with $p = g_\phi^2/M^2$, $s = g_\psi^2/M^2$ and $\kappa = -g_\phi g_\psi/M^2$, so the rank-one relation $\kappa^2 = ps$ is automatic: a single Proca mediator reaches only this subfamily of B4. The two relative-sign regimes are distinct. The first correction, $\Delta\mathcal L^{(M^{-4})} = (\partial_t J)^2/(2M^4)$, lies outside the span of the B4 operators by a derivative-counting argument. We classify it in an off-shell basis defined by integration by parts (IBP) alone, with no light-field equations of motion and no field redefinitions. The quotient spaces have dimensions 3, 12 and 3 in the $f^4$, $f^2g^2$ and $g^4$ sectors. Two algorithmically independent implementations of the same equivalence obtain exactly the same canonical bases and rational coordinates: one uses the Euler operator (SymPy), the other an explicit IBP relation matrix in exact rational arithmetic (pure Python). All scripts, outputs and hashes are public.

---

## 1. Introduction

The B4 model is a two-scalar EFT in 1+1 dimensions. Its departure from free propagation is carried by four-derivative contact operators with couplings $p$, $s$, $\kappa$, together with a quartic term $\Lambda f_x^4$. A natural question is whether such a contact structure can be the low-energy remnant of a simple heavy mediator. The candidate here is a massive vector (Proca) field $A_\mu$ of mass $M$, coupled linearly to a bilinear current of the light fields.

Two questions must be kept apart:

- **(T7a)** Does tree-level integration of the Proca field reproduce the B4 operator structure, including normalisations, signs and the rank-one constraint? Which operators appear at the next order in $1/M^2$, and how do they classify in a basis fixed without reference to the light-field dynamics?
- **(T7b)** Can the $\Lambda$-effect emerge as a correction instead of being put in by hand? In the Proca action used here, $(\Lambda/4) f_x^4$ is introduced explicitly. Tree-level integration therefore cannot speak to T7b. **This manuscript addresses T7a only.**

[AUTHOR TO COMPLETE: physical motivation of B4 and its role in the B4-BEC-Lambda-Bridge programme. The project-internal definition and status are summarised in §1.1.]

### 1.1 Context within the project

The B4 Lagrangian, as defined in the repository [B4-README; B4-FW], is

$$
\mathcal L_{\rm B4} = \tfrac12 w^2 - \tfrac12 q^2 - \tfrac{\Lambda}{4} q^4 + \tfrac12 v^2 - \tfrac12 r^2 - \tfrac12 m^2\psi^2 + \tfrac{p}{2} w^2 q^2 - \kappa\, w q v r + \tfrac{s}{2} v^2 r^2,
$$

with $w = \phi_t$, $q = \phi_x$, $v = \psi_t$, $r = \psi_x$; below $\phi \equiv f$ and $\psi \equiv g$. In B4, $p$, $s$ and $\kappa$ are independent parameters, to be treated independently unless a specific UV matching hypothesis relates them [B4-FW]. The Proca construction studied here is the fixed-EFT rank-one prescription introduced in the project's B4.10 stage and registered in the T3-fund Proca specification [SPEC]. It is a specific UV hypothesis of exactly this kind.

Three facts from the project set the scope of the present result:

- **The match is to a subfamily.** A single Proca mediator can reach only the rank-one subfamily $\kappa^2 = ps$ of B4 (§3.5). B4 points with $\kappa^2 \ne ps$ are not reproduced.
- **Not every B4 term comes from the mediator.** The registered Proca action contains the $\Lambda$ term by hand and no $\psi$ mass term. T7a compares only the four-derivative contact sector.
- **Existing verdicts on Proca are unchanged.** The repository documents the Proca mediator as a negative control for the UV-completion question [B4-README; UVCA]: its front speed exceeds the B4 unit cone on all sampled backgrounds. The repository README also states that no candidate reproduces the B4 operator structure. That statement refers to the full B4 structure: independent $p$, $s$, $\kappa$, the $\Lambda$ term and the mass term $m$. T7a is compatible with it. The present result is narrower: it establishes only that the rank-one contact subfamily $\kappa^2 = ps$ is reproduced at leading order, and it classifies the first correction. This is not a reproduction of the full B4 structure, so there is no contradiction between the partial match of T7a and the negative statement about B4 as a whole. T7a does not change any verdict on causality, stability or the UV-completion question, which are evaluated under separate registered protocols [ADD; SPEC; C4].

### 1.2 Method

The second half of T7a is a classification problem. Higher-order operators in EFTs are routinely reduced using the equations of motion or field redefinitions [Politzer 1980; Arzt 1995]. That is legitimate, but it mixes the operator content of the matching with the dynamics of the light fields. Here we deliberately work in an off-shell basis, related to the "Green's basis" of the matching literature [Gherardi et al. 2020]. Operators are identified only up to total derivatives, and equivalence is decided exactly, by two algorithmically independent methods (§6).

## 2. Model and conventions

We work in 1+1 dimensions with metric $\eta = \mathrm{diag}(-1, +1)$ and coordinates $(t, x)$. Subscripts denote partial derivatives. The action is that of the registered T3-fund Proca specification (§1 of `T3_FUND_PROCA_SPEC.md`):

$$
\mathcal L = \tfrac12 f_t^2 - \tfrac12 f_x^2 - \tfrac{\Lambda}{4} f_x^4 + \tfrac12 g_t^2 - \tfrac12 g_x^2 + \tfrac12 F^2 + \tfrac12 M^2 (A_0^2 - A_1^2) + A_1 J,
$$

with

$$
F \equiv F_{01} = \partial_t A_1 - \partial_x A_0, \qquad J = g_\phi\, f_t f_x + g_\psi\, g_t g_x, \qquad \Lambda \ge 0 .
$$

In these conventions $-\tfrac14 F_{\mu\nu}F^{\mu\nu} = \tfrac12 F_{01}^2$ and $-\tfrac12 M^2 A_\mu A^\mu = \tfrac12 M^2 (A_0^2 - A_1^2)$. The coupling $A_1 J$ equals $A_\mu J^\mu$ with $J^\mu = (0, J)$. The current has a single spatial component and is not conserved; the coupling is therefore not a covariant conserved-current coupling.

To avoid the sign ambiguity of the symbol $\Box$ in this signature, we write

$$
D \equiv \partial_t^2 - \partial_x^2 \quad(\text{so that } \Box = -D \text{ for } (-,+)).
$$

The B4 contact sector is written in the basis

$$
O_1 = (f_t f_x)^2, \qquad O_2 = (f_t f_x)(g_t g_x), \qquad O_3 = (g_t g_x)^2,
$$

$$
\mathcal L_{\rm B4} \supset \tfrac12\left[p\,O_1 - 2\kappa\,O_2 + s\,O_3\right].
$$

The coupling constants $g_\phi$, $g_\psi$ are real, with no restriction on their sign.

## 3. Tree-level integration of the Proca field

### 3.1 Gauss constraint

$A_0$ appears without time derivatives. Its Euler–Lagrange equation is

$$
F_x + M^2 A_0 = 0 \;\Rightarrow\; A_0 = -F_x/M^2, \qquad (M^2 - \partial_x^2) A_0 = -\partial_x\partial_t A_1 .
$$

In Fourier space ($\partial_t \to -i\omega$, $\partial_x \to ik$): $F = \frac{M^2}{M^2+k^2}(-i\omega) A_1$.

### 3.2 Equation for $A_1$ and the exact kernel

The $A_1$ equation is $F_t = J - M^2 A_1$. With $A_0$ eliminated, the quadratic form of $A_1$ has the coefficient

$$
c(\omega, k) = -\frac{M^2 (M^2 + k^2 - \omega^2)}{M^2 + k^2},
$$

so that $A_1 = K J$ with

$$
\boxed{K(\omega, k) = \frac{1 + k^2/M^2}{M^2 + k^2 - \omega^2}}, \qquad\text{i.e.}\qquad (D + M^2) A_1 = \left(1 - \frac{\partial_x^2}{M^2}\right) J .
$$

The factor $(1 - \partial_x^2/M^2)$ comes from the Gauss constraint. It is absent in a naive treatment that drops $A_0$.

### 3.3 Tree-level effective Lagrangian

At fixed source $J$ the $A$ sector is quadratic, and $J$ does not depend on $A$. The stationary value of $\tfrac12 c A_1^2 + A_1 J$ is $-J^2/(2c) = \tfrac12 J K J$. Hence, at tree level,

$$
\mathcal L_{\rm eff} = \tfrac12\, J\, A_1[J].
$$

### 3.4 Expansion in $1/M^2$

$$
M^2 K = 1 + \frac{\omega^2}{M^2} + \frac{\omega^4 - k^2\omega^2}{M^4} + \dots
\;\;\Rightarrow\;\;
A_1 = \frac{J}{M^2} - \frac{\partial_t^2 J}{M^4} + O(M^{-6}).
$$

Through order $M^{-4}$ the correction involves only time derivatives; spatial derivatives first enter at $O(M^{-6})$. Substituting and integrating by parts once:

$$
\mathcal L_{\rm eff} = \frac{J^2}{2M^2} + \frac{(\partial_t J)^2}{2M^4} + O(M^{-6}).
$$

### 3.5 Leading order: matching to B4

Expanding $J^2/(2M^2)$:

$$
\frac{J^2}{2M^2} = \tfrac12 \frac{g_\phi^2}{M^2} O_1 + \frac{g_\phi g_\psi}{M^2} O_2 + \tfrac12 \frac{g_\psi^2}{M^2} O_3 .
$$

Comparison with $\tfrac12[pO_1 - 2\kappa O_2 + sO_3]$ gives

$$
p = \frac{g_\phi^2}{M^2}, \qquad s = \frac{g_\psi^2}{M^2}, \qquad \kappa = -\frac{g_\phi g_\psi}{M^2}, \qquad \kappa^2 = ps .
$$

The rank-one relation is therefore automatic. Conversely, a single Proca mediator at tree level reaches only the rank-one subfamily $\kappa^2 = ps$ of the B4 parameter space; B4 with independent $p$, $s$, $\kappa$ is not reproduced. Since $O_3$ is even in $g_\psi$, the relative sign of the couplings is physical. It defines two distinct regimes:

| Regime | Sign condition | $\kappa$ |
|---|---|---|
| rank1_minus | $g_\phi g_\psi > 0$ | $-\sqrt{ps}$ |
| rank1_plus | $g_\phi g_\psi < 0$ | $+\sqrt{ps}$ |

**Remark on what is new.** The leading-order identification of $p$, $s$, $\kappa$ with $g_\phi$, $g_\psi$, $M$ is the defining prescription of the construction [SPEC; B4-FW], so the leading-order match is partly by construction. What T7a adds is the following: the derivation from the full action with the Gauss constraint retained; the exact kernel $K$; the sign of $\kappa$ traced through the action; and the classification of the first correction (§§4–6).

**Remark on the sign of $\kappa$.** The sign follows from the action together with the $-2\kappa$ convention of the B4 operator. It was not adjusted. An earlier internal draft carried the opposite sign; this was corrected against the registered specification before the checks reported here.

### 3.6 Symbolic verification

Every step of §§3.1–3.5 is checked by the script `t7a_step3_symbolic_check.py` (SymPy). It derives the equations from the full action, not from a pre-reduced system. The script evaluates 19 named expressions that must vanish identically:

- conventions: `conv_F2`, `conv_mass`;
- Gauss constraint and the form of $A_1$: `gauss_form`, `A1_form`, `F_in_A1`;
- the kernel: `K_claim`;
- couplings: `p`, `s`, `cO2`, `kappa_is_minus_gg_over_M2`, `kappa2_eq_ps`;
- the four sign quadrants: `regime_minus_pp`, `regime_minus_mm`, `regime_plus_pm`, `regime_plus_mp`;
- expansion coefficients through $M^{-4}$: `series_M0`, `series_M2`;
- on-shell identities: `K_eq_minus_inv_c`, `onshell_half_JKJ`.

All 19 vanish (`ALL_ZERO`). The verification is symbolic. No numerical verification under a protocol has been carried out.

## 4. Effective operators through $M^{-4}$

### 4.1 The correction is outside the B4 basis

Integration by parts preserves, monomial by monomial, the number of fields, the number of $\partial_t$ and the number of $\partial_x$. We write $(N_t, N_x)$ for the derivative class:

- $O_1$, $O_2$, $O_3$ are four-field monomials of class $(N_t, N_x) = (2, 2)$;
- every sector of

$$
(\partial_t J)^2 = g_\phi^2\,[\partial_t(f_tf_x)]^2 + 2 g_\phi g_\psi\, \partial_t(f_tf_x)\,\partial_t(g_tg_x) + g_\psi^2\,[\partial_t(g_tg_x)]^2
$$

is a four-field polynomial of class $(4, 2)$.

Consequently $(\partial_t J)^2$ is not IBP-equivalent to any element of $\mathrm{span}\{O_1, O_2, O_3\}$. It must be reported as a residual operator. Independently of this counting, the term is not a total derivative: for real fields with vanishing boundary terms, $\int (\partial_t J)^2 > 0$ whenever $\partial_t J \not\equiv 0$.

### 4.2 Sector decomposition

We define the coupling-free sector operators

$$
X_{f^4} = [\partial_t(f_tf_x)]^2, \qquad X_{f^2g^2} = \partial_t(f_tf_x)\,\partial_t(g_tg_x), \qquad X_{g^4} = [\partial_t(g_tg_x)]^2,
$$

$$
\Delta\mathcal L^{(M^{-4})} = \frac{g_\phi^2}{2M^4} X_{f^4} + \frac{g_\phi g_\psi}{M^4} X_{f^2g^2} + \frac{g_\psi^2}{2M^4} X_{g^4}.
$$

IBP preserves field content, so the three sectors are classified separately.

## 5. Off-shell operator classification

### 5.1 Equivalence rules (option A)

The classification uses the following rules, fixed before any computation (`T7_PROCA_PROTOCOL_DRAFT.md` §2, v0.6):

1. **Integration by parts is allowed.** Fields have compact support, or satisfy boundary conditions under which the relevant boundary terms vanish. Periodicity in $x$ is not assumed.
2. **Field redefinitions of the light fields are forbidden.** A local redefinition $f \to f + M^{-2} G[f, g]$ produces $\Delta\mathcal L = M^{-2} G\, E_f + O(M^{-4})$, with $E_f = \delta S/\delta f$. It is therefore equivalent to using the light-field equations at first order, and the two cannot be treated differently.
3. **Light-field equations of motion are forbidden** for removing or transforming operators. The heavy-field equation and the Gauss constraint are used; they constitute the tree-level integration itself.
4. **The basis is off-shell and IBP-independent.** Two operators are equivalent if and only if their difference is a total derivative under rule 1.

An on-shell reduction ("option B") is legitimate but is a different problem. It requires an explicitly defined on-shell equivalence, including the $\Lambda$-dependent terms generated through $E_f$. It is not mixed into this analysis.

### 5.2 Operator space

For each sector $S \in \{f^4, f^2g^2, g^4\}$, $V_S$ is the space of four-field monomials with field content $S$ and derivative class $(4, 2)$. Each factor has the form $\partial_t^a \partial_x^b \varphi$ with $a + b \ge 1$, $\varphi \in \{f, g\}$, and monomials are commutative products.

The restriction $a + b \ge 1$ reflects the shift symmetry $f \to f + c$, $g \to g + c$ of the UV action. Total derivatives that land in $V_S$ but come from preimages with underived fields are still accounted for: implicitly by the Euler criterion (§5.3), and explicitly in the second implementation (§6.2).

| Sector | $\lvert V_S\rvert$ |
|---|---|
| $f^4$ | 7 |
| $f^2g^2$ | 26 |
| $g^4$ | 7 |

### 5.3 Equivalence criterion

Two elements of $V_S$ are IBP-equivalent if and only if their difference is annihilated by the Euler operators

$$
E_\varphi = \sum_{a,b} (-D_t)^a (-D_x)^b \frac{\partial}{\partial(\partial_t^a\partial_x^b\varphi)}, \qquad \varphi \in \{f, g\}.
$$

This is exactness of the variational complex for polynomial Lagrangians without explicit $(t, x)$ dependence [Olver 1993]. Hence $V_S/\!\sim\; \cong\; \mathrm{Im}(E|_{V_S})$.

### 5.4 Canonical representatives

Monomials are ordered lexicographically by the key $(k_1, k_2, k_3)$:

- $k_1 = \max_i (a_i + b_i)$: the highest derivative order on a single factor (lower first);
- $k_2 = \sum_i (a_i + b_i)^2$: the concentration of derivatives (lower first);
- $k_3$: each factor is written as a triple $(\text{field}, a, b)$ with $f < g$. The four triples are sorted in ascending order, with repetitions kept, and the resulting lists are compared element by element.

The basis $B_S$ is built greedily. Monomials are traversed in key order, and a monomial is accepted if and only if its Euler image is linearly independent of the images already accepted. Linear independence is tested in a fixed common basis of jet monomials, with exact rational arithmetic. The coordinates $c_k$ of $X_S$ are the unique solution of

$$
E(X_S) = \sum_k c_k\, E(b_k),
$$

and the canonical representative is $R_S = \sum_k c_k b_k$.

## 6. Algorithmically independent cross-check

**What "independent" means here.** The two methods are independent as *algorithms and implementations*:
- they use different representations and different arithmetic libraries;
- they share no code;
- they decide equivalence by different mechanisms: the kernel of the Euler operator versus an explicitly constructed span of total derivatives.

They are **not** independent as *mathematical questions*. Both decide the same equivalence (difference is a total derivative), on the same operator spaces $V_S$, with the same ordering rule and the same greedy selection, all fixed in common approved specifications. Agreement therefore establishes that the classification is implemented correctly and is not an artefact of either algorithm. It does not independently validate the choice of equivalence (option A), of the space $V_S$, or of the canonical ordering. Those are conventions, fixed in advance and reported as such. Agreement is also not a numerical or physical verification in the sense of the draft T7 protocol.

### 6.1 Method I: Euler operator (SymPy)

Script `t7a_step4_classify.py`. Jet variables are plain commutative symbols, not SymPy `Function`/`Derivative` objects. The Euler matrix is built over $\mathbb Q$; the greedy selection uses an exact reduced echelon form.

The verdict is `ALL_CHECKS_PASS`. The controls are:

- **C1–C3:** known total derivatives ($D_t$, $D_x^2$, mixed sector) give zero coordinates.
- **C4:** each basis monomial gives its unit vector.
- **C5:** $f \leftrightarrow g$ symmetry, compared as IBP classes rather than raw coordinates.
- **C6:** the in-house Euler operator agrees with `sympy.calculus.euler.euler_equations`. A defect in the first version of this check was found and repaired: `euler_equations` omits the equation of a field absent from the Lagrangian. The repair was confirmed by a mutation test.
- **C7a:** $\int(\partial_t J)^2 > 0$ in both sign regimes. The integral is exact, over Gaussian-damped polynomial test fields.
- **C7b:** $\int X_S = \int R_S$ exactly for two test-field pairs. This is a necessary condition, not a sufficient one, and it is not an independent implementation.

### 6.2 Method II: explicit IBP relation matrix (pure Python, exact rationals)

Script `t7a_step4b_ibp_matrix.py`. It uses no SymPy, no Euler operator and no floating point, and it imports nothing from Method I. Monomials are represented combinatorially, and total derivatives act directly on the sorted factor tuples.

**Spaces.** $O_S$ is the space of class-$(4,2)$ monomials of sector $S$ with at least one underived factor, and $U_S = V_S \cup O_S$. Relations are $D_t Q$ for all $Q$ of class $(3,2)$ and $D_x Q$ for all $Q$ of class $(4,1)$, with field content $S$, **including preimages with underived fields**.

**Completeness (grading lemma).** $D_t$ and $D_x$ preserve field number and content and raise $N_t$, respectively $N_x$, by one. Any total divergence $Y = D_t A + D_x B \in U_S$ therefore receives contributions only from the class-$(3,2)$ component of $A$ and the class-$(4,1)$ component of $B$. All other components cancel among themselves. Hence the relations span all total derivatives in $U_S$.

**Isolation of $W_S$.** $W_S$ is the subspace of total derivatives lying in $V_S$. It is obtained from the reduced row echelon form of $T^{\mathsf T}$, where $T$ is the relation matrix with the $O_S$ coordinates ordered first. Rows whose pivot falls in the $V_S$ block have a vanishing $O_S$ part and form a basis of $W_S$. Then $\dim(V_S/\!\sim) = |V_S| - \dim W_S$.

**Basis and coordinates.** The basis is built with the same greedy rule, with independence now measured modulo $W_S$. This agrees with Method I if and only if $\ker E|_{V_S} = W_S$, which is precisely what the comparison tests.

The verdict is `ALL_CHECKS_PASS`. The controls are:

- **K1–K3:** known total derivatives reduce to zero modulo $W_S$.
- **K4:** each basis monomial gives its unit vector.
- **K5:** every basis vector of $W_S$ has a zero $O_S$ part.
- **K6a:** two independent enumerations of the preimages agree, and preimages with underived factors are present.
- **K7:** the primary criterion. The SHA-256 of the Method I output is verified first. Then the dimensions, the bases (same monomials, same order) and the exact rational coordinates are required to be **identical** in every sector.

| Sector | $\lvert O_S\rvert$ | relations | rank $T$ | $\dim W_S$ | $\dim(V_S/\!\sim)$ |
|---|---|---|---|---|---|
| $f^4$ | 18 | 26 | 19 | 4 | 3 |
| $f^2g^2$ | 81 | 106 | 80 | 14 | 12 |
| $g^4$ | 18 | 26 | 19 | 4 | 3 |

Agreement between the methods corresponds to $\ker E|_{V_S} = W_S$, i.e. to the statement that the Euler criterion and the explicit relation span define the same subspace. Here this equality is established by two facts. First, the dimensions agree in every sector: $|V_S| - \operatorname{rank} E|_{V_S} = \dim W_S$ (4, 14, 4). Second, the full greedy bases are identical (same monomials, same order), so every accept/reject decision in the ordered traversal of $V_S$ coincides. Together with the standard inclusion $W_S \subseteq \ker E|_{V_S}$ (the Euler operator annihilates total derivatives), equal dimensions give the equality of the subspaces. The agreement of the coordinates of $X_S$ alone would not establish it; it is reported as a consequence, not as the argument. The equality is the expected mathematical result (exactness of the variational complex in this polynomial setting). The comparison checks the implementations against it; it does not prove it in general.

Two further observations are reported only, without a pass/fail criterion. Most raw relations leave $V_S$: 21 of 26 in $f^4$ and 90 of 106 in $f^2g^2$. Only suitable combinations stay inside $V_S$. Dropping the preimages with underived fields (ablation K6b) does not change the dimensions here.

## 7. Results and limitations

### 7.1 Results

The two methods agree exactly.

| Sector | $\lvert V_S\rvert$ | $\dim(V_S/\!\sim)$ | Canonical basis $B_S$ (greedy order) | Coordinates of $X_S$ |
|---|---|---|---|---|
| $f^4$ | 7 | 3 | $f_x^2 f_{tt}^2$, $f_x f_t f_{tx} f_{tt}$, $f_{xx} f_t^2 f_{tt}$ | $(1, 2, 1)$ |
| $f^2g^2$ | 26 | 12 | the 12 monomials listed in Appendix A | $(0,0,1,1,0,1,0,0,0,1,0,0)$ |
| $g^4$ | 7 | 3 | as $f^4$ with $f \to g$ | $(1, 2, 1)$ |

The canonical representatives are

$$
R_{f^4} = f_x^2 f_{tt}^2 + 2 f_x f_t f_{tx} f_{tt} + f_{xx} f_t^2 f_{tt},
$$

$$
R_{f^2g^2} = f_x f_{tx} g_t g_{tt} + f_x f_{tt} g_x g_{tt} + f_{xx} f_t g_t g_{tt} + f_t f_{tx} g_x g_{tt},
$$

and $R_{g^4}$ is obtained from $R_{f^4}$ with $f \to g$. The full first correction in the off-shell basis is

$$
\Delta\mathcal L^{(M^{-4})} = \frac{g_\phi^2}{2M^4} R_{f^4} + \frac{g_\phi g_\psi}{M^4} R_{f^2g^2} + \frac{g_\psi^2}{2M^4} R_{g^4} \quad (\text{mod total derivatives}).
$$

**Sign of the mixed sector.** The mixed operator has no definite sign. For the first test-field pair of C7b, $\int X_{f^2g^2} = -2455075\,\pi/2097152 < 0$, while the full $\int(\partial_t J)^2$ remains positive in both sign regimes.

### 7.2 Limitations and scope of the claims

1. **Status.** These are verified draft results (verified in the sense of §6: two algorithmically independent implementations agree exactly); T7a is not registered. No criterion of the draft T7 protocol (PASS-FORMAL, PASS-VERIFIED, FAIL, INCONCLUSIVE) has been evaluated. The leading-order coefficient match of §3.5 is the content that PASS-FORMAL would assess, but no such verdict is claimed.
2. **Scope of `ALL_CHECKS_PASS` and `ALL_ZERO`.** These verdicts cover only the listed checks. They are not proofs beyond that scope.
3. **Arbiter.** The arbiter is the registered action and conventions; SymPy is a tool.
4. **Tree level, 1+1D.** Loop effects and higher dimensions are not addressed.
5. **Order.** The classification stops at $M^{-4}$. The operator structure at $M^{-6}$, where spatial derivatives of $J$ first enter, is not computed.
6. **Off-shell basis only.** The classification is in option A. An on-shell reduction may give a smaller basis, but it is a different question and is not addressed.
7. **Boundary conditions.** Equivalence is modulo total derivatives with vanishing boundary terms. Other domains require explicit boundary terms.
8. **No numerical verification** under a protocol. The invariant dispersion check foreseen in the draft protocol has not been performed.
9. **No statement about validity or unitarity.** In particular, nothing is concluded about the regime in which the $1/M^2$ expansion is reliable, about a unitarity bound, or about the backgrounds on which the EFT is applied. These are the subject of B4.12, which is open.
10. **T7b.** The $\Lambda$ term is introduced by hand in the UV action. The emergence of $\Lambda$ is not addressed.
11. **Rank-one subfamily only.** The match concerns the subfamily $\kappa^2 = ps$. B4 with independent $p$, $s$, $\kappa$, and the $\psi$ mass term, are outside the construction.
12. **Existing project verdicts.** Nothing here alters the documented status of the Proca mediator as a negative control for the UV-completion question, or the registered T3-fund results (§1.1).

## 8. Reproducibility

**Repository.** `https://github.com/Kretski/B4-BEC-Lambda-Bridge`, branch `master`, directory `t3fund_proca_b412_t7_drafts/`.

| Commit | Content |
|---|---|
| 61fefff | derivation and symbolic check (step 3, spec v0.5.2) |
| a72c323 | equivalence rules, option A (protocol draft v0.6, §2) |
| e7ab8d4 | canonical-representative rule v0.7.1; classification algorithm v0.8.1 |
| 98780c9 | Method I: Euler classification, outputs |
| 4074199 | Method II specification v0.9 |
| 2fc4d16 | Method II: IBP-matrix confirmation, outputs |

**SHA-256** (LF line endings, i.e. the repository blobs):

| File | SHA-256 |
|---|---|
| `t7a_step3_symbolic_check.py` | `83c4137ee9ccd0022ac29911b3369b62eabaec179cefa92d5980a72fe254613d` |
| `t7a_step4_classify.py` | `6b8b3d5c3ab5afd8bee16bfb4b0e764608dd3aa2a941cd80a02faf30b6399150` |
| `t7a_step4_classification.json` | `aaf07cc060d5508101bf66b711a1a5ba184b34bd0ddcb548a32c4f6e07f57269` |
| `t7a_step4b_ibp_matrix.py` | `2995149f19144f622f2f6e956b6f5e69dd23b076fdd24036a394f190df99633d` |
| `t7a_step4b_ibp.json` | `8dfc388ba9de2140aedf4f8d31921c91502ebf72c9cae91799fd645fae754020` |
| `T7A_STEP4_CLASSIFICATION_ALGORITHM_PROPOSAL.md` (v0.8.1) | `0df0221a16cfd1744b207e9f77e1a12bdf3f72e91a08757edb106abb64fe7eb2` |
| `T7A_STEP4B_IBP_MATRIX_PROPOSAL.md` (v0.9) | `fed152397d8edf6100ef574c69a12d3833d7d16e3678f1d6e3506d61ad13dcca` |

**Environment.** Python 3.13.18; SymPy 1.14.0 (Methods I and step 3); Method II uses the standard library only.

**Commands** (from `t3fund_proca_b412_t7_drafts/`):

```
python -I t7a_step3_symbolic_check.py > t7a_step3_output.txt
python -I t7a_step4_classify.py <outdir>
python -I t7a_step4b_ibp_matrix.py t7a_step4_classification.json <outdir>
```

**Re-execution check (2026-11-15).** The three scripts were re-run from a fresh clone at 2fc4d16 with the environment above. All five outputs are byte-identical to the committed files: `t7a_step3_output.txt`, `t7a_step4_classification.json`, `t7a_step4_output.txt`, `t7a_step4b_ibp.json` and `t7a_step4b_output.txt`.

---

## Appendix A. Canonical basis of the $f^2g^2$ sector (greedy order)

| $k$ | $b_k$ | $c_k$ |
|---|---|---|
| 0 | $f_x f_x g_{tt} g_{tt}$ | 0 |
| 1 | $f_x f_t g_{tx} g_{tt}$ | 0 |
| 2 | $f_x f_{tx} g_t g_{tt}$ | 1 |
| 3 | $f_x f_{tt} g_x g_{tt}$ | 1 |
| 4 | $f_x f_{tt} g_t g_{tx}$ | 0 |
| 5 | $f_{xx} f_t g_t g_{tt}$ | 1 |
| 6 | $f_{xx} f_{tt} g_t g_t$ | 0 |
| 7 | $f_t f_t g_{xx} g_{tt}$ | 0 |
| 8 | $f_t f_t g_{tx} g_{tx}$ | 0 |
| 9 | $f_t f_{tx} g_x g_{tt}$ | 1 |
| 10 | $f_{tx} f_{tt} g_x g_t$ | 0 |
| 11 | $f_{tt} f_{tt} g_x g_x$ | 0 |

## Appendix B. Exact integrals of control C7

The test fields are $f = P_f\, e^{-(t^2+x^2)}$ and $g = P_g\, e^{-(t^2+x^2)}$ with small integer polynomials, as recorded in `t7a_step4_classification.json`.

**C7a** ($\int(\partial_t J)^2$, test pair 1):

| Couplings | Regime | Integral |
|---|---|---|
| $g_\phi = g_\psi = 1$ | rank1_minus | $1342801\pi/131072$ |
| $g_\phi = 1$, $g_\psi = -1$ | rank1_plus | $7826279\pi/524288$ |

**C7b** ($\int X_S = \int R_S$; all values $\times\, \pi/2097152$):

| Pair | $f^4$ | $f^2g^2$ | $g^4$ |
|---|---|---|---|
| 1 | 5958329 | −2455075 | 20436637 |
| 2 | 3710429 | −15975 | 2892225 |

## References

**Literature** (bibliographic data checked 2026-11-15)

- [Politzer 1980] H. D. Politzer, "Power corrections at short distances," Nucl. Phys. B 172 (1980) 349–382. doi:10.1016/0550-3213(80)90172-8.
- [Arzt 1995] C. Arzt, "Reduced effective Lagrangians," Phys. Lett. B 342 (1995) 189–195. doi:10.1016/0370-2693(94)01419-D. arXiv:hep-ph/9304230.
- [Olver 1993] P. J. Olver, *Applications of Lie Groups to Differential Equations*, 2nd ed., Graduate Texts in Mathematics 107, Springer, New York (1993). doi:10.1007/978-1-4612-4350-2.
- [Gherardi et al. 2020] V. Gherardi, D. Marzocca, E. Venturini, "Matching scalar leptoquarks to the SMEFT at one loop," JHEP 07 (2020) 225. doi:10.1007/JHEP07(2020)225. arXiv:2003.12525. (Introduces the term "Green's basis": operators independent under integration by parts but not under the equations of motion.)

**Project documents** (repository `github.com/Kretski/B4-BEC-Lambda-Bridge`; commit of the cited version)

- [B4-README] `README.md`, at 05b9a6b: B4 Lagrangian; status of the UV candidates, including Proca as a negative control.
- [B4-FW] `docs/theoretical_framework.md`, at d31a842: B4 Lagrangian; independence of $p$, $s$, $\kappa$; Proca prescription.
- [UVCA] `docs/UV_CANDIDATE_ASSESSMENT.md`, at a3ce2a1: assessment of UV candidates, including the Proca front-speed results.
- [UVAP] `docs/UV_ACCEPTANCE_PROTOCOL.md`, at 05b9a6b: acceptance criteria T1–T8.
- [ADD] `docs/UV_ACCEPTANCE_PROTOCOL_ADDENDUM_T3_T6.md`: preferred-foliation causality convention registered at a3ce2a1; T3-fund protocol registered at 220bf76; Clarification 1 at d70c1e5.
- [SPEC] `t3fund_proca_spec_bundle/T3_FUND_PROCA_SPEC.md`, registered at e355af2: the Proca action and conventions used here (§1–2 of that document).
- [C4] `t3fund_proca_c4_result_bundle/T3_FUND_PROCA_C4_RESULT.md`, at d3de535: registered C4 result for the Proca mediator.
- [T7-DRAFTS] `t3fund_proca_b412_t7_drafts/`, at the commits listed in §8: the unregistered T7a drafts on which this manuscript is based.
- [ZEN] Archived software release v0.2.0 of the repository, DOI 10.5281/zenodo.24144641, as recorded in `CITATION.cff` (17b360d). It predates T7a and does not contain the results of this manuscript. The DOI resolution was not checked from the drafting environment.
