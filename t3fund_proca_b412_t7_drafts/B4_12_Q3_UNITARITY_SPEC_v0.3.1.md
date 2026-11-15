# B4.12 Q3 — Specification of the 1+1D perturbative-unitarity derivation (proposal v0.3.1)

**Status: PROPOSAL FOR REVIEW. Nothing has been derived or computed.** This specification implements Q3 of the registered B4.12 v0.4 (`B4_12_PROCA_PERTURBATIVITY_AND_UV_CLOSURE_v0.4.md`, 74c2a77; registration record 1e0ad68). The decisions of v0.4 §4 are used unchanged:
- vacuum and EFT contact as the primary case;
- the full coupled-channel matrix and its eigenvalues;
- derivation from S†S = 1 with no transferred constant;
- primary quantity (ii).

It does not change Q2 (registered, `641121c`–`36f17ce`), any registered document, or the B4.12 status (**OPEN**).

**Changes from v0.3 to v0.3.1** (review of v0.3, sha256 25df9c77…, "CONDITIONAL APPROVAL"):
1. §4 step 3.3: the reviewer's clarification is added verbatim, separating the exact bound on H(A) from the leading-order estimate for the tree-level eigenvalues. The v0.3 sentence claiming that a tree-level eigenvalue outside the bound can be reconciled with unitarity only if the next-order terms are at least as large as the tree-level term is withdrawn.
2. §4 step 3.2: the expected form of the second-order relation is stated, and the order of R_inel must be shown in the coupling counting, not assumed.
3. §5: Λ_unit is stated explicitly as a leading-order estimate, defined only up to perturbative accuracy.
4. Change-log wording as recommended by the reviewer.

**Changes from v0.2 to v0.3** (review of v0.2, sha256 7ee5acff…, "REVISION REQUIRED"). Only the substantive derivation in §4 step 3 and the dependent definition in §5 are changed; the change log and the corresponding §9 decision row are updated for traceability:
1. §4 step 3 is rewritten. The contraction condition S_el†S_el ≤ 1 is no longer applied to the first-order truncation S_el ≈ 1 + i·c·a, because there the linear terms cancel. The new route is:
   - the exact relation −i(T − T†) = T†T projected onto the elastic block (3.1);
   - its order-by-order expansion, including the tree-to-next-order relation (3.2);
   - a proposed admissibility criterion for the tree-level eigenvalues, built from exact consequences of step 2, with the approximation stated (3.3);
   - the outcome NO_UNIQUE_CRITERION if no unique criterion can be derived (3.3).
2. §5: Λ_unit is defined only if §4 step 3.3 gives a unique criterion; the other outcomes are listed.
3. Dependent wording only: the row "Reading of Q3" in §9.

All other sections are unchanged from v0.2.

**Changes from v0.1 to v0.2** (review of v0.1, sha256 fccc578c…, "REVISION REQUIRED"):
1. §4 step 3 is rewritten. The four-channel matrix is defined as the **elastic block** of the full S-matrix. Unitarity S†S = 1 is imposed on the full space of open channels. The v0.1 claim that the exact elastic block is a unitary matrix with phase eigenvalues is withdrawn. The absence of inelastic channels is not assumed beyond the order considered.
2. New §1.1: the parameter domain (p > 0, M > 0) and the treatment of the boundary p = 0.
3. New §1.2 and check Q3-C9: the normalisation of the f_R f_L channel, with the absence of a factor 1/√2 to be derived, not assumed.
4. Reviewer decisions D-a and D-b are recorded as accepted (§2, §6, §9).
5. §1 (two-particle sector) and Q3-C8 are reworded to match item 1.

---

## 1. Scope and conventions

- **Dimension and fields:** 1+1D, metric (−,+), D = ∂_t² − ∂_x². Light fields f, g with L_free = −½(∂f)² − ½(∂g)².
- **Vacuum:** f̄ = ḡ = 0 (primary case, v0.4 §4 item 4). Scattering about the B4 background is a separate sub-case and is **not** part of this specification.
- **Units and dimensions:** [f] = [g] = 0 in 1+1D; p, κ have dimension M⁻². The derivation confirms these dimensions; it does not assume them.
- **Kinematics:** massless external particles. One-particle states are right-movers R (ω = k > 0) and left-movers L (ω = −k > 0). A nontrivial two-particle collision requires one R and one L particle, with s_cm = 4E_R E_L > 0. RR and LL states have s_cm = 0 and are excluded (v0.4 §4 item 4).
- **Elastic two-particle subspace:**
  - Take fixed initial momenta (p_R, p_L). Energy–momentum conservation in 1+1D restricts final **two-particle** massless states to the same set of momenta {p_R, p_L}. The species of the R and L particles may change, which includes reflection.
  - The elastic two-particle subspace at fixed momenta is therefore finite-dimensional. The derivation must show this explicitly.
  - This subspace is **not** the full space of final states. States with four or more particles at the same total energy and momentum are kinematically open in general (§4 step 3).
- **Channels:** the basis |a_R b_L⟩ with a, b ∈ {f, g}, i.e. four channels: (f_R f_L), (f_R g_L), (g_R f_L), (g_R g_L). Their normalisation is fixed in §1.2.
- **State normalisation:** relativistic, ⟨k|k′⟩ = 2E (2π) δ(k − k′). Any other choice must be stated, and the result must be independent of it.

### 1.1 Parameter domain

- The Proca matching uses g_φ = M√p and g_ψ = −Mκ/√p, and the contact operator contains κ²/(2p).
- **Domain of the primary result:**
  - M > 0;
  - p > 0;
  - κ real; κ = 0 is allowed and is the decoupling check Q3-C6.
- In this domain the parametrisation (p, κ) ↔ (g_φ, g_ψ) is valid and one-to-one, up to the sign convention of g_ψ. The derivation states this explicitly.
- **The boundary p = 0** lies outside the domain:
  - It is **not** reached by substituting p = 0 into expressions that contain 1/p or 1/√p.
  - If it is discussed at all, it is treated as a limit along an explicitly stated path. Examples are fixed κ, where κ²/p diverges, or fixed κ²/p.
  - Such a limit is reported separately and is not part of Λ_unit.
- Every expression in the outputs carries its domain. Divisions by p, √p or M are listed with the condition under which they are taken, as in Q2 (checked denominators, not only reported).

### 1.2 Normalisation of the two-particle channels

R and L particles have disjoint momentum supports (k > 0 and k < 0). The channel f_R f_L nevertheless contains two particles of the same species, so the derivation must show from the Fock construction why this channel gets no extra factor 1/√2 relative to f_R g_L. A convention is not enough. The derivation must:

1. Define |a_R b_L⟩ = a_a†(p_R) a_b†(p_L)|0⟩ and compute its norm from the commutators, for a = b and for a ≠ b. It must show that with disjoint supports the exchange term vanishes, so the norm is the product of the one-particle norms in both cases.
2. Write the completeness relation on the two-particle sector for identical species, including the 1/2! symmetry factor that comes with integrating both momenta over the whole line.
3. Show that restricting to one R and one L particle (k₁ > 0, k₂ < 0) counts each state once, because the regions (k₁ > 0, k₂ < 0) and (k₁ < 0, k₂ > 0) describe the same state. The two-fold double counting therefore cancels the 1/2!.
4. Conclude whether an extra factor appears in the f_R f_L and g_R g_L channels relative to the mixed channels, and carry the result into the phase-space Jacobian (§4 step 1) and into N(s_cm).

This is check Q3-C9 (§7).

## 2. Contact interaction

The primary operator is the EFT contact from the Proca matching, written as in the registered v0.4 §4 item 2:

  L_contact = (p/2)(f_t f_x)² − κ (f_t f_x)(g_t g_x) + (κ²/2p)(g_t g_x)²,  for p > 0 (§1.1).

The last coefficient is **κ²/(2p)**. It is replaced by s/2 only where the rank-one condition s = κ²/p is imposed explicitly; this specification does not impose it.

**Other 4-point terms:**
- **The Λ term (decision D-a, accepted by the reviewer).** The light action contains −(Λ/4) f_x⁴, which in the vacuum is also a 4-point contact (f_R f_L → f_R f_L). It is a separate interaction of the light theory, not part of the Proca contact operator. It is derived and reported in a **separate subsection** of the outputs:
  - its amplitude;
  - its contribution to a(s_cm) and to the eigenvalues.

  It is **not** added into the primary matrix or into the primary Λ_unit.
- The M⁻⁴ residual operators (R_S from T7a step 4) are reported separately as a correction (v0.4 §4 item 2). They are not part of the primary result.

## 3. Channel construction

- Derive the tree-level 2→2 amplitude matrix 𝒯_ij(s_cm), with i, j over the four channels of §1, from L_contact.
- The **full** 4×4 matrix is built first, including ff ↔ gg transitions via κ and reflection channels, where kinematically allowed. Channels are not treated independently unless the derived matrix is block-diagonal.
- **Two routes for 𝒯 (required):**
  - (A) Feynman rules from L_contact: vertex factor with on-shell external legs;
  - (B) the first-order matrix element of the interaction Hamiltonian between mode-expanded two-particle states, with the normalisation of §1.2.

  They must agree identically.

## 4. Unitarity derivation (from S†S = 1)

1. **Elastic matrix element.**
   - Write S = 1 + iT on the full Hilbert space.
   - Relate the matrix elements of T between the elastic two-particle states of §1 to 𝒯, with the momentum-conserving δ-function.
   - The Jacobian of the 1+1D two-particle phase space for massless R + L must be **derived**, by two methods: direct δ-function integration, and the invariant form. It includes the identical-species counting of §1.2.

2. **Unitarity on the full space of open channels.**
   - S†S = 1 holds on the full space of states. Let P be the projector onto the elastic two-particle subspace at fixed momenta. The elastic block is S_el = P S P, the 4×4 matrix in the channel basis.
   - Inserting the completeness relation into ⟨i|S†S|j⟩ gives two parts:
     - the elastic intermediate states, which give S_el†S_el;
     - all other open states (four or more particles), which give a positive semidefinite remainder R_inel.
   - So S_el†S_el = 1 − R_inel, with R_inel ≥ 0. Hence S_el is a contraction: S_el†S_el ≤ 1.
   - The derivation must show this decomposition explicitly with the normalisation of step 1 and §1.2. The condition used in step 3 must be derived from it.

3. **From exact unitarity to a criterion on the tree-level eigenvalues (rewritten in v0.3).**

   **3.0 What is not claimed.**
   - The four-channel matrix is the **elastic block** of the full S-matrix. It is **not** claimed to be unitary on its own, and its exact eigenvalues are **not** claimed to be phases. At higher orders, inelastic channels (for example 2 → 4) can contribute, so R_inel ≠ 0 in general. It is **not** assumed that inelastic channels are absent at all orders.
   - The contraction condition S_el†S_el ≤ 1 is **not** applied to the first-order truncation S_el ≈ 1 + i·c·a. For Hermitian a, the linear terms cancel and S_el†S_el = 1 + c²a² + …. Applied to the truncation itself, the condition would demand c²a² ≤ 0, i.e. a = 0, which is not a criterion. The contraction condition bounds singular values of the exact elastic block. It gives **no** bound on the eigenvalues of the tree-level matrix directly.

   **3.1 Exact amplitude and the unitarity relation.**
   - Write the exact elastic block as S_el = 1 + i·c·A(s_cm), where A is the exact dimensionless elastic amplitude matrix, in general not Hermitian. The normalisation constant c follows from step 1 and §1.2. The derivation must show that c is real and must state its value; it is not inserted.
   - From S = 1 + iT and S†S = 1 on the full space: −i(T − T†) = T†T. Project this relation onto the elastic subspace with P and insert completeness, separating elastic and other open intermediate states. The result is the exact relation between the anti-Hermitian part of A, A†A and the inelastic remainder R_inel. The derivation must state this relation with all factors, in terms of c.

   **3.2 Perturbative expansion by orders.**
   - Expand A = A⁽¹⁾ + A⁽²⁾ + … in powers of the contact couplings. Tree level: A⁽¹⁾ = a = N·𝒯 from step 1.
   - The derivation must show from the real Lagrangian that a is Hermitian (check Q3-C7).
   - From 3.1, derive order by order:
     - at first order, the relation that fixes the anti-Hermitian part of A⁽¹⁾. It is expected to vanish; this must be derived.
     - at second order, the relation between the anti-Hermitian part of A⁽²⁾ (the leading imaginary part of the elastic amplitude) and a.
   - **Expected form, to be derived:**
     - The exact projected relation is expected to take the form −i(A − A†) = c·A†A + R_inel/c, with the normalisation of step 1.
     - With A = a + A⁽²⁾ + …, a = a†, and no contribution of R_inel at second order, this gives −i(A⁽²⁾ − A⁽²⁾†) = c·a².
     - The derivation must confirm both relations with all factors, or report the form it actually obtains.
   - **Order of R_inel:** it must be **shown** in the coupling counting used, not assumed. The counting must cover:
     - how many contact vertices a 2 → 4 amplitude requires;
     - the order of its square in R_inel;
     - that this order lies beyond the second order used above.
   - This tree-to-next-order relation is a mandatory output (§7). It is derived, not assumed.

   **3.3 Admissibility criterion for the tree-level eigenvalues (proposed route; to be derived and approved).**
   - The proposed route uses only exact consequences of step 2:
     - (i) For a contraction, the numerical range lies in the closed unit disk: |⟨v, S_el v⟩| ≤ 1 for every unit vector v in the elastic subspace. This holds without assuming that S_el is normal or diagonalisable.
     - (ii) Inserting S_el = 1 + i·c·A gives an exact region of the complex plane for ⟨v, A v⟩. Projecting this region onto the real axis bounds ⟨v, H(A) v⟩, where H(A) = (A + A†)/2. That gives a bound on all eigenvalues of the Hermitian part of the **exact** elastic amplitude. The bound depends only on the derived c.
     - (iii) At leading order H(A) = a + O(couplings²), so the extreme values of ⟨v, H(A) v⟩ are approximated by the extreme eigenvalues of a.
   - **Criterion (leading order):** an eigenvalue λ_i of the tree-level matrix a is *admissible* if it lies within the exact bound of (ii), i.e. |λ_i(a)| ≲ 1/|c| with c as derived.
   - **Exact bound versus tree-level estimate:** The exact bound applies to the Hermitian part of the exact amplitude A. Its application to the eigenvalues of the tree-level matrix a is only a leading-order self-consistency estimate, with corrections of order O(g²). A tree-level eigenvalue crossing the bound identifies the breakdown scale only up to perturbative accuracy; it does not establish an exact violation of unitarity or prove that the next-order contribution is as large as the tree-level term at the crossing point.
   - The derivation must establish steps (i)–(iii) explicitly, state the resulting inequality in terms of c, and state the approximation, namely which terms are neglected and at which order.
   - **Nature of the criterion:** it is a self-consistency criterion of perturbation theory, not an exact bound on the tree-level amplitude. The scale it defines is meaningful as an order-of-magnitude scale. Its numerical constant is the one produced by this criterion and no other.
   - **Alternatives and failure mode:**
     - If the derivation shows that steps (i)–(iii) do not hold as stated, this is reported. Examples: c is not real, or a is not Hermitian at tree level.
     - If more than one inequivalent criterion follows from the available relations with different constants, this is reported too.
     - In either case the outcome is recorded as **NO_UNIQUE_CRITERION**, with the reasons. No threshold is introduced by analogy with 3+1D, with partial-wave conventions, or by choosing among alternatives by convention.

4. **Order counting.** Inelastic 2 → 4 amplitudes do not occur at the order of L_contact considered (tree level, four-point). Their first contribution to R_inel is of higher order in the couplings. The derivation states the order explicitly. This is a statement about the order of the bound, not about the exact S-matrix.

## 5. Eigenvalue conditions

- Compute the eigenvalues λ_i(s_cm) of the dimensionless 4×4 tree-level matrix a(s_cm), and their dependence on s_cm.
- **Λ_unit** is defined only if §4 step 3.3 produces a unique admissibility criterion. It is then a **leading-order estimate**: the smallest √s_cm at which some tree-level eigenvalue λ_i(s_cm) first crosses the leading-order criterion. It is defined only up to perturbative accuracy, with corrections of order O(g²) (§4 step 3.3), and is reported with that qualifier. If different eigenvalues fail at different s_cm, all values are reported, and the smallest defines Λ_unit.
- **Outcomes other than a finite Λ_unit:**
  - If §4 step 3.3 gives **NO_UNIQUE_CRITERION**, Λ_unit is **not defined**. The result records this outcome and its reasons, and the Λ_unit field of the outputs (§7) holds that outcome instead of an expression.
  - If the criterion is unique but no eigenvalue ever fails within the domain of §1.1, this is reported as such. An example is eigenvalues that do not grow with s_cm. No finite Λ_unit is then given.
- **Approximation:** Λ_unit is an estimate of the scale at which the tree-level description stops being self-consistent (§4 step 3.3). It is not an exact bound, it does not establish an exact violation of unitarity, and it is not the scale of new physics.
- **Hypothesis H-Q3 (v0.4 §4), not assumed:** for the Proca contact, the coupling structure may give a single non-zero eigenvalue depending on a combination of p and κ. The derivation reports the rank and the eigenvalues as obtained, within the domain of §1.1.
- **Comparison with Λ_* ~ 1/√p:** diagnostic only (v0.4 §4). If Λ_unit depends on a combination such as p + κ²/p, the comparison is reformulated and not forced back to 1/√p.

## 6. Validity and exclusions

- The result concerns **vacuum** two-particle scattering only. It is **not** a bound on scattering about an arbitrary background, and not a statement about the B4 backgrounds.
- **Tree level, contact approximation:** valid only where the EFT expansion holds. The joint statement with Q2 (which scale binds first) is the separate joint Q2/Q3 check and is **not** part of this specification.
- **Full tree-level Proca exchange (decision D-b, accepted by the reviewer).** It contains momentum-dependent structure and is a different problem from the local contact approximation. It is a **separate later specification** and is not part of this one.
- **Elastic block only:** the bound is a perturbative statement about the elastic block (§4 step 3). It is not a statement about the exact S-matrix or about inelastic channels.
- **IR limitation:** massless scalars in 2D have IR problems (Coleman). They are no obstacle for derivative couplings at tree level, and are recorded as a limitation.

## 7. Mandatory checks and outputs

| # | Check | PASS condition |
|---|---|---|
| Q3-C1 | Kinematics | on-shell massless R + L: final two-particle momenta = {p_R, p_L}; RR/LL have s_cm = 0; derived, not assumed |
| Q3-C2 | Free limit | for L_contact = 0: S_el = 1 on the elastic subspace with the chosen normalisation |
| Q3-C3 | Phase space | the 1+1D Jacobian agrees between direct δ-function integration and the invariant form, with the counting of §1.2 |
| Q3-C4 | Amplitude, two routes | 𝒯 from Feynman rules (A) equals 𝒯 from the mode-expanded matrix element (B), identically |
| Q3-C5 | Dimension | a(s_cm) is dimensionless; 𝒯 has the derived dimension |
| Q3-C6 | Decoupling | for κ = 0 (p > 0), the f and g channels decouple (the g–g coefficient κ²/2p also vanishes) and a becomes block-diagonal |
| Q3-C7 | Hermiticity / T-symmetry | a(s_cm) has the symmetry implied by the real Lagrangian at tree level (derived and checked) |
| Q3-C8 | Perturbative consistency | with the derived normalisation, S_el = 1 + i·c·a satisfies S_el†S_el = 1 + O(a²); the O(a²) remainder is identified, not asserted to vanish |
| Q3-C9 | Channel normalisation | norms of \|f_R f_L⟩, \|f_R g_L⟩, \|g_R f_L⟩, \|g_R g_L⟩ computed from the commutators; the 1/2! factor and the R/L double counting cancel explicitly; the presence or absence of a 1/√2 factor is derived (§1.2) |

`ALL_CHECKS_PASS` covers only these checks.

**Outputs:**
- a standalone SymPy script `b412_q3_unitarity.py`;
- `b412_q3_output.txt`;
- `b412_q3_results.json`, containing:
  - 𝒯, N, c and a;
  - the eigenvalues of a and its rank;
  - the derived eigenvalue condition;
  - Λ_unit as a symbolic expression in p and κ, with its domain;
  - the separate Λ-term subsection (D-a);
  - Λ_unit/Λ_* as a symbolic diagnostic;
  - the list of checked denominators;
- a later report, after review.

## 8. Not included

- No numerical evaluation on registered backgrounds or registered parameter points.
- No joint Q2/Q3 conclusion.
- No full Proca exchange (D-b).
- No scattering about the B4 background.
- No claim about the exact S-matrix or about inelastic channels beyond the order considered.
- No change of the B4.12 status, no T7a registration, no Zenodo.

## 9. Recorded decisions

| Item | Decision (review of v0.1) |
|---|---|
| D-a: Λ f_x⁴ | Derived and reported separately; not in the primary matrix or Λ_unit |
| D-b: full Proca exchange | Separate later specification |
| Reading of Q3 | Perturbative unitarity of the tree-level eigenvalues of the elastic block, derived from S†S = 1 on the full space of open channels; no claim that the exact elastic block is unitary |
| Route to the threshold (v0.2 review) | Not directly from S_el†S_el ≤ 1 applied to the tree truncation; via −i(T − T†) = T†T expanded by orders, with an explicit criterion, or NO_UNIQUE_CRITERION |

**Requested in this version:** approval of the proposed route in §4 step 3.3 (numerical range of the contraction → bound on the Hermitian part of the exact elastic amplitude → tree-level approximation). If the route is not accepted, the specification can keep only §4 steps 3.1–3.2 and report the criterion as an open item. Execution remains forbidden until this specification is approved and registered.
