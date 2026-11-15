# B4.12 Q3 — Analytic memorandum on F3/F4 (proposal v0.2)

**Status: PROPOSAL FOR REVIEW.**
- This memorandum interprets the registered Q3 run. It introduces **no new computation**.
- Every number below is read from the registered outputs: `results_b412_q3/b412_q3_results.json`, commit `df4d8b9`, SHA-256 `8bbdd1d0…`, registered in `B4_12_Q3_RUN_REGISTRATION.md`, commit `17f53fe`.
- Where a line rearranges those numbers, the elementary algebra is shown and marked **(algebra)**.
- No option among O1, O2 and O3 is selected. Q3 and B4.12 remain **OPEN**.

**Changes from v0.1** (review of v0.1, sha256 ffca1b28…, not approved for commit). There are no new computations and no changes to registered results. The wording follows the reviewer's proposals:
1. A1 is restricted to a comparison of corresponding coupling structures.
2. A2 names p σ² as the parameter of the f → f f f channel. It introduces ε_f for an initial f packet and ε_el for the elastic block, and limits the "long-time" reading to the first-order setting.
3. A3 states that T_int and t correspond to each other but are not identified. It also states that the asymptotic-state reading motivates, but does not establish, a particular construction.
4. A4 is limited to the derived boost structure.
5. §5 (O1/O2/O3) is replaced by the reviewer's wording.

---

## 1. Setting

The results hold under all of the following:
- vacuum;
- massless f and g in 1+1D;
- free massless Fock asymptotic states;
- first order in the contact couplings.

L_contact = (p/2)(f_t f_x)² − κ(f_t f_x)(g_t g_x) + (κ²/2p)(g_t g_x)², with p > 0. The Λ term is reported separately. The normalisation convention is c = 1, a = 𝒯/(2s_cm).

## 2. What is established

**(E1) Operator split (F1).**
- L_contact = L_RL + L_ch exactly.
- L_RL = −p R_f²L_f² + κ(R_f²L_g² + R_g²L_f²) − (κ²/p) R_g²L_g² has boost weight 0.
- L_ch = (p/2)(R_f⁴ + L_f⁴) − κ(R_f²R_g² + L_f²L_g²) + (κ²/2p)(R_g⁴ + L_g⁴) has boost weights ±4.
- The operator is parity invariant.

**(E2) Elastic block.**
- Only L_RL contributes to the first-order elastic block (C4).
- The elastic matrix is diagonal in the channel basis, with eigenvalues of 𝒯 equal to −p s²/4, κ s²/4 (twice) and −κ² s²/(4p). Here s ≡ s_cm.

**(E3) First-order transitions (F3).**
- From |a_R b_L⟩ the only allowed first-order transitions are three:
  1. the elastic one;
  2. R → R R R with an L spectator;
  3. L → L L L with an R spectator.
- The full matrix elements with the spectator factorise as i·T₁₃·4π·T_int·δ(…)·δ_bb′·2|q|·2π·δ(q′ − q). This holds for all 10 nonzero species combinations.

**(E4) Splitting vertices (F3, route A = route B).** With K = k₁ + k₂ + k₃:

| Channel | T₁₃ / (K k₁k₂k₃) |
|---|---|
| f → f f f | −12 p |
| f → f g g | 4 κ |
| g → f f g | 4 κ |
| g → g g g | −12 κ²/p |
| Λ term: f → f f f (separate) | 6 Λ |

**(E5) Secular probabilities (F4).**
- For the packet φ(k) = √(8π/σ)·√k·e^{−k/σ}, at first order, P = C·t² with:

| Channel | C |
|---|---|
| f → f f f | 9 p² σ⁶ / (256 π²) |
| f → f g g | 3 κ² σ⁶ / (256 π²) |
| g → f f g | 3 κ² σ⁶ / (256 π²) |
| g → g g g | 9 κ⁴ σ⁶ / (256 π² p²) |
| Λ term: f → f f f (separate) | 9 Λ² σ⁶ / (1024 π²) |

- For the normalised product packet ψ_R × ψ_L, the sum over the final spectator equals 1. So P(R L → R R R L) = P(R → R R R) at first order.
- ESTABLISHED applies to generic parameter values. The κ channels vanish at κ = 0 and the Λ channel at Λ = 0.

## 3. Analytic reading

Each item states what follows from E1–E5 and how far it goes.

**(A1) Corresponding coupling structures.**
- Compare E2 with E4, sector by sector:
  - f sector: elastic −p, splitting −12p;
  - mixed: elastic κ, splitting 4κ;
  - g sector: elastic −κ²/p, splitting −12κ²/p.

  The factors are combinatorial, for example 12p = (p/2)·4! and 4κ = κ·2·2 **(algebra)**.
- For each corresponding coupling structure, the elastic entry and the chiral splitting vertex contain the same coupling coefficient, up to fixed combinatorial factors. There is no independent coupling parameter within L_contact that suppresses a given chiral vertex relative to its corresponding elastic interaction.
- This is a statement about coupling dependence, not about equality of the physical amplitudes or probabilities. For example, as κ → 0 the f-sector entries remain while the mixed and g-sector entries vanish. The comparison is channel by channel, not between arbitrary elements of the operator.

**(A2) Time scale and the relevant parameters (first-order setting).**
- **Parameter of the f → f f f channel:** ε_σ := p σ². Then P(f → f f f) = (9/256π²)·ε_σ²·(σt)² **(algebra)**. It is not a general EFT parameter of the whole operator, because κ is independent.
- **Initial f packet (both channels f → f f f and f → f g g):**
  - P_f(t) = (9p² + 3κ²) σ⁶ t²/(256π²).
  - With ε_f := σ² √(9p² + 3κ²), this gives P_f(t) = ε_f² (σt)²/(256π²) and σ t_f = 16π/ε_f, i.e. t_f = 16π / (σ³ √(9p² + 3κ²)) **(algebra)**.
- **Initial g packet:** t_g = 16π / (σ³ √(3κ² + 9κ⁴/p²)) **(algebra)**.
- **Elastic block:** at s_cm ~ σ² it is characterised by ε_el := σ² max(p, κ²/p). This also covers the mixed entry, since |κ| ≤ max(p, κ²/p) for p > 0. The conditional leading-order estimate Λ_RL corresponds to ε_el ≲ O(1). That is an energy condition.
- **Reading:** for small ε_f (or the analogous g parameter), the first-order splitting probability nevertheless reaches O(1) when σt ~ 16π/ε_f. Within the first-order setting and for small dimensionless coupling parameters, the secular limitation is therefore a long-time limitation, distinct from the energy condition of the elastic block.
- This does not show that the EFT has no separate energy limitation.
- The times t_f and t_g only mark where first-order perturbation theory around free Fock states stops being self-consistent for this packet. They are not physical lifetimes, and they say nothing about the exact evolution.

**(A3) Why it is secular.**
- Collinear massless states are exactly degenerate.
- An elastic R L collision of two packets lasts only a finite crossing time.
- The chiral term instead acts during the entire free propagation, because collinear particles never separate.
- The plane-wave T_int in F3 corresponds to the factor t in the finite-time packet calculation of F4 (registered `secular_link`). Both encode the same secular degeneracy in their respective formulations. Their precise identification requires a common regulator prescription and is not established here.
- This motivates examining the asymptotic-state construction. It does not establish that a particular alternative construction is necessary or sufficient. Other routes, such as a reorganisation of the perturbative expansion or a resummation, are not excluded.

**(A4) Boost structure.**
- L_RL has boost weight 0 (E1), and the elastic results and Λ_RL depend on s only.
- Since L_ch contains components with boost weights +4 and −4, the packet-dependent secular coefficients are not Lorentz scalars. The reported coefficients refer to the frame and packet parametrisation used in F4. No explicit transformation law for the packet probability under boosts is derived here.

**(A5) The conditional L_RL estimate.**
- Λ_RL² = 8/max(p, κ²/p) for κ ≠ 0 (registered).
- The two terms are equal at κ² = p², a special point. It is not the Q2 rank-one condition κ² = ps. The Q3 operator already has s = κ²/p, so equality needs s = p as well.
- The estimate is conditional on a prescription of type O1 or O2 that leaves the first-order RL amplitude unchanged.

## 4. What does not follow

- No statement about the exact interacting evolution, higher orders, or other initial states.
- No statement that no asymptotic construction or scattering description exists. No statement that the theory is non-unitary.
- No Λ_unit for the full L_contact. The registered outcome remains NO_UNIQUE_CRITERION, for the present setting.
- **Zero modes.** The zero-mode groups are only candidates for an additional chiral sector. F2 is a symbolic momentum-space test without the discrete lattice.
- **Relation to Q2.** It is not claimed that the chiral terms here and the background-dependent sound speeds of Q2 are the same physical effect. The memorandum notes only that both come from the same operator structure. A joint analysis is needed before any statement.
- **The Λ term.** It has its own chiral part and secular channel. No threshold is computed for it.

## 5. Consequences for the open decisions (no choice made)

- **O1:** A1–A3 motivate investigating whether the chiral interactions must be incorporated into the asymptotic Hamiltonian. The present calculation does not establish that a non-perturbative construction is necessary or sufficient.
- **O2:** A reduced S-matrix prescription would require a separate analysis of the relevant commutators, including terms of the type [Q_R, L_RL], and of their perturbative order. These terms are not calculated in the registered run, and A1 does not determine their magnitude.
- **O3:** Restricting the conclusion to the present free-Fock setting is the option directly supported by the registered result. It does not rule out O1 or O2.
- The zero-mode task (a full discrete operator construction in the periodic box) remains prior to a decision among the asymptotic prescriptions.

## 6. Main message

1. There are first-order collinear transitions with a spectator.
2. For the stated normalised packet they produce secular probabilities growing as t².
3. This limits the perturbative description around free massless Fock states at long times.
4. It does not prove that the theory is non-unitary, and it does not select an asymptotic setting.
5. The zero modes remain a separate, unresolved task.

## 7. Proposed use

If approved, this memorandum would be committed as an interpretation document next to the registered run. The run artefacts and registration records stay unchanged. It would form the basis for a short Q3 status report. The status report is not a closure of Q3.
