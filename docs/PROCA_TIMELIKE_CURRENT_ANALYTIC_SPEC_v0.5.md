# Proca with a timelike current: analytic review against B4.12 — specification v0.5

**Status:** PROPOSAL, NOT REGISTERED. Nothing is derived or computed before this specification is reviewed and registered.

This specification does not change:

- B4.12 (v0.4, registered; status OPEN);
- the common-M specification and result (3cba04e and follow-ups);
- the registered T3-fund Proca specification (`t3fund_proca_spec_bundle/T3_FUND_PROCA_SPEC.md`, e355af2);
- the T7a commits or the frozen T7a manuscript v0.3;
- `docs/UV_STRUCTURAL_CONSTRAINTS_NOTE_v0.3.md` (c7b92ba).

**Evidence level:** 1 (theoretical). The work is analytic, with symbolic checks only. There is no population, no numerical run and no new threshold.

**Changes from v0.4** (review of 2026-11-18, third round). Only two local corrections:

1. D3: an explicit rule for ‖K_lead‖₂ = 0 (R-UNDEFINED).
2. §3: the case "k_der = 0 on some applicable backgrounds and > M on the rest" is listed explicitly under OPEN, in the list and in the table.

Nothing else changed.

**Changes from v0.3** (review of 2026-11-18, second round; kept). Only two changes:

1. D3: K_exact and K_lead are now defined as the full quadratic matrices of L_light + L_exchange.
2. §3: the classes are made mutually exclusive by an ordered decision rule, with the k_der table of the review. To make "applicable backgrounds" decidable, the Stage 2 arithmetic now also covers evaluating the D3 formula (§4). No new computation type is added.

Nothing else changed.

**Changes from v0.2** (review of 2026-11-18, first round; kept). Only D3, the scale name in §3, and two small precision points requested in the same review were changed. The model, D0–D2, D5's procedure, the stages and the scope are unchanged.

1. D3: the remainder is the full K_rem = K_exact − K_lead. Δ_μν − η_μν/M² = (k_μk_ν − η_μν k²)/(M²(k² + M²)) is used. The pure k_μk_ν part is kept only as a separate diagnostic.
2. D3: the quadratic form is the full second variation of the effective action, including the background term ½J̄Δ δ²J.
3. D3 and §3(c): the scale is renamed k_der ("loss of dominance of the leading local term"). It is not called a strong-coupling scale. §3(c) is reworded accordingly.
4. D4: |γ_u| ≤ g_max, not γ_u.
5. D5: reduced-system hyperbolicity and physical hyperbolicity after the constraint are kept apart.

(Changes from v0.1, applied in v0.2 and kept: exact Lagrangian, sign and J_n in §0–§1; three-way IBP classification in D2; G′ = G in D4; coupled-system procedure in D5.)

## 0. Conventions (fixed)

- **Metric (−,+)**, η = diag(−1, 1), as in the registered T3-fund Proca specification §1 and T7a step 3.
- **Frame.** The frame of u is the lab frame of T7a: u^μ = (1, 0), u_μ = (−1, 0), n^μ = (0, 1). There X₁ = f_t, X₂ = f_x, X₃ = g_t, X₄ = g_x, with f ≡ φ and g ≡ ψ as in the registered specification.
- **Proca sector**, exactly as registered:

  L_Proca = ½F² + ½M²(A0² − A1²),   F ≡ A1_t − A0_x.

  This equals −¼F_{μν}F^{μν} − ½M²A_μA^μ in (−,+). Here A0 and A1 are the components A_μ with lower index. A0 has no time derivative (Gauss constraint).
- **Coupling:** +A_μJ^μ = A0·J^0 + A1·J^1.
- **B4 operator basis and normalisation:** as in `docs/notation_and_conventions.md` and B4.12 v0.4 §0. In particular L_B4 ⊃ −(Λ/4)f_x⁴.

## 1. Model (fixed before any derivation)

L = ½f_t² − ½f_x² + ½g_t² − ½g_x² + ½F² + ½M²(A0² − A1²) + A1·J_n + A0·J_u

- **J^1 = J_n = g_φ f_t f_x + g_ψ g_t g_x.** This is the T7a current, unchanged.
- **J^0 = J_u = γ_u f_x²**, i.e. J^μ ⊃ γ_u u^μ (n·∂φ)². This is the only addition. γ_u is real.
- **Bare Λ = 0.** The term −(Λ/4)f_x⁴ of the registered specification is removed. This is the T7b condition: any Λ term must come from the exchange.
- There is no other coupling, no contact term, and no coupling of J_u to g.
- **Background:** constant jets f = wt + qx, g = vt + rx, as registered. The Proca background is to be derived from the four Euler–Lagrange equations, not assumed. The expected values are Ā1 = J̄_n/M², Ā0 = −J̄_u/M², F̄ = 0.

## 2. Derivations

Each item gives the quantity to derive and a hypothesis to check. A hypothesis is not a result.

**D0 — sign at leading order.**

- From §1 with J_n = 0, eliminate A0 at leading order: M²A0 + J_u = 0 with the F terms dropped.
- Show that L_eff ⊃ −J_u²/(2M²) = −(γ_u²/2M²) f_x⁴. Comparing with −(Λ/4)f_x⁴ gives Λ_eff = 2γ_u²/M² > 0 for every real γ_u ≠ 0.
- Confirm that this is the same physical sign as Statement C2 in (+,−).

**D1 — leading local order, O(M⁻²).**

- Derive L_eff in full for J^μ = n^μJ_n + u^μJ_u.
- **H1:** the cross terms J_nJ_u vanish because u·n = 0. If so:
  - p, s and κ are exactly those of T7a, with the T7a normalisation (p = g_φ²/M², s = g_ψ²/M², κ = −g_φg_ψ/M²; rank one);
  - the only new operator is the Λ term;
  - no operator outside the B4 set is created at this order.
- **Consistency check for D1:** on the background, the coefficient of δf_x² in the second variation must equal that of B4 with Λ = Λ_eff, i.e. −(3/2)Λ_eff q² = −3γ_u²q²/M². It must be assembled from three sources: the background value Ā0 contracted with δ²J_u; the exchange of δA0; and the T7a terms.

**D2 — next order, O(∂²/M⁴).**

- Derive the O(∂²/M⁴) terms of L_eff, including:
  - derivatives on J_n;
  - derivatives on J_u;
  - mixed terms between J_n and J_u;
  - terms in ∂_μJ^μ = ∂_tJ_u + ∂_xJ_n.
- Do not assume that any of them can be removed by integration by parts.
- Use the approved T7a rules (variant A): IBP only; no field redefinitions; no light-field equations of motion. Use the approved T7a step-4 canonical basis and IBP matrix, unchanged.
- Sort every term into exactly one of three classes:
  1. **T7a-equivalent:** equal, after IBP, to a term already in the T7a M⁻⁴ classification;
  2. **new:** not contained in T7a;
  3. **undecided:** its effect on B4.12 cannot be decided analytically. It is listed with the reason.

**D3 — longitudinal structure and the scale k_der of the local expansion.**

- In 1+1D a Proca field has one physical degree of freedom, and A0 is fixed by the Gauss constraint.
- Derive the modified Gauss constraint with J_u, and show how J_u enters only through it.
- The Stückelberg form A_μ → A_μ + ∂_μπ/M, with coupling −(π/M)∂_μJ^μ up to a boundary term, is used **for interpretation only**. It is not used in any definition. All definitions use the exact tree-level elimination of (A0, A1). This needs no gauge fixing, because Proca has no gauge symmetry.

**Definitions (fixed):**

- **Background:** a fixed constant-gradient background Φ̄ = (w, q, v, r) with the Proca background of §1. It is homogeneous, so plane waves are exact. No WKB or local approximation is used.
- **Fluctuations:** light-field plane waves (δf, δg) ∝ e^{i(kx − ωt)}, with k the spatial wavenumber in the frame of u.
- **Exact quadratic form.** K_exact(ω, k; Φ̄) is the 2×2 matrix of the **full second variation of L_light + L_exchange** for (δf, δg). Here:
  - L_light = ½f_t² − ½f_x² + ½g_t² − ½g_x² is the free light part, with bare Λ = 0;
  - L_exchange is the tree-level effective action from the exact elimination of (A0, A1).

  Its free-light part is ½(ω² − k²) on each diagonal entry. With J = J̄ + δJ + ½δ²J + …, the exchange part contains both:
  - ½ δJ^μ Δ_μν(ω, k) δJ^ν, and
  - the background term ½ J̄^μ Δ_μν δ²J^ν.

  Δ_μν is the exact Proca propagator in (−,+):

  Δ_μν = (η_μν + k_μk_ν/M²)/(k_ρk^ρ + M²),   k_μ = (−ω, k).

  Because J̄ is constant, the background term picks Δ at zero momentum. The expectation is that it equals ½J̄·δ²J/M² exactly and contributes only to the leading local part. This is to be derived, not assumed.
- **Leading part.** K_lead(ω, k; Φ̄) is the 2×2 matrix of the full second variation of L_light + L_exchange^(lead), where L_exchange^(lead) is L_exchange with Δ_μν replaced by its leading local term η_μν/M². So K_lead contains:
  - the free light part;
  - the local exchange, i.e. the leading-order B4 quadratic form with Λ = Λ_eff.

  By construction, det K_lead(c_i k, k) = 0 on the leading-order branches. The norm used below is the spectral norm. It is usually not zero there, but it can be. For example, on a zero background with free light fields the branches are ω = ±k, and K_lead can vanish on them.
- **Rule for ‖K_lead‖₂ = 0 (fixed):** "zero" means exactly zero, decided symbolically, with the registered values rationalised as in Script A. No numerical tolerance is introduced.
  - **Identically zero on a branch.** If K_lead(c_i k, k; Φ̄) = 0 for all k > 0 on a real branch, that branch is marked **R-UNDEFINED**. If any branch of Φ̄ is R-UNDEFINED, the background is marked **R-UNDEFINED**, and k_der(Φ̄) is **undetermined**. The exception is when Stage 1 derives a separate analytic limit for that branch, stated explicitly and not assumed.
  - **Zero at isolated k.** If K_lead = 0 only at isolated k on a branch, those points are reported. R is not formed there, and no crossing is placed there. Elsewhere on the branch R is used as defined.
- **Remainder:** K_rem = K_exact − K_lead, using the exact identity

  Δ_μν − η_μν/M² = (k_μk_ν − η_μν k_ρk^ρ)/(M²(k_ρk^ρ + M²)).

  This includes both the k_μk_ν part and the −η_μν k² part. L_light is common to both matrices, so K_rem contains only exchange terms.
- **Separate diagnostic, not a criterion:** the pure longitudinal projection K_long = ½(k_μδJ^μ)²/(M²(k_ρk^ρ + M²)). It is reported alongside, labelled as a diagnostic of the linear current exchange only.
- **Ratio:** R(k; Φ̄, branch) = ‖K_rem‖₂/‖K_lead‖₂ (spectral norms). It is evaluated on each light branch ω = c_i k, where the c_i are the roots of the leading-order (B4, Λ = Λ_eff) characteristic polynomial F(c) on Φ̄. Only real c_i are used; a complex c_i is reported, and no ratio is formed for it.
- **k_der(Φ̄)** is the smallest k > 0 with max over branches of R = 1. It is the scale at which the leading local term stops dominating the exchange.
  - If R < 1 for all k > 0, then k_der = ∞ ("no crossing, below").
  - If R ≥ 1 for all k > 0, including k → 0, then k_der = 0 ("no crossing, above").
  - A singular point (k_ρk^ρ + M² = 0 on a branch) is reported, and no crossing is placed there.
- **Name and meaning.** k_der is **not** a strong-coupling or unitarity scale. A true strong-coupling scale k_sc would need an amplitude or perturbativity criterion. That is not defined here; it is B4.12 Q3 territory, and it is not derived in this review.
- **H3 (to check):** k_der with γ_u ≠ 0 is not smaller than k_der of the same background with γ_u = 0, i.e. T7a alone. D3 derives k_der/M analytically as a function of (γ_u, g_φ, g_ψ, Φ̄) where possible. If it is not possible in closed form, D3 states what remains, and H3 is OPEN.

**D4 — perturbativity in the B4.12 sense.**

- **Registered definitions (common-M specification 3cba04e):**
  - ρ_j = |u_j|/M for u_j ∈ {q, w, r, v}, with one common M;
  - couplings g_φ = M√p and |g_ψ|;
  - G = max(|q|, |w|, |r|, |v|);
  - c_max = max(√p, |κ|/√p);
  - criterion c_max·G ≤ g_max·ρ_max.
- **Extension (fixed):** the coupling γ_u is treated like g_φ and g_ψ, so |γ_u| ≤ g_max. From D1, |γ_u| = M√(Λ_eff/2). This gives c_Λ = √(Λ_eff/2) and c′_max = max(√p, |κ|/√p, √(Λ_eff/2)).
- **G′:** J_u = γ_u f_x² involves only q, which is already in G. So the set of gradient components does not change, and **G′ = G** exactly. This must be verified against the operator content of §1. It holds because no new gradient component enters.
- **H4:** c′_max ≥ c_max, by definition of the maximum. With G′ = G, c′_max·G′ ≥ c_max·G, so the criterion cannot become easier.

  D4 must confirm that the extension uses the same normalisation as the registered definitions:
  - g_j ρ_j = √X_j for the existing components;
  - the analogous identity for γ_u;
  - the same thresholds A, B, C.

  This shows that the counts of B4.12 cannot improve. It does not by itself give the per-background PASS/FAIL results. Those are Stage 2.

**D5 — principal symbol of the coupled system.**

- **System:** the four fields (δf, δg, δA0, δA1) linearised on the background of §1, with the full Lagrangian of §1. This includes the terms in which the background Proca field multiplies second variations of J_n and J_u.
- **Procedure (fixed):**
  1. Derive the modified Gauss constraint and eliminate δA0. Expected form: (M² − ∂_x²)A0 = −∂_xA1_t − J_u, to be derived.
  2. Write the reduced first-order system in the registered variables (a, b, c̃, d, A1, F) = (f_t, f_x, g_t, g_x, A1, F) of T3_FUND_PROCA_SPEC §3(iii), extended by the J_u terms. ∂_xJ_u contains b_x, so J_u can enter the principal part of the A1 equation. This must not be assumed absent.
  3. Linearise with the registered scaling (Ã = MδA1, F̃ = (k̂/M)δF, k̂ = √(M² + k²)). Identify the principal symbol P_∞ as the k → ∞ limit, as in the registered §3(iii).
     - If the J_u terms do not fit the pure k·P_s + k̂·P_k structure, report this as a finding.
     - Then use, and state, the weighting (Douglis–Nirenberg weights) needed to obtain a well-defined principal part.
  4. Hyperbolicity: P_∞ has real eigenvalues and a complete set of eigenvectors (strong hyperbolicity), as registered C2.
  5. Keep apart:
     - hyperbolicity of the reduced first-order system, and
     - physical hyperbolicity of the constrained system.

     The Gauss constraint must be preserved by the reduced evolution. Any eigenvalue or eigenvector of P_∞ that belongs only to constraint-violating directions is identified and reported separately. It is not counted for or against physical hyperbolicity.
- Report whether γ_u changes P_∞. Report its eigenvalues as functions of (γ_u, g_φ, g_ψ, Φ̄). This is the analytic prerequisite for any later C1/C2 check. No T3-fund numerics are done here.

## 3. Outcome classes (fixed before the derivation; mutually exclusive)

**Applicable backgrounds:** the 41 registered E3 backgrounds. A condition "on all applicable backgrounds" is established in one of two ways:
- analytically at Stage 1, for a parameter region that contains all 41;
- or by the Stage 2 arithmetic of the Stage 1 formula (§4).

Until it is established, it counts as undetermined.

**Decision rule, applied in this order; the first matching class is the outcome.**

**1. NOT-VIABLE** if at least one of the following holds:

- (a) D0 or D1 gives Λ_eff ≤ 0 for every real γ_u;
- (b) D1 shows that J_u changes p, s or κ, or creates a non-B4 operator at leading order, and this cannot be removed within the model of §1;
- (c) k_der = 0 ("no crossing, above") on all applicable backgrounds;
- (d) D5 shows that the current makes the **physical** (constrained) principal symbol non-hyperbolic: complex eigenvalues, or a defective eigenvalue that the analytic structure cannot resolve.

**2. VIABLE-FOR-NUMERICS** if NOT-VIABLE does not apply and **all** of the following hold:

- (V1) D0–D1 confirm H1: Λ_eff > 0, with p, s, κ and the operator set unchanged at leading order, and the D1 consistency check passed;
- (V2) k_der > M on all applicable backgrounds, from an explicit formula or an explicit bound;
- (V3) D5 gives a strongly hyperbolic physical principal symbol;
- (V4) D2 has no "undecided" terms.

**3. OPEN** in every other case, i.e. NOT-VIABLE does not apply and at least one of V1–V4 fails or is undetermined. In particular, OPEN results when:

- 0 < k_der ≤ M on at least one applicable background;
- k_der = 0 on some, but not all, applicable backgrounds, and k_der > M on the rest. (NOT-VIABLE (c) needs zero on **all** backgrounds, and V2 needs > M on **all**, so neither applies.)
- k_der is undetermined on some applicable background, including any background marked R-UNDEFINED;
- D2 has "undecided" terms;
- H1, the D1 check or D5 cannot be decided.

| Class | Condition on k_der (with the other conditions as above) |
|---|---|
| NOT-VIABLE | k_der = 0 on all applicable backgrounds |
| VIABLE-FOR-NUMERICS | k_der > M on all applicable backgrounds |
| OPEN | 0 < k_der ≤ M on at least one applicable background; **or** k_der = 0 on some but not all applicable backgrounds (with k_der > M on the rest); **or** k_der undetermined (including R-UNDEFINED) on at least one applicable background |

The comparison with M is a criterion for the range of the local approximation. It is not a claim that M is a strong-coupling scale.

**VIABLE-FOR-NUMERICS** only allows a separate specification for numerical checks: T3-fund C1–C3 and T5. It is **not** acceptance of a UV candidate.

**Reported separately, and not part of the class:**

- **H3.** Not lowering k_der relative to T7a alone is reported either way.
- **H4/D4.** If H4 holds, B4.12 cannot be closed positive through this mechanism with better counts than T7a. At best, Proca with a timelike current inherits **B4.12-OPEN**.

## 4. Stages and permissions

| Stage | Content | Status |
|---|---|---|
| 1 | D0–D5 analytically; symbolic checks in one SymPy script (pure algebra, no data); written report | after registration of this specification |
| 2 | Arithmetic only, on the 41 registered backgrounds with their registered Λ, p, s, κ and gradients: (i) the D4 condition; (ii) the D3 formula or bound for k_der from Stage 1, where Stage 1 did not settle it analytically. No new derivation and no new criterion | separate authorisation after the Stage 1 report |
| 3 | Numerical T3-fund / T5, only if the class is VIABLE-FOR-NUMERICS | separate specification |

## 5. What this specification does not do (unchanged)

- It changes no B4.12 threshold (A, B, C) and no T7a result.
- It does not claim that Λ > 0 is realised physically.
- It does not address loops, contact terms, or other UV classes.
- It does not use EXT-1, BEC data or GW data.
- Its results are not evidence for B4.

## 6. Reference correction (unchanged; approved)

B4.12 is not blocked by the definition of ρ. The source was found on 2026-11-15: `evidence/B4_10_proca_mass_sweep_v0_5/B4_12_REPORT.md` v0.5, with ρ_j = |component|/M and one common M. It was closed at document level then. B4.12 remains OPEN for Q2/Q3.

## 7. Requested from the review

1. Confirmation of the R-UNDEFINED rule in D3.
2. Confirmation that the mixed case (k_der = 0 on some backgrounds, > M on the rest) is now explicit under OPEN.
3. Approval for registration of v0.5. (The Stage 2 extension was approved in the review of v0.4.)
