# UV candidate: scalar mediator with a dynamical rapidity η — specification (v0.4.1)

**Changes from v0.1 (review of 2026-11-16, REVISION REQUIRED).**
1. The bare Λ term depends on η through n. The background equation for η̄, the existence region and the background counts are therefore derived anew (Parts C0, E, F4). The U1–U4 counts 17/41 and 24/41 are historical only.
2. The functional integral over σ is Gaussian but nonlocal. Its local form J_tot²/(2M²) is the leading term only (Part F1).
3. The effective mass of η includes the curvature of the Λ term. The mixed (δσ, δη) system is treated separately (Parts E5, F2).
4. F2 and F3 are aligned with these corrections, and the η-elimination condition is a question to be checked, not a presumed negative result.
5. A separate comparison document for dynamic η and dynamic u^μ is required by the reviewer. It is not part of this specification.

**Changes from v0.4 (review of 2026-11-16: APPROVED WITH ONE EDITORIAL CLARIFICATION).** Only the explanatory note on the historical wording in Part A is replaced by the reviewer's text. No other change.

**Changes from v0.3 (review of 2026-11-16, REVISION REQUIRED).** Only T6 items 3 and 4 are changed, plus a clarifying note on the historical wording in Part A. B1, T3-fund, T5, F1 and F2 are unchanged.
- T6.3: in 1+1D a smooth nonvanishing timelike field is always locally hypersurface-orthogonal (the Frobenius condition is identically satisfied). The item now separates local integrability from the global existence of a foliation.
- T6.4: for every finite real η, u·u = 1. The limit η̄ → ±∞ is an unbounded rapidity whose direction approaches the light cone. It is not a background with u·u = 0.

**Changes from v0.2 (review of 2026-11-16, REVISION REQUIRED).**
1. E5: the stability of η is decided on the full quadratic fluctuation system (δφ, δψ, δσ, δη), with all mixing blocks identified. A reduced (δσ, δη) block may be used only with an analytic proof that the reduction preserves the relevant spectrum.
2. D: the statement that k → 0 is the worst case of T5-static is removed as an assumption. It must be proved analytically, or positivity must be checked over the whole k range with analytic control of the asymptotics.

The registered T5 criteria are not changed. Only the way their fulfilment is established is made explicit.

**Status: PROPOSAL FOR REVIEW.** Specification only. No code, no run and no numerical result. No registered criterion, tolerance or background set is changed. B4.12 remains OPEN. T7a remains blocked. T7b (Λ-emergence) is not started. No Zenodo action.

**Guiding criterion (reviewer, 2026-11-16).** The task is not to show that a model with a dynamical η can be built. It is to test whether exactly this construction can generate B4 without losing:
- the independence of the operators;
- stability;
- the physical interpretation of u^μ.

## Registered framework this specification follows

| Item | Source |
|---|---|
| T1–T8 definitions | `docs/UV_ACCEPTANCE_PROTOCOL.md` (T1–T6 fixed before the scalar-mediator runs; T7–T8 added 2026-11-11) |
| Working interpretation: B4 is a Lorentz-violating EFT with a preferred observer field u^μ, c_B4 = c, Λ > 0 kept; preferred-foliation causality; T3 split into T3-B4 and T3-fund | `docs/UV_ACCEPTANCE_PROTOCOL_ADDENDUM_T3_T6.md`, `a3ce2a1` |
| T3-fund protocol (C0–C4, P1–P8, outcomes) | same addendum, registered `220bf76`; clarification 1 `d70c1e5` |
| Regime of this candidate: FUNDAMENTAL | regime table of the addendum (decision of 2026-11-13) |
| T6 requirements (items 1–5) | same addendum, "T6 as a stand-alone task" |
| Existing results (U1–U4) | `docs/SCALAR_MEDIATOR_RESULT.md`, `experiments/uv_candidates/scalar_mediator_dynamicU_check.py`, `docs/UV_CANDIDATE_ASSESSMENT.md` |

Under the registered convention a T3-B4 verdict is not a rejection. T3-B4 values are kept as historical results and are not re-run here.

## Part A — Existing baseline (U1–U4), stated without re-running

**Action used in U1–U4** (signature (+,−), 1+1D, M = 1 in the script):

    L = ½(∂φ)² + ½(∂ψ)² + ½(∂σ)² − ½M²σ² + ½f(∂η)² + σ Σ_i g_i (u·∂φ_i)(n·∂φ_i),
    u = (cosh η, sinh η),   n = (sinh η, cosh η),   φ₁ = φ, φ₂ = ψ.

The couplings are g₁ = √p and g₂ = −κ/√p, the rank-one assignment, as for Proca. η has **no potential**. The **Λ term is not included**.

**Results (first pass, as recorded):**
- **Stationary background.** σ̄ = J(η̄)/M², with J(η) = A cosh 2η + B sinh 2η. Stationarity J′(η̄) = 0 gives tanh 2η̄ = −B/A and requires |B| < |A|. Such an η̄ exists on **17 of 41** manifest backgrounds.
- σ̄ = 0 together with J′ = 0 requires A² = B², a degenerate, lightlike η̄. There is no generic solution. (Wording as recorded. For A² = B² with A ≠ 0, J′(η) = 0 requires tanh 2η = ±1, which has no solution at finite real η. This is a degenerate limiting configuration: the formal stationary solution lies at |η̄| → ∞, where the direction of u^μ approaches the light cone. It is not a finite-rapidity background and does not define a unit lightlike vector.)
- The principal symbol of the φ, ψ blocks in the u frame equals that of the scalar mediator with a_i = g_iσ̄. f and the η sector enter only at lower order.
- **T3-B4 (historical):** 0/17 with front ≤ 1; min 1.007, median 3.07, max 5.85; identical for f = 0.1, 1, 10. All 17 spectra real.
  - Method: a degree-8 fit of det K at X = k/M = 200, approximate.
  - Under the registered convention this is a historical T3-B4 value, not a rejection.

**What U1–U4 did not establish:**
- T1/T7 (the Λ term was absent);
- C0 on the 24 backgrounds without stationary η̄;
- C1 and C2 at the evidence level of T3-fund (the degree-8 fit is not an exact polynomial; P1 and P3 were not applied);
- C3 by P6, and C4;
- T5;
- T6.

Repeating U1–U4 is not authorised without a stated reason. The task builds on them.

## Part B — Action definition (fixed before any computation)

**B1. Candidate action.** It is the U1–U4 action plus the bare Λ term written covariantly in the frame of u:

    L_η = ½(∂φ)² + ½(∂ψ)² + ½(∂σ)² − ½M²σ² + ½f(∂η)² + σ Σ_i g_i (u·∂φ_i)(n·∂φ_i) − (Λ/4)(n·∂φ)⁴.

- The products are formed as in the U1–U4 script: u·∂φ = cosh η φ_t + sinh η φ_x and n·∂φ = sinh η φ_t + cosh η φ_x. At η = 0, (u·∂φ)(n·∂φ) = φ_tφ_x and (n·∂φ)⁴ = φ_x⁴ (lab-frame B4 form).
- **Signature.** This specification uses (+,−), as in U1–U4. B4.12 uses (−,+). Every result must state its signature, and the comparison with B4 must convert explicitly.

**B1a. Background functions (notation used below).** On a background with constant gradients (w, q) for φ and (v, r) for ψ:

    J(η) = A cosh 2η + B sinh 2η,   A = g₁wq + g₂vr,   B = ½[g₁(w² + q²) + g₂(v² + r²)],
    N(η) = sinh η w + cosh η q   (the value of n·∂φ),

with N′(η) = cosh η w + sinh η q. The constant-field Lagrangian is then σJ(η) − ½M²σ² − (Λ/4)N(η)⁴.

**B2. What is fixed and what is not.**
- **Fixed:**
  - the field content (φ, ψ, σ, η);
  - the rank-one couplings g₁ = √p, g₂ = −κ/√p;
  - Λ > 0 as a bare coupling;
  - the absence of an η potential;
  - f > 0.
- **Not added:**
  - any η potential, mass term or further coupling;
  - any second heavy field.

  A variant with such terms is a different candidate and needs its own specification.
- **Λ is inserted by hand.** This candidate therefore does **not** test whether −(Λ/4)φ_x⁴ emerges dynamically. That question is T7b, a separate task, not started.

**B3. Parameter values.**
- f ∈ {0.1, 1, 10}, as in U1–U4.
- M = 1.
- Λ from each background record.
- p, s, κ from each background record. Only (p, κ) enter the rank-one couplings. For general B4, s ≠ κ²/p, and this is recorded (Part F).

## Part C — T3-fund (registered protocol, FUNDAMENTAL regime)

**Regime.** FUNDAMENTAL, by the registered regime table. It is not changed.

**C0 — Background.** The background is constant φ, ψ gradients (w, q, v, r), constant σ̄ and constant η̄. With the action B1:
- the φ and ψ equations hold identically for constant fields (all terms are functions of constant first derivatives);
- the σ equation gives σ̄ = J(η̄)/M²;
- **the η equation** (constant fields, signature (+,−)) is

      σ̄ J′(η̄) − Λ N(η̄)³ N′(η̄) = 0,

  equivalently, after inserting σ̄, d/dη [ J(η)²/(2M²) − (Λ/4) N(η)⁴ ] = 0 at η = η̄.

Because of the Λ term, tanh 2η̄ = −B/A is **no longer** the stationarity condition, and |B| < |A| no longer defines the existence region by itself. To be derived before any root test:
- (a) the number of real solutions η̄ of the η equation on each background, including the possibility of none or of several;
- (b) for several solutions, which one is used, by a rule fixed here: every real solution is reported, each as a separate background, and none is discarded after the results;
- (c) the derivation that each solution is an exact solution of the full equations (analytic, or residual ≤ 10⁻¹⁰ by C0).

**Background counts.** The U1–U4 counts (17 with and 24 without a stationary η̄) are historical. They are **not** carried over. The counts for B1 are part of the result. Backgrounds without any real η̄ are INCONCLUSIVE for T3-fund (C0 not satisfied), are recorded under T6 item 4, and are not dropped. A time-dependent η background is not constructed in this task.

**C1–C3.**
- The exact characteristic polynomial det K(ω, k) is derived symbolically, and its exact factors are divided out (P1). No polynomial fit is used.
- The polynomial is evaluated on the P2 grid: 61 points in [10⁻³, 10³] plus 31 points in [10³, 10⁶].
- The roots are classified by P3, with 50-digit confirmation.
- The k = 0 sector is treated by P2a, including the ZERO MODE class.

**C1 — asymptotic statement.** Because the regime is FUNDAMENTAL, C1 PASS requires an analytic large-k argument, e.g. a strongly hyperbolic principal symbol for the full (φ, ψ, σ, η) system. Without it, C1 is INCONCLUSIVE ("numerically passed on [k_lo, k_hi]").

**C2 — principal symbol.** Documentation by P5:
- the first-order form;
- the normalisation of the variables;
- the reality of the principal roots;
- the condition number of the eigenvector matrix;
- for coinciding principal roots, an analytic argument on completeness.

The η block has the principal part f(ω² − k²). Its degeneracy with the luminal roots of the other fields must be addressed explicitly.

**C3.** By P6, with u(k) measured at 50 and 100 digits and the three-band rule.

**C4 — time domain.**
- **Background:** one designated background with stationary η̄.
- **Procedure:** P7 in full:
  - a Gaussian packet in the light fields, amplitude 10⁻³;
  - a norm criterion ≤ 3 over ≥ 1000/M;
  - an energy-drift target of 10⁻⁹, calibrated on the free system with the background values of σ and η re-chosen to solve the free equations;
  - convergence at (N, dt) and (2N, dt/2).

**Outcomes per background:** PASS, FAIL-ILLPOSED, FAIL-UNSTABLE or INCONCLUSIVE, as registered. They are reported for each value of f.

**Backgrounds (P8):** the B3 anchor plus the 40 records of `data/e3_manifest.json` (41 in total), with their recorded seeds.

## Part D — T5 (energy)

The T5 definitions of `UV_ACCEPTANCE_PROTOCOL.md` are used, with the heavy-scalar procedure of `HEAVY_SCALAR_POPULATION_REPORT.md` as template:
- **T5-static:** the quadratic energy form E₂(k) of the full fluctuation system (δφ, δψ, δσ, δη) about the stationary background must be positive for all k. Fulfilment is established in one of two ways, recorded with the result:
  1. **analytic:** a proof that the minimum eigenvalue of E₂(k) over k is attained at k → 0, after which the Schur-complement criterion at k → 0 decides. It is **not** assumed that k → 0 is the worst case: the nonlinear Λ term and the mixing of the fields can make a finite k critical. Or
  2. **over the whole range:** positivity checked on the P2 grid (k/M in [10⁻³, 10⁶]), together with an analytic large-k statement for E₂(k) and the separate analytic limit k → 0. A grid alone, without the asymptotic statement, gives INCONCLUSIVE for T5-static.
- **T5-Hamiltonian:** the number of negative eigenvalues of the Hermitian H₂(k) in (q, π) at k = 0.01, 0.1, 1 and 10. Tolerance −10⁻¹² for the minimum eigenvalue, as registered for the scalar mediator.
- **T5-Krein:** the sign of the mode energy for each real-frequency normal mode, reported separately.
- **Frame of the energy.** The energy is stated in the frame of u (preferred time), and also in the lab frame where the two differ. T5 is sufficient but not necessary for C3.

## Part E — T6 (origin, domain and stability of u^μ)

The five registered T6 items are answered for this construction:

1. **What u^μ is physically.** In this construction u is the timelike unit vector with rapidity η, a dynamical field with kinetic term ½f(∂η)² and no potential of its own. To be stated: what selects η̄. By the η equation of C0, η̄ is a function of the light-field gradients (w, q, v, r) and of Λ. That is, the frame is set by the matter backgrounds and by the Λ term, not externally. This is to be stated as a property of the construction, not as an interpretation of nature.
2. **How u is determined.** Dynamically, through the η equation f□η = σ ∂_ηJ − Λ(n·∂φ)³ ∂_η(n·∂φ) (signature and sign convention as in B1). Conditions to be derived from the full equation of C0, not from J′ = 0 alone:
   - existence of η̄;
   - its uniqueness;
   - whether it is a minimum or a maximum of the effective potential V_eff(η) = −[J(η)²/(2M²) − (Λ/4)N(η)⁴].
3. **Foliation.** In 1+1D the Frobenius condition for hypersurface orthogonality, which involves a totally antisymmetrised triple of indices, vanishes identically. A smooth, nonvanishing timelike field is therefore always locally hypersurface-orthogonal, for constant and for non-constant η alike. Local integrability is not a question here, and nothing "fails locally" for non-constant η. Two things are kept apart:
   - **local integrability:** automatic in 1+1D, stated as such;
   - **global existence of a foliation by the preferred time:**
     - for constant η̄ it is immediate: the leaves are the straight lines orthogonal to u;
     - for non-constant η(t, x) it can require additional conditions on the domain, the regularity of η and its behaviour at large |η| or at the boundary of the domain of Part E4. These conditions are to be stated, not assumed. Time-dependent η backgrounds are outside this task (Part C0), so only the conditions are recorded here.
4. **Where u is not defined.** The domain of definition of η̄ in background space is the set where the η equation of C0 has a real solution; it is to be derived and stated for B1. The limit η̄ → ±∞ is to be identified for B1. For every finite real η, u·u = cosh²η − sinh²η = 1, so u stays a unit timelike vector. As |η̄| → ∞ its components grow without bound and its normalised direction approaches a null direction. There is **no** finite lightlike limit of u as a unit vector, and this limit is **not** a background with u·u = 0. It is reported as unbounded rapidity at the edge of the domain, with the direction approaching the light cone. Outside the domain no stationary frame exists (the U1–U4 value 24/41 is historical and does not apply to B1). This is the analogue of the region where static observers do not exist.
5. **Stability and well-posedness of the η sector.**
   - **Low-energy effective mass** of η about η̄, after the leading local elimination of σ (F1), for signature (+,−) and the kinetic normalisation ½f(∂η)²:

         m²_η,eff = −(1/f) d²/dη² [ J(η)²/(2M²) − (Λ/4) N(η)⁴ ] at η = η̄.

     The Λ term contributes to the curvature, so the second derivative of σ̄J alone is not sufficient. The sign is to be derived before η is classified as light, heavy or unstable.
   - **Full fluctuation system.** The stability of η is decided on the complete quadratic fluctuation system (δφ, δψ, δσ, δη) about the background, before any elimination. To be done:
     - construct the full quadratic action;
     - identify every mixing block. On backgrounds with nonzero gradients the coupling σ Σ_i g_i (u·∂φ_i)(n·∂φ_i) mixes δσ and δη with δφ and δψ, and the Λ term mixes δη with δφ;
     - compute the physical eigenfrequencies of the full system, at k = 0 (by P2a) and as functions of k (by the T3-fund grid).
   - **Reduced block.** A reduced (δσ, δη) block may be used only together with an analytic proof that the reduction preserves the relevant spectrum, e.g. that the mixing blocks vanish on the background considered, or are of a stated order that cannot change the sign of any eigenvalue. Without such a proof the reduced block is a diagnostic and does not decide stability. Otherwise η could be declared stable while an instability sits in a mode mixed with the light fields.
   - m²_η,eff is a low-energy quantity. It does not replace the full system.
   - Whether η is light (|m_η,eff| ≪ M), heavy, or tachyonic (m²_η,eff < 0).
   - The η-sector contributions to C1–C3 and T5.

**T6 outcome classes:**
- **STATED-DYNAMICAL:** items 1–5 answered, η̄ exists, m²_η,eff ≥ 0 (confirmed in the full fluctuation system (δφ, δψ, δσ, δη)), and the η sector passes C1–C3 on the backgrounds with η̄.
- **PARTIAL:** η̄ exists on a stated subset only. The subset and its complement are both reported.
- **FAIL:** m²_η,eff < 0 at η̄ with an instability confirmed by C3 in the full fluctuation system, or the η sector fails C1–C3.
- **INCONCLUSIVE:** an item cannot be answered from the action.

## Part F — T1/T7 (matching to B4)

**Classical elimination and functional integration are reported separately.**

**F1. Elimination of σ.**
- **Classical elimination at tree level:** σ = J_tot/M² + O(∂²/M⁴), where J_tot = Σ g_i (u·∂φ_i)(n·∂φ_i). This gives L_eff ⊃ J_tot²/(2M²).
- **Functional integration:** σ enters quadratically, with kinetic term ½(∂σ)² and constant mass, and linearly through J_tot. The Gaussian integral is exact, but its result is **nonlocal**:

      S_eff ⊃ ½ ∫d²x J_tot (□ + M²)⁻¹ J_tot,

  with the operator sign fixed by the signature (+,−) and the boundary conditions, up to a field-independent constant. It is to be stated with its sign.
- **Low-energy expansion:** (□ + M²)⁻¹ = M⁻² − □/M⁴ + …. The leading local term is J_tot²/(2M²), and the corrections are O(M⁻⁴).
- **Consequence:** classical elimination and functional integration agree at leading order. It is **not** claimed that the exact functional integral gives only the local term. The O(M⁻⁴) terms are reported separately as remainders, as in v0.4 §4 item 2.

**F2. Status of η: the condition for its elimination.** η has no potential and its kinetic coefficient is f. To be derived:
- (i) The low-energy effective mass m²_η,eff of Part E5, which includes the Λ-term curvature, compared with M². The eigenfrequencies of the full fluctuation system (δφ, δψ, δσ, δη) of Part E5 are used to confirm the classification.
- (ii) **If η is light** (|m_η,eff| ≪ M), it **cannot be integrated out**. The low-energy theory then has three light fields (φ, ψ, η), not B4. This would be a negative T7 result for this construction at the stated order. It is not to be avoided by a choice of f after the result. Whether η is light is a question to be checked, not a presumed result. The Λ-term curvature may change the answer relative to U1–U4.
- (ii′) **If m²_η,eff < 0**, η̄ is a maximum of V_eff, and the elimination is not valid at that background. This is recorded under T6 and C3.
- (iii) **The limit f → 0** makes η an auxiliary field with the local constraint d/dη [J²/(2M²) − (Λ/4)N⁴] = 0, the field-dependent version of the C0 equation. Classically u is then fixed pointwise by the light-field gradients. To be derived:
  - whether this limit is regular;
  - which effective operator results;
  - whether the coefficients become functions of the background, i.e. coefficient functions rather than constants.
- (iv) Classical elimination of η (saddle point) and functional integration over η (not Gaussian: η enters through cosh and sinh) are not equivalent beyond tree level. Only tree level is in scope; the loop difference is recorded as not computed.

**F3. Structure to be checked (reviewer's list):**
1. **Sign and normalisation of −(Λ/4)φ_x⁴.** In this candidate the term is bare (B1). The check is that, at η̄ = 0, the term enters the lab-frame quadratic action as 1 + 3Λq² in d1, as in B4. For η̄ ≠ 0 the term is in the u frame.
2. **Independence of the coefficients.** With rank-one couplings the σ exchange gives the rank-one subfamily s = κ²/p only. This is known by construction, as for Proca and the scalar mediator. General B4 with independent (p, s, κ) is **not reachable** by this construction. This is declared here, before any computation, so that it is not reported later as a finding or hidden.
3. **Higher-order operators.** The M⁻⁴ remainders and all operators generated by η (Part F2) are listed with their order.
4. **Conditions for eliminating η:** Part F2, with the result stated as one of (ii), (iii) or not eliminable.
5. **Compatibility with T3-fund, T5 and T6:** reported per background set, not merged into one verdict.

**F4. The frame issue for matching.** B4 is defined in a fixed frame with coefficients (p, s, κ, Λ). In this construction the coupling structure is B4-like in the u frame. It agrees with the lab-frame B4 form only where η̄ = 0. By the C0 equation at η = 0 (J(0) = A, J′(0) = 2B, N(0) = q, N′(0) = w), η̄ = 0 is a solution if and only if

    2AB/M² − Λq³w = 0.

(In U1–U4, without the Λ term, this was B = 0.) To be stated:
- on which backgrounds η̄ = 0 (or |η̄| small) holds;
- what the comparison means where it does not hold, since then the u frame and the lab frame differ.

The comparison is made in the u frame. Any identification of the u frame with the lab frame of B4 must be justified, not assumed.

**T7 outcome classes:**
- **REPRODUCED-RANK-ONE:** at tree level, after the elimination established in F2, the effective operators equal B4 restricted to s = κ²/p in the u frame, with the Λ term of the correct sign and normalisation. Remainders are of stated higher order, and the field normalisation is identified.
- **NOT-REPRODUCED:** η cannot be eliminated (F2(ii)), or extra operators appear at the same order as the B4 operators.
- **CONDITIONAL:** reproduced only in the f → 0 limit or only on a stated background subset. The condition is part of the result.
- **INCONCLUSIVE:** a step of F1–F3 cannot be completed.

General B4 (independent s) is outside every class above by F3 item 2.

## Part G — Outcome rules (fixed before any computation)

- Each of T3-fund, T5, T6 and T7 is reported separately, with its own classes (Parts C–F), per background set and per value of f. They are not merged into a single verdict.
- **The candidate is not called accepted** unless T3-fund (registered form) and T6 are evaluated, by the registered acceptance rule. Even then the most it can be is a **rank-one** realisation with a **bare** Λ.
- **Not grounds for PASS:**
  - generating some nonlinear operator;
  - a real spectrum alone (T2);
  - agreement of spectral roots without operator identification (T7).
- **Not grounds for FAIL:**
  - a T3-B4 front above 1 (registered convention);
  - unbounded propagation speed (T3-fund does not bound speeds).
- **Backgrounds without η̄:** INCONCLUSIVE for T3-fund (C0), and recorded under T6 item 4. They are not dropped from the counts.
- **No parameter (f, M, couplings, background set, grid, tolerance) is changed after results are seen.**

## Out of scope

- T7b (dynamical emergence of −(Λ/4)φ_x⁴): a separate task, not started.
- T7a and the Proca B4.12 track.
- Dynamic u^μ with a different action (e.g. Einstein-aether type with a unit-norm constraint and further kinetic terms). Not merged with this candidate. A comparison of the two actions is a separate step, required before any merging.
- An η potential, an η mass term, or a second heavy field: different candidates.
- Time-dependent η backgrounds.
- Loop effects; positivity bounds (recalled from the literature, not re-derived).
- T4 (perturbativity) and T8 (observable). T8 is noted as open; it is not evaluated here.
- Any change of a registered protocol, criterion, tolerance or background set.

## Deliverables after approval (each with its own review)

1. A script specification: implementation details, two independent routes where the protocol requires them, and file names.
2. The script, then static review, then authorisation, then one run, then review of the results, then registration.

## Requested from the reviewer

1. Approval of Parts A–G as the specification of this candidate.
2. Confirmation of B1, the Λ term as −(Λ/4)(n·∂φ)⁴ (bare, in the u frame), and of the statement that T7b is not tested here.
3. Confirmation that the declared rank-one limitation (F3 item 2) and the η-elimination condition (F2) are the intended T7 framing.
4. Confirmation of the corrected η equation (C0) and of m²_η,eff (E5). The comparison document `UV_DYNAMIC_ETA_VS_DYNAMIC_U_COMPARISON.md`, required by the reviewer, is prepared separately.
