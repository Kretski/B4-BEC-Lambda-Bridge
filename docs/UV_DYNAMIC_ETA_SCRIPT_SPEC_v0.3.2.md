# UV candidate dynamic η — script specification (v0.3.2)

**Status: PROPOSAL FOR REVIEW.** This document specifies the implementation of the registered candidate specification. It contains no code and no result. No run is authorised. Order: this specification, then review, then the script, then static review, then explicit authorisation, then one run, then review of the results, then registration.

**Changes from v0.3.1 (review of 2026-11-16, one last correction).** Only the DC1 row of the table in §4 is changed. It now refers to the single comparison rule of S2.3, instead of stating a separate relative bound.

**Changes from v0.3 (review of 2026-11-16, approved with one correction).** Only the ZERO (agree) rule of S2.3 is changed. It now also requires the comparison rule, so two values below atol with opposite signs and a difference above the tolerance are classified disagree. Tolerances, method, steps and points are unchanged.

**Changes from v0.2 (review of 2026-11-16, REVISION REQUIRED).** Only two items are changed.
1. S6.2: mutually exclusive classes TACHYONIC (m² < 0), LIGHT (0 ≤ m² < 10⁻²M²) and HEAVY (m² ≥ 10⁻²M²), with a BOUNDARY flag.
2. S2.3: the claim of a guaranteed error bound is removed. A combined absolute and relative tolerance relative to the matrix scale is used, with explicit handling of zero entries and full reporting. Steps, precision and points are unchanged.

**Changes from v0.1 (review of 2026-11-16, REVISION REQUIRED).**
1. New section §1a, "Evidence scope and experimental validation", adopted verbatim from the review.
2. S1.2: the root counting is stated as exact for the rationalised decimal input values, not for the unknown exact physical parameters.
3. S2.3: route R2 now has a procedure, fixed in advance, for the precision, the finite-difference scheme, the steps and a convergence check. Non-convergence of R2 is reported separately and is not a FAIL of route R1.

**Registered basis.**
- `docs/UV_DYNAMIC_ETA_SPEC_v0.4.1.md` (`5fdd4f6`, SHA-256 1379e599…): the candidate, Parts A–G.
- `docs/UV_DYNAMIC_ETA_VS_DYNAMIC_U_COMPARISON_v0.3.1.md` (`5fdd4f6`, SHA-256 024aaae1…): the Einstein-aether generalisation.
- `docs/UV_ACCEPTANCE_PROTOCOL.md` and `docs/UV_ACCEPTANCE_PROTOCOL_ADDENDUM_T3_T6.md`: T1–T8, the T3-fund protocol (`220bf76`, `d70c1e5`) and T6.

Nothing in them is changed here. Where this document and the candidate specification differ, the candidate specification governs, and the difference is a defect of this document.

**Requirements of the review of 2026-11-16, carried into this specification:**
1. The choice between V1 and V2, with an independent comparison of both (S0.2).
2. A symbolic check of the reduction to (c₁₄, c₁₂₃) and of the variational equations (S0.1, S0.2).
3. For the Einstein-aether variant, the principal symbol of the full coupled system, not only of the isolated δη block (S3.2, S8).
4. m²_η,eff is normalised by c₁₄. The full dispersion of δη also contains c₁₂₃, and m²_η,eff is not presented as the complete elimination condition at arbitrary k (S6.3, S7.2).

## 1. Two scripts, not one

The work is split into two scripts with separate reviews, to keep each one reviewable:

| Script | Content | This document |
|---|---|---|
| `uv_dynamic_eta_check.py` (script A) | symbolic and spectral parts: S0–S7 below (C0–C3 of T3-fund, T5, T6, T7) | specified here |
| time-domain script (script B) | C4 of T3-fund (P7: wave packets, norm, calibrated energy drift, convergence) | **not specified here**; separate specification after script A |

Until script B has run, C4 is reported as NOT EVALUATED. By the registered outcome rules, "passes T3-fund on a set" therefore cannot be stated after script A alone.

## 1a. Evidence scope and experimental validation

(Text of the review of 2026-11-16, adopted verbatim as part of this specification.)

### Purpose and limits of Script A

Script A tests the internal theoretical consistency of the dynamic-η construction: background existence, variational equivalence, spectral properties, stability, hyperbolicity, preferred-foliation causality, and recovery of the rank-one B4 operator structure.

These checks are necessary theoretical tests. They are **not experimental validation of B4** and must not be described as evidence that B4 has been observed in nature.

Script B extends the analysis to time-domain evolution. Its successful completion would provide a further numerical consistency check, not experimental confirmation.

### Independent empirical validation programme

An independent empirical-validation line is required in parallel with the dynamic-η theoretical programme. It must use published experimental measurements, preserve the distinction between measured quantities and model-derived quantities, and report uncertainties and competing explanations.

**EXT-1 — Two-component BEC validation.**

The immediate empirical target is the two-component BEC model in the repository, compared quantitatively with published measurements of the density-mode and spin-mode sound velocities.

The analysis must:

- identify the exact experimental dataset, measured observables, uncertainties and calibration assumptions;
- derive predictions from the model before comparing them with the measurements;
- compare both density and spin branches, not only a selected branch;
- report residuals, uncertainty-aware goodness of fit and comparison with the relevant standard BEC description;
- distinguish fitted parameters from independently fixed parameters;
- document all exclusions, fitting choices and sensitivity tests.

A successful EXT-1 result would constitute an empirical test of the BEC-side model and its predictive adequacy for that physical platform. It would **not, by itself, demonstrate B4 or establish that the BEC system is a physical realisation of B4**. Any proposed mapping between the BEC model and B4 must be stated separately and tested independently.

### Other possible empirical line — GWOSC

An independent analysis of GWOSC gravitational-wave data may be conducted as a separate validation module. It must use the strain data and documented analysis conventions, rather than merely converting a published GWTC-4.0 bound.

Its result should be reported as an independent constraint, including the uncertainty and sensitivity to waveform and phase conventions. A constraint compatible with the B4 parameter range is not a detection or confirmation of B4.

### Evidence classification and reporting

All results must be labelled by evidence level:

1. **Theoretical consistency:** analytical derivations and internal mathematical checks.
2. **Numerical consistency:** reproducibility, symbolic/numerical cross-checks and time-domain simulations.
3. **Analogue-system validation:** comparison with independent measurements in a related physical platform, such as EXT-1.
4. **Independent empirical constraint:** analysis of observational data, such as GWOSC.
5. **B4-specific empirical support:** a measured observable for which B4 makes a quantitative, discriminating prediction that is tested against relevant competing models.

No result at levels 1–4 may be described as level 5 without a separate, explicit argument and test.

The dynamic-η programme and the empirical-validation programme are complementary but independent workstreams. Passing the theoretical checks does not substitute for empirical validation. A successful analogue-system comparison does not automatically validate B4.

**Evidence level of this specification.** All outputs of script A are at level 1 (theoretical consistency), together with the cross-checks of §4 (level 2). They are labelled so in the results JSON.

## 2. Inputs, parameters and environment (fixed before any run)

**Backgrounds.**
- The B3 anchor, `data/B3_max_speed_7_315096.json`.
- The 40 records of `data/e3_manifest.json` (key `selected`).
- Total 41, as in U1–U4 and in the candidate specification (P8).
- Each record supplies (p, s, κ, Λ, w, q, v, r). The script records the SHA-256 of both files and stops if p ≤ 0 for any record, because g₁ = √p.

**Couplings.** g₁ = √p and g₂ = −κ/√p (rank one, B2), with M = 1.

**f.** f ∈ {0.1, 1, 10} (B3).

**Einstein-aether variant.** It is implemented but **not run**. The values of (c₁₄, c₁₂₃) are not fixed in any registered document. The default run uses only dynamic η, i.e. c₁₄ = c₁₂₃ = f. A run of the variant requires a separate registered parameter document.

**Environment.** Python 3.13 and SymPy 1.14 as before, plus mpmath, which is bundled with SymPy. Invocation: `python3 -I uv_dynamic_eta_check.py <outdir>`. Inputs are read only from the two data files above. Signature (+,−) throughout. Every reported quantity states its frame: u frame or lab frame.

## 3. Tasks

### S0 — Symbolic checks of the comparison document (prerequisites)

**S0.1 Reduction of the Einstein-aether kinetic term.**
- Build K^{αβ}_{μν} with symbolic c₁…c₄ and u = (cosh η, sinh η).
- Show symbolically that −½K∂u∂u = ½c₁₄(D_uη)² − ½c₁₂₃(D_nη)².
- Show that the c₂ and c₃ terms are pointwise equal for the unit field.
- Show that c₃-term − c₂-term = ∂_α(u^μ∂_μu^α − u^α∂_μu^μ) for a general (non-unit) field.
- Show that c₁₄ = c₁₂₃ = f gives ½f(∂η)².

**S0.2 Variational equations, V1 and V2, two routes.**
- **V1:** the Euler–Lagrange equation of the full action with respect to η, after the substitution u(η), n(η).
- **V2:** n(u) = (u¹, u⁰)/√(u·u), with free variation in u^μ and λ. Then:
  - project the u equation on n and compare it with V1;
  - project it on u and solve for λ.
- Both procedures are implemented, as two independent routes. **Pass condition (DC0):** the tangential projection of V2 equals the V1 equation identically (symbolic difference 0, after `together`). Any difference is reported, not resolved by choice. λ is reported from V2 only.
- Note on SymPy: `euler_equations` drops identically zero equations. Each equation is therefore built explicitly from ∂L/∂(field) − D_α ∂L/∂(∂_α field).

### S1 — Background equation and C0

**S1.1** Derive symbolically the η equation for constant fields, d/dη[J²/(2M²) − (Λ/4)N⁴] = 0. Check that it agrees with σ̄J′ − ΛN³N′ = 0 after σ̄ = J/M².

**S1.2 Complete enumeration of η̄, exact.** With z = e^η > 0:

    J = ½[(A + B)z² + (A − B)z⁻²],   N = ½[(w + q)z + (q − w)z⁻¹].

The η equation is then a Laurent polynomial in z. Multiplied by the lowest needed power of z, it becomes a polynomial P(z) with coefficients in (A, B, w, q, Λ, M). Real η̄ ↔ positive real roots z of P. For each background:
- convert the floating-point record values to exact rationals (`sympy.Rational(str(x))`; the conversion is recorded);
- count the positive real roots exactly (Sturm sequence, `sympy.count_roots` on (0, ∞)). This is **exact root counting for the rationalised decimal input values**. It is not a statement about unknown exact values of the physical parameters, and it does not remove the uncertainty of the input data. Roots that are double, or whose separation is below 10⁻²⁰ relative, are flagged, because an infinitesimal change of the input could change the count;
- locate each root with `sympy.real_roots` / `nroots` at 50 digits;
- record η̄ = ln z, σ̄ = J(η̄)/M², and whether each root is simple.

Multiple roots are reported as separate backgrounds (C0(b)). Backgrounds with no positive root are counted and recorded (INCONCLUSIVE for T3-fund, T6 item 4). The U1–U4 counts are not used.

**S1.3 C0 residual.** For each (background, η̄) substitute into the full field equations (φ, ψ, σ, η, all constant): residual ≤ 10⁻¹⁰ at 50 digits. The φ and ψ equations hold identically for constant fields. This is stated analytically and confirmed numerically.

**S1.4 Lab-frame condition (F4).** For each background, report 2AB/M² − Λq³w and whether η̄ = 0 is one of the roots.

### S2 — Full quadratic fluctuation system

**S2.1** Expand L to second order in (δφ, δψ, δσ, δη) about each background (route R1: symbolic, SymPy). Identify every block of the 4 × 4 kinetic, mixed and potential structure. The required mixing blocks are recorded explicitly (E5):
- (δσ, δφ_i) from σΣg_i(u·∂φ_i)(n·∂φ_i);
- (δη, δφ_i) from the same term and from the Λ term;
- (δσ, δη).

**S2.2** Plane waves e^{i(kx − ωt)} give the matrix K(ω, k). det K(ω, k) is an exact polynomial in (ω, k) with background-dependent coefficients.

**S2.3 Route R2 (independent).** Build K(ω, k) numerically by finite-difference second derivatives of the Lagrangian density with respect to the fluctuation jets (fields and their first derivatives). It must not use R1 expressions. Procedure, fixed before any run:
- **precision:** mpmath with 100 significant digits;
- **scheme:** second-order central differences for each second derivative (four-point stencil for mixed derivatives, three-point for pure ones);
- **steps:** h₁ = 10⁻²⁰ and h₂ = 10⁻²⁵, with each step relative to max(1, |jet value|). For smooth entries of moderate size the truncation error O(h²) and the rounding error O(10⁻¹⁰⁰/h²) are **expected** to be far below the tolerance. This is not guaranteed: the actual error also depends on the scale of the function, the size of the relevant derivatives, cancellation between nearly equal values, and on entries that are zero or nearly zero. The tolerance below is therefore a combined absolute and relative one, and convergence is checked, not assumed;
- **comparison rule (used for both convergence and DC1).** For each point (ω, k) let S = max over all entries of |K_R1(ω, k)|, the scale of the matrix at that point. Two values a and b agree if

      |a − b| ≤ atol + rtol · max(|a|, |b|),   rtol = 10⁻³⁰,   atol = 10⁻³⁰ · S;

- **zero and near-zero entries.** An entry is classified **ZERO (agree)** only if both |a| ≤ atol and |b| ≤ atol **and** a and b satisfy the comparison rule above. For it the relative difference is not used as a criterion. If both are below atol but violate the rule (e.g. a = +0.9 atol, b = −0.9 atol, |a − b| = 1.8 atol), the entry is **disagree**. In every other case the comparison rule alone decides;
- **convergence check:** for every entry, the R2 values at h₁ and h₂ must agree by the comparison rule. Otherwise that entry is **R2 NOT CONVERGED**;
- **reporting:** for every entry, the absolute difference, the relative difference (where defined), the scale S, |a|, |b| and the classification (agree, ZERO (agree), disagree, R2 NOT CONVERGED);
- **points:** 5 fixed (ω, k) points per background, chosen before any run: (ω, k) = (0.3, 0.1), (1.1, 0.7), (2.0, 3.0), (0.5, 10.0), (7.0, 2.0), in units of M.

**Pass condition (DC1):** every converged entry agrees between R1 and R2 (at h₂) by the comparison rule, ZERO (agree) entries included. Any entry marked R2 NOT CONVERGED makes DC1 INCONCLUSIVE, not FAIL. Such entries are reported as a limitation of route R2, not as an error of route R1. Precision, scheme, steps, points, rtol and atol are not changed after results are seen.

### S3 — T3-fund C1–C3 (FUNDAMENTAL regime)

**S3.1** P1: divide exact factors out of det K symbolically before any root test.

**S3.2 Principal symbol, C2 (P5).** Derive the principal symbol of the full system (highest-derivative part), analytically. To be stated:
- which blocks enter it. Expected: the φψ block with the σ̄g_i and 3ΛN̄² contributions, the σ block, and the η block. The mixing terms are of lower differential order. **This expectation is to be derived, not assumed.**
- the first-order form;
- the normalisation of the variables;
- the reality of the principal roots;
- the condition number of the eigenvector matrix (indicator only, threshold 10⁸);
- for coinciding principal roots (relative 10⁻⁶), the analytic argument on a complete set of eigenvectors, without which the case is INCONCLUSIVE.

In the Einstein-aether variant the η block becomes c₁₄ω′² − c₁₂₃k′² in the ū frame. The principal symbol of the **whole** system is derived, not only that block.

**S3.3 C1.** An analytic large-k argument from S3.2 (strong hyperbolicity), if it can be given. Otherwise C1 is "numerically passed on [k_lo, k_hi]" = INCONCLUSIVE.

**S3.4 Grid and classification.**
- P2: k/M on 61 log points in [10⁻³, 10³] and 31 in [10³, 10⁶].
- P3: double-precision classification, r = |Im c|/max(1, |c|), with the bounds 10⁻⁷ and 10⁻⁵. Candidate and unresolved roots are recomputed at 50 digits, and a root is confirmed if r > 10⁻¹².
- P2a: k = 0 by det K(ω, 0) = 0 directly, including the ZERO MODE class.
- P4: fit γ ~ k^s on [10⁴, 10⁶] where confirmed complex roots exist (diagnostic).

**S3.5 C3 (P6).**
- γ(k) = max Im ω.
- u(k) = |γ₅₀ − γ₁₀₀| plus the refinement term (factor-4 refined grid inside any band with γ > u).
- The three-band rule with the margin 10.
- FAIL-UNSTABLE only after confirmation at both precisions and on ≥ 3 consecutive refined points.

**S3.6 Outcomes** per (background, η̄, f): PASS (C0–C3; C4 not evaluated), FAIL-ILLPOSED, FAIL-UNSTABLE, INCONCLUSIVE. Counts per f.

### S4 — T5

Following the candidate specification Part D:
- **S4.1** The Hermitian energy form E₂(k) of the full system in (q, π) variables, in the u frame. Also in the lab frame where η̄ ≠ 0.
- **S4.2 T5-static.** Either an analytic proof that min eig E₂(k) is attained at k → 0, or the full grid of S3.4 together with an analytic large-k statement and the analytic k → 0 limit. A grid alone gives INCONCLUSIVE.
- **S4.3 T5-Hamiltonian:** the number of negative eigenvalues of H₂(k) at k = 0.01, 0.1, 1, 10, with tolerance −10⁻¹².
- **S4.4 T5-Krein:** the sign of v†H₂v for each real-frequency normal mode, at the same k.

### S5 — T6

The five items of the candidate specification Part E:
- **S5.1** The domain of definition: the backgrounds with ≥ 1 positive root of P(z), from S1.2. The boundary behaviour: roots with z → 0 or z → ∞, i.e. unbounded rapidity, are flagged when the root moves beyond |η̄| > 10.
- **S5.2** Uniqueness: the number of roots per background.
- **S5.3** Minimum or maximum of V_eff(η) = −[J²/(2M²) − (Λ/4)N⁴] at each η̄: the sign of V_eff″.
- **S5.4 Foliation.** For constant η̄, the statement that the leaves are straight lines orthogonal to ū (analytic). Local integrability in 1+1D is stated as automatic. Non-constant η is out of scope.
- **S5.5** T6 outcome class per background: STATED-DYNAMICAL, PARTIAL, FAIL or INCONCLUSIVE, as defined in Part E.

### S6 — η mass and the low-energy η sector

- **S6.1** m²_η,eff = −(1/f)·d²/dη²[J²/(2M²) − (Λ/4)N⁴] at η̄, exact at 50 digits, with its sign. In the aether variant the prefactor is −1/c₁₄.
- **S6.2** Classification of η, mutually exclusive, with the threshold 10⁻² M² fixed before any run:
  1. **TACHYONIC:** m²_η,eff < 0;
  2. **LIGHT:** 0 ≤ m²_η,eff < 10⁻² M²;
  3. **HEAVY:** m²_η,eff ≥ 10⁻² M².

  A tachyonic mass term is never called light because its absolute value is small. m²_η,eff is computed at 50 digits. A value within 10⁻⁴⁰ of a class boundary (0 or 10⁻² M²) is flagged as BOUNDARY and is not assigned to a class.
- **S6.3** m²_η,eff is the mass term normalised by the time kinetic coefficient (f, or c₁₄ in the variant). The low-energy dispersion of the isolated δη block is c₁₄ω′² − c₁₂₃k′² − c₁₄m²_η,eff (dynamic η: c₁₄ = c₁₂₃ = f). The classification uses the eigenfrequencies of the full system (S3) to confirm it. m²_η,eff alone is not reported as the elimination condition.

### S7 — T1/T7

- **S7.1 Elimination of σ.**
  - Classical: σ = J_tot/M² at leading order.
  - Exact Gaussian: the nonlocal ½J_tot(□ + M²)⁻¹J_tot, expanded to O(M⁻⁴).
  - Verified symbolically: the leading local term J_tot²/(2M²), and the O(M⁻⁴) remainder written out and listed separately.
- **S7.2 Status of η.** Based on S3 and S6, each background is classified into F2 (ii), (ii′) or (iii):
  - The f → 0 limit is examined symbolically: whether the constraint d/dη[J²/(2M²) − (Λ/4)N⁴] = 0 can be solved locally for η as a function of the gradients (implicit function theorem: V_eff″ ≠ 0), and which effective operator results.
  - The regularity of the limit is a result, not an assumption.
  - The full k-dependent dispersion (c₁₄ and c₁₂₃) is used, not m² alone.
- **S7.3 Comparison with B4.** In the u frame, with η̄ fixed, the leading local effective theory for (φ, ψ) is compared with B4 restricted to s = κ²/p:
  - **operator by operator**, by the symbolic coefficients of the effective Lagrangian in (u·∂φ_i, n·∂φ_i);
  - **by the local characteristic polynomial** F of B4 (`docs/lambda_eff_derivation.md`).

  The Λ term is checked for sign and normalisation (1 + 3Λq² in d1 when η̄ = 0). General B4 (s ≠ κ²/p) is reported as outside the reach of the construction (F3 item 2), with the value of s − κ²/p per background.
- **S7.4 Outcome class** per background and f: REPRODUCED-RANK-ONE, NOT-REPRODUCED, CONDITIONAL or INCONCLUSIVE (Part F).

### S8 — Einstein-aether variant (implemented, not run)

All of S1–S7 are parametrised by (c₁₄, c₁₂₃), with c₁₄ = c₁₂₃ = f the default. In the variant:
- the δη kinetic block is c₁₄ω′² − c₁₂₃k′² (ū frame);
- the full principal symbol (S3.2) is derived for the whole coupled system;
- m²_η,eff uses 1/c₁₄.

The code path is covered by S0 and by a self-test at c₁₄ = c₁₂₃ = f, which must reproduce the default results exactly (**DC2**). No other values are run.

## 4. Internal checks (mandatory, distinct from the physics outcomes)

| Check | Content | Pass condition |
|---|---|---|
| DC0 | V1 versus the V2 tangential projection | identical (symbolic difference 0) |
| DC1 | K(ω, k): route R1 versus route R2 (procedure S2.3) | Every converged entry agrees under the combined absolute-relative comparison rule defined in S2.3; any R2 NOT CONVERGED entry gives INCONCLUSIVE |
| DC2 | aether path at c₁₄ = c₁₂₃ = f versus dynamic η | exactly identical outputs |
| DC3 | S0.1 reduction and identities | all symbolic differences 0 |
| DC4 | background enumeration: Sturm count = number of located roots; every root passes C0 | equality; C0 residual ≤ 10⁻¹⁰ |
| DC5 | the limit Λ = 0. There the η equation is J J′ = 0, so its roots are J′ = 0, which must reproduce tanh 2η̄ = −B/A where |B| < |A|, **and** J = 0 (σ̄ = 0, tanh 2η̄ = −A/B where |A| < |B|). The J = 0 roots also solve the U1–U4 equation σ̄J′ = 0 but are not counted in the recorded 17/41. They are reported as a note on the historical baseline, not as a correction of it | J′ = 0 roots identical to 10⁻³⁰; J = 0 roots listed |
| DC6 | no NaN; every reported number traceable to a named step | — |

`ALL_DC_PASS` is reported for DC0–DC6 only. It is not a physics verdict and not an independent proof.

## 5. Outputs

New folder `results_uv_dynamic_eta/` with:
- `uv_dynamic_eta_results.json`: all S0–S8 results, per background, per root η̄ and per f;
- `uv_dynamic_eta_output.txt`: the console log;
- `RUN_CONTEXT.txt`, `SHA256SUMS.txt`, `started_at.txt`, `finished_at.txt`, `exit_code.txt`, `console_output.txt`, as for the B4.12 runs;
- a copy of the executed script.

The physics outcomes (T3-fund without C4, T5, T6, T7) are reported **separately**, per background set and per f, with counts of every class. No single verdict. Degrees and counts are written with `int()` before JSON, and DiracDelta is not used.

## 6. Runtime and failure reporting

The 50- and 100-digit root work over 92 k-points × 41 backgrounds × up to several roots η̄ × 3 values of f may be slow. Its expected size is reported in the static-review package. If a run is slow or fails, the exact line and the partial output are reported. No grid, tolerance or background set is reduced to save time without a new reviewed version.

## 7. Out of scope

- Script B (C4, time domain): separate specification.
- Running the Einstein-aether variant: needs registered (c₁₄, c₁₂₃) values and criteria.
- T7b, T7a, B4.12.
- An η potential, a second heavy field, time-dependent η backgrounds.
- Loop effects. T4 and T8.
- Any change of a registered criterion, tolerance or background set.

## 8. Requested from the reviewer

1. Approval of the split into script A (this document) and script B (C4, later).
2. Approval of V1 + V2 as two routes (S0.2), and of R1 + R2 for K(ω, k) (S2.3).
3. Approval of the exact enumeration of η̄ by z = e^η and Sturm counting (S1.2).
4. The threshold 10⁻² M² and the three exclusive classes of S6.2 (adopted from the review of 2026-11-16).
5. Confirmation that the aether variant is implemented but not run (S8).

**B4.12: OPEN. T7a: blocked. T7b: not started. No run authorised. Zenodo: paused.**
