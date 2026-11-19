# BEC Stage 2: procedural annex to `STAGE2_PREREGISTRATION.md` (v2) — v0.7

**Status:** PROPOSAL FOR REVIEW. Specification. The Stage 2 script and its tests exist and are under static review (script v0.3). **None of them has been executed.** Only `py_compile` has been run, a syntax check. Nothing from Stage 2 has been run, searched, fitted or sampled. The only code executed in this work is the authorised symbolic check of A1 and A3. Not registered.

**Base:** 73347c3. **Supersedes:** v0.6 (5227003e…); v0.5 (fdaaab6c…); v0.4 (c59e1e55…); v0.3 (849ebffe…); v0.2 (33154dca…); v0.1 (ad965017…).

**Changes from v0.6 (follow the reviews of 2026-11-19 of annex v0.6 and script v0.2; no criterion or threshold of v2 changes):**

| § | Change |
|---|---|
| Status, §H | Updated: the script exists and is not executed; the stale requests are replaced. |
| J1 | "Not generated" is reported separately from "generated but not recovered". The threshold remains over the fixed number, 29. |
| J4 | The conclusion is limited to the **demonstrated sensitivity**. A NO MATCH is a result of the fixed procedure, not a proof that no match exists anywhere in the family. |
| J6 | States explicitly that the control does not guarantee equal sensitivity for every real B4 target. |
| §I R4 | Only the registered message of `lambda_eff.A_coefficients` counts as complex roots. Every other exception is numerical_error. |
| §I R10 | The three diagnostic states are kept distinct. |

**Changes from v0.5 to v0.6 (kept for the record; additions only):**

| § | Change |
|---|---|
| new §J | **Positive control: injection–recovery.** Synthetic targets are generated from known points of the family, and the same search is run on them. A recovery threshold and its consequence for reading a NO MATCH are fixed **before** any run. Without it, a NO MATCH on the B4 targets could not be distinguished from a failure of the search. |
| new §K | **The author's stated expectation before any run.** It is not a criterion and does not affect any verdict. |
| G | Order of work: the script must implement §J (script v0.3) and pass a further static review. |

**Changes from v0.4 to v0.5 (kept for the record). These follow the static review of script v0.1 on 2026-11-19, so that the annex and the script agree. No criterion or threshold of v2 changes.**

| § | Change |
|---|---|
| A2 | The zero-momentum sensitivity is **removed** (review decision on R2). It needs n₂/n₁, which the family does not fix. Only the symmetric frame is used. |
| B3 | SciPy is **pinned**. The B2a argument depends on the step-acceptance rule of this implementation. |
| new §I | The implementation decisions from the static review (R1, R3, R4, R5, R7, R8, R10, B5 outcome, C4 reasons) are recorded, so that the registered annex fixes them. |

**Changes from v0.3 to v0.4 (kept for the record) (§G only):**

| § | Change |
|---|---|
| G | The order of work now separates **preparing** the script, authorised after the annex was approved, from **running** it, which stays blocked until registration and a separate authorisation. The annex is not registered, and this version does not claim that it is. |

**Changes from v0.2 to v0.3 (kept for the record) (one correction, plus one implementation note from the same review):**

| § | Change |
|---|---|
| C1c | Before the rank-one check, each of the 29 targets is checked for p ≥ 0 and s ≥ 0. If p = 0 or s = 0, κ = 0 is required and reported explicitly. A negative p or s, a rank-one inconsistency, or p·s = 0 with κ ≠ 0 causes ABORT. No target is excluded individually and no other model is used. |
| B5 | The penalty counter wraps **every** call of the residual function, including the calls made for the finite-difference Jacobian. The penalty/ABORT rule for non-finite values applies to every call. |

**Changes from v0.1 to v0.2 (kept for the record) (only the three required clarifications, plus the record of the symbolic check):**

| § | Change |
|---|---|
| A1 | λ_min(Ẽ) is now the **operational definition** chosen for the v2 energy test, and it is stated as a concretisation of v2. Sylvester's law is cited only for the sign pattern. The reasons for the choice are given. An inertia-only control count is added. |
| A1, A3 | The symbolic check is recorded: report 995771de…, 15/15 PASS. "k = 0" means the limit k → 0⁺. |
| B2, B5 | Hard rejection is kept. It is now stated why the solver stays in the target region: every feasible cost is below the penalty cost, and the solver accepts only cost-decreasing steps. Non-finite residuals are mapped to the penalty. Penalty evaluations are counted per start. |
| C1 | Perturbed targets with p ≤ 0 or s ≤ 0 are rejected. Rules are added for κ = 0 and for verifying the rank-one property of each target. |

No criterion, threshold, class, matching rule, BEC family, f bound, seed, start count or sampling count of v2 is changed.

**What this annex does not change.** It does not change `docs/STAGE2_PREREGISTRATION.md` (git blob 39f66ce2…, sha256 8cf6dd9a…). In particular, it leaves unchanged:
- the four "match" criteria and their thresholds;
- the A/B/excluded/boundary classes and their thresholds;
- the branch-matching rule (primary = ordered by root value; the 24 permutations are a separately labelled sensitivity result);
- the BEC family and its parameter f ∈ [0, 0.99];
- seed 20261111;
- 80 starts per background;
- 2·10⁵ sampling points per region.

**What it does.** It closes the procedural gaps A–E named in the review of 2026-11-19. Where v2 is silent, the choice made here is marked **PROPOSED** and is open to change before registration. Where a choice follows from v2 or from registered code, the source is given.

## 0. History of the attempt (kept, not erased)

v2 discloses an earlier v1 run:
- it was started and stopped after 3 of 29 rows;
- best residuals were about 1.0 in region A;
- its outputs were deleted and not used.

This annex records the following in addition.
- **The author has seen these partial values.** Deleting the outputs does not remove that knowledge.
- **The v1 settings are not available.** Its bounds, objective and solver settings were not kept in the repository. So it cannot be stated whether the procedural choices below coincide with v1.
- **The choices below are made after that pilot.** None of them changes a criterion or a threshold of v2. A reviewer who considers any choice outcome-sensitive can fix it differently before registration.

## A. Energy diagnostic (classification into A and B)

**A1. Definition (checked symbolically: `STAGE2_ANNEX_A1_A3_SYMBOLIC_CHECK_REPORT_v0.1.md`, 995771de…, items A1.0–A1.6 PASS).**
- Take E − Σμ_iN_i of the two-component GP functional in the **lab frame**, to second order in the fluctuations (δn_i, δu_i). Here δu_i = (ħ/m_i)∂_xδθ_i.
- For a Fourier mode k, the result is the real symmetric 4×4 form
  E₂ = Σ_i [½m_in_i δu_i² + m_iV_i δn_i δu_i + ½(g_ii + ħ²k²/(4m_in_i)) δn_i²] + g₁₂ δn₁δn₂.
- **Normalisation.** Apply the congruence D = diag(g₁₁^{−1/2}, (m₁n₁)^{−1/2}, g₂₂^{−1/2}, (m₂n₂)^{−1/2}). In the order (δn₁, δu₁, δn₂, δu₂), the normalised matrix depends only on the variables of the registered module and family (a_i = g_iin_i/m_i, b_i = ħ²/(4m_i²), f, V_i):

  Ẽ(k) = [[1 + b₁k²/a₁, V₁/√a₁, f, 0], [V₁/√a₁, 1, 0, 0], [f, 0, 1 + b₂k²/a₂, V₂/√a₂], [0, 0, V₂/√a₂, 1]].

- **Operational definition (a concretisation of v2, not a consequence of a theorem).** For the v2 energy test, "the min eigenvalue of the lab-frame BdG energy form" is **defined** as λ_min(Ẽ(k)), the smallest eigenvalue of the normalised matrix above. The v2 thresholds (−10⁻⁹ and −10⁻⁶) are applied to this number.
- **What Sylvester's law gives, and what it does not.** The congruence D preserves the **inertia**: the numbers of positive, negative and zero eigenvalues. So "positive definite" means the same for M and for Ẽ, and for any other representation of the same quadratic form in invertible coordinates. It does **not** preserve eigenvalue magnitudes. The numerical thresholds therefore have a meaning only in a stated normalisation, and this annex states one.
- **Why a normalisation must be chosen at all, and why this one.**
  - (i) The unnormalised M mixes dimensions: g_ii in the δn entries and m_in_i in the δu entries. Its eigenvalues are not unit-free, so the dimensionless v2 thresholds cannot be applied to it as it stands.
  - (ii) Ẽ depends only on the family variables (a_i, b_i, f, V_i). It reduces to the identity at f = 0, V_i = 0, k → 0.
  - (iii) In the variables (δn_i, δu_i), the k-dependence enters only as a positive diagonal shift (A3). In the variables of the Bogoliubov amplitudes (δψ, δψ*), the phase direction costs energy ∝ k². There λ_min → 0 as k → 0 for every state, and a fixed threshold such as −10⁻⁹ would depend on the smallest k on the grid. This last point is a hand argument and is **not** covered by the symbolic check.
- **Where the choice matters.** For any normalisation, A versus "not positive semidefinite" is the same outside a neighbourhood of λ_min = 0. The choice affects only points whose λ_min lies near the thresholds: class A between −10⁻⁹ and 0, and class 4 in the excluded band.
- **Control count (reported, not used for classification).** For every classified point, the script also records the inertia-only verdict, i.e. whether Ẽ(0) is positive definite. It reports the number of points classified A with λ_min(Ẽ) in [−10⁻⁹, 0], where the inertia verdict is "not positive definite". This count measures how much the class-A result depends on the operational tolerance.
- **The v1 pilot.** Its implementation of the energy form is not available (§0). Whether v1 used this normalisation cannot be stated. The definition here is fixed now, after the pilot, and is recorded as such.

**A2. Lab frame (PROPOSED).**
- The family fixes only dV = (V₁ − V₂)/√a₁. The energy form depends on V₁ and V₂ separately.
- **Primary convention:** the symmetric frame, V₁ = +dV·√a₁/2 and V₂ = −dV·√a₁/2. This is the convention of Hamner et al., where u₁ = −v/2 and u₂ = +v/2, up to the sign of dV.
- **Zero-momentum sensitivity: removed in v0.5.** It would need the density ratio n₂/n₁, which the registered family does not fix. Fixing n₁ = n₂ without a physical basis was rejected in the static review (R2). The symmetric frame is the only frame used.

**A3. The k grid.**
- Proposed grid: k = 0 together with 200 points geometrically spaced in the dimensionless variable K = k√(b₁/a₁) ∈ [10⁻³, 10²].
- **Analytic result (symbolic check 995771de…, items A3.1–A3.6 PASS; monotonicity by the Courant–Fischer theorem, which is cited, not re-derived).** k² enters Ẽ only through the positive terms b_ik²/a_i on the diagonal. So Ẽ(k₂) − Ẽ(k₁) is positive semidefinite for k₂ > k₁, and λ_min(k) is non-decreasing in k². "k = 0" denotes the limit k → 0⁺ of the k ≠ 0 form, not a uniform change of the flow. The grid minimum is therefore attained at k = 0. The grid is kept as a consistency check: any grid point with λ_min(k) < λ_min(0) − 10⁻¹² is an implementation error and causes ABORT.
- **Consequence.** In this form, energetic stability reduces to positive definiteness at k = 0, a Landau-type condition. For one component it gives V² < a, the Landau criterion at k → 0.

**A4. Numerics (PROPOSED).**
- Eigenvalues by `numpy.linalg.eigvalsh` in float64 on the normalised 4×4.
- A point with any non-finite entry is classified "boundary / ill-conditioned" and counted.
- Values of λ_min between −10⁻⁶ and −10⁻⁹ are class 4, excluded, as in v2.

**A5. Dynamic class.**
- Uses `theory/bec_dispersion.py` (blob 16112643…, sha256 228086f4…) unchanged.
- `hydrodynamic_roots` gives four real simple roots, or raises an error.
- The raise reasons map to the v2 classes:
  - complex roots → class 3;
  - repeated or ill-conditioned → class 4.
- f² ≥ 1 cannot occur, since f ≤ 0.99 by v2.

## B. The numerical search

**B1. Parameter space and transformation (PROPOSED bounds; f as in v2).**

| parameter | search variable | bounds | starting distribution |
|---|---|---|---|
| ρ = a₂/a₁ | log₁₀ρ | [−2, 2] | uniform |
| f = g₁₂/√(g₁₁g₂₂) | f | [0, 0.99] (v2) | uniform |
| dV = (V₁ − V₂)/√a₁ | dV | [−4, 4] | uniform |
| κ_b = b₂/b₁ | log₁₀κ_b | [−2, 2] | uniform |

- **Fixed units:** a₁ = 1, b₁ = ¼, m₁ = ħ = 1. Ratios carry all the information.
- **Why both signs of dV.** The B4 roots are not symmetric under c → −c, and the unit-free x_j is invariant only under affine maps with positive slope.

**B2. Objective (PROPOSED).**
- Residual vector, with branches matched by the primary rule:
  r = [(x₂^B4 − x₂^BEC)/0.05, (x₃^B4 − x₃^BEC)/0.05, ε₁/0.25, ε₂/0.25, ε₃/0.25, ε₄/0.25],
  with x_j and ε_j as defined in v2. The scales are the v2 thresholds, so r ≤ 1 componentwise corresponds to criteria (2)–(3).
- The sign criterion (1) is **not** in the objective. It is evaluated afterwards.
- **Hard rejection.** The residual is the constant vector 10³·(1, …, 1) in any of these cases:
  - the point is not in the target region (A for the primary search, B for the control);
  - the BEC roots or β are not defined;
  - any residual component is non-finite.

  This rule is kept as the admissibility rule. It is not replaced by a smooth penalty.
- **Why the solver stays in the target region (B2a).** Two facts combine:
  1. **Every admissible cost is below the penalty cost.**
     - With x_j ∈ [0, 1], each x-residual is at most 1/0.05 = 20.
     - With b̂ ∈ [−1, 1], each ε_j is at most 2. The numerator is at most 2·max(|b̂^B4|, |b̂^BEC|), and the denominator is at least that maximum. So each ε-residual is at most 2/0.25 = 8.
     - Hence an admissible cost satisfies ½Σr² ≤ ½(2·20² + 4·8²) = 528.
     - The penalty cost is ½·6·10⁶ = 3·10⁶.
  2. **The solver is monotone.** In `scipy.optimize.least_squares(method="trf")`, a trial step replaces the current point only if the actual cost reduction is positive. This was read in the source of scipy 1.18.1, `scipy/optimize/_lsq/trf.py`: `if actual_reduction > 0: x = x_new`. The installed version is re-checked at static review.

  Starting from an admissible point, no accepted iterate can therefore be inadmissible. A trial step into an inadmissible point is rejected, and the trust region shrinks.
- **Known limitation (B2b).**
  - Near a region boundary, finite-difference Jacobian probes can land on the penalty. The Jacobian is then dominated by the jump, and the solver can stall. With `x_scale="jac"` this also distorts the scaling.
  - This affects efficiency, not admissibility.
  - It is made visible, not hidden: see B5.

**B3. Solver (PROPOSED).**
- `scipy.optimize.least_squares(method="trf")`, with bounds as in B1.
- `x_scale="jac"`, `ftol = xtol = gtol = 1e-10`, `max_nfev = 2000`.
- **SciPy pinned (v0.5):** version 1.18.1, and the sha256 of `scipy/optimize/_lsq/trf.py` = 0c7c29d5d3fd6915720e9547f051d47e14ca5a821e7872686aa8cf1d1c4e4893. Both are checked at start-up, before any computation; a mismatch is an ABORT. This is the implementation whose step-acceptance rule (`if actual_reduction > 0: x = x_new`) supports B2a. The versions are also recorded in RUN_CONTEXT.

**B4. Starts and seeds (PROPOSED).**
- 80 starts per background (v2), drawn from the B1 distributions.
- Generator: `numpy.random.Generator(PCG64(SeedSequence([20261111, region_code, background_index])))`, with region_code 0 for A and 1 for B.
- A start outside the target region is redrawn up to 10³ times. If no admissible start is found, the background is recorded as "no admissible start" for that region.

**B5. Convergence and reporting.**
- Per start, record: solver status, solver message, nfev, final cost, region membership of the final point, whether it was rejected, and **the number of residual evaluations that returned the penalty** (n_penalty_eval, counted by a wrapper around the residual function). The wrapper counts **every** call, including the finite-difference Jacobian calls of the solver. The non-finite → penalty rule of B2 is applied inside the wrapper, so it covers every call.
- If any final point is inadmissible although its start was admissible, B2a has been violated. This is reported as an implementation error and causes ABORT.
- **Converged** means status > 0 and the final point lies in the target region.
- The **best point** is the converged start with the lowest cost. Criteria (1)–(4) are evaluated on it.
- Per background, also report the number of non-converged and rejected starts, and the number of starts with n_penalty_eval > 0 (possible boundary stalling, B2b).
- A background with no converged start in the region is reported as **NO ADMISSIBLE FIT**. It counts as "no match".

## C. The robustness test (criterion 4)

**C1. What is perturbed (PROPOSED).**
- The B4 target inputs w, q, v, r, p, s and Λ are perturbed independently.
- κ is **not** perturbed independently. It is set to sign(κ)·√(ps) after the perturbation, with sign(κ) taken from the unperturbed target. The reason is that the registered kernel `theory/lambda_eff.py` is valid only for the rank-one family κ² = ps.
- **C1a. ps < 0, or p, s not positive.** In the rank-one Proca construction, p = g_φ²/M² and s = g_ψ²/M² are non-negative. A perturbed target with p ≤ 0 or s ≤ 0, where the unperturbed value was positive, is **inadmissible**. It is rejected and redrawn under the C4 rule, and it counts towards the 500-draw cap. √(ps) is never evaluated for ps < 0. No non-rank-one kernel is used.
- **C1b. κ = 0.** Under rank one, κ = 0 occurs only if p = 0 or s = 0.
  - For such a target, a coefficient that is exactly 0 is **not perturbed**, and κ stays 0. The other inputs are perturbed as in C2.
  - A target with κ = 0 but p > 0 and s > 0 contradicts rank one. It is treated under C1c.
- **C1c. Validity of the target inputs (before any search; in this order).**
  1. **Physical admissibility.** p ≥ 0 and s ≥ 0 for every one of the 29 targets, as required by p = g_φ²/M² and s = g_ψ²/M². A rank-one check alone would not exclude, for example, p = s = −1 with κ = 1.
  2. **Degenerate case.** If p = 0 or s = 0, then κ = 0 is required. The number of such targets is reported explicitly.
  3. **Rank-one consistency.** |κ² − ps| ≤ 10⁻¹²·max(κ², ps, 10⁻³⁰⁰).

  **Any failure of 1–3 causes ABORT** of the whole run. A failing target is not excluded individually, and no other kernel or model is substituted. This is a validity check of the registered input data, not a change to B4.

  Original wording of the rank-one step, kept for traceability: the selection script checks every one of the 29 targets for |κ² − ps| ≤ 10⁻¹²·max(κ², ps, 10⁻³⁰⁰). For comparison, the anchor's recorded difference is 3.5·10⁻¹⁸. A failure causes **ABORT**, because the kernel would not be valid for that target. The script also reports the number of targets with κ = 0. It is expected to be zero, but this is not assumed.

**C2. Size and distribution.**
- Each perturbation is δ = max(0.02·|x|, 0.005)·u, with u ~ U[−1, 1] independent per input.
- This is the reading of "relative 2 %, absolute 0.005": the larger of the two, with a uniform sign and magnitude.

**C3. No refit.** The best BEC point from B5 is held fixed. Criteria (1)–(3) are re-evaluated against each perturbed target.

**C4. Admissibility of perturbed targets.**
- A perturbed target is rejected and redrawn if it loses four real simple local roots, i.e. if `A_coefficients` raises an error.
- Draws continue until 50 admissible targets exist, up to 500 draws.
- If fewer than 50 admissible targets exist, criterion (4) is **INCONCLUSIVE** for that background, and the number is reported.
- Seed: SeedSequence([20261111, 2, background_index]).

**C5. Threshold.** Criterion (4) holds if (1)–(3) hold for at least 80 % of the admissible perturbed targets (v2).

## D. Sampling of sign tuples

- **Distribution:** the B1 distributions.
- **Seed:** SeedSequence([20261111, 4, region_code]).
- **Points:** draw until 2·10⁵ points have been classified into region A, and separately into region B. The cap is 2·10⁷ draws. The acceptance rate is reported.
- **Sign tuple:** the tuple of sign(b̂_j), ordered by root value. An entry with |b̂_j| < 0.05 is recorded as 0.
- **Status of the result:** empirical. A sign tuple not seen is "not observed in the sample", not "impossible".

## E. Identity and provenance of the targets

| item | path | git blob | sha256 |
|---|---|---|---|
| manifest | `data/e3_manifest.json` | f0a26e10… | 1bae892a08c3de1f928e4057ee0ad801294a90462478aaedf5b9c724e23d87e0 |
| anchor | `data/B3_max_speed_7_315096.json` | 9bb3faa4… | 99a6b04bc34821ddd5c5d58a78b3474d369aa9fc0ee264163096398b0028ccbd |
| B4 kernel and β | `theory/lambda_eff.py` | 6d9a8b37… | e07299762f3e6b89187f8dc8270b430c951e15816ea129983f040dd2f80c951f |
| BEC module | `theory/bec_dispersion.py` | 16112643… | 228086f4a1a3ab292c1798910dffcc4ba69cec9e17c88dd7046369f52b878865 |

- **Selection rule.** The targets are the anchor plus the manifest's `selected` rows, kept if `lambda_eff.A_coefficients` returns four real roots (|Im| ≤ 10⁻⁹, as coded there).
- **Expected count: 29** (v2). The selection script runs first and prints the indices. If the count is not 29, it ABORTS. The list of indices is recorded, with its hash, **before** any search.
- **β and M.** β_i = A_i/M² with A_i from `beta_M2`. Because b̂ = β/max|β|, M cancels. No value of M is needed.
- **Conditioning flags.**
  - `lambda_eff.conditioning` returns |F′(c_i0)|, the relative root gap, and near_zero_root (|c| < 0.05·max|c|).
  - **PROPOSED:** flags are reported per target. Targets are **not** excluded on them, since v2 does not say to exclude. A target with relative gap below 10⁻⁶ or |F′| below 10⁻² (the v2 BEC-side thresholds) is labelled ill-conditioned, and its result is shown separately as well as in the totals.
  - For example, the anchor root −0.175 is below 0.05·7.315, so the anchor carries near_zero_root.
- **Input verification.** All four files are verified by sha256 at start-up, as in Script A rev1.1. Any mismatch causes ABORT before any computation.

## F. Informative note: where Hamner et al. (2011) lies (no change to scope)

- **Hamner has f = √(σ₁σ₂) ≈ √0.98586 ≈ 0.9929 > 0.99.** It lies **outside** the v2 family bound f ≤ 0.99.
- Its other family coordinates, in its own units, are:
  - ρ = 0.3;
  - κ_b = 1, since the masses are equal;
  - dV = v/√ρ₁ ≈ 0.093 at v = 0.17, and ≈ 0.176 at v = 0.32.
- Stage 2 as registered tests B4 targets against the BEC family. It does not test the Hamner mixture.
- Whether to extend f is a change to v2 and is not proposed here. If the reviewer wants it, it must be registered as v3 **before** any run.

## G. Outputs and order of work

**Order:**
1. this annex;
2. review;
3. symbolic check of A1 and A3 (done: 995771de…, 15/15 PASS);
4. **writing** the script and its tests (v0.3 must implement §J). Authorised on 2026-11-19 after approval of the annex in substance. This is preparation only: no search, no fit, and no run on the targets;
5. static review of the complete code against v2 and this annex, with corrections where needed;
6. registration of v2 together with the annex, and of the reviewed script by hash. A separate decision; not done and not assumed;
7. a separate authorisation, then one run on a fresh clean LF clone after input-hash verification;
8. review of the results.

**Execution of the script remains blocked** until steps 6 and 7 are complete.

**Outputs:** `results_stage2_bec/`, containing:
- the target index list with its hash;
- per-background tables (best point, criteria (1)–(4), convergence counts, flags);
- the sign-tuple tables;
- RUN_CONTEXT.txt and SHA256SUMS.txt.

**Out of scope:** LHY; nonlinear operators; a common coupling family across backgrounds; non-rank-one B4; the ψ mass term; Hamner as a target (§F); EXT-1.

## I. Implementation decisions fixed by the static review (v0.5)

- **R1. The 24 permutations.**
  - Evaluated at the same best point of the primary search, with no re-search.
  - For an assignment π, B4 branch j is compared with BEC branch π(j) in (1)–(2). The two inner B4 roots are compared with their assigned BEC roots in (3).
  - Reported as the number of passing assignments and the best one, labelled "more freedom; not a match".
- **R3. background_index.** Anchor = 0, then the manifest `selected` rows in file order, after filtering.
- **R4. Finiteness on the B4 side, and classification of failures (v0.7).**
  - **complex_roots** only when `lambda_eff.A_coefficients` raises ValueError with exactly its registered message, "local roots are complex".
  - **numerical_error** for every other exception: any other ValueError, LinAlgError, RuntimeError, TypeError and so on. This applies to an exception raised anywhere in the parameter construction, the A_coefficients call or the array conversion, and also to output of the wrong shape.
  - **nonfinite** for non-finite c₀ or A_i.
  - **degenerate** for zero root span or all A_i = 0.
  - **Base targets:** a candidate lost for any reason other than complex roots (non-finite c₀ or A_i, a numerical error, or zero span or zero β) is an ABORT, never a silent exclusion.
  - **Perturbed targets:** each rejection reason is counted separately. The reasons are: C1a; complex roots; numerical error; non-finite; degenerate.
- **R5. Scope of the A1 control count.**
  - Counted over all class-A points of the D sampling.
  - Stored per best search point.
  - Not computed for every residual evaluation inside the searches.
- **R6. Region B.** It runs the same pipeline as region A and is reported separately as the control. It is not part of the primary result.
- **R7. Verdict.**
  - MATCH = (1) ∧ (2) ∧ (3) ∧ (4).
  - If (1)–(3) pass but (4) cannot be evaluated (fewer than 50 admissible perturbed targets), the verdict is **INCONCLUSIVE (criterion 4)**. It is never MATCH, and it is not a demonstrated NO MATCH.
  - Otherwise NO MATCH. If there is no converged admissible start, NO ADMISSIBLE FIT.
- **R8. H.** H is exactly `bec_dispersion.hydrodynamic_polynomial` (ascending coefficients of [(c−V₁)² − a₁][(c−V₂)² − a₂] − γ). |H′| is evaluated at the four roots with `polyder` of that polynomial, as an absolute value in units a₁ = 1. This is in addition to the module's own relative-gap and slope checks.
- **R10. Solver output checks, before the cost comparison.**
  - sol.cost, every component of sol.x, and (for an admissible final point) every component of the final residual must be finite.
  - The recomputed ½|r|² must equal sol.cost to relative 10⁻⁹.
  - Any violation is an ABORT.
  - **Three distinct diagnostic states (v0.7), never merged:**
    - (a) an **inadmissible final point**: per-start outcome `rejected_start_shifted_out` or `rejected`, not an ABORT unless B2a is violated;
    - (b) **non-finite** cost, x or final residual: ABORT, with an R10 non-finite message;
    - (c) a **cost mismatch** between the recomputed and the reported cost: ABORT, with an R10 cost-mismatch message.
- **B5. Explicit per-start outcome.** Each start records `rejected` (final point inadmissible) and exactly one `outcome`:
  - converged;
  - not_converged_in_region;
  - rejected_start_shifted_out (R9: the first counted evaluation was already a penalty);
  - rejected (excluded by the B2a ABORT).

  "no admissible start" is recorded separately.

## J. Positive control: injection–recovery (v0.6)

**Purpose.** To measure whether the fixed search (B1–B5) recovers a match when one is known to exist. A NO MATCH on the B4 targets is interpretable only if this control passes.

**J1. Synthetic targets.**
- For each region R ∈ {A, B}: 29 synthetic targets, the same number as the real ones. Synthetic target i = 0…28 is drawn from the B1 distributions with the generator SeedSequence([20261111, 6, region_code, i]), with redraws until the point is in region R, up to 10⁴ draws.
- Its target data are the BEC spectrum at that point: x_j from the sorted roots c₀, and b̂_j from β (`bec_dispersion.beta_coefficients`), in the same form as the B4 targets.
- The generating point y* is stored but is **never** passed to the search.
- If no point in region R is found within 10⁴ draws, synthetic target i is recorded as "not generated" and counted as **not recovered** for the threshold.
- **Separate reporting (v0.7).** Sampling difficulty and optimiser failure are different problems, so each region reports:
  - `n_target_fixed` (29);
  - `n_generated`;
  - `n_not_generated`;
  - `n_recovered`;
  - `n_generated_not_recovered`;
  - `recovery_rate` = n_recovered / 29.

  The threshold stays over the fixed number 29, so the rule does not depend on how many targets were generated.
- No real B4 target is used or looked at to make the synthetic targets.

**J2. Search.**
- Exactly B1–B5, with the same bounds, objective, solver, 80 starts and region-membership rule.
- The start seeds use background_index = 1000 + i, so they differ from the real targets' seeds (0…28) and from the generating stream.

**J3. Recovery.**
- A synthetic target is **recovered** if criteria (1)–(3) of v2 hold at the best converged point (primary branch matching).
- Criterion (4) is not part of the control, because the perturbation of C1–C2 acts on B4 inputs and a synthetic target has none.
- The recovered point need not equal y*: the map from family to spectrum need not be injective. |y − y*| is reported, not used.
- For each region, also report: the number of starts with each outcome, n_penalty_eval, and the cost at the best point.

**J4. Threshold, fixed before any run.**
- For each region separately: recovery rate ≥ 0.90, i.e. at least 27 of 29.
- **If the control passes in region R:** a NO MATCH on the B4 targets in region R is reported as a NO MATCH, with this reading: the fixed, registered procedure found no match, and the search has **demonstrated sufficient sensitivity on a pre-specified set of attainable synthetic targets**. This makes the negative result interpretable **within that demonstrated sensitivity**. It does **not** prove that no match exists anywhere in the admissible family. A local optimiser with 80 starts cannot prove non-existence.
- **If the control fails in region R:** every NO MATCH and NO ADMISSIBLE FIT in region R is reported as **NOT INTERPRETABLE (insufficient search sensitivity)**. A MATCH found in region R is still reported, since the control measures only sensitivity, with a note that the control failed.
- **No change after seeing the control.** The search, its settings and the threshold are not changed after the control has been seen. Any change requires a new registration and a new run.

**J5. Order inside the run.**
- The positive control runs **after** target selection and C1c, and **before** any real target is evaluated.
- Its results are written to `positive_control.json`, and its sha256 is recorded in RUN_CONTEXT, before the first real search starts.
- The run does not stop if the control fails. The real targets are still evaluated and reported under J4.

**J6. Limitation, stated in advance.**
- The synthetic targets lie exactly inside the family, so the control measures sensitivity to **attainable** targets.
- A B4 target near the family's boundary, near f = 0.99, or near the region boundaries may be harder to recover than the control suggests.
- **The control does not guarantee equal sensitivity for every real B4 target (v0.7).** It shows sensitivity on synthetic targets drawn from the family interior. A real target that is harder to reach may be missed even when the control passes.
- This limitation is reported together with the result. A passed control is never quoted without it.

## K. The author's stated expectation before any run (not a criterion)

This is recorded so that the result can later be compared with what was expected. It does not change any criterion or verdict.

**Expected outcome.** **NO MATCH** in region A for most or all of the 29 targets. The expectation for region B is less clear.

**Reasons:**
- **Root structure.** The registered B4 targets are strongly asymmetric non-perturbative backgrounds. The anchor's roots are −4.916, −0.175, 2.797 and 7.315, which give x ≈ (0, 0.388, 0.631, 1).
  - In the BEC family, the roots are Doppler-shifted sound pairs.
  - In region A, the Landau-type condition (A3) keeps the flow small compared with the sound speeds. That limits how asymmetric the root pattern can be.
  - Larger asymmetry requires larger counterflow, which tends to push points into region B or into dynamic instability.
- **Pilot evidence.** The v1 pilot (§0) gave best residuals of about 1.0 in region A on 3 of 29 rows, i.e. at the edge of the thresholds. Those settings are unknown, so this is a weak indication only.
- **Earlier records.** `UV_CANDIDATE_ASSESSMENT.md` records for the BEC mean field that it has "no κ operator" and Λ_BEC < 0. This is not tested by Stage 2, which compares roots and β only, but it does not point towards a match.

**What a NO MATCH would mean, if the positive control passes.** The rank-one B4 backgrounds are not reproduced, at the level of characteristic roots and β signs, by the mean-field two-component BEC family with f ≤ 0.99 in the symmetric frame without LHY. It would **not** mean that B4 is wrong in general, or that no other UV class or BEC regime realises it. It would close this particular identification.

**What a MATCH would mean.** A necessary condition holds, and nothing more (v2). The next step would be matching at the operator level with one fixed dictionary across backgrounds (`matching_protocol.md` §4–5).

## H. Requested from the reviewer

Decisions of 2026-11-19 on v0.1 (recorded):
- A2, B1, B3, B4, C3–C5, D and E were approved.
- A1, B2 and C1 were approved subject to the clarifications made in this version.
- The f bound stays as in v2 (§F).
- The symbolic check of A1 and A3 was authorised and has passed (995771de…).

Decisions of 2026-11-19 on v0.2 (recorded):
- A1 (operational definition and control count), A3, B2, B5, C1a and C1b were accepted.
- C1c required the p, s ≥ 0 check, which is added in v0.3.

Decisions of 2026-11-19 on v0.3–v0.6 (recorded):
- v0.3 was accepted in substance, and writing the script (Phase 3) was authorised. That authorisation did not cover any run, search, fit or registration.
- §G was corrected in v0.4.
- Script v0.1 had a static review; R2 led to the A2 change in v0.5.
- In v0.6, §J was accepted in principle, with threshold 27/29, subject to the J1, J4 and J6 clarifications made here. §K was accepted without change.

Requested now (v0.7):
1. Accept the J1, J4 and J6 clarifications and the §I R4 and R10 wording.
2. Do a limited static review of script v0.3, which implements §J and the R4 correction, together with its tests.
3. Then, separately, decide on registration and on authorisation to run the tests. Neither is requested or assumed here.
