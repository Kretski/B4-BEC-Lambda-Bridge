# T3-fund protocol: PROPOSAL, version 2 (not registered)

Status: version 2 of the proposal, awaiting review by the project lead. Version 1 (T3_FUND_PROTOCOL_PROPOSAL.md) is kept unchanged for traceability.
Nothing here is in force until it is registered in UV_ACCEPTANCE_PROTOCOL_ADDENDUM_T3_T6.md. No candidate is evaluated against T3-fund before registration.
Historical T3-B4 results are not changed by this protocol.

Changes from version 1: decisions of the project lead on the three open questions and four required clarifications (A-D), and three further decisions
(2026-11-13): energy tolerance calibrated on the free system with 1e-9 as target; no physical threshold on positive growth in C3; all five existing
constructions are tested as full models (k -> infinity). In detail:
two regimes (fundamental / EFT with declared k_max); FAIL-UNSTABLE counts as a T3-fund failure; C1 needs an analytic argument for k -> infinity (A);
C3 no longer borrows the root tolerance (B); C2 documentation requirements (C); C4 energy and norm definitions (D).

## What T3-fund tests
Under preferred-foliation causality T3-fund asks whether linear evolution in the preferred time T is well posed and stable, in the frame of the observer field.
It separates three properties: well-posed evolution, stability of the background, and bounded propagation speed. The third is NOT a T3-fund criterion.
Notation: omega_j(k) are the frequencies of the linearised system on a stationary background for real k; gamma(k) = max_j Im omega_j(k);
M is the heavy scale of the candidate; gamma_hat = gamma/M.

## Regime declaration (before the first test of a candidate)
Each candidate is declared as one of:
- FUNDAMENTAL: criteria are tested for all k, including k -> infinity.
- EFT: criteria are tested for 0 <= k <= k_max. Behaviour above k_max is not a ground for failure.
An EFT declaration is valid only if all of the following are recorded before the first test:
(1) the value of k_max; (2) its physical reason (mass of a heavy degree of freedom, threshold of new physics, or another justified scale); (3) k_max/M.
(4) Without a justified k_max the candidate is treated as FUNDAMENTAL and may not disregard bad asymptotic behaviour.
The declaration is not changed after results are seen.

## Criteria
C0 Background. Exact solution of the full equations (analytic statement, or residual <= 1e-10). Otherwise INCONCLUSIVE.
C1 Well-posedness. gamma(k) bounded on the tested range: for all k (FUNDAMENTAL) or for k <= k_max (EFT).
C2 Principal symbol (FUNDAMENTAL; for EFT only if the EFT is itself a closed evolution system up to k_max). Real characteristic roots and a complete set
   of eigenvectors. Real but defective roots do not guarantee C1.
C3 Mode stability. No growing mode on the tested range.
C4 Time-domain confirmation on at least one designated background per candidate.

Example (analytic). omega^2 = c^2 k^2 (1 + Lambda k^2):
- Lambda > 0: real for all k; C1-C3 hold in both regimes, although the front velocity is unbounded.
- Lambda < 0: gamma = c k sqrt(|Lambda| k^2 - 1) for k > k_* = 1/sqrt(|Lambda|). As a FUNDAMENTAL law it fails C1. As an EFT it fails only if k_max > k_*;
  if the declared k_max is below k_*, the sign alone is not a failure.
T3-fund therefore does not select the sign of Lambda in general; it rejects the negative sign of the uncut formula.

## Evidence levels for C1 (clarification A)
A grid cannot prove sup_k gamma(k) < infinity. For the FUNDAMENTAL regime:
- C1 holds [ANALYTIC] if an asymptotic argument is given (for example a strongly hyperbolic principal symbol, or an explicit large-k expansion of omega_j).
- Without it, a clean numerical result (P4) is recorded as "numerically passed on [k_lo, k_hi]" and the outcome for C1 is INCONCLUSIVE, not PASS.
- The fit gamma ~ k^s is a diagnostic; it can establish a failure (growth with k) but not a proof of boundedness.
For the EFT regime the range is finite; C1 is decided numerically on [0, k_max] with the grid of P2 restricted to that range.

## Numerical procedure and tolerances
P1 Exact factors of the characteristic polynomial are divided out symbolically before any root test.
P2 Grid: k/M on 61 logarithmic points in [1e-3, 1e3] and, for FUNDAMENTAL, 31 more in [1e3, 1e6]. For EFT: the points with k <= k_max, plus k_max itself.
P3 Root classification in double precision, r = |Im c| / max(1, |c|), c = omega/k:
   r <= 1e-7: numerically real;  1e-7 < r <= 1e-5: unresolved;  r > 1e-5: candidate complex root.
   A candidate complex root is NOT yet an instability: it must be confirmed by recomputing the roots of the same exact polynomial in 50-digit arithmetic
   (mpmath) after P1. Confirmed if r > 1e-12 there. Unresolved roots are also recomputed in 50 digits and reclassified; if still between the bounds: INCONCLUSIVE.
P4 Growth with k (diagnostic): where confirmed complex roots exist in [1e4, 1e6], fit gamma ~ k^s; s > 0.05 indicates ill-posedness.
P5 C2 documentation (clarification C). The verdict on C2 must record: the first-order form of the system used; the normalisation of the variables;
   the reality of the principal roots (P3); the condition number of the eigenvector matrix (indicator only, threshold 1e8, representation dependent);
   and, whenever two principal roots coincide within 1e-6 relative, an analytic argument for or against a complete set of eigenvectors.
   Without that argument a degenerate case is INCONCLUSIVE.
P6 C3 (clarification B; decision of 2026-11-13). The criterion is gamma(k) = max_j Im omega_j(k) <= 0 within numerical uncertainty. No physical threshold
   for an "acceptably small" positive growth rate is introduced: positive growth is a physical instability, not an error to be neglected.
   - Growth that exceeds the numerical uncertainty and is confirmed under refinement (roots of the exact polynomial recomputed at 50 and at 100 digits,
     and the k grid refined by a factor 4 inside the band): FAIL-UNSTABLE. Report gamma_hat_max, the band of k and the e-folding time 1/gamma.
   - Apparent growth within the numerical uncertainty (not confirmed, or changing sign under refinement): INCONCLUSIVE.
   - No positive growth on the tested range: C3 holds.
P7 C4 (clarification D).
   - Initial data: Gaussian wave packet(s) in the light fields, amplitude 1e-3 of the background gradient scale, heavy field unperturbed; fully specified in the script.
   - Norm: L2 norm over the box of all fluctuation fields (light and heavy), recorded at >= 100 times over a duration >= 1000/M; the time series is plotted.
     Test criterion: max over time of norm / initial norm <= 3. This is a criterion for this test, not a definition of stability.
   - Energy: the total conserved energy of the full system (all fields). For Lagrangian models without dissipation its conservation is a requirement on
     the integrator, not on the model. For a model with dissipation or sources the stated balance law replaces conservation.
     Target: relative drift of the total energy <= 1e-9 over the run. The target is not assumed to be attainable: it is calibrated on the free system
     (all couplings set to zero) with the same numerical method, time step and spatial grid, repeated at (2N, dt/2); the drift must decrease under refinement.
     If the free system does not reach the target, C4 is INCONCLUSIVE for the candidates, not a FAIL of the physical model.
     Implementation note (to be confirmed by the project lead): the calibration run uses the same background gradients and the same initial perturbation
     as the candidate run. A free system with no perturbation has constant fields and conserves energy trivially, so it would not test the integrator.
     Additional diagnostic, reported but not a pass criterion: drift of the fluctuation energy E_fl = E_total - E_background with the background part
     subtracted analytically. (The total energy is dominated by the background, so the total-energy drift alone is a weak test of the fluctuations.)
   - Convergence: all of the above repeated at (N, dt) and (2N, dt/2); results agree within 1 %.
P8 Backgrounds: the sets used for T2, plus the candidate's own stationary backgrounds; seeds fixed and recorded.

## Outcomes per background
- PASS: C0-C3 hold with the required evidence level (C4 on designated backgrounds).
- FAIL-ILLPOSED: C1 or C2 fails.
- FAIL-UNSTABLE: C1 and C2 hold, and a growing mode is confirmed in a finite band. Counts as a T3-fund failure.
- INCONCLUSIVE: C0 not satisfied; unresolved roots; no analytic argument for C1 (FUNDAMENTAL) or for a degenerate C2; C4 not converged.
Per candidate: the counts of the four outcomes on each background set, together with the declared regime. "Passes T3-fund on a set" is used only when
there is no FAIL and no INCONCLUSIVE on that set.

## Relation to the other tests
T2 overlaps with C2 and C3; T3-fund adds the regime declaration, the k -> infinity evidence level, completeness of eigenvectors, the separation of
ill-posed from unstable, the time-domain confirmation and the INCONCLUSIVE rules. T5 (positive energy) is sufficient but not necessary for C3.
T3-B4 keeps its historical values. T6 is independent.

## Regime of the existing candidates (decision of the project lead, 2026-11-13)
Rule: full UV constructions that contain their heavy fields explicitly are tested as full models, asymptotically to k -> infinity. Restricting the test to
k < M would test only the EFT approximation, not the construction that claims to complete it.
If inspection of the actual equations shows that a model is only an EFT with an explicit cutoff, its regime is changed BEFORE its test, with a justified k_max.
The regime is never changed after the fact to rescue a candidate with bad asymptotic behaviour.

| Candidate | Regime | Available input | Missing |
|---|---|---|---|
| Proca mediator | FUNDAMENTAL (k -> infinity) | exact characteristic polynomial; principal symbol | exact background check (C0); C2 eigenvector argument; high-precision confirmation; C4 |
| Scalar mediator | FUNDAMENTAL | principal symbol (1 - c^2) d1 d2, analytic | C1 asymptotic statement; C2 documentation; C4 |
| Scalar mediator + dynamical eta | FUNDAMENTAL | principal symbol as above; stationary eta on 17/41 backgrounds | C0 on the other 24; C1/C2; C4 |
| HeavyScalar, one light field | FUNDAMENTAL | principal symbol proportional to (omega^2 - k^2), analytic | C3 on a population; C4 |
| HeavyScalar, two light fields | FUNDAMENTAL | det K = (c^2 - 1) x quartic; Z > 0 criterion; population data; one time-domain run | C2 documentation; calibrated energy drift (P7); second resolution for C4 |

The existing results for these five candidates are NOT PASS or FAIL under this protocol. They are historical results, reused only where they cover the
new criteria; missing checks are marked as missing, not assumed.

## Decisions closed in this version
1. Energy tolerance: calibration on the free system; 1e-9 is the target (P7).
2. C3: no physical threshold on positive growth (P6).
3. Regime: all five existing constructions are FUNDAMENTAL.
Remaining point for the final review: the implementation note in P7 (calibration with a perturbation).
