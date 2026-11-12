# Addendum to the UV acceptance protocol: causal cone (T3) and preferred frame (T6)

Date: 2026-11-12. The original T1-T6 text in UV_ACCEPTANCE_PROTOCOL.md is not modified. Existing numerical results are not changed, only relabelled.
Working interpretation chosen by the project lead: option (b) of B4_CAUSALITY_ANALYSIS.md, i.e. B4 is a Lorentz-violating effective theory with a
preferred observer field (static observer / ZAMO) and Lambda > 0 is kept. This is a working hypothesis, not an established fact.

## Three speeds
- c_B4: the reference (unit) speed of the effective theory.
- c_fund: the fundamental causal limit.
- c_front: the front velocity of a given UV realisation (large-k limit of omega/k).

## T3 is split
- T3-B4: c_front <= c_B4. All T3 verdicts obtained so far (Proca, scalar mediator, scalar + dynamical eta, heavy scalar) are T3-B4 verdicts and keep their values.
- T3-fund: causality with respect to the fundamental structure. Its definition depends on one open decision (below). Until that decision is registered,
  T3-fund is "not evaluated" for every candidate, and no candidate is accepted.

## Registered decision (project lead, 2026-11-12): alternative (iii), preferred-foliation causality
> Causality Convention: The B4-Lambda framework adopts preferred-foliation causality defined by the observer field u^mu. No finite Lorentz-invariant
> upper bound on spatial front velocity is imposed at the effective-theory level. T3-fund evaluates the well-posedness and stability of evolution with
> respect to the preferred time. The physical origin, admissibility and domain of definition of the observer field are assessed independently under T6.

Accepted with it: c_B4 = c (speed of light in the normalisation used); no separate finite c_fund > c is introduced; unbounded front or group velocity is
allowed provided evolution is well posed in the preferred time. This is a working hypothesis for the interpretation, not an established fact.
Alternatives considered and not adopted: (i) finite cone c_fund = c_B4 (incompatible with Lambda > 0 for every k > 0); (ii) finite cone c_fund = R c_B4, R > 1
(requires a fundamental speed above the speed of light).

## T3-fund: what is checked
1. Real spectrum of linear excitations.
2. Well-posed Cauchy problem with respect to the preferred time.
3. Hyperbolicity, or another appropriate notion of well-posedness of the evolution system.
4. Bounded growth of solutions at fixed time (no growth rate increasing without bound in k).
5. No instability that makes the evolution physically unacceptable.
Unbounded propagation speed and unbounded amplitude growth are different properties; only the second fails T3-fund.
Note: items 1 and 3 overlap with T2. Under this convention T3-fund no longer bounds speeds, so the discriminating weight moves to T2, T5 and T6.
The operational protocol (criteria, tolerances, outcomes) is registered below under "T3-fund protocol (registered 2026-11-13)". Status for all candidates: not evaluated.

## T6 as a stand-alone task (now the central open physical problem)
It is not sufficient to specify u^mu mathematically. To be established:
1. What the observer field is physically.
2. How it is determined, dynamically or geometrically.
3. What its time foliation is (existence of a global preferred time; hypersurface orthogonality of the congruence).
4. What happens where a given congruence is not defined. Static observers do not exist inside the ergoregion of Kerr; ZAMOs remain defined outside the
   horizon. One frame cannot be assumed to cover the whole spacetime without an additional construction.
5. Its own stability and well-posedness if it is made dynamical. The dynamical-eta test already done is one negative data point for one construction.

## T3-fund protocol (registered 2026-11-13)
Registered by decision of the project lead after review of T3_FUND_PROTOCOL_PROPOSAL_v2.md. The registered text is version 2 with three operational
edits and no change of the physical criteria: (1) P7, the free calibration background must itself solve the free equations; (2) P6, an operational
definition of the numerical uncertainty; (3) P2a, a separate check at k = 0. From this point the criteria and tolerances are not changed in the light of results.
The existing results for all candidates are historical evidence, not PASS or FAIL under this protocol. Evaluation starts after this registration,
beginning with the calibration of P7 on the free system.

### What T3-fund tests
Under preferred-foliation causality T3-fund asks whether linear evolution in the preferred time T is well posed and stable, in the frame of the observer field.
It separates three properties: well-posed evolution, stability of the background, and bounded propagation speed. The third is NOT a T3-fund criterion.
Notation: omega_j(k) are the frequencies of the linearised system on a stationary background for real k; gamma(k) = max_j Im omega_j(k);
M is the heavy scale of the candidate; gamma_hat = gamma/M.

### Regime declaration (before the first test of a candidate)
Each candidate is declared as one of:
- FUNDAMENTAL: criteria are tested for all k, including k -> infinity.
- EFT: criteria are tested for 0 <= k <= k_max. Behaviour above k_max is not a ground for failure.
An EFT declaration is valid only if all of the following are recorded before the first test:
(1) the value of k_max; (2) its physical reason (mass of a heavy degree of freedom, threshold of new physics, or another justified scale); (3) k_max/M.
(4) Without a justified k_max the candidate is treated as FUNDAMENTAL and may not disregard bad asymptotic behaviour.
The declaration is not changed after results are seen.

### Criteria
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

### Evidence levels for C1 (clarification A)
A grid cannot prove sup_k gamma(k) < infinity. For the FUNDAMENTAL regime:
- C1 holds [ANALYTIC] if an asymptotic argument is given (for example a strongly hyperbolic principal symbol, or an explicit large-k expansion of omega_j).
- Without it, a clean numerical result (P4) is recorded as "numerically passed on [k_lo, k_hi]" and the outcome for C1 is INCONCLUSIVE, not PASS.
- The fit gamma ~ k^s is a diagnostic; it can establish a failure (growth with k) but not a proof of boundedness.
For the EFT regime the range is finite; C1 is decided numerically on [0, k_max] with the grid of P2 restricted to that range.

### Numerical procedure and tolerances
P1 Exact factors of the characteristic polynomial are divided out symbolically before any root test.
P2 Grid: k/M on 61 logarithmic points in [1e-3, 1e3] and, for FUNDAMENTAL, 31 more in [1e3, 1e6]. For EFT: the points with k <= k_max, plus k_max itself.
P3 Root classification in double precision, r = |Im c| / max(1, |c|), c = omega/k:
   r <= 1e-7: numerically real;  1e-7 < r <= 1e-5: unresolved;  r > 1e-5: candidate complex root.
   A candidate complex root is NOT yet an instability: it must be confirmed by recomputing the roots of the same exact polynomial in 50-digit arithmetic
   (mpmath) after P1. Confirmed if r > 1e-12 there. Unresolved roots are also recomputed in 50 digits and reclassified; if still between the bounds: INCONCLUSIVE.
P2a The point k = 0. The logarithmic grid starts at k/M = 1e-3, and the classification by c = omega/k does not apply at k = 0. The zero-wavenumber
   sector is therefore checked separately, by the analytic limit k -> 0 or by solving the characteristic system det K(omega, k = 0) = 0 directly for omega
   and classifying with r = |Im omega|/max(M, |omega|) and the same bounds as P3. Exponential growth at k = 0 fails C3 like any other mode.
   A double root omega = 0 of a massless field (solutions a + b t, a uniform shift of the background gradient) is not exponential growth; it is recorded
   as a zero mode, not as an instability.
P4 Growth with k (diagnostic): where confirmed complex roots exist in [1e4, 1e6], fit gamma ~ k^s; s > 0.05 indicates ill-posedness.
P5 C2 documentation (clarification C). The verdict on C2 must record: the first-order form of the system used; the normalisation of the variables;
   the reality of the principal roots (P3); the condition number of the eigenvector matrix (indicator only, threshold 1e8, representation dependent);
   and, whenever two principal roots coincide within 1e-6 relative, an analytic argument for or against a complete set of eigenvectors.
   Without that argument a degenerate case is INCONCLUSIVE.
P6 C3 (clarification B; decision of 2026-11-13). The criterion is gamma(k) = max_j Im omega_j(k) <= 0 within numerical uncertainty. No physical threshold
   for an "acceptably small" positive growth rate is introduced: positive growth is a physical instability, not an error to be neglected.
   - Numerical uncertainty u(k) is measured, not postulated: u(k) = |gamma computed at 50 digits - gamma computed at 100 digits| at the same k, and, on the
     refined k grid (factor 4 inside the band), the largest change of gamma between a coarse point and its neighbouring refined points beyond what the
     100-digit values themselves give. No threshold on gamma/M is introduced.
   - Growth is CONFIRMED if gamma > 10 u(k) at both precisions and gamma stays positive on at least three consecutive points of the refined grid:
     FAIL-UNSTABLE. Report gamma_hat_max, the band of k, u(k) and the e-folding time 1/gamma. (The factor 10 is a margin on the measured uncertainty.)
   - Apparent growth with gamma <= 10 u(k), or changing sign under refinement: INCONCLUSIVE.
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
     Calibration run (with a perturbation; a free system without perturbation conserves energy trivially and does not test the integrator):
       (a) the background gradients of the light fields are kept, to the extent that they also solve the free equations;
       (b) the background values of heavy and auxiliary fields are chosen so that they satisfy the FREE equations (they are not copied from the coupled
           background: with the coupling switched off, the coupled value of a sourced heavy field is in general no longer a solution);
       (c) the same profile, shape and relative amplitude of the initial perturbation as in the candidate run;
       (d) the same numerical method, time step and spatial grid, and then (2N, dt/2).
     Additional diagnostic, reported but not a pass criterion: drift of the fluctuation energy E_fl = E_total - E_background with the background part
     subtracted analytically. (The total energy is dominated by the background, so the total-energy drift alone is a weak test of the fluctuations.)
   - Convergence: all of the above repeated at (N, dt) and (2N, dt/2); results agree within 1 %.
P8 Backgrounds: the sets used for T2, plus the candidate's own stationary backgrounds; seeds fixed and recorded.

### Outcomes per background
- PASS: C0-C3 hold with the required evidence level (C4 on designated backgrounds).
- FAIL-ILLPOSED: C1 or C2 fails.
- FAIL-UNSTABLE: C1 and C2 hold, and a growing mode is confirmed in a finite band. Counts as a T3-fund failure.
- INCONCLUSIVE: C0 not satisfied; unresolved roots; no analytic argument for C1 (FUNDAMENTAL) or for a degenerate C2; C4 not converged.
Per candidate: the counts of the four outcomes on each background set, together with the declared regime. "Passes T3-fund on a set" is used only when
there is no FAIL and no INCONCLUSIVE on that set.

### Relation to the other tests
T2 overlaps with C2 and C3; T3-fund adds the regime declaration, the k -> infinity evidence level, completeness of eigenvectors, the separation of
ill-posed from unstable, the time-domain confirmation and the INCONCLUSIVE rules. T5 (positive energy) is sufficient but not necessary for C3.
T3-B4 keeps its historical values. T6 is independent.

### Regime of the existing candidates (decision of the project lead, 2026-11-13)
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

## Acceptance rule
No UV candidate is called accepted until T3-fund (in the registered form) and T6 are evaluated for it.
