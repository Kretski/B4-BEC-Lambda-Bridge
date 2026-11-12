# T3-fund protocol: PROPOSAL (not registered)

Status: proposal for the project lead. Nothing here is in force until it is confirmed and copied into UV_ACCEPTANCE_PROTOCOL_ADDENDUM_T3_T6.md.
No candidate has been evaluated against it. Criteria and tolerances must be fixed before the first run and not changed afterwards.

## What T3-fund tests
Under the registered convention (preferred-foliation causality) T3-fund asks whether linear evolution in the preferred time T is well posed and stable.
It places no bound on propagation speed. All quantities are evaluated in the frame of the observer field (the "lab" frame of the models).
Notation: for a stationary background and real wavenumber k the linearised system has frequencies omega_j(k); gamma(k) = max_j Im omega_j(k) is the
largest growth rate; M is the heavy scale of the candidate.

## Criteria
C0 Background. The background must be an exact solution of the full equations (stated analytically, or residual <= 1e-10 numerically).
   Otherwise: INCONCLUSIVE.
C1 Well-posedness (growth rate bounded in k). sup over real k of gamma(k) must be finite. An unbounded gamma(k) (for example gamma ~ k or k^2) means that
   arbitrarily short wavelengths grow arbitrarily fast: no continuous dependence on initial data, in any frame.
C2 Principal symbol. The highest-derivative part must have real characteristic roots and a complete set of eigenvectors (strong hyperbolicity for
   second-order systems; for higher-order dispersive systems: real leading-order frequencies). Real but defective roots (weak hyperbolicity) do not
   guarantee C1, because lower-order terms can then produce growth increasing with k.
C3 Mode stability. gamma(k) <= tolerance for all k: no exponentially growing mode at any wavelength. A background can satisfy C1 and C2 and fail C3
   (well posed but unstable in a finite band); this is reported as a separate outcome.
C4 Time-domain confirmation. For at least one background per candidate that passes C1-C3: full nonlinear evolution of a small perturbation, with bounded
   norm, conserved total energy, and unchanged result under refinement of the grid and the time step.

Examples that show what C1 separates (analytic):
- omega^2 = c^2 k^2 (1 + Lambda k^2), Lambda > 0: real for all k. Passes C1-C3 although the front velocity is unbounded.
- Same law with Lambda < 0: gamma = c k sqrt(|Lambda| k^2 - 1) for k > 1/sqrt(|Lambda|), growing like k^2. Fails C1 unless the theory is cut off below that k.
So T3-fund is not empty: it accepts the superluminal sign and rejects the subluminal sign of a bare k^4 term taken at face value.

## Numerical procedure and tolerances (proposed)
P1 Polynomial preparation. Exact factors of the characteristic polynomial (for example a luminal pair) are divided out symbolically before any root test.
   (Reason: near-degenerate roots of the unfactored polynomial gave spurious Im ~ 1e-7 in the heavy-scalar population test.)
P2 Grid. k/M on 61 logarithmic points in [1e-3, 1e3] (same as T2), plus an asymptotic segment of 31 logarithmic points in [1e3, 1e6].
P3 Reality tolerance. A root is real if |Im c| <= 1e-7 max(1, |c|), c = omega/k. Roots with 1e-7 < |Im c|/max(1,|c|) <= 1e-5 are "unresolved".
P4 C1 test. If gamma exceeds the tolerance in the asymptotic segment, fit gamma ~ k^s over [1e4, 1e6]: s > 0.05 fails C1.
   Where the principal symbol is known analytically, the analytic statement replaces the fit.
P5 C2 test. Roots of the principal symbol real within P3; eigenvector matrix of the principal symbol with condition number <= 1e8 (otherwise "defective").
P6 C4 test. Perturbation amplitude 1e-3 of the background scale, duration >= 1000/M; fluctuation norm at the end <= 3 times its initial value;
   relative drift of the total energy <= 1e-9; the same within 1 % when the time step is halved and the number of grid points doubled.
P7 Backgrounds. The same sets as for T2 (anchor + manifest, and the population cells), plus the candidate's own stationary backgrounds where applicable.
   Seeds are fixed and recorded.

## Outcomes per background
- PASS: C0-C3 hold (C4 on the designated backgrounds).
- FAIL-ILLPOSED: C1 or C2 fails (complex principal roots, or gamma growing with k).
- FAIL-UNSTABLE: C1 and C2 hold but gamma > tolerance in a finite band of k. Report the maximal gamma/M and the band.
- INCONCLUSIVE: C0 not satisfied; unresolved roots (P3); defective principal symbol without an analytic argument; C4 not converged.
Per candidate: report the counts of the four outcomes on each background set. A candidate "passes T3-fund on a set" only if there is no FAIL and no
INCONCLUSIVE on that set; otherwise the counts are reported without a verdict word.

## Relation to the other tests
- T2 (real spectrum, hyperbolic principal symbol on scanned backgrounds) overlaps with C2 and C3. T3-fund adds: the k -> infinity behaviour (C1),
  completeness of eigenvectors (C2), the separation of ill-posed from unstable, the time-domain confirmation (C4), and explicit INCONCLUSIVE rules.
- T5 (positive energy) is sufficient for C3 (a positive definite conserved quadratic energy excludes growth) but not necessary.
- T3-B4 (front <= c_B4) is unrelated and keeps its existing values.
- T6 is independent: T3-fund presupposes an observer field and says nothing about its origin.

## Open points for the project lead
1. Are the tolerances in P2-P6 acceptable (in particular the 1e-7 reality tolerance and the factor 3 on the norm)?
2. Should FAIL-UNSTABLE count as failing T3-fund, or be recorded under T2/T5 only? (Proposed: it fails T3-fund item 5, "no physically unacceptable instability".)
3. For theories that are only effective (valid below a cutoff k_max), should C1 be tested up to k_max instead of k -> infinity? If yes, k_max must be
   stated per candidate before the run.
