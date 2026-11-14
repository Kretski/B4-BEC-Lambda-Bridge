# T3-fund Proca: Stage 1 implementation and self-test, revision 3 (before any population evaluation)
Spec registered in commit e355af2 (SPEC_HASH 86c360326da7406f, MANIFEST_HASH 35a3c113d517f174). Nothing here was run on a manifest case (random backgrounds, seeds 20261115/20261116). Spec, thresholds, C3 grid, P6 criterion, Pred-4, populations, seeds and C4 rule are unchanged.

Files: `t3fund_proca_stage1.py`, `smeta_ref.py` (unchanged SM+η Stage 1 module: polynomial arithmetic, jet-Hessian determinant, grid, tolerances), `selftest_proca.py`, `selftest_proca_fail.py` + outputs, `selftest_proca_result.json`, and (added in this revision) `t3fund_proca_check.py` with `proca_check_symbolic_output.txt` — the symbolic check of the spec referred to below.

## Revision 2 corrections (kept): k=0 without rounding; C0 split into `abar_construction` and EL residuals; hidden-by-double; P6 completeness.
## Revision 3.1 (documentation only, after the final review)
Header comment of `t3fund_proca_stage1.py` updated (50-digit roots and classification at every grid point, not every 10th); docstring of `t3fund_proca_check.py` now states that only checks (1) C0 and (2) det K identity are executed there (checks (3)-(5) live in the spec bundle scripts; nonlinear C4 code does not exist yet). No code change; normal self-test repeated: worst deviations unchanged (Q 2.3e-51, det K 4.4e-16, 6x6 vs Q 3.7e-17, k=0 2.3e-17, C0 EL 0).

## Revision 3 corrections
1. Unified C3 classification by the 50-digit value, independent of the double-precision branch: r_m ≤ 1e-7 real; 1e-7 < r_m ≤ 1e-5 unresolved; r_m > 1e-5 complex. TOL_CONF (1e-12) no longer enters the classification (before, r_d>1e-5 with r_m between 1e-12 and 1e-7 was classified complex); it feeds only the diagnostic `real_with_tiny_imag` (real points with 1e-12<r_m≤1e-7). `hidden_by_double` counts points with r_d ≤ 1e-7 and r_m > 1e-7. New outputs `n_real`, `n_unresolved`.
2. A P6 PASS (or no complex point at all) can no longer mask an unevaluated or unresolved grid point: if c3 would be PASS and `nodata>0` or `unresolved_stuck>0`, C3 = INCONCLUSIVE. (The same rule also covers unresolved points, which the SM+η code handled only when no P6 verdict existed.) A P6 FAIL-UNSTABLE is not weakened by this rule.
3. P6 completeness now checks identity, not only the count: `complete` requires that the refinement points are unique, that the set of evaluated indices equals {0,…,len(fine)−1} (each exactly once, in order) and that every point was evaluated at both 50 and 100 digits; recorded `n_refined`, `n_unique_refined`, `n_rows`, `complete`.
4. Tests strengthened: hidden-by-double now requires exactly 92/92 hidden points, n_real=0, n_complex=92; new test: double screening forced complex on a real background → 92/92 'real', outcome PASS (classification follows r_m); new test: P6 forced to PASS on a complex window (dps-dependent imaginary part) with one unevaluated grid point → C3 INCONCLUSIVE (the old logic would have returned PASS); the vacuous assertion of the earlier draft was removed.
Note on C0: `EL_f, EL_g, EL_a0, EL_a1` are ∂L/∂(field value) on the constant background jets (all derivative terms of the EL equations vanish on constant jets); they use the same Abar and are not an independent check of the background. The independent symbolic check is `t3fund_proca_check.py` part (1): it builds the full Euler–Lagrange operators with `sympy.calculus.euler.euler_equations` for the four fields (including the t- and x-derivative terms), substitutes the background functions f=wt+qx, g=vt+rx, A0=0, A1=Ā, evaluates the derivatives (`.doit()`) and finds residuals [0,0,0,0] for general symbolic parameters (see output file).

## Self-test (independent SymPy EL operator, exact rationals), 8 random backgrounds × x∈{0.7,3,50}
Q coefficients 2.3e-51; det K/k⁸ (jet-Hessian) vs sympy 4.4e-16; ratio to Q vs −(1+x²)/x⁴ 9.1e-16; 6×6 vs k·roots(Q) 3.7e-17; k=0 vs closed form 2.3e-17; C0 EL residuals 0; 8/8 PASS (hyperbolic, μ-bound 0 violations, E true). FAIL path (4 random backgrounds): FAIL-UNSTABLE, P6 complete 365/365 unique, E false; ghost control FAIL-GHOST; unit tests (hidden 92/92; P6 100-digit failure → INCONCLUSIVE; classification by r_m; nodata not masked by P6 PASS) pass.
Not done (by instruction): nonlinear C4 code and its ≤1e-12 check against the SymPy EL operator (a later step; not claimed as passed); the 10 081 registered evaluations. Estimated cost ≈1 s per non-scope case, ≈5 s per scope case (1041), ≈2 h on 2 cores.
8966488f3de85510   selftest_proca.py
c3d6e0878daecf67   selftest_proca_fail.py
6f915e9b93c29e8e   smeta_ref.py
542e3fd4cb7969fc   t3fund_proca_check.py
7be1486a4fde2fe5   t3fund_proca_stage1.py
