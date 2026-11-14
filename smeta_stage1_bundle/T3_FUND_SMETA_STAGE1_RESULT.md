# T3-fund, SM+η: Stage 1 result (C0–C3 + C1 consistency)

Registered in commit d006f77 (SPEC_HASH b2fbbf16e6741ac7, MANIFEST_HASH 9e9eb854437015e0). No criterion, population, tolerance or prediction was changed after seeing results. Scope: tree level, this Lagrangian, the sampled backgrounds only.

Artefacts: `t3fund_smeta_stage1.py` (sha256[:16] 6f915e9b93c29e8e), `smeta_stage1_result.json` (39c77359831d231f), manifest sha256[:16] 46e0b19d4d8f0db8. 1429 evaluations = registered population.

## 0. Execution history (disclosure)
The first launch (background, 2 CPUs) was lost after ~2 h when the container restarted; no result had been written or read. `evaluate()` was not touched. Only the driver `__main__` was changed to write one JSON line per case and resume from it; the run was then completed in foreground chunks. Final `smeta_stage1_result.json` is assembled from all 1429 cases in registered order.

## 1. Outcomes per population (C3 outcome)
| Population | n | FAIL-UNSTABLE | FAIL-GHOST | DEGENERATE | PASS |
|---|---|---|---|---|---|
| S1 G1 (f=0.1, 1, 10) | 138 each | 138 each | 0 | 0 | 0 |
| S1 G2 (f=0.1, 1, 10) | 133 each | 133 each | 0 | 0 | 0 |
| S1 G3 (f=0.1, 1, 10) | 117 each | 117 each | 0 | 0 | 0 |
| Z0 (branch Z) G1/G2/G3 | 62/67/81 | all | 0 | 0 | 0 |
| H3 (f<0 control) | 40 | 0 | 40 | 0 | 0 |
| S2g | 2 | 0 | 2 | 0 | 0 |
| H3b (f=0 boundary) | 10 | 0 | 0 | 10 | 0 |
| Anchors (A_G1..A_G3) | 3 | 3 | 0 | 0 | 0 |

No S1 background (0 of 1164 cases, 388 backgrounds × 3 values of f) is C3 PASS. 1110 of 1164 S1 cases are complex already at k/M=1e-3; 1068 have complex k=0 roots. No case has a confirmed complex root in k/M ≥ 1e3 (UV band): the instabilities are at k/M<1e3. Largest grid growth rate in S1: 5.38 (units of M).

## 2. Predictions
- **Pred-1 (ghost): confirmed.** Classification follows the §3(ii) priority rule in all 1429 cases (0 mismatches); H3 40/40 FAIL-GHOST, S2g 2/2 FAIL-GHOST, H3b 10/10 DEGENERATE. Sign of the c⁸ coefficient (50 digits, k=1e3) equals sign(Z_tφ Z_tψ f) in all 1419 non-degenerate cases (0 mismatches; `pred1_c8.py`).
- **Pred-2 (k=0): confirmed.** Numerical k=0 class equals the closed form in 1164/1164 S1 cases; max relative error of the roots in x: 4.95e-11 (< 1e-9).
- **Pred-3 (UV speeds): confirmed.** 693 S1 cases have real roots at k=1e6; max deviation from {±1, ±1, c_i±}: 3.86e-5 (< 1e-4).
- **Pred-4 (no non-IR instability when the closed-form k=0 class is real): FALSIFIED.** There are 96 such cases (32 G3 draws × 3 f). All 96 are C3 FAIL-UNSTABLE: 8–43 of the 92 grid points complex (median 24.5); 50 of them already complex at k/M=1e-3, 46 first at larger k; none in the UV band k/M ≥ 1e3. So the real k=0 class does not imply stability at finite k; the instability is a finite-k band, not a UV one. Reported without changing any criterion.
- Diagnostic (not an outcome): the 693 S1 cases with real UV roots all have front speed max|c| > 1+1e-6 (consistent with the manifest: all S1 backgrounds are superluminal in the front sense). Under preferred-foliation causality this is allowed and not a verdict.

## 3. C0, C1 checks
- C0: max residual of the σ and η equations at 50 digits 2.1e-50 (all 1429).
- C1 first-order reduction (1377 cases that are not ghost/degenerate): A0 invertible; max |Im λ(P)| = 0 and coefficients real; cond(V) ≤ 2.17.
- C1 consistency check (b), eigenvalues of the first-order reduction vs roots of det K at k=0.7 to 1e-9: violated in 215 cases (210 Z0, 5 S1; max 4.2e-7). Likely cause (not separately verified): σ and η share the luminal speed 1, so double roots appear and double-precision root error of order √ε≈1e-8 is expected; the 50-digit recheck below was done for the check (c) only.
- C1 consistency check (c), max_k γ(k) ≤ γ_max(1+1e-9) in double precision: violated in 492 cases (281 S1, 210 Z0, 1 anchor), all at k/M ≳ 6e5 where double-precision roots of a degree-8 polynomial with near-degenerate roots are unreliable. By the registered rule such a violation means implementation inconsistency, not FAIL. Supplementary 50-digit recheck (post hoc, diagnostic only, `c1_mp_recheck.py`) at the points where the double value exceeded the bound: 492/492 satisfy the bound, max γ_mp/γ_max = 2.5e-12. The strict registered check in double therefore carries a numerical-precision flag on these 492 cases; this does not affect any C3 outcome (all are FAIL-UNSTABLE from the 50-digit-confirmed C3 criterion, with no unresolved or no-data points).

## 4. C4
The two designated backgrounds (G3 idx 6 and idx 7, f=1) are both C3 FAIL-UNSTABLE (idx 6: 19/92 grid points complex, max γ 0.36; idx 7: 29/92 complex, max γ 8.7e-3, complex at k/M=1e-3). By the registered rule C4 is **NOT-RUN**, and no replacement background is chosen. C4 for this candidate is therefore not available.

## 5. Conclusion (limited)
Within this model at tree level and for the registered samples, SM+η with a dynamical η has no background that reaches C3 PASS. Ghost cases (f<0, S2g) are excluded as registered, branch Z is unstable at k=0, and branch S is unstable either at k=0 (complex closed form) or at finite k (the 96 real-k=0 cases). This does not exclude other couplings, other populations or other mediators, and it does not establish anything about a UV completion of B4. Per the registered order, T7 for SM+η is not started (no population reached C3/C4 PASS).
