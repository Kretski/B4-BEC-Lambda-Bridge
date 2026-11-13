# T3-fund HS2 (HeavyScalar, two light fields) — C0–C3 result
Spec registered in commit 7b0e867 (SPEC_HASH 9d52c927a9acd3bd, MANIFEST_HASH 6c312456ac894009). Implementation: t3fund_hs2_stage1.py (no change of criteria). Predictions evaluated by t3fund_hs2_preds.py. Run: 954 backgrounds, 250 s.

## C0–C2 (numerical confirmation of the registered analytics)
- C0: |M²χ̄ − O| ≤ 5.6e-17 on all backgrounds.
- det K (explicit 3×3 matrix) vs k⁶·D·Q: max relative difference 6.0e-13 (954 backgrounds, 6 sample points each).
- First-order matrices (9×9): A real symmetric, cond(eigenvectors of A) = 1 to 1e-15; every root of D·Q lies within 3.3e-14 of an eigenvalue of kA+iB; λmax((B+Bᵀ)/2) ≤ 65.5 (largest on a H3 control background), finite everywhere.
- Hence C1 (finite growth bound for all k) and C2 (real symmetric A, complete eigenvectors) hold on every background with Z̄ > 0.

## C3 outcomes
| Set | n | PASS | FAIL-UNSTABLE | FAIL-GHOST |
|---|---|---|---|---|
| H1 K1 (rank-one) | 296 | 241 | 55 | – |
| H1 K3 | 300 | 282 | 18 | – |
| H1 K4 | 300 | 273 | 27 | – |
| H1 total | 896 | 796 | 100 | – |
| H2 (Z̄ not > 0) | 4 | – | – | 4 |
| H3 control | 50 | – | 50 | – |
| Anchors B2h, B2k, B2s, B2z | 4 | 4 | – | – |
No INCONCLUSIVE, no UNRESOLVED, no "hidden by double" case. P6: max u(k) = 2.7e-49, min γ/u = 5.1e46, every FAIL-UNSTABLE has ≥3 consecutive refined points above 10u.

## Registered predictions
- Pred-1 (FAIL-UNSTABLE ⇔ F has complex roots): 896/896 agree, 0 exceptions (K1 296/296, K3 300/300, K4 300/300).
- Pred-2: every unstable H1 background is already complex at k/M = 1e-3 (0 exceptions); no H1 background has a complex root in k/M ∈ [1e3,1e6] (0).
- Pred-3: H3 50/50 FAIL-UNSTABLE; H2 4/4 FAIL-GHOST; anchors 4/4 PASS.
- Pred-4 (K1 vs one-field HS, g=2a): MAX ROOT DISTANCE 1.5e-12 in double precision, i.e. slightly above the registered 1e-12. The worst case (K1 draw, k/M=30) recomputed in 50 digits gives distance 0.0; the excess is double-precision conditioning of near-degenerate roots, not a model difference. Reported as a nominal miss of the numerical threshold with this explanation; the criterion is not changed.

## Findings
- All 100 H1 instabilities are IR (k→0 local polynomial F has complex roots); the UV region is real for all 896 Z̄>0 backgrounds, as predicted by the k→∞ expansion ω = ±k ± √(s±ᵀZ̄⁻¹s±). With two light fields the Stage 1 conclusion holds again: UV completion adds no instability of its own on the sampled sets; instability fractions: K1 18.6 %, K3 6.0 %, K4 9.0 %.
- Designated C4 backgrounds (manifest, K3 idx 0 and K4 idx 0) are C3 PASS, so C4 uses them as registered.
- Not established: C4 for HS2, T7 (T7-a/T7-b not yet run), behaviour outside sampled ranges, UV completion of B4.
