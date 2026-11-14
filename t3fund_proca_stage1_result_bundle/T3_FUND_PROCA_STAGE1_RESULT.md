# T3-fund Proca Stage 1: population results (registered evaluation, no criterion changed)

Spec SPEC_HASH 86c360326da7406f, MANIFEST_HASH 35a3c113d517f174 (manifest file hash fd428e2d31195940); implementation committed in 0649531 (`t3fund_proca_stage1.py` sha256[:16] 7be1486a4fde2fe5), approved after review before any run. All 10 081 registered cases evaluated (PE 41, PN 10 000, H3 40); no case INCONCLUSIVE, no nodata, no unresolved or stuck grid point, no hidden-by-double, no real_with_tiny_imag.
Scope: tree level, this Lagrangian, linear C0-C3 only. C4 (time domain) was NOT run; T7 not started; B4.12 perturbativity not evaluated and not resolved by anything here.

## 0. Erratum (aggregation identifier, found in review before any commit)
The driver `t3fund_proca_stage1.py` (committed in 0649531, unchanged here) keys cases by `label`. The label `E3_stable_6052` occurs twice in the manifest (rank1_minus and rank1_plus), so the aggregated JSON written by the driver duplicated one and omitted the other, and a first draft of the cell table was wrong for the two PE cells. The raw JSONL was correct and complete (10 081 distinct cases, both collided cases evaluated). `aggregate_from_jsonl.py` regenerates `proca_stage1_result.json`, `cell_table.md` and `aggregate_stats.json` from the JSONL with the unique key (pop, cut, mode, label) and asserts that the 10 081 identities equal the manifest's. No case was recomputed; all totals and all prediction statistics are unchanged; only the split of PE between the two modes changed (14/7 and 14/6, not 14/7 and 15/6 as in the draft). The driver itself is not modified, so its hash stays 7be1486a4fde2fe5; its aggregated output file is replaced by the regenerated one.

## 1. Outcomes per population and cell
Total: PASS 9901, FAIL-UNSTABLE 140, FAIL-GHOST 40 (all H3), INCONCLUSIVE 0.

| population | cut | mode | N | PASS | FAIL-UNSTABLE | FAIL-GHOST | INCONCLUSIVE |
|---|---|---|---|---|---|---|---|
| PE | - | rank1_minus | 21 | 14 | 7 | 0 | 0 |
| PE | - | rank1_plus | 20 | 14 | 6 | 0 | 0 |
| PN | 0.01 | rank1_minus | 1000 | 1000 | 0 | 0 | 0 |
| PN | 0.01 | rank1_plus | 1000 | 1000 | 0 | 0 | 0 |
| PN | 0.04 | rank1_minus | 1000 | 1000 | 0 | 0 | 0 |
| PN | 0.04 | rank1_plus | 1000 | 1000 | 0 | 0 | 0 |
| PN | 0.09 | rank1_minus | 1000 | 1000 | 0 | 0 | 0 |
| PN | 0.09 | rank1_plus | 1000 | 1000 | 0 | 0 | 0 |
| PN | 0.25 | rank1_minus | 1000 | 1000 | 0 | 0 | 0 |
| PN | 0.25 | rank1_plus | 1000 | 1000 | 0 | 0 | 0 |
| PN | 0.95 | rank1_minus | 1000 | 959 | 41 | 0 | 0 |
| PN | 0.95 | rank1_plus | 1000 | 914 | 86 | 0 | 0 |
| H3 | 0.09 | rank1_minus | 40 | 0 | 0 | 40 | 0 |

PE: 28 PASS / 13 FAIL-UNSTABLE (rank1_minus 14/7 of 21, rank1_plus 14/6 of 20). PN: 9873 PASS / 127 FAIL-UNSTABLE; all failures are in the cut 0.95 cells (41/1000 rank1_minus, 86/1000 rank1_plus); cuts 0.01-0.25 are 6000/6000 PASS.

## 2. Failure mechanism (140 FAIL-UNSTABLE)
Every failure is non-strictly-hyperbolic in the UV (complex roots of Q(c;1), max|Im c| 0.009-0.47) and has complex grid points (49-92 of 92). In 115 the local EFT four-root polynomial F(c) already has complex roots (|Im c|>1e-3) and the first grid point k/M=1e-3 is complex (IR failure); in 25 F is real and the first point is real, the instability appears at higher k/M. Every PASS (9901) has 92/92 real points and strictly hyperbolic UV.

## 3. Registered predictions
- Pred-1: Hessian route vs Q max relative deviation 2.8e-12 on all 1041 scope cases (threshold 1e-9); 6x6 matrix vs Q 4.3e-14 on all 10 041 non-ghost cases. HOLDS.
- Pred-2: k=0 agrees with the closed form omega^2=M^2(1+pq^2+sr^2) in 10 041/10 041 (max rel. dev. 2.1e-16). HOLDS.
- Pred-3: UV speeds vs Q(c;1) max error 4.6e-11 on the 9901 real cases (threshold 1e-4). HOLDS.
- Pred-4 (scope 1041 = PE + first 100 PN per cell): E true in 1014, all PASS; E false in 27, all FAIL-UNSTABLE; counterexamples to "E => PASS" or "FAIL => not E": 0. HOLDS. (The converse was not predicted.)
- Pred-5: H3 40/40 FAIL-GHOST; kappa_kin=(1,1,1) in all 10 041 PE/PN cases; no GHOST/DEGENERATE. HOLDS.
- Pred-6: 115 cases with |Im c|>1e-3 for F: 115/115 first point complex and FAIL-UNSTABLE. HOLDS. (25 further failures have real F roots; not covered by the prediction.)
- C0: all residuals exactly 0. C1: mu-bound (50 digits, 92 points) defined in 1014 scope cases, 0 violations; undefined in 27 (not strictly hyperbolic, all FAIL).

## 4. UV front (max |c|) and C4 designation (rule of spec section 6, applied mechanically to the manifest order)
PN PASS: front min/median/max 1.04/1.30/7.95, 7899 of 9873 have front <= 2.0. Designated: PN_cut0.01_rank1_minus_0, PN_cut0.01_rank1_minus_1 and PE E3_stable_9632 (rank1_minus; first C3-PASS PE with front <= 2.0). B3 anchor is outside the C4 domain (front > 2). No replacement of designated backgrounds.

## 5. What this does and does not show
Within this model, a large part of the registered parameter region has an exact background and a linear spectrum that is real and strictly hyperbolic at all 92 grid points (C3 PASS), unlike the SM+eta result. This is linear mode stability on constant-gradient backgrounds only: it does not show nonlinear well-posedness or time-domain stability (C4), does not establish a UV completion of the B4 Lambda-effect (T7), and does not address the B4.12 perturbativity failure (0/41 for g<=1). The mediator-elimination idea remains a separate theoretical direction.

## 6. Next steps (not done)
Nonlinear C4 code with <=1e-12 check against the independent SymPy EL operator, then C4 on the three designated backgrounds; T7 for Proca only in a separate registered document after C4.
