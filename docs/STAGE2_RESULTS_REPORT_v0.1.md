# BEC Stage 2: results of the registered run — report v0.1

**Status:** report for review. Not registered. Criteria, thresholds and verdicts are those of v2 and annex v0.7. Nothing was changed after the result was seen.

## 1. Run identity and integrity

| item | value |
|---|---|
| code | registered script v0.4.3 (commit ad69257), sha256 3f5d5ea1…, as recorded in RUN_CONTEXT |
| machine | Lenovo, Windows 11; environment `stage2`: Python 3.13.18, numpy 2.5.5, scipy 1.18.1 |
| SciPy pin | trf.py LF-normalised 0c7c29d5… (pinned); raw a174a6e7…; CRLF |
| inputs | manifest 1bae892a…, anchor 99a6b04b…, lambda_eff e0729976…, bec_dispersion 228086f4… (all as registered) |
| time | 2026-11-20 00:28:18 – 01:37:18 UTC (69 min) |
| stderr | empty (no ABORT) |
| output integrity | all 6 entries of SHA256SUMS.txt verified on the received files |

**Hashes of the outputs:**
- targets.json: fa01175a…
- positive_control.json: e6db2798…
- per_target.json: de1795b6…
- sign_tuples.json: e78780a1…
- summary.json: 458fcac1…
- RUN_CONTEXT.txt: e5ca5bf3…

**Order of work (annex J5).** The run did C1c, then the positive control, then recorded its hash in RUN_CONTEXT, and only then evaluated the real targets.

**Targets.** 29 targets: the anchor plus 28 manifest rows. The 12 manifest rows not kept were all excluded for complex roots, i.e. the negative controls. There were no p·s = 0 cases, C1c passed, and no target is ill-conditioned.

`run_stdout.txt` was written by PowerShell redirection as UTF-16 LE. Its content is the console log. It is not covered by SHA256SUMS.

## 2. Positive control (annex J)

| region | recovered | generated | not generated | generated, not recovered | passed |
|---|---|---|---|---|---|
| A | 29/29 | 29 | 0 | 0 | yes |
| B | 29/29 | 29 | 0 | 0 | yes |

- Best cost on the synthetic targets: median ≈ 8·10⁻²⁷ in A and ≈ 3·10⁻²⁴ in B; at most 2·10⁻¹⁴.
- max ε = 0.000 throughout.

Per J4, a NO MATCH is interpretable **within the demonstrated sensitivity** in both regions.

## 3. Primary result: region A

**29/29 NO MATCH.** There are no INCONCLUSIVE, NO ADMISSIBLE FIT or NOT INTERPRETABLE verdicts.

| criterion | targets passing (of 29) |
|---|---|
| (1) sign | 8 |
| (2) magnitude, max ε ≤ 0.25 | **0** (max ε ranges from 1.011 to 1.979) |
| (3) inner roots, \|Δx\| ≤ 0.05 | 28 |

- **Best costs:** minimum 25.0, median 34.1. The positive control reached ≈ 10⁻²⁶.
- **Sensitivity over all 24 branch assignments:** 0 assignments pass (1)–(3) for any target.
- **Search outcomes over 29 × 80 starts:** 2236 converged, 84 not converged in the region, and no rejected or shifted-out starts. No start was without an admissible start point.
- **A1 control count:** 0 points in A had λ_min in [−10⁻⁹, 0] while not positive definite. So the class-A classification does not depend on the operational tolerance.

## 4. Control: region B (reported separately)

**29/29 NO MATCH.** Best costs: minimum 2.16 (target 13), median 30.0.

The closest case is target 13. It passes (1) and (3) and fails only (2), with max ε = 0.43. It is still no match. It is also in region B: energetically unstable in the lab frame, and not part of the primary result.

## 5. Structural finding from the sign-tuple sampling (annex D; empirical)

Both regions reached 2·10⁵ points (acceptance A 0.267, B 0.540).

- **B4 targets.** 11 distinct sign tuples of b̂, ordered by root value. In **every** B4 target, the branch at the **lowest** root has b̂ ≥ 0 (sign + or 0) and the branch at the **highest** root has b̂ ≤ 0 (sign − or 0).
- **BEC, region A.** 9 distinct tuples. In **every** one, the lowest-root branch has b̂ ≤ 0 and the highest-root branch has b̂ ≥ 0, i.e. the opposite sign at both outer branches. No B4 tuple was observed among the 2·10⁵ region-A points.
- **BEC, region B.** 15 tuples. One B4 tuple, (0, +, −, 0), appears 6 times in 2·10⁵. Three B4 targets carry it, target 13 among them.

**Reading (interpretation, not a registered criterion).**
- The β coefficients are the k² corrections to the characteristic speeds.
- In the stable BEC family, the outer characteristic branches bend **outwards**: the fastest right-moving and the fastest left-moving branch become faster at larger k. This is consistent with the positive quantum-pressure term k⁴/4.
- In the registered B4 backgrounds, the outer branches bend **inwards**, as the mediator corrections A_i give.
- The mismatch is therefore one of **sign of dispersion**, not of the characteristic roots: the roots themselves are reproduced in 28 of 29 targets.
- This is an observation from the sampling and the fits. An analytic statement, that the outer β signs of the stable BEC family are fixed, is a hand conjecture that remains to be proved.

## 6. Comparison with the stated expectation (annex §K)

- **Outcome: as expected.** NO MATCH in region A for all 29 targets.
- **Reason: not as expected.** §K expected the root structure (asymmetry and Landau limit) to be the main obstacle. It was not: criterion (3) passes for 28 of 29 targets in region A. The obstacle is criterion (2), the magnitude and sign pattern of β. §K's main structural reason is therefore **not supported** by the run. The expectation was right for a different reason than the one written.

## 7. Scope of the conclusion (as fixed in v2, J4 and J6)

**What the result says.** The rank-one B4 backgrounds are not reproduced, at the level of characteristic roots together with β signs and magnitudes, by the mean-field two-component BEC family with f ≤ 0.99, in the symmetric frame, without LHY, under the registered procedure. Its sensitivity was demonstrated on attainable synthetic targets.

**What it does not say:**
- that no match exists anywhere in the family;
- anything about other UV classes, beyond-mean-field corrections or other identifications;
- anything about operator-level matching, which a match would have required next.

## 8. Record note

The record of the aborted run attempt 1 was committed in ad69257 as `experiments/stage2_bec/run_attempts/attempt1_lenovo_v042/run_stderr.txt`.
- Git on Windows converted its CRLF line endings to LF at `git add`. The committed file is the LF form, sha256 eaf47bd5…, 506 bytes.
- It is not the uploaded CRLF form 3eaae0ab… (514 bytes) that note v0.2 lists.
- The content is otherwise byte-identical. The commit script's final blob check flagged this after the push.

## 9. Proposed next steps (for decision)

1. Register the results: the run outputs and this report, after review.
2. Optionally, an analytic note on the β-sign pattern of §5, to turn the empirical finding into a statement with stated premises. No new run would be needed.
