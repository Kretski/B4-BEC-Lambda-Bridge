# SM+η Stage 1: C1(b) full high-precision diagnostic (addendum)

Post-hoc and diagnostic only. No registered criterion, tolerance, population or C3 outcome is changed or recomputed; the 1429 Stage 1 results (smeta_stage1_result.json, sha256[:16] 39c77359831d231f) are used as they are. The registered check (b) in double precision remains recorded as violated in 215 cases; this addendum explains it, it does not overwrite it.

Correction to the Stage 1 bundle: `diag_c1.py` there is a 7-case sample for check (c), not a check (b) script, and the 50-digit recheck for (c) (`c1_mp_recheck.py`) covers only the grid points where the double value exceeded the bound, not the full grid.

## Method (`c1b_full.py`)
Cases: all 215 with registered max_dist > 1e-9 (210 Z0, 5 S1) + 5 random controls that passed (b) (seed 20261114). For k = 0.7, 3, 50: roots of det K (60-digit c-polynomial, mp.polyroots) × k compared in both directions with the eigenvalues of the 60-digit first-order reduction (ω = i·eig(ikP+C), nonzero modes); same comparison in double precision; cluster structure of the roots.

## Results
- 60-digit agreement: max multiset distance 3.3e-26 over all 220 cases and all three k (215 failing + 5 controls). The first-order reduction and det K are the same system; the registered failure is not an implementation inconsistency.
- 210 Z0 cases: det K has exactly double roots (cluster spread at the 60-digit √-precision floor, ≤1e-25) and the registered double-precision distance at k=0.7 (1.6e-8–4.2e-7) is of the order expected for double roots in double precision (√ε ≈ 1.5e-8 times scale). The double-precision discrepancy grows with k (observed roughly ∝ k²: e.g. 2.7e-8, 1.0e-6, 3.2e-4 at k=0.7, 3, 50 in the first case); only k=0.7 enters the flag.
- 5 S1 cases (G1 124 f=10; G1 146 f=1, 10; G1 165 f=10; G3 188 f=10): roots are simple but clusters of two to three roots lie within about 1e-3 of ω = ±k (near-luminal; min separation 7e-5–4.8e-4) and the double-precision first-order eigenproblem is moderately ill-conditioned (cond(V) 59–675). The discrepancies are 1.4e-9–1.4e-8. So these are near-degeneracy conditioning effects, not exact multiplicity (flag `sensitivity_explained` = False for them by construction); the conditioning explanation is supported by these numbers but is not a proof.
- Controls: all 5 pass at 60 digits (≤1e-26) and in double (≤4e-12).

## Status of C1
Check (b): double-precision violation in 215 cases, resolved as numerical (60-digit agreement 3.3e-26). Check (c): violation in 492 double-precision cases, 50-digit recheck at the flagged points (492/492 within the bound); the full-grid 50-digit check was not run. C1 is therefore recorded as: analytic well-posedness argument holds; registered double-precision consistency checks flagged 215 + 492 cases, all of which are numerically explained at high precision; not recorded as unconditionally passed. No C3 outcome depends on it.

Files: c1b_full.py, c1b_full.jsonl, c1b_s1_near_cluster.py, c1b_s1_near_cluster_output.txt.
