# B4.10 Experiment 3 — multi-background result (package v0.3.0)

Pre-registered manifest unchanged (`data/e3_manifest.json`, thresholds 0.10/0.30, quotas 8/6/6 per mode).
40 backgrounds + B3 anchor (reported separately, not in the quota). Grid: 241 log-spaced x in [1e-3,1e3] plus adaptive
refinement (reality transitions located to dlog 1e-7; overlap tracking refined to assigned overlap ≥0.98).
Run: `python run_experiment3.py results_e3`. Tests: `pytest` (68).

**One tolerance correction (disclosed):** the first run used a positivity threshold 1e-10·‖h‖, which wrongly flagged one
background at x≥631 (h_min/‖h‖=3.7e-11 > 0). Changed to 1e-13·‖h‖ (≈ eigensolver noise) before the reported results. No selection or
analysis criterion was changed.

## Outcome table (counts of backgrounds)
| group | n | local quartic real | full 6-root real at all x | Hamiltonian positive at all x | negative-energy branch | Krein sign change | composition exchange | exchange inside EFT window η=0.1 / 0.3 |
|---|---|---|---|---|---|---|---|---|
| rank1_minus stable | 8 | 8 | 8 | 7 | 1 | 0 | 8 | 0 / 0 |
| rank1_minus boundary | 6 | 6 | 5 (1 loses reality at x=0.1887) | 4 | 2 | 0 | 5 | 0 / 0 |
| rank1_minus negative ctrl | 6 | 0 | 0 (complex at every x) | 0 | – | – | – | – |
| rank1_plus stable | 8 | 8 | 8 | 8 | 0 | 0 | 8 | 0 / 1 |
| rank1_plus boundary | 6 | 6 | 6 | 4 | 2 | 0 | 6 | 0 / 0 |
| rank1_plus negative ctrl | 6 | 0 | 0 | 0 | – | – | – | – |
| B3 anchor | 1 | 1 | 1 | 1 | 0 | 0 | 1 (branches 0,1,3,5) | 0 / 0 |

## Findings
1. **Negative controls behave as controls**: all 12 have complex roots at every x (QZ and mpmath agree on reality at every checked node) and an indefinite h. The Proca completion does not stabilise them anywhere in 1e-3 ≤ x ≤ 1e3.
2. **Hamiltonian positivity and real roots are not the same thing.** h is never positive when roots are complex (0 inconsistencies), but 5 of 28 real-root stable/boundary backgrounds have real simple roots *and* an indefinite h with exactly one negative-energy slow scalar branch (Hamiltonian inertia = number of E<0 branches, verified at every node, 0 mismatches). Krein sign change along x: none in any background; the negative-energy branch exists already at x→0 (local EFT property, not created by the Proca sector).
3. **Closed-form criterion (exact, tested on 41/41 backgrounds):** at x→0 h is positive iff −C is positive definite, i.e. 1+3Λq²>Xw, 1>Xv, and (1+3Λq²)(1−Xv)>Xw. For Λ=0 this is Xw+Xv<1. The five negative-energy backgrounds have Xw+Xv between 1.15 and 1.89, but Xw+Xv>1 is not itself diagnostic: positive-h backgrounds reach 1.89 too because a large 3Λq² can compensate. The full criterion decides. (Positivity of h for all x was observed to follow the x→0 criterion in these 41 cases; not proven for all x.)
4. **The negative energy is the shift-symmetric one, δ²(H−μQ)** (naming corrected in v0.3.1: `energy.py` formerly called it the lab T00). For these five the lab-frame second-order T00 is positive for every branch and every sampled x. See `E3_ENERGY_NOTE.md` for the three separate results and the open physical question.
5. **Reality loss by Krein collision (1 case):** rank1_minus boundary #4487: real for x<0.18869, complex above (max relative Im 0.09). Just before the transition two adjacent simple roots with equal and opposite energies (±5.4e-8) merge — textbook Krein collision. The local quartic was real (rel_gap 0.068): the Proca completion creates a complex window here.
6. **Composition exchange is common but not universal in detail:** 27 of the 28 stable/boundary backgrounds (anchor excluded) show at least one crossing of the Proca energy share through 0.5; 101 crossings in total, x from 0.14 to 9.2 (median 0.65), |ω|/M median 1.66. Branches exchanged vary by background.
7. **Exchange vs EFT window:** only one crossing (rank1_plus stable #5461, branch 3, x=0.226, |ω|/M=0.26) lies inside the loose η=0.3 window; none inside η=0.1. Among the 28 real-spectrum stable/boundary backgrounds, the exchange lies outside the window for 27 at η=0.3 and for all 28 at η=0.1.
8. **Group velocities** (group, not front): every real-root background has max |v_g|>1 (min 1.37, median 2.4); the anchor 8.3. Superluminal group velocities persist in the full system — causality (front velocity) remains B4.11.
9. **Tracking**: no unresolved intervals; one background has a near-degenerate adjacent pair (relative gap 2.8e-5, the Krein-collision case) — expected.

## Limits
* Stratified diagnostic sample, not a random sample: no probability statements about backgrounds.
* One x-grid; B3 anchor provenance still unrecovered (generator not in T1.7.1–T1.7.3).
* EFT window is operational; g∝M_A, so UV perturbativity is not examined (B4.12). Front velocity not computed (B4.11).
