# HeavyScalar_TwoField: population report

Status: exploratory. The official criteria T1-T6 (docs/UV_ACCEPTANCE_PROTOCOL.md) are not modified. This report states what was checked and what stays open.

## Model
L = 1/2 X_phi + 1/2 X_psi + 1/2 (d chi)^2 - 1/2 M^2 chi^2 + chi (a X_phi + b X_psi + e Y), X_f = (d f)^2, Y = d phi . d psi, signature (+,-), M = 1.
Background: phi0 = w t + q x, psi0 = v t + r x, chi_bar = (a X_phi + b X_psi + e Y)/M^2 (exact constant solution).

## Analytic statements
1. det K(c, k) = (c^2 - 1) * P4(c; k) exactly (symbolic division, remainder 0, asserted in the script). One luminal pair on every background.
2. Principal symbol of (phi, psi, chi): diag(Z, 1) (omega^2 - k^2), Z = [[1 + 2a chi_bar, e chi_bar], [e chi_bar, 1 + 2b chi_bar]]. Front velocity 1 iff Z > 0.
3. Quadratic energy in velocities: H2 = 1/2 qdot^T A qdot + 1/2 y^T (-H_yy) y with A = diag(Z, 1). The potential form at wavenumber k is
   [[k^2 Z, i k u], [-i k u^T, k^2 + M^2]], u = (2a q + e r, 2b r + e q). Schur complement: positive for all k iff Z - u u^T / M^2 >= 0 (worst case k -> 0).
   T5-static: Z > 0 and Z - u u^T/M^2 >= 0.
4. A positive definite conserved quadratic energy excludes growing modes, so T5-static implies real frequencies (T2). The converse does not hold.

## Numerical criteria (fixed before the runs)
T2: quartic roots real on X = k/M in [1e-3, 1e3], 61 log points, |Im c| <= 1e-7 max(1, |c|). T3: min eig Z > 1e-12.
T4-proxy: max |eig(Z) - 1| <= 0.3. This is a diagnostic, NOT the official T4 (g <= 1, rho <= 0.3), which is not evaluated here.
T5-static: as above, tolerance 1e-12. T5-Hamiltonian: number of negative eigenvalues of the 6x6 Hermitian H2(k) in (q, pi) at k = 0.01, 0.1, 1, 10.
T5-Krein: sign of the mode energy v^dag H2 v for each real-frequency normal mode at the same k.

## Results (500 backgrounds per coupling set, uniform in [-1, 1]^4)
| (a, b, e) | seed | T2 | T3 | T4-proxy | T5-static | all four | T2 but not T5-static | H2 positive definite | backgrounds with a negative-energy mode |
|---|---|---|---|---|---|---|---|---|---|
| (0.2, 0.2, 0.4) | 21 | 474 | 500 | 432 | 465 | 428 | 9 | 465 | 10 |
| (0.1, 0.1, 0.2) | 21 | 500 | 500 | 500 | 500 | 500 | 0 | 500 | 0 |
| (0.3, -0.2, 0.1) | 21 | 500 | 500 | 500 | 500 | 500 | 0 | 500 | 0 |
| (0.2, 0.2, 0.4) | 22 | 478 | 500 | 450 | 475 | 449 | 3 | 475 | 5 |
| (0.1, 0.1, 0.2) | 22 | 500 | 500 | 500 | 500 | 500 | 0 | 500 | 0 |
| (0.3, -0.2, 0.1) | 22 | 500 | 500 | 500 | 500 | 500 | 0 | 500 | 0 |

- T5-static but not T2: 0 in every set (max relative Im c 2.9e-12 when T5-static holds), consistent with statement 4.
- H2 is positive definite on exactly the T5-static backgrounds in every set; no negative-energy mode occurs inside T5-static.
- Backgrounds with T2 but not T5-static are dynamically stable and energetically unstable (the analogue of BEC region B).
- Seed 21 reproduces the first exploratory run (474 / 432 / 465 / 428). Seed 22 is an independent repetition.

Wording: for two of the three coupling sets all 500 randomly generated backgrounds pass the four criteria used; for the third, 428 (seed 21) and 449 (seed 22) of 500 pass. This is a result for these populations, not a global theorem.

## Correction
A first pass tested the degree-6 polynomial without factoring out (c^2 - 1) and reported 61 / 107 / 12 backgrounds as T2 failures although T3 and T5 held.
These were numerical artefacts: near-degenerate roots give |Im c| ~ 1e-7. They are not results.

## Open
- Official T4 for this model (needs a definition of g and rho for the couplings a, b, e).
- T5: Krein analysis is done only at four wavenumbers and in the lab frame; no boost dependence.
- T1 and T7 are negative for this model (docs/SCALAR_MEDIATOR_RESULT.md): it does not reproduce the B4 operator structure. T8 is undefined.
- Other coupling sets, |gradients| > 1, quantum corrections.

## Reproduce
python experiments/uv_candidates/heavy_scalar_T2_T5_population.py --seed 21
Outputs: results/heavy_scalar_population/seed21_*.csv (one row per background) and seed21_summary.json (seed, tolerances, counts).
