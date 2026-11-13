# P7 calibration on the free system (T3-fund, C4)

Registered protocol: UV_ACCEPTANCE_PROTOCOL_ADDENDUM_T3_T6.md (P7, items a-d, Clarification 1). Script `experiments/uv_candidates/p7_free_calibration.py`
(specification block frozen before the first run, SPEC_HASH abdf4cc6180cc8ac; one syntax-only edit of two dictionary key names, no change of numerics).
Third resolution run by `p7_run_fine.py`, which calls the same frozen functions. Date 2026-11-13. Raw data in `results/p7_calibration/`.

## Setup (all couplings zero)
Periodic box L=1024, pseudo-spectral D, RK4, T=1000 (=1000/M, M=1), 201 records. Background: free packets in phi and psi (A=1, s=15, k0=0.3), chi_bg=0 (solves the free heavy equation).
Perturbation: Gaussian packets in phi, psi (s=15, k0=0.5), amplitude 1e-3 x max|d_x phi_bg|, heavy field unperturbed. Fluctuation = total run - background-only run.

## Result (protocol run)
| (N, dt) | max drift of total energy | drift of E_fl (diagnostic) | max norm / initial norm | target 1e-9 |
|---|---|---|---|---|
| (8192, 0.05) | 5.64e-9 | 2.6e-2 | 1.0 | not reached |
| (16384, 0.025) | 1.76e-10 | 8.0e-4 | 1.0 | reached |
| (32768, 0.0125) | 5.50e-12 | 2.5e-5 | 1.0 | reached |

Drift ratios between successive resolutions: 32.0 and 32.0 (RK4 energy dissipation of the resolved modes; decreases under refinement as required).
Fluctuation norm: final/initial = 0.99999996, 0.9999999988, 0.99999999996; max ratio 1.0 at all resolutions (agreement far better than 1 %).
Extra run with a perturbed chi packet (NOT part of the protocol, exercises the heavy-field integrator): drifts 5.64e-9 and 1.76e-10 (same ratio 32.0); fluctuation norm final/initial 0.9069 at both resolutions, i.e. the (E_fl-based) heavy-field behaviour is converged too.

## Decision of the project lead (reviewer text of 2026-11-13) and consequences
- Candidate base resolution: (16384, 0.025). Convergence check: (32768, 0.0125). (8192, 0.05) is not used for candidates because the free system does not reach 1e-9 there.
- The free system reaches the target at both resolutions to be used. The target 1e-9 is therefore attainable and C4 is not declared INCONCLUSIVE on this ground.
- E_fl stays a diagnostic. At the base resolution 1.76e-10 of the total energy corresponds to 8e-4 of the perturbation energy: the total-energy criterion is a weak test of the fluctuations, as the protocol notes. No criterion is changed.
- Wall time: about 15 min per resolution pair at base, about 70 min at the convergence resolution for the protocol pass alone (2 CPUs).

Status: calibration closed. Evaluation of candidates may start.
