# E3 manifest (for review — no Proca sweep has been run)

## Provenance
* Population reproduces the T1.7.2 scan: seed 172026, 10000 backgrounds, cut 0.95 (`generate_samples`, hash in SHA256_SOURCES.txt). Reproduction check: non-hyperbolic counts 286 (rank1_minus) and 661 (rank1_plus) equal the T1.7.2 `results/t172_robustness_cut_0.95.json` summary.
* **`eft_control` protocol definition (recovered from T1.7.2/T1.7.3 code):** `max(Xq,Xr,Xw,Xv) < cut`, cuts 0.5 / 0.7 / 0.95. B3 row: max X = 0.9106 → true only at 0.95, false at 0.5 and 0.7.
* **Still unrecovered:** the script that wrote the B3 row (`eft_domain`, `temporal_sum`, `max_abs_c`, row index 8899). It is not in T1.7.1/T1.7.2/T1.7.3, and the B3 row is NOT in the seed-172026 population (p searched over all 10000 backgrounds, also seeds 17303 and others). The B3 background is therefore kept as a separate anchor, not as a member of the E3 sample.

## Pre-registered rules (fixed before any E3 result)
* negative_control: max|Im c| > 1e-6; boundary: four real roots and rel_gap < 0.1; stable: four real roots and rel_gap >= 0.3; the band in between is a counted, unselected buffer. rel_gap = min root separation / max|root| of the local quartic.
* Per mode: 8 stable, 6 boundary, 6 negative controls (20); modes rank1_minus and rank1_plus; selection = evenly spaced quantiles in Xw+Xv within each class (deterministic).

## Class sizes in the population
| mode | stable | boundary | buffer | negative |
|---|---|---|---|---|
| rank1_minus | 3779 | 690 | 5245 | 286 |
| rank1_plus | 4729 | 739 | 3871 | 661 |

## Selected backgrounds (40)
| mode | class | bg index | p | s | Xq | Xr | Xw | Xv | Xw+Xv | rel_gap | max|Im c| | eft_ctrl 0.5/0.7/0.95 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| rank1_minus | stable | 6052 | 0.461 | 0.534 | 0.48 | 0.62 | 0.00 | 0.02 | 0.02 | 0.503 | 0.000 | 0/1/1 |
| rank1_minus | stable | 9632 | 0.065 | 0.116 | 0.06 | 0.75 | 0.38 | 0.06 | 0.44 | 0.345 | 0.000 | 0/0/1 |
| rank1_minus | stable | 1487 | 0.187 | 0.195 | 0.25 | 0.27 | 0.42 | 0.22 | 0.64 | 0.341 | 0.000 | 1/1/1 |
| rank1_minus | stable | 9582 | 0.099 | 0.499 | 0.76 | 0.26 | 0.77 | 0.04 | 0.81 | 0.337 | 0.000 | 0/0/1 |
| rank1_minus | stable | 5883 | 0.112 | 0.070 | 0.14 | 0.76 | 0.93 | 0.02 | 0.95 | 0.382 | 0.000 | 0/0/1 |
| rank1_minus | stable | 1506 | 0.466 | 0.127 | 0.84 | 0.03 | 0.53 | 0.57 | 1.10 | 0.349 | 0.000 | 0/0/1 |
| rank1_minus | stable | 7761 | 0.051 | 0.052 | 0.48 | 0.14 | 0.42 | 0.87 | 1.30 | 0.480 | 0.000 | 0/0/1 |
| rank1_minus | stable | 5713 | 0.602 | 0.124 | 0.94 | 0.87 | 0.94 | 0.95 | 1.89 | 0.437 | 0.000 | 0/0/1 |
| rank1_minus | boundary | 9745 | 0.339 | 0.342 | 0.31 | 0.02 | 0.03 | 0.01 | 0.03 | 0.082 | 0.000 | 1/1/1 |
| rank1_minus | boundary | 911 | 0.490 | 0.953 | 0.09 | 0.27 | 0.38 | 0.28 | 0.66 | 0.053 | 0.000 | 1/1/1 |
| rank1_minus | boundary | 2731 | 0.660 | 0.530 | 0.62 | 0.18 | 0.15 | 0.73 | 0.88 | 0.040 | 0.000 | 0/0/1 |
| rank1_minus | boundary | 5611 | 0.292 | 0.352 | 0.15 | 0.25 | 0.84 | 0.23 | 1.07 | 0.097 | 0.000 | 0/0/1 |
| rank1_minus | boundary | 4487 | 0.548 | 0.345 | 0.77 | 0.08 | 0.58 | 0.73 | 1.31 | 0.068 | 0.000 | 0/0/1 |
| rank1_minus | boundary | 2020 | 0.299 | 0.755 | 0.44 | 0.14 | 0.90 | 0.95 | 1.85 | 0.036 | 0.000 | 0/0/1 |
| rank1_minus | negative_control | 5433 | 0.393 | 0.113 | 0.01 | 0.06 | 0.88 | 0.18 | 1.06 | 0.362 | 0.192 | 0/0/1 |
| rank1_minus | negative_control | 1579 | 0.492 | 0.407 | 0.30 | 0.57 | 0.86 | 0.44 | 1.30 | 0.757 | 0.393 | 0/0/1 |
| rank1_minus | negative_control | 4026 | 0.954 | 0.639 | 0.49 | 0.82 | 0.74 | 0.68 | 1.42 | 0.363 | 0.235 | 0/0/1 |
| rank1_minus | negative_control | 6109 | 0.086 | 0.122 | 0.03 | 0.12 | 0.88 | 0.67 | 1.55 | 0.683 | 0.396 | 0/0/1 |
| rank1_minus | negative_control | 1945 | 0.927 | 0.979 | 0.88 | 0.25 | 0.74 | 0.95 | 1.69 | 0.330 | 0.354 | 0/0/1 |
| rank1_minus | negative_control | 6392 | 0.423 | 0.482 | 0.14 | 0.47 | 0.94 | 0.93 | 1.88 | 0.574 | 0.377 | 0/0/1 |
| rank1_plus | stable | 6052 | 0.461 | 0.534 | 0.48 | 0.62 | 0.00 | 0.02 | 0.02 | 0.401 | 0.000 | 0/1/1 |
| rank1_plus | stable | 522 | 0.177 | 0.100 | 0.43 | 0.57 | 0.19 | 0.26 | 0.45 | 0.355 | 0.000 | 0/1/1 |
| rank1_plus | stable | 8472 | 0.303 | 0.122 | 0.12 | 0.68 | 0.19 | 0.44 | 0.63 | 0.308 | 0.000 | 0/1/1 |
| rank1_plus | stable | 9056 | 0.440 | 0.898 | 0.11 | 0.04 | 0.34 | 0.44 | 0.78 | 0.490 | 0.000 | 1/1/1 |
| rank1_plus | stable | 5461 | 0.061 | 0.151 | 0.77 | 0.76 | 0.84 | 0.08 | 0.92 | 0.410 | 0.000 | 0/0/1 |
| rank1_plus | stable | 8325 | 0.253 | 0.534 | 0.12 | 0.16 | 0.42 | 0.65 | 1.06 | 0.397 | 0.000 | 0/1/1 |
| rank1_plus | stable | 7151 | 0.133 | 0.891 | 0.50 | 0.08 | 0.61 | 0.65 | 1.26 | 0.489 | 0.000 | 0/1/1 |
| rank1_plus | stable | 6465 | 0.058 | 0.448 | 0.72 | 0.34 | 0.95 | 0.94 | 1.89 | 0.424 | 0.000 | 0/0/1 |
| rank1_plus | boundary | 523 | 0.427 | 0.254 | 0.46 | 0.04 | 0.04 | 0.01 | 0.05 | 0.078 | 0.000 | 1/1/1 |
| rank1_plus | boundary | 7635 | 0.060 | 0.277 | 0.12 | 0.48 | 0.02 | 0.70 | 0.72 | 0.089 | 0.000 | 0/1/1 |
| rank1_plus | boundary | 825 | 0.081 | 0.759 | 0.09 | 0.88 | 0.14 | 0.84 | 0.98 | 0.095 | 0.000 | 0/0/1 |
| rank1_plus | boundary | 7900 | 0.925 | 0.234 | 0.92 | 0.35 | 0.69 | 0.46 | 1.15 | 0.012 | 0.000 | 0/0/1 |
| rank1_plus | boundary | 329 | 0.278 | 0.227 | 0.63 | 0.73 | 0.70 | 0.66 | 1.37 | 0.080 | 0.000 | 0/0/1 |
| rank1_plus | boundary | 3437 | 0.260 | 0.057 | 0.81 | 0.88 | 0.95 | 0.90 | 1.85 | 0.095 | 0.000 | 0/0/1 |
| rank1_plus | negative_control | 3202 | 0.178 | 0.538 | 0.19 | 0.01 | 0.49 | 0.52 | 1.01 | 0.047 | 0.040 | 0/1/1 |
| rank1_plus | negative_control | 6694 | 0.106 | 0.136 | 0.01 | 0.42 | 0.47 | 0.82 | 1.30 | 0.117 | 0.094 | 0/0/1 |
| rank1_plus | negative_control | 6891 | 0.408 | 0.079 | 0.10 | 0.64 | 0.59 | 0.82 | 1.41 | 0.049 | 0.102 | 0/0/1 |
| rank1_plus | negative_control | 4648 | 0.783 | 0.622 | 0.93 | 0.33 | 0.72 | 0.80 | 1.52 | 0.115 | 0.158 | 0/0/1 |
| rank1_plus | negative_control | 5924 | 0.523 | 0.872 | 0.34 | 0.05 | 0.85 | 0.79 | 1.65 | 0.164 | 0.188 | 0/0/1 |
| rank1_plus | negative_control | 5713 | 0.602 | 0.124 | 0.94 | 0.87 | 0.94 | 0.95 | 1.89 | 0.048 | 0.117 | 0/0/1 |

Selected with eft_control at cut 0.5: 5 of 40; at 0.95 all 40 by construction of the population.

## What E3 will compute per background (x=k/M, 1e-3..1e3, g_phi=M sqrt p, g_psi=-M kappa/sqrt p)
1 reality and simplicity of the six roots (QZ + mpmath); 2 positive definiteness of the reduced Hamiltonian h(k); 3 sign of canonical energy and Krein value; 4 Proca energy share and composition exchange (crossing x, |omega|/M there); 5 minimum branch gaps; 6 EFT window (|omega|/M<=eta, x<=eta) and protocol cuts, reported separately from numerical stability; 7 group velocities (reported as group, not front, velocities).
Negative controls: expected signature = complex roots and indefinite h; a diagnostic that reports positive h or real roots there is a failure and will be reported.
Pre-registered outcome table: exchange present/absent per background; Krein sign change yes/no; exchange inside EFT window yes/no. Failed backgrounds stay in all tables.
Open design question for you: for negative controls the local quartic is complex but the full 6-root system may be real at some x (Exp. 1 showed discriminant results only for the B3 row). Both outcomes will be recorded; no prior expectation is assumed.