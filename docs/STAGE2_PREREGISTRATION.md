# Stage 2 pre-registration, v2 (written BEFORE any Stage 2 result is used)
Disclosure: an earlier v1 run was started and stopped after 3 of 29 rows had printed (best residuals ~1.0 in region A; not used, deleted).
The criteria below follow the reviewer's plan and were fixed before any rerun.

## Classification of every BEC parameter set (energy diagnostic added, BdG core unchanged)
1. energetically stable: four real simple roots AND min eigenvalue of the lab-frame BdG energy form >= -1e-9 over the k grid  -> region A (primary)
2. dynamically stable, energetically unstable: four real simple roots, min eigenvalue < -1e-6                              -> region B (mandatory control)
3. dynamically unstable: complex roots (or immiscible g12^2 > g11 g22)                                                  -> excluded
4. boundary / ill-conditioned: repeated roots (gap < 1e-6), |H'| < 1e-2, or |lambda_min| between the two thresholds         -> excluded, counted
Energy form is frame dependent (V_i are lab velocities); stated with every result.

## B4 side
The 29 backgrounds with four real local roots (anchor + manifest). beta_i = A_i / M^2 (M^2 beta_i = A_i). Same conditioning flags as in lambda_eff.

## Unit-free comparison (beta is Galilean invariant; speed scale, length scale, shift are dictionary freedoms)
Per branch j (branches matched one-to-one):
  c-structure: x_j = (c_j - c_min)/(c_max - c_min);   magnitude: bhat_j = beta_j / max|beta|.
  eps_j = |bhat_j^B4 - bhat_j^BEC| / max(|bhat_j^B4|, |bhat_j^BEC|, eps0),   eps0 = 0.05.
Matching of branches: primary = ordered by root value (the only identification available on both sides without a dictionary);
sensitivity = best of all 24 one-to-one assignments, reported separately and labelled as having more freedom.
Four criteria, all required for "match":
  (1) sign: sign(bhat_j) equal on all branches with |bhat| >= 0.05 on both sides;
  (2) magnitude: max_j eps_j <= 0.25;
  (3) root structure: |x_j^B4 - x_j^BEC| <= 0.05 for the two inner roots;
  (4) robustness: (1)-(3) hold for >= 80% of 50 perturbed B4 targets (each invariant perturbed by relative 2%, absolute 0.005).
BEC family (4 params): rho = a2/a1, f = g12/sqrt(g11 g22) in [0, 0.99], dV = (V1 - V2)/sqrt(a1), kappa_b = b2/b1.
Search: 80 random starts per background, bounded least squares, seed 20261111; region membership enforced (hard rejection).
Also reported (no tuning): which sign tuples exist at all in A and B (random sampling, 2e5 points each).
Not tested: LHY (separate stage), nonlinear operators, common-coupling family across backgrounds. A match is necessary, not sufficient.
