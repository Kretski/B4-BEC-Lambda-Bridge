# B3 / B4.9 handoff package — inputs for B4.10

## Included
- `B3_max_speed_7_315096.json`: recovered original B3 maximum-speed row (rank1_minus, row index 8899), including the four B3 characteristic roots and scan metadata.
- `B4_9_proca_eigenvector_tracking_v0_1/`: B4.9 source, tests, README, and original result CSV/JSON outputs.
- In `proca_dispersion/eigenmodes.py`, the two corrections made in the current review are included:
  1. right singular vector reconstructed as `vh[-1].conj()`;
  2. `null_residual` computed directly as `||H @ v|| / max(1, sigma_max(H))`.

## B3 maximum-speed background

Metric/sign and Fourier conventions are documented in the B4.9 README and source. The exact scalar background and EFT coefficients are:

| Quantity | Value |
|---|---:|
| p | 0.05072130549385561 |
| s | 0.32639425996872784 |
| kappa | -0.12866679047568869 |
| Lambda | 0.9645378634411128 |
| w | -3.6323792917796456 |
| q | -4.237158273982973 |
| v | -0.47662815028331484 |
| r | -1.6064767477012714 |
| max |c| | 7.315096255616975 |

B3 roots (sorted here by numerical value):
`[-4.916360156290043, -0.17461942910556702, 2.7969453770704584, 7.315096255616975]`.

The row is marked `eft_domain=true`, `kinetic_positive=true`, and `hyperbolic_at_sample=true` in the recovered scan record. This is one sampled background, not a global theorem.

## What is held fixed in the proposed B4.10 mass sweep

**Chosen prescription: hold the low-energy EFT coefficients `p`, `s`, and `kappa`, as well as the scalar background `(w,q,v,r,Lambda)`, fixed.** For each positive mediator mass `M_A`, reconstruct the one-mediator couplings as

`g_phi(M_A) = M_A * sqrt(p)`

`g_psi(M_A) = -M_A * kappa / sqrt(p)`

under the convention `p=g_phi^2/M_A^2`, `s=g_psi^2/M_A^2`, `kappa=-g_phi*g_psi/M_A^2`.

This keeps the rank-one EFT Wilson coefficients fixed (up to floating-point rounding). It is **not** the B4.9 benchmark sweep, which held `g_phi`, `g_psi`, and the scalar background fixed. Those are different families of theories.

For the selected fixed-EFT prescription, `g_a` scales as `M_A`; `J0 = g_phi*w*q + g_psi*v*r` scales as `M_A`; hence `Abar=J0/M_A^2` scales as `1/M_A`. The products `Abar*g_a` in the sourced-background contact terms remain finite (constant for fixed scalar background and fixed EFT coefficients). The full theory should approach its local EFT description for fixed external `omega,k` as `M_A` becomes large, provided the derivative expansion remains controlled (`|omega|/M_A`, `k/M_A` small).

Do not describe this as the conventional fixed-coupling decoupling limit. If instead `g_phi,g_psi` were held fixed, `p,s,kappa` would scale as `1/M_A^2` and the EFT interactions would vanish as `M_A -> infinity`.

## Important caveat: domain check before B4.10

The recovered row has `Xw=0.669225999997443` and `Xv=0.07414841809678635`, so `Xw+Xv=0.7433744180942294`. This is above the earlier T1.7 cut `0.7`. Do not assume the T1.7 cut=.7 conclusions apply to this row. It is nevertheless marked `eft_control=true` in the recovered row. Recompute the precise definitions of all EFT control parameters from the original B3 code/protocol before making any domain claim. Also explicitly evaluate `p*w^2`, `s*v^2`, kinetic-matrix eigenvalues, and any other declared EFT expansion parameters; do not substitute these for the protocol's own definitions.

## Numerical validation priorities for B4.10
1. Verify `p,s,kappa` reconstructed from couplings at every `M_A`.
2. Test the fixed-EFT `M_A -> infinity` limit against the local B3 characteristic roots at fixed physical `k` (and separately document how the result changes if `k/M_A` is held fixed).
3. At each candidate complex root, compare polynomial roots with an independent generalized/linearized eigenvalue formulation and higher precision (`mpmath`).
4. Report `sigma_min(H)`, `sigma_2(H)`, root condition estimates, and normalized residuals; do not classify tiny imaginary parts as physical without precision/convergence checks.
5. Track modes by eigenvector overlap; use adaptive steps and invariant-subspace tracking near degeneracies.
6. Keep the original B4.9 illustrative benchmark distinct from this original B3 maximum-speed background.

## B4.9 scope and status

The B4.9 README documents a separate illustrative benchmark (`g_phi^2=0.05`, `g_psi=0`, `w=q=1`, `v=r=0`, `Lambda=0`) and its avoided crossing. Its fixed-coupling mediator-mass sweep is not the fixed-EFT prescription selected above. The original archive's 15-test result is recorded in `LOCAL_TEST_OUTPUT.txt`; after the two eigenvector/residual corrections above, rerun tests in the recipient environment. No B4.10 numerical result is claimed by this handoff package.
