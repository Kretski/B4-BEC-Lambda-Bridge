# T6.3: preferred foliation from source to detector and re-derivation of K(z)

Working definition (project lead, 2026-11-12): u_mu = -N d_mu T, a kinematic definition of the observer field, not a model of its origin.
T6 is split into T6.1 (timelike gradient of T), T6.2 (geometric choice of T), T6.3 (global compatibility source-detector), T6.4 (dynamical origin).
This note addresses T6.3 for propagation over cosmological distances. Script: experiments/uv_candidates/t63_Kz_derivation_check.py.

## Assumptions (each one is an input, not a result)
A1. Along the propagation path the foliation is the cosmological one: T = cosmic time, N = 1, flat FLRW.
A2. Local dispersion measured by u: omega_loc = c k_loc sqrt(1 + Lambda k_loc^2), with omega_loc = -u^mu k_mu and k_loc the proper wavenumber.
A3. Lambda is constant in proper units along the path.
A4. First order in Lambda k_loc^2.
A5. Stationary-phase convention h~(f) = A exp(i Psi) with Psi = 2 pi f t_c - ..., so dPsi/df = 2 pi t(f). The opposite Fourier convention flips the overall sign.

## Derivation
Comoving wavenumber is conserved, k_loc = k0 (1+z), k0 = 2 pi f / c with f the observed frequency.
Group velocity v_g = c (1 + 2 Lambda k^2)/sqrt(1 + Lambda k^2) = c (1 + 3/2 Lambda k^2 + ...). Equal comoving distance for all frequencies gives the arrival-time shift
    delta t(f) = -6 pi^2 Lambda f^2 K(z) / c^3,     K(z) = c * Int_0^z (1+z')^2 / H(z') dz'   (a length).
Integrating dPsi/df = 2 pi delta t:
    delta Psi(f) = -4 pi^3 Lambda K(z) f^3 / c^3   (plus terms absorbed in t_c and phi_c).
This reproduces the form of the Lambda-model formula; K(z) is fixed by A1-A4. For Lambda > 0 high frequencies arrive earlier.

## Numerical check
Exact group velocity integrated through the expansion (flat LCDM, Omega_m = 0.31) against the first-order formula:
ratio 1.00000 for Lambda k0^2 = 1e-6 and 0.9998 / 0.9996 for 1e-4 (z = 0.1 / 1.0). The deviation is the expected second-order term.

## Values (H0 = 67.7 km/s/Mpc, Omega_m = 0.31, illustrative)
| z | K(z) [Mpc] | D_L [Mpc] | K/D_L |
|---|---|---|---|
| 0.01 | 44.6 | 44.6 | 1.000 |
| 0.1 | 476.7 | 475.6 | 1.002 |
| 0.2 | 1020.0 | 1011.9 | 1.008 |
| 0.5 | 3025.2 | 2919.1 | 1.036 |
| 1.0 | 7432.3 | 6792.1 | 1.094 |
| 2.0 | 19141.3 | 15931.8 | 1.201 |
Scale: |delta Psi| = 1 rad at 100 Hz corresponds to Lambda = 1.5e-8 m^2 (z = 0.1) and 9.5e-10 m^2 (z = 1), i.e. sqrt(Lambda) ~ 1e-4 m and 3e-5 m.

## Global compatibility (what is and is not shown)
- The source region (Kerr, T = t_BL, ZAMO) and the cosmological region are not joined into one function T here. No such global T is constructed.
- Order-of-magnitude argument only: the near-zone path length is of order GM/c^2 (about 1e2 km for stellar-mass binaries) against K ~ 1e21-1e23 km, so
  the near-zone contribution to delta Psi is suppressed by ~1e-19 or more unless the wave is traced to alpha -> 0, where the low-k expansion fails anyway
  (T6_OBSERVER_FIELD_ZAMO.md). A mismatch between the asymptotic rest frame of the source and the comoving frame (peculiar velocity ~1e-3 c) changes
  k_loc at relative order v/c. Neither statement is a derivation.
- Open: an explicit matching of the two foliations (T6.3 proper), inhomogeneities along the path, and whether A3 holds.

## Comparison with the standard parametrisation (not verified)
The structure matches the alpha = 4 case of parametrised dispersion tests (E^2 = p^2 c^2 + A p^alpha c^alpha), whose distance measure also involves
Int (1+z')^(alpha-2)/H dz'. Conventions for the prefactor and sign differ between papers; this must be checked against the source before any comparison with published bounds.
