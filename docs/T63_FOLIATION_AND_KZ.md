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

## Comparison with the two historical formulas of the Lambda-model (2026-11-12)
Equivalent dimensionless forms, with [Lambda] = m^2 and f in Hz:
    delta Psi = -4 pi^3 Lambda K_new f^3 / c^3 = -4 pi^3 Lambda K_old f^3 / c^2,   K_old = Int_0^z (1+z')^2 / H dz' [s],  K_new = c K_old [m].
Group shift from the phase: delta t(f) = (1/2 pi) d(delta Psi)/df = -6 pi^2 Lambda K_old f^2 / c^2 = -(3 Lambda omega^2 / 2 c^2) K_old.
The phase coefficient 4 pi^3 and the time coefficient 3/2 are mutually consistent. (delta Psi = 2 pi f delta t is not valid for a frequency-dependent group delay.)

1. Matched-filter code of 2026-09-01: delta Psi = -4 pi^3 Lambda K_old f^3 / c^3 with K_old in seconds.
   Same f^3 dependence, same integrand (1+z')^2/H, same sign. With Lambda in m^2 the expression has dimension s/m: either Lambda is in other units there,
   or one factor of c is missing (a factor 3e8 in any numerical bound). OPEN until the code is inspected.
2. Formula sheet of 2026-08-22: Delta t = (3 Lambda omega_obs^2 / 2 c^2) Int_0^z (1+z')/H dz'.
   Coefficient: agrees. Sign: agrees if Delta t is defined as t_low - t_high (must be stated). Integrand: (1+z') is not what this dispersion law gives;
   the derivation and the numerical integration give (1+z')^2:

   | z | numeric delta t / (1.5 Lambda k0^2) | Int (1+z)^2/E | Int (1+z)/E |
   |---|---|---|---|
   | 0.1 | 0.107639 | 0.107640 | 0.102476 |
   | 0.5 | 0.683156 | 0.683158 | 0.544266 |
   | 1.0 | 1.678379 | 1.678386 | 1.113169 |
   | 2.0 | 4.322516 | 4.322545 | 2.175142 |
   (absolute values; dimensionless, c = H0 = 1, Omega_m = 0.31)
   The k_loc^2 dependence is not an extra assumption: v_g = d omega/dk = c (1 + 2 Lambda k^2)/sqrt(1 + Lambda k^2) follows from omega = c k sqrt(1 + Lambda k^2).
   This is the same cosmological dependence as the quadratic case (n = 2) of Jacob and Piran (2008): coefficient (1+n)/2, integrand (1+z)^n.
   Only the mathematical structure of the propagation is shared; this is not a claim that the Lambda-model is an LIV model, and the physical
   interpretation of Lambda is a separate question. (Reference checked by the second reviewer, not by this analysis.)

Sign convention: the negative cubic term means earlier arrival of high frequencies for Lambda > 0 provided h~(f) = A exp(i Psi) with Psi_GR containing +2 pi f t_c.
It must not be carried over to another Fourier convention without conversion.

## Inspection of the 2026-09-01 matched-filter code (waveform.py, fragments supplied by the project lead, 2026-11-12)
Code facts:
- cosmological_K_factor: K = Int_0^z (1+z')^2/H dz' with H in 1/s, so K is in seconds (H0 = 67.4, Om = 0.315, OL = 0.685).
- lambda_phase_correction: -(4 pi^3 Lambda K_z / C_SI^3) f^3.
- taylorf2_leading_order: Psi_GR = +2 pi f tc - phi_c - pi/4 + ..., and the waveform is A exp(+i Psi_total).
- Docstring derivation: Delta T(omega) = -(Lambda / 2 c^3) omega^2 K(z), then Delta Psi = 2 pi f Delta T.

T8.2 (Fourier convention): confirmed in the code itself. With +2 pi f tc and exp(+i Psi), the negative cubic term means earlier arrival of high
frequencies for Lambda > 0, as in the derivation above.

T8.1 (dimensions and normalisation), assuming Lambda is defined by omega = c k sqrt(1 + Lambda k^2) with [Lambda] = m^2:
- Correct phase with K in seconds: -(4 pi^3 Lambda K / c^2) f^3. The code divides by c^3. The code expression has dimension s/m.
- Therefore the parameter called Lambda in the code equals c * Lambda numerically in SI:  Lambda_true [m^2] = Lambda_code / c  (c = 2.998e8).
  The code underestimates the phase for a given Lambda by a factor c; any bound or sensitivity quoted for Lambda in m^2 from this code is too loose by that factor.
- Because the phase is linear in Lambda, existing results do not need to be re-run: the Lambda axis (values, grids, priors, bounds) is rescaled by exactly 1/c.
  Injection-recovery tests that use the same template for injection and recovery remain valid as pipeline tests.
- The docstring derivation contains two errors that cancel in the numerical coefficient: the time shift lacks the factor 3 (correct: -(3 Lambda/2c^2) omega^2 K),
  and Delta Psi = 2 pi f Delta T is used instead of Delta Psi = Int 2 pi Delta T df (which gives a factor 1/3). The coefficient 4 pi^3 is right; the power of c is not.
- Not checked: other files using the same formula (lalsim_waveform.py, likelihood.py, the gwosc_* scripts, Paper 2), and where numerical bounds were published.

## Status
- Cosmological propagation with T = cosmic time: derived and checked numerically (this note).
- T6.3 (global compatibility source-detector): NOT established. No single function T joining the Kerr/ZAMO region to the cosmological foliation is constructed.
- T8.1: resolved for waveform.py as stated above (factor c, under the stated definition of Lambda). Other files and published numbers: open.
- T8.2: resolved for waveform.py.
- T8 as a whole: formula and convention consistent after the 1/c rescaling; not confirmed as an observable until the remaining files and quoted bounds are checked.
