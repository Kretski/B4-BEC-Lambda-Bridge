# Causality analysis of B4 itself (before any UV completion)

Working choice of the project lead (2026-11-12): B4 is a Lorentz-violating effective theory; the unit speed is the fundamental speed of light;
the preferred frame is a geometric observer field (static observer / ZAMO), not dynamical so far. This is an interpretation choice, not a claim about nature.

## 1. Characteristic cones of B4 on a background
The B4 EFT is second order in derivatives, so its characteristic speeds are the roots of the local polynomial
F(c) = d1 d2 + p (w - q c)^2 d2 + s (v - r c)^2 d1, d1 = c^2 - 2 a_phi c - (1 + 3 Lambda q^2), d2 = c^2 - 2 a_psi c - 1.
For this local theory phase, group and front velocities coincide (no dispersion at this order).
- Purely time-dependent background (q = r = 0), analytic: roots +-1 and +-sqrt(1 - p w^2 - s v^2). All speeds <= 1 for any Lambda.
- Generic backgrounds, numerical (experiments/uv_candidates/b4_local_cone_scan.py, 4000 backgrounds per cell, kappa^2 = p s, cuts X <= 0.09 and X <= 0.95):

| Lambda | cut | hyperbolic | max speed <= 1 | min of max speed | median of max speed |
|---|---|---|---|---|---|
| 0 | 0.09 | 4000 | 0 | 1.0000003 | 1.031 |
| 0 | 0.95 | 3853 | 0 | 1.0000001 | 1.215 |
| 0.01 | 0.09 | 4000 | 0 | 1.00012 | 1.034 |
| 0.01 | 0.95 | 3854 | 0 | 1.00019 | 1.220 |
| 0.1 | 0.09 | 4000 | 0 | 1.00006 | 1.056 |
| 0.1 | 0.95 | 3870 | 0 | 1.00031 | 1.259 |
| 1 | 0.09 | 4000 | 0 | 1.00009 | 1.194 |
| 1 | 0.95 | 3946 | 0 | 1.0012 | 1.753 |

Result: on every sampled generic background the widest B4 characteristic cone exceeds the unit cone, already at Lambda = 0.
The excess is produced by the p, s, kappa operators as well as by Lambda; it is a property of the EFT, not of a UV completion.

## 2. Phase, group and front velocity of the Lambda-model dispersion
omega = c k sqrt(1 + Lambda k^2):  v_phase = c sqrt(1 + Lambda k^2),  v_group = c (1 + 2 Lambda k^2)/sqrt(1 + Lambda k^2),  v_front = lim omega/k.
- Lambda > 0: v_phase > c and v_group > c for all k > 0; v_front is unbounded (omega ~ k^2). No finite cone.
- Lambda < 0: v_phase < c, but omega becomes imaginary for k > 1/sqrt(|Lambda|): the description must be cut off there.
- This Lambda is the coefficient of a linear k^4 term, i.e. beta = c Lambda/2 in omega = c k + beta k^3. It is not the coefficient of the B4 operator phi_x^4.

## 3. Common causal cone
The local B4 system is hyperbolic on most sampled backgrounds and has its own widest cone (max |root|), larger than the unit cone.
With the k^4 dispersion there is no finite cone at all; evolution is still well posed with respect to the preferred time.

## 4. Region with Lambda > 0 and all modes inside the unit cone
Only q = r = 0 in the preferred frame (analytic). None of 32000 sampled generic backgrounds.

## 5. What could restore fundamental causality
- Not a UV mediator of the kinds tested (Proca, scalar, scalar + dynamical eta): the excess is already in the EFT.
- (a) Restrict the theory to backgrounds with q = r = 0 in the observer frame, and show that this restriction is preserved by the dynamics.
- (b) Change the causality criterion: in Lorentz-violating theories with a preferred time foliation, causality is defined by that foliation and
  superluminal dispersion is not by itself acausal. This is how the literature on Horava-type / Einstein-aether theories treats omega^2 = k^2 + k^4/M^2
  (recalled from the literature, not verified here). It replaces "front <= 1" by well-posed evolution in the preferred time plus T6 for the frame field.
- (c) Keep the Minkowski cone and take Lambda <= 0 with an explicit cutoff.
These are choices for the project lead; (b) changes T3 and must be registered in the protocol before use.

## Related observable (T8), to be verified against the literature
A k^4 correction gives a gravitational-wave phase shift proportional to f^3. This appears to be the alpha = 4 case of the standard parametrised
dispersion tests (omega^2 = k^2 + A k^alpha, phase shift ~ f^(alpha-1)) used in LIGO-Virgo-KAGRA analyses, where bounds are quoted for both signs of A.
Recalled, not checked.
