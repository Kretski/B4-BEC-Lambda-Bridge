# T6, first analysis: the observer field in Kerr (ZAMO vs static observers)

Scope: geometry only (Kerr, Boyer-Lindquist coordinates, signature -+++, G = c = 1). No dynamics for the observer field is proposed or tested here.
Symbolic checks: experiments/uv_candidates/t6_zamo_checks.py (all identities below marked [S] are verified there).

## ZAMO congruence
u_mu = -alpha d_mu t, i.e. u^mu = (1/alpha)(1, 0, 0, omega), with alpha^2 = Delta Sigma / A and omega = 2 M a r / A [S].
Sigma = r^2 + a^2 cos^2 theta, Delta = r^2 - 2 M r + a^2, A = (r^2 + a^2)^2 - a^2 Delta sin^2 theta.
- Unit timelike where alpha^2 > 0 [S], i.e. Delta > 0: everywhere outside the outer horizon r_+ = M + sqrt(M^2 - a^2), including the ergoregion.
- Zero angular momentum: u_phi = 0 [S].
- Hypersurface orthogonal by construction (u is proportional to dt, so u ^ du = 0): zero vorticity. The preferred time is the Boyer-Lindquist t
  and the leaves are t = const.
- Expansion zero (u^mu has only t and phi components and nothing depends on t, phi). Not geodesic: acceleration a_mu = d_mu ln alpha. Shear nonzero when a != 0.
- At r -> r_+: alpha -> 0. The congruence becomes null and the foliation t = const does not extend through the horizon.

## Static observers
u proportional to the Killing vector d/dt. Timelike only where g_tt < 0, i.e. r > r_E(theta) = M + sqrt(M^2 - a^2 cos^2 theta) [S]: not defined inside the ergoregion.
Hypersurface orthogonal iff g_tphi / g_tt is constant; d/dr (g_tphi/g_tt) = -2 M a (r^2 - a^2 cos^2 theta) sin^2 theta / (Sigma - 2 M r)^2 [S], nonzero for a != 0.
So for a != 0 the static congruence has vorticity and does not define a time foliation orthogonal to it. For a = 0 it coincides with the ZAMO congruence.

## Consequences for the causality convention
1. Preferred-foliation causality needs a hypersurface-orthogonal observer field. In Kerr the ZAMO field satisfies this; the static field does not (a != 0).
   Under the registered convention the observer field in Kerr should therefore be the ZAMO field, with preferred time t.
2. Domain: r > r_+ only. Inside the horizon this construction gives no observer field and no preferred time. This is a limit of the construction, to be
   stated wherever the model is applied near a horizon.
3. Local frequency of a mode with Killing energy E and axial angular momentum L: omega_loc = (E - omega L)/alpha [S]. It grows without bound as alpha -> 0.
   With omega = c k_loc sqrt(1 + Lambda k_loc^2), the correction Lambda k_loc^2 is of order one where alpha ~ sqrt(Lambda) (E - omega L)/c (order of magnitude):
   close enough to the horizon the low-k expansion of the Lambda-model is not valid, for any Lambda > 0.

## What this does not settle (T6 remains open)
- Physical origin: ZAMOs are defined by the stationarity and axisymmetry of Kerr. A generic dynamical spacetime (e.g. a binary merger) has no such Killing
  structure, so "ZAMO" is not a general rule for the observer field.
- A common form: in Kerr u_mu = -alpha d_mu t; for comoving observers in a cosmological (FRW) background u_mu = -d_mu(cosmic time). Both are
  u_mu = -N d_mu T for a time function T. Taking "u_mu = -N d_mu T with T a preferred time function" as the definition would cover both cases and makes
  hypersurface orthogonality automatic; what determines T in general (a field equation, or a geometric rule) is the open T6 question.
  This is the structure used in khronometric / Horava-type theories (recalled from the literature, not verified here).
- For gravitational-wave propagation over cosmological distances the relevant frame is presumably the cosmological one, not the ZAMO frame of the source;
  how the two are matched is not addressed.
- Stability and well-posedness of the observer field if it is made dynamical (T6 item 5): not addressed.
