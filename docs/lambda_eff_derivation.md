# Definition of A_i from the Proca action (B4.10 v0.5)

Where it was already written: README.md (step 3), `b410/fixed_eft.py::leading_shift_coefficients`, `b410/exact.py` docstring (Q(c;eps)).
None of them derived it from the action; `derive_Ai_from_action.py` (SymPy, independent of b410) now does.

## Action and background
L = 1/2 f_t^2 - 1/2 f_x^2 - (Lam/4) f_x^4 + 1/2 g_t^2 - 1/2 g_x^2 + 1/2 (A1_t - A0_x)^2 + 1/2 M^2 (A0^2 - A1^2) + A1 (g_phi f_t f_x + g_psi g_t g_x)
Background f = w t + q x, g = v t + r x, A0 = 0, A1 = Abar = (g_phi w q + g_psi v r)/M^2.

## Quadratic action (second variation; linear terms cancel on the background)
L2 = 1/2 f_t^2 - 1/2 (1+3 Lam q^2) f_x^2 + 1/2 g_t^2 - 1/2 g_x^2 + 1/2 (a1_t - a0_x)^2 + 1/2 M^2 (a0^2 - a1^2)
     + a1 [g_phi (q f_t + w f_x) + g_psi (r g_t + v g_x)] + Abar [g_phi f_t f_x + g_psi g_t g_x]
Coefficients entering the spectrum: 1+3 Lam q^2 (stiffness), contact terms
 a_phi = Abar g_phi = p w q - kappa v r,  a_psi = Abar g_psi = -kappa w q + s v r,   with p = g_phi^2/M^2, s = g_psi^2/M^2, kappa = -g_phi g_psi/M^2.

## Exact characteristic polynomial (verified symbolically)
Plane waves exp(i(kx - omega t)), c = omega/k, x = k/M:
 det K(omega, k) = - M^8 x^4 (1+x^2) * Q(c; x)  (ratio to Q equals 1 to 1e-16 numerically)
 Q(c; x) = F(c) - [x^2/(1+x^2)] c^2 d1(c) d2(c)
 d1 = c^2 - 2 a_phi c - (1 + 3 Lam q^2),  d2 = c^2 - 2 a_psi c - 1
 F  = d1 d2 + p (w - q c)^2 d2 + s (v - r c)^2 d1

## Definition of A_i
c_i(x) = c_i^0 + A_i x^2 + O(x^4), where F(c_i^0) = 0 (the four local B3 roots).
Implicit differentiation of Q(c; eps) = 0 with eps = x^2/(1+x^2) = x^2 + O(x^4):
 A_i = (dc/d eps)|_{eps=0} = c_i^2 d1(c_i) d2(c_i) / F'(c_i).
The two additional (heavy, Proca) roots are not given by this formula; they go like +-M/k at small x.

## Check, B3 anchor (M-independent)
roots c^0 = -4.9164, -0.1746, 2.7969, 7.3151; A_i = 17.8222, 0.001053, -7.7625, -11.4863 (SymPy and code agree to 1e-6; numeric spectrum at x = 0.01 matches c^0 + A x^2).
Units: A_i is in speed units (c = omega/k, c_B4 = 1); x = k/M. If a length-squared coefficient is wanted, Lam_eff,i = 2 A_i / (c_i^0 M^2) (definition by analogy to omega = c k sqrt(1+Lam k^2); it is root dependent and negative for the anchor).
Not covered here: which of these matches the nonlinear coefficient Lam (it enters only through 3 Lam q^2 in d1).
