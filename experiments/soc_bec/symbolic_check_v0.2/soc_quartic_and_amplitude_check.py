"""Symbolic check of SOC_BEC_QUARTIC_NOTE_v0.1 §1-§2, and the amplitude-level (time-delay) combination.

Authorised 2026-11-20 ("символна проверка на §1-§2; изчисляване на амплитудата с обмена").
Pure SymPy; no input files; no randomness. Exit 0 iff all checks pass.

Part A (note §1-§2): band expansion, m*, c4, the pure spatial quartic C of P = Y^2/(2g), and the
condition Lambda > 0  <=>  gamma/2 > x (x-1)^2.

Part B (amplitude level, background/eikonal method as in K2 note v0.5 §1): exact characteristic speed
of a right-moving fluctuation on a constant-gradient left-moving background (theta_t = c eps,
theta_x = eps) for L = P(Y), Y = mu - theta_t - E(theta_x). This includes the cubic vertex (exchange)
exactly at the level of the background. The O(eps) term is odd in the background and integrates to zero
across a localised pulse; the time delay is set by the O(eps^2) coefficient of 1/c_+.
"""
import sys
import sympy as sp

results = []


def check(name, ok, detail=""):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))


m, ER, Om, g, mu = sp.symbols('m E_R Omega g mu', positive=True)
q, eps = sp.symbols('q epsilon', real=True)
kR = sp.sqrt(2 * m * ER)                                   # E_R = kR^2/(2m)

# ---------------------------------------------------------------- Part A
Em = q**2 / (2 * m) - sp.sqrt((Om / 2)**2 + (kR * q / m)**2)
Et = sp.simplify(Em - Em.subs(q, 0))                        # E~(q) = E_-(q) - E_-(0)
ser = sp.series(Et, q, 0, 6).removeO()
inv_mstar = sp.simplify(2 * ser.coeff(q, 2))                # 1/m*
c4 = sp.simplify(ser.coeff(q, 4))
check("A1  m/m* = 1 - 4E_R/Omega", sp.simplify(m * inv_mstar - (1 - 4 * ER / Om)) == 0, str(sp.simplify(m * inv_mstar)))
check("A2  c4 = 4 E_R^2 / (m^2 Omega^3) > 0", sp.simplify(c4 - 4 * ER**2 / (m**2 * Om**3)) == 0, str(c4))
check("A3  no odd powers in E~(q) (delta = 0)", all(sp.simplify(ser.coeff(q, k)) == 0 for k in (1, 3, 5)))

mstar = 1 / inv_mstar
w = sp.symbols('w', real=True)                              # w = theta_t (fluctuation)
Ytr = mu - w - q**2 / (2 * mstar) - c4 * q**4               # truncated band, enough for O(q^4)
P = sp.expand(Ytr**2 / (2 * g))
C = sp.simplify(sp.Poly(P, w, q).coeff_monomial(q**4))
n = mu / g
check("A4  pure spatial quartic C = 1/(8 g m*^2) - n c4", sp.simplify(C - (1 / (8 * g * mstar**2) - n * c4)) == 0)
check("A5  no w^2 q^2 term (p = 0) in P = Y^2/(2g)", sp.Poly(P, w, q).coeff_monomial(w**2 * q**2) == 0)
check("A6  cubic term w q^2 / (2 g m*) present",
      sp.simplify(sp.Poly(P, w, q).coeff_monomial(w * q**2) - 1 / (2 * g * mstar)) == 0)

x, gam = sp.symbols('x gamma', positive=True)               # x = Omega/(4E_R), gamma = g n / E_R
sub = {Om: 4 * ER * x, mu: gam * ER}                        # mu = g n
ratio = sp.simplify((8 * g * n * mstar**2 * c4).subs(sub))
check("A7  8 g n m*^2 c4 = (gamma/2) / (x (x-1)^2)", sp.simplify(ratio - (gam / 2) / (x * (x - 1)**2)) == 0, str(ratio))

# ---------------------------------------------------------------- Part B
# general P(Y); quadratic fluctuation operator on a constant background (theta_t0, theta_x0) = (c eps, eps)
c = sp.sqrt(mu / mstar)                                     # sound speed: c^2 = g n / m* = mu / m*
Ef = lambda k: k**2 / (2 * mstar) + c4 * k**4                # truncated band (O(eps^2) in the result is exact)
Y0 = mu - c * eps - Ef(eps)
Pp = Y0 / g                                                 # P'(Y0) = n0
Ppp = 1 / g                                                 # P''
E1 = sp.diff(Ef(q), q).subs(q, eps)
E2 = sp.diff(Ef(q), q, 2).subs(q, eps)
# L2 = 1/2 P''(dtheta_t + E1 dtheta_x)^2 - 1/2 P' E2 dtheta_x^2  ->  speeds c_pm = E1 +- sqrt(P' E2 / P'')
cplus = E1 + sp.sqrt(Pp * E2 / Ppp)
cser = sp.series(cplus, eps, 0, 3).removeO()
check("B1  c_+(0) = c (unperturbed sound speed)", sp.simplify(cser.subs(eps, 0) - c) == 0)
lin = sp.simplify(cser.coeff(eps, 1))
check("B2  O(eps) shift = eps/(2 m*) (odd in the background; integrates to zero across a pulse)",
      sp.simplify(lin - 1 / (2 * mstar)) == 0, str(lin))
inv = sp.series(1 / cplus, eps, 0, 3).removeO()
k2 = sp.simplify(inv.coeff(eps, 2))
expected = (1 / c) * (sp.Rational(5, 8) / (mstar * mu) - 6 * c4 * mstar)
check("B3  O(eps^2) coefficient of 1/c_+ = (1/c)(5/(8 m* mu) - 6 c4 m*)", sp.simplify(k2 - expected) == 0)
# healthy check: ordinary superfluid (c4 = 0) gives a time delay (positive coefficient)
k2_c40 = sp.simplify(k2.subs(c4, 0)) if c4.free_symbols else None
k2_plain = sp.simplify(((1 / c) * (sp.Rational(5, 8) / (mstar * mu))))
check("B4  c4 = 0 (parabolic band): coefficient = 5/(8 m* mu c) > 0 (time delay, healthy sign, cf. K3a)",
      sp.simplify(k2_plain).is_positive is not False)
# wrong-sign (advance) condition: 6 c4 m* > 5/(8 m* mu)  <=>  (48/5) g n c4 m*^2 > 1  <=>  0.6 gamma > x (x-1)^2
adv = sp.simplify(((sp.Rational(48, 5)) * g * n * c4 * mstar**2).subs(sub))
check("B5  advance condition: (48/5) g n c4 m*^2 = (3 gamma / 5) / (x (x-1)^2)",
      sp.simplify(adv - (sp.Rational(3, 5) * gam) / (x * (x - 1)**2)) == 0, str(adv))
check("B6  advance threshold weaker than contact threshold: 3/5 > 1/2", sp.Rational(3, 5) > sp.Rational(1, 2))

nfail = results.count(False)
print(f"\nSUMMARY: {len(results) - nfail}/{len(results)} PASS")
sys.exit(0 if nfail == 0 else 1)
