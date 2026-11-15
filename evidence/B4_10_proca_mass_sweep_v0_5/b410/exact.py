"""Exact (rational-arithmetic) statements about the universal polynomial Q(c; eps).

With eps = x^2/(1+x^2) in (0,1), x = k/M_A, the fixed-EFT spectrum solves
    Q(c; eps) = F(c) - eps c^2 d1(c) d2(c) = 0           (degree six in c).
The discriminant of Q with respect to c is a polynomial in eps.  If it has no
real root in (0,1], no two of the six roots can ever collide for 0 < x < inf, so the
number of real roots (and hence of complex-conjugate pairs) is constant on the whole
range and equals its value at small eps.

The floating-point B3 parameters are converted to exact rationals (their exact binary
values); the result is exact for those numbers.
"""
from __future__ import annotations
import sympy as sp
from .fixed_eft import B3Point, contact_coefficients


def universal_polynomial(b3: B3Point):
    c, e = sp.symbols("c e")
    R = sp.Rational
    aphi, apsi = contact_coefficients(b3)
    d1 = c**2 - R(1 + 3*b3.Lam*b3.q**2) - 2*c*R(aphi)
    d2 = c**2 - 1 - 2*c*R(apsi)
    F = sp.expand(d1*d2 + R(b3.p)*(R(b3.w) - R(b3.q)*c)**2*d2
                  + R(b3.s)*(R(b3.v) - R(b3.r)*c)**2*d1)
    Q = sp.expand(F - e*c**2*d1*d2)
    return Q, c, e


def discriminant_report(b3: B3Point, eps_small=sp.Rational(1, 10**4)):
    Q, c, e = universal_polynomial(b3)
    D = sp.Poly(sp.expand(sp.discriminant(sp.Poly(Q, c))), e)
    real_roots = [float(r) for r in D.real_roots()]
    roots_in_0_1 = [r for r in real_roots if 0.0 < r <= 1.0]
    q_small = sp.Poly(Q.subs(e, eps_small), c)
    return {
        "discriminant_degree_in_eps": D.degree(),
        "discriminant_real_roots": real_roots,
        "discriminant_real_roots_in_(0,1]": roots_in_0_1,
        "discriminant_sign_at_eps_small": int(sp.sign(D.eval(eps_small))),
        "discriminant_sign_at_eps_1": int(sp.sign(D.eval(1))),
        "n_real_roots_of_Q_at_eps_small": int(q_small.count_roots()),
        "eps_small": float(eps_small),
        "conclusion_all_six_real_and_simple_for_all_x":
            bool(len(roots_in_0_1) == 0 and q_small.count_roots() == 6
                 and sp.Poly(Q.subs(e, eps_small), c).degree() == 6),
    }
