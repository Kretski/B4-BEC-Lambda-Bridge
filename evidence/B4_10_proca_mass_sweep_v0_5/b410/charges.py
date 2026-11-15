"""Lab energy, shift charges and the shift-symmetric energy, all from the NONLINEAR action.

L = 1/2 phi_t^2 - 1/2 phi_x^2 - Lam/4 phi_x^4 + 1/2 psi_t^2 - 1/2 psi_x^2 + 1/2 (A1_t-A0_x)^2
    + 1/2 M^2 (A0^2 - A1^2) + A1 (g_phi phi_t phi_x + g_psi psi_t psi_x).
Shift charges (phi -> phi+c, psi -> psi+c): Q_phi = dL/dphi_t, Q_psi = dL/dpsi_t.  Background phi0 = w t + q x,
psi0 = v t + r x, so mu_phi = w, mu_psi = v.  Expanding with Phi = bg + eps*delta (first-order waves only):

  T00 = sum_i (dL/dPhi_t,i) Phi_t,i - L           -> eps^2 coefficient  E_lab   (lab-frame second-order energy)
  Q_phi, Q_psi                                       -> eps^2 coefficients  Q2_phi, Q2_psi
  E_shift = E_lab - w Q2_phi - v Q2_psi              (= second variation of H - mu Q)

Claim tested numerically: E_shift equals the Legendre/Noether energy of the quadratic Lagrangian L2
(= 1/4 omega z^dag H' z = 1/4 u^dag h u).  Linear terms in second-order fields are total derivatives
with constant coefficients and are dropped (zero cycle average).
"""
from __future__ import annotations
import numpy as np
import sympy as sp

_cache = {}


def _build():
    if "f" in _cache:
        return _cache["f"]
    eps = sp.Symbol('eps')
    M, g1, g2, w, q, v, r, Lam = sp.symbols('M g1 g2 w q v r Lam')
    names = "dphi_t dphi_x dpsi_t dpsi_x da0 da0_x da1 da1_t".split()
    S = {n: sp.Symbol(n) for n in names}
    Abar = (g1*w*q + g2*v*r)/M**2
    pt, px = w + eps*S['dphi_t'], q + eps*S['dphi_x']
    st, sx = v + eps*S['dpsi_t'], r + eps*S['dpsi_x']
    a0, a0x = eps*S['da0'], eps*S['da0_x']
    a1, a1t = Abar + eps*S['da1'], eps*S['da1_t']
    P_t, P_x, S_t, S_x, A1t = sp.symbols('P_t P_x S_t S_x A1t')
    A1s = sp.Symbol('A1s'); A0s, A0xs = sp.symbols('A0s A0xs')
    L = (sp.Rational(1, 2)*P_t**2 - sp.Rational(1, 2)*P_x**2 - Lam/4*P_x**4 + sp.Rational(1, 2)*S_t**2 - sp.Rational(1, 2)*S_x**2
         + sp.Rational(1, 2)*(A1t - A0xs)**2 + M**2/2*(A0s**2 - A1s**2) + A1s*(g1*P_t*P_x + g2*S_t*S_x))
    pi_phi = sp.diff(L, P_t); pi_psi = sp.diff(L, S_t); pi_a = sp.diff(L, A1t)
    T00 = pi_phi*P_t + pi_psi*S_t + pi_a*A1t - L
    sub = {P_t: pt, P_x: px, S_t: st, S_x: sx, A1t: a1t, A1s: a1, A0s: a0, A0xs: a0x}
    def e2(expr):
        return sp.expand(expr.subs(sub)).coeff(eps, 2)
    exprs = {"E_lab": e2(T00), "Q2_phi": e2(pi_phi), "Q2_psi": e2(pi_psi)}
    params = (M, g1, g2, w, q, v, r, Lam)
    out = {}
    for k, ex in exprs.items():
        poly = sp.Poly(ex, *S.values())
        terms = []
        for mon, c in poly.terms():
            idx = [n for n, e in zip(S.keys(), mon) for _ in range(e)]
            terms.append((sp.lambdify(params, c, 'numpy'), tuple(idx)))
        out[k] = terms
    _cache["f"] = out
    return out


def cycle_average(kind, params, z, omega, k):
    """<eps^2 coefficient> for plane wave Re[z e^{i(kx-omega t)}], z = (phi, psi, a0, a1) true complex amplitudes."""
    phi, psi, a0, a1 = z
    amp = {'dphi_t': -1j*omega*phi, 'dphi_x': 1j*k*phi, 'dpsi_t': -1j*omega*psi, 'dpsi_x': 1j*k*psi,
           'da0': a0, 'da0_x': 1j*k*a0, 'da1': a1, 'da1_t': -1j*omega*a1}
    tot = 0.0
    for f, idx in _build()[kind]:
        c = float(f(*params))
        if len(idx) == 2:
            tot += c*0.5*np.real(amp[idx[0]]*np.conj(amp[idx[1]]))
        elif len(idx) == 0:
            tot += 0.0
        else:
            raise ValueError("unexpected monomial degree")
    return float(tot)


def energies(params, z, omega, k):
    """dict(E_lab, Q2_phi, Q2_psi, E_shift) for one mode; params=(M,g1,g2,w,q,v,r,Lam)."""
    M, g1, g2, w, q, v, r, Lam = params
    El = cycle_average("E_lab", params, z, omega, k)
    Qp = cycle_average("Q2_phi", params, z, omega, k)
    Qs = cycle_average("Q2_psi", params, z, omega, k)
    return dict(E_lab=El, Q2_phi=Qp, Q2_psi=Qs, E_shift=El - w*Qp - v*Qs)
