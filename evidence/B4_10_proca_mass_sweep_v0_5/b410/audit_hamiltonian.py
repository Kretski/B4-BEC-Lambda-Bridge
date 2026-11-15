"""Independent audit of Experiment 2 (does NOT use H_B49, H', null-vector Noether code).

Source of truth: the T1.7 action
    L = 1/2 phi_t^2 - 1/2 phi_x^2 - (Lam/4) phi_x^4 + 1/2 psi_t^2 - 1/2 psi_x^2
        + 1/2 (A1_t - A0_x)^2 + 1/2 M^2 (A0^2 - A1^2) + A1 (g_phi phi_t phi_x + g_psi psi_t psi_x)
(J_phi^1 = phi_t phi_x so that (1/2M^2) J.J = (p/2)(phi_t phi_x)^2 and kappa = -g_phi g_psi/M^2,
as in the T1.7.1 text.)  Background: phi0 = (w t + q x), psi0 = (v t + r x), A1 = Abar = J0/M^2.
Steps: (1) sympy second variation of L; (2) A0 constraint; (3) Legendre transform with A0
eliminated -> 6x6 Hamiltonian matrix h(k); (4) Hamilton-equation frequencies and positivity.
"""
from __future__ import annotations
import numpy as np
import sympy as sp


def second_variation_coefficients(M, gphi, gpsi, w, q, v, r, Lam):
    """sympy: expand L to O(eps^2) with Phi = bg + eps*delta; return dict of bilinear coefficients."""
    eps = sp.Symbol('eps')
    names = "dphi_t dphi_x dpsi_t dpsi_x da0 da0_x da1 da1_t".split()
    S = {n: sp.Symbol(n) for n in names}
    Mq, gp, gs, W, Q, V, R, Lm = (sp.nsimplify(float(x), rational=True) if False else sp.Float(float(x), 30)
                                  for x in (M, gphi, gpsi, w, q, v, r, Lam))
    Abar = (gp*W*Q + gs*V*R)/Mq**2
    pt, px = W + eps*S['dphi_t'], Q + eps*S['dphi_x']
    st, sx = V + eps*S['dpsi_t'], R + eps*S['dpsi_x']
    a0, a0x = eps*S['da0'], eps*S['da0_x']
    a1, a1t = Abar + eps*S['da1'], eps*S['da1_t']
    L = (sp.Rational(1, 2)*pt**2 - sp.Rational(1, 2)*px**2 - Lm/4*px**4
         + sp.Rational(1, 2)*st**2 - sp.Rational(1, 2)*sx**2
         + sp.Rational(1, 2)*(a1t - a0x)**2 + Mq**2/2*(a0**2 - a1**2) + a1*(gp*pt*px + gs*st*sx))
    ser = sp.expand(L)
    L1 = ser.coeff(eps, 1); L2 = ser.coeff(eps, 2)
    return S, sp.expand(L1), sp.expand(L2)


def terms_from_L2(S, L2):
    """Convert sympy L2 into (coef, field_i, d_i, field_j, d_j) tuples comparable with energy.terms."""
    poly = sp.Poly(L2, *S.values())
    out = {}
    for mon, c in poly.terms():
        idx = [n for n, e in zip(S.keys(), mon) for _ in range(e)]
        out[tuple(sorted(idx))] = float(c)
    return out


def hamiltonian_matrix(k, M, gphi, gpsi, w, q, v, r, Lam):
    """6x6 Hermitian h(k): H2 = (1/2) u^dag h u, u=(phi,psi,a1,P_phi,P_psi,P_a1), A0 eliminated.

    H_phi = 1/2 (P_phi - Yphi phi_x - g q a)^2 + 1/2(1+3 Lam q^2) phi_x^2 - g w a phi_x   (similarly psi)
    H_A   = 1/2 P_a^2 (1 + k^2/M^2) + 1/2 M^2 a^2   (A0 = -P_a'/M^2 from its constraint)
    """
    Abar = (gphi*w*q + gpsi*v*r)/M**2
    ik = 1j*k
    # linear maps: residual r_i = P_i - ik*Y_i*phi_i - g_i*b_i*a ;  u -> components
    def e(i):
        v_ = np.zeros(6, complex); v_[i] = 1.0; return v_
    rows = []
    h = np.zeros((6, 6), complex)
    for (g, a_, b_, Y, grad2, ph, pm) in ((gphi, w, q, Abar*gphi, 1 + 3*Lam*q*q, 0, 3),
                                          (gpsi, v, r, Abar*gpsi, 1.0, 1, 4)):
        rvec = e(pm) - ik*Y*e(ph) - g*b_*e(2)          # velocity = r_i
        h += np.outer(rvec.conj(), rvec)
        gradv = ik*e(ph)
        h += grad2*np.outer(gradv.conj(), gradv)
        # - g a_ a phi_x : hermitian cross term  -g a_ (a^* phi_x + c.c.)
        cross = -g*a_*(np.outer(e(2).conj(), gradv) + np.outer(gradv.conj(), e(2)))
        h += cross
    h += (1 + k*k/M**2)*np.outer(e(5).conj(), e(5)) + M**2*np.outer(e(2).conj(), e(2))
    return h


def hamilton_frequencies(h):
    """omega from  du/dt = J h u,  u ~ exp(-i omega t):  -i omega u = J h u."""
    J = np.block([[np.zeros((3, 3)), np.eye(3)], [-np.eye(3), np.zeros((3, 3))]])
    # u = (q,p): need J acting as (q,p)->(dH/dp, -dH/dq); h is 6x6 in order (q,p)
    lam = np.linalg.eigvals(J @ h)
    return np.sort((1j*lam).real), np.abs((1j*lam).imag).max()


def mode_energy_from_h(h, omega, k, M, gphi, gpsi, w, q, v, r, Lam, z_true):
    """cycle-averaged H2 = (1/4) u^dag h u for the mode with true complex amplitudes z_true=(phi,psi,a0,a1)."""
    Abar = (gphi*w*q + gpsi*v*r)/M**2
    phi, psi, a0, a1 = z_true
    Pphi = (-1j*omega + 1j*k*Abar*gphi)*phi + gphi*q*a1
    Ppsi = (-1j*omega + 1j*k*Abar*gpsi)*psi + gpsi*r*a1
    Pa = -1j*omega*a1 - 1j*k*a0
    u = np.array([phi, psi, a1, Pphi, Ppsi, Pa])
    return float(np.real(0.25*np.vdot(u, h @ u))), u


def energy_split_hamiltonian(u, k, M, gphi, gpsi, w, q, v, r, Lam):
    """Hamiltonian-representation split of the cycle-averaged H2 = (1/4) u^dag h u.

    scalar_i = (1/4)[|r_i|^2 + (1+3 Lam q^2)|k phi_i|^2],  A = (1/4)[(1+k^2/M^2)|P_a|^2 + M^2|a|^2],
    interaction = the -g w a phi_x cross terms.  Sums exactly to the total.
    NOTE: the split is representation dependent (differs from the Noether-term split).
    """
    Abar = (gphi*w*q + gpsi*v*r)/M**2
    ik = 1j*k
    out = {"phi": 0.0, "psi": 0.0, "int": 0.0}
    for name, (g, a_, b_, Y, grad2, ph, pm) in (("phi", (gphi, w, q, Abar*gphi, 1 + 3*Lam*q*q, 0, 3)),
                                                ("psi", (gpsi, v, r, Abar*gpsi, 1.0, 1, 4))):
        rv = u[pm] - ik*Y*u[ph] - g*b_*u[2]
        out[name] = 0.25*(abs(rv)**2 + grad2*abs(ik*u[ph])**2)
        out["int"] += 0.25*(-g*a_)*2*np.real(np.conj(u[2])*(ik*u[ph]))
    out["A"] = 0.25*((1 + k*k/M**2)*abs(u[5])**2 + M**2*abs(u[2])**2)
    return out
