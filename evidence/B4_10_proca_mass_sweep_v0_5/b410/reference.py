"""Independent reference computations (do not import B4.9 model/eigenmodes code).

* t17_roots:      characteristic roots c of det(c^2 A - c (B+B^T) + C), the T1.7
                  polynomial written from the Hessian formulas (no sympy).
* H_blocks:       the 4x4 Fourier matrix H(w,k) = w^2 H2 + w H1 + H0 written from
                  the B4.9 README formula.
* qz_roots:       the six finite roots of det H(w)=0 from the generalized eigenvalue
                  problem of the linearized pencil (scipy), two infinite roots dropped.
* mp_roots:       the same roots from the pole-cleared degree-six polynomial in
                  mpmath at high precision, built independently of B4.9.
* root_diagnostics: sigma_min, second smallest sigma, normalized residual and a
                  relative root condition number for each root.
"""
from __future__ import annotations
import numpy as np
import scipy.linalg as sla
import mpmath as mp


def t17_roots(p, s, kappa, Lam, w, q, v, r):
    A = np.array([[1 + p*q*q, -kappa*q*r], [-kappa*q*r, 1 + s*r*r]])
    B = np.array([[-kappa*r*v + 2*p*q*w, -kappa*q*v],
                  [-kappa*r*w, -kappa*q*w + 2*s*r*v]])
    C = np.array([[-1 - 3*Lam*q*q + p*w*w, -kappa*w*v],
                  [-kappa*w*v, -1 + s*v*v]])
    S = B + B.T
    P = np.polynomial.polynomial
    e = lambda i, j: np.array([C[i, j], -S[i, j], A[i, j]])
    poly = P.polysub(P.polymul(e(0, 0), e(1, 1)), P.polymul(e(0, 1), e(1, 0)))
    return np.sort_complex(P.polyroots(poly))


def H_blocks(M, k, gphi, gpsi, w, q, v, r, Lam):
    """(H0, H1, H2) with H(omega) = omega^2 H2 + omega H1 + H0 (README of B4.9)."""
    J0 = gphi*w*q + gpsi*v*r
    Abar = J0/M**2
    H0 = np.zeros((4, 4)); H1 = np.zeros((4, 4)); H2 = np.zeros((4, 4))
    H2[0, 0] = H2[1, 1] = H2[3, 3] = 1.0
    H0[0, 0] = -k**2*(1 + 3*Lam*q*q)
    H0[1, 1] = -k**2
    H1[0, 0] = -2*k*Abar*gphi
    H1[1, 1] = -2*k*Abar*gpsi
    # j_phi = gphi (w k - q omega), j_psi = gpsi (v k - r omega), symmetric entries
    for idx, (g, a, b) in enumerate(((gphi, w, q), (gpsi, v, r))):
        H0[idx, 3] = H0[3, idx] = g*a*k
        H1[idx, 3] = H1[3, idx] = -g*b
    H0[2, 2] = M**2 + k**2
    H1[2, 3] = H1[3, 2] = k
    H0[3, 3] = -M**2
    return H0, H1, H2


def H_of(omega, blocks):
    H0, H1, H2 = blocks
    return omega**2*H2 + omega*H1 + H0


def qz_roots(blocks):
    """Finite eigenvalues omega of the quadratic pencil (sorted by real part)."""
    H0, H1, H2 = blocks
    n = 4
    I = np.eye(n)
    Z = np.zeros((n, n))
    Aq = np.block([[Z, I], [-H0, -H1]])
    Bq = np.block([[I, Z], [Z, H2]])
    ev = sla.eig(Aq, Bq, right=False)
    ev = ev[np.isfinite(ev)]
    return ev[np.argsort(ev.real)]


def mp_roots(M, k, gphi, gpsi, w, q, v, r, Lam, dps=60):
    """Roots of the pole-cleared degree-six polynomial, evaluated in mpmath.

    P(om) = M^2 (M^2+k^2-om^2) D1 D2 + (M^2+k^2)(j1^2 D2 + j2^2 D1)
    """
    mp.mp.dps = dps
    M, k, gphi, gpsi, w, q, v, r, Lam = (mp.mpf(float(x)) for x in
                                         (M, k, gphi, gpsi, w, q, v, r, Lam))
    M2 = M*M
    Abar = (gphi*w*q + gpsi*v*r)/M2

    def pmul(a, b):
        out = [mp.mpf(0)]*(len(a) + len(b) - 1)
        for i, ai in enumerate(a):
            for j, bj in enumerate(b):
                out[i+j] += ai*bj
        return out

    def padd(a, b):
        n = max(len(a), len(b))
        a = a + [mp.mpf(0)]*(n - len(a)); b = b + [mp.mpf(0)]*(n - len(b))
        return [x + y for x, y in zip(a, b)]

    def pscale(a, c):
        return [c*x for x in a]

    D1 = [-k*k*(1 + 3*Lam*q*q), -2*k*Abar*gphi, mp.mpf(1)]
    D2 = [-k*k, -2*k*Abar*gpsi, mp.mpf(1)]
    j1 = [gphi*w*k, -gphi*q]
    j2 = [gpsi*v*k, -gpsi*r]
    pole = [M2 + k*k, mp.mpf(0), mp.mpf(-1)]
    t1 = pscale(pmul(pole, pmul(D1, D2)), M2)
    t2 = pscale(padd(pmul(pmul(j1, j1), D2), pmul(pmul(j2, j2), D1)), M2 + k*k)
    poly = padd(t1, t2)
    roots = mp.polyroots(poly[::-1], maxsteps=500, extraprec=2*dps)
    return sorted(roots, key=lambda z: (float(mp.re(z)), float(mp.im(z))))


def root_diagnostics(omega, blocks):
    """Per-root numerical health at a (possibly complex) root omega."""
    H0, H1, H2 = blocks
    H = H_of(complex(omega), blocks)
    U, s, Vh = np.linalg.svd(H)
    x = Vh[-1].conj()
    y = U[:, -1]
    resid = np.linalg.norm(H @ x)/max(1.0, s[0])
    dH = 2*omega*H2 + H1
    denom = abs(np.vdot(y, dH @ x))
    norms = [np.linalg.norm(Hj, 2) for Hj in (H0, H1, H2)]
    num = sum(abs(omega)**j*norms[j] for j in range(3))
    cond = num/(abs(omega)*denom) if abs(omega) > 0 and denom > 0 else np.inf
    return {"sigma_min": float(s[-1]), "sigma_2nd_smallest": float(s[-2]),
            "residual": float(resid), "root_condition": float(cond),
            "simple_root_ratio": float(s[-2]/max(s[-1], 1e-300))}
