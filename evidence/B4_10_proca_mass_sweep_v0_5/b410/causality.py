"""B4.11: front velocity of the full Proca system.

High-frequency (principal) part: with omega = c k, k -> infinity at fixed M, the homogeneous degree-6 part of the
pole-cleared polynomial P(omega,k) = M^2 (M^2+k^2-omega^2) D1 D2 + (M^2+k^2)(j1^2 D2 + j2^2 D1) is
        k^6 M^2 [ (1-c^2) d1 d2 + p (w-qc)^2 d2 + s (v-rc)^2 d1 ]  =  k^6 M^2 Q_inf(c)  (= -det H, so det H -> -M^2 k^6 Q_inf),
Q_inf = F - c^2 d1 d2 (F = local B3/T1.7 polynomial), independent of M (fixed-EFT prescription).
Front velocity of a hyperbolic system = max |root of the principal characteristic polynomial|.
Complex roots of Q_inf => Im omega ~ k at high k => the Cauchy problem is ill posed (no front velocity).
"""
from __future__ import annotations
import numpy as np
import mpmath as mp
from . import experiment2 as e2, reference as ref, audit_hamiltonian as ah
from .fixed_eft import B3Point, couplings, local_polynomial_coeffs, _polys


def principal_roots(b3: B3Point):
    c, im = e2.c_infinity(b3)
    P = np.polynomial.polynomial
    d1, d2, F = _polys(b3)
    Q = P.polysub(F, P.polymul(P.polymul([0, 0, 1.0], d1), d2))
    allr = P.polyroots(Q)
    return allr, bool(np.abs(allr.imag).max() < 1e-9)


def local_roots(b3: B3Point):
    r = np.polynomial.polynomial.polyroots(local_polynomial_coeffs(b3))
    return r, bool(np.abs(r.imag).max() < 1e-9)


def det_ratio(b3: B3Point, cs, k=1e5, dps=60):
    """det H(c k, k) / (k^6 Q_inf(c)) in mpmath: a single constant (-M^2 = -1) for every c proves that the
    homogeneous degree-6 part of det H is M^2 k^6 Q_inf(c)."""
    mp.mp.dps = dps
    g = couplings(b3, 1.0)
    d1, d2, F = _polys(b3)
    P = np.polynomial.polynomial
    Q = P.polysub(F, P.polymul(P.polymul([0, 0, 1.0], d1), d2))
    blocks = ref.H_blocks(1.0, k, *g, b3.w, b3.q, b3.v, b3.r, b3.Lam)
    out = []
    for c in cs:
        om = mp.mpf(c)*mp.mpf(k)
        H = mp.matrix(4, 4)
        for a in range(4):
            for b in range(4):
                H[a, b] = mp.mpf(blocks[2][a, b])*om**2 + mp.mpf(blocks[1][a, b])*om + mp.mpf(blocks[0][a, b])
        out.append(float(mp.det(H)/(mp.mpf(k)**6*P.polyval(c, Q))))
    return out


def evolve_front(b3: B3Point, sigma=0.05, L=200.0, N=2**15, times=(4.0, 8.0), thr=1e-7, M=1.0):
    """Pseudo-spectral time evolution of the reduced Hamiltonian system; returns front positions and max amplitude."""
    g = couplings(b3, M)
    P = (M, *g, b3.w, b3.q, b3.v, b3.r, b3.Lam)
    x = (np.arange(N) - N//2)*(L/N)
    ks = 2*np.pi*np.fft.fftfreq(N, d=L/N)
    u0 = np.zeros((6, N)); u0[0] = np.exp(-x**2/(2*sigma**2))
    u0k = np.fft.fft(np.fft.ifftshift(u0, axes=1), axis=1)
    J = np.block([[np.zeros((3, 3)), np.eye(3)], [-np.eye(3), np.zeros((3, 3))]])
    Gm = np.array([J @ ah.hamiltonian_matrix(k, *P) for k in ks])
    lam, V = np.linalg.eig(Gm)
    c0 = np.linalg.solve(V, u0k.T[..., None])[..., 0]
    out = {}
    for t in times:
        ck = c0*np.exp(lam*t)
        uk = np.einsum('nij,nj->ni', V, ck)
        u = np.real(np.fft.fftshift(np.fft.ifft(uk.T, axis=1), axes=1))
        amp = np.abs(u).max(axis=0)
        above = np.where(amp > thr)[0]
        out[t] = dict(right=float(x[above.max()]) if len(above) else None, left=float(x[above.min()]) if len(above) else None,
                      max_amp=float(amp.max()), growth=float(np.log(max(amp.max(), 1e-300))))
    return out, float(np.abs(lam.real).max())
