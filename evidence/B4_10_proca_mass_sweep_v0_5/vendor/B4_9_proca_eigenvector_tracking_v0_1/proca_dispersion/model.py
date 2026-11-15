"""Local EFT and full Proca dispersion relations (1+1 dimensions).

Conventions: metric (-,+), Fourier exp(-i omega t + i k x).
The scalar base symbol is diag(D_phi, D_psi), where
D_phi=omega^2-k^2(1+3 Lambda q0^2), D_psi=omega^2-k^2.
The background current perturbation is
delta J = i [g_phi(w0 k-q0 omega) Phi + g_psi(v0 k-r0 omega) Psi].
"""
from dataclasses import dataclass
import numpy as np
from numpy.polynomial import Polynomial as Poly

@dataclass(frozen=True)
class Background:
    M: float
    g_phi: float
    g_psi: float
    w: float
    q: float
    v: float
    r: float
    Lambda: float = 0.0

    def validate(self):
        vals = (self.M, self.g_phi, self.g_psi, self.w, self.q,
                self.v, self.r, self.Lambda)
        if not all(np.isfinite(x) for x in vals):
            raise ValueError("All background values must be finite.")
        if self.M <= 0 or self.Lambda < 0:
            raise ValueError("Require M>0 and Lambda>=0.")

def _trim(coeff, tol=1e-13):
    coeff = np.asarray(coeff, dtype=float)
    while len(coeff) > 1 and abs(coeff[-1]) <= tol * max(1.0, np.max(np.abs(coeff))):
        coeff = coeff[:-1]
    return coeff

def _real_roots(coeff, imag_tol=2e-7):
    coeff = _trim(coeff)
    roots = np.polynomial.polynomial.polyroots(coeff)
    roots = sorted(roots, key=lambda z: (float(np.real(z)), float(np.imag(z))))
    return np.asarray(roots, dtype=complex)

def local_polynomial(bg: Background):
    """Polynomial det M_EFT(c), ascending powers of c."""
    bg.validate()
    M2 = bg.M**2
    # a_phi(c)=g_phi(w-q*c), a_psi(c)=g_psi(v-r*c)
    c = Poly([0.0, 1.0])
    a1 = bg.g_phi * (bg.w - bg.q*c)
    a2 = bg.g_psi * (bg.v - bg.r*c)
    d1 = c*c - (1.0 + 3.0*bg.Lambda*bg.q**2)
    d2 = c*c - 1.0
    # The UV background Proca field generates a contact/seagull term
    # Abar_i * delta^2 J_i. In the characteristic symbol this is
    # -2*c*(J0/M^2)*g_i on each scalar diagonal.
    J0 = bg.g_phi*bg.w*bg.q + bg.g_psi*bg.v*bg.r
    Abar = J0/M2
    d1 = d1 - 2.0*c*Abar*bg.g_phi
    d2 = d2 - 2.0*c*Abar*bg.g_psi
    return d1*d2 + (a1*a1*d2 + a2*a2*d1)/M2

def local_roots(bg: Background):
    return _real_roots(local_polynomial(bg).coef)

def full_polynomial(bg: Background, k: float):
    """Polynomial in omega for the full coupled scalar-Proca spectrum.

    It is the pole-cleared determinant:
      M^2(M^2+k^2-omega^2) D1 D2
      +(M^2+k^2) [j1^2 D2+j2^2 D1] = 0,
    where j1=g_phi(w*k-q*omega), j2=g_psi(v*k-r*omega).
    It retains the Proca branch and is degree six generically.
    """
    bg.validate()
    if not np.isfinite(k) or k <= 0:
        raise ValueError("Require finite k>0.")
    om = Poly([0.0, 1.0])
    M2 = bg.M**2
    J0 = bg.g_phi*bg.w*bg.q + bg.g_psi*bg.v*bg.r
    Abar = J0/M2
    D1 = om*om - k*k*(1.0 + 3.0*bg.Lambda*bg.q**2) - 2.0*om*k*Abar*bg.g_phi
    D2 = om*om - k*k - 2.0*om*k*Abar*bg.g_psi
    j1 = bg.g_phi*(bg.w*k - bg.q*om)
    j2 = bg.g_psi*(bg.v*k - bg.r*om)
    return (M2*(M2+k*k-om*om)*D1*D2
            + (M2+k*k)*(j1*j1*D2 + j2*j2*D1))

def full_roots(bg: Background, k: float):
    return _real_roots(full_polynomial(bg, k).coef)

def local_symbol_matrix(bg: Background, c: complex):
    a = np.array([bg.g_phi*(bg.w-bg.q*c),
                  bg.g_psi*(bg.v-bg.r*c)], dtype=complex)
    J0 = bg.g_phi*bg.w*bg.q + bg.g_psi*bg.v*bg.r
    Abar = J0/bg.M**2
    base = np.diag([c*c-(1+3*bg.Lambda*bg.q**2)-2*c*Abar*bg.g_phi,
                    c*c-1-2*c*Abar*bg.g_psi]).astype(complex)
    return base + np.outer(a, a)/bg.M**2

def full_symbol_matrix(bg: Background, omega: complex, k: float):
    c = omega/k
    a = np.array([bg.g_phi*(bg.w-bg.q*c),
                  bg.g_psi*(bg.v-bg.r*c)], dtype=complex)
    J0 = bg.g_phi*bg.w*bg.q + bg.g_psi*bg.v*bg.r
    Abar = J0/bg.M**2
    base = np.diag([omega**2-k**2*(1+3*bg.Lambda*bg.q**2)-2*omega*k*Abar*bg.g_phi,
                    omega**2-k**2-2*omega*k*Abar*bg.g_psi]).astype(complex)
    denom = bg.M**2*(bg.M**2+k**2-omega**2)
    if abs(denom) < 1e-13:
        raise ZeroDivisionError("At the Proca pole use the pole-cleared polynomial.")
    return base + (k*k/denom)*(bg.M**2+k*k)*np.outer(a,a)

def pole_frequency(bg: Background, k: float):
    return float(np.sqrt(bg.M**2+k**2))
