"""Quadratic Lagrangian, Hessian symbol, Krein value and canonical energy (Exp. 2).

Reconstructed quadratic Lagrangian L2 (1+1D, metric (-,+), fields
Phi = (dphi, dpsi, a0, a1)) as a list of bilinear terms
    coef * (d_{di} Phi_i) (d_{dj} Phi_j),   d in {none, t, x}:

  1/2 phi_t^2 + 1/2 psi_t^2 + 1/2 a1_t^2          kinetic (A0 has no time derivative)
  -1/2 (1+3 Lam q^2) phi_x^2 - 1/2 psi_x^2 + 1/2 a0_x^2
  - a1_t a0_x                                      from 1/2 F01^2, F01 = a1_t - a0_x
  Abar g_phi phi_t phi_x + Abar g_psi psi_t psi_x  contact (seagull) terms
  g_phi a1 (w phi_x + q phi_t) + g_psi a1 (v psi_x + r psi_t)   mediator couplings
  1/2 M^2 a0^2 - 1/2 M^2 a1^2                      Proca mass

For plane waves Phi = Re[z e^{i(kx - omega t)}] the cycle average is
    <L2> = (1/4) z^dagger Hcal(omega,k) z,       Hcal the Hessian symbol.
B4.9's real matrix H is related by a CONSTANT unitary: H_B49 = U^dagger Hcal U,
U = diag(1,1,u,u), u = i.  (Verified numerically in tests; H_B49 was NOT
re-derived from an independent action: L2 is reconstructed from it and the two
agree to rounding.)

Diagnostics (kept separate):
  Krein value  K_n  = z^dagger (dH/domega) z       (unitary invariant)
  shift-symmetric energy  E_n = Noether energy of L2 = sum_terms coef (n_t - 1) <d Phi_i d Phi_j>
                    = (1/4) omega z^dagger H' z   on shell (since <L2> = 0)
so sign(E_n) = sign(omega) * sign(K_n); the two coincide only for omega > 0.
NAMING (corrected in v0.3.1): L2 is the expansion about a background that is LINEAR IN TIME
(phi0 = w t + q x).  Its Noether energy is the second variation of H - mu_phi Q_phi - mu_psi Q_psi
(mu_phi=w, mu_psi=v), NOT the lab-frame second-order T00; the two differ by w Q2_phi + v Q2_psi
(b410/charges.py, derived from the nonlinear action and checked numerically).
A0 is non-dynamical: it is slaved to (phi, psi, a1) by its own equation (row 2
of H z = 0), so evaluating E on the null vector includes the constraint.
"""
from __future__ import annotations
import numpy as np

T_, X_, N_ = 't', 'x', None
_D = {N_: lambda om, k: 1.0, T_: lambda om, k: -1j*om, X_: lambda om, k: 1j*k}
IDX = {'phi': 0, 'psi': 1, 'a0': 2, 'a1': 3}


def terms(M, gphi, gpsi, w, q, v, r, Lam):
    """List of (coef, i, di, j, dj, group) bilinear terms of L2."""
    Abar = (gphi*w*q + gpsi*v*r)/M**2
    return [
        (0.5, 0, T_, 0, T_, 'phi'), (0.5, 1, T_, 1, T_, 'psi'), (0.5, 3, T_, 3, T_, 'A'),
        (-0.5*(1 + 3*Lam*q*q), 0, X_, 0, X_, 'phi'), (-0.5, 1, X_, 1, X_, 'psi'),
        (0.5, 2, X_, 2, X_, 'A'),
        (-1.0, 3, T_, 2, X_, 'A'),
        (Abar*gphi, 0, T_, 0, X_, 'phi'), (Abar*gpsi, 1, T_, 1, X_, 'psi'),
        (gphi*w, 3, N_, 0, X_, 'phiA'), (gphi*q, 3, N_, 0, T_, 'phiA'),
        (gpsi*v, 3, N_, 1, X_, 'psiA'), (gpsi*r, 3, N_, 1, T_, 'psiA'),
        (0.5*M**2, 2, N_, 2, N_, 'A'), (-0.5*M**2, 3, N_, 3, N_, 'A'),
    ]


def hessian_symbol(tl, omega, k, group=None):
    """Hermitian 4x4 matrix Hcal with <L2> = (1/4) z^dag Hcal z (optionally one group)."""
    H = np.zeros((4, 4), complex)
    for coef, i, di, j, dj, g in tl:
        if group is not None and g != group:
            continue
        a, b = _D[di](omega, k), _D[dj](omega, k)
        H[j, i] += coef*np.conj(b)*a
        H[i, j] += coef*np.conj(a)*b
    return H


def u_matrix(u=1j):
    return np.diag([1, 1, u, u]).astype(complex)


def b49_from_symbol(Hc, u=1j):
    U = u_matrix(u)
    return U.conj().T @ Hc @ U


def energy_density_terms(tl, z, omega, k, group=None):
    """Cycle-averaged Noether energy density, sum_terms coef (n_t - 1) <dPhi dPhi>."""
    E = 0.0
    for coef, i, di, j, dj, g in tl:
        if group is not None and g != group:
            continue
        n_t = (di == T_) + (dj == T_)
        a = _D[di](omega, k)*z[i]
        b = _D[dj](omega, k)*z[j]
        E += coef*(n_t - 1)*0.5*np.real(a*np.conj(b))
    return float(E)


def lagrangian_average(tl, z, omega, k):
    return sum(coef*0.5*np.real(_D[di](omega, k)*z[i]*np.conj(_D[dj](omega, k)*z[j]))
               for coef, i, di, j, dj, g in tl)


def energy_density_time_domain(tl, z, omega, k, n=4096):
    """Brute-force check: sample real fields on a phase grid, form T^00 directly."""
    th = np.linspace(0, 2*np.pi, n, endpoint=False)
    ph = np.exp(1j*th)[None, :]
    def fld(i, d):
        return np.real(_D[d](omega, k)*z[i]*ph)[0]
    E = np.zeros(n)
    for coef, i, di, j, dj, g in tl:
        n_t = (di == T_) + (dj == T_)
        E += coef*(n_t - 1)*fld(i, di)*fld(j, dj)
    return float(E.mean())


def null_vector_b49(H):
    """Right singular vector of the smallest singular value (real roots: real vector)."""
    U, s, Vh = np.linalg.svd(H)
    return Vh[-1].conj(), float(s[-1]), float(s[-2])


def krein_value(Hprime, x):
    x = x/np.linalg.norm(x)
    return float(np.real(np.vdot(x, Hprime @ x)))


def canonical_energy_b49(omega, Hprime, x):
    """E = (1/4) omega z^dag H' z with |z|=1 (unit-norm null vector)."""
    return 0.25*omega*krein_value(Hprime, x)


def group_energies(M, gphi, gpsi, w, q, v, r, Lam, omega, k, x):
    """Energy density split by term group for the (unit-norm) B49 null vector x."""
    tl = terms(M, gphi, gpsi, w, q, v, r, Lam)
    z = u_matrix() @ (x/np.linalg.norm(x))
    return {g: energy_density_terms(tl, z, omega, k, group=g)
            for g in ('phi', 'psi', 'A', 'phiA', 'psiA')}
