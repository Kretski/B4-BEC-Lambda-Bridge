"""Leading k^2 corrections of the four scalar-like roots of the B4 Proca spectrum.

Exact spectrum (derived from the action in theory/derive_from_action.py, det K = -M^8 x^4 (1+x^2) Q):
    Q(c; x) = F(c) - x^2/(1+x^2) c^2 d1(c) d2(c),  x = k/M,  c = omega/k
    d1 = c^2 - 2 a_phi c - (1 + 3 Lam q^2),  d2 = c^2 - 2 a_psi c - 1
    F  = d1 d2 + p (w - q c)^2 d2 + s (v - r c)^2 d1
    a_phi = p w q - kappa v r,  a_psi = -kappa w q + s v r
Definitions:
    c_i(x) = c_i0 + A_i x^2 + O(x^4),            A_i = c_i0^2 d1 d2 / F'(c_i0)
    omega_i = c_i0 k + A_i k^3 / M^2 + ...;      M^2 Lam_eff_i = 2 A_i / c_i0   (a *definition* by analogy with
                                                  omega = c k sqrt(1+Lam k^2); root dependent; NOT bare Lam)
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import numpy as np
P = np.polynomial.polynomial


@dataclass(frozen=True)
class Params:
    p: float; s: float; kappa: float; Lam: float
    w: float; q: float; v: float; r: float

    @classmethod
    def from_row(cls, row):
        return cls(row["p"], row["s"], row["kappa"], row["Lambda"], row["w"], row["q"], row["v"], row["r"])

    @classmethod
    def from_json(cls, path):
        return cls.from_row(json.load(open(path)))


def polys(P_: Params):
    a_phi = P_.p*P_.w*P_.q - P_.kappa*P_.v*P_.r
    a_psi = -P_.kappa*P_.w*P_.q + P_.s*P_.v*P_.r
    d1 = np.array([-(1 + 3*P_.Lam*P_.q**2), -2*a_phi, 1.0])
    d2 = np.array([-1.0, -2*a_psi, 1.0])
    s1 = np.array([P_.w, -P_.q]); s2 = np.array([P_.v, -P_.r])
    F = P.polyadd(P.polymul(d1, d2), P.polyadd(P_.p*P.polymul(P.polymul(s1, s1), d2),
                                               P_.s*P.polymul(P.polymul(s2, s2), d1)))
    return d1, d2, F


def local_roots(P_: Params):
    return np.sort_complex(P.polyroots(polys(P_)[2]))


def A_coefficients(P_: Params):
    """(c0, A) for the four local roots; requires all roots real and simple."""
    d1, d2, F = polys(P_)
    c0 = P.polyroots(F)
    if np.abs(c0.imag).max() > 1e-9:
        raise ValueError("local roots are complex")
    c0 = np.sort(c0.real)
    A = c0**2*P.polyval(c0, d1)*P.polyval(c0, d2)/P.polyval(c0, P.polyder(F))
    return c0, A


def lambda_eff_M2(P_: Params):
    c0, A = A_coefficients(P_)
    return c0, A, 2*A/c0


def exact_roots(P_: Params, x: float):
    """Six roots c of Q(c; x) (degree 6): the four scalar-like ones and the two heavy ones."""
    d1, d2, F = polys(P_)
    eps = x*x/(1 + x*x)
    Q = P.polysub(F, eps*P.polymul(P.polymul([0, 0, 1.0], d1), d2))
    return P.polyroots(Q)


def tracked_scalar_roots(P_: Params, x: float):
    """Four exact roots assigned one-to-one to the local roots (Hungarian assignment on |c - c0|)."""
    from scipy.optimize import linear_sum_assignment
    c0, _ = A_coefficients(P_)
    r = exact_roots(P_, x)
    cost = np.abs(c0[:, None] - r[None, :])
    ri, ci = linear_sum_assignment(cost)
    return c0, r[ci[np.argsort(ri)]]


def beta_M2(P_: Params):
    """Primary coefficient: omega_i = c_i0 k + beta_i k^3, beta_i = A_i / M^2.  Returns A_i (= M^2 beta_i).
    beta is unchanged by a Galilean shift omega -> omega + u k; Lam_eff = 2 beta / c0 is not (c0 shifts)."""
    return A_coefficients(P_)[1]


def conditioning(P_: Params):
    """|F'(c_i0)| and the relative gap of each local root to its nearest neighbour (A_i ~ 1/F')."""
    d1, d2, F = polys(P_)
    c0 = A_coefficients(P_)[0]
    dF = np.abs(P.polyval(c0, P.polyder(F)))
    scale = max(np.abs(c0).max(), 1e-300)
    gap = np.array([np.min(np.abs(np.delete(c0, i) - c0[i])) for i in range(4)])/scale
    return dict(abs_Fprime=dF, rel_gap=gap, near_zero_root=np.abs(c0) < 0.05*scale)
