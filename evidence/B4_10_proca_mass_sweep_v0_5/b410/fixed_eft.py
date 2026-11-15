"""Fixed-EFT prescription for the B4.10 mass sweep.

Held fixed: p, s, kappa (low-energy EFT coefficients), Lambda and the scalar
background (w, q, v, r).  For each M_A > 0 the one-mediator couplings are

    g_phi(M) =  M sqrt(p),     g_psi(M) = -M kappa / sqrt(p),

with p = g_phi^2/M^2, s = g_psi^2/M^2, kappa = -g_phi g_psi/M^2  (rank one).
This is NOT the B4.9 fixed-coupling family.

Structural facts used (derived in README.md, verified numerically in tests):
  * the contact-term coefficients  a_phi = Abar g_phi = p w q - kappa v r  and
    a_psi = Abar g_psi = -kappa w q + s v r  do not depend on M;
  * dividing the pole-cleared polynomial by M^2 (M^2+k^2) k^4 gives
        Q(c; x) = F(c) - [c^2 x^2/(1+x^2)] d1(c) d2(c),   x = k/M,  c = omega/k,
    with F(c) = d1 d2 + p (w - q c)^2 d2 + s (v - r c)^2 d1, i.e. the T1.7 B3
    characteristic polynomial.  Hence the whole spectrum c_i = omega_i/k depends
    on x only.
"""
from __future__ import annotations
import json
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from ._b49 import Background

DATA_B3 = Path(__file__).resolve().parent.parent / "data" / "B3_max_speed_7_315096.json"


@dataclass(frozen=True)
class B3Point:
    p: float
    s: float
    kappa: float
    Lam: float
    w: float
    q: float
    v: float
    r: float
    roots: tuple
    record: dict

    def args(self):
        return (self.p, self.s, self.kappa, self.Lam, self.w, self.q, self.v, self.r)


def load_b3(path=DATA_B3) -> B3Point:
    rec = json.loads(Path(path).read_text())
    return B3Point(p=rec["p"], s=rec["s"], kappa=rec["kappa"], Lam=rec["Lambda"],
                   w=rec["w"], q=rec["q"], v=rec["v"], r=rec["r"],
                   roots=tuple(sorted(rec["roots_real"])), record=rec)


def couplings(b3: B3Point, M: float):
    if not (np.isfinite(M) and M > 0):
        raise ValueError("M must be finite and positive")
    return M*np.sqrt(b3.p), -M*b3.kappa/np.sqrt(b3.p)


def background(b3: B3Point, M: float) -> Background:
    gphi, gpsi = couplings(b3, M)
    return Background(M=M, g_phi=gphi, g_psi=gpsi, w=b3.w, q=b3.q,
                      v=b3.v, r=b3.r, Lambda=b3.Lam)


def eft_coefficients(bg: Background):
    """(p, s, kappa) reconstructed from the couplings of a Background."""
    M2 = bg.M**2
    return bg.g_phi**2/M2, bg.g_psi**2/M2, -bg.g_phi*bg.g_psi/M2


def contact_coefficients(b3: B3Point):
    """M-independent contact-term coefficients (Abar*g_phi, Abar*g_psi)."""
    return (b3.p*b3.w*b3.q - b3.kappa*b3.v*b3.r,
            -b3.kappa*b3.w*b3.q + b3.s*b3.v*b3.r)


def _polys(b3: B3Point):
    P = np.polynomial.polynomial
    aphi, apsi = contact_coefficients(b3)
    d1 = np.array([-(1 + 3*b3.Lam*b3.q**2), -2*aphi, 1.0])
    d2 = np.array([-1.0, -2*apsi, 1.0])
    s1 = np.array([b3.w, -b3.q])      # (w - q c)
    s2 = np.array([b3.v, -b3.r])      # (v - r c)
    F = P.polyadd(P.polymul(d1, d2),
                  P.polyadd(b3.p*P.polymul(P.polymul(s1, s1), d2),
                            b3.s*P.polymul(P.polymul(s2, s2), d1)))
    return d1, d2, F


def local_polynomial_coeffs(b3: B3Point):
    """F(c), ascending powers: the M-independent local (B3) characteristic polynomial."""
    return _polys(b3)[2]


def leading_shift_coefficients(b3: B3Point):
    """A_i in c_i(x) = c_i^(0) + A_i x^2 + O(x^4),  A_i = c_i^2 d1 d2 / F'(c_i)."""
    P = np.polynomial.polynomial
    d1, d2, F = _polys(b3)
    dF = P.polyder(F)
    c0 = np.array(sorted(b3.roots))
    A = c0**2*P.polyval(c0, d1)*P.polyval(c0, d2)/P.polyval(c0, dF)
    return c0, A


def domain_report(b3: B3Point):
    """The four dimensionless background invariants and related quantities.

    Recomputed from the parameters.  The T1.7.1 ``rank_one_grid`` text defines
    eft_domain as p q^2, s r^2 < max_normalized, p w^2 < max_w2, s v^2 < max_v2
    with defaults 0.95 (README: 0.5 for the strict example).  The original T1.7.3
    scan code that wrote the B3 row is not available, so applying the 0.95 cut is a
    reconstruction consistent with eft_domain=true, not a proof.
    """
    A = np.array([[1 + b3.p*b3.q**2, -b3.kappa*b3.q*b3.r],
                  [-b3.kappa*b3.q*b3.r, 1 + b3.s*b3.r**2]])
    aphi, apsi = contact_coefficients(b3)
    return {
        "Xq_p_q2": b3.p*b3.q**2, "Xr_s_r2": b3.s*b3.r**2,
        "Xw_p_w2": b3.p*b3.w**2, "Xv_s_v2": b3.s*b3.v**2,
        "Xw_plus_Xv": b3.p*b3.w**2 + b3.s*b3.v**2,
        "kinetic_matrix_A_eigenvalues": [float(x) for x in np.linalg.eigvalsh(A)],
        "rank_one_residual_ps_minus_kappa2": b3.p*b3.s - b3.kappa**2,
        "contact_coefficient_Abar_g_phi": aphi,
        "contact_coefficient_Abar_g_psi": apsi,
        "max_abs_c": max(abs(c) for c in b3.roots),
        "note": "eft_domain reconstructed from T1.7.1 rank_one_grid (cut 0.95); T1.7.3 scan code unavailable",
    }
