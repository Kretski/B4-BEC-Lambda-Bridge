"""B4.12: perturbativity / control parameters of the fixed-EFT Proca realisation (dimensional analysis + exact identities).

Fixed EFT coefficients p, s, kappa (rank one: kappa^2 = p s) and  g_phi = M sqrt(p), g_psi = M sqrt(s).  Exact identities:
    g_phi * |q|/M = sqrt(Xq),  g_phi * |w|/M = sqrt(Xw),  g_psi * |r|/M = sqrt(Xr),  g_psi * |v|/M = sqrt(Xv)   (M-independent)
i.e. (coupling) x (background gradient / Proca mass) = sqrt(X).  Equivalent: background gradient = sqrt(X) * Lambda_*, with the
Stueckelberg/operator strong-coupling scale Lambda_*phi = M/g_phi = 1/sqrt(p) (estimate by dimensional analysis, NOT a computed unitarity bound).
Fluctuation EFT parameter: x = k/M (corrections ~ X x^2).
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import linear_sum_assignment
from .fixed_eft import B3Point, couplings
from . import experiment2 as e2
from .fixed_eft import local_polynomial_coeffs

# thresholds proposed for B4.12 (set before looking at the tables)
SETS = {"A_g1_rho0.3": (1.0, 0.3), "B_g1_rho0.1": (1.0, 0.1), "C_gStrong_rho0.3": (float(np.sqrt(4*np.pi)), 0.3)}


def invariants(b3: B3Point, M: float):
    gphi, gpsi = couplings(b3, M)
    gpsi = abs(gpsi); gphi = abs(gphi)
    X = dict(Xq=b3.p*b3.q**2, Xr=b3.s*b3.r**2, Xw=b3.p*b3.w**2, Xv=b3.s*b3.v**2)
    prods = dict(Xq=gphi*abs(b3.q)/M, Xw=gphi*abs(b3.w)/M, Xr=gpsi*abs(b3.r)/M, Xv=gpsi*abs(b3.v)/M)
    return dict(g_phi=gphi, g_psi=gpsi, X=X, g_times_rho={k: prods[k] for k in prods},
                sqrtX={k: float(np.sqrt(v)) for k, v in X.items()},
                Lambda_star_phi=1/np.sqrt(b3.p), Lambda_star_psi=1/np.sqrt(b3.s))


def controlled(b3: B3Point):
    sx = max(np.sqrt(b3.p*b3.q**2), np.sqrt(b3.s*b3.r**2), np.sqrt(b3.p*b3.w**2), np.sqrt(b3.s*b3.v**2))
    out = dict(sqrtX_max=float(sx))
    for name, (gmax, rmax) in SETS.items():
        out[name] = bool(sx <= gmax*rmax)
    out["g_needed_for_rho0.3"] = float(sx/0.3); out["alpha_needed_for_rho0.3"] = float((sx/0.3)**2/(4*np.pi))
    return out


def x_window(b3: B3Point, tols=(0.01, 0.1), x_lo=1e-3, x_hi=1.0, n=400):
    """Smallest x where the four scalar-like speeds deviate from the local (B3) roots by more than tol (relative to max|c0|)."""
    c0 = np.sort(np.polynomial.polynomial.polyroots(local_polynomial_coeffs(b3)).real)
    out = {t: None for t in tols}
    for x in np.geomspace(x_lo, x_hi, n):
        nd = e2.evaluate(b3, x)
        c = np.array([b["c"] for b in nd["branches"]])
        cost = np.abs(c[None, :] - c0[:, None])
        ri, ci = linear_sum_assignment(cost)
        dev = cost[ri, ci].max()/np.abs(c0).max()
        for t in tols:
            if out[t] is None and dev > t:
                out[t] = float(x)
    return out
