"""Two-component mean-field BdG dispersion with quantum pressure.

The four signed branches are expanded as:
    omega_j(k) = c0_j*k + beta_j*k**3 + O(k**5).

No LHY correction is included in this module.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class BECParams:
    m1: float
    m2: float
    n1: float
    n2: float
    g11: float
    g22: float
    g12: float
    V1: float = 0.0
    V2: float = 0.0
    hbar: float = 1.0

    def validate(self):
        values = (
            self.m1, self.m2, self.n1, self.n2,
            self.g11, self.g22, self.g12,
            self.V1, self.V2, self.hbar,
        )
        if not all(np.isfinite(x) for x in values):
            raise ValueError("All parameters must be finite.")
        if min(self.m1, self.m2, self.n1, self.n2, self.hbar) <= 0:
            raise ValueError("Masses, densities and hbar must be positive.")
        if self.g11 <= 0 or self.g22 <= 0:
            raise ValueError("Intra-component couplings must be positive.")


def coefficients(params: BECParams):
    """Return a1, a2, b1, b2 and gamma."""
    params.validate()
    p = params

    a1 = p.g11 * p.n1 / p.m1
    a2 = p.g22 * p.n2 / p.m2
    b1 = p.hbar**2 / (4.0 * p.m1**2)
    b2 = p.hbar**2 / (4.0 * p.m2**2)
    gamma = p.g12**2 * p.n1 * p.n2 / (p.m1 * p.m2)

    return a1, a2, b1, b2, gamma


def hydrodynamic_polynomial(params: BECParams):
    """Polynomial H(c,0), ascending powers of c."""
    a1, a2, _, _, gamma = coefficients(params)
    p = params

    d1 = np.polynomial.polynomial.polymul(
        np.array([p.V1**2 - a1, -2*p.V1, 1.0]),
        np.array([p.V2**2 - a2, -2*p.V2, 1.0]),
    )
    d1[0] -= gamma
    return d1


def hydrodynamic_roots(params: BECParams, imag_tol=1e-9, gap_tol=1e-6):
    """Return four real, simple signed hydrodynamic roots."""
    roots = np.polynomial.polynomial.polyroots(
        hydrodynamic_polynomial(params)
    )

    if len(roots) != 4 or np.max(np.abs(roots.imag)) > imag_tol:
        raise ValueError("Hydrodynamic spectrum does not have four real roots.")

    roots = np.sort(roots.real)

    # repeated roots (e.g. the zero-sound-speed mode at g12^2 = g11*g22 and V=0) are split only by ~sqrt(eps) numerically
    scale_c = max(float(np.max(np.abs(roots))), 1e-300)
    if np.min(np.diff(roots)) < gap_tol * scale_c:
        raise ValueError("Repeated or nearly repeated hydrodynamic roots (relative gap below tolerance).")

    deriv = np.polynomial.polynomial.polyder(
        hydrodynamic_polynomial(params)
    )
    slopes = np.polynomial.polynomial.polyval(roots, deriv)

    scale = max(1.0, np.max(np.abs(slopes)))
    if np.any(np.abs(slopes) < 1e-10 * scale):
        raise ValueError("Repeated or ill-conditioned hydrodynamic root.")

    return roots


def beta_coefficients(params: BECParams):
    """Return (c0, beta), with omega = c0*k + beta*k**3 + O(k**5)."""
    p = params
    a1, a2, b1, b2, _ = coefficients(p)
    c0 = hydrodynamic_roots(p)

    d1 = (c0 - p.V1)**2 - a1
    d2 = (c0 - p.V2)**2 - a2

    deriv = np.polynomial.polynomial.polyder(
        hydrodynamic_polynomial(p)
    )
    h_c = np.polynomial.polynomial.polyval(c0, deriv)

    beta = (b1*d2 + b2*d1) / h_c
    return c0, beta


def characteristic(c, k, params: BECParams):
    """Evaluate the full finite-k characteristic polynomial H(c,k**2)."""
    a1, a2, b1, b2, gamma = coefficients(params)
    z = k*k
    d1 = (c - params.V1)**2 - a1 - b1*z
    d2 = (c - params.V2)**2 - a2 - b2*z
    return d1*d2 - gamma


def exact_phase_roots(k, params: BECParams):
    """Return the four signed finite-k phase velocities omega/k."""
    if k <= 0:
        raise ValueError("k must be positive.")

    a1, a2, b1, b2, gamma = coefficients(params)
    p = params
    z = k*k

    d1 = np.array([
        p.V1**2 - a1 - b1*z, -2*p.V1, 1.0
    ])
    d2 = np.array([
        p.V2**2 - a2 - b2*z, -2*p.V2, 1.0
    ])

    poly = np.polynomial.polynomial.polymul(d1, d2)
    poly[0] -= gamma
    roots = np.polynomial.polynomial.polyroots(poly)

    if np.max(np.abs(roots.imag)) > 1e-8:
        raise ValueError("Finite-k spectrum contains complex roots.")

    return np.sort(roots.real)


def galilean_shift(params: BECParams, u: float):
    """Apply a common Galilean boost to both background velocities."""
    from dataclasses import replace
    return replace(params, V1=params.V1 + u, V2=params.V2 + u)


def conditioning(params: BECParams):
    """Diagnostics for the leading expansion: |H'(c0)| per root and the relative gap to the nearest other root."""
    c0 = hydrodynamic_roots(params)
    deriv = np.polynomial.polynomial.polyder(hydrodynamic_polynomial(params))
    slopes = np.abs(np.polynomial.polynomial.polyval(c0, deriv))
    scale = max(float(np.max(np.abs(c0))), 1e-300)
    gaps = np.array([np.min(np.abs(np.delete(c0, i) - c0[i])) for i in range(4)]) / scale
    return {"abs_Hprime": slopes, "rel_gap": gaps}
