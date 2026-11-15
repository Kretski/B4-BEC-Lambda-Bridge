"""Scale-robust spectrum for the fixed-EFT family, plus guards around B4.9.

B4.9 ``full_roots`` trims the leading polynomial coefficient when
|a6| <= 1e-13 * max|a_i| and then silently returns fewer than six roots.  For the
fixed-EFT family a6 = -M^2 while a0 ~ M^4 k^4, so this happens for large M k^2
(e.g. M=40,k=80 or M=1000,k=140).  B4.9 is left unchanged; B4.10 (i) checks the
root count and (ii) uses the exact scaling omega/M = F(k/M) to compute the
spectrum at M=1, where all coefficients are O(1)-O(10^2).
"""
from __future__ import annotations
import numpy as np
from ._b49 import full_roots, full_polynomial
from . import reference as ref
from .fixed_eft import B3Point, background, couplings


class RootCountError(RuntimeError):
    pass


def b49_full_roots_checked(bg, k):
    roots = full_roots(bg, k)
    if len(roots) != 6:
        raise RootCountError(f"B4.9 full_roots returned {len(roots)} of 6 roots "
                             f"(M={bg.M:g}, k={k:g}); leading coefficient trimmed")
    return roots


def b49_leading_ratio(bg, k):
    co = full_polynomial(bg, k).coef
    return float(abs(co[-1])/np.abs(co).max())


def universal_omega_over_M(b3: B3Point, x: float):
    """The six roots omega/M at k/M = x, from the linearized pencil at M=1 (QZ)."""
    gphi, gpsi = couplings(b3, 1.0)
    blocks = ref.H_blocks(1.0, float(x), gphi, gpsi, b3.w, b3.q, b3.v, b3.r, b3.Lam)
    return ref.qz_roots(blocks)


def universal_c(b3: B3Point, x: float):
    """c = omega/k at k/M = x (six roots, sorted by real part)."""
    return universal_omega_over_M(b3, x)/x
