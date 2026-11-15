"""Experiment 1: M_A -> infinity at fixed EFT coefficients (original B3 maximum-speed point)."""
from __future__ import annotations
import numpy as np
from ._b49 import full_roots, local_roots
from .spectrum import (b49_full_roots_checked, b49_leading_ratio, universal_c,
                       universal_omega_over_M, RootCountError)
from . import reference as ref
from .fixed_eft import (B3Point, background, couplings, eft_coefficients,
                        leading_shift_coefficients)


def match_to_b3(b3: B3Point, c_roots):
    """The four roots closest to the B3 roots (one-to-one), plus the two others."""
    c = np.asarray(c_roots, complex)
    used, picked = set(), []
    for c0 in sorted(b3.roots):
        j = min((i for i in range(len(c)) if i not in used), key=lambda i: abs(c[i] - c0))
        used.add(j); picked.append(c[j])
    others = [c[i] for i in range(len(c)) if i not in used]
    return np.array(picked), np.array(others)


def local_limit_table(b3: B3Point, masses):
    """Local (B4.9 EFT) roots at fixed-EFT couplings vs the B3 roots, for each M."""
    rows = []
    b3r = np.array(sorted(b3.roots))
    t17 = ref.t17_roots(*b3.args()).real
    for M in masses:
        lr = np.sort(local_roots(background(b3, M)).real)
        rows.append({"M": float(M), "max_abs_diff_vs_B3": float(np.abs(lr - b3r).max()),
                     "max_abs_diff_vs_independent_T17": float(np.abs(lr - t17).max())})
    return rows


def coupling_reconstruction(b3: B3Point, masses):
    rows = []
    for M in masses:
        bg = background(b3, M)
        p, s, k = eft_coefficients(bg)
        rows.append({"M": float(M), "g_phi": bg.g_phi, "g_psi": bg.g_psi,
                     "rel_err_p": abs(p - b3.p)/b3.p, "rel_err_s": abs(s - b3.s)/b3.s,
                     "rel_err_kappa": abs(k - b3.kappa)/abs(b3.kappa)})
    return rows


def scaling_collapse(b3: B3Point, xs, masses):
    """max |c(M, k=xM) - c(M_ref, k=xM_ref)|: exact collapse in x = k/M expected."""
    rows = []
    for x in xs:
        ref_c = None
        worst = 0.0
        for M in masses:
            k = x*M
            c = np.array(sorted(b49_full_roots_checked(background(b3, M), k)/k,
                                key=lambda z: (z.real, z.imag)))
            if ref_c is None:
                ref_c = c
            worst = max(worst, float(np.abs(c - ref_c).max()))
        rows.append({"x": float(x), "max_abs_dc_over_masses": worst})
    return rows


def b49_validity_map(b3: B3Point, masses, xs):
    """Where B4.9 full_roots returns all six roots, and how well the universal
    (M=1, QZ) spectrum reproduces it / mpmath elsewhere."""
    rows = []
    for M in masses:
        for x in xs:
            k = x*M
            bg = background(b3, M)
            n = len(full_roots(bg, k))
            uni = universal_c(b3, x)
            gphi, gpsi = couplings(b3, M)
            mpz = np.array([complex(z)/k for z in ref.mp_roots(M, k, gphi, gpsi, b3.w, b3.q, b3.v, b3.r, b3.Lam)])
            err = float(np.abs(np.sort_complex(uni) - np.sort_complex(mpz)).max()/max(1.0, np.abs(mpz).max()))
            rows.append({"M": float(M), "x": float(x), "k": float(k), "b49_n_roots": int(n),
                         "b49_leading_ratio": b49_leading_ratio(bg, k),
                         "universal_vs_mpmath_rel_err": err})
    return rows


def shift_convergence(b3: B3Point, xs, M=1.0):
    """Full c_i(x) vs B3 root c_i^(0) and the leading prediction c_i^(0) + A_i x^2."""
    c0, A = leading_shift_coefficients(b3)
    rows = []
    for x in xs:
        k = x*M
        c = full_roots(background(b3, M), k)/k
        four, _ = match_to_b3(b3, c)
        four = four[np.argsort(four.real)]
        for i in range(4):
            shift = four[i].real - c0[i]
            pred = A[i]*x*x
            rows.append({"x": float(x), "branch": i, "c_B3": float(c0[i]),
                         "c_full": float(four[i].real), "shift": float(shift),
                         "shift_leading_pred": float(pred),
                         "residual": float(shift - pred),
                         "max_abs_imag": float(abs(four[i].imag))})
    return rows


def verification_table(b3: B3Point, points):
    """Per-root comparison: B4.9 polynomial, QZ, mpmath, and numerical health."""
    out = []
    for (M, k) in points:
        gphi, gpsi = couplings(b3, M)
        args = (M, k, gphi, gpsi, b3.w, b3.q, b3.v, b3.r, b3.Lam)
        blocks = ref.H_blocks(*args)
        om_poly = np.array(sorted(full_roots(background(b3, M), k), key=lambda z: (z.real, z.imag)))
        om_qz = ref.qz_roots(blocks)
        om_qz = np.array(sorted(om_qz, key=lambda z: (z.real, z.imag)))
        om_mp = ref.mp_roots(*args)
        for i in range(len(om_poly)):
            d = ref.root_diagnostics(om_poly[i], blocks)
            mpz = complex(om_mp[i])
            out.append({"M": float(M), "k": float(k), "root": i,
                        "omega_poly_re": float(om_poly[i].real), "omega_poly_im": float(om_poly[i].imag),
                        "abs_diff_qz": float(abs(om_qz[i] - om_poly[i])),
                        "abs_diff_mpmath": float(abs(mpz - om_poly[i])),
                        **d})
    return out
