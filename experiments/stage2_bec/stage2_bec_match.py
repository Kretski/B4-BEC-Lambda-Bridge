#!/usr/bin/env python3
"""BEC Stage 2: B4 targets versus the registered two-component BEC family.

Implements docs/STAGE2_PREREGISTRATION.md (v2; blob 39f66ce2, sha256 8cf6dd9a...) together with
STAGE2_PREREGISTRATION_PROCEDURAL_ANNEX_v0.7.md (sha256 71e1edc1...). Section labels in comments (A1..K, v2) refer to
those two documents. Script version v0.4.1:
  v0.2: R2 removed the zero-momentum sensitivity; B5 outcome field; R10 finiteness; R4/C4 separate
        rejection counters; SciPy pinned for B2a.
  v0.3: R4 exception classification by message, all exceptions caught inside b4_spectrum;
        positive control (annex J) run before the real targets; J4 applied to the verdicts.
  v0.4: R10 -- Residual.evaluate returns an explicit status; a non-finite final residual is an
        ABORT in check_final and is no longer merged with an ordinary inadmissible final point.
  v0.4.1: metadata only -- docstring version; start records without an admissible start carry
        final_status "no_admissible_start" (uniform report schema). No change in logic.

STATUS: written for static review only. Execution is NOT authorised (annex §G, steps 6-7).

Usage (only after registration and a separate authorisation):
    python3 -I stage2_bec_match.py --repo <fresh clean LF clone> --out <empty output dir>

Order of work inside one run (fixed):
    0. verify the sha256 of the four input files (E) and the pinned SciPy trf implementation (B2a);
       ABORT on mismatch, before any computation
    1. load theory/lambda_eff.py and theory/bec_dispersion.py from those verified files
    2. select the targets (E), check C1c on every target, write targets.json (+ hash); ABORT if count != 29
    2b. positive control (annex J), per region; write positive_control.json and record its hash
    3. per target and region (A primary, B control): search (B), criteria (v2), robustness (C), J4
    4. sign-tuple sampling (D)
    5. write results JSON, SHA256SUMS
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import os
import platform
import sys
import time

import numpy as np

# ----------------------------------------------------------------------------------------------
# Fixed constants (v2 and annex v0.4). Nothing below may be changed without a new registration.
# ----------------------------------------------------------------------------------------------
EXPECTED_SHA256 = {
    "data/e3_manifest.json": "1bae892a08c3de1f928e4057ee0ad801294a90462478aaedf5b9c724e23d87e0",
    "data/B3_max_speed_7_315096.json": "99a6b04bc34821ddd5c5d58a78b3474d369aa9fc0ee264163096398b0028ccbd",
    "theory/lambda_eff.py": "e07299762f3e6b89187f8dc8270b430c951e15816ea129983f040dd2f80c951f",
    "theory/bec_dispersion.py": "228086f4a1a3ab292c1798910dffcc4ba69cec9e17c88dd7046369f52b878865",
}
EXPECTED_N_TARGETS = 29                 # v2 / annex E

# B2a rests on the step-acceptance rule of scipy's trf ("if actual_reduction > 0: x = x_new"), read in
# this exact implementation. Version and file hash are pinned; ABORT otherwise (static review, item 5).
EXPECTED_SCIPY_VERSION = "1.18.1"
EXPECTED_TRF_SHA256 = "0c7c29d5d3fd6915720e9547f051d47e14ca5a821e7872686aa8cf1d1c4e4893"

SEED = 20261111                         # v2
N_STARTS = 80                           # v2
MAX_START_REDRAWS = 1000                # annex B4
REGION_CODE = {"A": 0, "B": 1}          # annex B4, D

# v2 classification thresholds (applied to lambda_min of the normalised form E~, annex A1)
LAM_STABLE = -1e-9                      # class A: lambda_min >= -1e-9
LAM_UNSTABLE = -1e-6                    # class B: lambda_min <  -1e-6
HPRIME_MIN = 1e-2                       # v2 class 4: |H'| < 1e-2 at a root (units a1 = 1)
# relative gap < 1e-6 and complex roots: handled by bec_dispersion.hydrodynamic_roots (annex A5)

# annex A3: k grid in K = k sqrt(b1/a1): 0 plus 200 geometric points in [1e-3, 1e2]
K_GRID = np.concatenate(([0.0], np.geomspace(1e-3, 1e2, 200)))
GRID_TOL = 1e-12                        # ABORT if lambda_min(k) < lambda_min(0) - 1e-12

# v2 criteria
EPS0 = 0.05
TH_SIGN = 0.05                          # (1): only branches with |bhat| >= 0.05 on both sides
TH_EPS = 0.25                           # (2): max_j eps_j <= 0.25
TH_X = 0.05                             # (3): |x_j^B4 - x_j^BEC| <= 0.05 for the two inner roots
N_PERT = 50                             # (4): 50 perturbed targets
MAX_PERT_DRAWS = 500                    # annex C4
TH_ROBUST = 0.80                        # (4): >= 80 %
PERT_REL, PERT_ABS = 0.02, 0.005        # annex C2
PERT_KEYS = ("w", "q", "v", "r", "p", "s", "Lambda")   # annex C1, fixed draw order

# annex B1: search variables y = (log10 rho, f, dV, log10 kappa_b)
Y_LO = np.array([-2.0, 0.0, -4.0, -2.0])
Y_HI = np.array([2.0, 0.99, 4.0, 2.0])

# annex B2
PENALTY_VALUE = 1e3
N_RES = 6
PENALTY_COST = 0.5 * N_RES * PENALTY_VALUE ** 2         # 3e6
FEASIBLE_COST_BOUND = 0.5 * (2 * (1 / TH_X) ** 2 + 4 * (2 / TH_EPS) ** 2)   # 528 (B2a)

# annex B3
LSQ_KW = dict(method="trf", x_scale="jac", ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=2000)

# annex J (positive control)
N_PC_TARGETS = 29                       # J1: same number as the real targets
MAX_PC_DRAWS = 10_000                   # J1
PC_INDEX_OFFSET = 1000                  # J2: start seeds use background_index = 1000 + i
PC_MIN_RECOVERED = 27                   # J4: >= 27 of 29 (rate >= 0.90), per region

# annex D
N_SAMPLE_PER_REGION = 200_000
MAX_SAMPLE_DRAWS = 20_000_000

# annex E conditioning labels (B4 side)
COND_REL_GAP = 1e-6
COND_FPRIME = 1e-2
RANK1_RTOL = 1e-12                      # annex C1c
COMPLEX_ROOTS_MESSAGE = "local roots are complex"   # theory/lambda_eff.py, A_coefficients (registered text)


class Abort(RuntimeError):
    """Any condition the annex classifies as ABORT."""


# ----------------------------------------------------------------------------------------------
# 0-1. Input verification and module loading (annex E)
# ----------------------------------------------------------------------------------------------
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_inputs(repo):
    seen = {}
    for rel, want in EXPECTED_SHA256.items():
        path = os.path.join(repo, rel)
        if not os.path.isfile(path):
            raise Abort(f"input missing: {rel}")
        got = sha256_file(path)
        if got != want:
            raise Abort(f"sha256 mismatch for {rel}: {got} != {want}")
        seen[rel] = got
    return seen


def verify_scipy():
    """Pin of the solver implementation used in the B2a argument. ABORT on mismatch."""
    import scipy
    import scipy.optimize._lsq.trf as trf
    if scipy.__version__ != EXPECTED_SCIPY_VERSION:
        raise Abort(f"scipy {scipy.__version__} != pinned {EXPECTED_SCIPY_VERSION} (B2a)")
    got = sha256_file(trf.__file__)
    if got != EXPECTED_TRF_SHA256:
        raise Abort(f"scipy trf.py sha256 {got} != pinned {EXPECTED_TRF_SHA256} (B2a)")
    return {"scipy": scipy.__version__, "trf_py": trf.__file__, "trf_sha256": got}


def load_module(name, path):
    """Load a module from the exact verified file (not from sys.path)."""
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ----------------------------------------------------------------------------------------------
# BEC side: family -> parameters (annex B1, A2), classification (A1, A3, A4, A5, v2)
# ----------------------------------------------------------------------------------------------
def family_to_bec(bec, y):
    """y = (log10 rho, f, dV, log10 kappa_b) -> BECParams.

    Fixed units (annex B1): a1 = 1, b1 = 1/4, m1 = hbar = 1, hence g11 n1 = 1.
    b2 = kappa_b * b1 -> m2 = 1/sqrt(kappa_b);  a2 = rho -> g22 n2 = rho m2;
    gamma = g12^2 n1 n2/(m1 m2) = f^2 a1 a2 with g12 = f sqrt(g11 g22).
    Densities: n1 = n2 = 1 is a bookkeeping choice only: a_i, b_i, gamma and the symmetric-frame
    V_i do not depend on it, and nothing computed here depends on n2/n1. (The zero-momentum
    sensitivity, which would depend on it, was removed: static review R2, annex v0.5 A2.)
    Frame (annex A2): symmetric only, V1 = +dV sqrt(a1)/2, V2 = -dV sqrt(a1)/2.
    """
    log_rho, f, dV, log_kb = (float(v) for v in y)
    rho, kb = 10.0 ** log_rho, 10.0 ** log_kb
    m1, m2 = 1.0, 1.0 / math.sqrt(kb)
    n1 = n2 = 1.0
    g11 = 1.0 * m1 / n1                 # a1 = 1
    g22 = rho * m2 / n2                 # a2 = rho
    g12 = f * math.sqrt(g11 * g22)
    a1 = 1.0
    delta = dV * math.sqrt(a1)          # V1 - V2
    V1, V2 = +delta / 2, -delta / 2     # annex A2: symmetric frame
    return bec.BECParams(m1=m1, m2=m2, n1=n1, n2=n2, g11=g11, g22=g22, g12=g12, V1=V1, V2=V2, hbar=1.0)


def energy_matrices(bec, P):
    """Normalised lab-frame energy form E~(k) on the K grid (annex A1), shape (len(K_GRID), 4, 4)."""
    a1, a2, b1, b2, _ = bec.coefficients(P)
    f = P.g12 / math.sqrt(P.g11 * P.g22)
    # K = k sqrt(b1/a1)  ->  b1 k^2/a1 = K^2,  b2 k^2/a2 = (b2/b1)(a1/a2) K^2
    K2 = K_GRID ** 2
    E = np.zeros((len(K_GRID), 4, 4))
    E[:, 0, 0] = 1.0 + K2
    E[:, 2, 2] = 1.0 + (b2 / b1) * (a1 / a2) * K2
    E[:, 1, 1] = 1.0
    E[:, 3, 3] = 1.0
    E[:, 0, 1] = E[:, 1, 0] = P.V1 / math.sqrt(a1)
    E[:, 2, 3] = E[:, 3, 2] = P.V2 / math.sqrt(a2)
    E[:, 0, 2] = E[:, 2, 0] = f
    return E


def energy_lambda(bec, P):
    """(lambda_min over the grid, lambda_min at K = 0, positive-definite verdict at K = 0).

    ABORT if the grid violates the monotonicity of annex A3 (implementation error)."""
    E = energy_matrices(bec, P)
    if not np.all(np.isfinite(E)):
        return None
    lam = np.linalg.eigvalsh(E)[:, 0]
    lam0 = float(lam[0])
    if np.any(lam < lam0 - GRID_TOL):
        raise Abort(f"A3 violated: grid lambda_min {lam.min():.3e} < lambda_min(0) {lam0:.3e}")
    pd0 = bool(lam0 > 0.0)              # inertia-only verdict (annex A1 control count)
    return float(lam.min()), lam0, pd0


def classify(bec, P):
    """v2 classification with annex A1/A4/A5. Returns dict with 'cls' in {1, 2, 3, 4} and diagnostics.

    1 = region A, 2 = region B, 3 = dynamically unstable (excluded), 4 = boundary/ill-conditioned (excluded).
    """
    out = {"cls": None, "lam_min": None, "lam0": None, "pd0": None, "reason": ""}
    try:
        roots = bec.hydrodynamic_roots(P)
    except ValueError as exc:
        msg = str(exc)
        if "four real roots" in msg:
            out.update(cls=3, reason="complex roots")
        elif "Repeated" in msg or "ill-conditioned" in msg:
            out.update(cls=4, reason="repeated/ill-conditioned root (module)")
        else:                           # e.g. validate(): non-finite parameters
            out.update(cls=4, reason=f"module error: {msg}")
        return out
    dH = np.polynomial.polynomial.polyder(bec.hydrodynamic_polynomial(P))
    hprime = np.abs(np.polynomial.polynomial.polyval(roots, dH))
    if not np.all(np.isfinite(hprime)) or np.min(hprime) < HPRIME_MIN:
        out.update(cls=4, reason="|H'| < 1e-2 (v2)")
        return out
    e = energy_lambda(bec, P)
    if e is None:
        out.update(cls=4, reason="non-finite energy form (A4)")
        return out
    lam_min, lam0, pd0 = e
    out.update(lam_min=lam_min, lam0=lam0, pd0=pd0)
    if lam_min >= LAM_STABLE:
        out.update(cls=1, reason="A")
    elif lam_min < LAM_UNSTABLE:
        out.update(cls=2, reason="B")
    else:
        out.update(cls=4, reason="lambda_min in excluded band (v2)")
    return out


def bec_spectrum(bec, P):
    """(c0 sorted, beta) or None."""
    try:
        c0, beta = bec.beta_coefficients(P)
    except ValueError:
        return None
    if not (np.all(np.isfinite(c0)) and np.all(np.isfinite(beta))):
        return None
    return np.asarray(c0, float), np.asarray(beta, float)


# ----------------------------------------------------------------------------------------------
# Unit-free comparison (v2)
# ----------------------------------------------------------------------------------------------
def x_structure(c):
    c = np.sort(np.asarray(c, float))
    span = c[-1] - c[0]
    if not span > 0:
        return None
    return (c - c[0]) / span


def bhat(beta):
    beta = np.asarray(beta, float)
    m = np.max(np.abs(beta))
    if not m > 0:
        return None
    return beta / m


def eps_vec(bh4, bhb):
    den = np.maximum(np.maximum(np.abs(bh4), np.abs(bhb)), EPS0)
    return np.abs(bh4 - bhb) / den


def criteria(x4, bh4, xb, bhb, perm=(0, 1, 2, 3)):
    """v2 criteria (1)-(3). perm[j] = BEC branch assigned to B4 branch j (primary: identity)."""
    perm = np.asarray(perm)
    xbp, bhbp = xb[perm], bhb[perm]
    both = (np.abs(bh4) >= TH_SIGN) & (np.abs(bhbp) >= TH_SIGN)
    c1 = bool(np.all(np.sign(bh4[both]) == np.sign(bhbp[both])))
    eps = eps_vec(bh4, bhbp)
    c2 = bool(np.max(eps) <= TH_EPS)
    dx = np.abs(x4[[1, 2]] - xbp[[1, 2]])
    c3 = bool(np.max(dx) <= TH_X)
    return {"c1_sign": c1, "c2_magnitude": c2, "c3_roots": c3,
            "max_eps": float(np.max(eps)), "max_dx_inner": float(np.max(dx))}


# ----------------------------------------------------------------------------------------------
# B4 side: targets (annex E, C1c)
# ----------------------------------------------------------------------------------------------
KEYS = ("p", "s", "kappa", "Lambda", "w", "q", "v", "r")


def target_record(row, label):
    return {"label": label, **{k: float(row[k]) for k in KEYS}}


def c1c_check(t):
    """Annex C1c, in order. Raises Abort on any failure."""
    p, s, k = t["p"], t["s"], t["kappa"]
    if not (p >= 0.0 and s >= 0.0):
        raise Abort(f"C1c-1 failed for {t['label']}: p = {p}, s = {s}")
    if (p == 0.0 or s == 0.0) and k != 0.0:
        raise Abort(f"C1c-2 failed for {t['label']}: p s = 0 but kappa = {k}")
    if abs(k * k - p * s) > RANK1_RTOL * max(k * k, p * s, 1e-300):
        raise Abort(f"C1c-3 failed for {t['label']}: |kappa^2 - p s| = {abs(k * k - p * s):.3e}")
    return bool(p == 0.0 or s == 0.0)


def b4_spectrum(lam, t):
    """Returns (spectrum, reason). spectrum = (c0 sorted, A = M^2 beta, x, bhat) or None.

    reason in {"ok", "complex_roots", "numerical_error", "nonfinite", "degenerate"}:
      complex_roots   : lambda_eff.A_coefficients raised ValueError with exactly its registered message
                        "local roots are complex" (|Im| > 1e-9)
      numerical_error : any other exception (any other ValueError, LinAlgError, RuntimeError, TypeError,
                        ...), raised anywhere in Params construction, A_coefficients or the array
                        conversion, or output of the wrong shape
      nonfinite       : c0 or A_i not finite (e.g. F'(c0) = 0 at a repeated root; static review R4)
      degenerate      : zero root span (x undefined) or all A_i = 0 (bhat undefined)
    """
    try:
        P_ = lam.Params(t["p"], t["s"], t["kappa"], t["Lambda"], t["w"], t["q"], t["v"], t["r"])
        c0, A = lam.A_coefficients(P_)
        c0, A = np.asarray(c0, float), np.asarray(A, float)
    except ValueError as exc:
        # only the registered message of lambda_eff.A_coefficients counts as complex roots (static review R4)
        if str(exc) == COMPLEX_ROOTS_MESSAGE:
            return None, "complex_roots"
        return None, "numerical_error"
    except Exception:                            # any other exception: recorded, never propagated silently
        return None, "numerical_error"
    if c0.shape != (4,) or A.shape != (4,):
        return None, "numerical_error"
    if not (np.all(np.isfinite(c0)) and np.all(np.isfinite(A))):
        return None, "nonfinite"
    x, bh = x_structure(c0), bhat(A)
    if x is None or bh is None:
        return None, "degenerate"
    return (c0, A, x, bh), "ok"


def select_targets(lam, repo):
    with open(os.path.join(repo, "data/B3_max_speed_7_315096.json")) as fh:
        anchor = json.load(fh)
    with open(os.path.join(repo, "data/e3_manifest.json")) as fh:
        manifest = json.load(fh)
    candidates = [target_record(anchor, f"anchor:{anchor.get('index')}")]
    for row in manifest["selected"]:
        candidates.append(target_record(row, f"manifest:{row['mode']}:{row['index']}"))
    targets, not_kept = [], {}
    for t in candidates:
        spec, reason = b4_spectrum(lam, t)
        if spec is not None:                     # four real roots (|Im| <= 1e-9 in A_coefficients)
            targets.append(t)
        else:
            not_kept[t["label"]] = reason
    n_not_complex = sum(1 for r in not_kept.values() if r != "complex_roots")
    if len(targets) != EXPECTED_N_TARGETS or n_not_complex:
        # R4: a base target lost for any reason other than complex roots is an ABORT, never a silent exclusion
        raise Abort(f"target count {len(targets)} != {EXPECTED_N_TARGETS} or non-root exclusions; "
                    f"not kept: {not_kept}")
    n_degenerate = 0
    for t in targets:                            # C1c on every target, before any search
        n_degenerate += c1c_check(t)
    for i, t in enumerate(targets):
        t["background_index"] = i                # order: anchor first, then manifest 'selected' order
    return targets, n_degenerate, not_kept


def conditioning_label(lam, t):
    P_ = lam.Params(t["p"], t["s"], t["kappa"], t["Lambda"], t["w"], t["q"], t["v"], t["r"])
    c = lam.conditioning(P_)
    ill = bool(np.min(c["rel_gap"]) < COND_REL_GAP or np.min(c["abs_Fprime"]) < COND_FPRIME)
    return {"abs_Fprime": [float(v) for v in c["abs_Fprime"]],
            "rel_gap": [float(v) for v in c["rel_gap"]],
            "near_zero_root": [bool(v) for v in c["near_zero_root"]],
            "ill_conditioned": ill}


# ----------------------------------------------------------------------------------------------
# B. Search
# ----------------------------------------------------------------------------------------------
EVAL_STATUSES = ("admissible", "out_of_region", "undefined_dispersion", "nonfinite_residual")


class Residual:
    """Residual function with hard rejection (B2) and a counter of every call (B5)."""

    def __init__(self, bec, x4, bh4, region_cls):
        self.bec, self.x4, self.bh4, self.region_cls = bec, x4, bh4, region_cls
        self.n_calls = 0
        self.n_penalty = 0
        self.first_is_penalty = None

    def evaluate(self, y):
        """(residual vector or None, status) without counting.

        status (R10, script v0.4), one of EVAL_STATUSES:
          "admissible"           : in the target region, dispersion defined, residual finite
          "out_of_region"        : the point is not in the target region (any other class)
          "undefined_dispersion" : BEC roots or beta not defined, or zero span / zero beta
          "nonfinite_residual"   : residual computed but not finite (returned for diagnosis)
        """
        P = family_to_bec(self.bec, y)
        cl = classify(self.bec, P)
        if cl["cls"] != self.region_cls:
            return None, "out_of_region"
        spec = bec_spectrum(self.bec, P)
        if spec is None:
            return None, "undefined_dispersion"
        xb, bhb = x_structure(spec[0]), bhat(spec[1])
        if xb is None or bhb is None:
            return None, "undefined_dispersion"
        r = np.concatenate(((self.x4[[1, 2]] - xb[[1, 2]]) / TH_X, eps_vec(self.bh4, bhb) / TH_EPS))
        if not np.all(np.isfinite(r)):
            return r, "nonfinite_residual"
        return r, "admissible"

    def __call__(self, y):
        self.n_calls += 1
        r, status = self.evaluate(y)
        ok = status == "admissible"
        if not ok:                                # B2: every non-admissible status gets the penalty during the search
            self.n_penalty += 1
            r = np.full(N_RES, PENALTY_VALUE)
        if self.first_is_penalty is None:
            self.first_is_penalty = not ok
        return r


def check_final(sol, evaluator, where):
    """R10: finiteness of the solver output, then the cost consistency check.

    Three distinct states (annex §I R10):
      (a) inadmissible final point (out_of_region / undefined_dispersion) -> returned, no ABORT here
      (b) non-finite cost, x or final residual                              -> ABORT
      (c) cost mismatch at an admissible final point                        -> ABORT
    Returns (final residual or None, final point admissible?, final status)."""
    if not np.isfinite(sol.cost):
        raise Abort(f"R10: non-finite solver cost at {where}")
    if not np.all(np.isfinite(sol.x)):
        raise Abort(f"R10: non-finite solver x at {where}")
    r_final, status = evaluator.evaluate(sol.x)
    if status not in EVAL_STATUSES:
        raise Abort(f"R10: unknown evaluation status {status!r} at {where}")
    if status == "nonfinite_residual":
        raise Abort(f"R10: non-finite final residual at {where}")
    if status == "admissible":
        if not np.all(np.isfinite(r_final)):          # defensive; excluded by the status above
            raise Abort(f"R10: non-finite final residual at {where}")
        if abs(0.5 * float(np.dot(r_final, r_final)) - float(sol.cost)) > 1e-9 * max(1.0, float(sol.cost)):
            raise Abort(f"R10: cost mismatch at final point at {where}")
        return r_final, True, status
    return None, False, status


def start_outcome(converged, final_ok, first_is_penalty):
    """B5: one explicit outcome per start (no_admissible_start is recorded separately)."""
    if converged:
        return "converged"
    if final_ok:
        return "not_converged_in_region"          # status <= 0, final point admissible
    if first_is_penalty:
        return "rejected_start_shifted_out"       # R9: first counted evaluation was already a penalty
    return "rejected"                             # unreachable: excluded by the B2a ABORT


def draw_y(rng):
    """One draw from the annex B1 starting distributions (fixed order of draws)."""
    return np.array([rng.uniform(-2.0, 2.0), rng.uniform(0.0, 0.99),
                     rng.uniform(-4.0, 4.0), rng.uniform(-2.0, 2.0)])


def search(bec, least_squares, x4, bh4, region, background_index):
    region_cls = 1 if region == "A" else 2
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED, REGION_CODE[region], background_index])))
    starts = []
    for s_i in range(N_STARTS):
        y0, admissible = None, False
        for _ in range(MAX_START_REDRAWS):
            y = draw_y(rng)
            if Residual(bec, x4, bh4, region_cls).evaluate(y)[1] == "admissible":
                y0, admissible = y, True
                break
        if not admissible:
            starts.append({"start": s_i, "no_admissible_start": True, "final_status": "no_admissible_start"})
            continue
        res_fun = Residual(bec, x4, bh4, region_cls)
        sol = least_squares(res_fun, y0, bounds=(Y_LO, Y_HI), **LSQ_KW)
        r_final, final_ok, final_status = check_final(sol, Residual(bec, x4, bh4, region_cls),
                                        f"background {background_index}, region {region}, start {s_i}")
        # B2a: an admissible first evaluation must not end at an inadmissible point
        if res_fun.first_is_penalty is False and not final_ok:
            raise Abort(f"B2a violated: background {background_index}, region {region}, start {s_i}")
        converged = bool(sol.status > 0 and final_ok)
        rec = {"start": s_i, "no_admissible_start": False, "y0": y0.tolist(), "y": sol.x.tolist(),
               "status": int(sol.status), "message": str(sol.message), "nfev": int(sol.nfev),
               "njev": None if sol.njev is None else int(sol.njev),
               "cost": float(sol.cost), "final_in_region": bool(final_ok), "final_status": final_status,
               "rejected": not final_ok,         # B5: final point inadmissible (only possible if first eval was a penalty)
               "n_calls": res_fun.n_calls, "n_penalty_eval": res_fun.n_penalty,
               "first_eval_penalty": bool(res_fun.first_is_penalty),
               "converged": converged,
               "outcome": start_outcome(converged, final_ok, res_fun.first_is_penalty)}
        starts.append(rec)
    conv = [s for s in starts if s.get("converged")]
    best = min(conv, key=lambda s: s["cost"]) if conv else None
    summary = {
        "n_starts": N_STARTS,
        "n_no_admissible_start": sum(1 for s in starts if s.get("no_admissible_start")),
        "n_converged": len(conv),
        "n_not_converged": sum(1 for s in starts if not s.get("no_admissible_start") and not s.get("converged")),
        "n_rejected": sum(1 for s in starts if s.get("rejected")),
        "outcomes": {o: sum(1 for s in starts if s.get("outcome") == o)
                     for o in ("converged", "not_converged_in_region", "rejected_start_shifted_out", "rejected")},
        "n_first_eval_penalty": sum(1 for s in starts if s.get("first_eval_penalty")),
        "n_starts_with_penalty_eval": sum(1 for s in starts if s.get("n_penalty_eval", 0) > 0),
        "no_admissible_fit": best is None,
    }
    return best, starts, summary


# ----------------------------------------------------------------------------------------------
# C. Robustness
# ----------------------------------------------------------------------------------------------
def perturb_target(t, rng):
    """One perturbed target (annex C1, C1a, C1b, C2). Returns (target or None if C1a-inadmissible)."""
    u = rng.uniform(-1.0, 1.0, size=len(PERT_KEYS))      # always 7 draws, fixed order
    new = dict(t)
    for key, ui in zip(PERT_KEYS, u):
        x = t[key]
        if key in ("p", "s") and x == 0.0:                # C1b: an exact zero is not perturbed
            continue
        new[key] = x + max(PERT_REL * abs(x), PERT_ABS) * ui
    if (t["p"] > 0.0 and new["p"] <= 0.0) or (t["s"] > 0.0 and new["s"] <= 0.0):
        return None                                       # C1a
    if t["kappa"] == 0.0:
        new["kappa"] = 0.0                                # C1b
    else:
        new["kappa"] = math.copysign(math.sqrt(new["p"] * new["s"]), t["kappa"])
    return new


def robustness(lam, t, xb, bhb):
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED, 2, t["background_index"]])))
    n_draws = n_rej_c1a = n_pass = 0
    n_rej = {"complex_roots": 0, "numerical_error": 0, "nonfinite": 0, "degenerate": 0}
    admissible = 0
    while admissible < N_PERT and n_draws < MAX_PERT_DRAWS:
        n_draws += 1
        tp = perturb_target(t, rng)
        if tp is None:
            n_rej_c1a += 1
            continue
        spec, reason = b4_spectrum(lam, tp)
        if spec is None:
            n_rej[reason] += 1                    # C4: counted separately by reason (static review, item 4)
            continue
        _, _, x4p, bh4p = spec
        admissible += 1
        c = criteria(x4p, bh4p, xb, bhb)
        n_pass += c["c1_sign"] and c["c2_magnitude"] and c["c3_roots"]
    out = {"n_draws": n_draws, "n_admissible": admissible, "n_rejected_C1a": n_rej_c1a,
           **{f"n_rejected_{k}": v for k, v in n_rej.items()}, "n_pass_1_3": int(n_pass)}
    if admissible < N_PERT:
        out.update(c4_robust=None, status="INCONCLUSIVE (fewer than 50 admissible perturbed targets)")
    else:
        frac = n_pass / admissible
        out.update(fraction=frac, c4_robust=bool(frac >= TH_ROBUST), status="evaluated")
    return out


# ----------------------------------------------------------------------------------------------
# Per-target evaluation
# ----------------------------------------------------------------------------------------------
def evaluate_target(lam, bec, least_squares, t, control_passed):
    spec, reason = b4_spectrum(lam, t)
    if spec is None:                              # cannot happen after select_targets; guarded anyway
        raise Abort(f"B4 target {t['label']} not evaluable: {reason}")
    c0, A, x4, bh4 = spec
    out = {"label": t["label"], "background_index": t["background_index"],
           "c0_B4": c0.tolist(), "bhat_B4": bh4.tolist(), "x_B4": x4.tolist(),
           "conditioning": conditioning_label(lam, t), "regions": {}}
    for region in ("A", "B"):
        best, starts, summary = search(bec, least_squares, x4, bh4, region, t["background_index"])
        reg = {"summary": summary, "starts": starts}
        if best is None:
            reg.update(match=False, verdict_raw="NO ADMISSIBLE FIT",
                       verdict=apply_j4("NO ADMISSIBLE FIT", control_passed[region]))
        else:
            P = family_to_bec(bec, best["y"])
            spec = bec_spectrum(bec, P)
            xb, bhb = x_structure(spec[0]), bhat(spec[1])
            crit = criteria(x4, bh4, xb, bhb)
            rob = robustness(lam, t, xb, bhb)
            verdict = final_verdict(crit, rob)
            match = verdict == "MATCH"
            # v2 sensitivity: best of the 24 one-to-one assignments, at the same best point (no re-search)
            perms = []
            for perm in itertools.permutations(range(4)):
                c = criteria(x4, bh4, xb, bhb, perm)
                perms.append({"perm": list(perm), **c, "pass_1_3": c["c1_sign"] and c["c2_magnitude"] and c["c3_roots"]})
            best_perm = min(perms, key=lambda d: (not d["pass_1_3"], d["max_eps"] / TH_EPS + d["max_dx_inner"] / TH_X))
            cl = classify(bec, P)
            reg.update(best=best, c0_BEC=spec[0].tolist(), bhat_BEC=bhb.tolist(), x_BEC=xb.tolist(),
                       criteria_primary=crit, robustness=rob, match=match,
                       verdict_raw=verdict, verdict=apply_j4(verdict, control_passed[region]),
                       sensitivity_24_permutations={"best": best_perm, "n_pass_1_3": sum(d["pass_1_3"] for d in perms),
                                                    "label": "more freedom than the primary rule; not a match"},
                       energy=cl)
        reg["role"] = "primary" if region == "A" else "control (reported separately; not part of the primary result)"
        out["regions"][region] = reg
    return out


def positive_control(bec, least_squares, region):
    """Annex J: injection-recovery on synthetic targets drawn from the family, per region."""
    target_cls = 1 if region == "A" else 2
    items = []
    for i in range(N_PC_TARGETS):
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED, 6, REGION_CODE[region], i])))
        y_star = x_t = bh_t = None
        n_draws = 0
        while n_draws < MAX_PC_DRAWS:                   # J1
            n_draws += 1
            y = draw_y(rng)
            P = family_to_bec(bec, y)
            if classify(bec, P)["cls"] != target_cls:
                continue
            spec = bec_spectrum(bec, P)
            if spec is None:
                continue
            xs, bs = x_structure(spec[0]), bhat(spec[1])
            if xs is None or bs is None:
                continue
            y_star, x_t, bh_t = y, xs, bs
            break
        if y_star is None:
            items.append({"i": i, "generated": False, "n_draws": n_draws, "recovered": False})
            continue
        # J2: same search; y_star is never passed; start seeds use background_index = 1000 + i
        best, starts, summary = search(bec, least_squares, x_t, bh_t, region, PC_INDEX_OFFSET + i)
        rec = {"i": i, "generated": True, "n_draws": n_draws, "y_star": y_star.tolist(),
               "x_target": x_t.tolist(), "bhat_target": bh_t.tolist(), "search_summary": summary}
        if best is None:
            rec.update(recovered=False, reason="no admissible fit")
        else:
            spec = bec_spectrum(bec, family_to_bec(bec, best["y"]))
            crit = criteria(x_t, bh_t, x_structure(spec[0]), bhat(spec[1]))        # J3: (1)-(3) only
            rec.update(recovered=bool(crit["c1_sign"] and crit["c2_magnitude"] and crit["c3_roots"]),
                       criteria=crit, best_cost=best["cost"], best_y=best["y"],
                       dist_y_to_y_star=float(np.linalg.norm(np.asarray(best["y"]) - y_star)))
        items.append(rec)
    n_gen = sum(1 for it in items if it["generated"])
    n_rec = sum(1 for it in items if it["recovered"])
    return {"region": region,
            "n_target_fixed": N_PC_TARGETS,
            "n_generated": n_gen,
            "n_not_generated": N_PC_TARGETS - n_gen,
            "n_recovered": n_rec,
            "n_generated_not_recovered": n_gen - n_rec,
            "recovery_rate": n_rec / N_PC_TARGETS,         # J4: over the fixed number, not over n_generated
            "threshold": f">= {PC_MIN_RECOVERED} of {N_PC_TARGETS}",
            "passed": n_rec >= PC_MIN_RECOVERED,
            "scope": "sensitivity demonstrated on attainable synthetic targets from the family interior; "
                     "not a guarantee of equal sensitivity for every real B4 target (J4, J6)",
            "items": items}


def apply_j4(verdict, control_passed):
    """Annex J4: if the positive control failed in a region, NO MATCH / NO ADMISSIBLE FIT become NOT INTERPRETABLE."""
    if control_passed:
        return verdict
    if verdict in ("NO MATCH", "NO ADMISSIBLE FIT"):
        return "NOT INTERPRETABLE (insufficient search sensitivity)"
    return verdict + " [positive control failed in this region]"


def final_verdict(crit, rob):
    """v2 match = (1) and (2) and (3) and (4). R7: (1)-(3) pass but (4) not evaluable -> INCONCLUSIVE."""
    pass_1_3 = bool(crit["c1_sign"] and crit["c2_magnitude"] and crit["c3_roots"])
    if pass_1_3 and rob["c4_robust"] is True:
        return "MATCH"
    if pass_1_3 and rob["c4_robust"] is None:
        return "INCONCLUSIVE (criterion 4: fewer than 50 admissible perturbed targets)"
    return "NO MATCH"


# ----------------------------------------------------------------------------------------------
# D. Sign-tuple sampling
# ----------------------------------------------------------------------------------------------
def sign_tuple(bhb):
    return tuple(int(0 if abs(v) < TH_SIGN else (1 if v > 0 else -1)) for v in bhb)


def sample_sign_tuples(bec, region):
    target_cls = 1 if region == "A" else 2
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED, 4, REGION_CODE[region]])))
    counts, n_in, n_draws = {}, 0, 0
    n_ctrl = 0                                   # annex A1 control count (class A only)
    cls_counts = {1: 0, 2: 0, 3: 0, 4: 0}
    while n_in < N_SAMPLE_PER_REGION and n_draws < MAX_SAMPLE_DRAWS:
        n_draws += 1
        y = draw_y(rng)
        P = family_to_bec(bec, y)
        cl = classify(bec, P)
        cls_counts[cl["cls"]] += 1
        if cl["cls"] != target_cls:
            continue
        spec = bec_spectrum(bec, P)
        bhb = None if spec is None else bhat(spec[1])
        if bhb is None:
            continue
        n_in += 1
        key = ",".join(str(v) for v in sign_tuple(bhb))
        counts[key] = counts.get(key, 0) + 1
        if target_cls == 1 and cl["lam_min"] <= 0.0 and not cl["pd0"]:
            n_ctrl += 1
    return {"region": region, "n_in_region": n_in, "n_draws": n_draws,
            "acceptance_rate": n_in / n_draws if n_draws else None,
            "reached_target": n_in >= N_SAMPLE_PER_REGION,
            "class_counts_all_draws": {str(k): v for k, v in cls_counts.items()},
            "sign_tuples": dict(sorted(counts.items())),
            "control_count_A_lam_in_[-1e-9,0]_not_pd": n_ctrl if target_cls == 1 else None,
            "control_count_scope": "all class-A points of this sampling stream; per-point verdicts are also "
                                   "stored for every best search point; NOT computed for every residual "
                                   "evaluation inside the searches (R5)",
            "status": "empirical: a tuple not seen is 'not observed in the sample', not 'impossible'"}


# ----------------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------------
def write_json(path, obj):
    with open(path, "w", newline="\n") as fh:
        json.dump(obj, fh, indent=1, sort_keys=False, default=float)
        fh.write("\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    repo, out = os.path.abspath(args.repo), os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    if os.listdir(out):
        raise Abort(f"output directory not empty: {out}")

    hashes = verify_inputs(repo)                                    # 0: before any computation
    scipy_pin = verify_scipy()
    lam = load_module("b4_lambda_eff", os.path.join(repo, "theory/lambda_eff.py"))
    bec = load_module("b4_bec_dispersion", os.path.join(repo, "theory/bec_dispersion.py"))
    import scipy
    from scipy.optimize import least_squares

    context = {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
               "scipy": scipy.__version__, "script_sha256": sha256_file(os.path.abspath(__file__)),
               "inputs_sha256": hashes, "scipy_pin": scipy_pin, "start_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    print("inputs verified:", json.dumps(hashes, indent=1))

    targets, n_degenerate, not_kept = select_targets(lam, repo)     # 2: selection + C1c
    tpath = os.path.join(out, "targets.json")
    write_json(tpath, {"n_targets": len(targets), "n_with_p_or_s_zero": n_degenerate,
                       "candidates_not_kept": not_kept, "targets": targets})
    context["targets_sha256"] = sha256_file(tpath)
    print(f"targets: {len(targets)} selected; C1c passed; p*s = 0 cases: {n_degenerate}; "
          f"targets.json sha256 {context['targets_sha256']}")
    print("indices:", [t["label"] for t in targets])

    # 2b. positive control (annex J5): before any real target; written and hashed first
    pc = {region: positive_control(bec, least_squares, region) for region in ("A", "B")}
    pcpath = os.path.join(out, "positive_control.json")
    write_json(pcpath, pc)
    context["positive_control_sha256"] = sha256_file(pcpath)
    control_passed = {region: pc[region]["passed"] for region in ("A", "B")}
    for region in ("A", "B"):
        print(f"positive control {region}: recovered {pc[region]['n_recovered']}/{N_PC_TARGETS} "
              f"(generated {pc[region]['n_generated']}), passed: {pc[region]['passed']}")
    print(f"positive_control.json sha256 {context['positive_control_sha256']}")

    results = []
    for t in targets:                                               # 3
        r = evaluate_target(lam, bec, least_squares, t, control_passed)
        results.append(r)
        print(f"[{t['background_index']:2d}] {t['label']:28s}  A: {r['regions']['A']['verdict']:17s}  "
              f"B: {r['regions']['B']['verdict']:17s}  ill-cond: {r['conditioning']['ill_conditioned']}")
    write_json(os.path.join(out, "per_target.json"), results)

    sampling = [sample_sign_tuples(bec, "A"), sample_sign_tuples(bec, "B")]   # 4
    write_json(os.path.join(out, "sign_tuples.json"), sampling)

    nA = sum(r["regions"]["A"]["match"] for r in results)
    nB = sum(r["regions"]["B"]["match"] for r in results)
    ill = [r for r in results if r["conditioning"]["ill_conditioned"]]
    summary = {"n_targets": len(results),
               "PRIMARY_region_A": "match_region_A, no_admissible_fit_A, inconclusive_c4_A",
               "CONTROL_region_B": "match_region_B, no_admissible_fit_B, inconclusive_c4_B (reported separately)",
               "match_region_A": nA, "match_region_B": nB,
               "inconclusive_c4_A": sum(r["regions"]["A"]["verdict"].startswith("INCONCLUSIVE") for r in results),
               "inconclusive_c4_B": sum(r["regions"]["B"]["verdict"].startswith("INCONCLUSIVE") for r in results),
               "no_admissible_fit_A": sum(r["regions"]["A"]["verdict_raw"] == "NO ADMISSIBLE FIT" for r in results),
               "no_admissible_fit_B": sum(r["regions"]["B"]["verdict_raw"] == "NO ADMISSIBLE FIT" for r in results),
               "not_interpretable_A": sum(r["regions"]["A"]["verdict"].startswith("NOT INTERPRETABLE") for r in results),
               "not_interpretable_B": sum(r["regions"]["B"]["verdict"].startswith("NOT INTERPRETABLE") for r in results),
               "positive_control_passed": control_passed,
               "reading_of_no_match": "NO MATCH means the fixed, registered procedure found no match; with a passed "
                                      "positive control it is interpretable within the demonstrated sensitivity. "
                                      "It is not a proof that no match exists anywhere in the family (J4).",
               "ill_conditioned_targets": [r["label"] for r in ill],
               "match_A_among_ill_conditioned": sum(r["regions"]["A"]["match"] for r in ill),
               "reading": "a match is necessary, not sufficient (v2); region B is the mandatory control"}
    write_json(os.path.join(out, "summary.json"), summary)
    context["end_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(out, "RUN_CONTEXT.txt"), "w", newline="\n") as fh:
        for k, v in context.items():
            fh.write(f"{k}: {json.dumps(v) if isinstance(v, dict) else v}\n")
    with open(os.path.join(out, "SHA256SUMS.txt"), "w", newline="\n") as fh:
        for name in sorted(os.listdir(out)):
            if name != "SHA256SUMS.txt":
                fh.write(f"{sha256_file(os.path.join(out, name))}  {name}\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Abort as exc:
        print(f"ABORT: {exc}", file=sys.stderr)
        sys.exit(2)
