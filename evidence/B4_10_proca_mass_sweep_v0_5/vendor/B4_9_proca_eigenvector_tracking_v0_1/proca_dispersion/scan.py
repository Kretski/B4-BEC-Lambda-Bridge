"""Continuation-based spectral scans for the finite-mass Proca dispersion."""
from itertools import permutations
import numpy as np
from .model import full_roots, local_roots, pole_frequency

def _assignment(previous, current):
    """Minimum-total-distance bijection; six roots make exhaustive assignment cheap."""
    best_perm = None
    best_cost = float("inf")
    for perm in permutations(range(len(current))):
        cost = sum(abs(previous[i]-current[perm[i]]) for i in range(len(previous)))
        if cost < best_cost:
            best_cost, best_perm = cost, perm
    return np.asarray([current[j] for j in best_perm], dtype=complex)

def track_roots(bg, k_values):
    """Track all roots along ordered positive k values by global nearest assignment."""
    ks = np.asarray(k_values, dtype=float)
    if ks.ndim != 1 or len(ks) < 2 or np.any(~np.isfinite(ks)) or np.any(ks <= 0):
        raise ValueError("k_values must contain at least two finite positive values.")
    if np.any(np.diff(ks) <= 0):
        raise ValueError("k_values must be strictly increasing.")
    roots0 = np.asarray(full_roots(bg, float(ks[0])), dtype=complex)
    tracks = np.empty((len(ks), len(roots0)), dtype=complex)
    tracks[0] = roots0[np.argsort(roots0.real)]
    for n in range(1, len(ks)):
        raw = np.asarray(full_roots(bg, float(ks[n])), dtype=complex)
        tracks[n] = _assignment(tracks[n-1], raw)
    return tracks

def summarize_scan(bg, k_values, tracks=None):
    ks = np.asarray(k_values, dtype=float)
    if tracks is None:
        tracks = track_roots(bg, ks)
    imag = np.abs(tracks.imag)
    # Relative imaginary component avoids classifying high-frequency modes by absolute scale.
    scale = np.maximum(1.0, np.abs(tracks.real))
    complex_mask = imag > 1e-6*scale
    pair_gaps = []
    positive_gaps = []
    for ik, row in enumerate(tracks):
        vals = np.sort_complex(row)
        pair_gaps.extend(abs(vals[i]-vals[j]) for i in range(len(vals)) for j in range(i+1,len(vals)))
        positive = np.sort(row.real[(row.real > 0) & (np.abs(row.imag) <= 1e-6*np.maximum(1.0,np.abs(row.real)))])
        if len(positive) >= 2:
            positive_gaps.extend((abs(positive[i+1]-positive[i]), float(ks[ik])) for i in range(len(positive)-1))
    local = local_roots(bg)
    low = tracks[np.argmin(np.abs(tracks), axis=1)]
    return {
        "n_k": int(len(ks)),
        "n_modes": int(tracks.shape[1]),
        "complex_points": int(np.count_nonzero(complex_mask)),
        "max_abs_imag_omega": float(np.max(imag)),
        "max_relative_imag": float(np.max(imag/scale)),
        "all_real_within_tolerance": bool(not np.any(complex_mask)),
        "minimum_pairwise_frequency_gap": float(min(pair_gaps)),
        "minimum_positive_frequency_gap": float(min(positive_gaps, key=lambda x:x[0])[0]) if positive_gaps else None,
        "k_at_minimum_positive_gap_over_M": float(min(positive_gaps, key=lambda x:x[0])[1]/bg.M) if positive_gaps else None,
        "local_characteristic_roots": [[float(x.real),float(x.imag)] for x in local],
        "k_min_over_M": float(ks[0]/bg.M),
        "k_max_over_M": float(ks[-1]/bg.M),
    }

def scan_spectrum(bg, k_values):
    tracks = track_roots(bg, k_values)
    return tracks, summarize_scan(bg, k_values, tracks)


def avoided_crossing_diagnostic(bg, k_min_over_M, k_max_over_M, n=401):
    """Locate the narrowest gap between the two highest positive real modes.

    This is a diagnostic for a candidate scalar/Proca avoided crossing, not a
    universal mode-identity classifier. The interval must bracket the crossing.
    """
    if not (0 < k_min_over_M < k_max_over_M) or n < 3:
        raise ValueError("Require 0<k_min_over_M<k_max_over_M and n>=3.")
    ratios = np.linspace(k_min_over_M, k_max_over_M, n)
    best = None
    for ratio in ratios:
        roots = np.asarray(full_roots(bg, float(ratio*bg.M)))
        real_positive = np.sort(roots.real[(roots.real > 0) &
            (np.abs(roots.imag) <= 1e-6*np.maximum(1.0,np.abs(roots.real)))])
        if len(real_positive) < 2:
            continue
        lower, upper = real_positive[-2], real_positive[-1]
        gap = float(upper-lower)
        if best is None or gap < best["gap_omega"]:
            best = {"k_over_M":float(ratio), "omega_lower":float(lower),
                    "omega_upper":float(upper), "gap_omega":gap}
    return best
