"""Stage 1: hardening of the mean-field BdG dispersion core (no LHY).  Fixed parameters, sorted k, scaled error R(k)."""
import numpy as np
import pytest
from dataclasses import replace

from theory.bec_dispersion import (BECParams, beta_coefficients, coefficients, conditioning, exact_phase_roots,
                                   galilean_shift, hydrodynamic_roots)

BASE = dict(m1=1.0, m2=1.3, n1=1.0, n2=0.8, g11=1.0, g22=1.2, g12=0.3, V1=0.15, V2=-0.1, hbar=1.0)
SETS = {
    "base": BASE,
    "weak_coupling": {**BASE, "g12": 0.01},
    "strong_stable": {**BASE, "g12": 1.0},                      # g12^2 < g11 g22 = 1.2 (miscible)
    "masses_densities": {**BASE, "m1": 0.7, "m2": 2.0, "n1": 2.5, "n2": 0.3},
    "counterflow_subsonic": {**BASE, "V1": 0.4, "V2": -0.4},
    "supersonic_counterflow": {**BASE, "V1": 3.0, "V2": -3.0},  # real roots, Landau (not dynamical) regime
}
KS = (0.01, 0.02, 0.04)                                          # sorted, fixed
FLOOR = 1e-12                                                    # truncation error must sit far above polyroots noise (~1e-15)


def mk(name):
    return BECParams(**SETS[name])


def scaled_errors(p, ks=KS):
    c0, beta = beta_coefficients(p)
    err = np.array([np.max(np.abs(exact_phase_roots(k, p) - (c0 + beta*k**2))) for k in ks])
    return err, err/np.array(ks)**4


@pytest.fixture
def params():
    return mk("base")


def test_four_real_signed_roots(params):
    roots = hydrodynamic_roots(params)
    assert len(roots) == 4 and np.all(np.isfinite(roots)) and np.all(np.diff(roots) > 0)


@pytest.mark.parametrize("name", list(SETS))
def test_scaled_error_is_constant(name):
    """R(k) = |c_exact - c0 - beta k^2| / k^4 stays ~constant for small k (leading omitted term is O(k^4))."""
    assert list(KS) == sorted(KS)
    err, R = scaled_errors(mk(name))
    assert err[0] > FLOOR, "error at the smallest k is at the numerical noise floor: increase k"
    assert np.allclose(R, R[0], rtol=2e-2)
    assert err[0] < err[1] < err[2]


@pytest.mark.parametrize("name", list(SETS))
def test_conditioning_is_reasonable(name):
    p = mk(name)
    cd = conditioning(p)
    assert np.all(cd["abs_Hprime"] > 1e-2)
    assert np.all(cd["rel_gap"] > 1e-2)


@pytest.mark.parametrize("name", list(SETS))
def test_beta_is_galilean_invariant_and_c0_shifts(name):
    p = mk(name)
    c0, beta = beta_coefficients(p)
    c1, b1 = beta_coefficients(galilean_shift(p, 3.7))
    assert np.allclose(c1, c0 + 3.7, atol=1e-10)
    assert np.allclose(b1, beta, rtol=1e-9, atol=1e-10)


def test_zero_intercomponent_coupling_gives_pairs_around_V(params):
    p = replace(params, g12=0.0)
    roots = hydrodynamic_roots(p)
    a1, a2, *_ = coefficients(p)
    expected = np.sort([p.V1 - np.sqrt(a1), p.V1 + np.sqrt(a1), p.V2 - np.sqrt(a2), p.V2 + np.sqrt(a2)])
    assert np.allclose(roots, expected, atol=1e-12)


def test_rank_one_boundary_is_flagged_not_silent():
    """g12^2 = g11 g22 (miscibility boundary): a zero-sound-speed mode; must raise, not return a number."""
    p = BECParams(**{**BASE, "g12": float(np.sqrt(BASE["g11"]*BASE["g22"])), "V1": 0.0, "V2": 0.0})
    with pytest.raises(ValueError):
        hydrodynamic_roots(p)


def test_immiscible_region_is_flagged():
    p = BECParams(**{**BASE, "g12": 1.4})
    with pytest.raises(ValueError):
        hydrodynamic_roots(p)


def test_invalid_mass_rejected(params):
    with pytest.raises(ValueError):
        hydrodynamic_roots(replace(params, m1=0.0))


@pytest.mark.parametrize("name", ["base", "strong_stable"])
def test_mpmath_high_precision_check(name):
    mp = pytest.importorskip("mpmath")
    mp.mp.dps = 50
    p = mk(name)
    c0, beta = beta_coefficients(p)
    a1, a2, b1, b2, g = [mp.mpf(x) for x in coefficients(p)]
    V1, V2 = mp.mpf(p.V1), mp.mpf(p.V2)
    for k in (0.01, 0.02):
        z = mp.mpf(k)**2
        f = lambda c: ((c - V1)**2 - a1 - b1*z)*((c - V2)**2 - a2 - b2*z) - g
        exact = np.array([float(mp.findroot(f, mp.mpf(r))) for r in c0])
        assert np.max(np.abs(exact_phase_roots(k, p) - exact)) < 1e-12     # float roots reproduce the 50-digit ones
        trunc = np.max(np.abs(exact - (c0 + beta*k**2)))
        assert trunc/k**4 > 0 and abs(trunc/k**4 - scaled_errors(p, (k,))[1][0]) < 1e-3*trunc/k**4
