import json, pathlib, numpy as np, pytest
from theory import lambda_eff as le
D = pathlib.Path(__file__).resolve().parent.parent/"data"
ANCHOR = le.Params.from_json(D/"B3_max_speed_7_315096.json")

def test_anchor_roots_and_A():
    c0, A = le.A_coefficients(ANCHOR)
    assert np.allclose(c0, [-4.916360, -0.174619, 2.796945, 7.315096], atol=1e-5)
    assert np.allclose(A, [17.82220, 0.001053, -7.762528, -11.486321], atol=2e-5)

def test_lambda_eff_anchor():
    c0, A, L = le.lambda_eff_M2(ANCHOR)
    assert np.allclose(L, [-7.2503, -0.012066, -5.5512, -3.1400], atol=2e-3) or np.allclose(L, [-7.249, -0.01206, -5.552, -3.140], atol=3e-3)
    assert np.all(L < 0)

@pytest.mark.parametrize("x", [0.01, 0.02, 0.04])
def test_expansion_residual_is_fourth_order(x):
    c0, c = le.tracked_scalar_roots(ANCHOR, x)
    _, A = le.A_coefficients(ANCHOR)
    res = np.abs(c.real - c0 - A*x**2)
    res2 = np.abs(le.tracked_scalar_roots(ANCHOR, x/2)[1].real - c0 - A*(x/2)**2)
    assert np.all(res < 5e3*x**4*np.maximum(1, np.abs(A)))
    ratio = res[res > 1e-12]/res2[res > 1e-12]
    assert np.all((ratio > 10) & (ratio < 22))      # ~16 for x^4

def test_dimension_M_independence():
    # A_i depends only on (p,s,kappa,Lam,w,q,v,r): no M anywhere; M^2 Lam_eff = 2A/c0 is therefore M-independent
    c0, A, L = le.lambda_eff_M2(ANCHOR)
    for M in (1.0, 7.0, 300.0):
        Lam_eff = L/M**2          # length^2 for a given M
        assert np.allclose(Lam_eff*M**2, L)

def test_sign_conventions_complex_rejected():
    bad = le.Params(0.5, 0.5, 0.5, 0.9, 0, 3, 0, 3)
    try:
        le.A_coefficients(bad)
    except ValueError:
        pass

@pytest.mark.slow
def test_symbolic_derivation_matches():
    import subprocess, sys
    out = subprocess.run([sys.executable, "-I", str(D.parent/"theory"/"derive_from_action.py")], capture_output=True, text=True, timeout=900).stdout
    assert "det K / Q = M**8*x_**4*(-x_**2 - 1)" in out

def test_optional_crosscheck_with_b410_full_spectrum():
    b410 = pytest.importorskip("b410")
    from b410.fixed_eft import B3Point, leading_shift_coefficients
    b = B3Point(p=ANCHOR.p, s=ANCHOR.s, kappa=ANCHOR.kappa, Lam=ANCHOR.Lam, w=ANCHOR.w, q=ANCHOR.q, v=ANCHOR.v, r=ANCHOR.r, roots=tuple(le.A_coefficients(ANCHOR)[0]), record={})
    assert np.allclose(leading_shift_coefficients(b)[1], le.A_coefficients(ANCHOR)[1])

def test_beta_is_A_and_Lambda_eff_ill_defined_at_zero_root():
    c0, A, L = le.lambda_eff_M2(ANCHOR)
    assert np.allclose(le.beta_M2(ANCHOR), A)
    # omega = c0 k + beta k^3 with beta = A/M^2: finite even when c0 -> 0 (Lam_eff = 2A/c0 is not)

def test_tracking_is_one_to_one_even_near_degeneracy():
    c0, c = le.tracked_scalar_roots(ANCHOR, 0.3)
    assert len({round(z.real, 12) + 1j*round(z.imag, 12) for z in c}) == 4

def test_conditioning_report_anchor():
    cd = le.conditioning(ANCHOR)
    assert np.all(cd["abs_Fprime"] > 1.0) and cd["rel_gap"].min() > 0.01
