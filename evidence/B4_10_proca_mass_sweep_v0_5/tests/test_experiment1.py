import numpy as np
import pytest
from b410 import experiment1 as e1
from b410 import reference as ref
from b410._b49 import Background, coupled_matrix, full_roots, local_roots, track_eigenmodes
from b410.fixed_eft import (load_b3, background, couplings, eft_coefficients,
                            leading_shift_coefficients, domain_report)

B3 = load_b3()
MASSES = [1e-3, 0.1, 1.0, 7.0, 1e3, 1e6]


def test_b3_record_is_internally_consistent():
    assert abs(B3.p*B3.s - B3.kappa**2) < 1e-15
    t17 = ref.t17_roots(*B3.args())
    assert np.abs(t17.imag).max() < 1e-12
    assert np.allclose(np.sort(t17.real), B3.roots, atol=1e-12)


def test_domain_invariants_are_recomputed_not_assumed():
    d = domain_report(B3)
    assert d["Xq_p_q2"] == pytest.approx(B3.record["Xq"], rel=1e-12)
    assert d["Xr_s_r2"] == pytest.approx(B3.record["Xr"], rel=1e-12)
    assert d["Xw_p_w2"] == pytest.approx(B3.record["Xw"], rel=1e-12)
    assert d["Xv_s_v2"] == pytest.approx(B3.record["Xv"], rel=1e-12)
    assert min(d["kinetic_matrix_A_eigenvalues"]) > 0


def test_couplings_reconstruct_fixed_eft_coefficients_at_every_mass():
    for row in e1.coupling_reconstruction(B3, MASSES):
        assert row["rel_err_p"] < 1e-13 and row["rel_err_s"] < 1e-13 and row["rel_err_kappa"] < 1e-13
    # and g scales linearly with M (fixed-EFT prescription)
    g1 = np.array(couplings(B3, 1.0)); g5 = np.array(couplings(B3, 5.0))
    assert np.allclose(g5, 5*g1, rtol=1e-14)


def test_local_limit_equals_b3_roots_for_every_mass():
    for row in e1.local_limit_table(B3, MASSES):
        assert row["max_abs_diff_vs_B3"] < 1e-12
        assert row["max_abs_diff_vs_independent_T17"] < 1e-12


def test_exact_collapse_of_the_spectrum_in_x_equals_k_over_M():
    for row in e1.scaling_collapse(B3, [0.01, 0.14, 0.5, 2.0], [0.3, 1.0, 7.0]):
        assert row["max_abs_dc_over_masses"] < 1e-9


def test_negative_control_fixed_coupling_family_does_not_collapse():
    # B4.9 family: g held fixed at their M=1 values while M varies.
    g = couplings(B3, 1.0)
    x = 0.3
    cs = []
    for M in (1.0, 3.0):
        bg = Background(M=M, g_phi=g[0], g_psi=g[1], w=B3.w, q=B3.q, v=B3.v, r=B3.r, Lambda=B3.Lam)
        k = x*M
        cs.append(np.array(sorted(full_roots(bg, k).real/k)))
    assert np.abs(cs[0] - cs[1]).max() > 1e-2


def test_leading_shift_prediction_has_fourth_order_residual():
    xs = [0.01, 0.02]
    rows = e1.shift_convergence(B3, xs)
    for branch in (0, 2, 3):
        r = {row["x"]: abs(row["residual"]) for row in rows if row["branch"] == branch}
        ratio = r[0.02]/r[0.01]
        assert 13.0 < ratio < 19.5          # x^4 scaling gives 16
    # the leading coefficients themselves
    c0, A = leading_shift_coefficients(B3)
    assert np.allclose(c0, B3.roots)
    assert A[0] > 0 and A[2] < 0 and A[3] < 0   # superluminal |c| decreases with x


def test_dressed_proca_frequency_at_k_to_zero():
    # omega_P(k->0) = M sqrt(1 + p q^2 + s r^2) = M sqrt(det A)
    pred = np.sqrt(1 + B3.p*B3.q**2 + B3.s*B3.r**2)
    for M in (1.0, 7.0):
        k = 1e-6*M
        om = np.sort(full_roots(background(B3, M), k).real)
        assert abs(om[-1]/M - pred) < 1e-5 and abs(om[0]/M + pred) < 1e-5


def test_readme_matrix_matches_b49_coupled_matrix():
    rng = np.random.default_rng(3)
    for _ in range(40):
        M = rng.uniform(0.2, 5); k = rng.uniform(0.01, 3)
        om = complex(rng.normal(), rng.normal())*2
        bg = background(B3, M)
        blocks = ref.H_blocks(M, k, bg.g_phi, bg.g_psi, bg.w, bg.q, bg.v, bg.r, bg.Lambda)
        assert np.abs(ref.H_of(om, blocks) - coupled_matrix(bg, om, k)).max() < 1e-10*max(1, abs(om)**2)


POINTS = [(1.0, 0.01), (1.0, 0.14), (1.0, 0.2077), (1.0, 0.5), (3.0, 1.0), (0.2, 0.05)]


def test_polynomial_roots_agree_with_qz_and_mpmath():
    tab = e1.verification_table(B3, POINTS)
    assert len(tab) == 6*len(POINTS)
    for row in tab:
        scale = max(1.0, abs(row["omega_poly_re"]))
        assert row["abs_diff_qz"] < 1e-9*scale
        assert row["abs_diff_mpmath"] < 1e-9*scale


def test_roots_are_simple_and_residuals_tiny():
    for row in e1.verification_table(B3, POINTS):
        assert row["residual"] < 1e-12
        assert row["simple_root_ratio"] > 1e6          # sigma_2nd_smallest >> sigma_min
        assert abs(row["omega_poly_im"]) < 1e-12*max(1.0, abs(row["omega_poly_re"]))
        assert np.isfinite(row["root_condition"])


def test_qz_has_exactly_six_finite_roots():
    bg = background(B3, 1.0)
    blocks = ref.H_blocks(1.0, 0.3, bg.g_phi, bg.g_psi, bg.w, bg.q, bg.v, bg.r, bg.Lambda)
    assert len(ref.qz_roots(blocks)) == 6


def test_exact_no_complex_frequencies_for_any_x():
    from b410.exact import discriminant_report
    rep = discriminant_report(B3)
    assert rep["discriminant_real_roots_in_(0,1]"] == []
    assert rep["n_real_roots_of_Q_at_eps_small"] == 6
    assert rep["conclusion_all_six_real_and_simple_for_all_x"] is True


def test_overlap_tracker_is_stable_along_x_sweep():
    res = track_eigenmodes(background(B3, 1.0), np.linspace(0.02, 1.5, 200))
    assert res["mean_overlap"][1:].min() > 0.999
    assert res["null_residual"].max() < 1e-12
    assert np.abs(res["omega"].imag).max() < 1e-9


def test_universal_spectrum_matches_b49_where_b49_is_valid():
    from b410.spectrum import universal_c
    for M, x in [(0.3, 0.14), (1.0, 0.01), (7.0, 0.5), (1.0, 2.0)]:
        k = x*M
        b49 = np.sort_complex(full_roots(background(B3, M), k)/k)
        uni = np.sort_complex(universal_c(B3, x))
        assert np.abs(b49 - uni).max() < 1e-8*max(1.0, np.abs(uni).max())


def test_b49_silently_drops_roots_at_large_scales_and_b410_detects_it():
    from b410.spectrum import b49_full_roots_checked, RootCountError, universal_c
    bg = background(B3, 1000.0)
    k = 140.0
    assert len(full_roots(bg, k)) < 6                  # B4.9 behaviour, documented
    with pytest.raises(RootCountError):
        b49_full_roots_checked(bg, k)
    # the universal solver still returns all six, in agreement with mpmath
    x = k/1000.0
    gphi, gpsi = couplings(B3, 1000.0)
    mpz = np.sort_complex(np.array([complex(z)/k for z in
                                    ref.mp_roots(1000.0, k, gphi, gpsi, B3.w, B3.q, B3.v, B3.r, B3.Lam)]))
    uni = np.sort_complex(universal_c(B3, x))
    assert len(uni) == 6
    assert np.abs(uni - mpz).max() < 1e-9*np.abs(mpz).max()


def test_universal_solver_stable_for_extreme_x():
    from b410.spectrum import universal_c
    for M, k in [(1e5, 1e-2), (1.0, 1e-4), (1.0, 1e3)]:
        gphi, gpsi = couplings(B3, M)
        mpz = np.sort_complex(np.array([complex(z)/k for z in
                                        ref.mp_roots(M, k, gphi, gpsi, B3.w, B3.q, B3.v, B3.r, B3.Lam)]))
        uni = np.sort_complex(universal_c(B3, k/M))
        assert np.abs(uni - mpz).max() < 1e-8*np.abs(mpz).max()
