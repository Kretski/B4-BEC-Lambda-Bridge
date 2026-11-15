import json
import numpy as np
import pytest
import b410._b49  # noqa: F401
from b410 import experiment2 as e2, reference as ref, energy as en
from b410.fixed_eft import load_b3, B3Point, leading_shift_coefficients, contact_coefficients

B3 = load_b3()


@pytest.fixture(scope="module")
def run():
    nodes, unres = e2.sweep(B3, n0=121)
    lab, st = e2.track(nodes)
    tbl = e2.tracked_table(nodes, lab)
    xs = np.array([n["x"] for n in nodes])
    return nodes, unres, lab, st, tbl, xs


def test_no_complex_roots_qz_and_mpmath(run):
    nodes, *_ = run
    assert max(n["max_imag"] for n in nodes) < 1e-12
    for x in (1e-3, 0.03, 0.2, 1.0, 30.0):
        m = e2.mp_check(B3, x)
        assert m["max_diff"] < 1e-9 and m["mp_max_abs_imag"] < 1e-12


def test_negative_control_complex_roots_are_detected():
    bad = B3Point(p=0.3394370989189766, s=0.5576612789314608, kappa=-0.43507577121685514, Lam=0.0244906774933632,
                  w=1.1617196379863794, q=1.642277396198334, v=1.2345796836260892, r=-0.5063847645833194, roots=(), record={})
    assert max(e2.evaluate(bad, x)["max_imag"] for x in np.geomspace(1e-3, 1e3, 25)) > 1.0


def test_six_branches_tracked_and_match_omega_order(run):
    nodes, unres, lab, st, tbl, xs = run
    assert not unres
    assert (lab == lab[0]).all()                      # overlap tracking == omega ordering
    assert min(s[0] for s in st) >= 0.98
    assert all(len(n["branches"]) == 6 for n in nodes)


def test_tracker_follows_a_diabatic_swap_that_omega_order_would_miss():
    """Synthetic: two Hermitian levels E=+-t with fixed diabatic vectors cross exactly at t=0."""
    def vecs(t):
        w = np.array([-abs(t), abs(t)]); V = np.eye(2)
        order = np.argsort([-t, t])                    # sort by eigenvalue
        return V[:, order]
    Pa, Pb = vecs(-0.3), vecs(0.3)
    perm, asg, mar = e2.match(Pa, Pb)
    assert list(perm) == [1, 0] and asg.min() == pytest.approx(1.0)


def test_subspace_overlap_tracks_exact_degenerate_cluster():
    rng = np.random.default_rng(0)
    Q, _ = np.linalg.qr(rng.normal(size=(5, 5)))
    def cluster(seed):
        R, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(2, 2)))
        return Q[:, :2] @ R                            # arbitrary basis of the same 2-dim eigenspace
    A, B = cluster(1), cluster(2)
    assert e2.subspace_overlap(A, B) == pytest.approx(1.0, abs=1e-12)
    perm, asg, mar = e2.match(A, B)
    assert asg.min() < 0.999 or mar.min() < 0.5       # individual vectors are NOT reliable
    other = Q[:, 2:4]
    assert e2.subspace_overlap(A, other) < 1e-12


def test_all_branches_have_positive_energy_and_krein_equals_sign_omega(run):
    *_, tbl, xs = run
    for rows in tbl:
        for r in rows:
            assert r["E"] > 0
            assert np.sign(r["K"]) == np.sign(r["omega"])
            assert r["E"] == pytest.approx(r["E_noether"], rel=1e-7)


def test_krein_theory_consistent_with_simple_real_roots(run):
    """Simple real roots => z^dag H' z != 0 (no Krein collision) and no sign change on a branch."""
    *_, tbl, xs = run
    for br in range(6):
        K = np.array([r[br]["K"] for r in tbl])
        assert (K != 0).all() and len(set(np.sign(K))) == 1
    assert min(abs(r["K"]) for rows in tbl for r in rows) > 1e-4


def test_conditioning_diagnostics_are_sane(run):
    *_, tbl, xs = run
    for rows in tbl:
        for r in rows:
            assert r["sigma_2nd"] > 1e8*max(r["sigma_min"], 1e-300) or r["sigma_min"] < 1e-6
            assert r["residual"] < 1e-12
            assert r["sigma_min"] < r["sigma_2nd"]


def test_small_x_limit_recovers_b3_and_leading_shift(run):
    nodes, *_ = run
    x = 1e-3
    nd = e2.evaluate(B3, x)
    c = np.array([b["c"] for b in nd["branches"]])[1:5]
    c0, A = leading_shift_coefficients(B3)
    assert np.allclose(c, c0 + A*x**2, atol=5e-9)
    # dressed Proca roots: omega = +-M sqrt(det A)
    Amat = np.array([[1 + B3.p*B3.q**2, -B3.kappa*B3.q*B3.r], [-B3.kappa*B3.q*B3.r, 1 + B3.s*B3.r**2]])
    omP = np.array([nd["branches"][0]["omega"], nd["branches"][5]["omega"]])
    # the pair is split at O(x) by the contact term (+-M sqrt(det A) - k*const), mean |omega| is O(x^2)-exact
    assert np.abs(omP).mean() == pytest.approx(np.sqrt(np.linalg.det(Amat)), rel=1e-5)


def test_large_x_limit_matches_exact_c_infinity():
    c_inf, im = e2.c_infinity(B3)
    assert im < 1e-12
    c = np.array([b["c"] for b in e2.evaluate(B3, 1e3)["branches"]])
    assert np.allclose(c, c_inf, atol=1e-5)


def test_composition_exchange_is_between_positive_energy_branches(run):
    nodes, unres, lab, st, tbl, xs = run
    ex = e2.exchange_summary(tbl, xs)
    exch = [e["branch"] for e in ex if e["exchanged_energy_metric"]]
    assert exch == [0, 1, 3, 5]
    for e in ex:
        assert e["krein_sign_changes"] == 0 and e["energy_sign_changes"] == 0
    cross = sorted(c for e in ex for c in e["crossing_x_energy_metric"])
    assert 0.1 < cross[0] and cross[-1] < 0.3


def test_shares_are_energy_partition(run):
    *_, tbl, xs = run
    for r in tbl[len(tbl)//2]:
        s = e2.shares(r)
        assert s["f_scalar"] + s["f_proca"] + s["f_int"] == pytest.approx(1.0, abs=1e-9)


def test_euclidean_amplitude_metric_is_not_a_justified_substitute(run):
    """Documented metric dependence: amplitude metric crosses 0.5 where energy metric does not."""
    *_, tbl, xs = run
    for br in (2, 4):
        fa = np.array([e2.shares(r[br])["amp_proca"] for r in tbl])
        fe = np.array([e2.shares(r[br])["f_proca"] for r in tbl])
        assert len(e2.crossings(xs, fa)) >= 1 and len(e2.crossings(xs, fe)) == 0


def test_eft_domain_recovery_and_marking(run):
    nodes, unres, lab, st, tbl, xs = run
    ok, X = e2.eft_domain_flags(B3)
    assert ok and X["Xw"] == pytest.approx(0.669226, abs=1e-6)
    assert not e2.eft_domain_flags(B3, e2.STRICT_CUTS)[0]
    mark = e2.eft_mark(B3, tbl)
    assert all(m["status"] == "unestablished" for m in mark)
    ex = e2.exchange_summary(tbl, xs)
    for e in ex:
        for cx in e["crossing_x_energy_metric"]:
            j = int(np.argmin(np.abs(np.log(xs/cx))))
            assert not mark[j]["per_branch_in_window_eta0.3"][e["branch"]]
            assert not mark[j]["per_branch_in_window_eta0.1"][e["branch"]]


def test_driver_writes_all_artifacts(tmp_path):
    import importlib.util, pathlib
    spec = importlib.util.spec_from_file_location("run_experiment2", pathlib.Path(__file__).resolve().parent.parent / "run_experiment2.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    mod.main(str(tmp_path), n0=41)
    for f in ("exp2_energy_tracking.csv", "exp2_krein_signature.csv", "exp2_canonical_energy.csv",
              "exp2_mixing_diagnostics.csv", "exp2_degeneracy_report.json", "exp2_summary.json"):
        assert (tmp_path/f).stat().st_size > 100
    s = json.loads((tmp_path/"exp2_summary.json").read_text())
    assert s["complex_roots_found"] is False and s["canonical_energy"]["all_positive"] is True
    assert s["eft"]["verdict"].startswith("every point")
