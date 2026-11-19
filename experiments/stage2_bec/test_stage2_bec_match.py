"""Tests for stage2_bec_match.py (script v0.4.2), for static review.

They use synthetic inputs only. None of them evaluates a registered target, runs a search on a
target, or samples the family at Stage 2 scale.

Not executed: running them needs a separate authorisation. When authorised:
    B4_REPO=<fresh clean LF clone> python3 -I -m pytest -q test_stage2_bec_match.py
(-I drops the script directory from sys.path, so the module is loaded by path below.)
"""
import importlib.util
import math
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get("B4_REPO")
pytestmark = pytest.mark.skipif(not REPO, reason="set B4_REPO to a clean clone")


def _load(name, path):
    # registered in sys.modules before execution, as in the script's load_module (v0.4.2)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    try:
        spec.loader.exec_module(mod)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    return mod


@pytest.fixture(scope="module")
def S():
    return _load("stage2_bec_match", os.path.join(HERE, "stage2_bec_match.py"))


@pytest.fixture(scope="module")
def mods(S):
    S.verify_inputs(REPO)
    lam = S.load_module("t_lambda_eff", os.path.join(REPO, "theory/lambda_eff.py"))
    bec = S.load_module("t_bec_dispersion", os.path.join(REPO, "theory/bec_dispersion.py"))
    return lam, bec


# ------------------------------------------------------------------ E: inputs
def test_T01_hash_mismatch_aborts(S, tmp_path):
    for rel in S.EXPECTED_SHA256:
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"x")
    with pytest.raises(S.Abort):
        S.verify_inputs(str(tmp_path))


def test_T02_hashes_of_clean_clone_pass(S):
    assert set(S.verify_inputs(REPO)) == set(S.EXPECTED_SHA256)


# ------------------------------------------------------------------ B1: units of the family
def test_T03_family_units(S, mods):
    _, bec = mods
    y = np.array([math.log10(0.3), 0.7, 0.4, math.log10(2.5)])
    P = S.family_to_bec(bec, y)
    a1, a2, b1, b2, gamma = bec.coefficients(P)
    assert a1 == pytest.approx(1.0) and b1 == pytest.approx(0.25)
    assert a2 == pytest.approx(0.3) and b2 / b1 == pytest.approx(2.5)
    assert gamma == pytest.approx(0.7 ** 2 * a1 * a2)
    assert (P.V1 - P.V2) / math.sqrt(a1) == pytest.approx(0.4)
    assert P.V1 == pytest.approx(-P.V2)                       # A2: symmetric frame


def test_T04_scipy_pin_aborts_on_mismatch(S, monkeypatch):
    # (former T04, zero-momentum frame, removed with that sensitivity: static review R2)
    monkeypatch.setattr(S, "EXPECTED_SCIPY_VERSION", "0.0.0")
    with pytest.raises(S.Abort):
        S.verify_scipy()
    monkeypatch.undo()
    # hash mismatch, independent of the installed version (script v0.4): the version check is made to
    # pass by setting the expected version to the installed one, so only the hash check can fail
    import scipy
    monkeypatch.setattr(S, "EXPECTED_SCIPY_VERSION", scipy.__version__)
    monkeypatch.setattr(S, "EXPECTED_TRF_SHA256", "0" * 64)
    with pytest.raises(S.Abort, match="trf.py sha256"):
        S.verify_scipy()


# ------------------------------------------------------------------ A1, A3: energy form
def test_T05_energy_matrix_entries(S, mods):
    _, bec = mods
    P = S.family_to_bec(bec, np.array([math.log10(0.5), 0.3, 0.2, 0.0]))
    E = S.energy_matrices(bec, P)
    assert E.shape == (len(S.K_GRID), 4, 4)
    assert np.allclose(E, np.transpose(E, (0, 2, 1)))
    E0 = E[0]
    assert E0[0, 0] == 1.0 and E0[2, 2] == 1.0 and E0[0, 2] == pytest.approx(0.3)
    assert E0[0, 1] == pytest.approx(P.V1) and E0[2, 3] == pytest.approx(P.V2 / math.sqrt(0.5))


def test_T06_grid_minimum_at_k0(S, mods):
    _, bec = mods
    for y in ([0.0, 0.2, 0.1, 0.0], [0.5, 0.9, 1.5, -1.0], [-1.0, 0.0, 3.0, 1.0]):
        P = S.family_to_bec(bec, np.array(y))
        lam_min, lam0, _ = S.energy_lambda(bec, P)
        assert lam_min == pytest.approx(lam0, abs=1e-12)


def test_T07_one_component_landau_limit(S, mods):
    # f = 0 decouples the blocks; block 1 is PD iff V1^2 < a1 = 1
    _, bec = mods
    for dV, pd in ((1.9, True), (2.1, False)):            # V1 = dV/2
        P = S.family_to_bec(bec, np.array([2.0, 0.0, dV, 0.0]))   # rho = 100: block 2 stays PD
        _, lam0, pd0 = S.energy_lambda(bec, P)
        assert pd0 is pd


# ------------------------------------------------------------------ A5, v2: classification
def test_T08_classify_static_miscible_is_A(S, mods):
    _, bec = mods
    cl = S.classify(bec, S.family_to_bec(bec, np.array([0.0, 0.5, 0.0, 0.0])))
    assert cl["cls"] == 1 and cl["pd0"] is True


def test_T09_classify_returns_valid_class(S, mods):
    _, bec = mods
    rng = np.random.Generator(np.random.PCG64(1))
    for _ in range(50):
        cl = S.classify(bec, S.family_to_bec(bec, S.draw_y(rng)))
        assert cl["cls"] in (1, 2, 3, 4)
        if cl["cls"] == 1:
            assert cl["lam_min"] >= S.LAM_STABLE
        if cl["cls"] == 2:
            assert cl["lam_min"] < S.LAM_UNSTABLE


# ------------------------------------------------------------------ v2: comparison quantities
def test_T10_x_and_bhat(S):
    assert np.allclose(S.x_structure([3.0, -1.0, 0.0, 1.0]), [0.0, 0.25, 0.5, 1.0])
    assert np.allclose(S.bhat([2.0, -4.0, 1.0, 0.0]), [0.5, -1.0, 0.25, 0.0])
    assert S.x_structure([1.0, 1.0, 1.0, 1.0]) is None and S.bhat([0, 0, 0, 0]) is None


def test_T11_eps_bounded_by_2(S):
    rng = np.random.Generator(np.random.PCG64(2))
    for _ in range(1000):
        a, b = rng.uniform(-1, 1, 4), rng.uniform(-1, 1, 4)
        assert np.all(S.eps_vec(a, b) <= 2.0 + 1e-15)


def test_T12_criteria_identity_pass_and_fail(S):
    x = np.array([0.0, 0.3, 0.6, 1.0]); b = np.array([1.0, -0.4, 0.02, 0.5])
    c = S.criteria(x, b, x, b)
    assert c["c1_sign"] and c["c2_magnitude"] and c["c3_roots"]
    c = S.criteria(x, b, x + np.array([0, 0.06, 0, 0]), b)
    assert not c["c3_roots"]
    c = S.criteria(x, b, x, b * np.array([1, -1, 1, 1]))
    assert not c["c1_sign"]


def test_T13_B2a_cost_bound(S):
    assert S.FEASIBLE_COST_BOUND == pytest.approx(528.0)
    assert S.PENALTY_COST == pytest.approx(3e6)
    assert S.FEASIBLE_COST_BOUND < S.PENALTY_COST


# ------------------------------------------------------------------ B2, B5: residual wrapper
def test_T14_wrapper_counts_every_call_and_penalises(S, mods):
    _, bec = mods
    y_ref = np.array([0.0, 0.5, 0.0, 0.0])               # class A point used as a synthetic target
    c0, beta = S.bec_spectrum(bec, S.family_to_bec(bec, y_ref))
    R = S.Residual(bec, S.x_structure(c0), S.bhat(beta), region_cls=1)
    r = R(y_ref)
    assert R.n_calls == 1 and R.n_penalty == 0 and R.first_is_penalty is False
    assert np.allclose(r, 0.0)
    R2 = S.Residual(bec, S.x_structure(c0), S.bhat(beta), region_cls=2)   # same point, wrong region
    r2 = R2(y_ref)
    assert R2.n_penalty == 1 and np.all(r2 == S.PENALTY_VALUE) and R2.first_is_penalty is True


def test_T15_jacobian_calls_are_counted(S, mods):
    # least_squares on a synthetic target; every call (including FD Jacobian) must reach the wrapper
    from scipy.optimize import least_squares
    _, bec = mods
    y_ref = np.array([0.2, 0.6, 0.3, 0.1])
    c0, beta = S.bec_spectrum(bec, S.family_to_bec(bec, y_ref))
    R = S.Residual(bec, S.x_structure(c0), S.bhat(beta), region_cls=S.classify(bec, S.family_to_bec(bec, y_ref))["cls"])
    calls = []
    def f(y):
        calls.append(1)
        return R(y)
    sol = least_squares(f, y_ref + 0.01, bounds=(S.Y_LO, S.Y_HI), **S.LSQ_KW)
    assert R.n_calls == len(calls) > 0
    # EXPECTATION, NOT YET CONFIRMED: function evaluations plus one extra probe per variable for each
    # 2-point finite-difference Jacobian (4 search variables), i.e. nfev + 4*njev, as nfev and njev are
    # counted by scipy 1.18.1 (the pinned version). If this assertion fails when the tests are first run
    # under the pinned version, it is reported as a failed expectation and NOT adjusted to make it pass;
    # any change goes through review. The decisive B5 check is n_calls == len(calls) above.
    assert R.n_calls == sol.nfev + 4 * sol.njev


# ------------------------------------------------------------------ C1c, C1a, C1b
def _t(**kw):
    base = dict(label="t", p=0.05, s=0.3, kappa=-math.sqrt(0.015), Lambda=0.5, w=-1.0, q=-1.2, v=-0.3, r=-0.8)
    base.update(kw)
    return base


def test_T16_c1c_rules(S):
    assert S.c1c_check(_t()) is False
    with pytest.raises(S.Abort):
        S.c1c_check(_t(p=-1.0, s=-1.0, kappa=1.0))      # rank-one consistent, physically inadmissible
    with pytest.raises(S.Abort):
        S.c1c_check(_t(p=0.0, kappa=0.1))
    with pytest.raises(S.Abort):
        S.c1c_check(_t(kappa=0.2))                      # rank-one inconsistent
    assert S.c1c_check(_t(p=0.0, kappa=0.0)) is True


def test_T17_perturbation_rules(S):
    rng = np.random.Generator(np.random.PCG64(3))
    t = _t()
    for _ in range(200):
        tp = S.perturb_target(t, rng)
        if tp is None:
            continue
        assert tp["kappa"] ** 2 == pytest.approx(tp["p"] * tp["s"])
        assert math.copysign(1.0, tp["kappa"]) == math.copysign(1.0, t["kappa"])
        for k in S.PERT_KEYS:
            assert abs(tp[k] - t[k]) <= max(0.02 * abs(t[k]), 0.005) + 1e-15
    t0 = _t(p=0.0, kappa=0.0)
    for _ in range(50):
        tp = S.perturb_target(t0, rng)
        if tp is not None:
            assert tp["p"] == 0.0 and tp["kappa"] == 0.0
    tiny = _t(p=0.001, kappa=-math.sqrt(0.001 * 0.3))   # absolute 0.005 dominates: p can cross 0
    n_none = sum(S.perturb_target(tiny, rng) is None for _ in range(400))
    assert n_none > 0                                   # C1a rejections occur


def test_T18_sign_tuple(S):
    assert S.sign_tuple([0.5, -0.049, -0.05, 0.0]) == (1, 0, -1, 0)


# ------------------------------------------------------------------ R7, R10, R4/C4 (script v0.2)
def test_T19_verdict_R7(S):
    ok = {"c1_sign": True, "c2_magnitude": True, "c3_roots": True}
    bad = {"c1_sign": True, "c2_magnitude": False, "c3_roots": True}
    assert S.final_verdict(ok, {"c4_robust": True}) == "MATCH"
    assert S.final_verdict(ok, {"c4_robust": False}) == "NO MATCH"
    assert S.final_verdict(ok, {"c4_robust": None}).startswith("INCONCLUSIVE")
    assert S.final_verdict(bad, {"c4_robust": None}) == "NO MATCH"
    assert S.final_verdict(bad, {"c4_robust": True}) == "NO MATCH"


class _Ev:
    def __init__(self, r, status):
        self.r, self.status = r, status

    def evaluate(self, y):
        return self.r, self.status


def test_T20_check_final_R10(S):
    x = np.zeros(4)
    r = np.ones(6)
    assert S.check_final(SimpleNamespace(cost=3.0, x=x), _Ev(r, "admissible"), "t")[1:] == (True, "admissible")
    with pytest.raises(S.Abort, match="non-finite solver cost"):
        S.check_final(SimpleNamespace(cost=float("nan"), x=x), _Ev(r, "admissible"), "t")
    with pytest.raises(S.Abort, match="non-finite solver x"):
        S.check_final(SimpleNamespace(cost=3.0, x=np.array([0, np.nan, 0, 0])), _Ev(r, "admissible"), "t")
    with pytest.raises(S.Abort, match="non-finite final residual"):
        S.check_final(SimpleNamespace(cost=3.0, x=x), _Ev(np.array([1, 1, 1, 1, 1, np.inf]), "nonfinite_residual"), "t")
    with pytest.raises(S.Abort, match="cost mismatch"):
        S.check_final(SimpleNamespace(cost=2.0, x=x), _Ev(r, "admissible"), "t")  # 0.5*|r|^2 = 3
    # an ordinary inadmissible final point is NOT an ABORT, and its status is kept distinct
    for st in ("out_of_region", "undefined_dispersion"):
        assert S.check_final(SimpleNamespace(cost=3e6, x=x), _Ev(None, st), "t")[1:] == (False, st)
    with pytest.raises(S.Abort, match="unknown evaluation status"):
        S.check_final(SimpleNamespace(cost=3.0, x=x), _Ev(r, True), "t")


def test_T20b_real_residual_nonfinite_aborts(S, mods):
    # R10 with the real Residual.evaluate (script v0.4): a non-finite target bhat gives a non-finite
    # residual at an admissible point; evaluate reports it distinctly, the search wrapper penalises it,
    # and check_final ABORTs instead of reporting an ordinary inadmissible final point
    _, bec = mods
    y_ref = np.array([0.0, 0.5, 0.0, 0.0])               # class A (T08)
    c0, beta = S.bec_spectrum(bec, S.family_to_bec(bec, y_ref))
    bh_bad = S.bhat(beta).copy()
    bh_bad[0] = np.nan
    R = S.Residual(bec, S.x_structure(c0), bh_bad, region_cls=1)
    r, status = R.evaluate(y_ref)
    assert status == "nonfinite_residual" and not np.all(np.isfinite(r))
    assert np.all(R(y_ref) == S.PENALTY_VALUE) and R.n_penalty == 1
    with pytest.raises(S.Abort, match="non-finite final residual"):
        S.check_final(SimpleNamespace(cost=1.0, x=y_ref), S.Residual(bec, S.x_structure(c0), bh_bad, 1), "t")
    # the same point with a finite target is admissible, and with the wrong region it is out_of_region
    R_ok = S.Residual(bec, S.x_structure(c0), S.bhat(beta), region_cls=1)
    assert R_ok.evaluate(y_ref)[1] == "admissible"
    assert S.Residual(bec, S.x_structure(c0), S.bhat(beta), region_cls=2).evaluate(y_ref)[1] == "out_of_region"


def test_T21_b4_rejection_reasons(S, mods):
    lam, _ = mods
    # complex roots from the registered kernel: p = s = kappa = 0, Lambda = -1, q = 1 gives
    # F = (c^2 + 2)(c^2 - 1), so two roots are +-i sqrt(2)
    spec, reason = S.b4_spectrum(lam, _t(p=0.0, s=0.0, kappa=0.0, Lambda=-1.0, w=0.0, q=1.0, v=0.0, r=0.0))
    assert spec is None and reason == "complex_roots"

    # the other reasons are tested on the mapping itself with a stand-in kernel (a repeated root would
    # be split by ~sqrt(eps) in polyroots and could be reported as complex instead)
    def fake(c0, A=None, exc=None):
        def A_coefficients(P_):
            if exc is not None:
                raise exc
            return np.asarray(c0, float), np.asarray(A, float)
        return SimpleNamespace(Params=lambda *a: a, A_coefficients=A_coefficients)
    c = [-1.0, -0.5, 0.5, 1.0]
    assert S.b4_spectrum(fake(c, [1.0, np.inf, 1.0, 1.0]), _t())[1] == "nonfinite"
    assert S.b4_spectrum(fake(c, exc=np.linalg.LinAlgError("x")), _t())[1] == "numerical_error"
    # R4 (script v0.3): only the registered message counts as complex roots; everything else is numerical_error
    assert S.b4_spectrum(fake(c, exc=ValueError("local roots are complex")), _t())[1] == "complex_roots"
    assert S.b4_spectrum(fake(c, exc=ValueError("some other problem")), _t())[1] == "numerical_error"
    assert S.b4_spectrum(fake(c, exc=RuntimeError("x")), _t())[1] == "numerical_error"
    assert S.b4_spectrum(fake(c, exc=TypeError("x")), _t())[1] == "numerical_error"
    assert S.b4_spectrum(fake(c, [1.0, 2.0, 3.0]), _t())[1] == "numerical_error"          # wrong shape
    assert S.b4_spectrum(fake(c, [0.0, 0.0, 0.0, 0.0]), _t())[1] == "degenerate"
    assert S.b4_spectrum(fake([1.0, 1.0, 1.0, 1.0], [1.0, 2.0, 3.0, 4.0]), _t())[1] == "degenerate"
    spec, reason = S.b4_spectrum(fake(c, [1.0, -2.0, 0.5, 4.0]), _t())
    assert reason == "ok" and np.allclose(spec[3], [0.25, -0.5, 0.125, 1.0])


def test_T22_start_outcome(S):
    assert S.start_outcome(True, True, False) == "converged"
    assert S.start_outcome(False, True, False) == "not_converged_in_region"
    assert S.start_outcome(False, False, True) == "rejected_start_shifted_out"


# ------------------------------------------------------------------ annex J (script v0.3)
def test_T23_apply_j4(S):
    nip = "NOT INTERPRETABLE (insufficient search sensitivity)"
    for v in ("MATCH", "NO MATCH", "NO ADMISSIBLE FIT", "INCONCLUSIVE (criterion 4: x)"):
        assert S.apply_j4(v, True) == v
    assert S.apply_j4("NO MATCH", False) == nip
    assert S.apply_j4("NO ADMISSIBLE FIT", False) == nip
    assert S.apply_j4("MATCH", False).startswith("MATCH") and "control failed" in S.apply_j4("MATCH", False)
    assert S.apply_j4("INCONCLUSIVE (criterion 4: x)", False).startswith("INCONCLUSIVE")


@pytest.mark.parametrize("region", ["A", "B"])
def test_T24_positive_control_bookkeeping(S, mods, monkeypatch, region):
    # reduced sizes for a test only; the registered constants are not changed by the script.
    # Both regions are covered (script v0.4).
    from scipy.optimize import least_squares
    _, bec = mods
    monkeypatch.setattr(S, "N_PC_TARGETS", 2)
    monkeypatch.setattr(S, "N_STARTS", 2)
    monkeypatch.setattr(S, "PC_MIN_RECOVERED", 2)
    out = S.positive_control(bec, least_squares, region)
    assert out["region"] == region
    assert out["n_target_fixed"] == 2
    assert out["n_generated"] + out["n_not_generated"] == 2
    assert out["n_recovered"] + out["n_generated_not_recovered"] == out["n_generated"]
    assert out["recovery_rate"] == out["n_recovered"] / 2
    assert out["passed"] == (out["n_recovered"] >= 2)
    for it in out["items"]:
        if it["generated"]:
            assert "y_star" in it and it["search_summary"]["n_starts"] == 2


def test_T25_registered_constants(S):
    assert (S.N_PC_TARGETS, S.PC_MIN_RECOVERED, S.MAX_PC_DRAWS, S.PC_INDEX_OFFSET) == (29, 27, 10_000, 1000)
    assert (S.SEED, S.N_STARTS, S.N_SAMPLE_PER_REGION) == (20261111, 80, 200_000)
    assert S.COMPLEX_ROOTS_MESSAGE == "local roots are complex"


def test_T26_positive_control_generated_not_recovered(S, mods, monkeypatch):
    # forces the case "generated but not recovered" deterministically: the search is replaced by a stub
    # that never finds an admissible fit; the bookkeeping must count every generated target as not
    # recovered and must not count it as "not generated" (annex J1, script v0.4)
    _, bec = mods
    monkeypatch.setattr(S, "N_PC_TARGETS", 3)
    calls = []

    def stub_search(bec_, ls_, x_t, bh_t, region, background_index):
        calls.append(background_index)
        return None, [], {"n_starts": 0, "stub": True}
    monkeypatch.setattr(S, "search", stub_search)
    out = S.positive_control(bec, None, "A")
    assert out["n_generated"] == len(calls)
    assert out["n_recovered"] == 0
    assert out["n_generated_not_recovered"] == out["n_generated"]
    assert out["n_generated"] + out["n_not_generated"] == 3
    assert out["passed"] is False
    assert calls == [S.PC_INDEX_OFFSET + it["i"] for it in out["items"] if it["generated"]]   # J2 seeds
    for it in out["items"]:
        if it["generated"]:
            assert it["recovered"] is False and it["reason"] == "no admissible fit"
            assert "y_star" in it                                                     # stored, J1


# ------------------------------------------------------------------ v0.4.2: module loading
def test_T27_load_module_registers_and_cleans_up(S, tmp_path):
    # a registered module with @dataclass loads, is registered, and its dataclass works
    name = "t27_bec_dispersion"
    mod = S.load_module(name, os.path.join(REPO, "theory/bec_dispersion.py"))
    assert sys.modules[name] is mod
    P = mod.BECParams(m1=1.0, m2=1.0, n1=1.0, n2=1.0, g11=1.0, g22=1.0, g12=0.5, V1=0.0, V2=0.0, hbar=1.0)
    assert P.g12 == 0.5
    # the same name cannot be loaded twice
    with pytest.raises(S.Abort):
        S.load_module(name, os.path.join(REPO, "theory/bec_dispersion.py"))
    sys.modules.pop(name, None)
    # a module that fails during execution is removed from sys.modules
    bad = tmp_path / "bad_module.py"
    bad.write_text("raise RuntimeError('boom')\n")
    with pytest.raises(RuntimeError):
        S.load_module("t27_bad", str(bad))
    assert "t27_bad" not in sys.modules
