import json, sys, importlib.util
from pathlib import Path
import numpy as np
import pytest
import b410._b49  # noqa: F401
from b410 import experiment3 as e3, experiment2 as e2, audit_hamiltonian as ah
from b410.fixed_eft import load_b3, couplings

ROOT = Path(__file__).resolve().parent.parent
MAN = json.loads((ROOT/"data"/"e3_manifest.json").read_text())
SEL = MAN["selected"]


def get(mode, idx):
    return [r for r in SEL if r["mode"] == mode and r["index"] == idx][0]


def test_manifest_quotas_and_rank_one_identity():
    from collections import Counter
    c = Counter((r["mode"], r["cls"]) for r in SEL)
    for mode in ("rank1_minus", "rank1_plus"):
        assert c[(mode, "stable")] == 8 and c[(mode, "boundary")] == 6 and c[(mode, "negative_control")] == 6
    assert len(SEL) == 40
    for r in SEL:
        assert abs(r["p"]*r["s"] - r["kappa"]**2) < 1e-15
        assert (r["kappa"] < 0) == (r["mode"] == "rank1_minus")


def test_manifest_backgrounds_regenerate_from_vendored_t173_generator():
    sys.path.insert(0, str(ROOT/"vendor_t173"))
    import robustness_scan as rs
    gs = {g["index"]: g for g in rs.generate_samples(10000, (0.05, 1.0), (0.05, 1.0), (0.0, 1.0), 0.95, 172026)}
    for r in SEL[::4]:
        g = gs[r["index"]]
        assert g["p"] == r["p"] and g["w"] == pytest.approx(r["w"], rel=1e-14) and g["r"] == pytest.approx(r["r"], rel=1e-14)


def test_anchor_reproduces_experiment2():
    b3 = load_b3(); res = e3.run_background(b3, mp_every=60); s = e3.summarize_background(b3, res)
    assert s["full_real_all_x"] and s["h_positive_all_nodes"] and s["n_negative_energy_samples"] == 0
    assert s["exchanged_branches_segment_ends"] == [0, 1, 3, 5]
    assert sorted(c["x"] for c in s["crossings"])[0] == pytest.approx(0.1399, rel=2e-3)
    assert not any(any(c["in_window"].values()) for c in s["crossings"])


def test_analytic_static_criterion_predicts_hamiltonian_positivity_at_small_x():
    """-C positive definite for the rank-one family: (1+3 Lam q^2 - Xw) > 0, 1 - Xv > 0, (1+3 Lam q^2)(1-Xv) > Xw."""
    pts = [load_b3()] + [e3.make_point(r) for r in SEL]
    for b3 in pts:
        Xw, Xv, a = b3.p*b3.w**2, b3.s*b3.v**2, 1 + 3*b3.Lam*b3.q**2
        crit = (a - Xw > 0) and (1 - Xv > 0) and (a*(1 - Xv) - Xw > 0)
        g = couplings(b3, 1.0); A = (1.0, *g, b3.w, b3.q, b3.v, b3.r, b3.Lam)
        ev = np.linalg.eigvalsh(ah.hamiltonian_matrix(1e-3, *A))
        assert crit == bool(ev[0] > 1e-13*max(1, ev[-1]))


def test_hamiltonian_inertia_equals_number_of_negative_energy_branches():
    for mode, idx in (("rank1_minus", 5713), ("rank1_plus", 7900)):
        b3 = e3.make_point(get(mode, idx)); res = e3.run_background(b3, n0=61, mp_every=10**6)
        s = e3.summarize_background(b3, res)
        assert s["inertia_mismatch_nodes"] == 0 and s["n_negative_energy_samples"] > 0
        assert s["krein_sign_changes"] == 0 and s["full_real_all_x"]


def test_reality_loss_is_a_krein_collision_of_opposite_energy_branches():
    b3 = e3.make_point(get("rank1_minus", 4487)); res = e3.run_background(b3, n0=121, mp_every=10**6)
    assert len(res["transitions"]) == 1 and not res["transitions"][0]["to_real"]
    assert res["transitions"][0]["x"] == pytest.approx(0.18869, rel=1e-3)
    last = [n for n in res["nodes"] if n["real"]][-1]
    om = np.array([b["omega"] for b in last["branches"]]); E = np.array([b["E"] for b in last["branches"]])
    j = int(np.argmin(np.diff(om)))
    assert E[j]*E[j+1] < 0 and abs(E[j]) == pytest.approx(abs(E[j+1]), rel=1e-3)


def test_positive_control_all_energies_positive():
    b3 = e3.make_point(get("rank1_plus", 5461)); res = e3.run_background(b3, n0=121, mp_every=40)
    s = e3.summarize_background(b3, res)
    assert s["h_positive_all_nodes"] and s["n_negative_energy_samples"] == 0 and s["mp_reality_agrees"]


def test_negative_controls_are_complex_everywhere_with_indefinite_hamiltonian():
    for r in [q for q in SEL if q["cls"] == "negative_control"][:4]:
        b3 = e3.make_point(r); res = e3.run_background(b3, n0=61, mp_every=10)
        assert not any(n["real"] for n in res["nodes"])
        assert not any(n["h_positive"] for n in res["nodes"])
        assert all(m["mp_complex"] and m["qz_complex"] for m in res["mp"])


def test_driver_runs_and_writes_artifacts(tmp_path):
    spec = importlib.util.spec_from_file_location("run_experiment3", ROOT/"run_experiment3.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    mod.main(str(tmp_path), limit=3)
    for f in ("e3_background_summary.csv", "e3_nodes.csv", "e3_summary.json"):
        assert (tmp_path/f).stat().st_size > 100
    assert json.loads((tmp_path/"e3_summary.json").read_text())["anchor_reported_separately"] is True
