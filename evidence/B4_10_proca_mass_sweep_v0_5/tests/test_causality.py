import json, importlib.util, sys
from pathlib import Path
import numpy as np
import pytest
import b410._b49  # noqa: F401
from b410 import causality as ca, experiment3 as e3
from b410.fixed_eft import load_b3

ROOT = Path(__file__).resolve().parent.parent
MAN = json.loads((ROOT/"data"/"e3_manifest.json").read_text())


def pt(mode, idx):
    return e3.make_point([r for r in MAN["selected"] if r["mode"] == mode and r["index"] == idx and r["cls"] != "negative_control"][0])


def test_homogeneous_degree_six_part_of_det_H_is_Q_inf():
    for b3 in (load_b3(), pt("rank1_minus", 9745), pt("rank1_plus", 5461)):
        ratios = ca.det_ratio(b3, [-3, -1, 0.5, 2, 5])
        assert np.allclose(ratios, -1.0, atol=1e-6)         # -M^2 with M = 1, same constant for every c


def test_principal_roots_match_x_to_infinity_limit_of_full_spectrum():
    from b410 import experiment2 as e2
    b3 = load_b3(); allr, real = ca.principal_roots(b3)
    c = np.array([b["c"] for b in e2.evaluate(b3, 1e3)["branches"]])
    assert real and np.allclose(np.sort(allr.real), c, atol=2e-6)


@pytest.mark.parametrize("mode,idx", [(None, None), ("rank1_minus", 9745), ("rank1_plus", 5461)])
def test_time_domain_front_equals_principal_characteristic_speed(mode, idx):
    b3 = load_b3() if mode is None else pt(mode, idx)
    allr, real = ca.principal_roots(b3); assert real
    td, maxre = ca.evolve_front(b3)
    vr = (td[8.0]["right"] - td[4.0]["right"])/4.0; vl = (td[8.0]["left"] - td[4.0]["left"])/4.0
    assert vr == pytest.approx(allr.real.max(), rel=3e-3) and vl == pytest.approx(allr.real.min(), rel=3e-3)
    assert maxre < 1e-6                                     # no growth: well posed


def test_complex_principal_symbol_gives_high_frequency_blowup_ill_posed():
    b3 = pt("rank1_minus", 4487)
    allr, real = ca.principal_roots(b3); assert not real
    td, maxre = ca.evolve_front(b3)
    kmax = np.pi*2**15/200.0
    assert maxre == pytest.approx(np.abs(allr.imag).max()*kmax, rel=0.1)   # growth rate ~ Im c * k
    assert td[8.0]["max_amp"] > 1e50


def test_table_41_backgrounds_all_real_ones_superluminal_controls_ill_posed(tmp_path):
    spec = importlib.util.spec_from_file_location("run_b411", ROOT/"run_b411.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    mod.main(str(tmp_path), time_domain=False)
    s = json.loads((tmp_path/"b411_summary.json").read_text())
    assert s["principal_real"] == 28 and s["principal_complex"] == ["rank1_minus_boundary_4487"]
    assert s["negative_controls_principal_complex"] == 12 and s["n_subluminal"] == 0
    assert s["front_velocity_min"] == pytest.approx(1.3686, abs=1e-3) and s["front_velocity_max"] == pytest.approx(8.2606, abs=1e-3)


def test_small_amplitude_population_still_superluminal():
    sys.path.insert(0, str(ROOT/"vendor_t173"))
    from robustness_scan import generate_samples
    fr = []
    for bg in generate_samples(300, (0.05, 1.0), (0.05, 1.0), (0.0, 1.0), 0.05, 172026):
        for sg in (-1.0, 1.0):
            b3 = e3.make_point(dict(p=bg["p"], s=bg["s"], kappa=float(sg*np.sqrt(bg["p"]*bg["s"])), Lambda=bg["Lambda"], w=float(bg["w"]), q=float(bg["q"]), v=float(bg["v"]), r=float(bg["r"])))
            allr, ok = ca.principal_roots(b3); assert ok
            fr.append(np.abs(allr.real).max())
    assert min(fr) > 1.0
