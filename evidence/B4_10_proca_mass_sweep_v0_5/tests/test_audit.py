"""Independent audit of Experiment 2: second variation of the ORIGINAL action (sympy),
A0 constraint, Legendre transform -> 6x6 Hamiltonian, Hamilton-equation frequencies, positivity."""
import numpy as np
import pytest
import b410._b49  # noqa: F401
from b410 import audit_hamiltonian as ah, experiment2 as e2, energy as en
from b410.fixed_eft import load_b3, B3Point, couplings

B3 = load_b3()
G = couplings(B3, 1.0)
A = (1.0, *G, B3.w, B3.q, B3.v, B3.r, B3.Lam)
NAME = {(0, None): 'dphi', (1, None): 'dpsi', (2, None): 'da0', (3, None): 'da1'}


def key(i, d):
    base = ['dphi', 'dpsi', 'da0', 'da1'][i]
    return base if d is None else f"{base}_{d}"


def test_sympy_second_variation_equals_reconstructed_L2():
    S, L1, L2 = ah.second_variation_coefficients(*A)
    d = ah.terms_from_L2(S, L2)
    mine = {}
    for coef, i, di, j, dj, g in en.terms(*A):
        kk = tuple(sorted((key(i, di), key(j, dj))))
        mine[kk] = mine.get(kk, 0.0) + coef
    assert set(d) == set(mine)
    for kk in d:
        assert d[kk] == pytest.approx(mine[kk], rel=1e-12, abs=1e-12)


def test_background_is_stationary_no_linear_term_in_undifferentiated_fields():
    S, L1, L2 = ah.second_variation_coefficients(*A)
    import sympy as sp
    assert sp.simplify(L1.coeff(S['da0'])) == 0 and abs(float(L1.coeff(S['da1']))) < 1e-25


@pytest.mark.parametrize("x", [1e-3, 0.05, 0.2, 1.0, 30.0])
def test_hamilton_equations_reproduce_the_six_roots(x):
    om, im = ah.hamilton_frequencies(ah.hamiltonian_matrix(x, *A))
    ref_om = np.array([b["omega"] for b in e2.evaluate(B3, x)["branches"]])
    assert im < 1e-9*np.abs(om).max()
    assert np.allclose(om, ref_om, rtol=1e-9, atol=1e-12*np.abs(ref_om).max())


def test_reduced_hamiltonian_is_positive_definite_on_dense_grid():
    xs = np.geomspace(1e-4, 1e4, 400)
    mins = [np.linalg.eigvalsh(ah.hamiltonian_matrix(x, *A)).min() for x in xs]
    assert min(mins) > 0
    # normalized by k^2 at small k the smallest eigenvalue stays bounded away from 0
    assert min(m/x**2 for m, x in zip(mins, xs) if x < 1e-2) > 0.5


def test_mode_energy_from_hamiltonian_equals_whitham_energy_and_noether():
    for x in (0.01, 0.2, 3.0):
        h = ah.hamiltonian_matrix(x, *A)
        for b in e2.evaluate(B3, x)["branches"]:
            z = en.u_matrix() @ b["vec"]
            E, u = ah.mode_energy_from_h(h, b["omega"], x, *A, z)
            assert E == pytest.approx(b["E"], rel=1e-8)
            sp = ah.energy_split_hamiltonian(u, x, *A)
            assert sum(sp.values()) == pytest.approx(E, rel=1e-9)
            # A-sector share coincides between the two independent representations
            assert sp["A"]/E == pytest.approx(e2.shares(b)["f_proca"], abs=1e-8)


def test_lab_frame_second_order_T00_is_also_positive():
    """T00^(2)_lab = E_shift-sym + w*g_phi<a1 phi_x> + v*g_psi<a1 psi_x>  (second-order charges)."""
    for x in np.geomspace(1e-3, 1e3, 49):
        for b in e2.evaluate(B3, x)["branches"]:
            phi, psi, a0, a1 = en.u_matrix() @ b["vec"]
            dE = (B3.w*G[0]*0.5*np.real(a1*np.conj(1j*x*phi)) + B3.v*G[1]*0.5*np.real(a1*np.conj(1j*x*psi)))
            assert b["E"] + dE > 0


def test_negative_control_complex_background_is_not_positive_definite():
    bad = B3Point(p=0.3394370989189766, s=0.5576612789314608, kappa=-0.43507577121685514, Lam=0.0244906774933632,
                  w=1.1617196379863794, q=1.642277396198334, v=1.2345796836260892, r=-0.5063847645833194, roots=(), record={})
    g = couplings(bad, 1.0); a = (1.0, *g, bad.w, bad.q, bad.v, bad.r, bad.Lam)
    xs = np.geomspace(1e-3, 1e3, 40)
    assert min(np.linalg.eigvalsh(ah.hamiltonian_matrix(x, *a)).min() for x in xs) < 0
    assert max(ah.hamilton_frequencies(ah.hamiltonian_matrix(x, *a))[1] for x in xs) > 1.0


def test_exchange_conclusion_is_stable_under_numerical_resolution():
    res = {}
    for n0 in (61, 241, 481):
        nodes, _ = e2.sweep(B3, n0=n0)
        lab, _ = e2.track(nodes)
        tbl = e2.tracked_table(nodes, lab); xs = np.array([n["x"] for n in nodes])
        ex = e2.exchange_summary(tbl, xs)
        res[n0] = {e["branch"]: e["crossing_x_energy_metric"][0] for e in ex if e["exchanged_energy_metric"]}
    assert set(res[61]) == set(res[241]) == set(res[481]) == {0, 1, 3, 5}
    for b in (0, 1, 3, 5):
        assert res[241][b] == pytest.approx(res[481][b], rel=1e-3)
        assert res[61][b] == pytest.approx(res[481][b], rel=3e-2)


def test_shift_symmetric_energy_equals_second_order_T00_minus_mu_Q_from_nonlinear_action():
    """Independent identity: E_Whitham = h-energy = [eps^2 T00 - w Q2_phi - v Q2_psi] of the NONLINEAR action."""
    from b410 import charges as ch
    import json
    from pathlib import Path
    man = json.loads((Path(__file__).resolve().parent.parent/"data"/"e3_manifest.json").read_text())
    from b410 import experiment3 as e3
    pts = [B3] + [e3.make_point(r) for r in man["selected"] if r["cls"] != "negative_control" and r["index"] in (5713, 5461)]
    for b3 in pts:
        g = couplings(b3, 1.0); P = (1.0, *g, b3.w, b3.q, b3.v, b3.r, b3.Lam)
        for x in (0.05, 0.5, 5.0):
            nd = e2.evaluate(b3, x)
            for b in nd["branches"]:
                z = en.u_matrix() @ b["vec"]
                e = ch.energies(P, z, b["omega"], x)
                assert e["E_shift"] == pytest.approx(b["E"], rel=1e-8, abs=1e-14)
                # lab energy differs by mu.Q2
                assert e["E_lab"] - b["E"] == pytest.approx(b3.w*e["Q2_phi"] + b3.v*e["Q2_psi"], rel=1e-8, abs=1e-14)


def test_negative_shift_energy_branch_has_positive_lab_energy():
    from b410 import charges as ch, experiment3 as e3
    import json
    from pathlib import Path
    man = json.loads((Path(__file__).resolve().parent.parent/"data"/"e3_manifest.json").read_text())
    b3 = e3.make_point([r for r in man["selected"] if r["index"] == 5713 and r["mode"] == "rank1_minus"][0])
    g = couplings(b3, 1.0); P = (1.0, *g, b3.w, b3.q, b3.v, b3.r, b3.Lam)
    found = False
    for x in (0.01, 0.1, 1.0):
        for b in e2.evaluate(b3, x)["branches"]:
            e = ch.energies(P, en.u_matrix() @ b["vec"], b["omega"], x)
            if e["E_shift"] < 0:
                found = True; assert e["E_lab"] > 0
    assert found
