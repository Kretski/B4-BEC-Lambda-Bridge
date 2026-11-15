"""Analytic / numerical checks of the quadratic Hamiltonian and energy metric."""
import numpy as np
import pytest
import b410._b49  # noqa: F401  (vendored B4.9 first on sys.path)
from b410 import energy as en, reference as ref
from b410.fixed_eft import load_b3, couplings

B3 = load_b3()


def setup(M=1.0, k=0.7):
    g = couplings(B3, M)
    a = (M, g[0], g[1], B3.w, B3.q, B3.v, B3.r, B3.Lam)
    return a, en.terms(*a), ref.H_blocks(M, k, *g, B3.w, B3.q, B3.v, B3.r, B3.Lam)


@pytest.mark.parametrize("M,k,om", [(1, .3, .4), (1, 1.7, -1.3), (3.0, .2, 5.0), (.2, 4.0, 11.0)])
def test_b49_H_is_hessian_symbol_up_to_constant_unitary(M, k, om):
    a, tl, bl = setup(M, k)
    Hb = en.b49_from_symbol(en.hessian_symbol(tl, om, k), 1j)
    assert np.abs(Hb - ref.H_of(om, bl)).max() < 1e-12*max(1, np.abs(Hb).max())
    # negative control: wrong phase must NOT match
    Hw = en.b49_from_symbol(en.hessian_symbol(tl, om, k), -1j)
    assert np.abs(Hw - ref.H_of(om, bl)).max() > 1e-3


def test_hessian_symbol_is_hermitian_and_a0_has_no_time_derivative():
    a, tl, bl = setup()
    H = en.hessian_symbol(tl, 0.9, 0.7)
    assert np.allclose(H, H.conj().T, atol=1e-14)
    # a0 row/col carries no omega^2 or omega^1 diagonal piece: H[2,2] omega independent
    assert abs(en.hessian_symbol(tl, 0.9, .7)[2, 2] - en.hessian_symbol(tl, 5.0, .7)[2, 2]) < 1e-14
    assert H[2, 2].real == pytest.approx(1.0 + 0.7**2)


@pytest.mark.parametrize("k", [0.05, 0.3, 1.0, 4.0])
def test_whitham_energy_equals_noether_equals_time_domain(k):
    a, tl, bl = setup(1.0, k)
    for om in ref.qz_roots(bl).real:
        x, s1, s2 = en.null_vector_b49(ref.H_of(om, bl))
        Hp = 2*om*bl[2] + bl[1]
        Ew = en.canonical_energy_b49(om, Hp, x)
        z = en.u_matrix() @ x
        assert abs(en.lagrangian_average(tl, z, om, k)) < 1e-12
        En = en.energy_density_terms(tl, z, om, k)
        Et = en.energy_density_time_domain(tl, z, om, k)
        assert En == pytest.approx(Ew, rel=1e-9) and Et == pytest.approx(Ew, rel=1e-9)
        # group split sums back to the total
        assert sum(en.group_energies(*a, om, k, x).values()) == pytest.approx(Ew, rel=1e-9)


@pytest.mark.parametrize("k", [0.3, 1.0])
def test_group_velocity_identity_uses_same_H_prime(k):
    a, tl, bl = setup(1.0, k)
    g = couplings(B3, 1.0)
    h = 1e-5
    rp = ref.qz_roots(ref.H_blocks(1.0, k + h, *g, B3.w, B3.q, B3.v, B3.r, B3.Lam)).real
    rm = ref.qz_roots(ref.H_blocks(1.0, k - h, *g, B3.w, B3.q, B3.v, B3.r, B3.Lam)).real
    fd = (rp - rm)/(2*h)
    for om, v_fd in zip(ref.qz_roots(bl).real, fd):
        x, _, _ = en.null_vector_b49(ref.H_of(om, bl))
        Hp = 2*om*bl[2] + bl[1]
        Hk = (ref.H_of(om, ref.H_blocks(1.0, k + h, *g, B3.w, B3.q, B3.v, B3.r, B3.Lam))
              - ref.H_of(om, ref.H_blocks(1.0, k - h, *g, B3.w, B3.q, B3.v, B3.r, B3.Lam)))/(2*h)
        assert -(x@Hk@x)/(x@Hp@x) == pytest.approx(v_fd, rel=1e-5, abs=1e-7)


def test_energy_sign_is_sign_omega_times_krein_sign():
    a, tl, bl = setup(1.0, 0.8)
    for om in ref.qz_roots(bl).real:
        x, _, _ = en.null_vector_b49(ref.H_of(om, bl))
        Hp = 2*om*bl[2] + bl[1]
        assert np.sign(en.canonical_energy_b49(om, Hp, x)) == np.sign(om)*np.sign(en.krein_value(Hp, x))


def test_a0_is_slaved_by_its_constraint():
    """Row 2 of H z = 0 gives a0 = -(H20 z0 + H21 z1 + H23 z3)/H22; the null vector obeys it."""
    a, tl, bl = setup(1.0, 0.8)
    for om in ref.qz_roots(bl).real:
        H = ref.H_of(om, bl)
        x, _, _ = en.null_vector_b49(H)
        a0 = -(H[2, 0]*x[0] + H[2, 1]*x[1] + H[2, 3]*x[3])/H[2, 2]
        assert a0 == pytest.approx(x[2], abs=1e-12)
