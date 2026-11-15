import numpy as np
from b410 import uv
from b410.fixed_eft import load_b3

def test_identity_g_rho_equals_sqrtX():
    b3 = load_b3()
    for M in (0.5, 7.0, 300.0):
        inv = uv.invariants(b3, M)
        for k in inv["sqrtX"]:
            assert abs(inv["g_times_rho"][k] - inv["sqrtX"][k]) < 1e-12

def test_anchor_controlled_flags():
    c = uv.controlled(load_b3())
    assert abs(c["sqrtX_max"] - np.sqrt(0.9106)) < 1e-3
    assert not c["A_g1_rho0.3"] and not c["B_g1_rho0.1"] and c["C_gStrong_rho0.3"]
    assert c["g_needed_for_rho0.3"] > 3.0

def test_x_window_monotone():
    w = uv.x_window(load_b3(), n=120)
    assert w[0.01] is not None and w[0.1] is not None and w[0.01] <= w[0.1]
