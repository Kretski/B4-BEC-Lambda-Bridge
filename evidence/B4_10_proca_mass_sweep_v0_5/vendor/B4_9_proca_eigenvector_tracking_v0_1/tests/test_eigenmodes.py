import numpy as np
from proca_dispersion import (Background, full_polynomial, coupled_matrix,
    track_eigenmodes, summarize_eigenmode_scan, avoided_crossing_diagnostic)

def b3():
    return Background(M=1.0,g_phi=np.sqrt(0.05),g_psi=0.0,
        w=1.0,q=1.0,v=0.0,r=0.0,Lambda=0.0)

def b1():
    return Background(M=1.0,g_phi=1.0,g_psi=1.0,
        w=np.sqrt(0.5),q=0.0,v=np.sqrt(0.5),r=0.0,Lambda=0.0)

def test_4x4_determinant_matches_pole_cleared_polynomial():
    bg=b3()
    for k,omega in [(0.13,0.41),(0.7,1.2),(3.0,3.14)]:
        det=np.linalg.det(coupled_matrix(bg,omega,k))
        poly=full_polynomial(bg,k)(omega)
        # Determinant equals - the pole-cleared dispersion polynomial.
        assert abs(det+poly) < 2e-9*max(1.0,abs(poly))

def test_eigenvectors_are_null_at_all_real_b3_roots():
    bg=b3()
    result=track_eigenmodes(bg,np.geomspace(1e-3,0.2,30))
    summary=summarize_eigenmode_scan(result)
    assert summary["max_null_residual"] < 1e-8
    assert summary["all_real_within_tolerance"]

def test_eigenvector_tracker_continues_through_avoided_crossing():
    bg=b3()
    ks=np.linspace(2.7,3.6,181)
    result=track_eigenmodes(bg,ks)
    # The two upper positive-frequency branches are scalar/Proca mixed.
    # Eigenvector continuation should remain continuous and exchange dominant
    # scalar-vs-Proca composition across the avoided crossing.
    upper=np.argsort(result["omega"][0].real)[-2:]
    f0=result["scalar_fraction"][0,upper]
    f1=result["scalar_fraction"][-1,upper]
    assert np.min(result["mean_overlap"][1:]) > 0.5
    # The branch that starts scalar-dominated becomes more Proca-dominated,
    # while the partner gains scalar character: an avoided-crossing exchange.
    assert f0[0] > 0.2 and f1[0] < 0.1
    assert f0[1] < 0.1 and f1[1] > 0.1

def test_b1_boundary_eigenmode_scan_remains_real_with_small_null_residual():
    result=track_eigenmodes(b1(),np.geomspace(1e-3,0.2,30))
    summary=summarize_eigenmode_scan(result)
    assert summary["all_real_within_tolerance"]
    assert summary["max_null_residual"] < 1e-7

def test_complex_stress_case_is_flagged_by_eigenmode_scan():
    bg=Background(M=1.0,g_phi=0.8108811308204898,g_psi=1.3198751273700329,
        w=-0.5549437639433696,q=0.39273626882885226,
        v=-1.7629934306179855,r=-0.4494727955570852,Lambda=0.16151817312910333)
    result=track_eigenmodes(bg,np.geomspace(0.4,0.8,20))
    summary=summarize_eigenmode_scan(result)
    assert not summary["all_real_within_tolerance"]


def test_complex_mode_null_vectors_use_full_complex_frequency():
    bg=Background(M=1.0,g_phi=0.8108811308204898,g_psi=1.3198751273700329,
        w=-0.5549437639433696,q=0.39273626882885226,
        v=-1.7629934306179855,r=-0.4494727955570852,Lambda=0.16151817312910333)
    result=track_eigenmodes(bg,[0.49,0.50,0.51])
    assert np.max(result["null_residual"]) < 1e-7
