import numpy as np
import pytest
from proca_dispersion import Background, track_roots, scan_spectrum, local_roots, full_roots, avoided_crossing_diagnostic

def b3():
    return Background(M=1.0, g_phi=np.sqrt(0.05), g_psi=0.0,
                      w=1.0, q=1.0, v=0.0, r=0.0, Lambda=0.0)

def b1():
    return Background(M=1.0, g_phi=1.0, g_psi=1.0,
                      w=np.sqrt(0.5), q=0.0, v=np.sqrt(0.5),
                      r=0.0, Lambda=0.0)

def test_tracks_six_roots_over_log_scan():
    ks=np.geomspace(1e-3,0.3,40)
    tr=track_roots(b3(),ks)
    assert tr.shape==(40,6)
    assert np.all(np.isfinite(tr))

def test_scan_detects_complex_frequency_diagnostics():
    ks=np.geomspace(1e-3,0.1,20)
    tr,summary=scan_spectrum(b3(),ks)
    assert summary["n_modes"]==6
    assert summary["complex_points"]==0
    assert summary["all_real_within_tolerance"]

def test_b1_boundary_scan_has_six_finite_roots():
    ks=np.geomspace(1e-3,0.2,30)
    tr,summary=scan_spectrum(b1(),ks)
    assert tr.shape==(30,6)
    assert summary["all_real_within_tolerance"]

def test_first_low_energy_roots_approach_local_eft_roots():
    bg=b3()
    k=1e-5
    roots=full_roots(bg,k)
    low=roots[np.argsort(np.abs(roots))[:4]]/k
    assert np.max(np.abs(np.sort_complex(low)-np.sort_complex(local_roots(bg))))<2e-4

def test_root_tracker_rejects_nonmonotonic_k():
    with pytest.raises(ValueError):
        track_roots(b3(),[0.01,0.02,0.015])

def test_root_tracker_rejects_zero_k():
    with pytest.raises(ValueError):
        track_roots(b3(),[0.0,0.01])


def test_avoided_crossing_is_resolved_by_dense_k_scan():
    bg=b3()
    ks=np.linspace(2.7,3.6,181)
    gaps=[]
    for k in ks:
        roots=full_roots(bg,float(k))
        pos=np.sort(roots.real[(roots.real>0)&(np.abs(roots.imag)<1e-8)])
        # The two upper positive branches are the Proca-like and
        # superluminal scalar-like modes in this interval.
        gaps.append(pos[-1]-pos[-2])
    idx=int(np.argmin(gaps))
    assert 2.9 < ks[idx] < 3.3
    assert gaps[idx] > 1e-4
    assert gaps[idx] < 0.06

def test_complex_frequency_detector_on_out_of_domain_stress_case():
    bg=Background(M=1.0,g_phi=0.8108811308204898,g_psi=1.3198751273700329,
        w=-0.5549437639433696,q=0.39273626882885226,
        v=-1.7629934306179855,r=-0.4494727955570852,Lambda=0.16151817312910333)
    roots=full_roots(bg,0.5)
    assert np.max(np.abs(roots.imag))>0.1


def test_avoided_crossing_diagnostic_locates_small_nonzero_gap():
    result=avoided_crossing_diagnostic(b3(),2.7,3.6,n=901)
    assert 2.9 < result["k_over_M"] < 3.3
    assert 1e-4 < result["gap_omega"] < 0.06
