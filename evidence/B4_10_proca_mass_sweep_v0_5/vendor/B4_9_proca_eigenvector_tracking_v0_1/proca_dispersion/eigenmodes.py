"""Full coupled 4x4 Fourier pencil and eigenvector-overlap mode tracking.

Field order: (delta_phi, delta_psi, delta_A0, delta_A1).
For real omega,k the matrix is real symmetric. Its determinant is proportional
to the pole-cleared scalar dispersion polynomial used in B4.7/B4.8.
"""
from itertools import permutations
import numpy as np
from .model import full_roots

def coupled_matrix(bg, omega, k):
    bg.validate()
    if k <= 0:
        raise ValueError("k must be positive.")
    M2=bg.M**2
    J0=bg.g_phi*bg.w*bg.q+bg.g_psi*bg.v*bg.r
    Abar=J0/M2
    D1=omega**2-k**2*(1+3*bg.Lambda*bg.q**2)-2*omega*k*Abar*bg.g_phi
    D2=omega**2-k**2-2*omega*k*Abar*bg.g_psi
    j1=bg.g_phi*(bg.w*k-bg.q*omega)
    j2=bg.g_psi*(bg.v*k-bg.r*omega)
    # Proca block in (A0,A1) order, signature (-,+).
    K00=M2+k**2
    K01=omega*k
    K11=omega**2-M2
    return np.array([
        [D1, 0.0, 0.0, j1],
        [0.0, D2, 0.0, j2],
        [0.0, 0.0, K00, K01],
        [j1, j2, K01, K11],
    ],dtype=complex)

def null_vector(bg, omega, k):
    """Unit right singular vector associated with the smallest singular value."""
    H=coupled_matrix(bg,complex(omega),k)
    _,s,vh=np.linalg.svd(H)
    vec=vh[-1].conj().astype(complex)
    # Deterministic phase convention: largest component real and nonnegative.
    pivot=int(np.argmax(np.abs(vec)))
    if abs(vec[pivot])>0:
        vec*=np.exp(-1j*np.angle(vec[pivot]))
    residual = np.linalg.norm(H @ vec) / max(1.0, float(s[0]))
    return vec, float(residual)

def mode_data(bg, omega, k):
    vec,residual=null_vector(bg,omega,k)
    scalar=float(np.sum(np.abs(vec[:2])**2))
    proca=float(np.sum(np.abs(vec[2:])**2))
    return {"vector":vec,"residual":residual,
            "scalar_fraction":scalar,"proca_fraction":proca}

def _best_overlap_assignment(prev_vecs, cur_vecs):
    n=len(prev_vecs)
    best=None; best_score=-1.0
    for perm in permutations(range(n)):
        score=sum(abs(np.vdot(prev_vecs[i],cur_vecs[perm[i]])) for i in range(n))
        if score>best_score:
            best_score=score; best=perm
    return best,float(best_score/n)

def track_eigenmodes(bg,k_values):
    """Track six roots by eigenvector overlap, not by frequency proximity.

    The returned arrays are aligned mode-by-mode along increasing k.
    """
    ks=np.asarray(k_values,float)
    if ks.ndim!=1 or len(ks)<2 or np.any(~np.isfinite(ks)) or np.any(ks<=0) or np.any(np.diff(ks)<=0):
        raise ValueError("k_values must be finite, positive and strictly increasing.")
    roots=np.asarray(full_roots(bg,float(ks[0])),complex)
    order=np.argsort(roots.real)
    roots=roots[order]
    data=[mode_data(bg,z,float(ks[0])) for z in roots]
    vecs=[d["vector"] for d in data]
    n=len(roots)
    freqs=np.empty((len(ks),n),complex); fractions=np.empty((len(ks),n))
    residuals=np.empty((len(ks),n))
    freqs[0]=roots
    for j,d in enumerate(data):
        fractions[0,j]=d["scalar_fraction"]
        residuals[0,j]=d["residual"]
    mean_overlaps=np.ones(len(ks))
    for ik in range(1,len(ks)):
        new_roots=np.asarray(full_roots(bg,float(ks[ik])),complex)
        new_data=[mode_data(bg,z,float(ks[ik])) for z in new_roots]
        new_vecs=[d["vector"] for d in new_data]
        perm,mean_ov=_best_overlap_assignment(vecs,new_vecs)
        mean_overlaps[ik]=mean_ov
        roots=new_roots[list(perm)]
        data=[new_data[j] for j in perm]
        vecs=[d["vector"] for d in data]
        freqs[ik]=roots
        for j,d in enumerate(data):
            fractions[ik,j]=d["scalar_fraction"]
            residuals[ik,j]=d["residual"]
    return {"omega":freqs,"scalar_fraction":fractions,
            "proca_fraction":1.0-fractions,"null_residual":residuals,
            "mean_overlap":mean_overlaps,"k":ks}

def summarize_eigenmode_scan(result):
    om=result["omega"]
    imag=np.abs(om.imag)
    scale=np.maximum(1.0,np.abs(om.real))
    complex_mask=imag>1e-6*scale
    return {
        "n_k":int(len(result["k"])),
        "n_modes":int(om.shape[1]),
        "complex_points":int(np.count_nonzero(complex_mask)),
        "max_abs_imag_omega":float(np.max(imag)),
        "all_real_within_tolerance":bool(not np.any(complex_mask)),
        "max_null_residual":float(np.max(result["null_residual"])),
        "minimum_mean_eigenvector_overlap":float(np.min(result["mean_overlap"][1:])),
        "minimum_mean_eigenvector_overlap_k_over_M":float(result["k"][1+int(np.argmin(result["mean_overlap"][1:]))]),
        "scalar_fraction_range":[float(np.min(result["scalar_fraction"])),float(np.max(result["scalar_fraction"]))],
    }
