# Full first-order reduction for C1 and C2 including the field amplitudes (review request of 2026-11-13).
# U_t + A U_x = B U  with  U = (a, a_t, a_x, [b, b_t, b_x,] s_or_c, ..._t, ..._x)   (Scalar mediator: 9 variables; HeavyScalar: 6 variables)
# Checks on every background of the Stage-1 sets: (i) A real symmetric (eigenvalues, condition number of the eigenvector matrix, incl. zero eigenvalues of the amplitude rows),
# (ii) consistency of the reduction: max Re eig(-ikA + B) equals gamma from det K at all grid points (extra eigenvalues are the 'constraint' modes u_x - ik u, lambda = 0),
# (iii) energy bound: for any k, d|U|^2/dt = 2 Re<U,(B)U> <= 2 lam_max((B+B^T)/2) |U|^2, so gamma(k) <= lam_max((B+B^T)/2) =: b for all k including k = 0.
# usage: python t3fund_stage1_firstorder.py --manifest e3_manifest.json --out DIR
import numpy as np, json, sys, argparse
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import t3fund_stage1 as T

def AB(model, th):
    if model == 'B':
        g, w, q = th; Z = 1+2*g*g*(w*w-q*q); A = np.zeros((6, 6)); B = np.zeros((6, 6))
        A[1, 2] = A[2, 1] = A[4, 5] = A[5, 4] = -1
        B[0, 1] = 1; B[1, 4] = -2*g*w/Z; B[1, 5] = 2*g*q/Z; B[3, 4] = 1; B[4, 3] = -1; B[4, 1] = 2*g*w; B[4, 2] = -2*g*q
        return A, B
    gp, gs, w, q, v, r = th; sb = gp*w*q+gs*v*r; A = np.zeros((9, 9)); B = np.zeros((9, 9))
    A[1, 1] = 2*gp*sb; A[1, 2] = A[2, 1] = -1; A[4, 4] = 2*gs*sb; A[4, 5] = A[5, 4] = -1; A[7, 8] = A[8, 7] = -1
    B[0, 1] = 1; B[1, 7] = -gp*q; B[1, 8] = -gp*w; B[3, 4] = 1; B[4, 7] = -gs*r; B[4, 8] = -gs*v
    B[6, 7] = 1; B[7, 6] = -1; B[7, 1] = gp*q; B[7, 2] = gp*w; B[7, 4] = gs*r; B[7, 5] = gs*v
    return A, B

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument('--manifest', required=True); ap.add_argument('--out', default='.'); a = ap.parse_args()
    P = T.build(a.manifest); res = {}
    ks = T.KGRID
    for name, L in P.items():
        cond = 0; sym = True; realeig = True; dmax_lo = 0; dmax_hi = 0; bound_viol = 0; bmax = 0; gmax_fo = 0; n0 = set()
        for model, tag, th in L:
            th = [float(x) for x in th]; A, B = AB(model, th)
            sym &= bool(np.allclose(A, A.T)); ev, V = np.linalg.eig(A); realeig &= bool(np.abs(ev.imag).max() < 1e-12)
            cond = max(cond, float(np.linalg.cond(V))); n0.add(int((np.abs(ev) < 1e-12).sum()))
            b = float(np.linalg.eigvalsh(0.5*(B+B.T)).max()); bmax = max(bmax, b)
            gfo = 0
            for k in ks:
                lam = np.linalg.eigvals(-1j*k*A+B); g_fo = max(0.0, lam.real.max())
                rd = T.roots_dbl(model, th, k); g_p = max(0.0, (rd.imag*k).max())
                d = abs(g_fo-g_p)
                if k <= 1e3: dmax_lo = max(dmax_lo, d)
                else: dmax_hi = max(dmax_hi, d/k)
                gfo = max(gfo, g_fo)
            gmax_fo = max(gmax_fo, gfo)
            if gfo > b*(1+1e-9)+1e-12: bound_viol += 1
        res[name] = dict(n=len(L), A_symmetric=sym, A_real_eigs=realeig, max_cond_eigvec=cond, zero_eigenvalue_counts=sorted(n0), max_abs_diff_gamma_k_le_1e3=dmax_lo,
                         max_rel_diff_gamma_over_k_above_1e3=dmax_hi, energy_bound_max=bmax, max_gamma_first_order=gmax_fo, bound_violations=bound_viol)
        print(name, json.dumps(res[name]), flush=True)
    json.dump(res, open(a.out.rstrip('/')+'/t3fund_stage1_firstorder_result.json', 'w'), indent=1)
