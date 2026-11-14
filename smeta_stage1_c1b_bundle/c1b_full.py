"""
Full C1(b) diagnostic for SM+eta Stage 1 (post-hoc, diagnostic only; does NOT change any registered criterion or C3 outcome).

Case set: every case in smeta_stage1_result.json with C1.max_dist_poly_root_to_fo_eig > 1e-9 (215 cases), plus --controls N
randomly chosen cases that PASSED the registered check (b), as a sanity control.

For each case and for k in (0.7, 3.0, 50.0):
  (i)   roots of det K (c-polynomial, 60 digits) times k  vs  eigenvalues of the first-order reduction (60 digits, mp.eig):
        omega = i*lambda, lambda = eig(i k P + C).  Reported: max_dist_mp (multiset distance, both directions).
  (ii)  same comparison in double precision (as registered): max_dist_dbl.
  (iii) root-cluster structure of det K at k: largest cluster of roots within 1e-6 (multiplicity m) and the spread inside
        clusters (<=1e-25 => exactly multiple at 60 digits); predicted double-precision error scale eps_mach^(1/m).
Verdict per case: 'MP-CONSISTENT' if max_dist_mp <= 1e-20 at all three k, else 'MP-MISMATCH'.
sensitivity_explained = MP-CONSISTENT and exact multiple roots (m>=2) and the registered double distance at k=0.7 <= 1e3*eps_mach^(1/m).
Note: the double-precision discrepancy grows with k (observed ~k^2), so it is reported per k; only k=0.7 (the registered k) enters the flag.
Resumable: results are appended to c1b_full.jsonl (one line per case); rerun the same command to continue.
Usage: python c1b_full.py --procs 2 [--controls 10] [--limit N]
"""
import json, sys, os, argparse, random
import numpy as np, mpmath as mp
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t3fund_smeta_stage1 as S
DPS = 60
KS = (0.7, 3.0, 50.0)
EPS = np.finfo(float).eps

def first_order_mp(H):
    n = 12; A0 = mp.zeros(n, n); A1 = mp.zeros(n, n); B = mp.zeros(n, n)
    for a in range(4):
        r1, r2, r3 = 3*a, 3*a+1, 3*a+2; v, t, x = 3*a, 3*a+1, 3*a+2
        A0[r1, v] = 1; B[r1, t] = -1
        A0[r2, x] = 1; A1[r2, t] = -1
        for b in range(4):
            vb, tb, xb = 3*b, 3*b+1, 3*b+2
            A0[r3, tb] += -H[3*a+1][tb]; A0[r3, xb] += -H[3*a+1][xb]
            A1[r3, tb] += -H[3*a+2][tb]; A1[r3, xb] += -H[3*a+2][xb]
            B[r3, vb] += H[3*a][vb]; B[r3, tb] += H[3*a][tb] - H[3*a+1][vb]; B[r3, xb] += H[3*a][xb] - H[3*a+2][vb]
    return A0, A1, B

def multiset_dist(a, b):
    """max over a of min over b, and vice versa (both directions)."""
    d1 = max(min(abs(x-y) for y in b) for x in a); d2 = max(min(abs(x-y) for x in a) for y in b)
    return max(d1, d2)

def clusters(roots, tol=1e-6):
    r = sorted(roots, key=lambda z: (float(mp.re(z)), float(mp.im(z)))); used = [False]*len(r); out = []
    for i in range(len(r)):
        if used[i]: continue
        cl = [i]; used[i] = True
        for j in range(i+1, len(r)):
            if not used[j] and abs(r[j]-r[i]) < tol: cl.append(j); used[j] = True
        out.append(cl)
    m = max(len(c) for c in out); sep = None
    for c in out:
        if len(c) >= 2:
            s = max(abs(r[a]-r[b]) for a in c for b in c)  # spread inside the cluster
            sep = s if sep is None else max(sep, s)
    return m, (float(sep) if sep is not None else None)

def work(item):
    case, cid, dbl_dist = item
    mp.mp.dps = DPS
    g, f, bg, Hmp, one = S.to_mp(case, DPS)
    A0, A1, B = first_order_mp(Hmp)
    P = -(A0**-1)*A1; C = -(A0**-1)*B
    Hd = [[float(x) for x in row] for row in Hmp]
    out = dict(_cid=cid, pop=case['pop'], dist_dbl_registered=dbl_dist, k={})
    worst_mp = 0; worst_dbl = 0; mmax = 1
    for k in KS:
        kk = mp.mpf(k)
        # det K roots (omega = c*k)
        d = S.poly_c(Hmp, kk, False); co = [mp.re(x) for x in d]
        while len(co) > 1 and abs(co[-1]) < mp.mpf(10)**(-DPS+10)*max(abs(x) for x in co): co = co[:-1]
        rts = mp.polyroots(co[::-1], maxsteps=2000, extraprec=8*DPS)
        wk = [r*kk for r in rts]
        # first-order eigenvalues: -i omega U = i k P U + C U  ->  omega = i * eig(i k P + C)
        Mx = mp.mpc(0, 1)*kk*P + C
        E, _ = mp.eig(Mx, left=False, right=False) if False else (mp.eig(Mx, left=False, right=False), None)
        wfo = [mp.mpc(0, 1)*e for e in E]
        # FO has extra zero-frequency modes; compare against nonzero det K roots both directions only on the det K side and on FO nonzero side
        wfo_nz = [w for w in wfo if abs(w) > mp.mpf(10)**(-20)]
        dm = multiset_dist(wk, wfo_nz) if wfo_nz else None
        # double precision as registered
        rd, _ = S.roots_dbl(Hd, k); wd = list(rd*k)
        c1 = S.c1_info(Hd, k); ev = list(1j*np.linalg.eigvals(1j*k*(-np.linalg.inv(np.array(S.first_order(Hd)[0]))@S.first_order(Hd)[1]) + (-np.linalg.inv(S.first_order(Hd)[0])@S.first_order(Hd)[2])))
        ev_nz = [e for e in ev if abs(e) > 1e-9]
        dd = max(max(min(abs(a-b) for b in ev_nz) for a in wd), max(min(abs(a-b) for a in wd) for b in ev_nz))
        m, spread = clusters(wk)
        out['k'][str(k)] = dict(n_detK_roots=len(wk), n_fo_nonzero=len(wfo_nz), dist_mp=float(dm) if dm is not None else None,
                                dist_dbl=float(dd), max_cluster_mult=m, cluster_spread=spread)
        worst_mp = max(worst_mp, float(dm) if dm is not None else 1e9); worst_dbl = max(worst_dbl, float(dd)); mmax = max(mmax, m)
    out['max_dist_mp'] = worst_mp; out['max_dist_dbl'] = worst_dbl; out['max_cluster_mult'] = mmax
    out['predicted_dbl_scale'] = float(EPS**(1.0/mmax))
    out['verdict'] = 'MP-CONSISTENT' if worst_mp <= 1e-20 else 'MP-MISMATCH'
    sp = [v['cluster_spread'] for v in out['k'].values() if v['cluster_spread'] is not None]
    out['exact_multiple_roots'] = bool(mmax >= 2 and sp and max(sp) <= 1e-25)   # spread at the 60-digit sqrt-precision floor => exactly multiple
    out['dist_dbl_k07'] = out['k']['0.7']['dist_dbl']
    out['sensitivity_explained'] = bool(out['verdict'] == 'MP-CONSISTENT' and out['exact_multiple_roots'] and out['dist_dbl_k07'] <= 1e3*out['predicted_dbl_scale'])
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--procs', type=int, default=2); ap.add_argument('--controls', type=int, default=10); ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--result', default='smeta_stage1_result.json'); ap.add_argument('--manifest', default='t3fund_smeta_manifest.json'); ap.add_argument('--out', default='c1b_full.jsonl')
    a = ap.parse_args()
    R = json.load(open(a.result)); man = json.load(open(a.manifest))
    P = {json.dumps([p['pop'], p['cs'], p.get('idx'), p['f'], p.get('label')]): p for p in S.build(man)}
    bad = [r for r in R if 'C1' in r and r['C1']['max_dist_poly_root_to_fo_eig'] > 1e-9]
    good = [r for r in R if 'C1' in r and r['C1']['max_dist_poly_root_to_fo_eig'] <= 1e-9]
    random.seed(20261114); ctrl = random.sample(good, min(a.controls, len(good)))
    items = [(P[r['_cid']], r['_cid'], r['C1']['max_dist_poly_root_to_fo_eig']) for r in bad + ctrl]
    if a.limit: items = items[:a.limit]
    done = set()
    if os.path.exists(a.out): done = {json.loads(l)['_cid'] for l in open(a.out)}
    todo = [it for it in items if it[1] not in done]
    print('cases', len(items), '(failing', len(bad), '+ controls', len(ctrl), ') todo', len(todo), flush=True)
    with Pool(a.procs) as pool, open(a.out, 'a') as fh:
        for i, x in enumerate(pool.imap_unordered(work, todo, chunksize=1)):
            fh.write(json.dumps(x) + '\n'); fh.flush()
            if (i+1) % 10 == 0: print('progress', i+1, len(todo), flush=True)
    X = [json.loads(l) for l in open(a.out)]
    import collections
    print('SUMMARY failing:', collections.Counter((x['verdict'], x['sensitivity_explained']) for x in X if x['dist_dbl_registered'] > 1e-9))
    print('SUMMARY controls:', collections.Counter((x['verdict'], x['sensitivity_explained']) for x in X if x['dist_dbl_registered'] <= 1e-9))
    print('max_dist_mp overall', max(x['max_dist_mp'] for x in X), 'max cluster mult', max(x['max_cluster_mult'] for x in X))
