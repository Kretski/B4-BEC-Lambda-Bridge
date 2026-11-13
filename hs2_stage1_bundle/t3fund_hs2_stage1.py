# T3-fund HS2 (HeavyScalar, two light fields): C0-C3.  Implements T3_FUND_HS2_SPEC.md, sections 2-4 (registered in commit 7b0e867, SPEC_HASH 9d52c927a9acd3bd, MANIFEST_HASH 6c312456ac894009).
# Implementation notes (not changes of the spec): polynomial Q built with pmul/padd in c; the luminal factor D=1-c^2 is classified separately (roots +-1); k=0 polynomial obtained
# by cofactor expansion of the explicit 3x3 matrix K(omega,0) with polynomial entries; first-order A,B built per the spec (U=(f1,f1_t,f1_x,f2,f2_t,f2_x,xi,xi_t,xi_x)).
import numpy as np, json, sys, time, argparse
import mpmath as mp
from multiprocessing import Pool
TOL_REAL, TOL_UNRES, TOL_CONF = 1e-7, 1e-5, 1e-12
KGRID = np.concatenate([np.geomspace(1e-3, 1e3, 61), np.geomspace(1e3, 1e6, 32)[1:]])
M = 1.0

def pmul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] = out[i+j] + x*y
    return out
def padd(a, b):
    n = max(len(a), len(b)); a = list(a)+[0]*(n-len(a)); b = list(b)+[0]*(n-len(b)); return [x+y for x, y in zip(a, b)]
def pscale(a, s): return [s*x for x in a]

def pieces(th, one=1):
    a, b, e, w, q, v, r = th
    Cp = [[a, e/2], [e/2, b]]; wv = [w, v]; qv = [q, r]
    chi = (a*(w*w-q*q)+b*(v*v-r*r)+e*(w*v-q*r))/(one*M*M)
    Z = [[1+2*chi*Cp[0][0], 2*chi*Cp[0][1]], [2*chi*Cp[1][0], 1+2*chi*Cp[1][1]]]
    return Cp, wv, qv, chi, Z

def Qpoly(th, kv, one=1):
    """Q(c,eps) = (D+eps M^2) D detZ - 4 eps s^T adj(Z) s, coefficient list in c (low->high)."""
    Cp, wv, qv, chi, Z = pieces(th, one)
    eps = one/(kv*kv); dZ = Z[0][0]*Z[1][1]-Z[0][1]*Z[1][0]
    adj = [[Z[1][1], -Z[0][1]], [-Z[1][0], Z[0][0]]]
    s = [[sum(Cp[i][j]*qv[j] for j in range(2)), sum(Cp[i][j]*wv[j] for j in range(2))] for i in range(2)]   # s_i = s_i0 + c s_i1
    ss = [0]
    for i in range(2):
        for j in range(2): ss = padd(ss, pscale(pmul(s[i], s[j]), adj[i][j]))
    D = [one, 0*one, -one]; De = [one+eps*M*M, 0*one, -one]
    return padd(pscale(pmul(De, D), dZ), pscale(ss, -4*eps))

def roots_dbl(th, kv):
    co = [complex(x) for x in Qpoly([float(t) for t in th], float(kv))]
    return np.roots(co[::-1])
def roots_mp(th, kv, dps):
    mp.mp.dps = dps
    t = [mp.mpf(float(x)) for x in th]; k = mp.mpf(float(kv))
    co = Qpoly(t, k, one=mp.mpf(1))
    try: return mp.polyroots(co[::-1], maxsteps=400, extraprec=4*dps)
    except mp.libmp.NoConvergence: return None
def gamma_mp(th, kv, dps):
    rt = roots_mp(th, kv, dps)
    if rt is None: return None
    return max(mp.mpf(0), max(mp.im(c) for c in rt))*mp.mpf(float(kv))      # luminal factor D has real roots; gamma from Q only

def Kmat(th, om, k):
    Cp, wv, qv, chi, Z = pieces(th)
    sg = [om*wv[j]+k*qv[j] for j in range(2)]
    K = np.zeros((3, 3), complex)
    for i in range(2):
        for j in range(2): K[i, j] = Z[i][j]*(k*k-om*om)
        K[i, 2] = -2j*sum(Cp[i][j]*sg[j] for j in range(2))
        K[2, i] = 2j*sum(Cp[j][i]*sg[j] for j in range(2))
    K[2, 2] = k*k-om*om+M*M
    return K
def poly_vs_K(th):
    """max relative difference between det K and k^6 D Q (double) at 6 sample points."""
    worst = 0
    for k, c in ((0.3, 0.4), (1.1, -0.7), (2.5, 1.6), (0.05, 3.0), (7.0, 0.2), (1.0, -1.3)):
        lhs = np.linalg.det(Kmat(th, c*k, k)); D = 1-c*c
        co = [complex(x) for x in Qpoly([float(t) for t in th], k)]; Qv = sum(co[i]*c**i for i in range(len(co)))
        rhs = k**6*D*Qv; worst = max(worst, abs(lhs-rhs)/max(abs(rhs), 1e-300))
    return worst

def first_order(th):
    Cp, wv, qv, chi, Z = pieces([float(t) for t in th]); Cp = np.array(Cp); Zi = np.linalg.inv(np.array(Z)); wv = np.array(wv); qv = np.array(qv)
    A = np.zeros((9, 9)); B = np.zeros((9, 9))
    for blk in (0, 1, 2):
        i0 = 3*blk; A[i0+1, i0+2] = -1.; A[i0+2, i0+1] = -1.
    for i in range(2):
        B[3*i, 3*i+1] = 1.
    zc = Zi@Cp
    for i in range(2):
        B[3*i+1, 7] = -2*(zc@wv)[i]; B[3*i+1, 8] = 2*(zc@qv)[i]
    B[7, 6] = -M*M; B[6, 7] = 1.
    for j in range(2):
        B[7, 3*j+1] = 2*(Cp@wv)[j]; B[7, 3*j+2] = -2*(Cp@qv)[j]
    return A, B
def fo_check(th, kv=0.7):
    A, B = first_order(th); ev = np.linalg.eigvals(kv*A+1j*B)
    cK = np.roots([complex(x) for x in Qpoly([float(t) for t in th], kv)][::-1]); target = list(cK*kv)+[kv, -kv]
    miss = max(min(abs(t-e) for e in ev) for t in target)
    return dict(max_dist_poly_root_to_fo_eig=float(miss), A_symmetric=bool(np.allclose(A, A.T)),
                lam_max_sym_B=float(np.linalg.eigvalsh((B+B.T)/2).max()), cond_eigvec_A=float(np.linalg.cond(np.linalg.eigh(A)[1])))

def k0_poly(th, dps=50):
    """det K(omega,0) by cofactor expansion with polynomial (in omega) entries, 50 digits; returns coefficient list low->high (real parts)."""
    mp.mp.dps = dps; t = [mp.mpf(float(x)) for x in th]; one = mp.mpf(1)
    Cp, wv, qv, chi, Z = pieces(t, one); zero = mp.mpf(0)
    E = [[None]*3 for _ in range(3)]
    for i in range(2):
        for j in range(2): E[i][j] = [zero, zero, -Z[i][j]]
        E[i][2] = [zero, -2j*mp.mpf(1)*sum(Cp[i][j]*wv[j] for j in range(2))] ; E[i][2] = [mp.mpc(x) for x in E[i][2]]
        E[2][i] = [zero, 2j*sum(Cp[j][i]*wv[j] for j in range(2))]; E[2][i] = [mp.mpc(x) for x in E[2][i]]
    E[2][2] = [M*M, zero, -one]
    def mul(*ps):
        r = [one]
        for p in ps: r = pmul(r, p)
        return r
    det = [0]
    for j in range(3):
        cols = [c for c in range(3) if c != j]
        minor = padd(mul(E[1][cols[0]], E[2][cols[1]]), pscale(mul(E[1][cols[1]], E[2][cols[0]]), -1))
        det = padd(det, pscale(mul(E[0][j], minor), (-1)**j))
    return [mp.re(c) for c in det], max(abs(mp.im(c)) for c in det)
def k0_class(th):
    co, imres = k0_poly(th)
    while len(co) > 1 and abs(co[-1]) < mp.mpf(10)**-40: co = co[:-1]
    nz = 0
    while abs(co[nz]) < mp.mpf(10)**-40: nz += 1
    rest = co[nz:]; rt = mp.polyroots(rest[::-1], maxsteps=300, extraprec=200)
    maxr = float(max(abs(mp.im(x))/max(1, abs(x)) for x in rt)) if rt else 0.0
    om2 = sorted(float(mp.re(x)**2) for x in rt)
    # formula cross-check: omega_h^2 = M^2 + 4 s_w^T Z^-1 s_w
    Cp, wv, qv, chi, Z = pieces([float(x) for x in th]); Cp = np.array(Cp); sw = Cp@np.array(wv); oh2 = M*M+4*sw@np.linalg.inv(np.array(Z))@sw
    return dict(zero_root_multiplicity=nz, max_r_nonzero=maxr, im_residual_of_coeffs=float(imres), heavy_omega2_formula=float(oh2),
                nonzero_roots=[complex(x) for x in rt][:4])

def evaluate(args):
    tag, th = args; th = [float(x) for x in th]; out = dict(set=tag, theta=th)
    Cp, wv, qv, chi, Z = pieces(th); ev = np.linalg.eigvalsh(np.array(Z))
    out['C0'] = dict(chibar=float(chi), residual_chi_eq=abs(M*M*chi-(th[0]*(th[3]**2-th[4]**2)+th[1]*(th[5]**2-th[6]**2)+th[2]*(th[3]*th[5]-th[4]*th[6]))), zmin=float(ev[0]))
    if ev[0] <= 0: out['outcome'] = 'FAIL-GHOST'; return out
    out['poly_vs_K_maxrelerr'] = poly_vs_K(th); out['first_order'] = fo_check(th)
    k0 = k0_class(th); out['k0'] = k0; k0_ok = k0['max_r_nonzero'] <= TOL_REAL and k0['nonzero_roots'] and True
    G = []; hidden = 0; conf = []; stuck = 0; nodata = 0
    for kv in KGRID:
        rd = roots_dbl(th, kv); r_d = max(abs(c.imag)/max(1, abs(c)) for c in rd)
        rm = roots_mp(th, kv, 50)
        if rm is None: nodata += 1; continue
        r_m = max(float(abs(mp.im(c)))/max(1.0, float(abs(c))) for c in rm)
        gm = max(0.0, float(max(mp.im(c) for c in rm)))*kv
        if r_d <= TOL_REAL:
            cls = 'real'
            if r_m > TOL_CONF: hidden += 1
        elif r_d <= TOL_UNRES:
            if r_m <= TOL_REAL: cls = 'real'
            elif r_m <= TOL_UNRES: stuck += 1; cls = 'unresolved_stuck'
            else: cls = 'complex' if r_m > TOL_CONF else 'real'
        else: cls = 'complex' if r_m > TOL_CONF else 'real'
        G.append((float(kv), r_d, r_m, gm, cls))
        if cls == 'complex': conf.append(float(kv))
    out.update(n_grid=len(G), n_nodata=nodata, hidden_by_double=hidden, unresolved_stuck=stuck, n_complex=len(conf), complex_at_kmin=bool(G and G[0][4] == 'complex'),
               complex_in_uv_band=bool(any(k >= 1e3 for k in conf)), gamma_max_grid=max((g[3] for g in G), default=0.0))
    verdict = None
    if conf:
        ks = [g[0] for g in G]; idx = [i for i, g in enumerate(G) if g[4] == 'complex']
        lo, hi = max(0, min(idx)-1), min(len(ks)-1, max(idx)+1); fine = []
        for i in range(lo, hi): fine += list(np.geomspace(ks[i], ks[i+1], 5)[:-1])
        fine.append(ks[hi]); rows = []
        for kv in fine:
            g50 = gamma_mp(th, kv, 50); g100 = gamma_mp(th, kv, 100)
            if g50 is None or g100 is None: continue
            mp.mp.dps = 120; rows.append((kv, g50, g100, abs(g50-g100)))
        flags = [0 if r[1] <= r[3] else (1 if r[1] <= 10*r[3] else 2) for r in rows]; run = best = 0
        for f in flags:
            run = run+1 if f == 2 else 0; best = max(best, run)
        signchg = any((a[1] > 0) != (b[1] > 0) for a, b in zip(rows, rows[1:])) and any(f == 1 for f in flags)
        if best >= 3: verdict = 'FAIL-UNSTABLE'
        elif any(f == 1 for f in flags) or signchg or best > 0: verdict = 'INCONCLUSIVE'
        else: verdict = 'PASS'
        mp.mp.dps = 120
        out['P6'] = dict(n_refined=len(fine), max_u_str=mp.nstr(max((r[3] for r in rows), default=mp.mpf(0)), 8), max_consecutive_gt10u=best,
                         min_ratio_gamma_over_u=(mp.nstr(min(r[1]/r[3] for r in rows if r[3] > 0 and r[1] > 0), 6) if any(r[3] > 0 and r[1] > 0 for r in rows) else None),
                         max_gamma=float(max(r[1] for r in rows)) if rows else None)
    c3 = verdict if verdict else ('INCONCLUSIVE' if (stuck or nodata) else 'PASS')
    out['C3_numeric'] = c3
    out['outcome'] = 'FAIL-UNSTABLE' if (c3 == 'FAIL-UNSTABLE' or (c3 == 'PASS' and not k0_ok)) else c3
    return out

def Fdisc(th):
    from t3fund_hs2_manifest_lib import Fcoef
    a, b, e, w, q, v, r = th; A2, A1, A0 = Fcoef(a, b, e, w, q, v, r); return A1*A1-4*A2*A0

def build(man):
    th = lambda d: (d['a'], d['b'], d['e'], d['w'], d['q'], d['v'], d['r'])
    P = []
    P += [('H1_'+d['cs'], th(d)) for d in man['H1']]; P += [('H2_'+d['cs'], th(d)) for d in man['H2']]
    P += [('H3_'+d['cs'], th(d)) for d in man['H3_control']]; P += [('ANCHOR_'+d['label'], th(d)) for d in man['anchors']]
    return P

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--manifest'); ap.add_argument('--out'); ap.add_argument('--only', default=None); ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--procs', type=int, default=4); a = ap.parse_args()
    man = json.load(open(a.manifest)); P = build(man)
    if a.only: P = [p for p in P if p[0].startswith(a.only)]
    if a.limit: P = P[:a.limit]
    t0 = time.time(); print('N', len(P), flush=True)
    with Pool(a.procs) as pool: R = pool.map(evaluate, P, chunksize=4)
    for r, p in zip(R, P):
        r['F_disc'] = None
    json.dump(R, open(a.out, 'w'), indent=1, default=lambda o: str(o)); print('done', time.time()-t0, flush=True)
