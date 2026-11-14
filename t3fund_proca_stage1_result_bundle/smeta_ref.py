import os
# T3-fund SM + dynamical eta: C0-C3 (and C1 consistency). Implements T3_FUND_SMETA_SPEC.md (registered in commit d006f77; SPEC_HASH b2fbbf16e6741ac7, MANIFEST_HASH 9e9eb854437015e0).
# Implementation notes (not changes of the spec): det K(ck,k)/k^8 is built by exact Leibniz expansion with polynomial entries in c (entries scaled by 1/k^2, eps=1/k);
# the 50-digit root computation (P3 confirmation) is performed at grid points where the double-precision ratio r_d > 1e-7 and at every 10th grid point (hidden-complex diagnostic); at the other
# points the class is 'real' from double precision (identical classification to computing it everywhere, since class 'real' for r_d<=1e-7 does not depend on the multiprecision result).
import numpy as np, json, sys, time, argparse
import mpmath as mp
from multiprocessing import Pool
TOL_REAL, TOL_UNRES, TOL_CONF, TAU = 1e-7, 1e-5, 1e-12, 1e-12
KGRID = np.concatenate([np.geomspace(1e-3, 1e3, 61), np.geomspace(1e3, 1e6, 32)[1:]])
M = 1.0
FIELDS = 4          # phi, psi, sigma, eta ; jets (v,t,x)

def pmul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] = out[i+j] + x*y
    return out
def padd(a, b):
    n = max(len(a), len(b)); a = list(a)+[0]*(n-len(a)); b = list(b)+[0]*(n-len(b)); return [x+y for x, y in zip(a, b)]
def pscale(a, s): return [s*x for x in a]

# ---------- background (branch S / Z) and Hessian on jets, generic number type via `lib` (math for float, mp for mpf) ----------
def background(g, w, q, v, r, lib, one):
    al = [w, v]; be = [q, r]
    A = g[0]*w*q + g[1]*v*r; B = one/2*(g[0]*(w*w+q*q) + g[1]*(v*v+r*r))
    if abs(B) < abs(A): br = 'S'; sig = (1 if A > 0 else -1)*lib.sqrt(A*A-B*B); eta = lib.atanh(-B/A)/2
    elif abs(A) < abs(B): br = 'Z'; sig = 0*one; eta = lib.atanh(-A/B)/2
    else: raise ValueError('boundary')
    return dict(al=al, be=be, A=A, B=B, br=br, sig=sig, eta=eta)
def hessian(g, f, bg, lib, one):
    """12x12 jet Hessian, index 3*a+j (a: phi,psi,sigma,eta ; j: v,t,x), closed form on the constant background."""
    sig, eta = bg['sig'], bg['eta']; s2 = lib.sinh(2*eta); c2 = lib.cosh(2*eta)
    H = [[0*one]*12 for _ in range(12)]
    def put(i, j, val): H[i][j] = val; H[j][i] = val
    Ueta_sum = 0*one; Uetaeta_sum = 0*one
    for i in range(2):
        a, b = bg['al'][i], bg['be'][i]; gi = g[i]
        H[3*i+1][3*i+1] = one + sig*gi*s2; H[3*i+2][3*i+2] = -one + sig*gi*s2; put(3*i+1, 3*i+2, sig*gi*c2)
        put(6, 3*i+1, gi*(s2*a + c2*b)); put(6, 3*i+2, gi*(s2*b + c2*a))                       # sigma value - phi_t, phi_x
        put(9, 3*i+1, sig*gi*(2*c2*a + 2*s2*b)); put(9, 3*i+2, sig*gi*(2*c2*b + 2*s2*a))        # eta value - phi_t, phi_x
        Ueta_sum += gi*(c2*(a*a+b*b) + 2*s2*a*b); Uetaeta_sum += gi*(2*s2*(a*a+b*b) + 4*c2*a*b)
    H[6][6] = -M*M*one; H[7][7] = one; H[8][8] = -one
    H[9][9] = sig*Uetaeta_sum; H[10][10] = f; H[11][11] = -f
    put(6, 9, Ueta_sum)
    return H
def kinetic_values(g, f, bg, lib, one):
    s2 = lib.sinh(2*bg['eta']); return [one + bg['sig']*g[0]*s2, one + bg['sig']*g[1]*s2, one, f]

# ---------- det K(ck,k)/k^8 as polynomial in c ----------
def Phi(eps, i_unit):
    """Phi[j][l](c) as coefficient lists (low->high in c) for jets j,l in (v,t,x), scaled by 1/k^2; eps=1/k."""
    z = 0*eps; I = i_unit; one = z+1
    return {(0, 0): [eps*eps], (0, 1): [z, -I*eps], (1, 0): [z, I*eps], (0, 2): [I*eps], (2, 0): [-I*eps],
            (1, 1): [z, z, one], (1, 2): [z, -one], (2, 1): [z, -one], (2, 2): [one]}
def Kentries(H, eps, i_unit):
    P = Phi(eps, i_unit); E = [[None]*4 for _ in range(4)]
    for a in range(4):
        for b in range(4):
            acc = [0]
            for j in range(3):
                for l in range(3):
                    h = H[3*a+j][3*b+l]
                    if h != 0: acc = padd(acc, pscale(P[(j, l)], h))
            E[a][b] = acc
    return E
PERMS = None
def det4(E):
    global PERMS
    if PERMS is None:
        import itertools
        def sgn(p):
            s = 1; p = list(p)
            for i in range(len(p)):
                for j in range(i+1, len(p)):
                    if p[i] > p[j]: s = -s
            return s
        PERMS = [(p, sgn(p)) for p in itertools.permutations(range(4))]
    det = [0]
    for p, s in PERMS:
        t = E[0][p[0]]
        for r in (1, 2, 3): t = pmul(t, E[r][p[r]])
        det = padd(det, pscale(t, s))
    return det
def poly_c(H, k, dbl=True):
    if dbl: eps = 1.0/k; I = 1j
    else: eps = mp.mpf(1)/k; I = mp.mpc(0, 1)
    d = det4(Kentries(H, eps, I)); return d

def roots_from(poly):
    re = [complex(mp.re(x)) if not isinstance(x, (int, float, complex)) else complex(x) for x in poly]
    return re
def roots_dbl(H, k):
    d = poly_c(H, k, True); co = np.array([complex(x) for x in d]); imres = np.abs(co.imag).max()/max(np.abs(co).max(), 1e-300)
    co = co.real
    while len(co) > 1 and abs(co[-1]) < 1e-14*np.abs(co).max(): co = co[:-1]
    return np.roots(co[::-1]), imres
def roots_mp(Hmp, k, dps):
    mp.mp.dps = dps; d = poly_c(Hmp, mp.mpf(k), False); co = [mp.re(x) for x in d]
    try: return mp.polyroots(co[::-1], maxsteps=500, extraprec=6*dps)
    except mp.libmp.NoConvergence: return None

def to_mp(case, dps):
    mp.mp.dps = dps; one = mp.mpf(1); g = [mp.mpf(case['gphi']), mp.mpf(case['gpsi'])]; f = mp.mpf(case['f'])
    w, q, v, r = [mp.mpf(case[x]) for x in 'wqvr']; bg = background(g, w, q, v, r, mp, one); return g, f, bg, hessian(g, f, bg, mp, one), one

# ---------- k = 0 ----------
def k0_info(Hmp, dps=50):
    mp.mp.dps = dps; one = mp.mpf(1); I = mp.mpc(0, 1); z = mp.mpf(0)
    P = {(0, 0): [one], (0, 1): [z, -I], (1, 0): [z, I], (1, 1): [z, z, one]}
    E = [[None]*4 for _ in range(4)]
    for a in range(4):
        for b in range(4):
            acc = [0]
            for j in range(2):
                for l in range(2):
                    h = Hmp[3*a+j][3*b+l]
                    if h != 0: acc = padd(acc, pscale(P[(j, l)], h))
            E[a][b] = acc
    d = det4(E); co = [mp.re(x) for x in d]; imres = max(abs(mp.im(x)) for x in d)
    scale = max(abs(x) for x in co)
    while len(co) > 1 and abs(co[-1]) < scale*mp.mpf(10)**-40: co = co[:-1]
    nz = 0
    while abs(co[nz]) < scale*mp.mpf(10)**-40: nz += 1
    rest = co[nz:]; rt = mp.polyroots(rest[::-1], maxsteps=500, extraprec=300)
    maxr = float(max(abs(mp.im(x))/max(1, abs(x)) for x in rt)) if rt else 0.0
    x2 = sorted(set(round(float(mp.re(x)**2), 10) for x in rt))
    return dict(zero_root_multiplicity=nz, max_r_nonzero=maxr, im_residual_of_coeffs=float(imres), omega2_unique=x2, n_nonzero_roots=len(rt))
def k0_closed_form(case, bg, H, f):
    """closed-form k=0 class (branch S, Z_t>0): returns ('real'|'complex', sorted x roots)."""
    g = [case['gphi'], case['gpsi']]; A, B, sig = bg['A'], bg['B'], bg['sig']; al = bg['al']; be = bg['be']
    Zt = [1-g[i]*B for i in range(2)]; N = [2*g[i]*(A*al[i]-B*be[i]) for i in range(2)]; m = [g[i]*(A*be[i]-B*al[i])/sig for i in range(2)]
    S1 = sum(m[i]**2/Zt[i] for i in range(2)); S2 = sum(N[i]**2/Zt[i] for i in range(2)); S3 = sum(m[i]*N[i]/Zt[i] for i in range(2))
    a2 = f; a1 = (4*sig**2-S2)-f*(1+S1); a0 = -(1+S1)*(4*sig**2-S2)-S3**2; disc = a1*a1-4*a2*a0
    if disc < 0: return 'complex', None
    x = sorted([(-a1-np.sqrt(disc))/(2*a2), (-a1+np.sqrt(disc))/(2*a2)])
    return ('real' if x[0] >= 0 else 'complex'), [float(t) for t in x]

# ---------- first-order reduction (C1) ----------
def first_order(H):
    """U=(Phi_a, Phi_a,t, Phi_a,x); A0 U_t + A1 U_x + B U = 0 (rows: Phi_t def, compatibility, linearised EL)."""
    H = np.array(H, dtype=float); n = 12; A0 = np.zeros((n, n)); A1 = np.zeros((n, n)); Bm = np.zeros((n, n))
    for a in range(4):
        r1, r2, r3 = 3*a, 3*a+1, 3*a+2; v, t, x = 3*a, 3*a+1, 3*a+2
        A0[r1, v] = 1.0; Bm[r1, t] = -1.0                       # d_t Phi - Phi_t = 0
        A0[r2, x] = 1.0; A1[r2, t] = -1.0                       # d_t Phi_x - d_x Phi_t = 0
        for b in range(4):
            vb, tb, xb = 3*b, 3*b+1, 3*b+2
            A0[r3, tb] += -H[3*a+1, tb]; A0[r3, xb] += -H[3*a+1, xb]
            A1[r3, tb] += -H[3*a+2, tb]; A1[r3, xb] += -H[3*a+2, xb]
            Bm[r3, vb] += H[3*a, vb]; Bm[r3, tb] += H[3*a, tb] - H[3*a+1, vb]; Bm[r3, xb] += H[3*a, xb] - H[3*a+2, vb]
    return A0, A1, Bm
def c1_info(H, kv=0.7):
    A0, A1, Bm = first_order(H); out = {}
    try:
        A0i = np.linalg.inv(A0)
    except np.linalg.LinAlgError: return dict(error='A0 singular')
    P = -A0i@A1; Cs = -A0i@Bm; lam, V = np.linalg.eig(P); cond = np.linalg.cond(V)
    Ct = np.linalg.solve(V, Cs@V); Hm = (Ct+Ct.conj().T)/2; gmax = float(np.linalg.eigvalsh(Hm).max())
    ev = np.linalg.eigvals(-1j*np.zeros((n:=12, n)) + 0)   # placeholder to keep names local
    # eigenvalues omega of (-i omega A0 + i k A1 + B) U = 0  <=>  omega = eig of (k*(-A0^-1 A1)*... ) ; U_t = P U_x + C U -> -i omega U = i k P U + C U
    Mx = 1j*kv*P + Cs; w = 1j*np.linalg.eigvals(Mx)   # omega = i * lambda
    return dict(max_imag_lambdaP=float(np.abs(lam.imag).max()), cond_V=float(cond), gamma_max=gmax, eig_omega_k07=w, k07=kv)

def evaluate(case):
    out = dict(case=dict(cs=case['cs'], idx=case.get('idx'), pop=case['pop'], f=case['f'], label=case.get('label')))
    g, f, bg, Hmp, one = to_mp(case, 50)
    out['branch'] = bg['br']; out['sigma'] = float(bg['sig']); out['eta'] = float(bg['eta'])
    # C0 residuals of EL (sigma and eta equations; phi,psi trivial) at 50 digits
    s2 = mp.sinh(2*bg['eta']); c2 = mp.cosh(2*bg['eta']); A, B, sig = bg['A'], bg['B'], bg['sig']
    J = A*c2+B*s2; Jp = 2*(B*c2+A*s2)
    out['C0'] = dict(res_sigma=float(abs(M*M*sig-J)), res_eta=float(abs(sig*Jp)))
    kap = kinetic_values(g, f, bg, mp, one); out['kappa'] = [float(x) for x in kap]
    if min(kap) < -TAU: out['outcome'] = 'FAIL-GHOST'; return out
    if min(abs(x) for x in kap) <= TAU: out['outcome'] = 'DEGENERATE'; return out
    Hd = [[float(x) for x in row] for row in Hmp]
    out['C1'] = c1_info(Hd); w1 = out['C1'].pop('eig_omega_k07'); k07 = out['C1'].pop('k07')
    rd, ir = roots_dbl(Hd, k07); target = list(rd*k07); miss = max(min(abs(t_-e) for e in w1) for t_ in target)
    out['C1']['max_dist_poly_root_to_fo_eig'] = float(miss); out['C1']['coef_imag_rel'] = float(ir)
    k0 = k0_info(Hmp); out['k0'] = k0; k0_ok = k0['max_r_nonzero'] <= TOL_REAL
    if bg['br'] == 'S':
        cls, xs = k0_closed_form(case, dict(A=float(A), B=float(B), sig=float(sig), al=[float(x) for x in bg['al']], be=[float(x) for x in bg['be']]), Hd, float(f))
        out['k0_closed'] = dict(cls=cls, x=xs, numeric_cls='real' if k0_ok else 'complex', agree=bool(cls == ('real' if k0_ok else 'complex')))
        if xs is not None and k0_ok:
            nx = k0['omega2_unique']; out['k0_closed']['max_x_err'] = float(max(min(abs(a-b)/max(1, abs(b)) for b in nx) for a in xs)) if nx else None
    G = []; hidden = 0; conf = []; stuck = 0; nodata = 0; gam_dbl_max = 0.0; speeds_uv = None
    for ik, kv in enumerate(KGRID):
        rd, ir = roots_dbl(Hd, kv); r_d = float(max(abs(c.imag)/max(1, abs(c)) for c in rd)); gd = float(max(0.0, rd.imag.max())*kv)
        if ik == len(KGRID)-1: speeds_uv = sorted(float(c.real) for c in rd) if r_d <= TOL_REAL else None
        need_mp = (r_d > TOL_REAL) or (ik % 10 == 0)
        if need_mp:
            rm = roots_mp(Hmp, kv, 50)
            if rm is None: nodata += 1; continue
            r_m = max(float(abs(mp.im(c)))/max(1.0, float(abs(c))) for c in rm); gm = max(0.0, float(max(mp.im(c) for c in rm)))*kv
        else: r_m, gm = r_d, gd
        if r_d <= TOL_REAL:
            cls = 'real'
            if r_m > TOL_CONF: hidden += 1
        elif r_d <= TOL_UNRES:
            if r_m <= TOL_REAL: cls = 'real'
            elif r_m <= TOL_UNRES: stuck += 1; cls = 'unresolved_stuck'
            else: cls = 'complex' if r_m > TOL_CONF else 'real'
        else: cls = 'complex' if r_m > TOL_CONF else 'real'
        G.append((float(kv), r_d, r_m, gm, cls)); gam_dbl_max = max(gam_dbl_max, gd)
        if cls == 'complex': conf.append(float(kv))
    out['C1']['max_gamma_grid_dbl'] = gam_dbl_max; out['C1']['gamma_within_bound'] = bool(gam_dbl_max <= out['C1'].get('gamma_max', 0)*(1+1e-9)+1e-9)
    # Pred-3 data: UV principal speeds (branch S, real roots at k=1e6)
    if bg['br'] == 'S' and speeds_uv is not None:
        gl = [float(x) for x in g]; Zt = [float(1-gl[i]*B) for i in range(2)]; Y = [float(gl[i]*A) for i in range(2)]; R = [np.sqrt(1+(gl[i]*float(sig))**2) for i in range(2)]
        pred = sorted([1., -1., 1., -1.]+[(Y[i]+s*R[i])/Zt[i] for i in range(2) for s in (1, -1)])
        out['uv_speeds'] = dict(observed=speeds_uv, predicted=[float(x) for x in pred], max_err=float(max(abs(a-b) for a, b in zip(speeds_uv, pred))), front=float(max(abs(x) for x in pred)))
    out.update(n_grid=len(G), n_nodata=nodata, hidden_by_double=hidden, unresolved_stuck=stuck, n_complex=len(conf), complex_at_kmin=bool(G and G[0][4] == 'complex'),
               complex_in_uv_band=bool(any(k >= 1e3 for k in conf)), gamma_max_grid=max((g_[3] for g_ in G), default=0.0))
    verdict = None
    if conf:
        ks = [g_[0] for g_ in G]; idx = [i for i, g_ in enumerate(G) if g_[4] == 'complex']
        lo, hi = max(0, min(idx)-1), min(len(ks)-1, max(idx)+1); fine = []
        for i in range(lo, hi): fine += list(np.geomspace(ks[i], ks[i+1], 5)[:-1])
        fine.append(ks[hi]); rows = []
        def gam(kv, dps):
            rt = roots_mp(Hmp, kv, dps)
            return None if rt is None else max(mp.mpf(0), max(mp.im(c) for c in rt))*mp.mpf(float(kv))
        for kv in fine:
            g50 = gam(kv, 50); g100 = gam(kv, 100)
            if g50 is None or g100 is None: continue
            mp.mp.dps = 120; rows.append((kv, g50, g100, abs(g50-g100)))
        flags = [0 if r_[1] <= r_[3] else (1 if r_[1] <= 10*r_[3] else 2) for r_ in rows]; run = best = 0
        for fl in flags:
            run = run+1 if fl == 2 else 0; best = max(best, run)
        signchg = any((a[1] > 0) != (b[1] > 0) for a, b in zip(rows, rows[1:])) and any(fl == 1 for fl in flags)
        if best >= 3: verdict = 'FAIL-UNSTABLE'
        elif any(fl == 1 for fl in flags) or signchg or best > 0: verdict = 'INCONCLUSIVE'
        else: verdict = 'PASS'
        mp.mp.dps = 120
        out['P6'] = dict(n_refined=len(fine), max_u_str=mp.nstr(max((r_[3] for r_ in rows), default=mp.mpf(0)), 8), max_consecutive_gt10u=best,
                         min_ratio_gamma_over_u=(mp.nstr(min(r_[1]/r_[3] for r_ in rows if r_[3] > 0 and r_[1] > 0), 6) if any(r_[3] > 0 and r_[1] > 0 for r_ in rows) else None),
                         max_gamma=float(max(r_[1] for r_ in rows)) if rows else None)
    c3 = verdict if verdict else ('INCONCLUSIVE' if (stuck or nodata) else 'PASS')
    out['C3_numeric'] = c3
    out['outcome'] = 'FAIL-UNSTABLE' if (c3 == 'FAIL-UNSTABLE' or (c3 == 'PASS' and not k0_ok)) else c3
    return out

def build(man):
    P = []
    S = {(d['cs'], d['idx']): d for d in man['draws']}
    def mk(d, f, pop, label=None):
        c = dict(cs=d['cs'], idx=d.get('idx'), gphi=d['gphi'], gpsi=d['gpsi'], w=d['w'], q=d['q'], v=d['v'], r=d['r'], f=f, pop=pop, label=label); return c
    for d in man['draws']:
        if d['pop'] == 'S1':
            for f in man['fmain']: P.append(mk(d, f, 'S1'))
        elif d['pop'] == 'S2g': P.append(mk(d, 1.0, 'S2g'))
        elif d['pop'] == 'S2x': P.append(mk(d, 1.0, 'S2x'))
        elif d['pop'] == 'Z0': P.append(mk(d, 1.0, 'Z0'))
    for e in man['H3_control_fneg']: P.append(mk(S[(e['cs'], e['idx'])], e['f'], 'H3'))
    for e in man['H3_boundary_f0']: P.append(mk(S[(e['cs'], e['idx'])], e['f'], 'H3b'))
    for a in man['anchors']: P.append(dict(cs=a['cs'], idx=None, gphi=a['gphi'], gpsi=a['gpsi'], w=a['w'], q=a['q'], v=a['v'], r=a['r'], f=a['f'], pop='ANCHOR', label=a['label']))
    return P

def _wrap(p):
    return p, evaluate(p)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--manifest'); ap.add_argument('--out'); ap.add_argument('--pop', default=None); ap.add_argument('--limit', type=int, default=0); ap.add_argument('--procs', type=int, default=4); a = ap.parse_args()
    man = json.load(open(a.manifest)); P = build(man)
    if a.pop: P = [p for p in P if p['pop'] in a.pop.split(',')]
    if a.limit: P = P[:a.limit]
    t0 = time.time()
    def cid(p): return json.dumps([p['pop'], p['cs'], p.get('idx'), p['f'], p.get('label')])
    jl = a.out + 'l'
    done = {}
    if os.path.exists(jl):
        for line in open(jl):
            try: r = json.loads(line); done[r['_cid']] = r
            except Exception: pass
    todo = [p for p in P if cid(p) not in done]
    print('N', len(P), 'already', len(done), 'todo', len(todo), flush=True)
    with Pool(a.procs) as pool, open(jl, 'a') as fh:
        for p, r in pool.imap_unordered(_wrap, todo, chunksize=1):
            r = dict(r); r['_cid'] = cid(p); done[r['_cid']] = r
            fh.write(json.dumps(r, default=lambda o: str(o)) + '\n'); fh.flush()
            if len(done) % 25 == 0: print('progress', len(done), len(P), round(time.time()-t0), flush=True)
    R = [done[cid(p)] for p in P]
    json.dump(R, open(a.out, 'w'), indent=1, default=lambda o: str(o)); print('done', time.time()-t0, flush=True)
