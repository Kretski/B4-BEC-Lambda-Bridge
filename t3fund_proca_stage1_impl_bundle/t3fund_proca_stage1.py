# T3-fund Proca mediator: C0-C3 + C1 checks. Implements T3_FUND_PROCA_SPEC.md (registered in commit e355af2; SPEC_HASH 86c360326da7406f, MANIFEST_HASH 35a3c113d517f174).
# Implementation notes (not changes of the spec): the six roots at each k/M=x are computed from Q(c;eps) built by exact polynomial arithmetic; the C3 classification loop (double screening at every grid point,
# 50-digit roots at EVERY grid point, classification by the 50-digit value r_m alone, P6 three-band rule with 50/100-digit evaluation of every refinement point) follows the SM+eta Stage 1 loop. M=1 (spectrum depends on x=k/M only).
import os, sys, json, time, argparse, itertools
import numpy as np, mpmath as mp, sympy as sp
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smeta_ref as R          # unchanged SM+eta Stage 1 module: pmul/padd/pscale, Phi/Kentries/det4 (jet-Hessian Leibniz determinant), KGRID, tolerances
TOL_REAL, TOL_UNRES, TOL_CONF = R.TOL_REAL, R.TOL_UNRES, R.TOL_CONF
KGRID = R.KGRID
SCOPE_PN = 100        # first PN draws per cell that get the extended checks (Hessian route, mu check, exact discriminant)

# ---------- parameters ----------
def pars(case, dps=50):
    mp.mp.dps = dps; m = mp.mpf
    P = dict(p=m(case['p']), s=m(case['s']), kap=m(case['kappa']), Lam=m(case['Lambda']), w=m(case['w']), q=m(case['q']), v=m(case['v']), r=m(case['r']), gp=m(case['gphi']), gs=m(case['gpsi']))
    # use the registered prescription (g_phi=sqrt p, g_psi=-kappa/sqrt p) as stored in the manifest; p,s,kappa as stored
    P['Abar'] = P['gp']*P['w']*P['q']+P['gs']*P['v']*P['r']        # M=1
    P['aphi'] = P['p']*P['w']*P['q']-P['kap']*P['v']*P['r']; P['apsi'] = -P['kap']*P['w']*P['q']+P['s']*P['v']*P['r']
    return P

# ---------- Q(c; eps) ----------
def Qcoeffs(P, eps):
    """ascending coefficients in c of Q = F - eps c^2 d1 d2."""
    one = mp.mpf(1)
    d1 = [-(1+3*P['Lam']*P['q']**2), -2*P['aphi'], one]; d2 = [-one, -2*P['apsi'], one]
    wq = [P['w'], -P['q']]; vr = [P['v'], -P['r']]
    t1 = R.pmul(R.pmul(wq, wq), d2); t2 = R.pmul(R.pmul(vr, vr), d1)
    F = R.padd(R.padd(R.pmul(d1, d2), R.pscale(t1, P['p'])), R.pscale(t2, P['s']))
    c2d = R.pmul([0, 0, one], R.pmul(d1, d2))
    return R.padd(F, R.pscale(c2d, -eps))
def eps_of(x): return x*x/(1+x*x)
def roots_Q_dbl(P, x):
    co = [float(mp.re(z)) for z in Qcoeffs(P, eps_of(mp.mpf(x)))]
    return np.roots(co[::-1])
def roots_Q_mp(P, x, dps):
    mp.mp.dps = dps; Pp = P if dps == 50 else None
    co = [mp.re(z) for z in Qcoeffs(P, eps_of(mp.mpf(x)))]
    try: return mp.polyroots(co[::-1], maxsteps=500, extraprec=6*dps)
    except mp.libmp.NoConvergence: return None

# ---------- 6x6 reduced matrix (spec 3(iii)) ----------
def lin_parts(P):
    I = mp.mpc(0, 1); gp, gs, Ab = P['gp'], P['gs'], P['Abar']; Lam, w, q, v, r = P['Lam'], P['w'], P['q'], P['v'], P['r']
    Ps = mp.zeros(6, 6); Pk = mp.zeros(6, 6)
    Ps[0, 1] = I*(1+3*Lam*q**2); Ps[0, 0] = -2*I*gp*Ab; Ps[0, 4] = -gp*I*w
    Ps[1, 0] = I; Ps[2, 3] = I; Ps[2, 2] = -2*I*gs*Ab; Ps[2, 4] = -gs*I*v; Ps[3, 2] = I
    Pk[0, 5] = -gp*q; Pk[2, 5] = -gs*r; Pk[4, 5] = 1
    Pk[5, 0] = gp*q; Pk[5, 1] = gp*w; Pk[5, 2] = gs*r; Pk[5, 3] = gs*v; Pk[5, 4] = -1
    return Ps, Pk
def Mlin(Ps, Pk, k):
    kk = mp.mpf(k); return kk*Ps+mp.sqrt(1+kk*kk)*Pk

# ---------- jet Hessian of the Lagrangian (symbolic once, lambdified to mpmath) ----------
_HF = None
def hessian_fun():
    global _HF
    if _HF is not None: return _HF
    names = [f'{fld}{j}' for fld in ('f', 'g', 'a0', 'a1') for j in ('v', 't', 'x')]
    S = {n: sp.Symbol(n, real=True) for n in names}; gp, gs, Lam = sp.symbols('gp gs Lam', real=True)
    half = sp.Rational(1, 2); M = 1
    ft, fx, gt, gx = S['ft'], S['fx'], S['gt'], S['gx']; a0, a0x, a1, a1t = S['a0v'], S['a0x'], S['a1v'], S['a1t']
    L = half*ft**2-half*fx**2-Lam/4*fx**4+half*gt**2-half*gx**2+half*(a1t-a0x)**2+half*M**2*(a0**2-a1**2)+a1*(gp*ft*fx+gs*gt*gx)
    syms = [S[n] for n in names]; H = sp.hessian(L, syms)
    fn = sp.lambdify([gp, gs, Lam]+syms, H, modules='mpmath'); _HF = fn; return fn
def hessian_mp(P):
    fn = hessian_fun(); mp.mp.dps = 50
    # order of jets: (fv,ft,fx),(gv,gt,gx),(a0v,a0t,a0x),(a1v,a1t,a1x); background values
    vals = [0, P['w'], P['q'], 0, P['v'], P['r'], 0, 0, 0, P['Abar'], 0, 0]
    return fn(P['gp'], P['gs'], P['Lam'], *vals)
def detK_poly(P, x):
    """det K(omega=c k,k)/k^8 as polynomial in c via the jet Hessian (Leibniz), eps_scale=1/k (M=1)."""
    H = hessian_mp(P); Hl = [[H[i, j] for j in range(12)] for i in range(12)]
    mp.mp.dps = 50; k = mp.mpf(x)
    d = R.det4(R.Kentries(Hl, mp.mpf(1)/k, mp.mpc(0, 1))); return [mp.re(z) for z in d], [mp.im(z) for z in d]

# ---------- C0 ----------
_GF = None
def c0_residuals(P):
    """Euler-Lagrange residuals on the constant-jet background.  On constant jets every derivative term of the EL equations vanishes, so the EL residual of a field is dL/d(value jet) evaluated
    at the background jets; the gradient of L is taken symbolically from the Lagrangian (independent of the definition of Abar) and lambdified.  'abar_construction' is the algebraic identity
    Abar = (gphi w q + gpsi v r)/M^2 used to build the background (it is the same expression used to define Abar, i.e. it is NOT an independent check)."""
    global _GF
    mp.mp.dps = 50
    if _GF is None:
        names = [f'{fld}{j}' for fld in ('f', 'g', 'a0', 'a1') for j in ('v', 't', 'x')]
        Sy = {n: sp.Symbol(n, real=True) for n in names}; gp, gs, Lam = sp.symbols('gp gs Lam', real=True); half = sp.Rational(1, 2)
        ft, fx, gt, gx = Sy['ft'], Sy['fx'], Sy['gt'], Sy['gx']; a0, a0x, a1, a1t = Sy['a0v'], Sy['a0x'], Sy['a1v'], Sy['a1t']
        L = half*ft**2-half*fx**2-Lam/4*fx**4+half*gt**2-half*gx**2+half*(a1t-a0x)**2+half*(a0**2-a1**2)+a1*(gp*ft*fx+gs*gt*gx)
        syms = [Sy[n] for n in names]; grad = [sp.diff(L, z) for z in syms]; _GF = sp.lambdify([gp, gs, Lam]+syms, grad, modules='mpmath')
    vals = [0, P['w'], P['q'], 0, P['v'], P['r'], 0, 0, 0, P['Abar'], 0, 0]
    g_ = _GF(P['gp'], P['gs'], P['Lam'], *vals)
    # EL_X = dL/dX - d_t(dL/dX_t) - d_x(dL/dX_x); on constant jets the derivative terms vanish -> residual = dL/dX_value (indices 0,3,6,9 = f,g,a0,a1 values)
    res = dict(EL_f=float(abs(g_[0])), EL_g=float(abs(g_[3])), EL_a0=float(abs(g_[6])), EL_a1=float(abs(g_[9])))
    res['abar_construction'] = float(abs(P['Abar']-(P['gp']*P['w']*P['q']+P['gs']*P['v']*P['r'])))
    return res

# ---------- exact discriminant criterion (Pred-4) ----------
def disc_criterion(case):
    """E = six real roots at eps=1e-4 and disc_c Q(.;eps) has no real root in (0,1]; exact rational arithmetic on the binary values of the manifest numbers."""
    R_ = sp.Rational; c, e = sp.symbols('c e')
    p, s, kap, Lam, w, q, v, r = [R_(case[k]) for k in ('p', 's', 'kappa', 'Lambda', 'w', 'q', 'v', 'r')]
    gp, gs = [R_(case[k]) for k in ('gphi', 'gpsi')]
    aphi = p*w*q-kap*v*r; apsi = -kap*w*q+s*v*r
    d1 = c**2-2*aphi*c-(1+3*Lam*q**2); d2 = c**2-2*apsi*c-1
    F = sp.expand(d1*d2+p*(w-q*c)**2*d2+s*(v-r*c)**2*d1); Qe = sp.expand(F-e*c**2*d1*d2)
    D = sp.Poly(sp.expand(sp.discriminant(sp.Poly(Qe, c))), e)
    nroots = D.count_roots(sp.Rational(0), sp.Rational(1)) if D.degree() > 0 else 0
    # count_roots counts roots in the closed interval [0,1] (multiplicity ignored); exclude a root at 0 (disc vanishes at e=0 when F has repeated roots)
    r0 = 1 if D.eval(0) == 0 else 0
    nreal_small = sp.Poly(Qe.subs(e, R_(1, 10**4)), c).count_roots()
    return dict(n_disc_roots_closed_interval=int(nroots), disc_zero_at_0=bool(r0), n_real_roots_small_eps=int(nreal_small), degree=int(D.degree()),
                E=bool(nreal_small == 6 and (nroots-r0) == 0))

# ---------- k = 0 ----------
def k0_info(P):
    Ps, Pk = lin_parts(P); mp.mp.dps = 50; M0 = Mlin(Ps, Pk, 0)
    ev = mp.eig(M0, left=False, right=False); om = [mp.mpc(0, 1)*e for e in ev]      # omega = i*lambda
    nz = [z for z in om if abs(z) > mp.mpf('1e-10')]
    maxr = float(max((abs(mp.im(z))/max(1, abs(z)) for z in nz), default=0))
    closed = 1+P['p']*P['q']**2+P['s']*P['r']**2          # mp, full precision
    devs = [abs(mp.re(z)**2-closed)/max(1, closed) for z in nz]; dev = max(devs) if devs else mp.inf
    # no rounding before the comparison; only the reported numbers below are converted to float
    return dict(n_zero=len(om)-len(nz), n_nonzero=len(nz), max_r_nonzero=maxr, omega2=sorted(float(mp.re(z)**2) for z in nz), omega2_closed=float(closed), max_rel_dev=float(dev),
                agree=bool(len(nz) == 2 and dev <= mp.mpf('1e-9')))

# ---------- evaluation ----------
def evaluate(case):
    out = dict(case=dict((k, case.get(k)) for k in ('pop', 'label', 'cut', 'mode')))
    scope = bool(case['pop'] == 'PE' or case.get('scope'))
    P = pars(case, 50)
    out['C0'] = c0_residuals(P)
    kap_kin = [1.0, case.get('kinetic_sign_psi', 1.0), 1.0]; out['kappa_kin'] = kap_kin; TAU = R.TAU
    if min(kap_kin) < -TAU: out['outcome'] = 'FAIL-GHOST'; return out
    if min(abs(z) for z in kap_kin) <= TAU: out['outcome'] = 'DEGENERATE'; return out
    # UV principal symbol: Q(c;1)
    rU = mp.polyroots([mp.re(z) for z in Qcoeffs(P, mp.mpf(1))][::-1], maxsteps=500, extraprec=300)
    maxim = float(max(abs(mp.im(z)) for z in rU)); cs = sorted(float(mp.re(z)) for z in rU)
    sepU = min(abs(rU[i]-rU[j]) for i in range(6) for j in range(i))
    out['uv'] = dict(max_imag=maxim, speeds=cs, min_sep=float(sepU), front=float(max(abs(z) for z in rU)), strictly_hyperbolic=bool(maxim <= 1e-9 and sepU > 1e-9))
    # k=0
    k0 = k0_info(P); out['k0'] = k0; k0_ok = k0['max_r_nonzero'] <= R.TOL_REAL
    # C1 consistency (b): eigenvalues of M(k) vs k*roots(Q) at x in {0.7,3,50}, 50 digits
    Ps, Pk = lin_parts(P); dev = 0.0
    for x in (0.7, 3.0, 50.0):
        ev = mp.eig(Mlin(Ps, Pk, x), left=False, right=False); om = sorted([mp.mpc(0, 1)*e for e in ev], key=lambda z: (float(mp.re(z)), float(mp.im(z))))
        rt = sorted([z*mp.mpf(x) for z in roots_Q_mp(P, x, 50)], key=lambda z: (float(mp.re(z)), float(mp.im(z))))
        dev = max(dev, float(max(abs(a-b)/max(1, abs(b)) for a, b in zip(om, rt))))
    out['C1b'] = dict(max_rel_dev_M_vs_Q=dev, ok=bool(dev <= 1e-9))
    # C3 grid
    G = []; hidden = 0; hidden_small = 0; conf = []; stuck = 0; nodata = 0; speeds_uv = None
    for ik, kv in enumerate(KGRID):
        rd = roots_Q_dbl(P, kv); r_d = float(max(abs(z.imag)/max(1, abs(z)) for z in rd))
        if ik == len(KGRID)-1: speeds_uv = sorted(float(z.real) for z in rd) if r_d <= TOL_REAL else None
        rm = roots_Q_mp(P, kv, 50)            # 50-digit roots at EVERY grid point (superset of the registered confirmation points)
        if rm is None: nodata += 1; continue
        r_m = max(float(abs(mp.im(z)))/max(1.0, float(abs(z))) for z in rm)
        # unified classification by the 50-digit value r_m (independent of the double-precision screening branch): r_m<=1e-7 real; 1e-7<r_m<=1e-5 unresolved; r_m>1e-5 complex.
        # TOL_CONF (1e-12) is not used for the classification; it only feeds the diagnostic counter real_with_tiny_imag.
        if r_m <= TOL_REAL: cls = 'real'
        elif r_m <= TOL_UNRES: stuck += 1; cls = 'unresolved_stuck'
        else: cls = 'complex'
        if r_d <= TOL_REAL and r_m > TOL_REAL: hidden += 1               # double precision said real, 50 digits say not (never left as 'real')
        if cls == 'real' and r_m > R.TOL_CONF: hidden_small += 1         # 1e-12 < r_m <= 1e-7: real by the 1e-7 tolerance, recorded separately
        G.append((float(kv), r_d, r_m, cls))
        if cls == 'complex': conf.append(float(kv))
    out['n_real'] = sum(1 for g_ in G if g_[3] == 'real'); out['n_unresolved'] = sum(1 for g_ in G if g_[3] == 'unresolved_stuck')
    # Pred-3
    if speeds_uv is not None: out['uv']['pred3_max_err'] = float(max(abs(a-b) for a, b in zip(speeds_uv, cs)))
    out.update(n_grid=len(G), n_nodata=nodata, hidden_by_double=hidden, real_with_tiny_imag=hidden_small, unresolved_stuck=stuck, n_complex=len(conf), complex_at_kmin=bool(G and G[0][3] == 'complex'), complex_in_uv_band=bool(any(k >= 1e3 for k in conf)))
    verdict = None
    if conf:
        ks = [g_[0] for g_ in G]; idx = [i for i, g_ in enumerate(G) if g_[3] == 'complex']
        lo, hi = max(0, min(idx)-1), min(len(ks)-1, max(idx)+1); fine = []
        for i in range(lo, hi): fine += list(np.geomspace(ks[i], ks[i+1], 5)[:-1])
        fine.append(ks[hi]); rows = []
        def gam(kv, dps):
            rt = roots_Q_mp(P, kv, dps)
            return None if rt is None else max(mp.mpf(0), max(mp.im(z) for z in rt))*mp.mpf(float(kv))
        fidx = []                                              # indices (into `fine`) of the refinement points evaluated at BOTH 50 and 100 digits
        for ifine, kv in enumerate(fine):
            g50 = gam(kv, 50); g100 = gam(kv, 100)
            if g50 is None or g100 is None: continue
            mp.mp.dps = 120; rows.append((kv, g50, g100, abs(g50-g100))); fidx.append(ifine)
        flags = [0 if r_[1] <= r_[3] else (1 if r_[1] <= 10*r_[3] else 2) for r_ in rows]; run = best = 0
        for fl in flags:
            run = run+1 if fl == 2 else 0; best = max(best, run)
        signchg = any((a[1] > 0) != (b[1] > 0) for a, b in zip(rows, rows[1:])) and any(fl == 1 for fl in flags)
        if best >= 3: verdict = 'FAIL-UNSTABLE'
        elif any(fl == 1 for fl in flags) or signchg or best > 0: verdict = 'INCONCLUSIVE'
        else: verdict = 'PASS'
        complete = bool(len(set(fine)) == len(fine) and fidx == list(range(len(fine))) and len(rows) == len(fine))   # every refinement point present exactly once, in order, no duplicates, all evaluated at 50 and 100 digits
        if not complete: verdict = 'INCONCLUSIVE'          # a refinement point that could not be evaluated at 50/100 digits is never PASS (nor an unqualified FAIL)
        mp.mp.dps = 120
        out['P6'] = dict(n_refined=len(fine), n_unique_refined=len(set(fine)), n_rows=len(rows), complete=complete, max_u_str=mp.nstr(max((r_[3] for r_ in rows), default=mp.mpf(0)), 8), max_consecutive_gt10u=best, max_gamma=float(max((r_[1] for r_ in rows), default=0)) if rows else None)
    c3 = verdict if verdict else 'PASS'
    if c3 == 'PASS' and (stuck or nodata): c3 = 'INCONCLUSIVE'            # an unresolved or unevaluated grid point never allows a positive conclusion (also after a P6 PASS)
    out['C3_numeric'] = c3; out['outcome'] = 'FAIL-UNSTABLE' if (c3 == 'FAIL-UNSTABLE' or (c3 == 'PASS' and not k0_ok)) else c3
    # Pred-6 data
    Fr = np.roots([float(mp.re(z)) for z in Qcoeffs(P, mp.mpf(0))][::-1]); out['pred6'] = dict(F_max_imag=float(np.abs(Fr.imag).max()), first_point_complex=bool(G and G[0][3] == 'complex'))
    # extended checks
    if scope:
        mp.mp.dps = 50
        # Hessian route: ratio of det K/k^8 poly to Q poly must be -(1+x^2)/x^4 ; roots to 1e-9
        hd = 0.0; rr = 0.0
        for x in (0.7, 3.0, 50.0):
            dre, dim_ = detK_poly(P, x); qco = [mp.re(z) for z in Qcoeffs(P, eps_of(mp.mpf(x)))]
            exp_ratio = -(1+mp.mpf(x)**2)/mp.mpf(x)**4
            rr = max(rr, float(max(abs(dre[i]/qco[i]-exp_ratio)/abs(exp_ratio) for i in range(len(qco)) if abs(qco[i]) > mp.mpf('1e-30'))))
            hd = max(hd, float(max(abs(z) for z in dim_))/float(max(abs(z) for z in dre)))
        out['hessian_route'] = dict(max_rel_ratio_dev=rr, max_imag_rel=hd, ok=bool(rr <= 1e-9))
        # mu check at all 92 grid points
        out['mu'] = mu_check(P, out['uv']['strictly_hyperbolic'])
        out['pred4'] = disc_criterion(case)
    return out

def mu_check(P, hyper):
    if not hyper: return dict(defined=False)
    mp.mp.dps = 50; Ps, Pk = lin_parts(P); Pinf = Ps+Pk
    lam, V = mp.eig(Pinf); Vi = mp.inverse(V); viol = 0; maxgam = mp.mpf(0)
    for k in KGRID:
        kk = mp.mpf(float(k)); kh = mp.sqrt(1+kk*kk); ev = mp.eig(kk*Ps+kh*Pk, left=False, right=False); gam = max(mp.re(e) for e in ev)
        Rr = Vi*((kh-kk)*Pk)*V; H = (Rr+Rr.H)/2; mu = max(mp.eigh(H, eigvals_only=True))
        if gam > mu*(1+mp.mpf('1e-9'))+mp.mpf('1e-9'): viol += 1
        maxgam = max(maxgam, gam)
    cond = mp.mnorm(V, 1)*mp.mnorm(Vi, 1)
    return dict(defined=True, violations=int(viol), max_gamma=float(maxgam), cond_V_inf=float(cond))

def build(man):
    P = []
    for e in man['PE']: P.append(dict(e, scope=True))
    cnt = {}
    for e in man['PN']:
        key = (e['cut'], e['mode']); cnt[key] = cnt.get(key, 0)+1; P.append(dict(e, scope=bool(cnt[key] <= SCOPE_PN)))
    for e in man['H3']: P.append(dict(e, scope=False))
    return P
def _wrap(p): return p, evaluate(p)
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--manifest'); ap.add_argument('--out'); ap.add_argument('--pop', default=None); ap.add_argument('--limit', type=int, default=0); ap.add_argument('--procs', type=int, default=2); a = ap.parse_args()
    man = json.load(open(a.manifest)); P = build(man)
    if a.pop: P = [p for p in P if p['pop'] in a.pop.split(',')]
    if a.limit: P = P[:a.limit]
    cid = lambda p: p['label']
    jl = a.out + 'l'; done = {}
    if os.path.exists(jl):
        for line in open(jl):
            try: r = json.loads(line); done[r['_cid']] = r
            except Exception: pass
    todo = [p for p in P if cid(p) not in done]; t0 = time.time()
    print('N', len(P), 'already', len(done), 'todo', len(todo), flush=True)
    with Pool(a.procs) as pool, open(jl, 'a') as fh:
        for p, r in pool.imap_unordered(_wrap, todo, chunksize=1):
            r = dict(r); r['_cid'] = cid(p); done[r['_cid']] = r
            fh.write(json.dumps(r, default=lambda o: str(o))+'\n'); fh.flush()
            if len(done) % 25 == 0: print('progress', len(done), len(P), round(time.time()-t0), flush=True)
    if all(cid(p) in done for p in P):
        json.dump([done[cid(p)] for p in P], open(a.out, 'w'), indent=1, default=lambda o: str(o)); print('done', time.time()-t0, flush=True)
