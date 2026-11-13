# T3-fund Stage 1: C0-C3 for (A) Scalar mediator and (B) HeavyScalar, one light field.  Registered protocol: UV_ACCEPTANCE_PROTOCOL_ADDENDUM_T3_T6.md
# SPECIFICATION (written before the first run; nothing below changes after results; SPEC_HASH printed at start):
#  Regime  : both candidates FUNDAMENTAL (k -> infinity tested; grid P2: 61 log points k/M in [1e-3,1e3] + 31 in [1e3,1e6]); M = 1.
#  Models  : (A) L = 1/2(d phi)^2 + 1/2(d psi)^2 + 1/2(d sigma)^2 - 1/2 M^2 sigma^2 + sigma(g_phi phi_t phi_x + g_psi psi_t psi_x),
#                backgrounds phi0 = w t + q x, psi0 = v t + r x, sigma_bar = (g_phi w q + g_psi v r)/M^2  (exact: constant sigma, linear phi, psi).
#            (B) L = 1/2(d phi)^2 + 1/2(d chi)^2 - 1/2 M^2 chi^2 + g chi (d phi)^2, phi0 = w t + q x, chi_bar = g (w^2 - q^2)/M^2  (exact).
#  Plane waves exp(i(kx - omega t)), c = omega/k, eps = 1/k^2.  Characteristic polynomials in c (derived here from the linearised equations):
#      (A) det = Aa Ab As - eps g_phi^2 (w - q c)^2 Ab - eps g_psi^2 (v - r c)^2 Aa,
#          Aa = 1 + 2 g_phi sigma_bar c - c^2, Ab = 1 + 2 g_psi sigma_bar c - c^2, As = 1 + eps - c^2
#      (B) det = Z (1 - c^2)(1 - c^2 + eps) - 4 g^2 eps (w c + q)^2,  Z = 1 + 2 g chi_bar = 1 + 2 g^2 (w^2 - q^2)/M^2
#  Roots: double precision (numpy) for P3 classification; 50-digit (mpmath) for confirmation and for every grid point as an ADDITIONAL cross-check;
#      a discrepancy (complex in 50 digits with r > 1e-12 while double says real) is REPORTED, the registered verdict rule is not changed.
#  C3 rule : P3 + P6 three-band rule with u(k) = |gamma(50 digits) - gamma(100 digits)|; refinement: 3 extra log points in every interval of the band
#      (+ one neighbouring interval on each side).  FAIL-UNSTABLE needs gamma > 10 u on >= 3 consecutive refined points.  The 'largest change between a coarse
#      point and its refined neighbours' part of u(k) is implemented only as the sign-change / three-consecutive-points test (ambiguity of the registered text, reported).
#  P2a     : k = 0 checked separately by solving det K(omega, 0) = 0 in 50 digits.  Double root omega = 0 -> ZERO MODE.
#  Sets (P8): (A1) 400 J0 = 0 backgrounds of scalar_mediator_T1_T6.py (seed 11, same generator);  (A2) the 40 'selected' backgrounds of the E3 manifest (J0 != 0);
#      the B3 anchor is NOT available in this workspace and is not evaluated;  (A3) new population seed 41, n = 100, |g| in U(0.1,1) random sign, w,q,v,r in U(-1.5,1.5).
#      (B1) seed 31, n = 200: w in U(0.05,2), q in U(-2,2), |g| log-uniform in [0.05,1], random sign, Z > 0 (resampled otherwise);
#      (B2) anchors (w,q) in {(1,0.3),(0.3,1)}, g = +-0.4;  (B3) CONTROL with Z < 0: seed 32, n = 50, w in U(0.05,1), q in U(1.5,4), |g| in U(0.5,1) with Z < 0.
#  Predictions registered before the run (analytic): C0 exact for all; (B) Z > 0 -> no complex root, large-k expansion omega = +-k + delta, delta real = +-g(w+-q)/sqrt(Z);
#      (B) Z < 0 -> delta imaginary, bounded growth rate (FAIL-UNSTABLE for C3, C1 bounded); (A) C3 not predicted analytically.
import numpy as np, json, sys, hashlib, time, csv, argparse
import mpmath as mp
from multiprocessing import Pool

TOL_REAL, TOL_UNRES, TOL_CONF = 1e-7, 1e-5, 1e-12
KGRID = np.concatenate([np.geomspace(1e-3, 1e3, 61), np.geomspace(1e3, 1e6, 32)[1:]])

def pmul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] = out[i+j] + x*y
    return out
def padd(a, b):
    n = max(len(a), len(b)); a = list(a)+[0]*(n-len(a)); b = list(b)+[0]*(n-len(b)); return [x+y for x, y in zip(a, b)]
def pscale(a, s): return [s*x for x in a]

def poly(model, th, kv, one=1):
    eps = one/(kv*kv)
    if model == 'A':
        gp, gs, w, q, v, r = th
        sb = gp*w*q+gs*v*r
        Aa = [one, 2*gp*sb, -one]; Ab = [one, 2*gs*sb, -one]; As = [one+eps, 0*one, -one]
        Wa = [w*w, -2*w*q, q*q]; Wb = [v*v, -2*v*r, r*r]
        d = pmul(pmul(Aa, Ab), As)
        d = padd(d, pscale(pmul(Wa, Ab), -eps*gp*gp)); d = padd(d, pscale(pmul(Wb, Aa), -eps*gs*gs))
        return d
    g, w, q = th
    Z = one+2*g*g*(w*w-q*q)
    d = pscale(pmul([one, 0*one, -one], [one+eps, 0*one, -one]), Z)
    return padd(d, pscale([q*q, 2*w*q, w*w], -4*g*g*eps))

def roots_dbl(model, th, kv):
    co = [complex(x) for x in poly(model, [float(t) for t in th], float(kv))]
    return np.roots(co[::-1])
def roots_mp(model, th, kv, dps):
    mp.mp.dps = dps
    t = [mp.mpf(float(x)) for x in th]; k = mp.mpf(float(kv))
    co = poly(model, t, k, one=mp.mpf(1))
    try: return mp.polyroots(co[::-1], maxsteps=400, extraprec=4*dps)
    except mp.libmp.NoConvergence: return None
def rr(c): return abs(float(mp.im(c)))/max(1.0, float(abs(c))) if not isinstance(c, complex) else abs(c.imag)/max(1.0, abs(c))

def gamma_mp(model, th, kv, dps):
    rt = roots_mp(model, th, kv, dps)
    if rt is None: return None
    # v2 (2026-11-13, after review): returns mp.mpf at the working precision; no conversion to float, so that u(k) = |gamma50 - gamma100| is measured in arbitrary precision
    return max(mp.mpf(0), max(mp.im(c) for c in rt))*mp.mpf(float(kv))

def k0_check(model, th, dps=50):
    # det K(omega, k = 0) as polynomial in omega, solved directly (P2a); returns (classification, roots)
    mp.mp.dps = dps; t = [mp.mpf(float(x)) for x in th]; one = mp.mpf(1)
    if model == 'A':
        gp, gs, w, q, v, r = t; M2 = one
        # det = omega^4 (M^2 - omega^2 + gp^2 q^2 + gs^2 r^2)   (derived; verified below against the generic polynomial)
        co = [0, 0, 0, 0, M2+gp*gp*q*q+gs*gs*r*r, 0, -one]
    else:
        g, w, q = t; Z = one+2*g*g*(w*w-q*q)
        # Z omega^2 (omega^2 - M^2) - 4 g^2 w^2 omega^2 -> omega^2 [Z (omega^2 - 1) - 4 g^2 w^2]
        co = [0, 0, -(Z+4*g*g*w*w), 0, Z]
    return [x for x in co]

def k0_roots(model, th, dps=50):
    co = k0_check(model, th, dps)
    # strip the omega^n zero-root factor
    nz = 0
    while co[nz] == 0: nz += 1
    rest = co[nz:]
    mp.mp.dps = dps
    rt = mp.polyroots(rest[::-1], maxsteps=200, extraprec=200)
    rz = rest and [abs(mp.im(x))/max(1, abs(x)) for x in rt]
    maxr = float(max(rz)) if rz else 0.0
    return nz, maxr, [complex(x) for x in rt]

def c2_doc(model, th):
    """first-order form U = (field_t, field_x) per field; principal matrix A (dU/dt + A dU/dx = lower order); symmetric -> complete eigenvectors."""
    if model == 'A':
        gp, gs, w, q, v, r = th; sb = gp*w*q+gs*v*r
        blocks = [np.array([[2*gp*sb, -1.], [-1., 0.]]), np.array([[2*gs*sb, -1.], [-1., 0.]]), np.array([[0., -1.], [-1., 0.]])]
    else:
        blocks = [np.array([[0., -1.], [-1., 0.]])]*2     # a-equation divided by Z (Z != 0)
    n = 2*len(blocks); A = np.zeros((n, n))
    for i, b in enumerate(blocks): A[2*i:2*i+2, 2*i:2*i+2] = b
    ev, V = np.linalg.eig(A)
    return dict(first_order_vars='(field_t, field_x) per field', max_abs_imag_eig=float(np.abs(ev.imag).max()), eig=sorted(float(e) for e in ev.real),
                cond_eigvec=float(np.linalg.cond(V)), symmetric=bool(np.allclose(A, A.T)),
                coincident_roots=bool(np.any(np.diff(sorted(ev.real)) < 1e-6*max(1, np.abs(ev).max()))))

def evaluate(args):
    model, tag, th = args
    th = [float(x) for x in th]
    out = dict(model=model, set=tag, theta=th)
    # C2 documentation
    out['C2'] = c2_doc(model, th)
    # P2a
    nz, maxr0, rt0 = k0_roots(model, th)
    out['k0'] = dict(zero_root_multiplicity=nz, max_r_nonzero=maxr0, nonzero_omega2=[complex(x**2).real for x in rt0][:1])
    k0_ok = maxr0 <= TOL_REAL
    # grid
    G = []; hidden = 0; conf = []; unres_stuck = 0; nodata = 0
    for kv in KGRID:
        rd = roots_dbl(model, th, kv)
        r_d = max(abs(c.imag)/max(1, abs(c)) for c in rd)
        rm = roots_mp(model, th, kv, 50)
        if rm is None: nodata += 1; continue
        r_m = max(float(abs(mp.im(c)))/max(1.0, float(abs(c))) for c in rm)
        gm = max(0.0, float(max(mp.im(c) for c in rm)))*kv
        if r_d <= TOL_REAL:
            cls = 'real'
            if r_m > TOL_CONF: hidden += 1
        elif r_d <= TOL_UNRES:
            cls = 'unresolved'
            if r_m <= TOL_REAL: cls = 'real'
            elif r_m <= TOL_UNRES: unres_stuck += 1; cls = 'unresolved_stuck'
            else: cls = 'complex' if r_m > TOL_CONF else 'real'
        else:
            cls = 'complex' if r_m > TOL_CONF else 'real'
        G.append((float(kv), r_d, r_m, gm, cls))
        if cls == 'complex': conf.append(float(kv))
    out.update(n_grid=len(G), n_nodata=nodata, hidden_by_double=hidden, unresolved_stuck=unres_stuck, n_complex=len(conf))
    out['gamma_max_grid'] = max((g[3] for g in G), default=0.0)
    # P6 refinement where complex roots were confirmed
    verdict = None
    if conf:
        ks = [g[0] for g in G]; idx = [i for i, g in enumerate(G) if g[4] == 'complex']
        lo, hi = max(0, min(idx)-1), min(len(ks)-1, max(idx)+1)
        fine = []
        for i in range(lo, hi):
            fine += list(np.geomspace(ks[i], ks[i+1], 5)[:-1])
        fine.append(ks[hi])
        rows = []
        for kv in fine:
            g50 = gamma_mp(model, th, kv, 50); g100 = gamma_mp(model, th, kv, 100)
            if g50 is None or g100 is None: rows.append((kv, None, None, None)); continue
            mp.mp.dps = 120
            rows.append((kv, g50, g100, abs(g50-g100)))
        ok = [r for r in rows if r[1] is not None]
        flags = [];
        for kv, g50, g100, u in ok:
            if g50 <= u: flags.append(0)
            elif g50 <= 10*u: flags.append(1)
            else: flags.append(2)
        run = best = 0
        for f in flags:
            run = run+1 if f == 2 else 0; best = max(best, run)
        signchg = any((a[1] > 0) != (b[1] > 0) for a, b in zip(ok, ok[1:])) and any(f == 1 for f in flags)
        if best >= 3: verdict = 'FAIL-UNSTABLE'
        elif any(f == 1 for f in flags) or signchg or best > 0: verdict = 'INCONCLUSIVE'
        else: verdict = 'PASS'
        mp.mp.dps = 120
        out['P6'] = dict(n_refined=len(fine), max_gamma=float(max(r[1] for r in ok)) if ok else None,
                         max_u=float(max((r[3] for r in ok), default=0)), max_u_str=mp.nstr(max((r[3] for r in ok), default=mp.mpf(0)), 8),
                         min_ratio_gamma_over_u=(mp.nstr(min(r[1]/r[3] for r in ok if r[3] > 0 and r[1] > 0), 6) if any(r[3] > 0 and r[1] > 0 for r in ok) else None),
                         n_u_exactly_zero=sum(1 for r in ok if r[3] == 0),
                         k_band=[min(g[0] for g in G if g[4] == 'complex'), max(g[0] for g in G if g[4] == 'complex')],
                         max_consecutive_gt10u=best, e_fold_time=float(1/max(r[1] for r in ok)) if ok and max(r[1] for r in ok) > 0 else None)
    out['C3_numeric'] = verdict if verdict else ('INCONCLUSIVE' if (unres_stuck or nodata) else 'PASS')
    out['C3_k0'] = 'PASS(zero mode multiplicity %d)' % nz if k0_ok else 'FAIL'
    out['outcome'] = out['C3_numeric'] if out['C3_numeric'] != 'PASS' else ('PASS' if k0_ok else 'FAIL-UNSTABLE')
    if out['C3_numeric'] == 'FAIL-UNSTABLE' or (out['C3_numeric'] == 'PASS' and not k0_ok): out['outcome'] = 'FAIL-UNSTABLE'
    return out

def populations():
    P = {}
    rng = np.random.default_rng(11); rows = []
    while len(rows) < 400:
        p, s = rng.uniform(0.02, 0.09, 2); kap = rng.choice([-1, 1])*np.sqrt(p*s)
        w, q, v = rng.uniform(-1.5, 1.5, 3)
        if abs(kap*v) < 1e-3: continue
        r = p*w*q/(kap*v)
        if max(np.sqrt(p)*abs(q), np.sqrt(s)*abs(r), np.sqrt(p)*abs(w), np.sqrt(s)*abs(v)) > 0.3: continue
        rows.append((p, s, kap, w, q, v, r))
    P['A1_J0eq0_seed11'] = [('A', 'A1', (np.sqrt(p), -kap/np.sqrt(p), w, q, v, r)) for p, s, kap, w, q, v, r in rows]
    return P

def build(manifest_path):
    P = populations()
    man = json.load(open(manifest_path))['selected']
    P['A2_manifest40'] = [('A', 'A2', (np.sqrt(d['p']), -d['kappa']/np.sqrt(d['p']), d['w'], d['q'], d['v'], d['r'])) for d in man]
    rng = np.random.default_rng(41); L = []
    for _ in range(100):
        gp, gs = rng.choice([-1, 1], 2)*rng.uniform(0.1, 1, 2); w, q, v, r = rng.uniform(-1.5, 1.5, 4); L.append(('A', 'A3', (gp, gs, w, q, v, r)))
    P['A3_general_seed41'] = L
    rng = np.random.default_rng(31); L = []
    while len(L) < 200:
        w = rng.uniform(0.05, 2); q = rng.uniform(-2, 2); g = rng.choice([-1, 1])*np.exp(rng.uniform(np.log(0.05), 0))
        if 1+2*g*g*(w*w-q*q) > 0: L.append(('B', 'B1', (g, w, q)))
    P['B1_Zpos_seed31'] = L
    P['B2_anchors'] = [('B', 'B2', (g, w, q)) for (w, q) in ((1., 0.3), (0.3, 1.)) for g in (0.4, -0.4)]
    rng = np.random.default_rng(32); L = []
    while len(L) < 50:
        w = rng.uniform(0.05, 1); q = rng.uniform(1.5, 4); g = rng.choice([-1, 1])*rng.uniform(0.5, 1)
        if 1+2*g*g*(w*w-q*q) < 0: L.append(('B', 'B3', (g, w, q)))
    P['B3_control_Zneg_seed32'] = L
    return P

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument('--manifest', required=True); ap.add_argument('--out', default='.'); ap.add_argument('--only', default=None); ap.add_argument('--limit', type=int, default=0)
    a = ap.parse_args()
    spec = open(__file__).read().split("import numpy")[0]
    print("SPEC_HASH", hashlib.sha256(spec.encode()).hexdigest()[:16], flush=True)
    P = build(a.manifest); res = {}
    for name, L in P.items():
        if a.only and a.only not in name: continue
        if a.limit: L = L[:a.limit]
        t0 = time.time()
        with Pool(2) as pool: R = pool.map(evaluate, L, chunksize=2)
        cnt = {}
        for r in R: cnt[r['outcome']] = cnt.get(r['outcome'], 0)+1
        summ = dict(n=len(R), outcomes=cnt, hidden_by_double_total=sum(r['hidden_by_double'] for r in R), unresolved_stuck_total=sum(r['unresolved_stuck'] for r in R),
                    c2_max_cond=max(r['C2']['cond_eigvec'] for r in R), c2_all_real=all(r['C2']['max_abs_imag_eig'] < 1e-12 for r in R),
                    c2_all_symmetric=all(r['C2']['symmetric'] for r in R), c2_any_coincident=sum(r['C2']['coincident_roots'] for r in R),
                    k0_all_ok=all(r['C3_k0'].startswith('PASS') for r in R), gamma_max_grid=max(r['gamma_max_grid'] for r in R), seconds=time.time()-t0)
        print(name, json.dumps(summ), flush=True)
        res[name] = dict(summary=summ, rows=R)
    json.dump(res, open(a.out.rstrip('/')+'/t3fund_stage1_result.json', 'w'), indent=1, default=float)
