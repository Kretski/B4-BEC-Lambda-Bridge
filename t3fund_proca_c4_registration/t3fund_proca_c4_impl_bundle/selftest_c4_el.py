"""Check of the ACTUAL C4 right-hand side (t3fund_proca_c4.rhs, pseudo-spectral) against the independent SymPy Euler-Lagrange operator.
Random smooth time-dependent fields f,g,A1 (trigonometric polynomials with integer-wavenumber modes on a periodic box, plus constant-gradient backgrounds f=w t+q x + ..., A1 = const + ...),
A0 from the Gauss constraint (mode by mode); the state (a,b,c,d,A1,F) is built from them; rhs(state) is compared with
  da = f_tt + EL_f,  db = f_xt,  dc = g_tt + EL_g,  dd = g_xt,  dA1 = A1_t,  dF = F_t + EL_A1   (sympy euler_equations returns EL = -(q_tt - RHS)),  and EL_a0 = 0 (Gauss).
Several (M-independent here: M=1 as in the evaluations; one extra case M=1.3 for the sympy/algebra part via rhs's M argument), random couplings and Lambda."""
import numpy as np, sympy as sp, json, sys
sys.path.insert(0, '.')
import t3fund_proca_c4 as C
t, x = sp.symbols('t x', real=True)
f = sp.Function('f')(t, x); g = sp.Function('g')(t, x); a0 = sp.Function('a0')(t, x); a1 = sp.Function('a1')(t, x)
half = sp.Rational(1, 2)
def run_case(seed, M):
    rr = np.random.default_rng(seed)
    gp, gs, Lam = float(rr.uniform(-1, 1)), float(rr.uniform(-1, 1)), float(rr.uniform(0, 1))
    w, q, v, r = [float(z) for z in rr.uniform(-1.2, 1.2, 4)]
    L = half*f.diff(t)**2 - half*f.diff(x)**2 - sp.Float(Lam)/4*f.diff(x)**4 + half*g.diff(t)**2 - half*g.diff(x)**2 + half*(a1.diff(t) - a0.diff(x))**2 + half*sp.Float(M)**2*(a0**2 - a1**2) + a1*(sp.Float(gp)*f.diff(t)*f.diff(x) + sp.Float(gs)*g.diff(t)*g.diff(x))
    EL = sp.calculus.euler.euler_equations(L, [f, g, a0, a1], [t, x])
    def modes():
        return [(int(rr.integers(1, 4)), float(rr.uniform(-2, 2)), float(rr.uniform(-.5, .5)), float(rr.uniform(0, 6))) for _ in range(3)]
    fm, gm, am = modes(), modes(), modes()
    fx_ = sum(A*sp.cos(m*x/4 - n*t + p) for m, n, A, p in fm) + w*t + q*x
    gx_ = sum(A*sp.cos(m*x/4 - n*t + p) for m, n, A, p in gm) + v*t + r*x
    A1x = sum(A*sp.cos(m*x/4 - n*t + p) for m, n, A, p in am) + 0.3
    A0x = 0
    for m, n, A, p in am:
        th = m*x/4 - n*t + p
        A0x += -sp.diff(A*sp.cos(th), t, x)/(M**2 + (m/4)**2)            # (M^2 - d_xx) A0 = -A1_tx for each mode
    sub = {f: fx_, g: gx_, a0: A0x, a1: A1x}
    ELn = [e.lhs.subs(sub).doit() - e.rhs.subs(sub).doit() for e in EL]
    Fv = A1x.diff(t) - A0x.diff(x)
    exprs = dict(a=fx_.diff(t), b=fx_.diff(x), c=gx_.diff(t), d=gx_.diff(x), A=A1x, F=Fv,
                 da=fx_.diff(t, 2) + ELn[0], db=fx_.diff(x, t), dc=gx_.diff(t, 2) + ELn[1], dd=gx_.diff(x, t), dA=A1x.diff(t), dF=Fv.diff(t) + ELn[3], gauss=ELn[2])
    N = 64; Lb = 8*np.pi; xs = np.arange(N)*Lb/N; kk = 2*np.pi*np.fft.rfftfreq(N, Lb/N); worst = 0.0; per = {}
    for t0 in (0.0, 0.37, 1.1):
        fn = {k: sp.lambdify(x, e.subs(t, t0), 'numpy') for k, e in exprs.items()}
        ev = {k: np.broadcast_to(np.array(fn[k](xs), float), xs.shape) for k in exprs}
        # linear-in-x background parts (q x, r x) do not appear in the state (only a,b,c,d,A1,F)
        S = np.array([ev['a'], ev['b'], ev['c'], ev['d'], ev['A'], ev['F'], 0*xs, 0*xs])
        th = dict(gphi=gp, gpsi=gs, Lambda=Lam, w=w, v=v)
        R = C.rhs(S, th, kk, M)
        for name, idx in (('da', 0), ('db', 1), ('dc', 2), ('dd', 3), ('dA', 4), ('dF', 5)):
            e = np.max(np.abs(R[idx] - ev[name])); per[name] = max(per.get(name, 0), e); worst = max(worst, e)
        per['gauss'] = max(per.get('gauss', 0), float(np.max(np.abs(ev['gauss']))))
    return worst, per, dict(gphi=gp, gpsi=gs, Lam=Lam, M=M)
out = []
for seed, M in [(1, 1.0), (2, 1.0), (3, 1.0), (4, 1.0), (5, 1.3), (6, 0.8)]:
    w_, per, par = run_case(seed, M); out.append(dict(seed=seed, worst=w_, per=per, par=par)); print(seed, M, 'worst', w_, per, flush=True)
W = max(o['worst'] for o in out); G = max(o['per']['gauss'] for o in out)
print('MAX residual (rhs vs sympy EL):', W, ' Gauss EL_a0:', G, 'PASS' if W <= 1e-12 and G <= 1e-12 else 'FAIL')
json.dump(dict(max_residual=W, max_gauss=G, cases=out, tolerance=1e-12), open('selftest_c4_el_result.json', 'w'), indent=1)
