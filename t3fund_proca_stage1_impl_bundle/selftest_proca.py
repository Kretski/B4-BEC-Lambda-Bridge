"""Self-test of t3fund_proca_stage1.py against the independent SymPy Euler-Lagrange operator on RANDOM backgrounds (not manifest cases; no population evaluation)."""
import sys, json, numpy as np, sympy as sp, mpmath as mp
sys.path.insert(0, '.')
import t3fund_proca_stage1 as S
rng = np.random.default_rng(20261115)
t, x = sp.symbols('t x', real=True); k, om, c = sp.symbols('k omega c')
f = sp.Function('f')(t, x); g = sp.Function('g')(t, x); a0 = sp.Function('a0')(t, x); a1 = sp.Function('a1')(t, x); e_ = sp.Symbol('e'); half = sp.Rational(1, 2)
report = []; worst = dict(Q=0, detK=0, hess_ratio=0, M_vs_Q=0, k0=0)
for trial in range(8):
    gp, gs = [float(z) for z in rng.uniform(-1, 1, 2)]; Lam = float(rng.uniform(0, 1)); w, q, v, r = [float(z) for z in rng.uniform(-1.2, 1.2, 4)]
    p, s, kap = gp**2, gs**2, -gp*gs
    case = dict(pop='TEST', scope=True, label='t%d' % trial, p=p, s=s, kappa=kap, Lambda=Lam, w=w, q=q, v=v, r=r, gphi=gp, gpsi=gs)
    P = S.pars(case, 50)
    # independent SymPy linearisation (exact rationals of the binary values), M=1
    Rr = sp.Rational; gpv, gsv, Lv, wv, qv, vv, rv = [Rr(z) for z in (gp, gs, Lam, w, q, v, r)]
    Abv = gpv*wv*qv+gsv*vv*rv
    Lb = (half*sp.diff(wv*t+qv*x+e_*f, t)**2-half*sp.diff(wv*t+qv*x+e_*f, x)**2-Lv/4*sp.diff(wv*t+qv*x+e_*f, x)**4+half*sp.diff(vv*t+rv*x+e_*g, t)**2-half*sp.diff(vv*t+rv*x+e_*g, x)**2
          +half*(sp.diff(Abv+e_*a1, t)-sp.diff(e_*a0, x))**2+half*((e_*a0)**2-(Abv+e_*a1)**2)+(Abv+e_*a1)*(gpv*sp.diff(wv*t+qv*x+e_*f, t)*sp.diff(wv*t+qv*x+e_*f, x)+gsv*sp.diff(vv*t+rv*x+e_*g, t)*sp.diff(vv*t+rv*x+e_*g, x)))
    L2 = sp.expand(sp.diff(Lb, e_, 2).subs(e_, 0)/2)
    eq2 = sp.calculus.euler.euler_equations(L2, [f, g, a0, a1], [t, x])
    amp = sp.symbols('Uf Ug Ua0 Ua1'); ph = sp.exp(sp.I*(k*x-om*t)); sub = {f: amp[0]*ph, g: amp[1]*ph, a0: amp[2]*ph, a1: amp[3]*ph}
    K = sp.Matrix([[sp.expand(eq.lhs.subs(sub).doit()/ph).coeff(a) for a in amp] for eq in eq2]); detK = sp.expand(K.det())
    for xv in (0.7, 3.0, 50.0):
        # (a) Q coefficients vs sympy
        eps = sp.Rational(xv)**2/(1+sp.Rational(xv)**2)
        aphi = sp.Rational(p)*wv*qv-sp.Rational(kap)*vv*rv; apsi = -sp.Rational(kap)*wv*qv+sp.Rational(s)*vv*rv
        d1 = c**2-2*aphi*c-(1+3*Lv*qv**2); d2 = c**2-2*apsi*c-1
        Qs = sp.Poly(sp.expand(d1*d2+sp.Rational(p)*(wv-qv*c)**2*d2+sp.Rational(s)*(vv-rv*c)**2*d1-eps*c**2*d1*d2), c).all_coeffs()[::-1]
        mine = S.Qcoeffs(P, S.eps_of(mp.mpf(xv)))
        # note: sympy uses p=gp^2 etc. as exact rationals of floats; mine uses the same floats converted to mp -> equal up to 50-digit rounding
        qdev = max(abs(mp.mpf(sp.N(a, 50))-b)/max(1, abs(b)) for a, b in zip(Qs, mine)); worst['Q'] = max(worst['Q'], float(qdev))
        # (b) detK implementation (jet Hessian Leibniz) vs sympy plane-wave determinant at (omega,k)=(c0 k, k)
        c0 = 0.37  # real test value
        dre, dim_ = S.detK_poly(P, xv)
        mine_val = sum(dre[i]*mp.mpf(c0)**i for i in range(len(dre)))          # det K/k^8 as polynomial in c
        ref_val = sp.N(detK.subs({om: sp.Rational(c0)*xv, k: sp.Rational(xv)}), 40)/sp.Rational(xv)**8
        worst['detK'] = max(worst['detK'], float(abs(mine_val-mp.mpf(str(sp.re(ref_val))))/max(1e-30, abs(mp.mpf(str(sp.re(ref_val)))))))
    out = S.evaluate(case)
    worst['hess_ratio'] = max(worst['hess_ratio'], out['hessian_route']['max_rel_ratio_dev']); worst['M_vs_Q'] = max(worst['M_vs_Q'], out['C1b']['max_rel_dev_M_vs_Q'])
    worst['k0'] = max(worst['k0'], out['k0']['max_rel_dev'])
    assert out['k0']['n_nonzero'] == 2 and out['k0']['agree']
    worst['C0_EL'] = max(worst.get('C0_EL', 0), max(out['C0'][k] for k in ('EL_f', 'EL_g', 'EL_a0', 'EL_a1')))
    report.append(dict(trial=trial, outcome=out['outcome'], hyper=out['uv']['strictly_hyperbolic'], front=out['uv']['front'], pred3=out['uv'].get('pred3_max_err'), mu=out['mu'], pred4=out['pred4'], n_complex=out['n_complex'], k0_agree=out['k0']['agree'], hidden=out['hidden_by_double'], tiny_imag=out['real_with_tiny_imag'], c0=out['C0']))
    print(trial, report[-1], flush=True)
print('WORST deviations:', worst)
json.dump(dict(worst=worst, trials=report), open('selftest_proca_result.json', 'w'), indent=1, default=str)
