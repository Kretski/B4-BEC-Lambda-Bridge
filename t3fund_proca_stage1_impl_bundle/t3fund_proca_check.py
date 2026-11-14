"""Verification of the analytic statements of T3_FUND_PROCA_SPEC.md (run before registration; no population evaluation).
Model (1+1D, metric (-,+), M=1 unless stated):
 L = 1/2 f_t^2 - 1/2 f_x^2 - (Lam/4) f_x^4 + 1/2 g_t^2 - 1/2 g_x^2 + 1/2 (A1_t - A0_x)^2 + 1/2 M^2 (A0^2 - A1^2) + A1 (gphi f_t f_x + gpsi g_t g_x)
Checks (IMPLEMENTED in THIS script, executed in full):
 (1) C0: Euler-Lagrange residuals of the four fields on the background f=wt+qx, g=vt+rx, A0=0, A1=Abar=(gphi w q+gpsi v r)/M^2 vanish (symbolic, full EL operators incl. t- and x-derivative terms).
 (2) det K(omega=c k,k=M x) = -M^8 x^4 (1+x^2) Q(c;eps), eps=x^2/(1+x^2) (exact rational instances, independent linearisation by plane waves).
Checks NOT in this script (done elsewhere, listed for reference):
 (3) reduced evolution system vs Euler-Lagrange equations at random smooth states (<=1e-12), (4) conservation of the energy E by the reduced system, (5) 6x6 linearised matrix vs roots of Q (50 digits):
     these are in t3fund_proca_check_numeric.py / t3fund_proca_check_mu_mp.py of the registered spec bundle (commit e355af2), not here. The nonlinear C4 code and its <=1e-12 check against the SymPy EL operator do not exist yet.
"""
import sympy as sp, numpy as np, mpmath as mp, sys, json
t, x = sp.symbols('t x', real=True)
M, gphi, gpsi, w, q, v, r, Lam = sp.symbols('M gphi gpsi w q v r Lam', real=True)
k, om, c = sp.symbols('k omega c')
f = sp.Function('f')(t, x); g = sp.Function('g')(t, x); a0 = sp.Function('a0')(t, x); a1 = sp.Function('a1')(t, x)
half = sp.Rational(1, 2)
def Lag(F_, G_, A0_, A1_):
    return (half*sp.diff(F_, t)**2 - half*sp.diff(F_, x)**2 - Lam/4*sp.diff(F_, x)**4 + half*sp.diff(G_, t)**2 - half*sp.diff(G_, x)**2
            + half*(sp.diff(A1_, t) - sp.diff(A0_, x))**2 + half*M**2*(A0_**2 - A1_**2)
            + A1_*(gphi*sp.diff(F_, t)*sp.diff(F_, x) + gpsi*sp.diff(G_, t)*sp.diff(G_, x)))
res = {}
# ---- (1) C0 (symbolic, expand) ----
Abar = (gphi*w*q + gpsi*v*r)/M**2
Lfull = Lag(f, g, a0, a1)
eqs = sp.calculus.euler.euler_equations(Lfull, [f, g, a0, a1], [t, x])
bg = {f: w*t + q*x, g: v*t + r*x, a0: sp.Integer(0), a1: Abar}
resid = [sp.expand(e.lhs.subs(bg).doit() - e.rhs.subs(bg).doit()) for e in eqs]
res['C0_residuals_symbolic'] = [str(z) for z in resid]
print('C0 residuals (symbolic):', resid, flush=True)
# ---- (2) det K identity, exact rational instances ----
e_ = sp.Symbol('e'); X = sp.Symbol('X', positive=True)
import random
random.seed(7); worst = 0
for trial in range(8):
    vals = {M: sp.Rational(random.randint(1, 5), 2), gphi: sp.Rational(random.randint(-5, 5), 3), gpsi: sp.Rational(random.randint(-5, 5), 3), w: sp.Rational(random.randint(-6, 6), 4),
            q: sp.Rational(random.randint(-6, 6), 4), v: sp.Rational(random.randint(-6, 6), 4), r: sp.Rational(random.randint(-6, 6), 4), Lam: sp.Rational(random.randint(0, 4), 3)}
    wv, qv, vv, rv, Lv, Mv_, gpv, gsv = [vals[z] for z in (w, q, v, r, Lam, M, gphi, gpsi)]
    Abv = (gpv*wv*qv + gsv*vv*rv)/Mv_**2
    def Li(F_, G_, A0_, A1_):
        return (half*sp.diff(F_, t)**2 - half*sp.diff(F_, x)**2 - Lv/4*sp.diff(F_, x)**4 + half*sp.diff(G_, t)**2 - half*sp.diff(G_, x)**2
                + half*(sp.diff(A1_, t) - sp.diff(A0_, x))**2 + half*Mv_**2*(A0_**2 - A1_**2) + A1_*(gpv*sp.diff(F_, t)*sp.diff(F_, x) + gsv*sp.diff(G_, t)*sp.diff(G_, x)))
    Lb = Li(wv*t + qv*x + e_*f, vv*t + rv*x + e_*g, e_*a0, Abv + e_*a1)
    L1 = sp.expand(sp.diff(Lb, e_).subs(e_, 0)); L2 = sp.expand(sp.diff(Lb, e_, 2).subs(e_, 0)/2)
    eq2 = sp.calculus.euler.euler_equations(L2, [f, g, a0, a1], [t, x])
    amp = sp.symbols('Uf Ug Ua0 Ua1'); ph = sp.exp(sp.I*(k*x - om*t))
    sub = {f: amp[0]*ph, g: amp[1]*ph, a0: amp[2]*ph, a1: amp[3]*ph}
    K = sp.Matrix([[sp.expand(eq.lhs.subs(sub).doit()/ph).coeff(a) for a in amp] for eq in eq2])
    dK = sp.expand(K.det().subs(om, c*k).subs(k, Mv_*X))
    pp, ss, kap = gpv**2/Mv_**2, gsv**2/Mv_**2, -gpv*gsv/Mv_**2
    aphi = pp*wv*qv - kap*vv*rv; apsi = -kap*wv*qv + ss*vv*rv
    d1 = c**2 - 2*aphi*c - (1 + 3*Lv*qv**2); d2 = c**2 - 2*apsi*c - 1
    Fc = d1*d2 + pp*(wv - qv*c)**2*d2 + ss*(vv - rv*c)**2*d1
    Qi = Fc - (X**2/(1 + X**2))*c**2*d1*d2
    resid2 = sp.expand(sp.together(dK - (-Mv_**8*X**4*(1 + X**2))*Qi))
    num_ = sp.numer(sp.together(resid2))
    print('trial', trial, 'L1 (linear terms on bg) =', L1, ' residual =', sp.simplify(num_), flush=True)
    worst = max(worst, abs(sp.simplify(num_)) if sp.simplify(num_).is_number else 1)
res['detK_identity_rational_instances_max_residual'] = str(worst)
json.dump(res, open('t3fund_proca_check_symbolic.json', 'w'))
print('done')
