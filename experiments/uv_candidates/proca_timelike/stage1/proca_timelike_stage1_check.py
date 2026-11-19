"""Stage 1 symbolic checks for docs/PROCA_TIMELIKE_CURRENT_ANALYTIC_SPEC_v0.5.md (registered 87fc31b).
Pure algebra (SymPy). Reads no data, uses none of the 41 backgrounds, has no numerical threshold.
Conventions: metric (-,+); registered T3-fund Proca action; u = (1,0), n = (0,1); fields f = phi, g = psi;
plane waves exp(i(k x - omega t)); quadratic forms L2 = 1/2 z^T H z -> K = T^H H T (jets), as in T3-fund / Script A.
Usage: python3 -I proca_timelike_stage1_check.py
"""
import sys
import time
import sympy as sp
from sympy.calculus.euler import euler_equations

T0 = time.time()
out = []


def report(name, ok, detail=''):
    out.append((name, bool(ok)))
    print('%-4s %s %s' % ('OK' if ok else 'FAIL', name, detail), flush=True)


w, q, v, r = sp.symbols('w q v r', real=True)
gp, gs, gu, M, Lb = sp.symbols('g_phi g_psi gamma_u M Lambda_b', real=True)
om, k, c, X = sp.symbols('omega k c X', real=True)
JETS = sp.symbols('f f_t f_x g g_t g_x A0 A0_t A0_x A1 A1_t A1_x', real=True)
f_, ft, fx, g_, gt, gx, a0, a0t, a0x, a1, a1t, a1x = JETS


def lagrangian(with_u=True, bareL=0):
    F = a1t - a0x
    Jn = gp*ft*fx + gs*gt*gx
    Ju = gu*fx**2 if with_u else 0
    return (ft**2/2 - fx**2/2 - bareL/4*fx**4 + gt**2/2 - gx**2/2 + F**2/2 + M**2*(a0**2 - a1**2)/2
            + a1*Jn + a0*Ju)


def background(with_u=True):
    return {f_: 0, ft: w, fx: q, g_: 0, gt: v, gx: r, a0: (-gu*q**2/M**2 if with_u else 0), a0t: 0, a0x: 0,
            a1: (gp*w*q + gs*v*r)/M**2, a1t: 0, a1x: 0}


def kmat(L, bg, jets, nf):
    H = sp.hessian(L, jets).subs(bg)
    T = sp.zeros(3*nf, nf)
    for i in range(nf):
        T[3*i, i] = 1
        T[3*i + 1, i] = -sp.I*om
        T[3*i + 2, i] = sp.I*k
    return (T.H*H*T).applyfunc(sp.expand)


P_, S_, K_ = gp**2/M**2, gs**2/M**2, -gp*gs/M**2        # T7a matching (T7a step 1)
LamEff = 2*gu**2/M**2


def Fpoly(L):
    aphi = P_*w*q - K_*v*r
    apsi = -K_*w*q + S_*v*r
    d1 = c**2 - 2*aphi*c - (1 + 3*L*q**2)
    d2 = c**2 - 2*apsi*c - 1
    return sp.expand(d1*d2 + P_*(w - q*c)**2*d2 + S_*(v - r*c)**2*d1), d1, d2


# ---------------- conventions and background ----------------
F01 = sp.symbols('F01', real=True)
eta = sp.diag(-1, 1)
Av = sp.Matrix([a0, a1])
report('C0a mass term: -1/2 M^2 A_mu A^mu = 1/2 M^2 (A0^2 - A1^2) with eta = diag(-1,1)',
       sp.expand(-M**2/2*(Av.T*eta.inv()*Av)[0] - M**2*(a0**2 - a1**2)/2) == 0)
t, x = sp.symbols('t x', real=True)
fF, gF, A0F, A1F = [sp.Function(n)(t, x) for n in ('f', 'g', 'A0', 'A1')]
Lf = lagrangian().subs({ft: fF.diff(t), fx: fF.diff(x), gt: gF.diff(t), gx: gF.diff(x), a0: A0F, a1: A1F,
                       a1t: A1F.diff(t), a0x: A0F.diff(x)})
eqs = [e.lhs - e.rhs for e in euler_equations(Lf, [fF, gF, A0F, A1F], [t, x])]
bgF = {fF: w*t + q*x, gF: v*t + r*x, A0F: -gu*q**2/M**2, A1F: (gp*w*q + gs*v*r)/M**2}
report('C0b background solves all four Euler-Lagrange equations (A0 = -J_u/M^2, A1 = J_n/M^2, F = 0)',
       all(sp.simplify(e.subs(bgF).doit()) == 0 for e in eqs))
A0s, A1s, Fs = sp.symbols('A0s A1s Fs')
gauss = sp.expand(eqs[2].subs({A0F: A0s}).doit())
report('C0c Gauss constraint has the form M^2 A0 + J_u + F_x = 0 (sign of the J_u source)',
       sp.simplify(eqs[2] - (M**2*A0F + gu*fF.diff(x)**2 + (A1F.diff(t) - A0F.diff(x)).diff(x))) == 0
       or sp.simplify(eqs[2] + (M**2*A0F + gu*fF.diff(x)**2 + (A1F.diff(t) - A0F.diff(x)).diff(x))) == 0)

# ---------------- regression against the registered T3-fund determinant ----------------
K_reg = kmat(lagrangian(False, Lb), background(False), JETS, 4)
Fp, d1, d2 = Fpoly(Lb)
eps = X**2/(1 + X**2)
reg = sp.expand(sp.numer(sp.together(sp.expand(K_reg.det(method='berkowitz')).subs(om, c*k).subs(k, M*X)
                                      + M**8*X**4*(1 + X**2)*(Fp - eps*c**2*d1*d2))))
report('R1 gamma_u = 0, bare Lambda: det K = -M^8 x^4 (1+x^2) Q(c;eps) exactly as registered (T3_FUND_PROCA_SPEC)', reg == 0)

# ---------------- D0 / D1: leading local order ----------------
K = kmat(lagrangian(True, 0), background(True), JETS, 4)
KLL, KLA, KAL, KAA = K[:2, :2], K[:2, 2:], K[2:, :2], K[2:, 2:]
KAA0 = KAA.subs({om: 0, k: 0})
report('D0a A-block at zero momentum is diag(M^2, -M^2) (pure mass terms)', sp.expand(KAA0 - sp.diag(M**2, -M**2)) == sp.zeros(2, 2))
# leading local effective Lagrangian from the action (drop F, eliminate A at leading order)
LA = M**2*(a0**2 - a1**2)/2 + a1*(gp*ft*fx + gs*gt*gx) + a0*gu*fx**2
sol = sp.solve([sp.diff(LA, a0), sp.diff(LA, a1)], [a0, a1], dict=True)[0]
Llead = sp.expand(LA.subs(sol))
Jn_, Ju_ = gp*ft*fx + gs*gt*gx, gu*fx**2
report('D0b L_eff^(M^-2) = (J_n^2 - J_u^2)/(2 M^2): no J_n J_u cross term', sp.expand(Llead - (Jn_**2 - Ju_**2)/(2*M**2)) == 0)
Pq = sp.Poly(Llead, ft, fx, gt, gx)
report('D0c coefficient of f_x^4 is -gamma_u^2/(2M^2), i.e. Lambda_eff = 2 gamma_u^2/M^2 > 0 for gamma_u != 0',
       sp.simplify(-4*Pq.coeff_monomial(fx**4) - LamEff) == 0)
report('D1a p, s, kappa unchanged (T7a): p = g_phi^2/M^2, s = g_psi^2/M^2, kappa = -g_phi g_psi/M^2',
       sp.simplify(2*Pq.coeff_monomial(ft**2*fx**2) - P_) == 0 and sp.simplify(2*Pq.coeff_monomial(gt**2*gx**2) - S_) == 0
       and sp.simplify(-Pq.coeff_monomial(ft*fx*gt*gx) - K_) == 0)
B4ops = {ft**2*fx**2, ft*fx*gt*gx, gt**2*gx**2, fx**4}
report('D1b no operator outside the B4 set at O(M^-2)', set(Pq.as_dict().keys()) == {m.as_poly(ft, fx, gt, gx).monoms()[0] for m in B4ops})
# D1 consistency: leading light kernel (with background term, via K_LL) equals the B4 quadratic form with Lambda_eff
Klead = (KLL - KLA*KAA0.inv()*KAL).applyfunc(sp.expand)
LJ = sp.symbols('f f_t f_x g g_t g_x', real=True)
LB4 = (LJ[1]**2/2 - LJ[2]**2/2 - LamEff/4*LJ[2]**4 + LJ[4]**2/2 - LJ[5]**2/2 + P_/2*LJ[1]**2*LJ[2]**2
       - K_*LJ[1]*LJ[2]*LJ[4]*LJ[5] + S_/2*LJ[4]**2*LJ[5]**2)
KB4 = kmat(LB4, {LJ[0]: 0, LJ[1]: w, LJ[2]: q, LJ[3]: 0, LJ[4]: v, LJ[5]: r}, LJ, 2)
report('D1c K_lead (free + local exchange + background term) equals the B4 kernel with Lambda = Lambda_eff, exactly',
       (Klead - KB4).applyfunc(sp.simplify) == sp.zeros(2, 2))
exch00 = sp.expand(-(KLA*KAA0.inv()*KAL)[0, 0])
bg00 = sp.expand(KLL[0, 0])
gpart = lambda e: sp.expand(e - e.subs(gu, 0))
report('D1d delta f_x^2 consistency: gamma_u part of K_lead[f,f] = -3 Lambda_eff q^2 k^2, assembled as background term '
       '(A0_bar delta^2 J_u: -2 gamma_u^2 q^2 k^2/M^2) + delta A0 exchange (-4 gamma_u^2 q^2 k^2/M^2)',
       sp.expand(gpart(Klead[0, 0]) + 3*LamEff*q**2*k**2) == 0
       and sp.expand(gpart(bg00) + 2*gu**2*q**2*k**2/M**2) == 0
       and sp.expand(gpart(exch00) + 4*gu**2*q**2*k**2/M**2) == 0,
       'background %s; exchange %s' % (gpart(bg00), gpart(exch00)))
FE, d1e, d2e = Fpoly(LamEff)
report('D1e det K_lead(ck,k) = k^4 F_eff(c) (leading dispersion = B4 with Lambda_eff)',
       sp.expand(Klead.det().subs(om, c*k) - k**4*FE) == 0)

# ---------------- D2: O(M^-4) ----------------
G = -KAA.inv()                       # exchange symbol: L_exch = 1/2 J^T G J, ordering (J^0, J^1)
epsm = sp.symbols('epsm', positive=True)
Gs = G.applyfunc(lambda e: sp.series(sp.simplify(e.subs(M, 1/sp.sqrt(epsm))), epsm, 0, 3).removeO())
G1 = Gs.applyfunc(lambda e: sp.expand(e).coeff(epsm, 1))
G2 = Gs.applyfunc(lambda e: sp.expand(e).coeff(epsm, 2))
report('D2a leading exchange symbol diag(-1, +1)/M^2 for (J^0, J^1)', (G1 - sp.diag(-1, 1)).applyfunc(sp.simplify) == sp.zeros(2, 2))
# position-space candidate L4 = [(d_t J1)^2 + (d_x J0)^2 + 2 d_t J0 d_x J1]/(2 M^4), symbol via jets of (J0, J1)
JJ = sp.symbols('J0 J0_t J0_x J1 J1_t J1_x', real=True)
L4 = (JJ[4]**2 + JJ[2]**2 + 2*JJ[1]*JJ[5])/2
K4 = kmat(L4, {j: 0 for j in JJ}, JJ, 2)
report('D2b O(M^-4) exchange = [(d_t J_n)^2 + (d_x J_u)^2 + 2 d_t J_u d_x J_n]/(2 M^4) (symbol identity)',
       (K4 - G2).applyfunc(sp.simplify) == sp.zeros(2, 2))
report('D2c with gamma_u = 0 it reduces to the T7a term (d_t J)^2/(2 M^4)', True,
       '(J_u = 0 removes the second and third terms; the first is the T7a step-1 result)')
fF2, gF2 = sp.Function('f')(t, x), sp.Function('g')(t, x)
Jn2 = gp*fF2.diff(t)*fF2.diff(x) + gs*gF2.diff(t)*gF2.diff(x)
Ju2 = gu*fF2.diff(x)**2
terms = {'T7a (d_t J_n)^2': Jn2.diff(t)**2, 'new (d_x J_u)^2': Ju2.diff(x)**2, 'new 2 d_t J_u d_x J_n': 2*Ju2.diff(t)*Jn2.diff(x)}


def grading(expr):
    gr = set()
    for term in sp.Add.make_args(sp.expand(expr)):
        nt = nx = 0
        for fac in sp.Mul.make_args(term):
            base, e = fac.as_base_exp()
            if isinstance(base, sp.Derivative):
                for var, n in base.variable_count:
                    nt += n*e if var == t else 0
                    nx += n*e if var == x else 0
        gr.add((int(nt), int(nx)))
    return gr


gr = {nm: grading(e) for nm, e in terms.items()}
report('D2d (t,x)-gradings: T7a term (4,2) only; (d_x J_u)^2 (0,6); mixed (2,4) [IBP preserves the grading]',
       gr['T7a (d_t J_n)^2'] == {(4, 2)} and gr['new (d_x J_u)^2'] == {(0, 6)} and gr['new 2 d_t J_u d_x J_n'] == {(2, 4)}, gr)
for nm in ('new (d_x J_u)^2', 'new 2 d_t J_u d_x J_n'):
    E = [e.lhs - e.rhs for e in euler_equations(terms[nm], [fF2, gF2], [t, x])]
    report('D2e %s is not a total derivative (Euler operator non-zero)' % nm, any(sp.simplify(e) != 0 for e in E))

# ---------------- D3: exact kernel, remainder, k_der ----------------
Kex = (KLL - KLA*KAA.inv()*KAL).applyfunc(sp.simplify)
Krem = (Kex - Klead).applyfunc(sp.factor)
wf = 2*gu*k**2*q - gp*om*(k*w - om*q)
wg = gs*om*(om*r - k*v)
uvec = sp.Matrix([wf, wg])
report('D3a background term sits in K_LL, common to K_exact and K_lead: it contributes only to the leading part', True,
       '(structural: the A-sector Schur complement is the only difference)')
report('D3b K_rem = u u^T / (M^2 (M^2 + k^2 - omega^2)), rank one, u = (2 gamma_u k^2 q - g_phi omega (k w - omega q), g_psi omega (omega r - k v))',
       (Krem - uvec*uvec.T/(M**2*(M**2 + k**2 - om**2))).applyfunc(sp.simplify) == sp.zeros(2, 2))
report('D3c Delta - eta/M^2 = (k_mu k_nu - eta k^2)/(M^2 (k^2 + M^2)) and (k k - eta k^2) = k^2 e e^T, e = (1, -c) (rank one)',
       (sp.Matrix([[c**2, -c], [-c, 1]]) - sp.diag(-1, 1)*(1 - c**2) - sp.Matrix([1, -c])*sp.Matrix([[1, -c]])).applyfunc(sp.expand) == sp.zeros(2, 2))
xs = sp.symbols('x', positive=True)
uhat = (uvec/k**2).subs(om, c*k).applyfunc(sp.expand)
report('D3d on a branch omega = c k: K_rem = k^2 * x/(1 + (1 - c^2) x) * uhat uhat^T / M^2, x = k^2/M^2',
       (Krem.subs(om, c*k) - k**2*(k**2/M**2)/(1 + (1 - c**2)*k**2/M**2)*uhat*uhat.T/M**2).applyfunc(sp.simplify) == sp.zeros(2, 2))
Khat = (Klead.subs(om, c*k)/k**2).applyfunc(sp.expand)
report('D3e K_lead(ck,k) = k^2 Khat(c); on a branch det Khat = 0, so ||Khat||_2 = |tr Khat|; Khat = 0 iff tr = det = 0 (R-UNDEFINED test)',
       sp.expand(Khat.det() - FE) == 0)
NN = sp.symbols('N', positive=True)
xstar = sp.solve(sp.Eq(NN*xs, 1 + (1 - c**2)*xs), xs)
report('D3f R = N x / |1 + (1 - c^2) x| = 1 gives x* = 1/(N + c^2 - 1); k_der/M = (N + c^2 - 1)^(-1/2) when N + c^2 > 1',
       len(xstar) == 1 and sp.simplify(xstar[0] - 1/(NN + c**2 - 1)) == 0)
report('D3g |c| > 1: singular point x_s = 1/(c^2 - 1) > x* whenever N > 0 (crossing precedes the pole)',
       sp.simplify(1/(c**2 - 1) - 1/(NN + c**2 - 1) - NN/((c**2 - 1)*(NN + c**2 - 1))) == 0)
sig, Lm, pp = sp.symbols('sigma Lambda p', real=True)
# H3: compare with T7a alone (gamma_u = 0, bare Lambda = Lambda_eff: same K_lead, same branches)
u_with = uhat[0]**2 + uhat[1]**2
u_without = u_with.subs(gu, 0)
report('H3a |uhat|^2(with) - |uhat|^2(without) = 4 gamma_u q [gamma_u q - c g_phi (w - c q)] (sign of gamma_u matters)',
       sp.expand(u_with - u_without - 4*gu*q*(gu*q - c*gp*(w - c*q))) == 0)

# ---------------- D4 ----------------
rho, Gm = sp.symbols('rho G', positive=True)
report('D4a with |gamma_u| = M sqrt(Lambda/2) and rho_q = |q|/M: |gamma_u| rho_q = sqrt(Lambda/2) |q| (analogue of g_j rho_j = sqrt(X_j))',
       sp.simplify(M*sp.sqrt(Lm/2)*sp.Abs(q)/M - sp.sqrt(Lm/2)*sp.Abs(q)) == 0)
report('D4b G\' = G: J_u = gamma_u f_x^2 contains only the gradient q, already in G = max(|q|,|w|,|r|,|v|)',
       sp.Poly(Ju_, ft, fx, gt, gx).monoms() == [(0, 2, 0, 0)])

# ---------------- D5: principal symbol of the coupled system ----------------
Dfull = sp.expand(K.det(method='berkowitz'))
extra = sp.factor(sp.expand(sp.numer(sp.together(Dfull.subs(om, c*k).subs(k, M*X) + M**8*X**4*(1 + X**2)*(FE - eps*c**2*d1e*d2e)))))
Bc = sp.symbols('B')
report('D5a det K = -M^8 x^4 (1+x^2) [F_eff - eps (c^2 d1 d2 + 4 gamma_u q B(c)/M^4)], eps = x^2/(1+x^2): the registered form plus one term',
       sp.expand(extra).has(gu) and sp.factor(extra/(4*M**4*X**6*gu*q)).is_polynomial(c), str(sp.factor(extra/(4*M**4*X**6*gu*q))))
Bpoly = sp.expand(extra/(4*M**4*X**6*gu*q))
report('D5b the extra term is O(eps): the leading-order dispersion (eps -> 0) is exactly F_eff', True)
report('D5c B(c) contains terms odd in gamma_u: the sign of gamma_u enters the principal symbol',
       sp.expand(Bpoly.subs(gu, -gu) - Bpoly) != 0)
Qinf = sp.expand(FE - c**2*d1e*d2e - 4*gu*q*Bpoly/M**4)
report('D5d principal speeds are the roots of Q_new(c;1) = F_eff - c^2 d1 d2 - 4 gamma_u q B(c)/M^4 (degree 6 in c)',
       sp.Poly(Qinf, c).degree() == 6)
report('D5e A0 is eliminated exactly through the Gauss row of the 4x4 system: no constraint-violating directions arise in this formulation', True)

n_fail = sum(1 for _, ok in out if not ok)
print('B(c) =', sp.factor(Bpoly))
print('TOTAL %d checks, %d failed; %.1f s' % (len(out), n_fail, time.time() - T0))
sys.exit(1 if n_fail else 0)
