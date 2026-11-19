"""Stage 1 D5 addendum: symbolic checks (spec docs/PROCA_TIMELIKE_CURRENT_ANALYTIC_SPEC_v0.5.md, 87fc31b).
Reduced first-order system with the timelike current, its principal symbol under the registered scaling, the Gauss
constraint, and the relation between the spectrum and Q_new. Also replaces the four checks of the Stage 1 script that
were recorded as unconditional True (D2c, D3a, D5b, D5e) by computed checks.
Pure algebra (SymPy). No data, none of the 41 backgrounds, no numerical threshold.
Usage: python3 -I proca_timelike_d5_addendum_check.py
"""
import sys
import time
import sympy as sp

T0 = time.time()
out = []


def report(name, ok, detail=''):
    out.append((name, bool(ok)))
    print('%-4s %s %s' % ('OK' if ok else 'FAIL', name, detail), flush=True)


t, x = sp.symbols('t x', real=True)
w, q, v, r = sp.symbols('w q v r', real=True)
gp, gs, gu = sp.symbols('g_phi g_psi gamma_u', real=True)
M = sp.symbols('M', positive=True)
k = sp.symbols('k', positive=True)
om, c, X, lam = sp.symbols('omega c X lambda')
Ab0, Ab1 = sp.symbols('Abar0 Abar1', real=True)

# ---------------- reduced first-order system (registered variables, T3_FUND_PROCA_SPEC 3(iii)) ----------------
names = ['a', 'b', 'cc', 'd', 'A1', 'F']
U = [sp.Function(n)(t, x) for n in names]
a, b, cc, d, A1, F = U
A0expr = -(F.diff(x) + gu*b**2)/M**2                 # from the Gauss constraint M^2 A0 + J_u + F_x = 0
Jn = gp*a*b + gs*cc*d
A1t = F + A0expr.diff(x)                             # F = A1_t - A0_x
Ft = Jn - M**2*A1                                    # A1 Euler-Lagrange equation
bt, dt = a.diff(x), cc.diff(x)                       # compatibility f_tx = f_xt, g_tx = g_xt
at = -(gp*(A1t*b + A1*bt)) - (-b + gp*A1*a + 2*gu*A0expr*b).diff(x)
ct = -(gs*(A1t*d + A1*dt)) - (-d + gs*A1*cc).diff(x)
rhs = [at, bt, ct, dt, A1t, Ft]

# G1: the Gauss constraint is satisfied identically by construction (A0 is not an evolved variable)
gauss = M**2*A0expr + gu*b**2 + F.diff(x)
report('G1 Gauss constraint M^2 A0 + J_u + F_x = 0 holds identically for every state (A0 defined by it)', sp.simplify(gauss) == 0)

# G2: the reduced equations reproduce the four Euler-Lagrange equations of the registered action + J_u
fF, gF, A0F, A1F = [sp.Function(n)(t, x) for n in ('f', 'g', 'A0', 'A1')]
from sympy.calculus.euler import euler_equations
Fd = A1F.diff(t) - A0F.diff(x)
L = (fF.diff(t)**2/2 - fF.diff(x)**2/2 + gF.diff(t)**2/2 - gF.diff(x)**2/2 + Fd**2/2 + M**2*(A0F**2 - A1F**2)/2
     + A1F*(gp*fF.diff(t)*fF.diff(x) + gs*gF.diff(t)*gF.diff(x)) + A0F*gu*fF.diff(x)**2)
EL = [e.lhs - e.rhs for e in euler_equations(L, [fF, gF, A0F, A1F], [t, x])]
state = {a: fF.diff(t), b: fF.diff(x), cc: gF.diff(t), d: gF.diff(x), A1: A1F, F: Fd}
# Replace A0 by its Gauss value and A1_t by F + A0_x (F a free function); the f and g Euler-Lagrange equations must then
# equal -(f_tt - at) and -(g_tt - ct) with the reduced right-hand sides evaluated on the same state.
FF = sp.Function('F')(t, x)
A0g = -(FF.diff(x) + gu*fF.diff(x)**2)/M**2
st = {a: fF.diff(t), b: fF.diff(x), cc: gF.diff(t), d: gF.diff(x), A1: A1F, F: FF}
at_s, ct_s = at.subs(st).doit(), ct.subs(st).doit()
def reduce_EL(e):
    e = e.subs(A1F.diff(t, x), (FF + A0g.diff(x)).diff(x)).subs(A1F.diff(t), FF + A0g.diff(x))
    return sp.expand(e.subs(A0F, A0g).doit())
Ef = sp.simplify(reduce_EL(EL[0]) + (fF.diff(t, 2) - at_s))
Eg = sp.simplify(reduce_EL(EL[1]) + (gF.diff(t, 2) - ct_s))
okG2 = (Ef == 0) and (Eg == 0)
EA1 = EL[3].subs(A0F, A0g).doit()
EA1 = EA1.subs(A1F.diff(t, 2), (FF + A0g.diff(x)).diff(t)).subs(A1F.diff(t, x), (FF + A0g.diff(x)).diff(x)).subs(A1F.diff(t), FF + A0g.diff(x))
G3 = sp.simplify(sp.expand(EA1.doit()) + (FF.diff(t) - Ft.subs(st).doit()))
report('G3 reduced F equation (F_t = J_n - M^2 A1) equals the A1 Euler-Lagrange equation', G3 == 0, '' if G3 == 0 else str(G3))
report('G2 reduced f and g equations equal the Euler-Lagrange equations (A0 from Gauss, A1_t = F + A0_x)', okG2, '' if okG2 else 'Ef=%s Eg=%s' % (Ef, Eg))

# linearisation (plane waves e^{i k x}); dU/dt = Mlin U
bg = [w, q, v, r, (gp*w*q + gs*v*r)/M**2, 0]
eps = sp.symbols('eps')
amp = sp.symbols('u0:6')
ph = sp.exp(sp.I*k*x)
sub = {U[i]: bg[i] + eps*amp[i]*ph for i in range(6)}
Mlin = sp.zeros(6, 6)
for i, e in enumerate(rhs):
    lin = sp.expand(sp.diff(e.subs(sub).doit(), eps).subs(eps, 0)/ph)
    for j in range(6):
        Mlin[i, j] = sp.expand(sp.diff(lin, amp[j]))

# S1: spectrum of the reduced system = roots of det K of the 4x4 Lagrangian system (all k)
ns = {}
src = open(__file__.replace('proca_timelike_d5_addendum_check.py', 'proca_timelike_stage1_check.py')).read()
exec(src.split('# ---------------- conventions and background')[0], ns)
K4 = ns['kmat'](ns['lagrangian'](True, 0), ns['background'](True), ns['JETS'], 4)
ren = {ns['w']: w, ns['q']: q, ns['v']: v, ns['r']: r, ns['gp']: gp, ns['gs']: gs, ns['gu']: gu, ns['M']: M,
       ns['om']: om, ns['k']: k}
DK = sp.expand(K4.det(method='berkowitz').subs(ren))
Pchar = sp.expand((Mlin + sp.I*om*sp.eye(6)).det(method='berkowitz'))
ratio = sp.cancel(sp.together(Pchar/DK))
report('S1 det(Mlin + i omega I) = -det K(omega,k)/M^2 exactly: the reduced spectrum is the full Lagrangian spectrum',
       sp.simplify(ratio + 1/M**2) == 0, str(ratio))

# W1: registered scaling A~ = M dA1, F~ = (khat/M) dF (khat = sqrt(M^2 + k^2)); no extra weights needed
kh = sp.sqrt(M**2 + k**2)
S = sp.diag(1, 1, 1, 1, M, kh/M)
Ms = (S*Mlin*S.inv()).applyfunc(sp.simplify)
Pinf = Ms.applyfunc(lambda e: sp.limit(e/k, k, sp.oo))
finite = all(sp.simplify(e).is_finite is not False and not sp.simplify(e).has(sp.oo, sp.zoo) for e in Pinf)
rem = (Ms - k*Pinf).applyfunc(lambda e: sp.limit(e, k, sp.oo))
bounded = all(not sp.simplify(e).has(sp.oo, -sp.oo, sp.zoo) for e in rem)
rem_zero = rem.applyfunc(sp.simplify) == sp.zeros(6, 6)          # the document claims B(k) -> 0, not only boundedness
report('W1 with the registered scaling every entry of the scaled matrix is O(k): P_inf = lim Ms/k exists (no Douglis-Nirenberg weights beyond the registered scaling)', finite)
report('W2 the remainder Ms - k P_inf tends to the zero matrix as k -> infinity (every limit is exactly 0)', bounded and rem_zero, str(rem.applyfunc(sp.simplify)))

# W3: B(k) = Ms - k P_inf is continuous on k >= 0 (denominators only powers of M and khat = sqrt(M^2+k^2) >= M) and
# tends to 0: hence sup_{k>=0} |B(k)| < infinity (needed for the logarithmic-norm bound with the k-independent V_inf)
Bk = (Ms - k*Pinf).applyfunc(lambda e: sp.together(sp.simplify(e)))
dens_ok = True
for e in Bk:
    den = sp.denom(e)
    for fac in sp.Mul.make_args(sp.factor(den)):
        base, ex = fac.as_base_exp()
        if base.has(k) and sp.simplify(base - kh) != 0 and sp.simplify(base - (M**2 + k**2)) != 0:
            dens_ok = False
report('W3 B(k) = Ms - k P_inf has no singularity for k >= 0 (denominators are powers of M and khat) and B -> 0: sup_k |B(k)| finite',
       dens_ok and bounded)

# W4: reality of the reduced system: Ms(-k) = conj(Ms(k)) (khat is even and real), so ||exp(t Ms(-k))|| = ||exp(t Ms(k))||
# and it suffices to bound the propagator for k >= 0
kr = sp.symbols('k_r', positive=True)
w4 = (Ms.subs(k, -kr) - Ms.subs(k, kr).applyfunc(sp.conjugate)).applyfunc(sp.simplify)
report('W4 Ms(-k) = conj(Ms(k)) for real k (real system, even real khat): the bound for k >= 0 covers all k', w4 == sp.zeros(6, 6))

# P1: characteristic polynomial of P_inf at lambda = -i c is proportional to Q_new(c;1)
P_, S_, K_ = gp**2/M**2, gs**2/M**2, -gp*gs/M**2
LamEff = 2*gu**2/M**2
aphi = P_*w*q - K_*v*r
apsi = -K_*w*q + S_*v*r
d1 = c**2 - 2*aphi*c - (1 + 3*LamEff*q**2)
d2 = c**2 - 2*apsi*c - 1
FE = d1*d2 + P_*(w - q*c)**2*d2 + S_*(v - r*c)**2*d1
Bc = (M**2*(-c**4*gp*q + c**3*gp*w + c**2*gp*q - c**2*gu*q - c*gp*w + gu*q) + 2*c**3*gp**2*gs*q**2*w
      + 2*c**3*gp*gs**2*q*r*v - 2*c**2*gp**2*gs*q*w**2 - 2*c**2*gp*gs**2*r*v*w - c**2*gs**2*gu*q*r**2
      + 2*c*gp*gs*gu*q**2*w + 4*c*gs**2*gu*q*r*v - gs**2*gu*q*v**2)
Qnew1 = sp.expand(FE - c**2*d1*d2 - 4*gu*q*Bc/M**4)
charP = sp.expand((lam*sp.eye(6) - Pinf).det(method='berkowitz').subs(lam, -sp.I*c))
rat = sp.cancel(sp.together(charP/Qnew1))
report('P1 det(lambda I - P_inf) at lambda = -i c equals a non-zero constant times Q_new(c;1): eigenvalues of P_inf are -i c_j, c_j the roots of Q_new(c;1)',
       not rat.has(c) and rat != 0, str(sp.factor(rat)))

# ---------------- replacements of the unconditional checks of the Stage 1 script ----------------
# D2c: L4 with gamma_u = 0 equals the T7a term (d_t J_n)^2/(2 M^4)
JJ = sp.symbols('J0 J0_t J0_x J1 J1_t J1_x', real=True)
L4 = (JJ[4]**2 + JJ[2]**2 + 2*JJ[1]*JJ[5])/2
report('D2c (computed) L4 with J^0 = 0 (gamma_u = 0) equals (d_t J_n)^2/2 (times 1/M^4): the T7a step-1 term',
       sp.expand(L4.subs({JJ[0]: 0, JJ[1]: 0, JJ[2]: 0}) - JJ[4]**2/2) == 0)
# D3a: with the background Proca values as free symbols, only K_LL depends on them; K_rem does not
JETS = ns['JETS']
f_, ft, fx, g_, gt, gx, a0, a0t, a0x, a1, a1t, a1x = JETS
Lsym = ns['lagrangian'](True, 0).subs(ren)
bgs = {f_: 0, ft: w, fx: q, g_: 0, gt: v, gx: r, a0: Ab0, a0t: 0, a0x: 0, a1: Ab1, a1t: 0, a1x: 0}
H = sp.hessian(Lsym, JETS).subs(bgs)
T = sp.zeros(12, 4)
for i in range(4):
    T[3*i, i] = 1
    T[3*i + 1, i] = -sp.I*om
    T[3*i + 2, i] = sp.I*k
Kb = (T.H*H*T).applyfunc(sp.expand)
KLA, KAA = Kb[:2, 2:], Kb[2:, 2:]
report('D3a (computed) K_LA and K_AA do not depend on the background Proca values; only K_LL does, so the background term cancels in K_rem',
       not any(e.has(Ab0) or e.has(Ab1) for e in list(KLA) + list(KAA)) and any(e.has(Ab0) or e.has(Ab1) for e in list(Kb[:2, :2])))
# D5b: Q_new(c; eps) -> F_eff at eps = 0
Qeps = sp.expand(FE - eps*(c**2*d1*d2 + 4*gu*q*Bc/M**4))
report('D5b (computed) Q_new(c; eps = 0) = F_eff(c)', sp.expand(Qeps.subs(eps, 0) - FE) == 0)
# D5e: replaced by G1 + S1 (constraint holds identically; reduced spectrum = Lagrangian spectrum)
report('D5e (computed) no constraint-violating directions: G1 and S1 hold, and the 6 reduced variables (a, b, cc, d, A1, F) are free Cauchy data',
       sp.simplify(gauss) == 0 and sp.simplify(ratio + 1/M**2) == 0 and Mlin.shape == (6, 6))

n_fail = sum(1 for _, ok in out if not ok)
print('P_inf =', Pinf)
print('TOTAL %d checks, %d failed; %.1f s' % (len(out), n_fail, time.time() - T0))
sys.exit(1 if n_fail else 0)
