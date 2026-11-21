"""
SOC_PULSE_CROSSING_DELAY_NOTE_v0.2, section 5: one run (approved 2026-11-21).
Part A (symbolic A1-A5), Part B (B0 synthetic kinematics, B1 exact simple wave, B2 preparation time,
B3 SOC scan), positive controls PC1-PC3, decision rows D1-D4.
Part B builds its solutions only from exact lambda_pm(u, y(u)) with y(u) integrated numerically;
symbolic coefficients are used only in the final comparisons.
Units for SOC: E_R = k_R = 1 (m = 1/2, b = 2, a = 2x), mu = gamma.
"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

FAIL = []
REC = {}


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)
    REC[name] = ok


def zero(e):
    return sp.simplify(sp.expand(e)) == 0


def close(num, ref, rel=1e-3, absf=1e-9):
    return abs(num - ref) <= rel*abs(ref) + absf


# =============================================================== Part A
print("=== Part A: symbolic ===", flush=True)
E1, E2, yv, lam = sp.symbols('E1 E2 y lambda', real=True)
A = sp.Matrix([[E1, 1], [yv*E2, E1]])
ev = sp.solve(sp.Eq((A - lam*sp.eye(2)).det(), 0), lam)
forms = [E1 + sp.sqrt(yv*E2), E1 - sp.sqrt(yv*E2)]
okA1 = len(ev) == 2 and all(any(zero(e - f) for f in forms) for e in ev) and \
    all(any(zero(e - f) for e in ev) for f in forms)
rminus = sp.Matrix([1, -sp.sqrt(yv*E2)])
okA1 &= zero(((A - (E1 - sp.sqrt(yv*E2))*sp.eye(2))*rminus)[0]) and \
    zero(((A - (E1 - sp.sqrt(yv*E2))*sp.eye(2))*rminus)[1])
check("A1 eigenvalues E' +- sqrt(y E''); lambda_- right eigenvector (1, -sqrt(y E'')) -> dy/du = -sqrt(y E'')", okA1)

u, M, mu = sp.symbols('u M mu', positive=True)
c3, c4 = sp.symbols('c3 c4', real=True)
c = sp.sqrt(mu/M)
# A2 Galilean
yG = M*(c - u/(2*M))**2
EppG = 1/M
okA2 = zero(sp.diff(yG, u)**2 - yG*EppG) and sp.simplify(sp.diff(yG, u).subs(u, 0) + c) == 0
lpG = u/M + (c - u/(2*M))          # sqrt(y E'') = c - u/2M for u < 2Mc
lmG = u/M - (c - u/(2*M))
okA2 &= zero(sp.expand(lpG - (c + u/(2*M)))) and zero(sp.expand(lmG - (-c + 3*u/(2*M))))
check("A2 Galilean: y = M(c - u/2M)^2 solves A1; f1 = 1/2M, f2 = 0, g1 = 3/2M, g2 = 0", okA2)

# A3 series y(u) to O(u^3)
y2, y3 = sp.symbols('y2 y3')
Et = u**2/(2*M) + c3*u**3 + c4*u**4
Epp = sp.diff(Et, u, 2)
Ep = sp.diff(Et, u)
ys = mu - c*u + y2*u**2 + y3*u**3
resid = sp.expand(sp.diff(ys, u)**2 - ys*Epp)
sol = sp.solve([resid.coeff(u, 1), resid.coeff(u, 2)], [y2, y3], dict=True)
okA3 = len(sol) == 1 and zero(sol[0][y2] - (1/(4*M) - sp.Rational(3, 2)*c3*sp.sqrt(mu*M)))
check("A3 simple wave y = mu - c u + y2 u^2 + y3 u^3, y2 = 1/(4M) - (3/2) c3 sqrt(mu M)", okA3,
      f"y3 = {sp.simplify(sol[0][y3])}")
ys = ys.subs(sol[0])

# A4 f1, f2, g1, g2
root = sp.sqrt(ys*Epp)
lp = sp.series(Ep + root, u, 0, 3).removeO()
lm = sp.series(Ep - root, u, 0, 3).removeO()
f1, f2 = sp.simplify(lp.coeff(u, 1)), sp.simplify(lp.coeff(u, 2))
g1, g2 = sp.simplify(lm.coeff(u, 1)), sp.simplify(lm.coeff(u, 2))
okA4 = zero(lp.coeff(u, 0) - c) and zero(lm.coeff(u, 0) + c)
galsub = {c3: 0, c4: 0}
okA4 &= zero(f1.subs(galsub) - 1/(2*M)) and zero(f2.subs(galsub)) and \
    zero(g1.subs(galsub) - 3/(2*M)) and zero(g2.subs(galsub))
BR = (2*f1**2 - f1*g1)/(8*c**3) - f2/(2*c**2)          # section 4 bracket
BR = sp.simplify(BR)
BRgal = sp.simplify(BR.subs(galsub))
okA4 &= zero(BRgal + 1/(32*c**3*M**2))
print(f"   f1 = {f1}\n   f2 = {f2}\n   g1 = {g1}\n   g2 = {g2}", flush=True)
print(f"   bracket = {BR}\n   Galilean bracket = {BRgal}", flush=True)
check("A4 f1, f2, g1, g2 derived; Galilean limit; Galilean bracket = -1/(32 c^3 M^2)", okA4)
# v0.3 prescription, for comparison only: K2 = (s1^2 - c s2)/c^3 per unit of int u^2 dx
cplus_v03 = Ep.subs(u, sp.Symbol('e')) + sp.sqrt((mu - c*sp.Symbol('e') - Et.subs(u, sp.Symbol('e')))*Epp.subs(u, sp.Symbol('e')))
K2v03 = sp.simplify(sp.series(1/cplus_v03, sp.Symbol('e'), 0, 3).removeO().coeff(sp.Symbol('e'), 2))

# A5 independent symbolic re-derivation with generic f1,f2,g1,g2 and polynomial profile
F1, F2, G1, G2, cc, t1, x0, eps = sp.symbols('F1 F2 G1 G2 cc t1 x0 eps', real=True)
P = x0*(1 - x0)*(1 - 2*x0)
uu = eps*P
lpp = cc + F1*uu + F2*uu**2
lmm = -cc + G1*uu + G2*uu**2
lmm_prime = sp.diff(lmm, x0)
T1s, T2s, sv = sp.symbols('T1s T2s s_dummy', real=True)
tt = t1 + x0/(2*cc) + eps*T1s + eps**2*T2s
rhs = sp.series((1 + lmm_prime*tt)/(lpp - lmm), eps, 0, 3).removeO()
d0 = sp.simplify(rhs.coeff(eps, 0) - 1/(2*cc))
d1 = sp.expand(rhs.coeff(eps, 1))
okA5 = zero(d0) and not d1.has(T1s) and not d1.has(T2s)
tau1 = sp.integrate(d1.subs(x0, sv), (sv, 0, x0))
d2 = sp.expand(rhs.coeff(eps, 2))
okA5 &= not d2.has(T2s)
d2 = sp.expand(d2.subs(T1s, tau1))
tau2_end = sp.integrate(d2, (x0, 0, 1))
dT1 = sp.simplify(2*tau1.subs(x0, 1))
dT2 = sp.simplify(2*tau2_end)
intP2 = sp.integrate(P**2, (x0, 0, 1))
pred = ((2*F1**2 - F1*G1)/(8*cc**3) - F2/(2*cc**2))*intP2
okA5 &= zero(dT1) and zero(dT2 - pred) and zero(sp.diff(dT2, G2)) and zero(sp.diff(dT2, t1))
check("A5 generic expansion: O(eps) = 0, O(eps^2) = section 4 bracket * int P^2, no g2, no t1", okA5,
      f"coefficient of g2 = {sp.simplify(sp.diff(dT2, G2))}")

# =============================================================== Part B machinery
print("\n=== Part B machinery: exact simple wave + probe ODE ===", flush=True)
W = 1.0
S = lambda x: np.sin(2*np.pi*x)*np.sin(np.pi*x)**2                  # int S = 0
Sp = lambda x: 2*np.pi*np.cos(2*np.pi*x)*np.sin(np.pi*x)**2 + np.sin(2*np.pi*x)*2*np.pi*np.sin(np.pi*x)*np.cos(np.pi*x)
S2 = lambda x: np.sin(np.pi*x)**2                                   # int S2 = 1/2 (PC2)
S2p = lambda x: 2*np.pi*np.sin(np.pi*x)*np.cos(np.pi*x)
intS2sq = quad(lambda x: S(x)**2, 0, 1, epsabs=1e-15, epsrel=1e-13)[0]


class Physical:
    """exact lambda_pm(u) from a band (E', E'', E''') and the numerically integrated simple wave."""
    def __init__(self, Ep, Epp, Eppp, muv, U=0.06):
        self.Ep, self.Epp, self.Eppp = Ep, Epp, Eppp
        f = lambda uu, yy: [-np.sqrt(yy[0]*Epp(uu))]
        self.up = solve_ivp(f, (0, U), [muv], method='DOP853', rtol=1e-13, atol=1e-15, dense_output=True)
        self.um = solve_ivp(f, (0, -U), [muv], method='DOP853', rtol=1e-13, atol=1e-15, dense_output=True)
        self.c = np.sqrt(muv*Epp(0.0))

    def y(self, uu):
        return self.up.sol(uu)[0] if uu >= 0 else self.um.sol(uu)[0]

    def lams(self, uu):
        yy = self.y(uu)
        r = np.sqrt(yy*self.Epp(uu))
        lp, lm = self.Ep(uu) + r, self.Ep(uu) - r
        dlm = 1.5*self.Epp(uu) - yy*self.Eppp(uu)/(2*r)
        return lp, lm, dlm


class Synthetic:
    def __init__(self, cv, F1v, F2v, G1v, G2v):
        self.c, self.p = cv, (F1v, F2v, G1v, G2v)

    def lams(self, uu):
        F1v, F2v, G1v, G2v = self.p
        return self.c + F1v*uu + F2v*uu**2, -self.c + G1v*uu + G2v*uu**2, G1v + 2*G2v*uu


def probe_dT(model, ampl, prof=S, profp=Sp, t1v=0.0):
    """delta T = 2 tau(W); also returns min over pulse labels of 1 + lambda_-'(x0) t at t = 0 and t_exit."""
    cv = model.c

    def rhs(x, tau):
        uu = ampl*prof(x)
        lp, lm, dlm = model.lams(uu)
        t = t1v + x/(2*cv) + tau[0]
        return [(1 + dlm*ampl*profp(x)*t)/(lp - lm) - 1/(2*cv)]
    s = solve_ivp(rhs, (0, W), [0.0], method='DOP853', rtol=1e-13, atol=1e-17)
    tau_end = s.y[0, -1]
    t_exit = t1v + W/(2*cv) + tau_end
    xs = np.linspace(0, W, 4001)
    lin = np.array([model.lams(ampl*prof(x))[2]*ampl*profp(x) for x in xs])
    pc3 = min(1.0, float(np.min(1 + lin*t_exit)))      # linear in t: endpoints t = 0 (value 1) and t_exit
    return 2*tau_end, pc3


def richardson(vals, order_start):
    """vals at eps, eps/2, eps/4 of R(eps) = a0 + a1 eps + a2 eps^2 ...; eliminates eps and eps^2."""
    r0, r1, r2 = vals
    q1 = [2*r1 - r0, 2*r2 - r1]
    return (4*q1[1] - q1[0])/3


EPS = [1e-2, 5e-3, 2.5e-3]


def coeff2(model, **kw):
    out, pcs = [], []
    for e in EPS:
        d, pc = probe_dT(model, e, **kw)
        out.append(d/e**2)
        pcs.append(pc)
    return richardson(out, 2), min(pcs)


# SOC band
def soc_band(xv, delv, mirror=False):
    av = 2.0*xv
    Epn = lambda k: 2*k - 2*(2*k - delv/2)/np.sqrt((2*k - delv/2)**2 + av**2)
    k0 = brentq(Epn, -10, 10, xtol=1e-15, rtol=1e-15)

    def parts(uu):
        k = k0 + uu
        Dd = 2*k - delv/2
        R = np.sqrt(Dd**2 + av**2)
        return 2*k - 2*Dd/R, 2 - 4*av**2/R**3, 24*av**2*Dd/R**5
    s = -1.0 if mirror else 1.0
    Ep = lambda uu: s*parts(s*uu)[0]
    Epp = lambda uu: parts(s*uu)[1]
    Eppp = lambda uu: s*parts(s*uu)[2]
    D0 = 2*k0 - delv/2
    R0 = np.sqrt(D0**2 + av**2)
    Mv = 1/(2 - 4*av**2/R0**3)
    c3v = s*4*av**2*D0/R0**5
    c4v = 2*av**2*(av**2 - 4*D0**2)/R0**7
    return Ep, Epp, Eppp, Mv, c3v, c4v


BRf = sp.lambdify((M, mu, c3, c4), BR, 'numpy')
f1f = sp.lambdify((M, mu, c3, c4), f1, 'numpy')
K2f = sp.lambdify((M, mu, c3, c4), K2v03, 'numpy')

# =============================================================== B0 synthetic kinematics
print("\n=== B0: synthetic kinematics (decides the g2 question) ===", flush=True)
okB0i = True
for g1v, g2v in [(1.5, 0.0), (1.5, 2.0)]:
    v, pc = coeff2(Synthetic(1.0, 0.0, 0.0, g1v, g2v))
    okB0i &= abs(v) < 1e-6
    print(f"   f = 0, g1 = {g1v}, g2 = {g2v}: dT/eps^2 -> {v:.3e}   (PC3 min {pc:.4f})", flush=True)
check("B0(i) f = 0 gives dT = O(eps^3) (|dT/eps^2| < 1e-6)", okB0i)
F1v, F2v, G1v = 0.5, 0.3, 1.5
predB0 = ((2*F1v**2 - F1v*G1v)/8 - F2v/2)*intS2sq
vals = {}
for g2v in (-2.0, 0.0, 2.0):
    vals[g2v], pc = coeff2(Synthetic(1.0, F1v, F2v, G1v, g2v))
    print(f"   (f1,f2,g1) = (0.5,0.3,1.5), g2 = {g2v:+.1f}: dT/eps^2 = {vals[g2v]:.8f}  "
          f"(section 4: {predB0:.8f}; PC3 min {pc:.4f})", flush=True)
slope = np.polyfit(list(vals), list(vals.values()), 1)[0]
rev_slope = 0.5*intS2sq
print(f"   fitted slope d(dT/eps^2)/dg2 = {slope:.3e};  review's proposed slope = {rev_slope:.6f}", flush=True)
okB0ii = all(close(v, predB0) for v in vals.values()) and abs(slope) <= 1e-3*abs(predB0)
check("B0(ii) independent of g2 and equal to section 4 bracket (rel 1e-3)", okB0ii)
D4 = abs(slope) > 1e-3*abs(predB0)

# =============================================================== B1 / B2 / PCs
print("\n=== B1: exact simple wave vs section 4 formula ===", flush=True)
cases = []
gal = Physical(lambda uu: uu, lambda uu: 1.0, lambda uu: 0.0, 1.0)
cases.append(("Galilean M=1 mu=1", gal, (1.0, 1.0, 0.0, 0.0)))
for (gam, xv, delv, mir, lab) in [(0.2, 1.2, 0.0, False, "SOC g=0.2 x=1.20 d=0"),
                                  (0.2, 1.1, 0.0, False, "SOC g=0.2 x=1.10 d=0"),
                                  (0.2, 1.2, 0.05, False, "SOC g=0.2 x=1.20 d=0.05 dir +1"),
                                  (0.2, 1.2, 0.05, True, "SOC g=0.2 x=1.20 d=0.05 dir -1")]:
    Ep, Epp, Eppp, Mv, c3v, c4v = soc_band(xv, delv, mir)
    cases.append((lab, Physical(Ep, Epp, Eppp, gam), (Mv, gam, c3v, c4v)))
okB1, okPC3, gal_num = True, True, None
for lab, mod, pars in cases:
    num, pc = coeff2(mod)
    ref = BRf(*pars)*intS2sq
    k2 = K2f(*pars)*intS2sq
    ok = close(num, ref)
    okB1 &= ok
    okPC3 &= pc > 0
    if lab.startswith("Galilean"):
        gal_num = num
    print(f"   {lab:32s}: numeric {num:+.8e}  formula {ref:+.8e}  rel.err {abs(num/ref - 1):.1e}  "
          f"[v0.3 prescription {k2:+.4e}]  PC3 min {pc:.4f}", flush=True)
check("B1 exact simple wave agrees with section 4 formula (rel 1e-3) in all cases", okB1)

print("\n=== B2: preparation time t1 = 0 vs t1 = W/(4c) ===", flush=True)
okB2 = True
for lab, mod, pars in cases:
    a0, p0 = coeff2(mod, t1v=0.0)
    a1, p1 = coeff2(mod, t1v=W/(4*mod.c))
    okB2 &= close(a1, a0)
    okPC3 &= min(p0, p1) > 0
    print(f"   {lab:32s}: {a0:+.8e} vs {a1:+.8e}  (PC3 min {min(p0, p1):.4f})", flush=True)
check("B2 coefficient independent of preparation time (rel 1e-3)", okB2)

print("\n=== Positive controls ===", flush=True)
okPC1 = True
for lab, mod, pars in cases:
    d, _ = probe_dT(mod, 0.0)
    okPC1 &= d == 0.0
check("PC1 no pulse -> dT = 0 exactly", okPC1)
okPC2 = True
for lab, mod, pars in cases[:2]:
    vals2 = []
    for e in EPS:
        d, pc = probe_dT(mod, e, prof=S2, profp=S2p)
        vals2.append(d/e)
        okPC3 &= pc > 0
    lin = richardson(vals2, 1)
    ref = -(f1f(*pars)/(2*mod.c**2))*0.5
    okPC2 &= close(lin, ref)
    print(f"   {lab:32s}: linear coeff {lin:+.8e}  expected {ref:+.8e}", flush=True)
check("PC2 pulse with int u != 0: linear coefficient -(f1/2c^2) int(u/eps) (rel 1e-3)", okPC2)
check("PC3 no characteristic crossing over the whole pulse up to probe exit (all runs)", okPC3)

# =============================================================== B3 scan (formula; c3 = 0 at delta = 0)
print("\n=== B3: SOC scan at delta = 0 (per unit int u^2 dx0) ===", flush=True)
for gam in (0.1, 0.2, 0.5):
    print(f"\n gamma = {gam}\n   x    | bracket SOC   | Galilean same M,mu | SOC - Galilean | v0.3 prescription", flush=True)
    for xv in np.round(np.arange(1.02, 1.5001, 0.02), 2):
        _, _, _, Mv, c3v, c4v = soc_band(xv, 0.0)
        b = BRf(Mv, gam, 0.0, c4v)
        bg = BRf(Mv, gam, 0.0, 0.0)
        print(f"   {xv:.2f} | {b:+.6e} | {bg:+.6e}      | {b - bg:+.6e}  | {K2f(Mv, gam, 0.0, c4v):+.6e}", flush=True)

# =============================================================== decision
print("\n=== Decision (pre-registered table, v0.2 section 6) ===", flush=True)
partA = all(REC[k] for k in REC if k.startswith("A"))
partB = REC.get("B1 exact simple wave agrees with section 4 formula (rel 1e-3) in all cases", False)
pcs = all(REC[k] for k in REC if k.startswith("PC"))
if not pcs:
    print("positive control failed -> run FAIL, not interpreted")
elif D4:
    print("D4: non-zero g2 slope in B0")
elif not (partA and partB):
    print("D3: A and B disagree -> FAIL, not interpreted")
else:
    galB = float(BRgal.subs({M: 1, mu: 1}))
    print(f"Galilean bracket: symbolic {galB:+.6e}, numeric (B1) {gal_num/intS2sq:+.6e}")
    print("D1: Galilean advance confirmed" if galB < 0 and gal_num < 0 else "D2: Galilean delay")
    print("D4 not triggered: B0 slope consistent with zero")
print("\nALL CHECKS:", "PASS" if not FAIL else f"FAIL {FAIL}", flush=True)
sys.exit(0 if not FAIL else 1)
