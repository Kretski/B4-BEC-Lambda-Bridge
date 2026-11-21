"""
Symbolic check of the hand rewrite in SOC_DETUNING_AND_PULSE_CROSSING_REPORT_v0.1 section II.4.
Requested 2026-11-21 14:03 before registration. One run, python3 -I.
  R1  printed A4 bracket (text of registered stdout d00758a4)  ==  -(1 + 12Q - 5 tau^2)/(32 c^3 M^2)
  R2  bracket recomputed from scratch (simple wave, f1, f2, g1)  ==  same rewrite
  R3  c3 = 0: bracket / Galilean bracket (same M, mu)  ==  1 + 12Q
  R4  numeric: rewrite * int S^2 reproduces the B1 'formula' values of stdout d00758a4 (rel 1e-8)
      Q = 8 mu M^2 c4,  tau = 6 c3 M sqrt(mu M),  c = sqrt(mu/M),  mu = g n
"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq

FAIL = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


M, mu = sp.symbols('M mu', positive=True)
c3, c4, u, y2, y3 = sp.symbols('c3 c4 u y2 y3', real=True)
c = sp.sqrt(mu/M)
Q = 8*mu*M**2*c4
tau = 6*c3*M*sp.sqrt(mu*M)
rewrite = -(1 + 12*Q - 5*tau**2)/(32*c**3*M**2)

printed = sp.sympify("(180*M**3*c3**2*mu - 96*M**2*c4*mu - 1)/(32*sqrt(M)*mu**(3/2))",
                     locals={'M': M, 'mu': mu, 'c3': c3, 'c4': c4})
check("R1 printed A4 bracket == -(1 + 12Q - 5 tau^2)/(32 c^3 M^2)", sp.simplify(printed - rewrite) == 0)

Et = u**2/(2*M) + c3*u**3 + c4*u**4
Epp, Ep = sp.diff(Et, u, 2), sp.diff(Et, u)
ys = mu - c*u + y2*u**2 + y3*u**3
res = sp.expand(sp.diff(ys, u)**2 - ys*Epp)
sol = sp.solve([res.coeff(u, 1), res.coeff(u, 2)], [y2, y3], dict=True)[0]
ys = ys.subs(sol)
lp = sp.series(Ep + sp.sqrt(ys*Epp), u, 0, 3).removeO()
lm = sp.series(Ep - sp.sqrt(ys*Epp), u, 0, 3).removeO()
f1, f2, g1 = lp.coeff(u, 1), lp.coeff(u, 2), lm.coeff(u, 1)
bracket = (2*f1**2 - f1*g1)/(8*c**3) - f2/(2*c**2)
check("R2 bracket recomputed from the simple wave == rewrite", sp.simplify(bracket - rewrite) == 0)

ratio = sp.simplify(bracket.subs(c3, 0)/bracket.subs({c3: 0, c4: 0}))
check("R3 c3 = 0: bracket / Galilean bracket (same M, mu) == 1 + 12Q", sp.simplify(ratio - (1 + 12*Q)) == 0,
      f"ratio = {ratio}")

# R4 numeric against registered B1 'formula' values (stdout d00758a4)
S = lambda x: np.sin(2*np.pi*x)*np.sin(np.pi*x)**2
intS2 = quad(lambda x: S(x)**2, 0, 1, epsabs=1e-15, epsrel=1e-13)[0]
rw = sp.lambdify((M, mu, c3, c4), rewrite, 'numpy')


def soc(xv, delv):
    av = 2.0*xv
    Epn = lambda k: 2*k - 2*(2*k - delv/2)/np.sqrt((2*k - delv/2)**2 + av**2)
    k0 = brentq(Epn, -10, 10, xtol=1e-15, rtol=1e-15)
    D0 = 2*k0 - delv/2
    R0 = np.sqrt(D0**2 + av**2)
    return 1/(2 - 4*av**2/R0**3), 4*av**2*D0/R0**5, 2*av**2*(av**2 - 4*D0**2)/R0**7


registered = [("SOC g=0.2 x=1.20 d=0", 0.2, 1.2, 0.0, -8.19478898e-01),
              ("SOC g=0.2 x=1.10 d=0", 0.2, 1.1, 0.0, -2.56268465e+00),
              ("SOC g=0.2 x=1.20 d=0.05", 0.2, 1.2, 0.05, -7.14043305e-01),
              ("Galilean M=1 mu=1", None, None, None, -4.88281250e-03)]
ok4 = True
for lab, gam, xv, dv, reg in registered:
    if gam is None:
        val = rw(1.0, 1.0, 0.0, 0.0)*intS2
    else:
        Mv, c3v, c4v = soc(xv, dv)
        val = rw(Mv, gam, c3v, c4v)*intS2
    # registered values are printed with 9 significant digits -> compare to 1e-8 relative
    ok = abs(val/reg - 1) < 1e-8
    ok4 &= ok
    print(f"   {lab:26s}: rewrite {val:+.9e}  registered {reg:+.8e}  rel {abs(val/reg - 1):.1e}", flush=True)
check("R4 rewrite reproduces registered B1 formula values (rel 1e-8)", ok4)

print("\nALL CHECKS:", "PASS" if not FAIL else f"FAIL {FAIL}", flush=True)
sys.exit(0 if not FAIL else 1)
