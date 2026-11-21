"""
SOC_DETUNING_TOLERANCE_NOTE_v0.1, section 5: symbolic (Part A) and exact numeric (Part B) check,
with positive controls PC1-PC3. One run, python3 -I. No fitting, no free parameters.

Conventions (as in the note and in experiments/soc_bec/dipole_check_v0.1):
  lower band E_-(p) = p^2/2m + E_R - R(p),  R = sqrt(D^2 + a^2),  D = b p - delta/2,  a = Omega/2,  b = k_R/m.
  numeric units: E_R = 1, k_R = 1 (so m = 1/2, b = 2, a = 2x, x = Omega/4E_R, mu = g n = gamma).
  Amplitude level: eikonal background method of SOC_BEC_QUARTIC_NOTE_v0.3 section 7
  (right-mover on the left-moving constant-gradient background theta_t = c*eps, theta_x = eps).
"""
import sys
import numpy as np
import sympy as sp
from scipy.optimize import brentq

FAIL = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def is_zero(expr):
    return sp.simplify(sp.expand(expr)) == 0


# ---------------------------------------------------------------- Part A (symbolic)
print("=== Part A: symbolic ===", flush=True)
p, dl, a, b, m, ER, x = sp.symbols('p delta a b m E_R x', positive=True)
D = b*p - dl/2
R = sp.sqrt(D**2 + a**2)

# A1: derivatives of R
forms = {1: b*D/R, 2: b**2*a**2/R**3, 3: -3*b**3*a**2*D/R**5, 4: -3*b**4*a**2*(a**2 - 4*D**2)/R**7}
check("A1 derivatives of R (n=1..4, incl. a^2-4D^2 form)",
      all(is_zero(sp.diff(R, p, n) - f) for n, f in forms.items()))

# A2: minimum condition and leading-order d = D0/a ~ -delta/(4 E_R (x-1))
# substitute m b^2 = 2 E_R and a = 2 E_R x
bb = sp.sqrt(2*ER/m)
aa = 2*ER*x
p1 = sp.symbols('p1')
Ep = p/m - b*D/R                                   # E_-'(p)
Ep_sub = Ep.subs({b: bb, a: aa}).subs(p, p1*dl)
lin = sp.series(Ep_sub, dl, 0, 2).removeO().coeff(dl, 1)
p1_sol = sp.solve(sp.Eq(lin, 0), p1)
ok_a2 = len(p1_sol) == 1
if ok_a2:
    D0_lin = (bb*p1_sol[0]*dl - dl/2)
    d_lin = sp.simplify(D0_lin/aa)
    ok_a2 = is_zero(d_lin - (-dl/(4*ER*(x - 1))))
check("A2 minimum condition -> d = D0/a = -delta/(4 E_R (x-1)) + O(delta^3)", ok_a2,
      f"d_lin = {d_lin if len(p1_sol)==1 else p1_sol}")
# oddness of D0(delta): next correction is O(delta^3) -> check quadratic coefficient of p0(delta) is zero
p2 = sp.symbols('p2')
Ep_sub2 = Ep.subs({b: bb, a: aa}).subs(p, p1_sol[0]*dl + p2*dl**2)
quad = sp.series(Ep_sub2, dl, 0, 3).removeO().coeff(dl, 2)
check("A2b p0(delta) has no delta^2 term (D0 odd in delta)", is_zero(sp.solve(sp.Eq(quad, 0), p2)[0]))

# A3: O(d^2) forms in terms of d = D0/a (exact in d otherwise)
d = sp.symbols('d', real=True)
R0 = a*sp.sqrt(1 + d**2)
invM = 1/m - b**2*a**2/R0**3
c3 = b**3*a**2*(a*d)/(2*R0**5)
c4p = b**4*a**2*(a**2 - 4*a**2*d**2)/(8*R0**7)
invMstar = 1/m - b**2/a
# x-1 expressed through m b^2/a = 1/x  ->  x = a/(m b^2)
xx = a/(m*b**2)
s_invM = sp.series(invM, d, 0, 4).removeO()
s_c3 = sp.series(c3, d, 0, 3).removeO()
s_c4 = sp.series(c4p, d, 0, 4).removeO()
check("A3a 1/M = (1/m*)(1 + 3 d^2/(2(x-1))) + O(d^4)",
      is_zero(s_invM - invMstar*(1 + 3*d**2/(2*(xx - 1)))))
check("A3b c3 = b^3 d/(2 a^2) + O(d^3)", is_zero(s_c3 - b**3*d/(2*a**2)))
check("A3c c4' = c4 (1 - 15 d^2/2) + O(d^4), c4 = b^4/(8 a^3)",
      is_zero(s_c4 - b**4/(8*a**3)*(1 - sp.Rational(15, 2)*d**2)))

# A4/A5/A6: eikonal amplitude coefficient with a general non-even band
eps, mu, M = sp.symbols('epsilon mu M', positive=True)
C3, C4, C5 = sp.symbols('c3 c4 c5', real=True)
u = sp.symbols('u')
Et = u**2/(2*M) + C3*u**3 + C4*u**4 + C5*u**5
c = sp.sqrt(mu/M)
E0, E1, E2 = [sp.diff(Et, u, n).subs(u, eps) for n in range(3)]
tau = 6*C3*M*sp.sqrt(mu*M)
Q = 8*mu*M**2*C4

# right-mover on left-moving background theta_t = +c eps
cplus = E1 + sp.sqrt((mu - c*eps - E0)*E2)
ser = sp.series(cplus, eps, 0, 3).removeO()
s1, s2 = ser.coeff(eps, 1), ser.coeff(eps, 2)
check("A4a s1 = 1/(2M) + t/2, t = 6 c3 mu/c", is_zero(s1 - (1/(2*M) + 3*C3*mu/c)))
check("A4b c s2 = 6 c4 mu - 1/(4M^2) - (t - 1/M)^2/8",
      is_zero(c*s2 - (6*C4*mu - 1/(4*M**2) - (6*C3*mu/c - 1/M)**2/8)))
K2 = sp.series(1/cplus, eps, 0, 3).removeO().coeff(eps, 2)
target = sp.Rational(5, 8) + tau/4 + 3*tau**2/8 - sp.Rational(3, 4)*Q
check("A4c K2 c^3 M^2 = 5/8 + tau/4 + 3 tau^2/8 - (3/4) Q  (c5 absent)",
      is_zero(K2*c**3*M**2 - target))
K2_c30 = K2.subs(C3, 0)
check("A5 c3 = 0 reduces to registered B3: K2 = (1/c)(5/(8 M mu) - 6 c4 M)",
      is_zero(K2_c30 - (1/c)*(sp.Rational(5, 8)/(M*mu) - 6*C4*M)))
# A6: left-mover on right-moving background theta_t = -c eps, computed directly
cminus_mag = sp.sqrt((mu + c*eps - E0)*E2) - E1
K2L = sp.series(1/cminus_mag, eps, 0, 3).removeO().coeff(eps, 2)
check("A6 left-mover computed directly: K2L c^3 M^2 = 5/8 - tau/4 + 3 tau^2/8 - (3/4) Q",
      is_zero(K2L*c**3*M**2 - (sp.Rational(5, 8) - tau/4 + 3*tau**2/8 - sp.Rational(3, 4)*Q)))


# ---------------------------------------------------------------- numeric pipeline (exact in delta)
def band_coeffs(xv, delv):
    """Exact local coefficients at the minimum; units E_R = k_R = 1 (m = 1/2, b = 2, a = 2x)."""
    av = 2.0*xv
    Ep_n = lambda k: 2*k - 2*(2*k - delv/2)/np.sqrt((2*k - delv/2)**2 + av**2)
    k0 = brentq(Ep_n, -10, 10, xtol=1e-15, rtol=1e-15)
    D0 = 2*k0 - delv/2
    R0 = np.sqrt(D0**2 + av**2)
    invM = 2 - 4*av**2/R0**3
    c3n = 4*av**2*D0/R0**5
    c4n = 2*av**2*(av**2 - 4*D0**2)/R0**7
    return k0, D0, av, invM, c3n, c4n


def margins(xv, gam, delv):
    k0, D0, av, invM, c3n, c4n = band_coeffs(xv, delv)
    Mv = 1/invM
    Qv = 8*gam*Mv**2*c4n
    tv = 6*c3n*Mv*np.sqrt(gam*Mv)
    G = Qv - 1                                             # contact: Lambda > 0 iff G > 0
    F = {s: 0.75*Qv - 0.625 - s*tv/4 - 3*tv**2/8 for s in (+1, -1)}   # advance iff F > 0
    return G, F, Qv, tv, D0, av


# uniqueness of the minimum: E'' >= 2 - 2/x > 0 for x > 1; also verified on a grid
def unique_min_on_grid(xv, delv):
    av = 2.0*xv
    k = np.linspace(-10, 10, 200001)
    Ep_n = 2*k - 2*(2*k - delv/2)/np.sqrt((2*k - delv/2)**2 + av**2)
    return int(np.sum(np.diff(np.sign(Ep_n)) != 0)) == 1


# ---------------------------------------------------------------- positive controls
print("\n=== Positive controls ===", flush=True)
reg_contact = {0.1: 1.204, 0.2: 1.280, 0.5: 1.420}
reg_amp = {0.1: 1.222, 0.2: 1.303, 0.5: 1.454}
ok1 = True
for gam in (0.1, 0.2, 0.5):
    xc = brentq(lambda xv: margins(xv, gam, 0.0)[0], 1.0001, 3.0, xtol=1e-12)
    xa = brentq(lambda xv: margins(xv, gam, 0.0)[1][+1], 1.0001, 3.0, xtol=1e-12)
    okg = abs(xc - reg_contact[gam]) < 1e-3 and abs(xa - reg_amp[gam]) < 1e-3
    ok1 &= okg
    print(f"  gamma={gam}: contact x* = {xc:.5f} (reg {reg_contact[gam]}), "
          f"amplitude x* = {xa:.5f} (reg {reg_amp[gam]})", flush=True)
check("PC1 delta = 0 thresholds reproduce v0.3 to 1e-3", ok1)

k0, D0, av, invM, c3n, c4n = band_coeffs(1.2, 1.0)        # Omega = 4.8 E_R, delta = 1 E_R
reg = dict(k0=-0.558996, e2=1.049877, e3=-1.100962, e4=-0.765348)
got = dict(k0=k0, e2=invM, e3=6*c3n, e4=24*c4n)
ok2 = all(abs(got[k] - reg[k]) < 1e-4 for k in reg)
print("  " + ", ".join(f"{k}: {got[k]:.6f} (reg {reg[k]})" for k in reg), flush=True)
check("PC2 Omega=4.8, delta=1 reproduces registered dipole check (stdout 3584bd7f) to 1e-4", ok2)

ok3 = True
gam, delv = 0.2, 1e-3
for xv in (1.1, 1.2, 1.3, 1.4):
    k0, D0, av, invM, c3n, c4n = band_coeffs(xv, delv)
    dex = D0/av
    invMstar = 2*(1 - 1/xv)
    c4_0 = 16/(8*av**3)
    approx = dict(d=-delv/(4*(xv - 1)), invM=invMstar*(1 + 3*dex**2/(2*(xv - 1))),
                  c3=8*dex/(2*av**2), c4=c4_0*(1 - 7.5*dex**2),
                  tau=-(3/(4*np.sqrt(2)))*delv*np.sqrt(gam/xv)*(xv - 1)**-2.5)
    Mv = 1/invM
    exact = dict(d=dex, invM=invM, c3=c3n, c4=c4n, tau=6*c3n*Mv*np.sqrt(gam*Mv))
    rel = {k: abs(exact[k]/approx[k] - 1) for k in exact}
    ok3 &= all(v < 1e-3 for v in rel.values())
    print(f"  x={xv}: " + ", ".join(f"{k} rel.err {v:.1e}" for k, v in rel.items()), flush=True)
check("PC3 exact values agree with section 2-3 O(d^2) forms at delta = 1e-3 (rel < 1e-3)", ok3)

if FAIL:
    print(f"\nFAILED: {FAIL}\nPart B is not interpreted (pre-registered rule).", flush=True)
    sys.exit(1)

# ---------------------------------------------------------------- Part B (exact numeric scan)
print("\n=== Part B: critical detuning (E_R units), exact in delta, delta in [0, 0.5] ===", flush=True)
dgrid = np.linspace(0.0, 0.5, 2001)


def critical(fun):
    """first delta in (0, 0.5] where fun changes sign from > 0; also count sign changes."""
    v = np.array([fun(dv) for dv in dgrid])
    if v[0] <= 0:
        return "no window at delta=0", 0
    sc = int(np.sum(np.diff(np.sign(v)) != 0))
    idx = np.where(v <= 0)[0]
    if len(idx) == 0:
        return "none <= 0.5", sc
    i = idx[0]
    return brentq(fun, dgrid[i - 1], dgrid[i], xtol=1e-13), sc


def fmt(r):
    v, sc = r
    s = f"{v:.5f}" if isinstance(v, float) else v
    return s + (f" [{sc} sign changes]" if sc > 1 else "")


uniq_ok = True
xs = np.round(np.arange(1.02, 1.5001, 0.02), 2)
results = {}
for gam in (0.1, 0.2, 0.5):
    print(f"\n gamma = {gam}", flush=True)
    print("   x     | delta_c contact | delta_c amp s=+1 | delta_c amp s=-1 | min amp  | |D0|=a/2 at", flush=True)
    for xv in xs:
        uniq_ok &= unique_min_on_grid(xv, 0.5)
        rc = critical(lambda dv: margins(xv, gam, dv)[0])
        rp = critical(lambda dv: margins(xv, gam, dv)[1][+1])
        rm = critical(lambda dv: margins(xv, gam, dv)[1][-1])
        mins = [r[0] for r in (rp, rm) if isinstance(r[0], float)]
        mn = f"{min(mins):.5f}" if mins else "-"
        dq = brentq(lambda dv: abs(margins(xv, gam, dv)[4]) - margins(xv, gam, dv)[5]/2, 0.0, 50.0, xtol=1e-12)
        results[(gam, xv)] = (rc, rp, rm, dq)
        print(f"   {xv:.2f}  | {fmt(rc):15s} | {fmt(rp):16s} | {fmt(rm):16s} | {mn:8s} | {dq:.5f}", flush=True)
check("B uniqueness of the minimum on a grid (delta = 0.5, all x)", uniq_ok)

# ---------------------------------------------------------------- pre-registered category
rc, rp, rm, dq = results[(0.2, 1.2)]
vals = [r[0] for r in (rp, rm) if isinstance(r[0], float)]
print("\n=== Pre-registered category (gamma = 0.2, x = 1.2, amplitude level, worse direction) ===", flush=True)
if not vals:
    print("delta_tol > 0.5 E_R -> robust", flush=True)
else:
    dt = min(vals)
    cat = "robust" if dt >= 0.1 else ("sensitive" if dt >= 0.01 else "fragile")
    print(f"delta_tol = {dt:.5f} E_R  ->  {cat}", flush=True)
print(f"contact-level critical delta at the same point: {fmt(rc)}", flush=True)
print("\nALL CHECKS:", "PASS" if not FAIL else f"FAIL {FAIL}", flush=True)
sys.exit(0 if not FAIL else 1)
