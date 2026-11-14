"""Numeric parts (3)-(5) of the Proca verification (see t3fund_proca_check.py for (1)-(2), symbolic)."""
import numpy as np, sympy as sp, mpmath as mp, json
rng = np.random.default_rng(20261114)
out = {}
# ---------------- (3) reduced system vs Euler-Lagrange at random smooth states ----------------
t, x = sp.symbols('t x', real=True)
Mv, gphi_v, gpsi_v, Lam_v = 1.3, 0.7, -0.9, 0.6
f = sp.Function('f')(t, x); g = sp.Function('g')(t, x); a0 = sp.Function('a0')(t, x); a1 = sp.Function('a1')(t, x)
half = sp.Rational(1, 2)
L = (half*f.diff(t)**2 - half*f.diff(x)**2 - sp.Float(Lam_v)/4*f.diff(x)**4 + half*g.diff(t)**2 - half*g.diff(x)**2
     + half*(a1.diff(t) - a0.diff(x))**2 + half*Mv**2*(a0**2 - a1**2) + a1*(gphi_v*f.diff(t)*f.diff(x) + gpsi_v*g.diff(t)*g.diff(x)))
EL = sp.calculus.euler.euler_equations(L, [f, g, a0, a1], [t, x])
def trig(nm, seed):
    rr = np.random.default_rng(seed); e = 0
    for j in range(3):
        m = int(rr.integers(1, 4)); n = float(rr.uniform(-2, 2)); amp = float(rr.uniform(-.5, .5)); ph = float(rr.uniform(0, 6))
        e += amp*sp.cos(m*x - n*t + ph)
    return e
def mode_list(seed):
    rr = np.random.default_rng(seed); out = []
    for j in range(3):
        m = int(rr.integers(1, 4)); n = float(rr.uniform(-2, 2)); amp = float(rr.uniform(-.5, .5)); ph = float(rr.uniform(0, 6)); out.append((m, n, amp, ph))
    return out
fm, gm, am = mode_list(1), mode_list(2), mode_list(3)
fx_ = sum(a*sp.cos(m*x - n*t + p) for m, n, a, p in fm) + 0.8*t + 0.4*x      # include a gradient background
gx_ = sum(a*sp.cos(m*x - n*t + p) for m, n, a, p in gm) - 0.5*t + 0.9*x
A1x = sum(a*sp.cos(m*x - n*t + p) for m, n, a, p in am) + 0.3
# Gauss constraint (M^2 - d_xx) A0 = -A1_tx, solved mode by mode: for cos(m x - n t + p): A1_tx = n m cos(...)?  d_t d_x cos(theta)= (-(-n))... compute symbolically
A0x = 0
for m, n, a, p in am:
    th = m*x - n*t + p
    A1tx = sp.diff(a*sp.cos(th), t, x)         # = -a*n*m*cos(th)?? (sympy decides)
    # for A0 = beta*cos(th):  (M^2 + m^2) beta cos(th) = -A1tx  -> beta cos(th) = -A1tx/(M^2+m^2)
    A0x += -A1tx/(Mv**2 + m**2)
sub = {f: fx_, g: gx_, a0: A0x, a1: A1x}
Fv = A1x.diff(t) - A0x.diff(x)
res3 = {}
pts = [(0.3, 0.7), (1.1, -0.4), (2.0, 1.9)]
def num(expr, pt): return float(sp.N(expr.subs({t: pt[0], x: pt[1]})))
ELn = [e.lhs.subs(sub).doit() - e.rhs.subs(sub).doit() for e in EL]
a_, b_ = fx_.diff(t), fx_.diff(x); c_, d_ = gx_.diff(t), gx_.diff(x)
A1t_red = Fv - Fv.diff(x, 2)/Mv**2
Ja = gphi_v*a_*b_ + gpsi_v*c_*d_
rhs_a = (1 + 3*Lam_v*b_**2)*b_.diff(x) - gphi_v*(A1t_red*b_ + A1x.diff(x)*a_ + 2*A1x*a_.diff(x))
rhs_c = d_.diff(x) - gpsi_v*(A1t_red*d_ + A1x.diff(x)*c_ + 2*A1x*c_.diff(x))
rhs_F = Ja - Mv**2*A1x
chk = dict(
    gauss_EL_a0=max(abs(num(ELn[2], p)) for p in pts),
    A1t_identity=max(abs(num(A1x.diff(t) - A1t_red, p)) for p in pts),
    EL_f=max(abs(num(ELn[0] + (fx_.diff(t, 2) - rhs_a), p)) for p in pts),
    EL_g=max(abs(num(ELn[1] + (gx_.diff(t, 2) - rhs_c), p)) for p in pts),
    EL_A1=max(abs(num(ELn[3] + (Fv.diff(t) - rhs_F), p)) for p in pts))
# sympy euler_equations returns dL/dq - d_t dL/dq_t - d_x dL/dq_x, i.e. minus the evolution form: EL = -(q_tt - RHS)
print('(3) EL vs reduced system (max abs residual over 3 points):', chk)
out['reduced_vs_EL'] = chk
# ---------------- (4) energy conservation of the reduced system (spectral) ----------------
Lbox = 2*np.pi*4; N = 128; xg = np.arange(N)*Lbox/N; kg = 2*np.pi*np.fft.fftfreq(N, d=Lbox/N)
D = lambda u: np.real(np.fft.ifft(1j*kg*np.fft.fft(u)))
def rhs(S, g_phi, g_psi, M, Lam):
    a, b, c, d, A, F = S
    Fxx = D(D(F)); A1t = F - Fxx/M**2
    da = (1 + 3*Lam*b**2)*D(b) - g_phi*(A1t*b + D(A)*a + 2*A*D(a))
    dc = D(d) - g_psi*(A1t*d + D(A)*c + 2*A*D(c))
    db = D(a); dd = D(c)
    dA = A1t
    dF = g_phi*a*b + g_psi*c*d - M**2*A
    return np.array([da, db, dc, dd, dA, dF])
def energy(S, M, Lam):
    a, b, c, d, A, F = S
    return Lbox/N*np.sum(0.5*a**2 + 0.5*b**2 + Lam/4*b**4 + 0.5*c**2 + 0.5*d**2 + 0.5*F**2 + D(F)**2/(2*M**2) + 0.5*M**2*A**2)
S0 = np.array([0.3*np.cos(xg*(j+1)/4 + j) + 0.1*np.sin(2*xg + j) + (0.5 if j < 4 else 0.0) for j in range(6)])
# make a,b consistent smooth fields (b=f_x of some f not required: reduced system conserves E for any (a,b,c,d,A,F) with db=D(a) only if b,d are gradients; use a=f_t,b=f_x from f)
fgrid = 0.4*np.sin(xg/4) + 0.2*np.cos(xg/2) + 0.3*np.sin(xg); ggrid = 0.3*np.cos(xg/4 + 1) + 0.2*np.sin(xg*3/4)
S0[1] = D(fgrid) + 0.4; S0[3] = D(ggrid) + 0.9
Mv2, gp, gs, Lm = 1.3, 0.7, -0.9, 0.6
R = rhs(S0, gp, gs, Mv2, Lm)
def dEdt(h):
    return (energy(S0 + h*R, Mv2, Lm) - energy(S0 - h*R, Mv2, Lm))/(2*h)
h = 1e-3; Dh = (4*dEdt(h/2) - dEdt(h))/3; Dh2 = (4*dEdt(h/4) - dEdt(h/2))/3; Dr = (16*Dh2 - Dh)/15
scale = abs(energy(S0, Mv2, Lm))
print('(4) dE/dt (Richardson) =', Dr, ' relative to E:', Dr/scale)
out['dEdt_relative'] = float(Dr/scale)
# ---------------- (5) linearised matrix vs Q (50 digits) ----------------
mp.mp.dps = 50
def Qpoly(p, s, kap, Lam, w, q, v, r, eps):
    c = mp.mpf  # placeholder
    return None
def lin_matrix(g_phi, g_psi, Lam, w, q, v, r, kk, M=1):
    mp.mp.dps = 50
    I = mp.mpc(0, 1); kh = mp.sqrt(M**2 + kk**2); Abar = (g_phi*w*q + g_psi*v*r)/M**2
    Mx = mp.zeros(6, 6)   # order: da, db, dc, dd, At, Ft  (M=1)
    Mx[0, 1] = I*kk*(1 + 3*Lam*q**2); Mx[0, 0] = -2*I*kk*g_phi*Abar; Mx[0, 5] = -g_phi*q*kh; Mx[0, 4] = -g_phi*I*kk*w
    Mx[1, 0] = I*kk
    Mx[2, 3] = I*kk; Mx[2, 2] = -2*I*kk*g_psi*Abar; Mx[2, 5] = -g_psi*r*kh; Mx[2, 4] = -g_psi*I*kk*v
    Mx[3, 2] = I*kk
    Mx[4, 5] = kh
    Mx[5, 0] = kh*g_phi*q; Mx[5, 1] = kh*g_phi*w; Mx[5, 2] = kh*g_psi*r; Mx[5, 3] = kh*g_psi*v; Mx[5, 4] = -kh
    return Mx
def Q_roots(g_phi, g_psi, Lam, w, q, v, r, kk, M=1):
    mp.mp.dps = 50
    p, s, kap = g_phi**2/M**2, g_psi**2/M**2, -g_phi*g_psi/M**2
    aphi = p*w*q - kap*v*r; apsi = -kap*w*q + s*v*r; xx = mp.mpf(kk)/M; eps = xx**2/(1 + xx**2)
    c = sp.Symbol('c')
    d1 = c**2 - 2*aphi*c - (1 + 3*Lam*q**2); d2 = c**2 - 2*apsi*c - 1
    Fc = sp.expand(d1*d2 + p*(w - q*c)**2*d2 + s*(v - r*c)**2*d1)
    Qs = sp.Poly(sp.expand(Fc - sp.Float(str(eps), 50)*c**2*d1*d2), c)
    co = [mp.mpf(str(z)) for z in Qs.all_coeffs()]
    return mp.polyroots(co, maxsteps=500, extraprec=300)
worst = 0; worst_xinf = 0
for trial in range(12):
    gp_, gs_ = mp.mpf(float(rng.uniform(-1, 1))), mp.mpf(float(rng.uniform(-1, 1))); Lm_ = mp.mpf(float(rng.uniform(0, 1)))
    w_, q_, v_, r_ = [mp.mpf(float(z)) for z in rng.uniform(-1.2, 1.2, 4)]
    for kk in (0.05, 0.7, 3.0, 40.0):
        Mx = lin_matrix(gp_, gs_, Lm_, w_, q_, v_, r_, mp.mpf(kk))
        ev = mp.eig(Mx, left=False, right=False)
        om = sorted([(mp.mpc(0, 1)*e) for e in ev], key=lambda z: (mp.re(z), mp.im(z)))         # omega = i*lambda  (U ~ e^{-i omega t}, d_t = -i omega)
        rts = sorted([z*kk for z in Q_roots(gp_, gs_, Lm_, w_, q_, v_, r_, kk)], key=lambda z: (mp.re(z), mp.im(z)))
        worst = max(worst, max(abs(a - b) for a, b in zip(om, rts)))
print('(5) max |eig(i*M_lin) - k*roots(Q)| over 48 (background,k) pairs, 50 digits:', mp.nstr(worst, 5))
out['lin_matrix_vs_Q'] = mp.nstr(worst, 5)
json.dump(out, open('t3fund_proca_check_numeric.json', 'w'), indent=1)
