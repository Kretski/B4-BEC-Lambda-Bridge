"""T7a step 3 (DRAFT, unregistered): independent SymPy check of the Proca tree-level matching.
Arbiter = registered action and conventions (T3_FUND_PROCA_SPEC.md section 1, metric (-,+)); SymPy only checks the derivation."""
import sympy as sp
from sympy.calculus.euler import euler_equations
t, x = sp.symbols('t x', real=True)
M, gphi, gpsi, Lam = sp.symbols('M g_phi g_psi Lambda', positive=True)
f, g, A0, A1 = [sp.Function(n)(t, x) for n in ('f', 'g', 'A0', 'A1')]
ft, fx, gt, gx = f.diff(t), f.diff(x), g.diff(t), g.diff(x)
F = A1.diff(t) - A0.diff(x)                       # F = F_01 = d_t A1 - d_x A0, lower indices, metric (-,+)
J = gphi*ft*fx + gpsi*gt*gx
L = ft**2/2 - fx**2/2 - Lam*fx**4/4 + gt**2/2 - gx**2/2 + F**2/2 + M**2*(A0**2 - A1**2)/2 + A1*J
res = {}
# (0) convention check: -1/4 F_{mu nu}F^{mu nu} with eta=diag(-1,1) equals 1/2 F_01^2 ; -1/2 M^2 A_mu A^mu = 1/2 M^2 (A0^2 - A1^2)
eta = sp.diag(-1, 1); Fmn = sp.Matrix([[0, F], [-F, 0]]); Fup = eta*Fmn*eta
res['conv_F2'] = sp.simplify(-sp.Rational(1, 4)*sum(Fmn[i, j]*Fup[i, j] for i in range(2) for j in range(2)) - F**2/2)
Av = sp.Matrix([A0, A1]); res['conv_mass'] = sp.simplify(-M**2/2*(Av.T*eta*Av)[0] - M**2*(A0**2 - A1**2)/2)
# (1) EL equations from the full action
eqs = euler_equations(L, [f, g, A0, A1], [t, x])
EA0 = eqs[2].lhs - eqs[2].rhs; EA1 = eqs[3].lhs - eqs[3].rhs
# Gauss: A0 eq must equal F_x + M^2 A0 = 0 (up to sign)
res['gauss_form'] = sp.simplify(EA0 - (F.diff(x) + M**2*A0))
# A1 eq: -F_t - M^2 A1 + J = 0 (up to sign)
res['A1_form'] = sp.simplify(EA1 - (-F.diff(t) - M**2*A1 + J))
# (2) eliminate A0 and F in Fourier space (plane wave exp(i(kx - w t))): check (D + ...) identity
w, k = sp.symbols('omega k', real=True)
a0, a1, j = sp.symbols('a0 a1 j')   # amplitudes
dt, dx = -sp.I*w, sp.I*k
Fa = dt*a1 - dx*a0
a0_sol = sp.solve(sp.Eq(dx*Fa + M**2*a0, 0), a0)[0]      # Gauss
Fa_sol = sp.simplify(Fa.subs(a0, a0_sol))
res['F_in_A1'] = sp.simplify(Fa_sol - M**2/(M**2 + k**2)*dt*a1)   # F = M^2 (M^2 - d_x^2)^-1 d_t A1
a1_eq = sp.simplify((-dt*Fa_sol - M**2*a1 + j).subs(a0, a0_sol))   # A1 equation after eliminating A0
a1_sol = sp.simplify(sp.solve(sp.Eq(a1_eq, 0), a1)[0]/j)          # A1 = K j
K = sp.simplify(a1_sol)
claim = (1 + k**2/M**2)/(M**2 + k**2 - w**2)                       # ((M^2-d_x^2)/M^2)/(d_t^2 + M^2 - d_x^2)
res['K_claim'] = sp.simplify(K - claim)
# (3) series in 1/M^2 (eps = 1/M^2)
eps = sp.symbols('eps', positive=True)
Ks = sp.series((claim*M**2).subs(M, 1/sp.sqrt(eps)), eps, 0, 3).removeO()
print('K*M^2 series in 1/M^2:', sp.expand(Ks))
# (4) kappa sign from leading order
O1, O2, O3 = sp.symbols('O1 O2 O3')
p, s, kap = sp.symbols('p s kappa')
Jf = gphi*sp.Symbol('Ofx') + gpsi*sp.Symbol('Ogx')   # placeholders for f_t f_x, g_t g_x
Lead = sp.expand(Jf**2/(2*M**2))
cO1 = Lead.coeff(sp.Symbol('Ofx'), 2); cO3 = Lead.coeff(sp.Symbol('Ogx'), 2)
cO2 = Lead.coeff(sp.Symbol('Ofx'), 1).coeff(sp.Symbol('Ogx'), 1)
res['p'] = sp.simplify(2*cO1 - gphi**2/M**2); res['s'] = sp.simplify(2*cO3 - gpsi**2/M**2)
res['kappa_from_cross'] = sp.simplify(-cO2 + gphi*gpsi/M**2*(-1)*(-1) - (gphi*gpsi/M**2 - gphi*gpsi/M**2))
kappa = sp.simplify(-cO2)
print('cross coefficient of O2 =', sp.simplify(cO2), ' => kappa = -cross =', kappa)
res['kappa_is_minus_gg_over_M2'] = sp.simplify(kappa + gphi*gpsi/M**2)
res['kappa2_eq_ps'] = sp.simplify(kappa**2 - (gphi**2/M**2)*(gpsi**2/M**2))
for n, v in res.items(): print(n, '->', v)
print('ALL_ZERO' if all(sp.simplify(v) == 0 for v in res.values()) else 'NONZERO_FOUND')

# (5) on-shell action: quadratic form after substituting Gauss solution; stationary value with source equals 1/2 K |j|^2
a1c, a0c = sp.symbols('a1c a0c')                   # complex conjugate amplitudes
dtc, dxc = sp.I*w, -sp.I*k
Fa_c = dtc*a1c - dxc*a0c
a0c_sol = sp.solve(sp.Eq(dxc*Fa_c + M**2*a0c, 0), a0c)[0]
Q = sp.Rational(1, 2)*(Fa_c.subs(a0c, a0c_sol)*Fa.subs(a0, a0_sol) + M**2*(a0c_sol*a0_sol - a1c*a1) )
c = sp.simplify(sp.diff(Q, a1c, a1))*2/1            # S = (1/2) c a1c a1
c = sp.simplify(c)
print('quadratic coefficient c =', c)
print('check K = -1/c :', sp.simplify(K + 1/c))
print('check on-shell: stationary value factor -1/(2c) vs K/2 :', sp.simplify(-1/(2*c) - K/2))
