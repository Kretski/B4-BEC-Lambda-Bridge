"""T6: symbolic checks of the ZAMO and static observer congruences in Kerr (Boyer-Lindquist, signature -+++, G = c = 1)."""
import sympy as sp
t, r, th, ph, M, a = sp.symbols('t r theta phi M a', real=True)
Sig = r**2 + a**2*sp.cos(th)**2
Del = r**2 - 2*M*r + a**2
A = (r**2 + a**2)**2 - a**2*Del*sp.sin(th)**2
gtt = -(1 - 2*M*r/Sig); gtp = -2*M*a*r*sp.sin(th)**2/Sig; gpp = A*sp.sin(th)**2/Sig
# inverse of the (t, phi) block
det = sp.simplify(gtt*gpp - gtp**2)
print("det of (t,phi) block + Delta sin^2:", sp.simplify(det + Del*sp.sin(th)**2))          # expect 0
gTT = sp.simplify(gpp/det)                                                               # g^{tt}
alpha2 = sp.simplify(-1/gTT)
print("alpha^2 - Delta*Sigma/A:", sp.simplify(alpha2 - Del*Sig/A))                         # expect 0
omega = sp.simplify(-gtp/gpp)
print("omega - 2 M a r / A:", sp.simplify(omega - 2*M*a*r/A))                              # expect 0
# ZAMO: u_mu = -alpha dt  =>  u^mu = (1/alpha)(1,0,0,omega); norm
al = sp.sqrt(alpha2)
ut, up = 1/al, omega/al
print("ZAMO norm + 1:", sp.simplify(gtt*ut**2 + 2*gtp*ut*up + gpp*up**2 + 1))               # expect 0
print("ZAMO u_phi (angular momentum):", sp.simplify(gtp*ut + gpp*up))                      # expect 0
# static observer: xi = d/dt; hypersurface orthogonal iff g_tphi/g_tt is constant
ratio = sp.simplify(gtp/gtt)
print("d/dr (g_tphi/g_tt):", sp.factor(sp.simplify(sp.diff(ratio, r))))                    # nonzero for a != 0
print("same at a = 0:", sp.simplify(sp.diff(ratio, r).subs(a, 0)))                         # 0
# domains: g_tt < 0 (static timelike) vs Delta > 0 (ZAMO timelike)
print("g_tt = 0 at r =", sp.solve(sp.Eq(sp.numer(sp.together(gtt)), 0), r))
print("Delta = 0 at r =", sp.solve(sp.Eq(Del, 0), r))
# local frequency seen by ZAMO for a mode with Killing energy E and axial angular momentum L: omega_loc = -u^mu k_mu, k_t = -E, k_phi = L
E, L = sp.symbols('E L', real=True)
print("omega_loc * alpha - (E - omega L):", sp.simplify((-(ut*(-E) + up*L))*al - (E - omega*L)))   # expect 0
