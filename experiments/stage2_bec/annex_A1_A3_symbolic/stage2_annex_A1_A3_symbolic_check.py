"""Symbolic check of STAGE2_PREREGISTRATION_PROCEDURAL_ANNEX_v0.1, items A1 and A3 only.

Authorised scope (review of 2026-11-19): A1 (normalised 4x4 energy form from the quadratic GP
functional) and A3 (E(k) - E(0) positive semidefinite, hence lambda_min non-decreasing in k^2,
minimum at k = 0; one-component limit V^2 < a). No numerical Stage 2 work, no search, no fit.

Pure SymPy, no input files, no randomness. Exit code 0 iff all checks pass.
"""
import sys
import sympy as sp

results = []


def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))


# ---------------------------------------------------------------- symbols
hb, k, L = sp.symbols('hbar k L', positive=True)
x = sp.symbols('x', real=True)
eps = sp.symbols('epsilon', real=True)
m, n0, g, V, mu = {}, {}, {}, {}, {}
for i in (1, 2):
    m[i] = sp.symbols(f'm{i}', positive=True)
    n0[i] = sp.symbols(f'n{i}', positive=True)
    g[i] = sp.symbols(f'g{i}{i}', positive=True)   # g_ii > 0 (validated in bec_dispersion.BECParams)
    V[i] = sp.symbols(f'V{i}', real=True)
    mu[i] = sp.symbols(f'mu{i}', real=True)
g12 = sp.symbols('g12', real=True)

# fluctuation amplitudes of one real Fourier mode: cos and sin parts
Nc, Ns, Uc, Us = {}, {}, {}, {}
for i in (1, 2):
    Nc[i], Ns[i], Uc[i], Us[i] = sp.symbols(f'Nc{i} Ns{i} Uc{i} Us{i}', real=True)

# ---------------------------------------------------------------- A1: GP functional, hydrodynamic form
# Psi_i = sqrt(n_i) exp(i theta_i); u_i = (hbar/m_i) d_x theta_i.
# (hbar^2/2m)|d_x Psi|^2 = hbar^2 (d_x n)^2 / (8 m n) + (m/2) n u^2.
# Energy density minus mu_i n_i (grand-canonical, lab frame):
#   e = sum_i [ hbar^2 (n_i')^2/(8 m_i n_i) + (m_i/2) n_i u_i^2 + g_ii n_i^2/2 - mu_i n_i ] + g12 n1 n2
n, u = {}, {}
for i in (1, 2):
    n[i] = n0[i] + eps*(Nc[i]*sp.cos(k*x) + Ns[i]*sp.sin(k*x))
    u[i] = V[i] + eps*(Uc[i]*sp.cos(k*x) + Us[i]*sp.sin(k*x))

e = g12*n[1]*n[2]
for i in (1, 2):
    e += (hb**2*sp.diff(n[i], x)**2/(8*m[i]*n[i]) + m[i]/2*n[i]*u[i]**2
          + g[i]*n[i]**2/2 - mu[i]*n[i])

# expand in eps to second order, then average over one period 2*pi/k
ser = sp.series(e, eps, 0, 3).removeO()
e1 = sp.expand(ser.coeff(eps, 1))
e2 = sp.expand(ser.coeff(eps, 2))
P = 2*sp.pi/k


def avg(expr):
    return sp.simplify(sp.integrate(sp.expand(expr), (x, 0, P))/P)


# first order: the mode average of a single cos/sin is zero, so the linear term vanishes for k != 0
# irrespective of mu; check it, and record the uniform (k = 0) stationarity condition separately.
e1_avg = avg(e1)
check("A1.0  first-order term averages to zero for a k != 0 mode", sp.simplify(e1_avg) == 0, str(e1_avg))
mu_stat = {1: g[1]*n0[1] + g12*n0[2] + m[1]*V[1]**2/2,
           2: g[2]*n0[2] + g12*n0[1] + m[2]*V[2]**2/2}
for i in (1, 2):
    d_uniform = sp.diff(e.subs(eps, 0).subs({sp.cos(k*x): 1, sp.sin(k*x): 0}), n0[i])
    check(f"A1.0' uniform stationarity in n{i} gives mu{i} = g{i}{i} n{i} + g12 n_j + m{i} V{i}^2/2",
          sp.simplify(d_uniform.subs(mu[i], mu_stat[i])) == 0)

e2_avg = avg(e2)

# claimed form (annex A1): E2 = sum_i [ 1/2 m n du^2 + m V dn du + 1/2 (g + hbar^2 k^2/(4 m n)) dn^2 ] + g12 dn1 dn2
# for a mode with cos and sin parts the period average is 1/2 * [form(c) + form(s)]


def claimed(Na1, Ua1, Na2, Ua2):
    out = g12*Na1*Na2
    for i, (Ni, Ui) in zip((1, 2), ((Na1, Ua1), (Na2, Ua2))):
        out += (m[i]*n0[i]*Ui**2/2 + m[i]*V[i]*Ni*Ui
                + (g[i] + hb**2*k**2/(4*m[i]*n0[i]))*Ni**2/2)
    return out


claimed_avg = sp.Rational(1, 2)*(claimed(Nc[1], Uc[1], Nc[2], Uc[2]) + claimed(Ns[1], Us[1], Ns[2], Us[2]))
check("A1.1  second-order period average equals the annex E2 form (cos and sin parts decouple)",
      sp.simplify(sp.expand(e2_avg - claimed_avg)) == 0)

# matrix of the form, E2 = 1/2 z^T M z, z = (dn1, du1, dn2, du2)
z = sp.Matrix([Nc[1], Uc[1], Nc[2], Uc[2]])
form_c = claimed(Nc[1], Uc[1], Nc[2], Uc[2])
M = sp.hessian(form_c, z)
M_expected = sp.Matrix([
    [g[1] + hb**2*k**2/(4*m[1]*n0[1]), m[1]*V[1], g12, 0],
    [m[1]*V[1], m[1]*n0[1], 0, 0],
    [g12, 0, g[2] + hb**2*k**2/(4*m[2]*n0[2]), m[2]*V[2]],
    [0, 0, m[2]*V[2], m[2]*n0[2]]])
check("A1.2  matrix M of the form (Hessian) as stated", sp.simplify(M - M_expected) == sp.zeros(4, 4))
check("A1.3  M is symmetric", sp.simplify(M - M.T) == sp.zeros(4, 4))

# congruence D M D with D = diag(g11^-1/2, (m1 n1)^-1/2, g22^-1/2, (m2 n2)^-1/2)
D = sp.diag(1/sp.sqrt(g[1]), 1/sp.sqrt(m[1]*n0[1]), 1/sp.sqrt(g[2]), 1/sp.sqrt(m[2]*n0[2]))
Et = sp.simplify(D*M*D)
a = {i: g[i]*n0[i]/m[i] for i in (1, 2)}
b = {i: hb**2/(4*m[i]**2) for i in (1, 2)}
f = g12/sp.sqrt(g[1]*g[2])
Et_expected = sp.Matrix([
    [1 + b[1]*k**2/a[1], V[1]/sp.sqrt(a[1]), f, 0],
    [V[1]/sp.sqrt(a[1]), 1, 0, 0],
    [f, 0, 1 + b[2]*k**2/a[2], V[2]/sp.sqrt(a[2])],
    [0, 0, V[2]/sp.sqrt(a[2]), 1]])
check("A1.4  D M D equals the annex normalised matrix E~(k) in a_i, b_i, f, V_i",
      sp.simplify(Et - Et_expected) == sp.zeros(4, 4))
check("A1.5  D is invertible (all diagonal entries finite and nonzero for positive g_ii, m_i, n_i)",
      all(sp.simplify(D[j, j]).is_positive for j in range(4)))
# note: Sylvester -> same inertia (signs) of M and E~; NOT the same eigenvalue magnitudes (recorded in report)

# a_i, b_i positive under the stated physical conditions
check("A1.6  a_i = g_ii n_i / m_i > 0 and b_i = hbar^2/(4 m_i^2) > 0 for g_ii, n_i, m_i, hbar > 0",
      all(a[i].is_positive and b[i].is_positive for i in (1, 2)))

# ---------------------------------------------------------------- A3: monotonicity in k^2
Pdiff = sp.simplify(Et - Et.subs(k, 0))
P_expected = sp.diag(b[1]*k**2/a[1], 0, b[2]*k**2/a[2], 0)
check("A3.1  E~(k) - E~(0) = diag(b1 k^2/a1, 0, b2 k^2/a2, 0)", sp.simplify(Pdiff - P_expected) == sp.zeros(4, 4))
check("A3.2  that difference is positive semidefinite (diagonal, entries >= 0)",
      all((sp.simplify(Pdiff[j, j]).is_nonnegative) for j in range(4)))
# monotone in k^2: for 0 <= k1 < k2, E~(k2) - E~(k1) = diag(b1 (k2^2 - k1^2)/a1, 0, ...) >= 0
k1, k2 = sp.symbols('k1 k2', positive=True)
dP = sp.simplify(Et.subs(k, k2) - Et.subs(k, k1))
dP_expected = sp.diag(b[1]*(k2**2 - k1**2)/a[1], 0, b[2]*(k2**2 - k1**2)/a[2], 0)
check("A3.3  E~(k2) - E~(k1) = diag(b_i (k2^2 - k1^2)/a_i) on the n-entries, zero elsewhere",
      sp.simplify(dP - dP_expected) == sp.zeros(4, 4))
# Courant-Fischer / Weyl: lambda_min(A + P) >= lambda_min(A) for P >= 0 (theorem, cited, not re-proved here)
# hence lambda_min(k) non-decreasing in k^2, and min over any k grid containing k = 0 is attained at k = 0.

# ---------------------------------------------------------------- one-component limit
a1s = sp.symbols('a', positive=True)
V1s = sp.symbols('V', real=True)
E1 = sp.Matrix([[1, V1s/sp.sqrt(a1s)], [V1s/sp.sqrt(a1s), 1]])
det1 = sp.simplify(E1.det())
check("A3.4  one component at k = 0: det = 1 - V^2/a, trace = 2", sp.simplify(det1 - (1 - V1s**2/a1s)) == 0 and E1.trace() == 2)
# 2x2 symmetric is positive definite iff trace > 0 and det > 0, i.e. iff V^2 < a
check("A3.5  one component: positive definite  <=>  V^2 < a (Landau at k -> 0)",
      sp.solve_univariate_inequality(det1 > 0, V1s, relational=False) == sp.Interval.open(-sp.sqrt(a1s), sp.sqrt(a1s)))
# the one-component restriction of E~ (drop component 2, set f = 0) reproduces E1
E_one = Et_expected.subs({g12: 0}).extract([0, 1], [0, 1]).subs(k, 0)
check("A3.6  the component-1 block of E~(0) at f = 0 equals the one-component matrix",
      sp.simplify(E_one - E1.subs({a1s: a[1], V1s: V[1]})) == sp.zeros(2, 2))

# ---------------------------------------------------------------- summary
n_fail = sum(1 for _, ok in results if not ok)
print(f"\nSUMMARY: {len(results) - n_fail}/{len(results)} PASS")
sys.exit(0 if n_fail == 0 else 1)
