"""Symbolic checks for UV_STRUCTURAL_CONSTRAINTS_NOTE_v0.2.md (sections A, B, C).
Pure algebra (SymPy); reads no data and has no numerical thresholds.
Usage: python3 -I uv_structural_constraints_check.py
Conventions: signature (+,-), 1+1D, u frame: X1 = u.dphi, X2 = n.dphi, X3 = u.dpsi, X4 = n.dpsi;
B4 (docs/notation_and_conventions.md, m = 0):
  L_B4 = 1/2(X1^2 - X2^2) - Lam/4 X2^4 + 1/2(X3^2 - X4^2) + p/2 X1^2 X2^2 - kappa X1 X2 X3 X4 + s/2 X3^2 X4^2
"""
import sys
import sympy as sp

X1, X2, X3, X4 = sp.symbols('X1 X2 X3 X4', real=True)
p, s, kap, Lam = sp.symbols('p s kappa Lambda', real=True)
out = []


def report(name, ok, detail=''):
    out.append((name, bool(ok)))
    print('%-5s %s %s' % ('OK' if ok else 'FAIL', name, detail))


def quartic(L):
    """quartic part of a polynomial Lagrangian in X1..X4."""
    P = sp.Poly(sp.expand(L), X1, X2, X3, X4)
    return sum(cf*sp.Mul(*[x**e for x, e in zip((X1, X2, X3, X4), m)]) for m, cf in P.terms() if sum(m) == 4)


B4q = -Lam/4*X2**4 + p/2*X1**2*X2**2 - kap*X1*X2*X3*X4 + s/2*X3**2*X4**2

# ---------------- A: no B4 quartic is Lorentz invariant ----------------
# Lorentz-invariant local quartics in first derivatives of two scalars: products of
# Xp = X1^2 - X2^2, Xq = X3^2 - X4^2, Y = X1 X3 - X2 X4 (the invariants are frame independent; written in the u frame).
Xp, Xq, Y = X1**2 - X2**2, X3**2 - X4**2, X1*X3 - X2*X4
basis = [Xp**2, Xq**2, Y**2, Xp*Xq, Xp*Y, Xq*Y]
a = sp.symbols('a0:6', real=True)
eqs = sp.Poly(sp.expand(sum(ai*bi for ai, bi in zip(a, basis)) - B4q), X1, X2, X3, X4).coeffs()
sol = sp.solve(eqs, list(a) + [p, s, kap, Lam], dict=True)
onlyfree = all(all(sp.simplify(d.get(v, v)) == 0 for v in (p, s, kap, Lam)) for d in sol) and len(sol) > 0
report('A1 B4 quartic in the span of Lorentz-invariant quartics only for p = s = kappa = Lambda = 0', onlyfree, sol)
for nm, op, cv in (('Lambda X2^4', -Lam/4*X2**4, Lam), ('p X1^2 X2^2', p/2*X1**2*X2**2, p),
                   ('kappa X1X2X3X4', -kap*X1*X2*X3*X4, kap), ('s X3^2 X4^2', s/2*X3**2*X4**2, s)):
    e2 = sp.Poly(sp.expand(sum(ai*bi for ai, bi in zip(a, basis)) - op), X1, X2, X3, X4).coeffs()
    s2 = sp.solve(e2, list(a) + [cv], dict=True)          # unknowns: the six invariant coefficients and this operator's own coefficient
    report('A2 single operator %s is not Lorentz invariant unless its coefficient is 0' % nm,
           len(s2) > 0 and all(sp.simplify(d.get(cv, cv)) == 0 for d in s2), s2)

# ---------------- B: single mediator -> rank one; healthy scalars -> kappa^2 <= p s ----------------
M = sp.symbols('M', positive=True)
g1, g2 = sp.symbols('g1 g2', real=True)
O1, O2 = X1*X2, X3*X4                               # (u.dphi)(n.dphi), (u.dpsi)(n.dpsi)
chi = sp.symbols('chi', real=True)
# leading local order: constant fields, L_chi = -1/2 M^2 chi^2 + chi J  (gradient terms of chi are O(d^2/M^4))
J = g1*O1 + g2*O2
Lchi = -M**2/2*chi**2 + chi*J
chis = sp.solve(sp.diff(Lchi, chi), chi)[0]
Leff = sp.expand(Lchi.subs(chi, chis))
report('B1 scalar elimination: L_eff = J^2/(2 M^2)', sp.expand(Leff - J**2/(2*M**2)) == 0)
q = quartic(Leff)
P = sp.Poly(q, X1, X2, X3, X4)
p_e = 2*P.coeff_monomial(X1**2*X2**2)
s_e = 2*P.coeff_monomial(X3**2*X4**2)
k_e = -P.coeff_monomial(X1*X2*X3*X4)
report('B2 matching p = g1^2/M^2, s = g2^2/M^2, kappa = -g1 g2/M^2',
       sp.simplify(p_e - g1**2/M**2) == 0 and sp.simplify(s_e - g2**2/M**2) == 0 and sp.simplify(k_e + g1*g2/M**2) == 0)
report('B3 single mediator: kappa^2 = p s exactly', sp.simplify(k_e**2 - p_e*s_e) == 0)
# N healthy scalar mediators (N = 3 written out): p s - kappa^2 = Cauchy-Binet sum of squares >= 0
N = 3
G1 = sp.symbols('g1_0:%d' % N, real=True)
G2 = sp.symbols('g2_0:%d' % N, real=True)
Ms = sp.symbols('M_0:%d' % N, positive=True)
pN = sum(G1[i]**2/Ms[i]**2 for i in range(N))
sN = sum(G2[i]**2/Ms[i]**2 for i in range(N))
kN = -sum(G1[i]*G2[i]/Ms[i]**2 for i in range(N))
cb = sum((G1[i]*G2[j] - G1[j]*G2[i])**2/(Ms[i]**2*Ms[j]**2) for i in range(N) for j in range(i + 1, N))
report('B4 N = 3 healthy scalars: p s - kappa^2 = sum_{a<b} (g1a g2b - g1b g2a)^2/(Ma^2 Mb^2) >= 0',
       sp.simplify(pN*sN - kN**2 - cb) == 0)
# a ghost-like mediator (wrong-sign mass term in this leading-order sense) flips the sign of its contribution
Lgh = +M**2/2*chi**2 + chi*J
chig = sp.solve(sp.diff(Lgh, chi), chi)[0]
report('B5 wrong-sign mediator gives L_eff = -J^2/(2 M^2) (shows the PSD statement needs healthy mediators)',
       sp.expand(Lgh.subs(chi, chig) + J**2/(2*M**2)) == 0)
# extra operator: coupling also to (n.dphi)^2 creates X1 X2^3, which is not a B4 operator
g3 = sp.symbols('g3', real=True)
J3 = g1*O1 + g2*O2 + g3*X2**2
q3 = sp.Poly(sp.expand(J3**2/(2*M**2)), X1, X2, X3, X4)
report('B6 coupling to (n.dphi)^2 as well generates X1 X2^3 (not in B4) with coefficient g1 g3/M^2',
       sp.simplify(q3.coeff_monomial(X1*X2**3) - g1*g3/M**2) == 0)

# vector mediator with currents in independent directions: the coefficient matrix is a Gram matrix, not rank one
A0_, A1_ = sp.symbols('A0_ A1_', real=True)
dot2 = lambda U, V: U[0]*V[0] - U[1]*V[1]
Jvec = (g1*O1, g2*O2)                                   # J^mu = g1 u^mu O1 + g2 n^mu O2 (u frame components)
LAv = M**2/2*dot2((A0_, A1_), (A0_, A1_)) + dot2((A0_, A1_), Jvec)
LeffV = sp.expand(LAv.subs(sp.solve([sp.diff(LAv, A0_), sp.diff(LAv, A1_)], [A0_, A1_], dict=True)[0]))
PV = sp.Poly(quartic(LeffV), X1, X2, X3, X4)
pV, sV, kV = 2*PV.coeff_monomial(X1**2*X2**2), 2*PV.coeff_monomial(X3**2*X4**2), -PV.coeff_monomial(X1*X2*X3*X4)
report('B7 one Proca mediator, currents along u and n: p = -g1^2/M^2, s = g2^2/M^2, kappa = 0, so kappa^2 != p s (rank two)',
       sp.simplify(pV + g1**2/M**2) == 0 and sp.simplify(sV - g2**2/M**2) == 0 and kV == 0 and sp.simplify(kV**2 - pV*sV) != 0,
       'p=%s s=%s kappa=%s' % (pV, sV, kV))
# common polarisation n^mu for both operators: rank one again
Jcom = (0, g1*O1 + g2*O2)
LAc = M**2/2*dot2((A0_, A1_), (A0_, A1_)) + dot2((A0_, A1_), Jcom)
LeffC = sp.expand(LAc.subs(sp.solve([sp.diff(LAc, A0_), sp.diff(LAc, A1_)], [A0_, A1_], dict=True)[0]))
PC = sp.Poly(quartic(LeffC), X1, X2, X3, X4)
pC, sC, kC = 2*PC.coeff_monomial(X1**2*X2**2), 2*PC.coeff_monomial(X3**2*X4**2), -PC.coeff_monomial(X1*X2*X3*X4)
report('B8 one Proca mediator with a common polarisation: kappa^2 = p s (rank one)', sp.simplify(kC**2 - pC*sC) == 0,
       'p=%s s=%s kappa=%s' % (pC, sC, kC))

# ---------------- C: tree-level sign of Lambda ----------------
h, al, be = sp.symbols('h alpha beta', real=True)
# scalar mediator coupled to a general bilinear of dphi: O = alpha X1^2 + beta X1 X2 + h X2^2
Os = al*X1**2 + be*X1*X2 + h*X2**2
qs = sp.Poly(sp.expand(Os**2/(2*M**2)), X1, X2)
Lam_s = -4*qs.coeff_monomial(X2**4)
report('C1 healthy scalar, any bilinear coupling: Lambda_eff = -2 h^2/M^2 <= 0', sp.simplify(Lam_s + 2*h**2/M**2) == 0)
# N healthy scalars: Lambda_eff = -2 sum h_a^2/M_a^2 <= 0
H = sp.symbols('h_0:%d' % N, real=True)
LamN = -4*sum(H[i]**2/(2*Ms[i]**2) for i in range(N))
report('C2 N healthy scalars: Lambda_eff = -2 sum h_a^2/M_a^2 <= 0', sp.simplify(LamN + 2*sum(H[i]**2/Ms[i]**2 for i in range(N))) == 0)
# Proca mediator at leading local order: L = 1/2 M^2 A.A + A.J (F^2 terms are O(d^2)); metric (+,-); u = (1,0), n = (0,1)
A0, A1 = sp.symbols('A0 A1', real=True)
gu, gn = sp.symbols('gamma_u gamma_n', real=True)
Ov = X2**2
Jv = (gu*Ov, gn*Ov)                                     # J^mu = gamma_u u^mu O + gamma_n n^mu O (u frame components)
dot = lambda U, V: U[0]*V[0] - U[1]*V[1]
LA = M**2/2*dot((A0, A1), (A0, A1)) + dot((A0, A1), Jv)
solA = sp.solve([sp.diff(LA, A0), sp.diff(LA, A1)], [A0, A1], dict=True)[0]
LeffA = sp.expand(LA.subs(solA))
report('C3 Proca elimination: L_eff = -J.J/(2 M^2)', sp.expand(LeffA + dot(Jv, Jv)/(2*M**2)) == 0)
qa = sp.Poly(LeffA, X2)
Lam_v = sp.simplify(-4*qa.coeff_monomial(X2**4))
report('C4 Proca: Lambda_eff = 2 (gamma_u^2 - gamma_n^2)/M^2 (positive for a current along u)',
       sp.simplify(Lam_v - 2*(gu**2 - gn**2)/M**2) == 0, 'Lambda_eff = %s' % Lam_v)

n_fail = sum(1 for _, ok in out if not ok)
print('TOTAL %d checks, %d failed' % (len(out), n_fail))
sys.exit(1 if n_fail else 0)
