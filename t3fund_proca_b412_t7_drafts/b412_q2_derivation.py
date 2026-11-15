"""B4.12 Q2 — analytic derivation of H-Q2 and epsilon_EFT (specification: B4_12_Q2_DERIVATION_SPEC_v0.2.md, commit 641121c).

Usage (only after authorisation):
    python -I b412_q2_derivation.py <outdir>

Symbolic only (SymPy). No floating-point evaluation, no registered background.
Setting: registered T3-fund Proca action, 1+1D, metric (-,+), constant-gradient background (nu = 0).
Arbiter: the registered action and conventions; SymPy is a tool.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import sympy as sp
from sympy.calculus.euler import euler_equations

SPEC = "B4_12_Q2_DERIVATION_SPEC_v0.2.md"

# ---------------------------------------------------------------- symbols
t, x = sp.symbols("t x", real=True)
om, kk = sp.symbols("omega k", real=True)                    # Fourier: fields ~ exp(i(k x - omega t))
w, q, v, r = sp.symbols("w q v r", real=True)                # background: f = w t + q x, g = v t + r x
M, p = sp.symbols("M p", positive=True)
kap = sp.symbols("kappa", real=True)
Lam = sp.symbols("Lambda", nonnegative=True)
gph, gps = sp.symbols("g_phi g_psi", real=True)              # generic couplings (fixed-EFT substitution below)
eps, cc = sp.symbols("epsilon c_loc")
A0b, A1b = sp.symbols("A0bar A1bar")
c, u, c0 = sp.symbols("c u c0")                              # c = omega/k, u = (k/M)^2, c0 = M^-2 light root
s_b4 = sp.symbols("s", positive=True)                        # independent B4 coefficient (Q2-C3)

df, dg, a0, a1, S = [sp.Function(n)(t, x) for n in ("df", "dg", "a0", "a1", "S")]
U1, U2, U3, U4, Sa = sp.symbols("U1 U2 U3 U4 Sa")
LIGHT, HEAVY = [df, dg], [a0, a1]
AMPS = {df: U1, dg: U2, a0: U3, a1: U4, S: Sa}

FIXED_EFT = {gph: M * sp.sqrt(p), gps: -M * kap / sp.sqrt(p)}   # g_psi^2/M^2 = kappa^2/p (= s only at rank one)
K_T7A = (1 + kk**2 / M**2) / (M**2 + kk**2 - om**2)             # T7a step 3 (61fefff)


def D(e, var):
    return sp.diff(e, var)


def fourier(expr):
    """Replace every derivative of a field by (-i omega)^a (i k)^b times its amplitude (constant-coefficient linear)."""
    rep = {}
    for d in expr.atoms(sp.Derivative):
        if d.expr in AMPS:
            fac = sp.Integer(1)
            for var, n in d.variable_count:
                fac *= (-sp.I * om if var == t else sp.I * kk) ** n
            rep[d] = fac * AMPS[d.expr]
    expr = expr.xreplace(rep).xreplace(AMPS)
    return sp.expand(expr)


# Revision 5 (performance only): division-free Berkowitz determinant instead of the default Bareiss method, whose
# polynomial gcd steps did not finish within 30 min on the symbolic 4x4 matrix. The determinant is the same polynomial.
DET = "berkowitz"
IDENTICALLY_ZERO_EL = []      # fields whose Euler-Lagrange expression was identically zero (reported)


def el_matrix(L, fields, label):
    """Euler-Lagrange equations of L for `fields`, in Fourier space, as a matrix (rows: equations, cols: amplitudes)."""
    # One field at a time. sympy's euler_equations silently drops an equation that is identically 0 = 0 (it evaluates
    # to True). Two cases are distinguished:
    #   field ABSENT from L                -> error (lesson of T7a step 4, C6);
    #   field PRESENT, equation dropped    -> its Euler-Lagrange expression is identically zero; recorded as 0.
    amps = [AMPS[f] for f in fields]
    rows = []
    for f in fields:
        if not L.has(f):
            raise RuntimeError("field %s is absent from the Lagrangian" % f)
        eq = euler_equations(L, [f], [t, x])
        if len(eq) == 0:
            IDENTICALLY_ZERO_EL.append("%s in %s" % (f, label))
            ex = sp.Integer(0)
        elif len(eq) == 1:
            ex = fourier(eq[0].lhs - eq[0].rhs)
        else:
            raise RuntimeError("unexpected number of equations for field %s" % f)
        rows.append([ex.coeff(a) for a in amps])
        rest = sp.expand(ex - sum(ex.coeff(a) * a for a in amps))
        rows[-1].append(rest)          # inhomogeneous remainder (must be 0 for the fluctuation system)
    return rows


def is_zero(e):
    return sp.simplify(sp.cancel(sp.together(e))) == 0


def mat_zero(Mx):
    return all(is_zero(e) for e in Mx)


def main(outdir: str) -> int:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    T0 = time.time()
    log, res, checks = [], {}, {}

    def note(s):
        log.append(s)
        print(s, flush=True)

    # ------------------------------------------------ action (registered T3-fund Proca, section 1)
    fe = w * t + q * x + eps * df
    ge = v * t + r * x + eps * dg
    A0e = A0b + eps * a0
    A1e = A1b + eps * a1
    J = gph * D(fe, t) * D(fe, x) + gps * D(ge, t) * D(ge, x)
    F = D(A1e, t) - D(A0e, x)
    L_light = (sp.Rational(1, 2) * D(fe, t)**2 - sp.Rational(1, 2) * D(fe, x)**2 - Lam / 4 * D(fe, x)**4
               + sp.Rational(1, 2) * D(ge, t)**2 - sp.Rational(1, 2) * D(ge, x)**2)
    L = sp.expand(L_light + sp.Rational(1, 2) * F**2 + sp.Rational(1, 2) * M**2 * (A0e**2 - A1e**2) + A1e * J)

    # ------------------------------------------------ D1: background (linear order) and quadratic action
    L1 = L.coeff(eps, 1)
    bg_rows = el_matrix(L1, LIGHT + HEAVY, "L1_background")
    bg_eqs = [row[-1] for row in bg_rows]                     # linear order: only the inhomogeneous parts
    sol = sp.solve([bg_eqs[2], bg_eqs[3]], [A0b, A1b], dict=True)
    if len(sol) != 1:
        raise RuntimeError("background A not unique")
    sol = sol[0]
    res["D1_background_A"] = {str(kv): str(vv) for kv, vv in sol.items()}
    checks["D1_light_background_EOM_zero"] = all(is_zero(e.subs(sol)) for e in bg_eqs[:2])
    note("D1 background A: %s ; light EOM zero: %s" % (sol, checks["D1_light_background_EOM_zero"]))

    L2 = sp.expand(L.coeff(eps, 2).subs(sol))
    full = el_matrix(L2, LIGHT + HEAVY, "L2_full")
    checks["D1_homogeneous_fluctuation_system"] = all(is_zero(row[-1]) for row in full)
    Mfull = sp.Matrix([row[:-1] for row in full])
    Mll, Mla, Mal, Maa = Mfull[:2, :2], Mfull[:2, 2:], Mfull[2:, :2], Mfull[2:, 2:]
    Mred = (Mll - Mla * Maa.inv() * Mal).applyfunc(lambda e: sp.cancel(sp.together(e)))

    # ------------------------------------------------ Q2-C1: kernel from the Proca block with an external source
    a0e_, a1e_ = a0, a1
    Fs = D(a1e_, t) - D(a0e_, x)
    Lsrc = sp.Rational(1, 2) * Fs**2 + sp.Rational(1, 2) * M**2 * (a0e_**2 - a1e_**2) + a1e_ * S
    src = el_matrix(Lsrc, HEAVY, "Lsrc_kernel")
    eqs_src = [row[0] * U3 + row[1] * U4 + row[2] for row in src]
    ssol = sp.solve(eqs_src, [U3, U4], dict=True)[0]
    K_derived = sp.cancel(ssol[U4] / Sa)
    checks["Q2_C1_kernel"] = is_zero(K_derived - K_T7A)
    res["Q2_C1_K_derived"] = sp.srepr(K_derived)
    note("Q2-C1 kernel equals T7a K: %s" % checks["Q2_C1_kernel"])

    # ------------------------------------------------ Q2-C2: H-Q2  (Mred == M0 + K N), same EL convention
    J1 = sp.expand(J).coeff(eps, 1)
    L2_light = sp.expand(L2.subs({a0: 0, a1: 0}).doit())          # light kinetic + Lambda + Abar1*J2
    Lloc = sp.expand(L2_light + cc / 2 * J1**2)                     # local stand-in: K -> constant c_loc
    loc = el_matrix(Lloc, LIGHT, "Lloc_HQ2")
    Mloc = sp.Matrix([row[:-1] for row in loc])
    M0, N = Mloc.subs(cc, 0), Mloc.diff(cc)
    remainder = (Mred - (M0 + K_derived * N)).applyfunc(lambda e: sp.cancel(sp.together(e)))
    checks["Q2_C2_HQ2"] = mat_zero(remainder)
    res["Q2_C2_remainder_srepr"] = [sp.srepr(e) for e in remainder]
    note("Q2-C2 H-Q2 identity (remainder zero): %s" % checks["Q2_C2_HQ2"])

    # ------------------------------------------------ Q2-C3: M^-2 truncation vs B4 (rank-one subfamily; general B4)
    Mred2 = (M0 + N / M**2).subs(FIXED_EFT).applyfunc(sp.cancel)
    fB, gB = w * t + q * x + eps * df, v * t + r * x + eps * dg
    wB, qB, vB, rB = D(fB, t), D(fB, x), D(gB, t), D(gB, x)
    L_B4 = sp.expand(sp.Rational(1, 2) * wB**2 - sp.Rational(1, 2) * qB**2 - Lam / 4 * qB**4
                     + sp.Rational(1, 2) * vB**2 - sp.Rational(1, 2) * rB**2           # m = 0: no psi mass in the Proca action
                     + p / 2 * wB**2 * qB**2 - kap * wB * qB * vB * rB + s_b4 / 2 * vB**2 * rB**2)
    b4 = el_matrix(sp.expand(L_B4).coeff(eps, 2), LIGHT, "L_B4")
    MB4 = sp.Matrix([row[:-1] for row in b4])
    checks["Q2_C3_rank_one_match"] = mat_zero(Mred2 - MB4.subs(s_b4, kap**2 / p))
    diff_general = (MB4 - Mred2).applyfunc(sp.cancel)
    ratio = diff_general.applyfunc(lambda e: sp.cancel(e / (s_b4 - kap**2 / p)))
    checks["Q2_C3_general_B4_differs_by_(s-kappa^2/p)"] = (not mat_zero(diff_general)) and all(
        sp.diff(e, s_b4) == 0 for e in ratio)
    res["Q2_C3_general_B4_difference_over_(s-kappa^2/p)"] = [str(e) for e in ratio]
    note("Q2-C3 rank-one match: %s ; general B4 differs only via (s - kappa^2/p): %s"
         % (checks["Q2_C3_rank_one_match"], checks["Q2_C3_general_B4_differs_by_(s-kappa^2/p)"]))

    # ------------------------------------------------ Q2-C4: decoupling of the FULL system at g_phi = g_psi = 0
    M_dec = Mfull.subs({gph: 0, gps: 0})
    off = list(M_dec[:2, 2:]) + list(M_dec[2:, :2])
    roots_f = sp.solve(sp.factor(M_dec[0, 0]), om)
    roots_g = sp.solve(sp.factor(M_dec[1, 1]), om)
    roots_A = sp.solve(sp.factor(M_dec[2:, 2:].det(method=DET)), om)
    checks["Q2_C4_decoupling"] = (all(is_zero(e) for e in off) and is_zero(M_dec[0, 1]) and is_zero(M_dec[1, 0]))
    res["Q2_C4_branches"] = {"light_f": [str(z) for z in roots_f], "light_g": [str(z) for z in roots_g],
                             "proca": [str(z) for z in roots_A]}
    note("Q2-C4 decoupled: %s ; light f: %s ; light g: %s ; Proca: %s"
         % (checks["Q2_C4_decoupling"], roots_f, roots_g, roots_A))

    # ------------------------------------------------ D2(b): full system, Schur identity, pole roots
    Mfe = Mfull.subs(FIXED_EFT)
    Mrede = Mred.subs(FIXED_EFT).applyfunc(sp.cancel)
    Pfull = sp.expand(Mfe.det(method=DET))
    Daa = sp.cancel(Mfe[2:, 2:].det(method=DET))
    checks["D2b_schur_identity"] = is_zero(Pfull - Daa * Mrede.det(method=DET))
    pole = M**2 + kk**2 - om**2
    num_red = sp.numer(sp.together(Mrede.det(method=DET)))
    rem_red = sp.rem(sp.Poly(num_red, om), sp.Poly(pole, om)).as_expr()
    res["D2b_full_polynomial_degree_in_omega"] = int(sp.degree(Pfull, om))
    res["D2b_Daa"] = str(sp.factor(Daa))
    res["D2b_reduced_numerator_divisible_by_pole"] = bool(sp.simplify(rem_red) == 0)
    # both poles separately: the background may break omega -> -omega
    om_pole = {"plus": sp.sqrt(M**2 + kk**2), "minus": -sp.sqrt(M**2 + kk**2)}
    res["D2b_Pfull_at_pole_is_zero"] = {sg: bool(sp.expand(Pfull.subs(om, val)) == 0) for sg, val in om_pole.items()}
    note("D2(b) Schur identity: %s ; deg_omega Pfull = %s ; Pfull(+pole)=0: %s ; Pfull(-pole)=0: %s ; "
         "reduced numerator divisible by pole: %s"
         % (checks["D2b_schur_identity"], res["D2b_full_polynomial_degree_in_omega"],
            res["D2b_Pfull_at_pole_is_zero"]["plus"], res["D2b_Pfull_at_pole_is_zero"]["minus"],
            res["D2b_reduced_numerator_divisible_by_pole"]))

    # heavy branch: continuation from the decoupled Proca branch (implicit function theorem in the coupling scale eta)
    eta = sp.symbols("eta", real=True)
    P_eta = sp.expand(Mfull.subs({gph: eta * gph, gps: eta * gps}, simultaneous=True).det(method=DET))
    heavy = {}
    N_ETA = 4                                         # highest eta-derivative examined for the leading correction
    for sg, wh in om_pole.items():
        at_h = {om: wh, eta: 0}
        P0 = sp.expand(P_eta.subs(at_h))
        Pw = sp.factor(sp.expand(sp.diff(P_eta, om).subs(at_h)))
        generic_simple = (P0 == 0) and (Pw != 0)      # Pw not identically zero: simple for GENERIC parameters only
        entry = dict(decoupled_root=str(wh), is_root_at_eta0=bool(P0 == 0),
                     P_omega_at_root=str(Pw),
                     simple_root_generic=bool(generic_simple),
                     exceptional_set="parameters where P_omega_at_root = 0 (resonance of the heavy root with a light "
                                     "root, or a degenerate Proca factor); continuation NOT established there")
        if generic_simple:
            # leading correction: first n with d^n P/d eta^n != 0 at (omega_h, 0); then
            # omega - omega_h = -(d^n P/d eta^n / n!) / P_omega * eta^n + O(eta^(n+1))
            lead = None
            for n_ in range(1, N_ETA + 1):
                Pn = sp.expand(sp.diff(P_eta, eta, n_).subs(at_h))
                if Pn != 0:
                    lead = (n_, sp.cancel(-Pn / sp.factorial(n_) / Pw))
                    break
            if lead is None:
                entry["leading_correction"] = ("not determined: d^n P/d eta^n = 0 identically at the decoupled root for "
                                               "n = 1..%d; no non-zero correction through this order" % N_ETA)
            else:
                entry["leading_correction"] = ("omega = omega_h + a_n*eta**%d + O(eta**%d) (a_n not identically zero; "
                                               "generic parameters)" % (lead[0], lead[0] + 1))
                entry["a_n"] = str(lead[1])
        entry["status"] = ("implicitly defined branch for generic (non-resonant) parameters: the unique root of Pfull "
                           "continuously connected to the decoupled Proca root (implicit function theorem); no closed "
                           "form is claimed; not established on the exceptional set") if generic_simple else \
            "continuation not established (decoupled root not simple, or not a root)"
        heavy[sg] = entry
    res["D2b_heavy_branch"] = heavy
    checks["D2b_heavy_continuation_defined_generic"] = all(e["simple_root_generic"] for e in heavy.values())
    note("D2(b) heavy branch by continuation (simple for generic parameters): %s"
         % checks["D2b_heavy_continuation_defined_generic"])

    # ------------------------------------------------ light relation in c = omega/k, u = (k/M)^2
    to_cu = {om: c * kk, M: kk / sp.sqrt(u)}
    detc = sp.cancel(sp.together(Mrede.det(method=DET).subs(to_cu) / kk**4))
    if kk in detc.free_symbols or M in detc.free_symbols:
        raise RuntimeError("light relation not a function of c and u only")
    Fcu = sp.expand(sp.numer(sp.together(detc)))
    F0 = sp.expand(Fcu.subs(u, 0))
    res["light_relation_F_c_u"] = sp.srepr(Fcu)
    res["F0_c_(M^-2 light polynomial)"] = str(sp.factor(F0))
    res["F0_degree_in_c"] = int(sp.degree(F0, c))

    # ------------------------------------------------ Q2-C6: M^2 K - 1 on the light variables
    Kh = sp.cancel((M**2 * K_derived).subs(to_cu))
    checks["Q2_C6_M2K_minus_1"] = is_zero(Kh - 1 - u * c**2 / (1 + u * (1 - c**2)))

    # ------------------------------------------------ D3 / Q2-C5: two independent routes for the u-expansion
    Fc, Fu = sp.diff(Fcu, c), sp.diff(Fcu, u)
    at0 = {c: c0, u: 0}
    # Route A: implicit-function derivatives of the exact light root
    # Revision 6 (performance only): sp.together instead of sp.cancel. No gcd reduction; the rational expression is
    # the same, only not reduced to lowest terms. All comparisons below act on numerators via pseudo-remainder.
    cA1 = sp.together(-(Fu / Fc).subs(at0))
    cA2 = sp.together((-(sp.diff(Fcu, c, 2) * cA1**2 + 2 * sp.diff(Fcu, c, u) * cA1 + sp.diff(Fcu, u, 2))
                       / (2 * Fc)).subs(at0))
    # Route B: perturbation theory, ansatz c = c0 + b1 u + b2 u^2 in the dispersion polynomial
    b1, b2 = sp.symbols("b1 b2")
    ser = sp.expand(Fcu.subs(c, c0 + b1 * u + b2 * u**2))      # exact polynomial; only u^1, u^2 coefficients used
    eq1, eq2 = ser.coeff(u, 1), ser.coeff(u, 2)
    # eq1 is linear in b1, eq2 is linear in b2: the unique solutions are written explicitly (same equations as
    # sp.solve, without its simplification step)
    eq1p, eq2p = sp.Poly(eq1, b1), sp.Poly(eq2, b2)
    if eq1p.degree() != 1 or eq2p.degree() != 1:
        raise RuntimeError("route B coefficient equations are not linear in b1, b2")
    cB1 = sp.together(-eq1p.coeff_monomial(1) / eq1p.coeff_monomial(b1))
    cB2 = sp.together((-eq2p.coeff_monomial(1) / eq2p.coeff_monomial(b2)).subs(b1, cB1))
    F0c = sp.Poly(F0.subs(c, c0), c0)

    def equal_mod_F0(e1, e2):
        num = sp.numer(sp.together(e1 - e2))
        return sp.prem(sp.Poly(sp.expand(num), c0), F0c).is_zero

    checks["Q2_C5_order_u1"] = equal_mod_F0(cA1, cB1)
    checks["Q2_C5_order_u2"] = equal_mod_F0(cA2, cB2)
    note("Q2-C5 two routes agree: u^1 %s, u^2 %s" % (checks["Q2_C5_order_u1"], checks["Q2_C5_order_u2"]))

    # signed coefficients; the error metrics of v0.2 are their absolute values
    eEFT1 = sp.together(2 * c0 * cA1 / (c0**2 - 1))   # eps_EFT = |eEFT1| u + O(u^2), c0**2 != 1, c0 != 0
    eTOT1 = sp.together(2 * cA1 / c0)                   # eps_tot = |eTOT1| u + O(u^2), c0 != 0
    eABS2 = sp.together(2 * c0 * cA1)                   # eps_abs = |eABS2| u^2 + O(u^3), c0 != 0
    res["D3_c1"], res["D3_c2"] = sp.srepr(cA1), sp.srepr(cA2)
    res["D3_eps_EFT_leading"] = dict(order="u", signed_coeff=sp.srepr(eEFT1), metric=sp.srepr(sp.Abs(eEFT1)))
    res["D3_eps_tot_leading"] = dict(order="u", signed_coeff=sp.srepr(eTOT1), metric=sp.srepr(sp.Abs(eTOT1)))
    res["D3_eps_abs_leading"] = dict(order="u**2", signed_coeff=sp.srepr(eABS2), metric=sp.srepr(sp.Abs(eABS2)))
    res["D3_special_cases"] = {
        "c0**2 == 1": "eps_EFT undefined on that branch (no B4 effect); only eps_abs reported. "
                      "F0(1) = %s ; F0(-1) = %s" % (sp.factor(F0.subs(c, 1)), sp.factor(F0.subs(c, -1))),
        "c0 == 0": "branch exists iff F0(0) = 0, F0(0) = %s. There omega_2**2 = 0, so eps_tot = |0 - c**2|/c**2 = 1 "
                   "whenever c(u) != 0 (exact, not a leading-order formula); eps_tot undefined if c(u) == 0. "
                   "CONDITIONAL on c1 != 0, with c1 = c1(c0 = 0) = %s (valid where F_c(0,0) != 0): "
                   "eps_EFT = c**2/|c**2 - 1| = c1**2 u**2 + O(u**3) and eps_abs = c**2 u = c1**2 u**3 + O(u**4). "
                   "If c1 = 0 the first non-zero coefficient of c(u) is not derived here: UNRESOLVED"
                   % (sp.factor(F0.subs(c, 0)), sp.together(cA1.subs(c0, 0))),
    }

    # ------------------------------------------------ D3 (iii): rho dependence at fixed X (scaling test)
    lam = sp.symbols("lambda_s", positive=True)
    scale = {w: lam * w, q: lam * q, v: lam * v, r: lam * r, p: p / lam**2, kap: kap / lam**2, Lam: Lam / lam**2}
    F_scaled = Fcu.subs(scale, simultaneous=True)
    ratio_s = sp.cancel(F_scaled / Fcu)                       # invariant up to an overall factor free of c and u
    checks["D3_iii_F_invariant_at_fixed_X"] = is_zero(sp.diff(ratio_s, c)) and is_zero(sp.diff(ratio_s, u))
    res["D3_iii_overall_scaling_factor"] = str(ratio_s)
    note("D3(iii) light relation depends on the background only through X-invariant combinations: %s"
         % checks["D3_iii_F_invariant_at_fixed_X"])

    # ------------------------------------------------ D4: admissibility (pole condition along a light branch)
    res["D4b_no_pole_crossing"] = ("necessary condition along a light branch: 1 + u*(1 - c_b**2) > 0, i.e. for "
                                   "c_b**2 > 1: (k/M)**2 < 1/(c_b**2 - 1); evaluated at leading order c_b -> c0")
    res["D4a_reality"] = ("NOT DERIVED for finite k/M. Only the limit k/M -> 0 is covered: real light branches there "
                          "iff all roots of F0(c) are real. Reality at finite k/M is open; T3-fund C3 is cited for the "
                          "full system on its registered grid, not as a proof here")
    res["D4c_separation_from_heavy"] = ("NOT DERIVED. Separation of light and heavy branches is a separate spectral "
                                        "property; the pole condition D4b is not claimed to be sufficient for it. Open")

    # ------------------------------------------------ D5: structural classification (no eps_max)
    nz = not sp.prem(sp.Poly(sp.expand(sp.numer(sp.together(eEFT1))), c0), F0c).is_zero
    res["D5_scope"] = ("GENERIC light branches only: simple roots c0 of F0 with c0 != 0 and c0**2 != 1. The special "
                       "branches of D3 (c0 = 0, c0**2 = 1) and degenerate roots are excluded from this classification")
    res["D5_classification"] = ("(ii) on generic branches: eps_EFT depends on k/M (leading ~ (k/M)^2, coefficient not "
                                "identically zero); at fixed X-invariants not on rho separately"
                                if (nz and checks["D3_iii_F_invariant_at_fixed_X"]) else
                                "(iii) not determined by this derivation")

    mandatory = ["D1_light_background_EOM_zero", "D1_homogeneous_fluctuation_system", "Q2_C1_kernel", "Q2_C2_HQ2",
                 "Q2_C3_rank_one_match", "Q2_C3_general_B4_differs_by_(s-kappa^2/p)", "Q2_C4_decoupling",
                 "Q2_C5_order_u1", "Q2_C5_order_u2", "Q2_C6_M2K_minus_1", "D2b_schur_identity",
                 "D2b_heavy_continuation_defined_generic"]
    verdict = "ALL_CHECKS_PASS" if all(checks[k] for k in mandatory) else "CHECKS_FAILED"
    res["identically_zero_EL_equations"] = list(IDENTICALLY_ZERO_EL)
    summary = dict(spec=SPEC, sympy=sp.__version__, seconds=round(time.time() - T0, 1),
                   checks={k: bool(vv) for k, vv in checks.items()}, mandatory=mandatory, results=res, verdict=verdict)
    (out / "b412_q2_results.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8", newline="\n")
    for k in mandatory:
        note("%s %s" % (k, "ZERO/PASS" if checks[k] else "FAIL"))
    note(verdict)
    (out / "b412_q2_output.txt").write_text("\n".join(log) + "\n", encoding="utf-8", newline="\n")
    return 0 if verdict == "ALL_CHECKS_PASS" else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(64)
    sys.exit(main(sys.argv[1]))
