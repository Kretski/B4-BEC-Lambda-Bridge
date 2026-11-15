"""B4.12 Q3 - 1+1D perturbative-unitarity derivation for the Proca EFT contact, revision 2.1.

Specification: B4_12_Q3_UNITARITY_SPEC_v0.3.1.md (commit 2a411d0, sha256 0c02465d...).
Conditions:    B4_12_Q3_SPEC_REGISTRATION.md (commit 5169551).
Scope of this revision: B4_12_Q3_CHIRAL_SECTOR_NOTE_v0.2.md (commit ee1e0c7, sha256 5b8ff30d...),
approved for revision 2 with seven conditions:
  1. F2b/F2c: operator commutators in a periodic box, zero-mode terms shown explicitly, nothing assumed zero;
  2. F3: programmatic enumeration of transitions, full contractions and normalisations, L_contact and the
     Lambda term checked separately;
  3. F4: secular term for a concrete normalised wave packet, first order only, no statement on the exact evolution;
  4. no choice among O1/O2/O3;
  5. no Lambda_unit for the full L_contact;
  6. no one-loop computation;
  7. the result for L_RL is conditional.
Normalisation convention (reviewer decision): c = 1, N = 1/(2 s_cm), a = T/(2 s_cm).
Revision 2.1 (static review of revision 2, sha256 27008abe...): F2 commutator test renamed and limited to a
symbolic momentum-space test; F3 adds the full spectator matrix elements R L -> R R R L and R L -> R L L L and
solves the energy balance explicitly; F4 checks that every probability integral is evaluated, finite and nonzero
(otherwise NOT_ESTABLISHED) and derives the relation between the packeted R -> R R R and the full two-particle
transition with a spectator.

Usage (only after authorisation):
    python -I b412_q3_unitarity.py <outdir>

Symbolic only (SymPy). No floating-point evaluation, no registered parameter point, no background.
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from math import factorial
from pathlib import Path

import sympy as sp
from sympy.physics.quantum import Dagger

SPEC = "B4_12_Q3_UNITARITY_SPEC_v0.3.1.md"
REGISTRATION = "B4_12_Q3_SPEC_REGISTRATION.md"
NOTE = "B4_12_Q3_CHIRAL_SECTOR_NOTE_v0.2.md"
REVISION = "2.1"

# ---------------------------------------------------------------- symbols and domain (spec 1.1)
p = sp.Symbol("p", positive=True)                 # domain: p > 0 (boundary p = 0 not reached by substitution)
kap = sp.Symbol("kappa", real=True)               # kappa real, kappa = 0 allowed (Q3-C6)
Lam = sp.Symbol("Lambda", nonnegative=True)       # separate light-theory term (D-a), as in the registered Q2 script
ER, EL = sp.symbols("E_R E_L", positive=True)     # energies of the incoming right- and left-mover
s = sp.Symbol("s_cm", positive=True)              # s_cm = 4 E_R E_L (derived in C1)
C_NORM = sp.Integer(1)                            # c = 1 (reviewer decision); invariance shown with a generic c
t, x = sp.symbols("t x", real=True)
SPECIES = ("f", "g")
CHANNELS = [(a, b) for a in SPECIES for b in SPECIES]   # |a_R b_L>: ff, fg, gf, gg
DENOMINATORS = []                                 # (expression, condition) for every division taken


def div(num, den, cond):
    DENOMINATORS.append((str(den), cond))
    return num / den


def L_contact(ft, fx, gt, gx):
    J1 = ft * fx
    J2 = gt * gx
    return p / 2 * J1**2 - kap * J1 * J2 + div(kap**2, 2 * p, "p > 0 (spec 1.1)") * J2**2


def L_lambda(ft, fx, gt, gx):
    return -Lam / 4 * fx**4


OPERATORS = (("L_contact", L_contact), ("L_lambda", L_lambda))


# ---------------------------------------------------------------- helpers
def int_delta(expr, var, arg, lo, hi):
    """Integrate expr * DiracDelta(arg) over var in (lo, hi) for arg linear in var.

    Returns (value, root, condition that the root lies inside). The Jacobian 1/|d arg/d var| is explicit.
    """
    slope = sp.diff(arg, var)
    assert sp.diff(slope, var) == 0, "argument not linear"
    root = sp.solve(sp.Eq(arg, 0), var)[0]
    inside = sp.And(sp.Gt(root, lo) if lo is not None else True, sp.Lt(root, hi) if hi is not None else True)
    return sp.simplify(expr.subs(var, root) / sp.Abs(slope)), root, sp.simplify(inside)


def leg(species, mover, inout, om=None):
    """Plane-wave data of an external leg: omega, k, sign (+1 incoming, -1 outgoing)."""
    if om is None:
        om = ER if mover == "R" else EL
    k = om if mover == "R" else -om
    return dict(sp=species, om=om, k=k, sg=1 if inout == "in" else -1)


def route_A(L, legs):
    """Route (A): Feynman rule with on-shell external legs, no R/L split.

    Each field is replaced by sum_j e_j exp(i sg_j (k_j x - om_j t)); d_t gives -i sg om, d_x gives +i sg k.
    With <0|f(x)|k> = exp(i k.x) (relativistic normalisation) and S - 1 = i Int L_int at first order, the
    amplitude is the coefficient of e_1 ... e_n (this sums over all assignments of legs to fields).
    """
    es = sp.symbols("e1:%d" % (len(legs) + 1))
    der = {("f", "t"): 0, ("f", "x"): 0, ("g", "t"): 0, ("g", "x"): 0}
    for e, lg in zip(es, legs):
        der[(lg["sp"], "t")] += e * (-sp.I * lg["sg"] * lg["om"])
        der[(lg["sp"], "x")] += e * (sp.I * lg["sg"] * lg["k"])
    expr = sp.expand(L(der[("f", "t")], der[("f", "x")], der[("g", "t")], der[("g", "x")]))
    return sp.expand(sp.Poly(expr, *es).coeff_monomial(sp.Mul(*es)))


# light-cone variables: on shell f = f_R(x - t) + f_L(x + t); R* = d_t f_R (so d_x f_R = -R*), L* = d_t f_L (d_x f_L = +L*)
Rf, Rg, Lf, Lg = sp.symbols("Rf Rg Lf Lg")
P0f, P0g = sp.symbols("P0f P0g")                 # box zero mode of d_t f, d_t g (pi_0 / L_box); absent in the continuum
CHIRAL_VARS = (Rf, Rg, Lf, Lg)


def chiral_decomposition(L, zero_modes=False):
    """Expand L in R/L variables; group monomials by (deg R, deg L) or, with zero modes, by (deg R, deg L, deg P0)."""
    zf, zg = (P0f, P0g) if zero_modes else (0, 0)
    expr = sp.expand(L(zf + Lf + Rf, Lf - Rf, zg + Lg + Rg, Lg - Rg))
    gens = CHIRAL_VARS + ((P0f, P0g) if zero_modes else ())
    poly = sp.Poly(expr, *gens)
    groups = {}
    for mon, coef in poly.terms():
        dR, dL = mon[0] + mon[1], mon[2] + mon[3]
        key = (dR, dL, mon[4] + mon[5]) if zero_modes else (dR, dL)
        groups.setdefault(key, []).append((mon, coef))
    return groups


def fields_of(mon):
    """List of (direction, species) for the oscillator part of a chiral monomial."""
    out = []
    for (d, sp_), n in zip((("R", "f"), ("R", "g"), ("L", "f"), ("L", "g")), mon[:4]):
        out += [(d, sp_)] * n
    return out


def main(outdir: str) -> int:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    T0 = time.time()
    log = []
    checks = {}
    res = {}
    findings = []

    def note(msg):
        line = "[%7.1fs] %s" % (time.time() - T0, msg)
        log.append(line)
        print(line, flush=True)

    note("spec %s; conditions %s; scope %s; revision %s; sympy %s" % (SPEC, REGISTRATION, NOTE, REVISION,
                                                                     sp.__version__))

    # ============================================================ Q3-C1 kinematics
    k3, k4 = sp.symbols("k3 k4", real=True)
    Ein, Pin = ER + EL, ER - EL
    sol_RL = sp.solve([sp.Eq(k3 + (-k4), Ein), sp.Eq(k3 + k4, Pin)], [k3, k4], dict=True)   # k3 > 0 (R), k4 < 0 (L)
    c1_rl = len(sol_RL) == 1 and sp.simplify(sol_RL[0][k3] - ER) == 0 and sp.simplify(sol_RL[0][k4] + EL) == 0
    rr_gap = sp.simplify(Ein - Pin)          # RR final needs this = 0
    ll_gap = sp.simplify(Ein + Pin)          # LL final needs this = 0
    s_cm = sp.expand(Ein**2 - Pin**2)
    c1 = c1_rl and rr_gap == 2 * EL and ll_gap == 2 * ER and sp.simplify(s_cm - 4 * ER * EL) == 0
    checks["Q3_C1_kinematics"] = bool(c1)
    res["C1"] = dict(RL_final_momenta={"k3": str(sol_RL[0][k3]), "k4": str(sol_RL[0][k4])},
                     RR_final_requires_zero=str(rr_gap), LL_final_requires_zero=str(ll_gap), s_cm=str(s_cm))
    note("C1: RL final momenta %s; RR needs %s = 0, LL needs %s = 0; s_cm = %s -> %s"
         % (sol_RL[0], rr_gap, ll_gap, s_cm, "PASS" if c1 else "FAIL"))

    # ============================================================ Q3-C9 channel normalisation (spec 1.2)
    p1, p1p = sp.symbols("p1 p1p", positive=True)    # R momenta (k > 0)
    q2, q2p = sp.symbols("q2 q2p", negative=True)    # L momenta (k < 0)

    def comm(k, kp):                                 # [a(k), a^dag(k')] = 2|k| (2 pi) delta(k - k')
        return 2 * sp.Abs(k) * 2 * sp.pi * sp.DiracDelta(k - kp)

    norm_same = comm(p1p, p1) * comm(q2p, q2) + comm(p1p, q2) * comm(q2p, p1)
    norm_diff = comm(p1p, p1) * comm(q2p, q2)
    exchange = sp.simplify(comm(p1p, q2) * comm(q2p, p1))
    c9_norm = exchange == 0 and sp.simplify(norm_same - norm_diff) == 0
    k1, k2 = sp.symbols("k1 k2", real=True)
    F = sp.exp(-(k1**2 + k2**2)) * (k1 - k2)**2
    reg1 = sp.integrate(F, (k1, 0, sp.oo), (k2, -sp.oo, 0))
    reg2 = sp.integrate(F, (k1, -sp.oo, 0), (k2, 0, sp.oo))
    double_count = sp.simplify(reg1 - reg2) == 0
    c9 = bool(c9_norm and double_count)
    checks["Q3_C9_channel_normalisation"] = c9
    res["C9"] = dict(exchange_term=str(exchange), norm_identical_minus_distinct=str(sp.simplify(norm_same - norm_diff)),
                     regions_equal=bool(double_count),
                     conclusion=("no extra 1/sqrt(2) for f_R f_L or g_R g_L: the exchange term vanishes on disjoint "
                                 "supports, and the 1/2! of the identical-boson completeness relation is cancelled by "
                                 "the two equivalent regions (k1>0,k2<0), (k1<0,k2>0)") if c9 else "NOT ESTABLISHED")
    note("C9: exchange term = %s; regions equal = %s -> %s" % (exchange, double_count, "PASS" if c9 else "FAIL"))

    # ============================================================ elastic fibre (condition 3) and Q3-C3 phase space
    P3, P4 = sp.symbols("P3 P4", real=True)
    J = sp.Abs(sp.Matrix([[sp.diff(P3 - P4, P3), sp.diff(P3 - P4, P4)],
                          [sp.diff(P3 + P4, P3), sp.diff(P3 + P4, P4)]]).det())
    fibre_norm = sp.simplify(2 * ER * 2 * EL * J)
    q = sp.Symbol("q", positive=True)
    E, P = Ein, Pin
    integrand = 1 / (2 * sp.pi * 2 * P3) / (2 * sp.pi * 2 * q) * (2 * sp.pi)**2
    v1, r1, in1 = int_delta(integrand, P3, P3 - q - P, 0, None)
    v1, r2, in2 = int_delta(v1, q, sp.expand((P3 + q - E).subs(P3, r1)), 0, None)
    phi1 = sp.simplify(v1)
    kp, km = sp.symbols("kplus kminus", real=True)
    jac_lc = sp.Abs(sp.Matrix([[sp.Rational(1, 2), sp.Rational(1, 2)], [sp.Rational(1, 2), -sp.Rational(1, 2)]]).det())
    ksq = -kp * km
    w_a = 1 / sp.Abs(sp.diff(ksq, km))
    rs = sp.sqrt(s)
    PmK_sq = sp.expand(-(rs - kp) * (rs - km)).subs(km, 0)
    w_b, kroot, _ = int_delta(w_a.subs(km, 0), kp, PmK_sq, 0, None)
    phi2 = sp.simplify(jac_lc * w_b)
    phi1_s = sp.simplify(phi1.subs(EL, s / (4 * ER)))
    fn = sp.simplify(fibre_norm.subs(EL, s / (4 * ER)))
    c3 = sp.simplify(phi1_s - phi2) == 0 and sp.simplify(phi2 * fn - 1) == 0
    checks["Q3_C3_phase_space"] = bool(c3)
    res["fibre"] = dict(definition=("direct-integral decomposition of the state space over total (E, P); fibre over "
                                    "(E, P) with E > |P|; the two-particle part of the fibre is C^4 over species "
                                    "(channels |a_R b_L>), the fibre inner product is obtained by removing "
                                    "(2 pi)^2 delta^2(P_f - P_i); this C^4 is not the continuum of momentum states"),
                        jacobian_p_to_EP=str(J), fibre_norm=str(fn),
                        root_conditions=[str(in1), str(in2)])
    res["C3"] = dict(method1_direct=str(phi1), method1_in_s=str(phi1_s), method2_invariant_CM=str(phi2),
                     product_with_fibre_norm=str(sp.simplify(phi2 * fn)))
    note("C3: Phi_2 direct = %s = %s; invariant (CM) = %s; fibre norm = %s -> %s"
         % (phi1, phi1_s, phi2, fn, "PASS" if c3 else "FAIL"))

    # ============================================================ route B ingredients
    eps = sp.Symbol("epsilon")
    pf, pg, fx_, gx_, ft_, gt_ = sp.symbols("pi_f pi_g fx gx ft gt", real=True)
    Lint = L_contact(ft_, fx_, gt_, gx_)
    uf = -sp.diff(Lint, ft_).subs({ft_: pf, gt_: pg})
    ug = -sp.diff(Lint, gt_).subs({ft_: pf, gt_: pg})
    Ltot = sp.Rational(1, 2) * (ft_**2 - fx_**2 + gt_**2 - gx_**2) + eps * Lint
    sub = {ft_: pf + eps * uf, gt_: pg + eps * ug}

    def upto1(e_):
        e_ = sp.expand(e_)
        return sp.expand(e_.coeff(eps, 0) + eps * e_.coeff(eps, 1))

    pi_f_check = upto1(sp.diff(Ltot, ft_).subs(sub, simultaneous=True) - pf)
    pi_g_check = upto1(sp.diff(Ltot, gt_).subs(sub, simultaneous=True) - pg)
    H = (pf * ft_ + pg * gt_ - Ltot).subs(sub, simultaneous=True)
    H_target = sp.Rational(1, 2) * (pf**2 + fx_**2 + pg**2 + gx_**2) - eps * Lint.subs({ft_: pf, gt_: pg})
    H_diff = upto1(H - H_target)
    b0 = pi_f_check == 0 and pi_g_check == 0 and H_diff == 0
    res["B0_Hint_equals_minus_Lint_first_order"] = bool(b0)
    res["B0_scope"] = "first order in the couplings only; higher-order Legendre terms are not derived (note F2)"
    note("B0: H_int = -L_int(pi) + O(eps^2): %s" % b0)

    kk = sp.Symbol("kk", positive=True)
    pin, pout = sp.symbols("p_in p_out", positive=True)
    me_ann, _, _ = int_delta(1 / (2 * sp.pi * 2 * kk) * (-sp.I * kk) * 2 * kk * 2 * sp.pi, kk, kk - pin, 0, None)
    me_cre, _, _ = int_delta(1 / (2 * sp.pi * 2 * kk) * (sp.I * kk) * 2 * kk * 2 * sp.pi, kk, kk - pout, 0, None)
    res["B1_one_particle"] = dict(annihilation=str(me_ann), creation=str(me_cre))
    note("B1: <0|d_t f_R|p> = %s * phase;  <p|d_t f_R|0> = %s * phase" % (me_ann, me_cre))

    def me_pair(x_, y_, a_in, c_out, w_in, w_out):
        A = me_ann.subs(pin, w_in)
        C = me_cre.subs(pout, w_out)
        return C * A * (int(x_ == c_out and y_ == a_in) + int(y_ == c_out and x_ == a_in))

    def route_B(L, a, b, cc_, d):
        total = 0
        for mon, coef in chiral_decomposition(L).get((2, 2), []):
            rx = [s_ for d_, s_ in fields_of(mon) if d_ == "R"]
            lx = [s_ for d_, s_ in fields_of(mon) if d_ == "L"]
            total += coef * me_pair(rx[0], rx[1], a, cc_, ER, ER) * me_pair(lx[0], lx[1], b, d, EL, EL)
        return sp.expand(total)

    # ============================================================ amplitude matrix, Q3-C4 two routes
    def T_matrix(L, route):
        Tm = sp.zeros(4, 4)
        for i, (c_, d) in enumerate(CHANNELS):
            for j, (a, b) in enumerate(CHANNELS):
                if route == "A":
                    Tm[i, j] = route_A(L, [leg(a, "R", "in"), leg(b, "L", "in"), leg(c_, "R", "out"), leg(d, "L", "out")])
                else:
                    Tm[i, j] = route_B(L, a, b, c_, d)
        return Tm

    TA = T_matrix(L_contact, "A")
    TB = T_matrix(L_contact, "B")
    TA_lam = T_matrix(L_lambda, "A")
    TB_lam = T_matrix(L_lambda, "B")
    c4 = (sp.simplify(TA - TB) == sp.zeros(4, 4) and sp.simplify(TA_lam - TB_lam) == sp.zeros(4, 4) and b0)
    checks["Q3_C4_two_routes"] = bool(c4)
    res["channels"] = ["%s_R %s_L" % ch for ch in CHANNELS]
    res["T_contact"] = [[str(TA[i, j]) for j in range(4)] for i in range(4)]
    res["T_lambda_separate"] = [[str(TA_lam[i, j]) for j in range(4)] for i in range(4)]
    res["C4_scope"] = ("route A (full Lagrangian, plane waves) equals route B (L_RL part only, Fock matrix elements) "
                       "for the four two-particle channels; this shows that only L_RL contributes to the elastic "
                       "block at first order. Transitions to other states are in F3")
    note("C4: T(A) = T(B) for L_contact and for the Lambda term; B0 = %s -> %s" % (b0, "PASS" if c4 else "FAIL"))
    note("T_contact = %s" % TA)
    note("T_lambda (separate) = %s" % TA_lam)

    TA_s = sp.simplify(TA.subs(EL, s / (4 * ER)))
    TL_s = sp.simplify(TA_lam.subs(EL, s / (4 * ER)))
    boost_inv = all(sp.diff(e, ER) == 0 for e in list(TA_s) + list(TL_s))
    res["T_contact_in_s"] = [[str(TA_s[i, j]) for j in range(4)] for i in range(4)]
    res["T_lambda_in_s"] = [[str(TL_s[i, j]) for j in range(4)] for i in range(4)]
    note("T depends on s_cm only: %s" % boost_inv)

    # ============================================================ dimensionless matrix, c = 1; Q3-C2, C5, C6, C7, C8
    N = div(1, C_NORM * fn, "s_cm > 0")
    a = sp.simplify(N * TA_s)
    a_lam = sp.simplify(N * TL_s)
    S_el = sp.eye(4) + sp.I * C_NORM * a
    cg = sp.Symbol("c_generic", positive=True)
    lamT = sp.Symbol("lambda_T", real=True)
    lam_a_generic = lamT * div(1, cg * fn, "c_generic > 0, s_cm > 0")
    conv_inv = sp.simplify(lam_a_generic * cg - lamT / fn) == 0
    res["normalisation"] = dict(c=str(C_NORM), N=str(N), a="T/(2 s_cm)",
                               convention_independence=("lambda(a_c) = lambda(T)/(c * 2 s_cm) for any c > 0, so "
                                                        "|lambda(a_c)| <= 1/c  <=>  |lambda(T)| <= 2 s_cm: %s" % conv_inv))
    T_free = T_matrix(lambda *z: 0, "A")
    S_free = sp.eye(4) + sp.I * C_NORM * N * T_free
    checks["Q3_C2_free_limit"] = bool(T_free == sp.zeros(4, 4) and S_free == sp.eye(4))
    lam_ = sp.Symbol("lam", positive=True)
    scl = {ER: lam_ * ER, EL: lam_ * EL, p: p / lam_**2, kap: kap / lam_**2, Lam: Lam / lam_**2}
    dimT = sp.simplify(TA.subs(scl, simultaneous=True) - lam_**2 * TA) == sp.zeros(4, 4)
    a_E = sp.simplify(N.subs(s, 4 * ER * EL) * TA)
    dima = sp.simplify(a_E.subs(scl, simultaneous=True) - a_E) == sp.zeros(4, 4)
    checks["Q3_C5_dimension"] = bool(dimT and dima and boost_inv and conv_inv)
    a0 = a.subs(kap, 0)
    f_idx = [CHANNELS.index(("f", "f"))]
    g_idx = [i for i in range(4) if i not in f_idx]
    block = all(a0[i, j] == 0 for i in f_idx for j in g_idx) and all(a0[j, i] == 0 for i in f_idx for j in g_idx)
    gg0 = a0[CHANNELS.index(("g", "g")), CHANNELS.index(("g", "g"))] == 0
    checks["Q3_C6_decoupling"] = bool(block and gg0)
    herm = sp.simplify(a - a.H) == sp.zeros(4, 4)
    checks["Q3_C7_hermiticity"] = bool(herm)
    rem = sp.simplify(S_el.H * S_el - sp.eye(4) - a * a)
    checks["Q3_C8_perturbative_consistency"] = bool(rem == sp.zeros(4, 4) and herm)
    res["a"] = [[str(a[i, j]) for j in range(4)] for i in range(4)]
    res["a_lambda_separate"] = [[str(a_lam[i, j]) for j in range(4)] for i in range(4)]
    res["C8_remainder"] = ("S_el^dag S_el - 1 = a^2 (c = 1; identified, nonzero for a != 0). This is a statement about "
                           "the first-order truncation, not a proof of unitarity")
    for k_ in ("Q3_C2_free_limit", "Q3_C5_dimension", "Q3_C6_decoupling", "Q3_C7_hermiticity",
               "Q3_C8_perturbative_consistency"):
        note("%s: %s" % (k_, "PASS" if checks[k_] else "FAIL"))

    # ============================================================ spec 4.2 / 3.1: exact relation (generic matrices)
    cs = sp.Symbol("c", positive=True)
    xs = sp.Matrix(4, 4, lambda i, j: sp.Symbol("x%d%d" % (i, j), real=True))
    ys = sp.Matrix(4, 4, lambda i, j: sp.Symbol("y%d%d" % (i, j), real=True))
    Ag = xs + sp.I * ys
    Sg = sp.eye(4) + sp.I * cs * Ag
    Rg_ = sp.eye(4) - Sg.H * Sg
    rel_ok = sp.expand(-sp.I * (Ag - Ag.H) - cs * Ag.H * Ag - Rg_ / cs) == sp.zeros(4, 4)
    res["relation_3_1"] = ("-i(A - A^dag) = c A^dag A + R_inel/c holds identically for S_el = 1 + i c A with R_inel "
                           "defined by S_el^dag S_el = 1 - R_inel (an identity for the defined R_inel): %s" % rel_ok)
    note("3.1 exact projected relation: %s" % rel_ok)

    # ============================================================ spec 3.2: order by order, R_inel^(2) kept
    ep = sp.Symbol("varepsilon", positive=True)
    an, A2, R2, R4 = sp.symbols("a A2 R2 R4", commutative=False)
    Aser = ep * an + ep**2 * A2
    Rser = ep**2 * R2 + ep**4 * R4
    ser = sp.expand(-sp.I * (Aser - Dagger(Aser)) - cs * Dagger(Aser) * Aser - Rser / cs).subs(Dagger(an), an)
    o1 = sp.expand(ser.coeff(ep, 1))
    o2 = sp.expand(ser.coeff(ep, 2))
    res["order_1"] = "%s = 0" % o1
    res["order_2_general"] = "%s = 0, i.e. -i(A2 - A2^dag) = c a^2 + R_inel^(2)/c" % o2
    res["order_2_status"] = ("the form without R_inel^(2) is NOT claimed for L_contact: connected 2->4 enters at "
                             "fourth order, but the disconnected collinear splitting of F3/F4 is first order and its "
                             "contribution to R_inel^(2) is secular in the free-Fock construction (F4)")
    note("3.2 order 1: %s = 0; order 2 (general): %s = 0" % (o1, o2))

    V, I_, Ext = sp.symbols("V I E_ext")
    tree = sp.solve([sp.Eq(Ext, 4 * V - 2 * I_), sp.Eq(I_ - V + 1, 0)], [V, I_], dict=True)[0]
    V_2to4 = tree[V].subs(Ext, 6)
    res["coupling_counting"] = dict(external_legs=str(sp.Eq(Ext, 4 * V - 2 * I_)), tree_vertices=str(tree[V]),
                                    connected_2to4_vertices=int(V_2to4),
                                    connected_contribution_to_R_inel_order=int(2 * V_2to4),
                                    disconnected="first-order spectator transitions exist (F3); not covered by this count")
    note("counting: connected 2->4 needs V = %s, R_inel at order %s; disconnected splitting: see F3"
         % (V_2to4, 2 * V_2to4))

    # ============================================================ F1: chiral split, boost weights, parity
    F1 = {}
    for name, L in OPERATORS:
        groups = chiral_decomposition(L)
        F1[name] = {"%d_R_%d_L" % k_: dict(boost_weight=int(k_[0] - k_[1]),
                                          terms=[str(sp.Mul(c_, *[v**e for v, e in zip(CHIRAL_VARS, m)]))
                                                 for m, c_ in v_])
                    for k_, v_ in sorted(groups.items())}
    g_c = chiral_decomposition(L_contact)
    L_RL_expr = sum(c_ * sp.Mul(*[v**e for v, e in zip(CHIRAL_VARS, m)]) for m, c_ in g_c.get((2, 2), []))
    L_ch_expr = sum(c_ * sp.Mul(*[v**e for v, e in zip(CHIRAL_VARS, m)])
                    for k_ in ((4, 0), (0, 4)) for m, c_ in g_c.get(k_, []))
    split_complete = set(g_c) <= {(2, 2), (4, 0), (0, 4)}
    parity = all(sp.expand(L(sp.Symbol("A"), -sp.Symbol("B"), sp.Symbol("C"), -sp.Symbol("D"))
                           - L(sp.Symbol("A"), sp.Symbol("B"), sp.Symbol("C"), sp.Symbol("D"))) == 0
                 for _, L in OPERATORS)
    checks["F1_chiral_split"] = bool(split_complete and parity)
    res["F1"] = dict(groups=F1, L_RL=str(L_RL_expr), L_ch=str(L_ch_expr), split_complete=bool(split_complete),
                     parity_invariant=bool(parity),
                     parity_use="x -> -x maps R <-> L; the L-sector results follow from the R-sector ones")
    note("F1: L_RL = %s" % L_RL_expr)
    note("F1: L_ch = %s" % L_ch_expr)
    note("F1: split complete %s; parity invariant %s -> %s" % (split_complete, parity,
                                                              "PASS" if checks["F1_chiral_split"] else "FAIL"))

    # ============================================================ F2: periodic box, zero modes, operator commutators
    Lb = sp.Symbol("L_box", positive=True)
    q0, pi0 = sp.symbols("q0 pi0", real=True)
    zm = q0 + pi0 * t / Lb                                       # zero-mode part of f in a periodic box (no winding)
    zm_dt, zm_dx = sp.diff(zm, t), sp.diff(zm, x)
    F2 = dict(box=("periodic box of length L_box, modes k_n = 2 pi n / L_box; n > 0: R oscillators, n < 0: L "
                   "oscillators, [a_n, a_m^dag] = delta_nm; field zero mode q0 + pi0 t / L_box, non-compact, no "
                   "winding; normal ordering with respect to the oscillators"),
              zero_mode_field_drops=("derivative couplings only: q0 never appears"),
              zero_mode_in_d_t=str(zm_dt), zero_mode_in_d_x=str(zm_dx))
    # decomposition with the zero mode of the canonical momentum
    zgroups = {}
    for name, L in OPERATORS:
        gz = chiral_decomposition(L, zero_modes=True)
        zgroups[name] = gz
        F2["zero_mode_terms_" + name] = {
            "%d_R_%d_L_%d_P0" % k_: [str(sp.Mul(c_, *[v**e for v, e in zip(CHIRAL_VARS + (P0f, P0g), m)]))
                                    for m, c_ in v_] for k_, v_ in sorted(gz.items()) if k_[2] > 0}
    # the P0-free part coincides with the continuum decomposition
    same = True
    for name, L in OPERATORS:
        cont = {k_: sorted(str(c_) + str(m) for m, c_ in v_) for k_, v_ in chiral_decomposition(L).items()}
        box0 = {k_[:2]: sorted(str(c_) + str(m[:4]) for m, c_ in v_) for k_, v_ in zgroups[name].items() if k_[2] == 0}
        same = same and cont == box0
    F2["P0_free_part_equals_continuum"] = bool(same)

    # [H0, M] for a normal-ordered momentum-conserving oscillator monomial: (sum E_created - sum E_annihilated) M,
    # with momentum conservation sum k_created = sum k_annihilated from Int dx over the box; zero modes carry no energy
    # change ([pi0, H0] = 0). E = k for R oscillators, E = -k for L oscillators.
    def commutes_with_H0(nR, nL):
        rows = []
        for cR in range(nR + 1):
            for cL in range(nL + 1):
                kR = sp.symbols("kR1:%d" % (nR + 1), positive=True)
                kL = sp.symbols("kL1:%d" % (nL + 1), negative=True)
                mom = sum(kR[:cR]) - sum(kR[cR:]) + sum(kL[:cL]) - sum(kL[cL:])
                ener = sum(kR[:cR]) - sum(kR[cR:]) - (sum(kL[:cL]) - sum(kL[cL:]))
                if nR + nL == 0:
                    rows.append(True)
                    continue
                var = (kR + kL)[0]
                sol = sp.solve(sp.Eq(mom, 0), var)
                if not sol:                                     # momentum conservation impossible: term vanishes
                    rows.append(True)
                    continue
                rows.append(sp.simplify(ener.subs(var, sol[0])) == 0)
        return all(rows)

    comm_table = {}
    for name in zgroups:
        for k_ in sorted(zgroups[name]):
            nR, nL, nP = k_
            comm_table["%s %d_R_%d_L_%d_P0" % (name, nR, nL, nP)] = bool(commutes_with_H0(nR, nL))
    F2["commutes_with_H0_by_group"] = comm_table
    qR_ok = comm_table.get("L_contact 4_R_0_L_0_P0", False) and comm_table.get("L_contact 0_R_4_L_0_P0", False)
    rl_not = comm_table.get("L_contact 2_R_2_L_0_P0", True) is False
    # [Q_R, Q_L]: Q_R built from n > 0 oscillators only, Q_L from n < 0 only (definition; no zero mode in either)
    F2["Q_definitions"] = ("Q_R = -Int dx :L_ch^R: from the (4,0,0) group (n > 0 oscillators only); Q_L likewise from "
                           "(0,4,0) (n < 0 only). The zero-mode groups (deg P0 > 0) are not assigned to Q_R or Q_L: "
                           "they contain the shared zero mode and their assignment is not unique; they are listed "
                           "separately")
    F2["Q_R_Q_L_commutator"] = ("zero: normal-ordered products of two disjoint sets of mutually commuting oscillators "
                                "(n > 0 and n < 0); this holds for the definitions above, which exclude the zero mode")
    # Wick reordering of the unordered chiral monomials
    def wick_coeffs(n):
        return {n - 2 * j: sp.Integer(factorial(n)) / (factorial(j) * 2**j * factorial(n - 2 * j)) for j in range(n // 2 + 1)}
    F2["wick_single_field"] = {"X^%d" % n: {":X^%d:" % m: "%s C^%d" % (cf, (n - m) // 2) for m, cf in wick_coeffs(n).items()}
                               for n in (2, 4)}
    F2["wick_remainder_status"] = ("the reordering produces c-numbers and bilinear :R^2: terms times "
                                   "C = <R*^2> = Sum_{n>0} k_n / (2 L_box) (UV-divergent); Int dx :R*^2: is a "
                                   "momentum-conserving bilinear with E = k, so it commutes with H0 (row 2_R_0_L); it "
                                   "rescales the R-sector Hamiltonian and is reported, not dropped")
    # box-normalised scaling of matrix elements: Int dx -> L_box, each oscillator field L_box^(-1/2), each P0 L_box^(-1)
    scal = {}
    for name in zgroups:
        for k_ in sorted(zgroups[name]):
            nR, nL, nP = k_
            scal["%s %d_R_%d_L_%d_P0" % (name, nR, nL, nP)] = str(Lb**(1 - sp.Rational(nR + nL, 2) - nP))
    F2["box_scaling_of_matrix_elements"] = scal
    F2["box_scaling_assumption"] = "fixed pi0 eigenvalue; box-normalised Fock states; leading power only"
    F2["commutator_test_scope"] = ("SYMBOLIC MOMENTUM-SPACE COMMUTATOR TEST: applies [H0, M] = (sum E_created - sum "
                                   "E_annihilated) M to normal-ordered monomials under formal momentum conservation, "
                                   "with continuous signed momentum symbols. It does NOT construct the Fock operators "
                                   "in the periodic box and does NOT impose the discrete lattice k_n = 2 pi n / L_box, "
                                   "n != 0; it is not a full operator construction in the box")
    checks["F2_symbolic_commutator_test"] = bool(same and qR_ok and rl_not)
    res["F2"] = F2
    note("F2: zero mode d_t = %s, d_x = %s; P0-free part = continuum %s; [Q_R,H0] = 0 %s; L_RL does not commute %s"
         % (zm_dt, zm_dx, same, qR_ok, rl_not))
    note("F2: zero-mode groups (L_contact): %s" % F2["zero_mode_terms_L_contact"])

    # ============================================================ F3: programmatic first-order transition table
    # reduction of energy-momentum conservation to separate R and L energy sums (derived, not assumed)
    SRc, SLc, SRa, SLa = sp.symbols("S_R_cre S_L_cre S_R_ann S_L_ann")
    red = sp.solve([sp.Eq(SRc + SLc, SRa + SLa), sp.Eq(SRc - SLc, SRa - SLa)], [SRc, SLc], dict=True)[0]
    res["F3_balance_reduction"] = ("E: S_R_cre + S_L_cre = S_R_ann + S_L_ann; P: S_R_cre - S_L_cre = S_R_ann - S_L_ann "
                                   "(P = +E for R, -E for L)  =>  %s" % red)

    def solve_sum(n, target):
        """Explicitly solve e_1 + ... + e_n = target for e_i > 0.

        Returns (solvable, certificate). Feasible: explicit witness e_i = target/n, checked by substitution and
        positivity. Infeasible: n = 0 with target != 0, or the solved e_1 = target - (e_2 + ... + e_n) is
        nonpositive for all positive e_2..e_n.
        """
        if n == 0:
            ok = sp.simplify(target) == 0
            return ok, ("empty sum = %s" % target) + (" holds" if ok else " cannot hold")
        es_ = sp.symbols("e1:%d" % (n + 1), positive=True)
        e1_sol = sp.solve(sp.Eq(sp.Add(*es_), target), es_[0])
        e1_val = e1_sol[0] if e1_sol else target - sp.Add(*es_[1:])
        wit = {e_: target / n for e_ in es_}
        wit_ok = (sp.simplify(sp.Add(*es_).subs(wit) - target) == 0) and bool((target / n).is_positive)
        if wit_ok:
            return True, "witness e_i = %s" % (target / n)
        infeasible = bool(sp.simplify(e1_val).is_nonpositive)
        return (False, "solved e1 = %s <= 0 for all positive e2..e%d" % (e1_val, n)) if infeasible else \
               (None, "undetermined: e1 = %s" % e1_val)

    def allowed(ann_R, ann_L, cre_R, cre_L):
        """Explicit solution of S_R_cre = E_R * ann_R and S_L_cre = E_L * ann_L with positive created energies."""
        okR, certR = solve_sum(cre_R, ER * ann_R)
        okL, certL = solve_sum(cre_L, EL * ann_L)
        if okR is None or okL is None:
            return None, ["R: " + certR, "L: " + certL]
        return bool(okR and okL), ["R: " + certR, "L: " + certL]

    table = []
    for name, L in OPERATORS:
        for k_, terms in sorted(chiral_decomposition(L).items()):
            dR, dL = k_
            for ann_R in (0, 1):
                for ann_L in (0, 1):
                    if ann_R > dR or ann_L > dL:
                        table.append(dict(operator=name, part="%d_R_%d_L" % k_, annihilated="R" * ann_R + "L" * ann_L,
                                          verdict="impossible", reason="part has no field of that direction"))
                        continue
                    cre_R, cre_L = dR - ann_R, dL - ann_L
                    ok, bad = allowed(ann_R, ann_L, cre_R, cre_L)
                    final = "R" * (cre_R + (1 - ann_R)) + "L" * (cre_L + (1 - ann_L))
                    table.append(dict(operator=name, part="%d_R_%d_L" % k_, annihilated="R" * ann_R + "L" * ann_L or "-",
                                      created="R" * cre_R + "L" * cre_L, final=final,
                                      verdict={True: "allowed", False: "forbidden", None: "UNDETERMINED"}[ok],
                                      certificate=bad))
    allowed_rows = sorted({(r["operator"], r["part"], r["annihilated"], r["final"]) for r in table
                           if r["verdict"] == "allowed"})
    expected = {("L_contact", "2_R_2_L", "RL", "RL"), ("L_contact", "4_R_0_L", "R", "RRRL"),
                ("L_contact", "0_R_4_L", "L", "RLLL")}
    contact_allowed = {r for r in allowed_rows if r[0] == "L_contact"}
    res["F3_table"] = table
    res["F3_allowed"] = [list(r) for r in allowed_rows]
    note("F3: allowed first-order transitions: %s" % allowed_rows)

    # splitting amplitudes R -> R R R by two routes, with state normalisations and symmetry factors
    kv = sp.symbols("k1:4", positive=True)
    Ktot = sp.Add(*kv)

    def route_B_split(L, a, finals):
        total = 0
        for mon, coef in chiral_decomposition(L).get((4, 0), []):
            slots = [s_ for d_, s_ in fields_of(mon)]
            acc = 0
            for j in range(4):
                if slots[j] != a:
                    continue
                rest = slots[:j] + slots[j + 1:]
                for perm in itertools.permutations(range(3)):
                    if all(rest[perm[i]] == finals[i] for i in range(3)):
                        acc += me_ann.subs(pin, Ktot) * sp.Mul(*[me_cre.subs(pout, kv[i]) for i in range(3)])
            total += coef * acc
        return sp.expand(total)

    def sym_factor(finals):
        f_ = 1
        for s_ in SPECIES:
            f_ *= factorial(finals.count(s_))
        return f_

    split = {}
    split_ok = True
    for name, L in OPERATORS:
        for a_ in SPECIES:
            for finals in itertools.combinations_with_replacement(SPECIES, 3):
                legs = [leg(a_, "R", "in", om=Ktot)] + [leg(finals[i], "R", "out", om=kv[i]) for i in range(3)]
                TAs = sp.factor(route_A(L, legs))
                TBs = sp.factor(route_B_split(L, a_, finals))
                eq = sp.simplify(TAs - TBs) == 0
                split_ok = split_ok and eq
                split["%s %s_R -> %s" % (name, a_, "".join(f_ + "_R " for f_ in finals).strip())] = dict(
                    route_A=str(TAs), route_B=str(TBs), equal=bool(eq),
                    final_state_norm_terms=sym_factor(list(finals)),
                    completeness_factor="1/%d" % sym_factor(list(finals)),
                    _T=TAs, _L=L, _a=a_, _finals=finals)
    checks["F3_transitions"] = bool(contact_allowed == expected and split_ok)
    res["F3_splitting"] = {k_: {kk_: vv for kk_, vv in v_.items() if not kk_.startswith("_")}
                           for k_, v_ in split.items()}
    res["F3_normalisation"] = ("one-particle <k|k'> = 2E (2pi) delta(k - k'); three-particle <k1 k2 k3|k1' k2' k3'> = "
                               "sum over species-preserving permutations of products of one-particle factors "
                               "(number of terms = product of m_s!), completeness factor 1/product of m_s!")
    note("F3: splitting amplitudes route A = route B: %s; table matches note v0.2 for L_contact: %s"
         % (split_ok, contact_allowed == expected))

    # ------------------------------------------------ F3 (2.1): full first-order matrix elements with a spectator
    # <out| i Int d^2x :L: |in> by explicit Wick enumeration: every external particle is either attached to a vertex
    # field (matching direction and species) or contracted with an external particle of the other side; all vertex
    # fields are attached (normal ordering). Spacetime integral: phase exp(i(dP x - dE t)); with u = dE + dP,
    # v = dE - dP: (2pi)^2 delta(dE) delta(dP) = 2 (2pi)^2 delta(u) delta(v) if neither vanishes identically;
    # if v == 0 identically (only R attached): Int dt dx exp(i (u/2)(x - t)) = T_int 2pi delta(u/2) = 4pi T_int delta(u);
    # u == 0 identically (only L attached) likewise.
    T_int = sp.Symbol("T_int", positive=True)        # Int dt over the full time axis (secular, formally delta(0))

    def dd(e_):
        return sp.DiracDelta(sp.expand(e_))

    def fock_me(L, ins, outs):
        """First-order <outs| i Int d^2x :L: |ins>. Particles: (direction, species, momentum symbol)."""
        def energy(pt):
            return pt[2] if pt[0] == "R" else -pt[2]
        total = 0
        for mon, coef in [mc for grp in chiral_decomposition(L).values() for mc in grp]:
            slots = fields_of(mon)
            nv = len(slots)
            for na in range(len(ins) + 1):
                for ann in itertools.combinations(range(len(ins)), na):
                    nc = nv - na
                    if nc < 0 or nc > len(outs):
                        continue
                    for cre in itertools.combinations(range(len(outs)), nc):
                        rest_in = [ins[i] for i in range(len(ins)) if i not in ann]
                        rest_out = [outs[i] for i in range(len(outs)) if i not in cre]
                        if len(rest_in) != len(rest_out):
                            continue
                        contr = 0
                        for perm in itertools.permutations(range(len(rest_out))):
                            if all(rest_in[i][:2] == rest_out[perm[i]][:2] for i in range(len(rest_in))):
                                contr += sp.Mul(*[2 * energy(rest_in[i]) * 2 * sp.pi
                                                  * dd(rest_out[perm[i]][2] - rest_in[i][2])
                                                  for i in range(len(rest_in))])
                        if contr == 0:
                            continue
                        attached = [("ann", ins[i]) for i in ann] + [("cre", outs[i]) for i in cre]
                        vert = 0
                        for perm in itertools.permutations(range(nv)):
                            if all(slots[perm[i]] == attached[i][1][:2] for i in range(nv)):
                                vert += sp.Mul(*[(me_ann.subs(pin, energy(pt)) if kind == "ann" else
                                                  me_cre.subs(pout, energy(pt))) for kind, pt in attached])
                        if vert == 0:
                            continue
                        dEv = sum(energy(pt) * (1 if kind == "cre" else -1) for kind, pt in attached)
                        dPv = sum(pt[2] * (1 if kind == "cre" else -1) for kind, pt in attached)
                        u_, v_ = sp.expand(dEv + dPv), sp.expand(dEv - dPv)
                        if v_ == 0 and u_ != 0:
                            stx = 4 * sp.pi * T_int * dd(u_)
                        elif u_ == 0 and v_ != 0:
                            stx = 4 * sp.pi * T_int * dd(v_)
                        elif u_ == 0 and v_ == 0:
                            stx = sp.Symbol("Vol_spacetime")
                        else:
                            stx = 2 * (2 * sp.pi)**2 * dd(u_) * dd(v_)
                        total += sp.I * coef * vert * stx * contr
        return sp.expand(total)

    Kin = sp.Symbol("K", positive=True)              # incoming R momentum
    qin = sp.Symbol("q", negative=True)              # incoming L momentum
    qout = sp.Symbol("q_out", negative=True)         # outgoing L momentum
    pout_R = sp.Symbol("p_out", positive=True)       # outgoing R momentum (L-splitting case)
    lv = sp.symbols("l1:4", negative=True)           # outgoing L momenta (L-splitting case)
    spect = {}
    spect_ok = True
    for name, L in OPERATORS:
        for a_ in SPECIES:
            for b_ in SPECIES:
                for finals in itertools.combinations_with_replacement(SPECIES, 3):
                    for b2 in SPECIES:
                        # R -> R R R with L spectator
                        ins = [("R", a_, Kin), ("L", b_, qin)]
                        outs = [("R", finals[i], kv[i]) for i in range(3)] + [("L", b2, qout)]
                        full = fock_me(L, ins, outs)
                        legs_k = [leg(a_, "R", "in", om=Kin)] + [leg(finals[i], "R", "out", om=kv[i]) for i in range(3)]
                        T13K = route_A(L, legs_k)
                        pred = (sp.I * T13K * 4 * sp.pi * T_int * dd(2 * (Ktot - Kin))       # u = dE + dP for R-only vertex
                                * int(b_ == b2) * 2 * (-qin) * 2 * sp.pi * dd(qout - qin))
                        eqR = sp.simplify(full - sp.expand(pred)) == 0
                        # L -> L L L with R spectator (mirror, computed explicitly)
                        insL = [("R", b2, Kin), ("L", a_, qin)]
                        outsL = [("R", b_, pout_R)] + [("L", finals[i], lv[i]) for i in range(3)]
                        fullL = fock_me(L, insL, outsL)
                        legs_L = [leg(a_, "L", "in", om=-qin)] + [leg(finals[i], "L", "out", om=-lv[i]) for i in range(3)]
                        T13L = route_A(L, legs_L)
                        predL = (sp.I * T13L * 4 * sp.pi * T_int * dd(2 * (qin - sp.Add(*lv)))   # v = dE - dP for L-only vertex
                                 * int(b_ == b2) * 2 * Kin * 2 * sp.pi * dd(pout_R - Kin))
                        eqL = sp.simplify(fullL - sp.expand(predL)) == 0
                        spect_ok = spect_ok and eqR and eqL
                        if full != 0 or fullL != 0 or not (eqR and eqL):
                            spect["%s %s_R %s_L -> %s + %s_L" % (name, a_, b_, " ".join(f_ + "_R" for f_ in finals), b2)] = dict(
                                full_R_split=str(sp.factor(full)), factorised_R=bool(eqR),
                                full_L_split=str(sp.factor(fullL)), factorised_L=bool(eqL))
    checks["F3_spectator_matrix_elements"] = bool(spect_ok)
    res["F3_spectator"] = dict(
        definition=("<c_R d_R e_R b'_L| i Int d^2x :L: |a_R b_L> and the mirror, by explicit Wick enumeration over all "
                    "monomials of the operator (not only the (4,0) part)"),
        factorisation=("full = i T_{a->cde} * 4 pi T_int delta(2(k1+k2+k3-K)) * delta_{b b'} 2|q| (2pi) delta(q'-q): "
                       "%s" % spect_ok),
        spectator_delta=("the spectator contributes 2|q| (2pi) delta(q' - q); the total (2pi)^2 delta^2(P_f - P_i) of "
                         "the four-particle state reduces to delta(0) * delta(dP_R) on the collinear set, represented by "
                         "T_int"),
        four_particle_norm=("<RRR L|RRR L> = (sum over species-preserving permutations of the three R factors) x (L "
                            "factor); completeness factor 1/(product of m_s! over R species)"),
        elements=spect)
    note("F3 (2.1): spectator matrix elements factorise as derived: %s" % spect_ok)

    # ============================================================ F4: secular term, concrete normalised packet, first order
    tp, tt, sig = sp.symbols("t_prime t sigma", positive=True)
    kin = sp.Symbol("k_in", positive=True)
    dE = sp.Add(*kv) - kin                                       # E_f - E_i for an R initial and R R R final state
    on_shell = dE.subs(kin, Ktot)
    step1 = sp.simplify(sp.integrate(sp.exp(sp.I * on_shell * tp), (tp, 0, tt)))
    Dl = sp.Symbol("Delta", real=True, nonzero=True)
    off = sp.simplify(sp.integrate(sp.exp(sp.I * Dl * tp), (tp, 0, tt)))
    off_abs2 = sp.simplify(sp.expand_complex(off * sp.conjugate(off)).rewrite(sp.cos))
    Nn = sp.sqrt(8 * sp.pi / sig)
    phi = Nn * sp.sqrt(kin) * sp.exp(-kin / sig)                 # phi(k) on k > 0
    norm_packet = sp.simplify(sp.integrate(phi**2 / (2 * sp.pi * 2 * kin), (kin, 0, sp.oo)))
    F4 = dict(packet="phi(k) = sqrt(8 pi / sigma) sqrt(k) exp(-k/sigma), k > 0, measure dk/(2pi 2k)",
              packet_norm=str(norm_packet),
              step1_fixed_momenta=("Int_0^t exp(i (E_f - E_i) t') dt' with E_f - E_i = %s identically on the momentum "
                                   "shell -> %s (secular, no oscillation)" % (on_shell, step1)),
              step1_contrast_off_shell="|Int_0^t exp(i Delta t') dt'|^2 = %s (bounded for Delta != 0)" % off_abs2)
    secular = {}
    for key, d_ in split.items():
        if d_["_T"] == 0:
            continue
        name = key.split()[0]
        finals = list(d_["_finals"])
        # vertex with a free incoming energy k_in; momentum conservation is imposed by the delta below
        legs_k = [leg(d_["_a"], "R", "in", om=kin)] + [leg(finals[i], "R", "out", om=kv[i]) for i in range(3)]
        T13k = route_A(d_["_L"], legs_k)
        # step 2: packet amplitude. <f|V|k> = -<f|Int dx L|k> = -T13 (2pi) delta(K - k) at fixed time (H_int = -L_int)
        amp_k = -sp.I * tt * (1 / (2 * sp.pi * 2 * kin)) * phi * (-T13k) * 2 * sp.pi
        amp, root, cond = int_delta(amp_k, kin, kin - Ktot, 0, None)
        amp = sp.simplify(amp)
        # step 3: probability, integrated over the final phase space with the completeness factor
        dens = sp.simplify(amp * sp.conjugate(amp) / sym_factor(finals))
        meas = sp.Mul(*[1 / (2 * sp.pi * 2 * kv[i]) for i in range(3)])
        integrand = sp.expand(sp.simplify(dens * meas))
        prob = integrand
        for i in range(3):
            prob = sp.integrate(prob, (kv[i], 0, sp.oo))
        prob = sp.simplify(prob)
        coef_t2 = sp.simplify(prob / tt**2)
        evaluated = not prob.has(sp.Integral)
        finite = not (prob.has(sp.Piecewise) or prob.has(sp.oo) or prob.has(-sp.oo) or prob.has(sp.zoo)
                      or prob.has(sp.nan))
        t2_form = sp.diff(coef_t2, tt) == 0
        nz = bool(sp.simplify(coef_t2) != 0) and coef_t2.is_zero is not True
        status = "ESTABLISHED" if (evaluated and finite and t2_form and nz) else "NOT_ESTABLISHED"
        secular[key] = dict(packet_amplitude=str(amp), probability_first_order=str(prob),
                            coefficient_of_t2=str(coef_t2), operator=name, evaluated=bool(evaluated),
                            finite=bool(finite), t2_form=bool(t2_form), nonzero=bool(nz), status=status)
        note("F4 %s: P = %s [%s]" % (key, prob, status))
    F4["channels"] = secular
    F4["scope"] = ("first order in the couplings, for this packet only. Not established: behaviour under the exact "
                   "interacting evolution, higher orders, other packets, or the non-existence of any asymptotic "
                   "construction")
    # relation to the full two-particle transition with a spectator (product packet, first order)
    sigL = sp.Symbol("sigma_L", positive=True)
    psiL = sp.sqrt(8 * sp.pi / sigL) * sp.sqrt(-qin) * sp.exp(qin / sigL)          # L packet on q < 0
    normL = sp.simplify(sp.integrate(psiL**2 / (2 * sp.pi * 2 * (-qin)), (qin, -sp.oo, 0)))
    # spectator overlap: Int dq/(2pi 2|q|) psi_L(q) * 2|q| (2pi) delta(q' - q) = psi_L(q')
    ov, _, _ = int_delta(psiL / (2 * sp.pi * 2 * (-qin)) * 2 * (-qin) * 2 * sp.pi, qin, qout - qin, None, 0)
    ov_ok = sp.simplify(ov - psiL.subs(qin, qout)) == 0
    sumL = sp.simplify(sp.integrate((ov * sp.conjugate(ov)) / (2 * sp.pi * 2 * (-qout)), (qout, -sp.oo, 0)))
    F4["spectator_relation"] = dict(
        L_packet="psi_L(q) = sqrt(8 pi / sigma_L) sqrt(-q) exp(q/sigma_L), q < 0", L_packet_norm=str(normL),
        overlap_equals_psi_L=bool(ov_ok), sum_over_final_L=str(sumL),
        statement=("for the normalised product packet psi_R x psi_L, the first-order probability of R L -> R R R L, "
                   "summed over the final L momentum, equals P(R -> R R R) x (sum over final L = %s); the L -> L L L "
                   "channel with an R spectator has different final states and adds separately. This is derived for "
                   "the product packet at first order; it is not assumed for other initial states" % sumL),
        secular_link="the plane-wave T_int of F3 corresponds to the factor t of step 1 for the packet")
    nonzero = [k_ for k_, v_ in secular.items() if v_["operator"] == "L_contact" and v_["status"] == "ESTABLISHED"]
    all_est = all(v_["status"] == "ESTABLISHED" for v_ in secular.values())
    spect_rel = bool(sp.simplify(normL - 1) == 0 and ov_ok and sp.simplify(sumL - 1) == 0)
    checks["F4_secular_first_order"] = bool(sp.simplify(norm_packet - 1) == 0 and step1 == tt and all_est
                                            and len(nonzero) > 0 and spect_rel)
    F4["status_rule"] = ("a channel is ESTABLISHED only if its probability has no unevaluated Integral, no Piecewise, "
                         "no infinity, the form coefficient * t^2 and a nonzero coefficient; otherwise NOT_ESTABLISHED "
                         "and the check fails")
    res["F4"] = F4
    note("F4: packet norm = %s; fixed-momenta time integral = %s; nonzero secular channels (L_contact): %s"
         % (norm_packet, step1, nonzero))

    findings.append(dict(
        id="CHIRAL_COLLINEAR_SPLITTING",
        statement=("first-order transitions R L -> R R R L and R L -> R L L L exist through the chiral part L_ch (F3); "
                   "in first-order perturbation theory around free Fock states they give a secular term (F4) in "
                   "channels %s" % nonzero),
        consequence=("the standard perturbative construction around free massless Fock states cannot be completed "
                     "for L_contact; within Q3 the full-operator outcome is NO_UNIQUE_CRITERION for the present "
                     "setting. This is not a statement that no unitary theory or scattering description exists"),
        status="DERIVED AT FIRST ORDER; no option O1/O2/O3 selected (note v0.2, condition 4)"))

    # ============================================================ spec 3.3 criterion: outcomes
    xx, yy = sp.symbols("xx yy", real=True)
    ymin = sp.solve(sp.diff(cs**2 * (xx**2 + yy**2) - 2 * cs * yy, yy), yy)[0]
    proj = sp.solve(sp.Eq((cs**2 * (xx**2 + yy**2) - 2 * cs * yy).subs(yy, ymin), 0), xx)
    res["criterion_steps"] = dict(
        i="contraction: |<v, S_el v>| <= 1 for unit v (exact, no normality needed)",
        ii="S_el = 1 + i c A: |Re<v,A v>| <= 1/|c|, i.e. ||H(A)||_op <= 1/|c| (exact); projection endpoints %s"
           % [str(r) for r in proj],
        iii="H(A) = a + O(g^2) => leading-order estimate |lambda_i(a)| <~ 1 for c = 1 (not an exact inequality)",
        presupposes="a well-defined elastic block S_el; blocked for L_contact in the free-Fock setting (F4)")
    res["outcome_full_L_contact"] = ("NO_UNIQUE_CRITERION (present free-Fock setting; F3/F4)" if nonzero else
                                     "F4 NOT ESTABLISHED; outcome undetermined")
    res["Lambda_unit_full_L_contact"] = "NOT COMPUTED (condition 5)"

    # ============================================================ conditional result for L_RL (condition 7)
    diag = a.is_diagonal()
    eig = [a[i, i] for i in range(4)] if diag else list(a.eigenvals(multiple=True))
    eig_T = [sp.simplify(e_ * C_NORM * fn) for e_ in eig]
    cond = dict(label=("CONDITIONAL: elastic block of L_RL only. Valid only if a prescription of type O1 or O2 is "
                       "later shown to leave the first-order RL amplitude unchanged. Not a result for the full "
                       "contact operator and not a result for the full theory"),
                a_is_diagonal=bool(diag), rank_a=int(sp.Matrix(a).rank(simplify=True)),
                eigenvalues_a=[str(e_) for e_ in eig], eigenvalues_T=[str(e_) for e_ in eig_T])
    thr = []
    for e_ in eig_T:
        alpha = sp.simplify(e_ / s**2)
        if alpha == 0:
            thr.append(dict(eig=str(e_), s_cross="none (eigenvalue zero)"))
        elif sp.diff(alpha, s) != 0:
            thr.append(dict(eig=str(e_), s_cross="NOT OF FORM alpha*s^2; not solved"))
        else:
            s_i = sp.simplify(div(fn / s, sp.Abs(alpha), "alpha != 0 (kappa != 0 for kappa-dependent entries)"))
            thr.append(dict(eig=str(e_), alpha=str(alpha), s_cross=str(s_i), sqrt_s_cross=str(sp.sqrt(s_i)), _s=s_i))
    finite = [d_["_s"] for d_ in thr if "_s" in d_]
    lam_RL = sp.sqrt(sp.Min(*finite)) if finite else None
    cond["thresholds"] = [{k_: v_ for k_, v_ in d_.items() if k_ != "_s"} for d_ in thr]
    cond["Lambda_RL_conditional"] = str(lam_RL) if lam_RL is not None else "no finite crossing"
    cond["qualifier"] = ("leading-order estimate up to O(g^2); not an exact bound, not an exact violation of "
                         "unitarity, not the scale of new physics; kappa-dependent thresholds for kappa != 0 only")
    cond["diagnostic_times_sqrt_p"] = str(sp.simplify(lam_RL * sp.sqrt(p))) if lam_RL is not None else "n/a"
    res["L_RL_conditional"] = cond
    res["D_a_lambda_subsection"] = dict(T_lambda=res["T_lambda_in_s"], a_lambda=res["a_lambda_separate"],
                                        chiral_parts=res["F1"]["groups"]["L_lambda"],
                                        note=("reported separately; the Lambda term also has chiral parts, so no "
                                              "threshold is computed for it"))
    note("L_RL conditional: a diagonal %s; eigenvalues of T %s; Lambda_RL (conditional) = %s"
         % (diag, eig_T, cond["Lambda_RL_conditional"]))

    # ============================================================ summary
    mandatory = ["Q3_C1_kinematics", "Q3_C2_free_limit", "Q3_C3_phase_space", "Q3_C4_two_routes",
                 "Q3_C5_dimension", "Q3_C6_decoupling", "Q3_C7_hermiticity", "Q3_C8_perturbative_consistency",
                 "Q3_C9_channel_normalisation", "F1_chiral_split", "F2_symbolic_commutator_test",
                 "F3_transitions", "F3_spectator_matrix_elements", "F4_secular_first_order"]
    verdict = "ALL_CHECKS_PASS" if all(checks[k_] for k_ in mandatory) else "CHECKS_FAILED"
    res["denominators_checked"] = sorted(set(DENOMINATORS))
    summary = dict(spec=SPEC, registration=REGISTRATION, scope_note=NOTE, revision=REVISION, sympy=sp.__version__,
                   seconds=round(time.time() - T0, 1), checks={k_: bool(v_) for k_, v_ in checks.items()},
                   mandatory=mandatory, results=res, findings=findings, verdict=verdict,
                   scope=("ALL_CHECKS_PASS covers the listed checks only; it is not an independent proof, it does "
                          "not select O1/O2/O3, and the L_RL result is conditional"))
    (out / "b412_q3_results.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8", newline="\n")
    for k_ in mandatory:
        note("%s %s" % (k_, "PASS" if checks[k_] else "FAIL"))
    note("findings: %s" % [f_["id"] for f_ in findings])
    note(verdict)
    (out / "b412_q3_output.txt").write_text("\n".join(log) + "\n", encoding="utf-8", newline="\n")
    return 0 if verdict == "ALL_CHECKS_PASS" else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(64)
    sys.exit(main(sys.argv[1]))
