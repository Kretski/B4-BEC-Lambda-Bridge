"""B4.12 Q3 - zero modes and discrete operator construction in a periodic box, revision 1.1.

Revision 1.1 (static review of revision 1, sha256 3280f9b0...): ZC6 now also requires scal_ok (Z5 L-scaling
equals the registered table). No other change.

Specification: B4_12_Q3_ZERO_MODE_SPEC_v0.2.md (commit fcd7ace, sha256 2a198f7f1e74360c...).
Registered inputs referred to (not read at run time; the needed registered expressions are embedded below):
  - B4_12_Q3_UNITARITY_SPEC_v0.3.1.md (commit 2a411d0): operators, domain p > 0, kappa real, Lambda term separate (D-a);
  - run df4d8b9, results_b412_q3/b412_q3_results.json (sha256 8bbdd1d0...): F1 split, F2 zero-mode table,
    F2 commutes_with_H0_by_group, F2 box scaling, T_contact, T_lambda, F3 splitting amplitudes (route B).

Tasks (spec section 3): Z1 operator form of H_int on the lattice; Z2 zero-mode momenta separately; Z3 exact
lattice commutators on (N_R, N_L) sectors with a cutoff-independence test; Z4 classification of the pi-dependent
groups with the full free energy E0; Z5 zero-mode vacuum limit with commensurate external momenta; Z6 Wick
remainders in closed form. Checks ZC1-ZC7 (spec section 5).

Two routes (spec section 4, D-z3):
  route 1 - symbolic mode algebra: chiral R/L/pi decomposition, explicit integer labels, Kronecker momentum
            conservation, commutators by normal-ordering algebra of monomials;
  route 2 - explicit finite matrices: a and a^dag act on an occupation-number basis (sqrt factors); the
            interaction is built from the derivative fields d_t phi, d_x phi themselves (no R/L rewrite), with all
            creation parts left of all annihilation parts and the x-integral by lattice orthogonality. Route 2
            never reads the monomial table, coefficients or intermediate expressions of route 1. Any state reached
            outside the constructed basis is reported as NOT_COVERED, never as zero.

Order of limits (spec section 2, D-z2): all primary results are exact at finite L and finite N_max; only then
N_max -> infinity at fixed L where a limit exists; then L -> infinity at fixed pi_{f,0}, pi_{g,0}. The reverse
order is not used. Divergent quantities (Z6) are reported at finite N_max and L only.

Out of scope (spec section 7): O1/O2/O3, [Q_R, L_RL], probabilities and secular coefficients, the limit at fixed
pi/L = w, renormalisation, Lambda_unit, numerical evaluation on registered parameter points, T7a, Zenodo.

Usage (only after authorisation):
    python3 -I b412_q3_zero_modes.py <outdir>

Symbolic only (SymPy); concrete integers only for the labels (N_R, N_L, N_max). No floating-point evaluation.
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from collections import Counter, defaultdict
from math import comb, factorial
from pathlib import Path

import sympy as sp

SPEC = "B4_12_Q3_ZERO_MODE_SPEC_v0.2.md"
SPEC_COMMIT = "fcd7ace"
REGISTERED_RESULTS = "results_b412_q3/b412_q3_results.json (run df4d8b9, sha256 8bbdd1d0...)"
REVISION = "1.1"

# ---------------------------------------------------------------- run parameters (integer labels only)
N_LIST = (1, 2, 3)                                   # cutoffs N_max for Z3/Z4/ZC5
SECTORS = [(a, b) for a in range(3) for b in range(3)]   # initial sectors (N_R, N_L), 0 <= N_R, N_L <= 2
ZC1_N_LIST = (1, 2, 3, 4)                            # cutoffs for the mode-algebra check ZC1
Z6_VERIFY_N = (1, 2)                                 # cutoffs for the explicit Wick verification (route 1)

# ---------------------------------------------------------------- symbols and domain
p = sp.Symbol("p", positive=True)                    # p > 0 (spec v0.3.1 section 1.1)
kap = sp.Symbol("kappa", real=True)                  # kappa real
Lam = sp.Symbol("Lambda", nonnegative=True)          # separate light-theory term (D-a)
Lb = sp.Symbol("L", positive=True)                   # box length L
Pi_f, Pi_g = sp.symbols("pi_f0 pi_g0", real=True)    # zero-mode momenta pi_{f,0}, pi_{g,0}
PI = {"f": Pi_f, "g": Pi_g}
x, y, t = sp.symbols("x y t", real=True)
Nsym = sp.Symbol("N_max", positive=True, integer=True)
FIELDS = ("f", "g")
DENOMINATORS = []


def div(num, den, cond):
    DENOMINATORS.append((str(den), cond))
    return num / den


def L_contact(ft, fx, gt, gx):
    """Registered primary operator (spec v0.3.1 section 2)."""
    J1 = ft * fx
    J2 = gt * gx
    return p / 2 * J1**2 - kap * J1 * J2 + div(kap**2, 2 * p, "p > 0 (spec v0.3.1 1.1)") * J2**2


def L_lambda(ft, fx, gt, gx):
    """Separate Lambda term (decision D-a)."""
    return -Lam / 4 * fx**4


OPERATORS = (("L_contact", L_contact), ("L_lambda", L_lambda))
OPFUN = dict(OPERATORS)

# ---------------------------------------------------------------- registered expressions (embedded, run df4d8b9)
REG_F1 = {
    "L_contact": {
        (0, 4): ["Lf**4*p/2", "-Lf**2*Lg**2*kappa", "Lg**4*kappa**2/(2*p)"],
        (2, 2): ["-Lf**2*Rf**2*p", "Lg**2*Rf**2*kappa", "Lf**2*Rg**2*kappa", "-Lg**2*Rg**2*kappa**2/p"],
        (4, 0): ["Rf**4*p/2", "-Rf**2*Rg**2*kappa", "Rg**4*kappa**2/(2*p)"],
    },
    "L_lambda": {
        (0, 4): ["-Lambda*Lf**4/4"], (1, 3): ["Lambda*Lf**3*Rf"], (2, 2): ["-3*Lambda*Lf**2*Rf**2/2"],
        (3, 1): ["Lambda*Lf*Rf**3"], (4, 0): ["-Lambda*Rf**4/4"],
    },
}
REG_F2_ZERO_MODE = {   # registered F2 zero_mode_terms_L_contact, key (deg R, deg L, total deg P0); P0 = pi_0 / L
    (0, 2, 2): ["Lf**2*P0f**2*p/2", "-Lf*Lg*P0f*P0g*kappa", "Lg**2*P0g**2*kappa**2/(2*p)"],
    (0, 3, 1): ["Lf**3*P0f*p", "-Lf**2*Lg*P0g*kappa", "-Lf*Lg**2*P0f*kappa", "Lg**3*P0g*kappa**2/p"],
    (1, 1, 2): ["-Lf*P0f**2*Rf*p", "Lg*P0f*P0g*Rf*kappa", "Lf*P0f*P0g*Rg*kappa", "-Lg*P0g**2*Rg*kappa**2/p"],
    (1, 2, 1): ["-Lf**2*P0f*Rf*p", "Lg**2*P0f*Rf*kappa", "Lf**2*P0g*Rg*kappa", "-Lg**2*P0g*Rg*kappa**2/p"],
    (2, 0, 2): ["P0f**2*Rf**2*p/2", "-P0f*P0g*Rf*Rg*kappa", "P0g**2*Rg**2*kappa**2/(2*p)"],
    (2, 1, 1): ["-Lf*P0f*Rf**2*p", "Lg*P0g*Rf**2*kappa", "Lf*P0f*Rg**2*kappa", "-Lg*P0g*Rg**2*kappa**2/p"],
    (3, 0, 1): ["P0f*Rf**3*p", "-P0g*Rf**2*Rg*kappa", "-P0f*Rf*Rg**2*kappa", "P0g*Rg**3*kappa**2/p"],
}
REG_F2_ZERO_MODE_LAMBDA = {}   # registered: empty
REG_F2_COMMUTES = {
    "L_contact 0_R_2_L_2_P0": True, "L_contact 0_R_3_L_1_P0": True, "L_contact 0_R_4_L_0_P0": True,
    "L_contact 1_R_1_L_2_P0": False, "L_contact 1_R_2_L_1_P0": False, "L_contact 2_R_0_L_2_P0": True,
    "L_contact 2_R_1_L_1_P0": False, "L_contact 2_R_2_L_0_P0": False, "L_contact 3_R_0_L_1_P0": True,
    "L_contact 4_R_0_L_0_P0": True, "L_lambda 0_R_4_L_0_P0": True, "L_lambda 1_R_3_L_0_P0": False,
    "L_lambda 2_R_2_L_0_P0": False, "L_lambda 3_R_1_L_0_P0": False, "L_lambda 4_R_0_L_0_P0": True,
}
REG_F2_SCALING = {
    "L_contact 0_R_2_L_2_P0": "L_box**(-2)", "L_contact 0_R_3_L_1_P0": "L_box**(-3/2)",
    "L_contact 0_R_4_L_0_P0": "1/L_box", "L_contact 1_R_1_L_2_P0": "L_box**(-2)",
    "L_contact 1_R_2_L_1_P0": "L_box**(-3/2)", "L_contact 2_R_0_L_2_P0": "L_box**(-2)",
    "L_contact 2_R_1_L_1_P0": "L_box**(-3/2)", "L_contact 2_R_2_L_0_P0": "1/L_box",
    "L_contact 3_R_0_L_1_P0": "L_box**(-3/2)", "L_contact 4_R_0_L_0_P0": "1/L_box",
    "L_lambda 0_R_4_L_0_P0": "1/L_box", "L_lambda 1_R_3_L_0_P0": "1/L_box", "L_lambda 2_R_2_L_0_P0": "1/L_box",
    "L_lambda 3_R_1_L_0_P0": "1/L_box", "L_lambda 4_R_0_L_0_P0": "1/L_box",
}
REG_CHANNELS = [("f", "f"), ("f", "g"), ("g", "f"), ("g", "g")]   # |a_R b_L>
REG_T = {   # registered T[i][j], i = out channel, j = in channel (E_R, E_L of the incoming movers)
    "L_contact": [["-4*E_L**2*E_R**2*p", "0", "0", "0"], ["0", "4*E_L**2*E_R**2*kappa", "0", "0"],
                  ["0", "0", "4*E_L**2*E_R**2*kappa", "0"], ["0", "0", "0", "-4*E_L**2*E_R**2*kappa**2/p"]],
    "L_lambda": [["-6*E_L**2*E_R**2*Lambda", "0", "0", "0"], ["0", "0", "0", "0"], ["0", "0", "0", "0"],
                 ["0", "0", "0", "0"]],
}
REG_SPLIT = {   # registered F3 route B: a_R(k1+k2+k3) -> c1_R(k1) c2_R(k2) c3_R(k3)
    ("L_contact", "f", ("f", "f", "f")): "-12*k1*k2*k3*p*(k1 + k2 + k3)",
    ("L_contact", "f", ("f", "f", "g")): "0",
    ("L_contact", "f", ("f", "g", "g")): "4*k1*k2*k3*kappa*(k1 + k2 + k3)",
    ("L_contact", "f", ("g", "g", "g")): "0",
    ("L_contact", "g", ("f", "f", "f")): "0",
    ("L_contact", "g", ("f", "f", "g")): "4*k1*k2*k3*kappa*(k1 + k2 + k3)",
    ("L_contact", "g", ("f", "g", "g")): "0",
    ("L_contact", "g", ("g", "g", "g")): "-12*k1*k2*k3*kappa**2*(k1 + k2 + k3)/p",
    ("L_lambda", "f", ("f", "f", "f")): "6*Lambda*k1*k2*k3*(k1 + k2 + k3)",
    ("L_lambda", "f", ("f", "f", "g")): "0", ("L_lambda", "f", ("f", "g", "g")): "0",
    ("L_lambda", "f", ("g", "g", "g")): "0", ("L_lambda", "g", ("f", "f", "f")): "0",
    ("L_lambda", "g", ("f", "f", "g")): "0", ("L_lambda", "g", ("f", "g", "g")): "0",
    ("L_lambda", "g", ("g", "g", "g")): "0",
}

# chiral variables (route 1): R_phi, L_phi are the n > 0 and n < 0 parts of d_t phi; P0_phi = pi_{phi,0} / L
Rf, Rg, Lf, Lg, P0f, P0g = sp.symbols("Rf Rg Lf Lg P0f P0g")
CH = (Rf, Rg, Lf, Lg, P0f, P0g)
FAM = (("f", "R"), ("g", "R"), ("f", "L"), ("g", "L"))       # order of the first four chiral variables
P0_TO_PI = {P0f: Pi_f / Lb, P0g: Pi_g / Lb}
REG_LOCALS = dict(p=p, kappa=kap, Lambda=Lam, Rf=Rf, Rg=Rg, Lf=Lf, Lg=Lg, P0f=P0f, P0g=P0g,
                  L_box=Lb, E_R=sp.Symbol("E_R", positive=True), E_L=sp.Symbol("E_L", positive=True),
                  k1=sp.Symbol("k1", positive=True), k2=sp.Symbol("k2", positive=True),
                  k3=sp.Symbol("k3", positive=True))

# derivative-field symbols (route 2): the interaction is expanded in these, not in R/L
FT, FX, GT, GX = sp.symbols("f_t f_x g_t g_x")
FIELD_FACTORS = ((FT, ("f", "t")), (FX, ("f", "x")), (GT, ("g", "t")), (GX, ("g", "x")))


# ---------------------------------------------------------------- generic helpers
def is_zero(e):
    e = sp.expand(e)
    if e == 0:
        return True
    return sp.simplify(e) == 0


def direction(n):
    return "R" if n > 0 else "L"


def canon(occ):
    return tuple(sorted((m, o) for m, o in occ.items() if o > 0))


def sector_of(st):
    nr = sum(m[1] * o for m, o in st if m[1] > 0)
    nl = sum(-m[1] * o for m, o in st if m[1] < 0)
    return nr, nl


def level(st):
    """Largest mode label that a state of this sector may contain."""
    return max(sector_of(st))


def state_str(st):
    return "|" + " ".join("%s%+d^%d" % (m[0], m[1], o) for m, o in st) + ">" if st else "|0>"


def mono_str(cre, ann):
    return " ".join(["a+(%s%+d)" % m for m in cre] + ["a(%s%+d)" % m for m in ann]) or "1"


def tag_str(name, tag):
    return "%s %d_R_%d_L_%d_Pf_%d_Pg" % (name, *tag)


def tag3_str(name, tag):
    return "%s %d_R_%d_L_%d_P0" % (name, tag[0], tag[1], tag[2] + tag[3])


def states_in_sector(NR, NL, N):
    """All occupation states over modes 1 <= |n| <= N of f and g with the given (N_R, N_L)."""
    def fill(modes, target):
        if not modes:
            if target == 0:
                yield {}
            return
        m = modes[0]
        w = abs(m[1])
        for o in range(target // w + 1):
            for rest in fill(modes[1:], target - o * w):
                d = dict(rest)
                if o:
                    d[m] = o
                yield d
    R = [(f, n) for f in FIELDS for n in range(1, N + 1)]
    Lm = [(f, -n) for f in FIELDS for n in range(1, N + 1)]
    out = []
    for a in fill(R, NR):
        for b in fill(Lm, NL):
            out.append(canon({**a, **b}))
    return out


def clean(d):
    out = {}
    for k_, v_ in d.items():
        v_ = sp.expand(v_)
        if v_ != 0 and not is_zero(v_):
            out[k_] = v_
    return out


def vec_diff_zero(u, v):
    return all(is_zero(u.get(k_, 0) - v.get(k_, 0)) for k_ in set(u) | set(v))


# ================================================================ ZC2: derivative fields from the mode expansion
def derive_fields():
    """Differentiate phi = q0 + pi0 t / L + sum_n (2|k_n| L)^(-1/2) [a_n e^{i(k_n x - |k_n| t)} + h.c.].

    Returns the sign table SGN[(deriv, direction, part)] defined by
        coefficient of a_n (resp. a_n^dag) in d_deriv phi at x = t = 0  =  SGN * i * sqrt(pi |n|) / L,
    the zero-mode parts, the Fourier keys and the flags of ZC2. Everything is derived symbolically.
    """
    nu = sp.Symbol("nu", positive=True, integer=True)
    unit = sp.I * sp.sqrt(sp.pi * nu) / Lb
    a_op, ad_op = sp.symbols("a_op adag_op")
    sgn = {}
    flags = {}
    detail = {}
    for fld in FIELDS:
        q0, P = sp.symbols("q_%s0 Pz_%s" % (fld, fld), real=True)
        for dr, nval in (("R", nu), ("L", -nu)):
            k = 2 * sp.pi * nval / Lb
            om = sp.Abs(k)
            ph = sp.exp(sp.I * (k * x - om * t))
            phc = sp.exp(-sp.I * (k * x - om * t))
            u = (a_op * ph + ad_op * phc) / sp.sqrt(2 * om * Lb)
            phi = q0 + P * t / Lb + u
            d = {"t": sp.expand(sp.diff(phi, t)), "x": sp.expand(sp.diff(phi, x))}
            for dv in "tx":
                e = d[dv]
                zero_part = sp.simplify(e.subs({a_op: 0, ad_op: 0}))
                want = P / Lb if dv == "t" else 0
                flags["%s_%s_%s_zero_part" % (fld, dr, dv)] = bool(is_zero(zero_part - want))
                flags["%s_%s_%s_q0_absent" % (fld, dr, dv)] = q0 not in e.free_symbols
                for part, op_, expo in (("ann", a_op, ph), ("cre", ad_op, phc)):
                    cfull = e.coeff(op_)
                    c0 = sp.simplify(cfull.subs({x: 0, t: 0}))
                    xt_ok = is_zero(sp.powsimp(cfull - c0 * expo, combine="exp"))
                    sg = sp.simplify(c0 / unit)
                    ok = xt_ok and sg in (1, -1)
                    flags["%s_%s_%s_%s_form" % (fld, dr, dv, part)] = bool(ok)
                    key = (dv, dr, part)
                    if ok:
                        if key in sgn and sgn[key] != int(sg):
                            flags["%s_%s_%s_%s_f_g_equal" % (fld, dr, dv, part)] = False
                        sgn[key] = int(sg)
                    detail["%s %s d_%s %s" % (fld, dr, dv, part)] = str(c0)
                # Fourier key at t = 0: a_n carries exp(+i k_n x), a_n^dag carries exp(-i k_n x)
                flags["%s_%s_fourier_keys" % (fld, dr)] = bool(
                    is_zero(ph.subs(t, 0) - sp.exp(sp.I * 2 * sp.pi * nval * x / Lb)))
            osc_t = d["t"] - d["t"].subs({a_op: 0, ad_op: 0})
            osc_x = d["x"] - d["x"].subs({a_op: 0, ad_op: 0})
            if dr == "R":
                flags["%s_dx_R_equals_minus_dt_R" % fld] = bool(is_zero(osc_x + osc_t))
            else:
                flags["%s_dx_L_equals_plus_dt_L" % fld] = bool(is_zero(osc_x - osc_t))
    return sgn, flags, detail


# ================================================================ route 1: symbolic mode algebra
def c1(sgn, dr, part, n):
    """Route-1 mode coefficient of R_phi / L_phi (the d_t parts): SGN[t, dir, part] * i * sqrt(pi |n|) / L."""
    return sgn[("t", dr, part)] * sp.I * sp.sqrt(sp.pi * sp.Abs(n)) / Lb


def chiral_groups(Lfun, dt_rep, dx_rep):
    """Group the monomials of L_int(d_t f, d_x f, d_t g, d_x g) by (deg R, deg L, deg P0f, deg P0g)."""
    expr = sp.expand(Lfun(dt_rep["f"], dx_rep["f"], dt_rep["g"], dx_rep["g"]))
    groups = defaultdict(list)
    if expr == 0:
        return {}
    for mon, coef in sp.Poly(expr, *CH).terms():
        groups[(mon[0] + mon[1], mon[2] + mon[3], mon[4], mon[5])].append((mon, coef))
    return dict(groups)


def mono_expr(mon, coef):
    return sp.Mul(coef, *[v**e for v, e in zip(CH, mon)])


def r1_build(sgn, terms, N, sign=-1):
    """sign * Int_0^L dx :prod: for chiral monomials, as normal-ordered mode monomials.

    Key: (sorted creation modes, sorted annihilation modes); a mode is (field, n). Momentum conservation is the
    Kronecker delta of spec section 2: annihilation +n, creation -n, zero modes 0; the integral gives L.
    """
    C1 = {}
    out = defaultdict(lambda: sp.Integer(0))
    for mon, coef in terms:
        factors = []
        for fam, e in zip(FAM, mon[:4]):
            factors += [fam] * e
        pref = sign * coef * Lb * (Pi_f / Lb)**mon[4] * (Pi_g / Lb)**mon[5]
        lists = []
        for fld, dr in factors:
            s_ = 1 if dr == "R" else -1
            lists.append([(part, (fld, s_ * n)) for n in range(1, N + 1) for part in ("cre", "ann")])
        for combo in itertools.product(*lists):
            if sum(m[1] if part == "ann" else -m[1] for part, m in combo) != 0:
                continue
            c = pref
            for part, m in combo:
                key = (part, m[1])
                if key not in C1:
                    C1[key] = c1(sgn, direction(m[1]), part, m[1])
                c = c * C1[key]
            cre = tuple(sorted(m for part, m in combo if part == "cre"))
            ann = tuple(sorted(m for part, m in combo if part == "ann"))
            out[(cre, ann)] += c
    return clean(out)


def r1_H0(N):
    """Oscillator part of H0 = sum |k_n| a^dag a (normal-ordered, E_c = 0). The pi^2/(2L) part is treated in Z2."""
    return {(((f, n),), ((f, n),)): 2 * sp.pi * abs(n) / Lb
            for f in FIELDS for n in list(range(1, N + 1)) + list(range(-N, 0))}


def no_product(A, B):
    """Normal-ordered product of two normal-ordered operators (bosonic Wick algebra, mode by mode):
    a^al (a^dag)^be = sum_j j! C(al, j) C(be, j) (a^dag)^(be-j) a^(al-j)."""
    out = defaultdict(lambda: sp.Integer(0))
    for (c1_, a1), x1 in A.items():
        ca1 = Counter(a1)
        for (c2, a2), x2 in B.items():
            cc2 = Counter(c2)
            common = sorted(m for m in ca1 if m in cc2)
            for ks in itertools.product(*[range(min(ca1[m], cc2[m]) + 1) for m in common]):
                w = 1
                ca = Counter(a1)
                cc = Counter(c2)
                for m, j in zip(common, ks):
                    w *= factorial(j) * comb(ca1[m], j) * comb(cc2[m], j)
                    ca[m] -= j
                    cc[m] -= j
                cre = tuple(sorted(list(c1_) + list(cc.elements())))
                ann = tuple(sorted(list(ca.elements()) + list(a2)))
                out[(cre, ann)] += w * x1 * x2
    return out


def r1_commutator(A, B):
    ab = no_product(A, B)
    ba = no_product(B, A)
    return clean({k_: ab.get(k_, 0) - ba.get(k_, 0) for k_ in set(ab) | set(ba)})


def r1_apply(op, st):
    """Action of a normal-ordered operator on an occupation state (route 1)."""
    out = defaultdict(lambda: sp.Integer(0))
    for (cre, ann), c in op.items():
        occ = dict(st)
        rad = 1
        ok = True
        for m, k_ in Counter(ann).items():
            o = occ.get(m, 0)
            if o < k_:
                ok = False
                break
            rad *= factorial(o) // factorial(o - k_)
            occ[m] = o - k_
        if not ok:
            continue
        for m, k_ in Counter(cre).items():
            o = occ.get(m, 0)
            rad *= factorial(o + k_) // factorial(o)
            occ[m] = o + k_
        out[canon(occ)] += c * sp.sqrt(rad)
    return clean(out)


def energy_change(cre, ann):
    return sum(2 * sp.pi * abs(m[1]) / Lb for m in cre) - sum(2 * sp.pi * abs(m[1]) / Lb for m in ann)


def sector_change(cre, ann):
    dR = sum(m[1] for m in cre if m[1] > 0) - sum(m[1] for m in ann if m[1] > 0)
    dL = sum(-m[1] for m in cre if m[1] < 0) - sum(-m[1] for m in ann if m[1] < 0)
    return dR, dL


def r1_direct_me(sgn, terms, outs, ins):
    """<outs| Int_0^L dx :L: |ins> for pi-free chiral terms, distinct operators, box-normalised states.

    ins/outs: lists of (field, direction, label expression); labels may be symbolic. Distinct labels are assumed
    (the state norms are then 1 and every matching assignment of fields to operators is a distinct term).
    Momentum conservation is the caller's responsibility and is checked there.
    """
    targets = [("cre", f, d, n) for f, d, n in outs] + [("ann", f, d, n) for f, d, n in ins]
    total = sp.Integer(0)
    prod_c = sp.Mul(*[c1(sgn, d, part, n) for part, f, d, n in targets])
    for mon, coef in terms:
        if mon[4] or mon[5]:
            continue
        factors = []
        for fam, e in zip(FAM, mon[:4]):
            factors += [fam] * e
        if len(factors) != len(targets):
            continue
        cnt = sum(1 for perm in itertools.permutations(range(len(targets)))
                  if all(factors[i] == (targets[perm[i]][1], targets[perm[i]][2]) for i in range(len(targets))))
        total += coef * Lb * cnt * prod_c
    return sp.expand(total)


# ================================================================ route 2: explicit finite matrices
class NotCovered(Exception):
    pass


class Basis:
    """Occupation-number basis over all oscillator modes 1 <= |n| <= N of f and g, restricted to
    N_R <= BR and N_L <= BL. Membership is tested by this predicate; a state outside raises NotCovered."""

    def __init__(self, N, BR, BL):
        self.N, self.BR, self.BL = N, BR, BL

    def contains(self, st):
        if any(not (1 <= abs(m[1]) <= self.N) for m, o in st):
            return False
        nr, nl = sector_of(st)
        return nr <= self.BR and nl <= self.BL

    def describe(self):
        return dict(N_max=self.N, N_R_max=self.BR, N_L_max=self.BL)


def ladder(kind, mode, st, basis):
    """Matrix action of a (kind 'a') or a^dag (kind 'ad') of one mode on a basis state: (state, N factor)."""
    occ = dict(st)
    o = occ.get(mode, 0)
    if kind == "a":
        if o == 0:
            return None
        occ[mode] = o - 1
        return canon(occ), o
    occ[mode] = o + 1
    new = canon(occ)
    if not basis.contains(new):
        raise NotCovered("%s on %s leaves the basis %s" % (mode, state_str(st), basis.describe()))
    return new, o + 1


GAUSS = {0: (1, 0), 1: (0, 1), 2: (-1, 0), 3: (0, -1)}


class Route2:
    """Route 2. Interaction operators from the derivative fields, by their action on the basis."""

    def __init__(self, sgn, basis):
        self.sgn = sgn
        self.basis = basis
        self.cache = {}
        self.polys = {}
        for name, Lfun in OPERATORS:
            expr = sp.expand(Lfun(FT, FX, GT, GX))
            self.polys[name] = sp.Poly(expr, FT, FX, GT, GX).terms()
        self.degree = max(sum(m) for terms in self.polys.values() for m, c in terms)

    def _apply_A(self, A, st):
        if not A:
            yield st, 0, 0, 0, 0, 1, 1
            return
        fld, dv = A[0]
        for st1, key, nR, nL, ph, sg, rad in self._apply_A(A[1:], st):
            for mode, o in st1:
                if mode[0] != fld:
                    continue
                n = mode[1]
                st2, f_ = ladder("a", mode, st1, self.basis)
                yield (st2, key + n, nR + (n > 0), nL + (n < 0), ph + 1,
                       sg * self.sgn[(dv, direction(n), "ann")], rad * abs(n) * f_)

    def _apply_C(self, C, st, key, nR, nL, ph, sg, rad, N):
        if not C:
            if key == 0:
                yield st, nR, nL, ph, sg, rad
            return
        fld, dv = C[0]
        if len(C) == 1:
            cand = [key] if 1 <= abs(key) <= N else []        # the last creation factor closes the Kronecker delta
        else:
            cand = [n for n in range(-N, N + 1) if n != 0]
        for n in cand:
            st2, f_ = ladder("ad", (fld, n), st, self.basis)
            yield from self._apply_C(C[1:], st2, key - n, nR + (n > 0), nL + (n < 0), ph + 1,
                                     sg * self.sgn[(dv, direction(n), "cre")], rad * abs(n) * f_, N)

    def column(self, name, st):
        """{tag: {out_state: coefficient}} of H_int = -Int dx :L_name: acting on st; tag = (n_R, n_L, n_pi_f, n_pi_g)
        counted from the factors actually chosen (Z -> zero mode, n > 0 -> R, n < 0 -> L)."""
        ck = (name, st)
        if ck in self.cache:
            return self.cache[ck]
        N = self.basis.N
        terms = self.polys[name]
        acc = {}
        for mi, (exps, cmon) in enumerate(terms):
            factors = []
            for (sym, fd), e in zip(FIELD_FACTORS, exps):
                factors += [fd] * e
            for types in itertools.product("ZCA", repeat=len(factors)):
                if any(ty == "Z" and fd[1] == "x" for fd, ty in zip(factors, types)):
                    continue          # derived (ZC2): d_x phi has no zero-mode part
                zf = sum(1 for fd, ty in zip(factors, types) if ty == "Z" and fd[0] == "f")
                zg = sum(1 for fd, ty in zip(factors, types) if ty == "Z" and fd[0] == "g")
                A = [fd for fd, ty in zip(factors, types) if ty == "A"]
                Cc = [fd for fd, ty in zip(factors, types) if ty == "C"]
                for st1, key, nR, nL, ph, sg, rad in self._apply_A(A, st):
                    for st2, nR2, nL2, ph2, sg2, rad2 in self._apply_C(Cc, st1, key, nR, nL, ph, sg, rad, N):
                        d = acc.setdefault((mi, (nR2, nL2, zf, zg)), {}).setdefault(st2, {})
                        re, im = GAUSS[ph2 % 4]
                        cur = d.get(rad2, (0, 0))
                        d[rad2] = (cur[0] + sg2 * re, cur[1] + sg2 * im)
        out = defaultdict(lambda: defaultdict(lambda: sp.Integer(0)))
        for (mi, tag), d in acc.items():
            cmon = terms[mi][1]
            pref = -cmon * Lb * (sp.sqrt(sp.pi) / Lb)**(tag[0] + tag[1]) * (Pi_f / Lb)**tag[2] * (Pi_g / Lb)**tag[3]
            for st2, rads in d.items():
                val = sum((re + sp.I * im) * sp.sqrt(r) for r, (re, im) in rads.items())
                out[tag][st2] += pref * val
        res = {tag: clean(v_) for tag, v_ in out.items()}
        res = {tag: v_ for tag, v_ in res.items() if v_}
        self.cache[ck] = res
        return res

    def group_apply(self, name, tag, vec):
        out = defaultdict(lambda: sp.Integer(0))
        for st, c in vec.items():
            for st2, v_ in self.column(name, st).get(tag, {}).items():
                out[st2] += c * v_
        return clean(out)

    def H0(self, st):
        """Diagonal of H0 from the number operators a^dag a built from the ladder matrices, plus pi^2/(2L)."""
        e = (Pi_f**2 + Pi_g**2) / (2 * Lb)
        for mode, o in st:
            r = ladder("a", mode, st, self.basis)
            st1, f1 = r
            st2, f2 = ladder("ad", mode, st1, self.basis)
            assert st2 == st and f1 == f2 == o
            e += 2 * sp.pi * abs(mode[1]) / Lb * o
        return e


# ================================================================ main
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

    note("spec %s (%s); revision %s; sympy %s" % (SPEC, SPEC_COMMIT, REVISION, sp.__version__))
    note("registered expressions embedded from %s" % REGISTERED_RESULTS)

    # ============================================================ ZC2 derivative fields
    sgn, zc2_flags, zc2_detail = derive_fields()
    zc2 = all(zc2_flags.values()) and len(sgn) == 8
    checks["ZC2_derivative_fields"] = bool(zc2)
    res["ZC2"] = dict(flags=zc2_flags, coefficients_at_x_t_0=zc2_detail,
                      sign_table={"%s %s %s" % k_: v_ for k_, v_ in sorted(sgn.items())},
                      convention=("coefficient = SGN * i * sqrt(pi |n|) / L; a_n carries exp(+i k_n x), a_n^dag "
                                  "exp(-i k_n x) at t = 0; d_t phi = pi_0/L + R_phi + L_phi, d_x phi = L_phi - R_phi"))
    note("ZC2: %s; sign table %s" % ("PASS" if zc2 else "FAIL", res["ZC2"]["sign_table"]))
    if not zc2:
        note("ZC2 failed: the derivative fields are not of the derived form; stopping")
        return finish(out, T0, log, checks, res, findings, note)
    dt_rep = {"f": P0f + Rf + Lf, "g": P0g + Rg + Lg}
    dx_rep = {"f": Lf - Rf, "g": Lg - Rg}

    # ============================================================ ZC1 mode algebra
    def zc1_fourier(N, fa, fb):
        """Fourier coefficients F_m of [phi_fa(x), d_t phi_fb(y)] at t = 0, keyed by m (dependence exp(i k_m (x-y)))."""
        # generators: ('q', f), ('pi', f), ('a', mode), ('ad', mode); coefficient and Fourier index in x or y
        phi = {("q", fa): (sp.Integer(1), 0)}
        dphi = {("pi", fb): (1 / Lb, 0)}
        for n in list(range(1, N + 1)) + list(range(-N, 0)):
            w = 1 / sp.sqrt(2 * (2 * sp.pi * abs(n) / Lb) * Lb)
            phi[("a", (fa, n))] = (w, n)
            phi[("ad", (fa, n))] = (w, -n)
            dphi[("a", (fb, n))] = (sgn[("t", direction(n), "ann")] * sp.I * sp.sqrt(sp.pi * abs(n)) / Lb, n)
            dphi[("ad", (fb, n))] = (sgn[("t", direction(n), "cre")] * sp.I * sp.sqrt(sp.pi * abs(n)) / Lb, -n)

        def ccr(g1, g2):
            if g1[0] == "a" and g2[0] == "ad" and g1[1] == g2[1]:
                return 1
            if g1[0] == "ad" and g2[0] == "a" and g1[1] == g2[1]:
                return -1
            if g1[0] == "q" and g2[0] == "pi" and g1[1] == g2[1]:
                return sp.I
            if g1[0] == "pi" and g2[0] == "q" and g1[1] == g2[1]:
                return -sp.I
            return 0

        F = defaultdict(lambda: sp.Integer(0))
        translation = True
        for g1, (cx, jx) in phi.items():
            for g2, (cy, jy) in dphi.items():
                c_ = ccr(g1, g2)
                if c_ == 0:
                    continue
                if jx + jy != 0:          # exp(2 pi i (jx x + jy y) / L) must depend on x - y only
                    translation = False
                F[jx] += c_ * cx * cy
        return {m: sp.simplify(v_) for m, v_ in F.items()}, translation

    zc1_rows = {}
    zc1_ok = True
    for fld in FIELDS:
        for N in ZC1_N_LIST:
            F, tr = zc1_fourier(N, fld, fld)
            inside = all(is_zero(F.get(m, 0) - sp.I / Lb) for m in range(-N, N + 1))
            outside = all(is_zero(F.get(m, 0)) for m in list(range(N + 1, N + 3)) + list(range(-N - 2, -N)))
            extra = [m for m in F if abs(m) > N and not is_zero(F[m])]
            ok = tr and inside and outside and not extra
            zc1_ok = zc1_ok and ok
            zc1_rows["%s N_max=%d" % (fld, N)] = dict(fourier={str(m): str(F.get(m, 0)) for m in range(-N - 2, N + 3)},
                                                      translation_invariant=tr, truncated_delta=ok)
        Fc, trc = zc1_fourier(ZC1_N_LIST[-1], fld, "g" if fld == "f" else "f")
        cross_zero = all(is_zero(v_) for v_ in Fc.values())
        zc1_ok = zc1_ok and cross_zero
        zc1_rows["%s cross-field" % fld] = dict(all_zero=bool(cross_zero))
    # symbolic single label: the coefficient of exp(i k_m (x - y)) for |m| <= N_max is independent of N_max
    nu = sp.Symbol("nu", positive=True, integer=True)
    w_nu = 1 / sp.sqrt(2 * (2 * sp.pi * nu / Lb) * Lb)
    # m = +nu: [a_nu, a_nu^dag] from phi(x) a-part and d_t phi(y) a^dag-part; and [a_-nu^dag, a_-nu] (= -1) with
    # phi(x) a^dag_-nu e^{-i k_-nu x} = e^{+i k_nu x} and d_t phi(y) a_-nu e^{i k_-nu y} = e^{-i k_nu y}
    F_sym_pos = sp.simplify(w_nu * sgn[("t", "R", "cre")] * sp.I * sp.sqrt(sp.pi * nu) / Lb * 1
                            + w_nu * sgn[("t", "L", "ann")] * sp.I * sp.sqrt(sp.pi * nu) / Lb * (-1))
    F_sym_neg = sp.simplify(w_nu * sgn[("t", "L", "cre")] * sp.I * sp.sqrt(sp.pi * nu) / Lb * 1
                            + w_nu * sgn[("t", "R", "ann")] * sp.I * sp.sqrt(sp.pi * nu) / Lb * (-1))
    sym_ok = is_zero(F_sym_pos - sp.I / Lb) and is_zero(F_sym_neg - sp.I / Lb)
    zc1_ok = zc1_ok and sym_ok
    checks["ZC1_mode_algebra"] = bool(zc1_ok)
    res["ZC1"] = dict(rows=zc1_rows, symbolic_label=dict(m_positive=str(F_sym_pos), m_negative=str(F_sym_neg),
                                                         zero_mode="i/L from [q0, pi0/L]"),
                      finite_N=("[phi(x), d_t phi(y)] = i (1/L) sum_{|m| <= N_max} exp(i k_m (x-y)): every Fourier "
                                "coefficient with |m| <= N_max equals i/L (m = 0 from the zero mode), all others 0"),
                      distributional_limit=("each Fourier coefficient i/L is independent of N_max once N_max >= |m|; "
                                            "as N_max -> infinity at fixed L every coefficient equals i/L, i.e. "
                                            "i delta_L(x - y) in the sense of distributions"))
    note("ZC1: truncated periodic delta for f and g at N_max in %s, cross-field zero, symbolic coefficient %s -> %s"
         % (ZC1_N_LIST, F_sym_pos, "PASS" if zc1_ok else "FAIL"))

    # ============================================================ Z1 field level and comparison with F1
    groups = {name: chiral_groups(Lfun, dt_rep, dx_rep) for name, Lfun in OPERATORS}
    z1 = {}
    f1_ok = True
    f1_rows = {}
    for name in groups:
        mine = defaultdict(lambda: sp.Integer(0))
        for tag, terms in groups[name].items():
            if tag[2] == 0 and tag[3] == 0:
                mine[tag[:2]] += sp.Add(*[mono_expr(m, c) for m, c in terms])
        for key in set(mine) | set(REG_F1[name]):
            reg = sp.Add(*[sp.sympify(s_, locals=REG_LOCALS) for s_ in REG_F1[name].get(key, [])])
            mp = sp.Poly(mine.get(key, 0), *CH[:4]).terms() if mine.get(key, 0) != 0 else []
            rp = sp.Poly(reg, *CH[:4]).terms() if reg != 0 else []
            md, rd = dict(mp), dict(rp)
            same = set(md) == set(rd) and all(is_zero(md[k_] - rd[k_]) for k_ in md)
            f1_ok = f1_ok and same
            f1_rows["%s %d_R_%d_L" % (name, *key)] = bool(same)
    # pi-dependent groups versus the registered (symbolic) candidate list
    pi_vs_reg = {}
    pi_same = True
    for name, reg_tab in (("L_contact", REG_F2_ZERO_MODE), ("L_lambda", REG_F2_ZERO_MODE_LAMBDA)):
        mine = defaultdict(lambda: sp.Integer(0))
        for tag, terms in groups[name].items():
            if tag[2] + tag[3] > 0:
                mine[(tag[0], tag[1], tag[2] + tag[3])] += sp.Add(*[mono_expr(m, c) for m, c in terms])
        for key in set(mine) | set(reg_tab):
            reg = sp.Add(*[sp.sympify(s_, locals=REG_LOCALS) for s_ in reg_tab.get(key, [])])
            same = is_zero(mine.get(key, 0) - reg)
            pi_same = pi_same and same
            pi_vs_reg["%s %d_R_%d_L_%d_P0" % (name, *key)] = bool(same)
    for name in groups:
        for tag, terms in sorted(groups[name].items()):
            z1[tag_str(name, tag)] = dict(
                field_level=[str(mono_expr(m, c)) for m, c in terms],
                field_level_in_pi=[str(sp.simplify(mono_expr(m, c).subs(P0_TO_PI))) for m, c in terms],
                pi_dependent=bool(tag[2] + tag[3] > 0))
    res["Z1"] = dict(groups=z1, F1_field_level=f1_rows, pi_groups_vs_registered_F2=pi_vs_reg,
                     mode_level_rule=("H_int = -Int_0^L dx :L_int:; each R_phi (L_phi) factor contributes "
                                      "+i sqrt(pi|n|)/L a^dag_{phi,n} or -i sqrt(pi|n|)/L a_{phi,n} with n > 0 (n < 0); "
                                      "each P0_phi factor pi_{phi,0}/L; the integral gives L delta(sum_cre n - "
                                      "sum_ann n, 0)"))
    note("Z1: field-level pi-free part = registered F1 split: %s; pi groups = registered F2 list: %s"
         % (f1_ok, pi_same))
    if not pi_same:
        findings.append(dict(id="Z1-pi-groups", text="box pi-dependent groups differ from the registered F2 list",
                             detail=pi_vs_reg))

    # ============================================================ ZC3 / ZC6 mode level and the box-to-continuum map
    nA = sp.Symbol("n_a", integer=True)
    box_comm = sp.sqrt(2 * sp.Abs(2 * sp.pi * nA / Lb) * Lb) * sp.sqrt(2 * sp.Abs(2 * sp.pi * nA / Lb) * Lb)
    box_comm_target = 2 * sp.Abs(2 * sp.pi * nA / Lb) * Lb
    map_ok = is_zero(box_comm - box_comm_target)
    map_ratio = sp.simplify(box_comm_target / (2 * sp.Abs(2 * sp.pi * nA / Lb)))   # = L, the image of (2 pi) delta
    map_ok = map_ok and is_zero(map_ratio - Lb)
    res["map"] = dict(box="[a(k_n), a^dag(k_m)] = %s delta_nm with a(k_n) = sqrt(2|k_n| L) a_n" % box_comm_target,
                      continuum="[a(k), a^dag(k')] = 2|k| (2 pi) delta(k - k') (registered C9)",
                      correspondence="(2 pi) delta(k - k') <-> %s delta_nm" % map_ratio,
                      amplitude=("T = <out| Int dx :L: |in>_box * prod_j sqrt(2 |k_j| L) / L, the 1/L replacing "
                                 "(2 pi) delta(P_out - P_in)"))

    pifree = {name: [(m, c) for tag, terms in groups[name].items() if tag[2] == 0 and tag[3] == 0 for m, c in terms]
              for name in groups}
    n1 = sp.Symbol("n1", positive=True, integer=True)
    n2 = sp.Symbol("n2", negative=True, integer=True)
    ERs, ELs = REG_LOCALS["E_R"], REG_LOCALS["E_L"]
    mode_rows = {}
    mode_ok = True

    def T_box(name, outs, ins):
        me = r1_direct_me(sgn, pifree[name], outs, ins)
        norm = sp.Mul(*[sp.sqrt(2 * sp.Abs(2 * sp.pi * n_ / Lb) * Lb) for f_, d_, n_ in outs + ins])
        return sp.simplify(me * norm / Lb)

    elastic_T = {}
    for name in groups:
        for i, (c_, d_) in enumerate(REG_CHANNELS):
            for j, (a_, b_) in enumerate(REG_CHANNELS):
                tb = T_box(name, [(c_, "R", n1), (d_, "L", n2)], [(a_, "R", n1), (b_, "L", n2)])
                reg = sp.sympify(REG_T[name][i][j], locals=REG_LOCALS).subs(
                    {ERs: 2 * sp.pi * n1 / Lb, ELs: -2 * sp.pi * n2 / Lb})
                ok = is_zero(tb - reg)
                mode_ok = mode_ok and ok
                elastic_T[(name, i, j)] = tb
                mode_rows["%s elastic %s_R %s_L -> %s_R %s_L" % (name, a_, b_, c_, d_)] = dict(
                    T_box=str(tb), registered_at_lattice=str(reg), equal=bool(ok))
    m1, m2, m3 = sp.symbols("m1 m2 m3", positive=True, integer=True)
    split_T = {}
    for (name, a_, finals), sreg in REG_SPLIT.items():
        kk = {REG_LOCALS["k1"]: 2 * sp.pi * m1 / Lb, REG_LOCALS["k2"]: 2 * sp.pi * m2 / Lb,
              REG_LOCALS["k3"]: 2 * sp.pi * m3 / Lb}
        reg = sp.sympify(sreg, locals=REG_LOCALS).subs(kk)
        for dr, s_ in (("R", 1), ("L", -1)):
            ins = [(a_, dr, s_ * (m1 + m2 + m3))]
            outs = [(finals[0], dr, s_ * m1), (finals[1], dr, s_ * m2), (finals[2], dr, s_ * m3)]
            mom = sum(n_ for f_, d_, n_ in ins) - sum(n_ for f_, d_, n_ in outs)
            tb = T_box(name, outs, ins)
            ok = is_zero(mom) and is_zero(tb - reg)          # L -> L L L compared by parity (registered F1)
            mode_ok = mode_ok and ok
            split_T[(name, a_, finals, dr)] = tb
            mode_rows["%s %s_%s -> %s (distinct labels m1, m2, m3)" % (name, a_, dr, " ".join(f_ + "_" + dr for f_ in finals))] = \
                dict(T_box=str(tb), registered_at_lattice=str(reg), equal=bool(ok))
    # route-1 self-consistency: direct symbolic-label elements = enumerated operator at concrete labels (N_max = 3)
    selfc = True
    Nc = N_LIST[-1]
    for name in groups:
        op = r1_build(sgn, pifree[name], Nc, sign=-1)
        for i, (c_, d_) in enumerate(REG_CHANNELS):
            for j, (a_, b_) in enumerate(REG_CHANNELS):
                st_in = canon({(a_, 1): 1, (b_, -2): 1})
                st_out = canon({(c_, 1): 1, (d_, -2): 1})
                me_op = -r1_apply(op, st_in).get(st_out, 0)      # <out| Int :L: |in> = -<out| H_int |in>
                me_dir = r1_direct_me(sgn, pifree[name], [(c_, "R", 1), (d_, "L", -2)], [(a_, "R", 1), (b_, "L", -2)])
                selfc = selfc and is_zero(me_op - me_dir)
    checks["ZC3_continuum_consistency"] = bool(f1_ok and mode_ok and map_ok and selfc)
    res["ZC3"] = dict(field_level=f1_rows, mode_level=mode_rows, map_derived=bool(map_ok),
                      route1_direct_equals_enumerated=bool(selfc),
                      scope=("mode level uses the pi-free part only; distinct labels assumed for the splitting "
                             "elements; L -> L L L compared with the registered R -> R R R by parity"))
    note("ZC3: field level %s; mode level %s; map %s; route-1 self-consistency %s -> %s"
         % (f1_ok, mode_ok, map_ok, selfc, "PASS" if checks["ZC3_continuum_consistency"] else "FAIL"))

    # ============================================================ Z2 / ZC4 zero-mode momenta
    gens_seen = set()
    op_by = {}
    for N in N_LIST:
        for name in groups:
            for tag, terms in groups[name].items():
                op_by[(name, tag, N)] = r1_build(sgn, terms, N, sign=-1)
                for (cre, ann), c in op_by[(name, tag, N)].items():
                    gens_seen |= {("ad", m) for m in cre} | {("a", m) for m in ann}
                    gens_seen |= {("pi", f_) for f_, s_ in (("f", Pi_f), ("g", Pi_g)) if s_ in c.free_symbols}
    q_symbols_in_H = any(str(s_).startswith("q_") for name in groups for tag in groups[name]
                         for m, c in groups[name][tag] for s_ in c.free_symbols)
    q_gen_in_H = any(g[0] == "q" for g in gens_seen)

    def ccr_gen(g1, g2):
        if g1[0] == "pi" and g2[0] == "q" and g1[1] == g2[1]:
            return -sp.I
        return 0                                         # pi commutes with pi and with every oscillator

    z2 = {}
    zc4 = not q_symbols_in_H and not q_gen_in_H and all("q0_absent" not in k_ or v_ for k_, v_ in zc2_flags.items())
    for fld in FIELDS:
        nonzero = [g for g in sorted(gens_seen | {("pi", "f"), ("pi", "g")}) if ccr_gen(("pi", fld), g) != 0]
        z2["[pi_%s0, H0 + H_int]" % fld] = ("0: H0 + H_int is a polynomial in {pi_f0, pi_g0, a, a^dag} only "
                                            "(q_%s0 absent), and pi_%s0 commutes with each generator" % (fld, fld)
                                            if not nonzero else "NONZERO via generators %s" % nonzero)
        zc4 = zc4 and not nonzero
    checks["ZC4_zero_mode_momenta"] = bool(zc4)
    res["Z2"] = dict(q_absent_in_derivative_fields=all(v_ for k_, v_ in zc2_flags.items() if "q0_absent" in k_),
                     q_absent_in_H=not (q_symbols_in_H or q_gen_in_H), commutators=z2,
                     decomposition=("both commutators vanish (derived): the problem decomposes into joint "
                                    "(pi_f0, pi_g0) eigensectors, in which pi_phi0/L are c-numbers"
                                    if zc4 else "decomposition NOT used"),
                     H0_zero_mode_part="pi_f0^2/(2L) + pi_g0^2/(2L) commutes with H_int for the same reason")
    note("Z2/ZC4: q0 absent %s; [pi_f0, H] and [pi_g0, H] = %s -> %s"
         % (res["Z2"]["q_absent_in_H"], list(z2.values()), "PASS" if zc4 else "FAIL"))

    # ============================================================ Z3 / ZC5 lattice commutators, two routes
    QR, QL = ("L_contact", (4, 0, 0, 0)), ("L_contact", (0, 4, 0, 0))
    comm_names = []
    for name in groups:
        for tag in sorted(groups[name]):
            comm_names.append(("XH0", name, tag))
    comm_names.append(("QRQL", None, None))
    z3 = {}
    blocks = defaultdict(dict)            # blocks[cname][(N, in_state)] = {out: value}
    inter = defaultdict(dict)             # intermediate states for [Q_R, Q_L]
    zc5_ok = True
    not_covered = []
    sectors_reached = defaultdict(set)
    structure = {}
    for N in N_LIST:
        H0r1 = r1_H0(N)
        r1c = {}
        for kind, name, tag in comm_names:
            if kind == "XH0":
                X = op_by[(name, tag, N)]
                cm = r1_commutator(X, H0r1)
                # structure: [M, H0] = -(E_cre - E_ann) M for every normal-ordered momentum-conserving monomial M
                expect = clean({k_: -energy_change(*k_) * v_ for k_, v_ in X.items()})
                struct_ok = vec_diff_zero(cm, expect)
                dsec = sorted({sector_change(*k_) for k_ in cm})
                structure["%s N_max=%d" % (tag_str(name, tag), N)] = dict(
                    exactly_zero=not cm, structure_minus_dE_times_M=bool(struct_ok), n_monomials=len(cm),
                    sector_changes_dNR_dNL=[list(map(int, d_)) for d_ in dsec])
                zc5_ok = zc5_ok and struct_ok
                r1c[(kind, name, tag)] = cm
            else:
                cm = r1_commutator(op_by[(QR[0], QR[1], N)], op_by[(QL[0], QL[1], N)])
                structure["[Q_R^osc, Q_L^osc] N_max=%d" % N] = dict(exactly_zero=not cm, n_monomials=len(cm))
                r1c[(kind, name, tag)] = cm
        reach = 2 * 4 * N                       # degree 4 of L_int, at most N per created label, two factors
        basis = Basis(N, max(s_[0] for s_ in SECTORS) + reach, max(s_[1] for s_ in SECTORS) + reach)
        r2 = Route2(sgn, basis)
        assert r2.degree == 4
        for S in SECTORS:
            for st in states_in_sector(S[0], S[1], N):
                try:
                    cols = {name: r2.column(name, st) for name in groups}
                    for kind, name, tag in comm_names:
                        cname = "[%s, H0]" % tag_str(name, tag) if kind == "XH0" else "[Q_R^osc, Q_L^osc]"
                        v1 = r1_apply(r1c[(kind, name, tag)], st)
                        if kind == "XH0":
                            gx_ = cols[name].get(tag, {})
                            e_in = r2.H0(st)
                            v2 = clean({o: gx_[o] * (e_in - r2.H0(o)) for o in gx_})
                            for o in gx_:
                                sectors_reached[cname].add(sector_of(o))
                            mids = []
                        else:
                            qlv = r2.group_apply(QL[0], QL[1], {st: 1})
                            qrv = r2.group_apply(QR[0], QR[1], {st: 1})
                            a_ = r2.group_apply(QR[0], QR[1], qlv)
                            b_ = r2.group_apply(QL[0], QL[1], qrv)
                            v2 = clean({o: a_.get(o, 0) - b_.get(o, 0) for o in set(a_) | set(b_)})
                            mids = list(qlv) + list(qrv)
                            for o in mids:
                                sectors_reached[cname].add(sector_of(o))
                        agree = vec_diff_zero(v1, v2)
                        zc5_ok = zc5_ok and agree
                        if not agree:
                            findings.append(dict(id="ZC5-mismatch", text="%s N_max=%d %s" % (cname, N, state_str(st))))
                        blocks[cname][(N, st)] = v1
                        inter[cname][(N, st)] = mids
                except NotCovered as e_:
                    zc5_ok = False
                    not_covered.append(dict(N_max=N, state=state_str(st), reason=str(e_)))
        # group actions themselves (route 1 = route 2), as an additional cross-check of the construction
        for S in SECTORS:
            for st in states_in_sector(S[0], S[1], N):
                try:
                    for name in groups:
                        col = r2.column(name, st)
                        for tag in set(groups[name]) | set(col):
                            v1 = r1_apply(op_by[(name, tag, N)], st) if tag in groups[name] else {}
                            v2 = col.get(tag, {})
                            ok = vec_diff_zero(v1, v2)
                            zc5_ok = zc5_ok and ok
                            if not ok:
                                findings.append(dict(id="ZC5-group-mismatch",
                                                     text="%s N_max=%d %s" % (tag_str(name, tag), N, state_str(st))))
                except NotCovered as e_:
                    zc5_ok = False
                    not_covered.append(dict(N_max=N, state=state_str(st), reason=str(e_)))
        note("Z3: N_max=%d evaluated on sectors %s (basis %s); ZC5 so far %s"
             % (N, SECTORS, basis.describe(), zc5_ok))

    # cutoff independence against the complete set of external and intermediate sectors
    cutoff = {}
    for cname, bl in blocks.items():
        indep = True
        dependent_image = False
        compared = 0
        for N, N2 in itertools.combinations(N_LIST, 2):
            for (Nb, st), v_ in bl.items():
                if Nb != N or level(st) > N:
                    continue
                v2 = bl.get((N2, st))
                if v2 is None:
                    continue
                mids = inter[cname].get((N2, st), [])
                mid_level = max([level(o) for o in mids], default=0)
                for o in set(v_) | set(v2):
                    req = max(level(st), level(o), mid_level)
                    if req <= N:
                        compared += 1
                        if not is_zero(v_.get(o, 0) - v2.get(o, 0)):
                            indep = False
                    elif not is_zero(v2.get(o, 0)):
                        dependent_image = True
        cutoff[cname] = dict(blocks_independent_of_N_max=bool(indep), blocks_compared=compared,
                             image_of_sector_depends_on_N_max=bool(dependent_image),
                             rule=("a block (in, out) is compared between N_max and N_max' > N_max only if N_max >= "
                                   "N_req = max(level(in), level(out), largest level of the intermediate sectors "
                                   "reached at N_max'); level = max(N_R, N_L)"))
    z3 = dict(structure=structure, cutoff_independence=cutoff,
              sectors_reached={k_: sorted(map(list, v_)) for k_, v_ in sectors_reached.items()},
              not_covered=not_covered,
              Q_definitions=("Q_R^osc = H_int restricted to the (4,0,0,0) group of L_contact (n > 0 oscillators only), "
                             "Q_L^osc likewise (0,4,0,0), as in the registered F2"),
              pi_sq_part="[pi^2/(2L), X] = 0 for every X here (Z2); route 2 keeps it in H0 and it cancels")
    # registered symbolic F2 table versus the box
    f2_cmp = {}
    for name in groups:
        agg = defaultdict(lambda: True)
        for tag in groups[name]:
            for N in N_LIST:
                agg[(tag[0], tag[1], tag[2] + tag[3])] &= structure["%s N_max=%d" % (tag_str(name, tag), N)]["exactly_zero"]
        for key, box in agg.items():
            k3 = "%s %d_R_%d_L_%d_P0" % (name, *key)
            regv = REG_F2_COMMUTES.get(k3)
            f2_cmp[k3] = dict(box_commutes_with_H0=bool(box), registered=regv,
                              status="confirmed" if regv == box else ("not in registered table" if regv is None
                                                                       else "contradicted"))
    for k3, regv in REG_F2_COMMUTES.items():
        if k3 not in f2_cmp:
            f2_cmp[k3] = dict(box_commutes_with_H0=None, registered=regv, status="no box group (refined)")
    z3["registered_F2_table_vs_box"] = f2_cmp
    res["Z3"] = z3
    zc5_ok = zc5_ok and not not_covered
    checks["ZC5_lattice_commutators_two_routes"] = bool(zc5_ok)
    note("Z3: [Q_R,H0], [Q_L,H0], [Q_R,Q_L] exactly zero at all N_max: %s; registered F2 table: %s"
         % ([structure["%s N_max=%d" % (tag_str(*QR), N)]["exactly_zero"] for N in N_LIST]
            + [structure["%s N_max=%d" % (tag_str(*QL), N)]["exactly_zero"] for N in N_LIST]
            + [structure["[Q_R^osc, Q_L^osc] N_max=%d" % N]["exactly_zero"] for N in N_LIST],
            {k_: v_["status"] for k_, v_ in f2_cmp.items()}))
    note("ZC5: route 1 = route 2 on every evaluated sector, NOT_COVERED entries %d -> %s"
         % (len(not_covered), "PASS" if zc5_ok else "FAIL"))

    # ============================================================ Z4 classification of the pi-dependent groups
    NR_, NL_, NRp, NLp = sp.symbols("N_R N_L N_Rp N_Lp", nonnegative=True, integer=True)

    def E0(pf, pg, nr, nl):
        return (pf**2 + pg**2) / (2 * Lb) + 2 * sp.pi / Lb * (nr + nl)

    cancel = sp.simplify(E0(Pi_f, Pi_g, NRp, NLp) - E0(Pi_f, Pi_g, NR_, NL_))
    cancel_ok = is_zero(cancel - 2 * sp.pi / Lb * (NRp + NLp - NR_ - NL_))
    z4 = {}
    Nz = N_LIST[-1]
    for name in groups:
        for tag in sorted(groups[name]):
            if tag[2] + tag[3] == 0:
                continue
            op = op_by[(name, tag, Nz)]
            rows = sorted({(len(c_) - len(a_), *sector_change(c_, a_)) for (c_, a_) in op})
            mom_ok = all(r_[1] == r_[2] for r_ in rows)
            equal_E = [(c_, a_) for (c_, a_) in op if sum(sector_change(c_, a_)) == 0]
            example = None
            if equal_E:
                c_, a_ = equal_E[0]
                st_in = canon(Counter(a_))
                v_ = r1_apply({(c_, a_): op[(c_, a_)]}, st_in)
                example = dict(monomial=mono_str(c_, a_), in_state=state_str(st_in),
                               out={state_str(k_): str(w_) for k_, w_ in v_.items()})
            z4[tag_str(name, tag)] = dict(
                particle_number_changes=sorted({int(r_[0]) for r_ in rows}),
                changes_particle_number=any(r_[0] != 0 for r_ in rows),
                sector_changes_dNR_dNL=sorted({(int(r_[1]), int(r_[2])) for r_ in rows}),
                conserves_N_R_and_N_L_separately=all(r_[1] == 0 and r_[2] == 0 for r_ in rows),
                momentum_dNR_equals_dNL=bool(mom_ok),
                nonzero_equal_E0_elements=bool(equal_E),
                equal_E0_example=example, N_max=Nz)
    res["Z4"] = dict(E0_difference_same_pi_sector=str(cancel), cancellation_shown=bool(cancel_ok), groups=z4,
                     scope="classification only; no probabilities or secular coefficients (spec section 7)")
    note("Z4: E0(out) - E0(in) in the same pi sector = %s; groups with equal-E0 elements: %s"
         % (cancel, [k_ for k_, v_ in z4.items() if v_["nonzero_equal_E0_elements"]]))
    if any(v_["nonzero_equal_E0_elements"] and v_["changes_particle_number"] for v_ in z4.values()):
        findings.append(dict(id="Z4-equal-energy-number-changing",
                             text=("some pi-dependent groups have nonzero first-order elements between free states "
                                   "of equal E0 that change the particle number; classification only"),
                             groups=[k_ for k_, v_ in z4.items()
                                     if v_["nonzero_equal_E0_elements"] and v_["changes_particle_number"]]))

    # ============================================================ Z5 / ZC6 zero-mode vacuum limit
    kphys = sp.symbols("q1:5", positive=True)
    z5 = {}
    scal_ok = True
    for name in groups:
        for tag in sorted(groups[name]):
            nosc = tag[0] + tag[1]
            expr = Lb * sp.Mul(*[sp.sqrt(sp.pi * (kphys[j] * Lb / (2 * sp.pi))) / Lb for j in range(nosc)]) \
                * (Pi_f / Lb)**tag[2] * (Pi_g / Lb)**tag[3]
            ratio = sp.simplify(expr.subs(Lb, 2 * Lb) / expr)
            alpha = sp.nsimplify(sp.simplify(sp.expand_log(sp.log(ratio), force=True) / sp.log(2)))
            k3 = tag3_str(name, tag)
            reg_alpha = None
            if k3 in REG_F2_SCALING:
                b_, e_ = sp.sympify(REG_F2_SCALING[k3], locals=REG_LOCALS).as_base_exp()
                reg_alpha = e_ if b_ == Lb else None
                scal_ok = scal_ok and reg_alpha is not None and is_zero(alpha - reg_alpha)
            rel = sp.nsimplify(alpha - (-1))
            z5[tag_str(name, tag)] = dict(
                L_power=str(alpha), registered=str(reg_alpha), relative_to_elastic_and_splitting=str(rel),
                vanishes_relative=bool(rel < 0), pi_dependent=bool(tag[2] + tag[3] > 0))
    # primary commensurate sequences: labels n_j = k_j L / (2 pi) exactly; secondary floor sequences: error only
    kR_, kL_ = sp.symbols("kR kL", positive=True)
    th1, th2, th3 = sp.symbols("theta1:4", nonnegative=True)
    prim_ok = True
    sec_rows = {}
    for (name, i, j), tb in elastic_T.items():
        reg = sp.sympify(REG_T[name][i][j], locals=REG_LOCALS).subs({ERs: kR_, ELs: kL_})
        prim = sp.simplify(tb.subs({n1: kR_ * Lb / (2 * sp.pi), n2: -kL_ * Lb / (2 * sp.pi)}) - reg)
        prim_ok = prim_ok and is_zero(prim)
        if reg != 0:
            sec = sp.expand(tb.subs({n1: (kR_ - 2 * sp.pi * th1 / Lb) * Lb / (2 * sp.pi),
                                     n2: -(kL_ - 2 * sp.pi * th2 / Lb) * Lb / (2 * sp.pi)}) - reg)
            sec_rows["%s elastic %d<-%d" % (name, i, j)] = dict(
                limit=str(sp.limit(sec, Lb, sp.oo)), L_times_error_limit=str(sp.limit(Lb * sec, Lb, sp.oo)))
    for (name, a_, finals, dr), tb in split_T.items():
        kk3 = sp.symbols("kk1:4", positive=True)
        reg = sp.sympify(REG_SPLIT[(name, a_, finals)], locals=REG_LOCALS).subs(
            {REG_LOCALS["k1"]: kk3[0], REG_LOCALS["k2"]: kk3[1], REG_LOCALS["k3"]: kk3[2]})
        prim = sp.simplify(tb.subs({m1: kk3[0] * Lb / (2 * sp.pi), m2: kk3[1] * Lb / (2 * sp.pi),
                                    m3: kk3[2] * Lb / (2 * sp.pi)}) - reg)
        prim_ok = prim_ok and is_zero(prim)
        if reg != 0:
            sub = {m1: (kk3[0] - 2 * sp.pi * th1 / Lb) * Lb / (2 * sp.pi),
                   m2: (kk3[1] - 2 * sp.pi * th2 / Lb) * Lb / (2 * sp.pi),
                   m3: (kk3[2] - 2 * sp.pi * th3 / Lb) * Lb / (2 * sp.pi)}
            sec = sp.expand(tb.subs(sub) - reg)
            sec_rows["%s %s_%s -> %s" % (name, a_, dr, "".join(finals))] = dict(
                limit=str(sp.limit(sec, Lb, sp.oo)), L_times_error_limit=str(sp.limit(Lb * sec, Lb, sp.oo)))
    sec_ok = all(r_["limit"] == "0" for r_ in sec_rows.values())
    checks["ZC6_box_to_continuum_limit"] = bool(
        map_ok and prim_ok and mode_ok and scal_ok
    )
    res["Z5"] = dict(groups=z5, scaling_equals_registered_table=bool(scal_ok),
                     scaling_assumptions=("box-normalised Fock states with O(1) occupations; fixed physical momenta "
                                          "(labels n = k L / (2 pi)); fixed pi_f0, pi_g0 (zero-mode vacuum limit)"),
                     primary_sequences=("commensurate L_j with k = 2 pi n_j / L_j exactly: box elastic block and "
                                        "R -> R R R / L -> L L L elements equal the registered continuum results at "
                                        "every L_j: %s" % prim_ok),
                     secondary_floor_sequences=sec_rows, secondary_error_vanishes=bool(sec_ok),
                     not_analysed="limit at fixed pi/L = w (background; setting of Q2)")
    note("Z5: L powers equal registered table %s; primary sequences exact %s; floor-sequence error -> 0 %s"
         % (scal_ok, prim_ok, sec_ok))
    note("ZC6: %s" % ("PASS" if checks["ZC6_box_to_continuum_limit"] else "FAIL"))

    # ============================================================ Z6 / ZC7 Wick remainders
    nn = sp.Symbol("n", positive=True, integer=True)
    contr = sp.simplify(c1(sgn, "R", "ann", nn) * c1(sgn, "R", "cre", nn))      # <0| X X |0> per mode, X = R
    contrL = sp.simplify(c1(sgn, "L", "ann", -nn) * c1(sgn, "L", "cre", -nn))
    C_closed = sp.factor(sp.summation(contr, (nn, 1, Nsym)))
    C_def = sp.factor(sp.summation((2 * sp.pi * nn / Lb) / (2 * Lb), (nn, 1, Nsym)))
    C_ok = is_zero(C_closed - C_def) and is_zero(contr - contrL)

    def wick_w(e, j):
        return sp.Integer(factorial(e)) / (factorial(j) * 2**j * factorial(e - 2 * j))

    z6 = {}
    zc7_ok = C_ok
    for name in groups:
        for tag, terms in sorted(groups[name].items()):
            entries = []
            for mon, coef in terms:
                for js in itertools.product(*[range(e // 2 + 1) for e in mon[:4]]):
                    if sum(js) == 0:
                        continue
                    w = sp.Mul(*[wick_w(e, j) for e, j in zip(mon[:4], js)])
                    red = tuple(e - 2 * j for e, j in zip(mon[:4], js)) + mon[4:]
                    cf = coef * w * C_closed**sum(js)
                    deg = sum(red[:4])
                    if deg == 0:
                        kind_ = "c-number: Int dx = L"
                    elif deg == 1:
                        kind_ = "single field: Int dx :X: = 0 (no zero Fourier component)"
                    elif deg == 2 and sum(1 for e in red[:4] if e) == 1:
                        kind_ = "same-family bilinear: Int dx :X^2: = sum_n |k_n| a^dag a over that field and direction"
                    else:
                        kind_ = "other (operator computed at finite N_max)"
                    entries.append(dict(term=str(mono_expr(mon, coef)), reduced=str(mono_expr(red, 1)),
                                        coefficient_closed_form=str(sp.factor(cf.subs(P0_TO_PI))),
                                        contractions=int(sum(js)), kind=kind_, _red=red, _cf=cf))
            # explicit verification (route 1): unordered Int dx prod = Int dx :prod: + predicted remainders
            vrows = {}
            for N in Z6_VERIFY_N:
                ok_all = True
                for mon, coef in terms:
                    prod_ = {0: {((), ()): coef * (Pi_f / Lb)**mon[4] * (Pi_g / Lb)**mon[5]}}
                    for fam, e in zip(FAM, mon[:4]):
                        fld, dr = fam
                        s_ = 1 if dr == "R" else -1
                        fd = defaultdict(dict)
                        for n_ in range(1, N + 1):
                            lab = s_ * n_
                            fd[lab][((), ((fld, lab),))] = c1(sgn, dr, "ann", lab)
                            fd[-lab][(((fld, lab),), ())] = c1(sgn, dr, "cre", lab)
                        for _ in range(e):
                            new = defaultdict(lambda: defaultdict(lambda: sp.Integer(0)))
                            for k1_, A_ in prod_.items():
                                for k2_, B_ in fd.items():
                                    for kk_, vv in no_product(A_, B_).items():
                                        new[k1_ + k2_][kk_] += vv
                            prod_ = {k_: dict(v_) for k_, v_ in new.items()}
                    unordered = clean({k_: Lb * v_ for k_, v_ in prod_.get(0, {}).items()})
                    normal = r1_build(sgn, [(mon, coef)], N, sign=1)
                    pred = defaultdict(lambda: sp.Integer(0))
                    for k_, v_ in normal.items():
                        pred[k_] += v_
                    for js in itertools.product(*[range(e // 2 + 1) for e in mon[:4]]):
                        if sum(js) == 0:
                            continue
                        w = sp.Mul(*[wick_w(e, j) for e, j in zip(mon[:4], js)])
                        red = tuple(e - 2 * j for e, j in zip(mon[:4], js)) + mon[4:]
                        cf = coef * w * C_closed.subs(Nsym, N)**sum(js)
                        for k_, v_ in r1_build(sgn, [(red, cf)], N, sign=1).items():
                            pred[k_] += v_
                    ok = vec_diff_zero(unordered, clean(pred))
                    ok_all = ok_all and ok
                vrows["N_max=%d" % N] = bool(ok_all)
                zc7_ok = zc7_ok and ok_all
            # commutators of the integrated remainder operators with H0 (route 1), finite N_max
            comm_rows = {}
            for N in N_LIST:
                zero_all = True
                for en in entries:
                    rop = r1_build(sgn, [(en["_red"], en["_cf"].subs(Nsym, N))], N, sign=-1)
                    cm = r1_commutator(rop, r1_H0(N))
                    zero_all = zero_all and not cm
                comm_rows["N_max=%d" % N] = "all zero" if zero_all else "NONZERO for some remainder"
            z6[tag_str(name, tag)] = dict(
                remainders=[{k_: v_ for k_, v_ in en.items() if not k_.startswith("_")} for en in entries],
                explicit_verification=vrows, commutators_with_H0=comm_rows,
                sign_in_H_int="H_int remainder = -(Int dx of the listed remainder)")
    checks["ZC7_wick_remainders"] = bool(zc7_ok)
    res["Z6"] = dict(C_closed_form=str(C_closed), C_definition="sum_{n=1}^{N_max} k_n / (2 L)", C_equal=bool(C_ok),
                     contraction_per_mode_R=str(contr), contraction_per_mode_L=str(contrL),
                     cross_family_contractions=("zero: R and L of one field use disjoint modes (n > 0, n < 0); f and g "
                                                "modes commute"),
                     k_max="2 pi N_max / L (reported alongside; no relation to a physical cutoff is chosen)",
                     groups=z6, no_limit_taken="no N_max -> infinity limit for divergent quantities; no renormalisation")
    note("Z6: C(N_max, L) = %s (= sum k_n/(2L): %s); explicit verification and closed forms -> %s"
         % (C_closed, C_ok, "PASS" if zc7_ok else "FAIL"))

    return finish(out, T0, log, checks, res, findings, note)


MANDATORY = ["ZC1_mode_algebra", "ZC2_derivative_fields", "ZC3_continuum_consistency", "ZC4_zero_mode_momenta",
             "ZC5_lattice_commutators_two_routes", "ZC6_box_to_continuum_limit", "ZC7_wick_remainders"]


def finish(out, T0, log, checks, res, findings, note):
    verdict = "ALL_CHECKS_PASS" if all(checks.get(k_, False) for k_ in MANDATORY) else "CHECKS_FAILED"
    res["denominators_checked"] = sorted(set(DENOMINATORS))
    summary = dict(spec=SPEC, spec_commit=SPEC_COMMIT, registered=REGISTERED_RESULTS, revision=REVISION,
                   sympy=sp.__version__, parameters=dict(N_LIST=list(N_LIST), SECTORS=[list(s_) for s_ in SECTORS],
                                                         ZC1_N_LIST=list(ZC1_N_LIST), Z6_VERIFY_N=list(Z6_VERIFY_N)),
                   seconds=round(time.time() - T0, 1),
                   checks={k_: ("PASS" if checks.get(k_) else ("FAIL" if k_ in checks else "NOT_COVERED"))
                           for k_ in MANDATORY},
                   mandatory=MANDATORY, results=res, findings=findings, verdict=verdict,
                   scope=("ALL_CHECKS_PASS covers ZC1-ZC7 only; it is not an independent proof and does not select "
                          "O1, O2 or O3. Route 2 covers only the evaluated sectors and N_max values; statements "
                          "beyond them rest on route 1. The candidate status of the zero-mode groups changes only in "
                          "a later reviewed report."))
    (out / "b412_q3_zero_modes_results.json").write_text(json.dumps(summary, indent=1, default=str) + "\n",
                                                         encoding="utf-8", newline="\n")
    for k_ in MANDATORY:
        note("%s %s" % (k_, summary["checks"][k_]))
    note("findings: %s" % [f_["id"] for f_ in findings])
    note(verdict)
    (out / "b412_q3_zero_modes_output.txt").write_text("\n".join(log) + "\n", encoding="utf-8", newline="\n")
    return 0 if verdict == "ALL_CHECKS_PASS" else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(64)
    sys.exit(main(sys.argv[1]))
