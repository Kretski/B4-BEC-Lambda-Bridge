"""T7a step 4 (DRAFT v0.8.1, unregistered): symbolic classification of Delta L^(M^-4) = (d_t J)^2/(2 M^4)
in the off-shell IBP basis.

Specification: T7A_STEP4_CLASSIFICATION_ALGORITHM_PROPOSAL.md v0.8.1 (commit e7ab8d4), canonical rule
T7A_STEP4_CANONICAL_BASIS_PROPOSAL.md v0.7.1, option A rules T7_PROCA_PROTOCOL_DRAFT.md section 2 (v0.6).
Status: DRAFT - unregistered, not a T7a result. No numerical verification, no M^-6, no B4.12 conclusions.
The independent IBP-matrix implementation (mandatory before T7a registration) is NOT part of this script.

Usage: python -I t7a_step4_classify.py <output_dir>
"""
import hashlib
import itertools
import json
import platform
import sys
from fractions import Fraction

import sympy as sp
from sympy.calculus.euler import euler_equations
from sympy.polys.matrices import DomainMatrix

ALGO_VERSION = "v0.8.1"
FIELDS = ("f", "g")
NT, NX = 4, 2


# ---------------------------------------------------------------- section 1: jet symbols
def jsym(phi, a, b):
    return sp.Symbol(f"{phi}_{a}_{b}")


def parse(s):
    phi, a, b = s.name.split("_")
    return phi, int(a), int(b)


def is_jet(s):
    p = s.name.split("_")
    return len(p) == 3 and p[0] in FIELDS


def jets(P):
    return sorted((s for s in P.free_symbols if is_jet(s)), key=lambda s: parse(s))


def Dt(P):
    return sp.expand(sum((sp.diff(P, s) * jsym(*(lambda q: (q[0], q[1] + 1, q[2]))(parse(s))) for s in jets(P)), sp.S(0)))


def Dx(P):
    return sp.expand(sum((sp.diff(P, s) * jsym(*(lambda q: (q[0], q[1], q[2] + 1))(parse(s))) for s in jets(P)), sp.S(0)))


def euler(P, phi):
    out = sp.S(0)
    for s in jets(P):
        q = parse(s)
        if q[0] != phi:
            continue
        term = sp.diff(P, s)
        for _ in range(q[1]):
            term = -Dt(term)
        for _ in range(q[2]):
            term = -Dx(term)
        out += term
    return sp.expand(out)


def label(triples):
    return " ".join(f"{p}_{'t' * a}{'x' * b}" for p, a, b in triples)


# ---------------------------------------------------------------- section 2: V_S
CONTENT = {"f4": ("f", "f", "f", "f"), "f2g2": ("f", "f", "g", "g"), "g4": ("g", "g", "g", "g")}
SECTORS = ("f4", "f2g2", "g4")
FACTORS = sorted((p, a, b) for p in FIELDS for a in range(NT + 1) for b in range(NX + 1) if a + b >= 1)


def gen_VS(sector):
    want = tuple(sorted(CONTENT[sector]))
    out = []
    for combo in itertools.combinations_with_replacement(FACTORS, 4):  # non-decreasing => each multiset once
        if sum(t[1] for t in combo) == NT and sum(t[2] for t in combo) == NX and tuple(sorted(t[0] for t in combo)) == want:
            out.append(tuple(combo))
    return out


def mono_poly(m):
    return sp.Mul(*[jsym(*t) for t in m])


def key(m):  # v0.7.1: k1, k2, k3
    orders = [a + b for _, a, b in m]
    return (max(orders), sum(o * o for o in orders), tuple(sorted(m)))


def decompose(P):
    """Expanded polynomial -> {sorted triple tuple: Fraction}."""
    P = sp.expand(P)
    out = {}
    if P == 0:
        return out
    for term in sp.Add.make_args(P):
        coeff, rest = term.as_coeff_Mul()
        trip = []
        for f_ in sp.Mul.make_args(rest):
            base, e = f_.as_base_exp()
            trip += [parse(base)] * int(e)
        k = tuple(sorted(trip))
        out[k] = out.get(k, Fraction(0)) + Fraction(int(sp.Rational(coeff).p), int(sp.Rational(coeff).q))
    return {k: v for k, v in out.items() if v != 0}


# ---------------------------------------------------------------- section 3: Euler matrix
def euler_row(P):
    row = {}
    for comp in FIELDS:
        for k, v in decompose(euler(P, comp)).items():
            row[(comp, k)] = v
    return row


def frac_to_str(q):
    return f"{q.numerator}/{q.denominator}" if q.denominator != 1 else f"{q.numerator}"


class Echelon:
    """Incremental exact row reduction over Q (own implementation, used for the greedy choice)."""

    def __init__(self, cols):
        self.cols = cols
        self.idx = {c: i for i, c in enumerate(cols)}
        self.piv = {}  # pivot col -> reduced row (dict col->Fraction)

    def reduce(self, row):
        r = {self.idx[c]: v for c, v in row.items()}
        for pc in sorted(self.piv):
            if pc in r and r[pc] != 0:
                f = r[pc]
                for c, v in self.piv[pc].items():
                    r[c] = r.get(c, Fraction(0)) - f * v
                r = {c: v for c, v in r.items() if v != 0}
        return r

    def add(self, row):
        r = self.reduce(row)
        if not r:
            return False
        pc = min(r)
        inv = 1 / r[pc]
        r = {c: v * inv for c, v in r.items()}
        for k in list(self.piv):  # keep fully reduced
            if pc in self.piv[k]:
                f = self.piv[k][pc]
                nr = dict(self.piv[k])
                for c, v in r.items():
                    nr[c] = nr.get(c, Fraction(0)) - f * v
                self.piv[k] = {c: v for c, v in nr.items() if v != 0}
        self.piv[pc] = r
        return True


def dm_from_rows(rows, cols):
    idx = {c: i for i, c in enumerate(cols)}
    M = [[sp.Rational(0)] * len(cols) for _ in rows]
    for i, r in enumerate(rows):
        for c, v in r.items():
            M[i][idx[c]] = sp.Rational(v.numerator, v.denominator)
    return DomainMatrix.from_list_sympy(len(rows), len(cols), M).convert_to(sp.QQ) if rows else None


def dm_rank(rows, cols):
    return 0 if not rows else dm_from_rows(rows, cols).rank()


def solve_coords(basis_rows, target_row, cols):
    """Exact solve target = sum c_k basis_k (DomainMatrix rref over QQ). Raises if inconsistent."""
    n = len(basis_rows)
    idx = {c: i for i, c in enumerate(cols)}
    for c in target_row:
        if c not in idx:
            raise RuntimeError("target has Euler column outside the fixed column set")
    A = [[sp.Rational(0)] * (n + 1) for _ in cols]
    for k, r in enumerate(basis_rows):
        for c, v in r.items():
            A[idx[c]][k] = sp.Rational(v.numerator, v.denominator)
    for c, v in target_row.items():
        A[idx[c]][n] = sp.Rational(v.numerator, v.denominator)
    R, pivots = DomainMatrix.from_list_sympy(len(cols), n + 1, A).convert_to(sp.QQ).rref()
    if n in pivots:
        raise RuntimeError("inconsistent system: target not in span of basis Euler images")
    if len(pivots) != n:
        raise RuntimeError("basis rows not independent")
    Rl = R.to_Matrix()
    coords = [Fraction(0)] * n
    for i, pc in enumerate(pivots):
        v = sp.Rational(Rl[i, n])
        coords[pc] = Fraction(int(v.p), int(v.q))
    return coords


# ---------------------------------------------------------------- section 0: sector operators
def X_sector(sector):
    A = Dt(jsym("f", 1, 0) * jsym("f", 0, 1))
    B = Dt(jsym("g", 1, 0) * jsym("g", 0, 1))
    return sp.expand({"f4": A * A, "f2g2": A * B, "g4": B * B}[sector])


SWAP = {jsym(p, a, b): jsym("g" if p == "f" else "f", a, b) for p in FIELDS for a in range(9) for b in range(9)}


def swap(P):
    return sp.expand(P.xreplace(SWAP))


# ---------------------------------------------------------------- C6: sympy euler_equations
t_, x_ = sp.symbols("t x", real=True)
FF = {"f": sp.Function("f")(t_, x_), "g": sp.Function("g")(t_, x_)}


def to_function(P):
    rep = {}
    for s in jets(P):
        p, a, b = parse(s)
        d = FF[p]
        if a:
            d = sp.diff(d, t_, a)
        if b:
            d = sp.diff(d, x_, b)
        rep[s] = d
    return P.xreplace(rep)


def c6_check(P):
    # euler_equations omits the equation of a field the Lagrangian does not contain, so each component is
    # referenced separately; for an absent field the reference is 0 (and the own E_phi must vanish).
    L = to_function(P)
    ok = True
    for comp in FIELDS:
        if L.has(FF[comp]):
            (eq,) = euler_equations(L, [FF[comp]], [t_, x_])
            ref = sp.expand((eq.lhs - eq.rhs).doit())
        else:
            ref = sp.S(0)
        mine = sp.expand(to_function(euler(P, comp)).doit())
        ok = ok and sp.simplify(ref - mine) == 0
    return ok


# ---------------------------------------------------------------- C7: exact Gaussian integrals
T, Xs = sp.symbols("T X", real=True)
TEST_PAIRS = [
    {"P_f": 1 + T + Xs + T * Xs**2, "P_g": 2 - T + T**2 * Xs},
    {"P_f": T - Xs**2 + 3 * T**2 * Xs, "P_g": 1 + Xs + T**3},
]


def jet_value(P0, a, b):
    """d_t^a d_x^b (P0 e^{-(T^2+X^2)}) = Q e^{-(T^2+X^2)}; returns Q."""
    Q = P0
    for _ in range(a):
        Q = sp.expand(sp.diff(Q, T) - 2 * T * Q)
    for _ in range(b):
        Q = sp.expand(sp.diff(Q, Xs) - 2 * Xs * Q)
    return Q


def gauss_moment_1d(m, c):
    if m % 2:
        return sp.S(0)
    return sp.gamma(sp.Rational(m + 1, 2)) / sp.Integer(c) ** sp.Rational(m + 1, 2)


def integrate_4field(P, pair):
    """Exact integral over R^2 of a 4-field jet polynomial on the test fields; all terms carry e^{-4(T^2+X^2)}."""
    vals = {}
    for s in jets(P):
        p, a, b = parse(s)
        vals[s] = jet_value(pair["P_" + p], a, b)
    integrand = sp.Poly(sp.expand(P.xreplace(vals)), T, Xs)
    tot = sp.S(0)
    for (m, n), cf in integrand.terms():
        tot += cf * gauss_moment_1d(m, 4) * gauss_moment_1d(n, 4)
    return sp.nsimplify(sp.simplify(tot))


def dtJ_squared_poly(gphi, gpsi, pair):
    """Return (Q, integral) for (d_t J)^2 = Q e^{-4(T^2+X^2)} with concrete couplings."""
    ft, fx = jet_value(pair["P_f"], 1, 0), jet_value(pair["P_f"], 0, 1)
    gt, gx = jet_value(pair["P_g"], 1, 0), jet_value(pair["P_g"], 0, 1)
    E2 = sp.exp(-(T**2 + Xs**2))
    J = (gphi * ft * fx + gpsi * gt * gx) * E2**2
    dJ = sp.expand(sp.simplify(sp.diff(J, T) / E2**2))  # polynomial: d_t J = dJ * e^{-2r^2}
    Q = sp.Poly(sp.expand(dJ**2), T, Xs)
    tot = sum((cf * gauss_moment_1d(m, 4) * gauss_moment_1d(n, 4) for (m, n), cf in Q.terms()), sp.S(0))
    return dJ, sp.simplify(tot)


# ---------------------------------------------------------------- main
def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    with open(__file__, "rb") as fh:
        script_sha = hashlib.sha256(fh.read()).hexdigest()

    sectors = {}
    ok_all = True
    failures = []

    for S in SECTORS:
        VS = sorted(gen_VS(S), key=key)
        X = X_sector(S)
        Xd = decompose(X)
        VSset = set(VS)
        if not all(k in VSset for k in Xd):  # section 2.4
            raise RuntimeError(f"X_{S} not contained in V_S")
        rows = {m: euler_row(mono_poly(m)) for m in VS}
        xrow = euler_row(X)
        cols = sorted({c for r in rows.values() for c in r} | set(xrow))
        rank_all = dm_rank([rows[m] for m in VS], cols)
        ech = Echelon(cols)
        B = [m for m in VS if ech.add(rows[m])]
        rank_B = dm_rank([rows[m] for m in B], cols)
        rank_ok = (len(B) == rank_all) and (rank_B == len(B))
        c = solve_coords([rows[m] for m in B], xrow, cols)
        R = sp.expand(sum((sp.Rational(ci.numerator, ci.denominator) * mono_poly(m) for ci, m in zip(c, B)), sp.S(0)))
        resid_ok = euler(sp.expand(X - R), "f") == 0 and euler(sp.expand(X - R), "g") == 0
        sectors[S] = dict(VS=VS, B=B, c=c, R=R, X=X, cols=cols, rows=rows, rank_ok=rank_ok, resid_ok=resid_ok,
                          n_VS=len(VS), dim=rank_all)
        if not (rank_ok and resid_ok):
            ok_all = False
            failures.append(f"sector {S}: rank_ok={rank_ok} residual_zero={resid_ok}")

    def coords_of(S, P):
        d = sectors[S]
        return solve_coords([d["rows"][m] for m in d["B"]], euler_row(P), d["cols"])

    controls = {}

    def ctrl(name, passed, desc, extra=None):
        nonlocal ok_all
        controls[name] = {"result": "PASS" if passed else "FAIL", "description": desc}
        if extra:
            controls[name].update(extra)
        if not passed:
            ok_all = False
            failures.append(name)

    f = lambda a, b: jsym("f", a, b)
    g = lambda a, b: jsym("g", a, b)
    # C1-C3 total derivatives
    Y1 = Dt(f(1, 0) * f(0, 1) * f(1, 0) * f(1, 1))
    Y2 = Dx(Dx(f(1, 0) ** 4))
    Y3 = Dt(f(1, 0) * f(0, 1) * g(1, 0) * g(1, 1))
    ctrl("C1", all(v == 0 for v in coords_of("f4", Y1)) and set(decompose(Y1)) <= set(sectors["f4"]["VS"]),
         "D_t(f_t f_x f_t f_tx) in sector f4 -> c = 0")
    ctrl("C2", sp.expand(Y2 - (12 * f(1, 0) ** 2 * f(1, 1) ** 2 + 4 * f(1, 0) ** 3 * f(1, 2))) == 0
         and all(v == 0 for v in coords_of("f4", Y2)),
         "D_x^2(f_t^4) = 12 f_t^2 f_tx^2 + 4 f_t^3 f_txx in sector f4 -> c = 0")
    ctrl("C3", all(v == 0 for v in coords_of("f2g2", Y3)) and set(decompose(Y3)) <= set(sectors["f2g2"]["VS"]),
         "D_t(f_t f_x g_t g_tx) in sector f2g2 -> c = 0")
    # C4 unit vectors
    c4 = True
    for S in SECTORS:
        B = sectors[S]["B"]
        for k, m in enumerate(B):
            ck = coords_of(S, mono_poly(m))
            c4 = c4 and all((v == 1) if i == k else (v == 0) for i, v in enumerate(ck))
    ctrl("C4", c4, "every basis monomial b_k gives the unit vector e_k (all sectors)")
    # C5 f<->g as IBP classes
    a_ = sectors["f4"]["dim"] == sectors["g4"]["dim"]
    D1 = sp.expand(swap(sectors["g4"]["R"]) - sectors["f4"]["X"])
    b_ = euler(D1, "f") == 0 and euler(D1, "g") == 0
    Xm = sectors["f2g2"]["X"]
    inv = sp.expand(swap(Xm) - Xm) == 0
    D2 = sp.expand(swap(sectors["f2g2"]["R"]) - Xm)
    c_ = inv and euler(D2, "f") == 0 and euler(D2, "g") == 0
    recanon = coords_of("f2g2", swap(sectors["f2g2"]["R"])) == sectors["f2g2"]["c"]
    ctrl("C5", a_ and b_ and c_ and recanon,
         "f<->g compared as IBP classes: (a) dim f4 = dim g4; (b) E(sigma(R_g4) - X_f4) = 0; "
         "(c) X_f2g2 sigma-invariant, E(sigma(R_f2g2) - X_f2g2) = 0, recanonization gives the same c",
         {"parts": {"a": a_, "b": b_, "c_invariance_and_class": c_, "c_recanonization": recanon}})
    # C6 Euler operator vs sympy euler_equations
    ctrl("C6", all(c6_check(sectors[S]["X"]) for S in SECTORS),
         "own E_f, E_g on X_S equal sympy.calculus.euler.euler_equations (Function representation), all sectors")
    # C7 exact integrals
    c7a_details, c7a = [], True
    for gphi, gpsi in [(1, 1), (1, -1)]:
        dJ, I = dtJ_squared_poly(gphi, gpsi, TEST_PAIRS[0])
        okk = (dJ != 0) and bool(I > 0)
        c7a = c7a and okk
        c7a_details.append({"g_phi": gphi, "g_psi": gpsi, "regime": "g_phi g_psi > 0" if gphi * gpsi > 0 else "g_phi g_psi < 0",
                            "dtJ_nonzero": dJ != 0, "integral": str(I)})
    ctrl("C7a", c7a, "positivity of the full expression: int (d_t J)^2 > 0 with d_t J != 0, test pair 1, both sign regimes",
         {"cases": c7a_details})
    c7b_details, c7b = [], True
    for ip, pair in enumerate(TEST_PAIRS):
        for S in SECTORS:
            IX = integrate_4field(sectors[S]["X"], pair)
            IR = integrate_4field(sectors[S]["R"], pair)
            eq = sp.simplify(IX - IR) == 0
            c7b = c7b and eq
            c7b_details.append({"pair": ip + 1, "sector": S, "int_X": str(IX), "int_R": str(IR), "equal": eq})
    ctrl("C7b", c7b, "supplementary exact-integral check int X_S = int R_S, two test pairs; necessary, not sufficient; "
         "NOT an independent implementation", {"cases": c7b_details})
    controls["C7_test_fields"] = {"P_f/P_g (times exp(-(t^2+x^2)))":
                                  [{k: str(v).replace("T", "t").replace("X", "x") for k, v in p.items()} for p in TEST_PAIRS]}

    # ---------------------------------------------------------------- report
    gphi, gpsi, M = sp.symbols("g_phi g_psi M")
    pref = {"f4": gphi**2 / (2 * M**4), "f2g2": gphi * gpsi / M**4, "g4": gpsi**2 / (2 * M**4)}
    total_terms = []
    out = {"meta": {"algorithm_version": ALGO_VERSION, "python": platform.python_version(), "sympy": sp.__version__,
                    "script_sha256": script_sha,
                    "status": "DRAFT - unregistered, not a T7a result",
                    "spec": "T7A_STEP4_CLASSIFICATION_ALGORITHM_PROPOSAL.md v0.8.1 (e7ab8d4); canonical rule v0.7.1; option A (v0.6)",
                    "scope_exclusions": "no independent IBP-matrix implementation, no numerical verification, no M^-6, no B4.12 conclusions"},
           "sectors": {}, "controls": controls}
    for S in SECTORS:
        d = sectors[S]
        out["sectors"][S] = {
            "X_S": label_poly(d["X"]),
            "prefactor": str(pref[S]),
            "n_VS": d["n_VS"], "dim": d["dim"],
            "basis": [{"index": k, "triples": [list(t) for t in m], "label": label(m)} for k, m in enumerate(d["B"])],
            "coords": [frac_to_str(v) for v in d["c"]],
            "canonical_representative": label_poly(d["R"]),
            "residual_zero": d["resid_ok"], "rank_checks": d["rank_ok"]}
        for v, m in zip(d["c"], d["B"]):
            if v != 0:
                total_terms.append(f"({sp.sstr(pref[S] * sp.Rational(v.numerator, v.denominator))}) * [{label(m)}]")
    out["total"] = " + ".join(total_terms)
    out["verdict"] = "ALL_CHECKS_PASS" if ok_all else "CHECK_FAILED: " + ", ".join(failures)

    with open(f"{outdir}/t7a_step4_classification.json", "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=1, sort_keys=True, ensure_ascii=False)
        fh.write("\n")
    lines = [f"T7a step 4 symbolic classification - {out['meta']['status']}",
             f"algorithm {ALGO_VERSION}; python {out['meta']['python']}; sympy {out['meta']['sympy']}",
             f"script sha256 {script_sha}", ""]
    for S in SECTORS:
        s = out["sectors"][S]
        lines += [f"[sector {S}] prefactor {s['prefactor']}", f"  X_S = {s['X_S']}",
                  f"  |V_S| = {s['n_VS']}, dim(V_S/~) = {s['dim']}, rank_checks = {s['rank_checks']}, residual_zero = {s['residual_zero']}",
                  "  basis (greedy order):"]
        lines += [f"    b{b['index']}: {b['label']}" for b in s["basis"]]
        lines += [f"  coords c = [{', '.join(s['coords'])}]", f"  canonical representative R_S = {s['canonical_representative']}", ""]
    lines += ["Delta L^(M^-4) in the off-shell basis:", "  " + out["total"], "", "controls:"]
    for k in sorted(controls):
        if "result" in controls[k]:
            lines.append(f"  {k}: {controls[k]['result']} - {controls[k]['description']}")
    lines += ["", out["verdict"]]
    with open(f"{outdir}/t7a_step4_output.txt", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def label_poly(P):
    d = decompose(P)
    if not d:
        return "0"
    parts = []
    for k in sorted(d, key=key):
        parts.append(f"{frac_to_str(d[k])}*[{label(k)}]")
    return " + ".join(parts)


if __name__ == "__main__":
    main()
