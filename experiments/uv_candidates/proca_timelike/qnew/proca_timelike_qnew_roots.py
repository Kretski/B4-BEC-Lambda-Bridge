"""Q_new(c; 1) root arithmetic (V3 input) for docs/PROCA_TIMELIKE_CURRENT_ANALYTIC_SPEC_v0.5.md (87fc31b),
on the 41 registered E3 backgrounds, for sigma = +1 and -1 separately. Base commit 53e2227.

What is computed (no new criterion):
  Q_new(c; 1) = F_eff - c^2 d1 d2 - 4 gamma_u q B(c)/M^4, taken VERBATIM from the registered D5 checker
  (experiments/uv_candidates/proca_timelike/stage1/proca_timelike_d5_addendum_check.py, variable Qnew1; its P1 check
  shows det(lambda I - P_inf) at lambda = -i c equals Q_new(c; 1)). The checker is executed (stdout captured) and must
  report 0 failed checks.
  Conventions (as Stage 2): M = 1; g_phi = +sqrt(p), g_psi = -kappa/sqrt(p); gamma_u = sigma sqrt(Lambda/2).
  The coefficients of Q_new lie in Q(a, b), a = sqrt(p), b = sqrt(Lambda/2); they are computed exactly.

(1) REGISTERED RULE (T3_FUND_PROCA_SPEC 3(iv): "six simple real roots of Q(c;1), tolerance 1e-9 relative"):
    six roots by mpmath.polyroots at 50 digits (maxsteps 500, extraprec 300).
    PRIMARY class (specification text, reviewer decision 2026-11-19): relative tolerance 1e-9, with "relative" read as
    in every other relative 1e-9 test of the registered T3-fund code (C1 consistency, k = 0 check:
    |x - y| / max(1, |y|)):  max_j |Im c_j| / max(1, |c_j|) <= 1e-9 and min_{i<j} |c_i - c_j| / max(1, |c_i|, |c_j|)
    > 1e-9 -> SIX-SIMPLE-REAL; first condition violated -> COMPLEX-ROOT; otherwise REPEATED-REAL-INCONCLUSIVE.
    DIAGNOSTICS (reported, never used for the class): (a) absolute 1e-9, as in the registered t3fund_proca_stage1.py
    line for this rule; (b) pure relative (no floor): |Im c_j| / |c_j| and |c_i - c_j| / max(|c_i|, |c_j|).
    Every background where a diagnostic class differs from the primary class is listed.
(2) CERTIFIED CROSS-CHECK (decides nothing beyond the registered rule; reported separately):
    every real root of Q_new is a root of the exact rational norm polynomial N(c) = prod over (+-a, +-b) of Q_new
    (degree 24). For each rational isolating interval of the square-free part of N, interval arithmetic (outward
    rounding, adaptive precision) certifies either a sign change of Q_new at the endpoints together with Q_new' != 0 on
    the interval (one simple real root), or Q_new != 0 on the interval (no root). Rational roots of N are tested
    exactly (minimal polynomial). If every interval is decided: n simple real roots; n = 6 -> SIX-SIMPLE-REAL,
    n < 6 -> COMPLEX-ROOT (6 - n non-real roots). Otherwise UNDETERMINED (a repeated root cannot be certified this way).
(3) CONTROL: the eigenvalues of P_inf (from the same checker, 50 digits) equal -i times the six roots, relative 1e-9.

NOT computed: V3 PASS/FAIL, condition number of V_inf, Stage 3. A class SIX-SIMPLE-REAL makes the analytic argument
of D5 addendum v0.3 section 6 applicable for that background and sign; it does not by itself classify the model.
Usage (from the repository root): python3 -I experiments/uv_candidates/proca_timelike/qnew/proca_timelike_qnew_roots.py <outdir>
"""
import contextlib
import hashlib
import importlib.util
import io
import json
import platform
import sys
import time
from pathlib import Path

import mpmath as mp
import sympy as sp

STAGE2 = "experiments/uv_candidates/proca_timelike/stage2/proca_timelike_stage2.py"
CHECKER = "experiments/uv_candidates/proca_timelike/stage1/proca_timelike_d5_addendum_check.py"
STAGE1 = "experiments/uv_candidates/proca_timelike/stage1/proca_timelike_stage1_check.py"
EXPECTED_SHA256 = {
    STAGE2: "585143fa50d390b27e40b6aa1d8a2114a32e2b881600d6ee21d0cc22e27a5ff6",
    CHECKER: "2064fa845bde0c4e9494f2806ff30812642d41fe549d1fd98f0fc4903e5fe1ca",
    STAGE1: "b59393d2a6510e4b910d0ea58f7a286921603697ff5de5f4e4992f531133c29f",
    "data/B3_max_speed_7_315096.json": "99a6b04bc34821ddd5c5d58a78b3474d369aa9fc0ee264163096398b0028ccbd",
    "data/e3_manifest.json": "1bae892a08c3de1f928e4057ee0ad801294a90462478aaedf5b9c724e23d87e0",
}
DPS = 50
TOL = mp.mpf('1e-9')                 # registered tolerance (T3_FUND_PROCA_SPEC 3(iv))
SIGNS = (+1, -1)
# (interval digits, root-enclosure width 10^-n): the width is kept well above the working precision so that the
# endpoint values of Q_new (of size ~ width * |Q_new'|) stay resolvable; precision and width are raised together
LEVELS = ((DPS + 20, 20), (2*DPS + 20, 40), (4*DPS + 20, 80))
A_, B_ = sp.symbols('a_sqrt_p b_sqrt_Lambda_half', positive=True)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_hashes():
    for pth, h in EXPECTED_SHA256.items():
        if sha(pth) != h:
            raise SystemExit("ABORT: sha256 mismatch for %s" % pth)


def load_stage2():
    spec = importlib.util.spec_from_file_location("proca_stage2", STAGE2)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def load_checker():
    ns = {'__file__': str(Path(CHECKER)), '__name__': 'proca_d5_checker'}
    buf = io.StringIO()
    code = None
    with contextlib.redirect_stdout(buf):
        try:
            exec(compile(Path(CHECKER).read_text(encoding='utf-8'), CHECKER, 'exec'), ns)
        except SystemExit as e_:          # the checker ends with sys.exit(0 if all checks pass else 1)
            code = e_.code
    txt = buf.getvalue()
    if code not in (None, 0):
        raise SystemExit("ABORT: registered D5 checker exited with code %r" % (code,))
    if 'TOTAL 13 checks, 0 failed' not in txt:
        raise SystemExit("ABORT: registered D5 checker did not report 13 checks, 0 failed")
    return ns


# ---------------- exact coefficients in Q(a, b) ----------------
def reduce_ab(e, p, L2):
    """Reduce a polynomial expression in a, b with a^2 = p, b^2 = L2 to r0 + r1 a + r2 b + r3 a b (rationals)."""
    P = sp.Poly(sp.expand(e), A_, B_)
    r = [sp.Integer(0)]*4
    for (i, j), co in P.terms():
        f = sp.Rational(co)*p**(i//2)*L2**(j//2)
        r[(i % 2) + 2*(j % 2)] += f
    return r          # index 0: 1, 1: a, 2: b, 3: a b


def qnew_exact(ns, rec, sg):
    """Coefficients (highest degree first) of Q_new(c;1) as 4-tuples of rationals over (1, a, b, ab)."""
    p, kap, L = rec['p'], rec['kappa'], rec['Lambda']
    L2 = L/2
    subs = {ns['M']: 1, ns['gp']: A_, ns['gs']: -kap*A_/p, ns['gu']: sg*B_,
            ns['w']: rec['w'], ns['q']: rec['q'], ns['v']: rec['v'], ns['r']: rec['r']}
    Q = sp.Poly(sp.expand(ns['Qnew1'].subs(subs)), ns['c'])
    if Q.degree() != 6:
        raise SystemExit("ABORT: Q_new is not of degree 6")
    if sp.simplify(Q.LC() + 1) != 0:
        raise SystemExit("ABORT: Q_new leading coefficient is not exactly -1")
    return [reduce_ab(co, p, L2) for co in Q.all_coeffs()], p, L2


def to_sym(r4, p, L2):
    return r4[0] + r4[1]*sp.sqrt(p) + r4[2]*sp.sqrt(L2) + r4[3]*sp.sqrt(p)*sp.sqrt(L2)


def norm_poly(coeffs, p, L2, c):
    """N(c) = prod over the four conjugates (a -> +-a, b -> +-b): an exact rational polynomial."""
    Qab = sum((r4[0] + r4[1]*A_ + r4[2]*B_ + r4[3]*A_*B_)*c**(len(coeffs) - 1 - i) for i, r4 in enumerate(coeffs))
    Nab = sp.expand(Qab*Qab.subs(A_, -A_)*Qab.subs(B_, -B_)*Qab.subs({A_: -A_, B_: -B_}, simultaneous=True))
    P = sp.Poly(Nab, c)
    out = []
    for co in P.all_coeffs():
        r = reduce_ab(co, p, L2)
        if any(x_ != 0 for x_ in r[1:]):
            raise SystemExit("ABORT: norm polynomial is not rational")
        out.append(r[0])
    N = sp.Poly(out, c)
    if N.is_zero or N.degree() != 24:
        raise SystemExit("ABORT: norm polynomial is not of degree 24")
    return N


# ---------------- registered numerical rule ----------------
def registered_rule(coeffs, p, L2):
    with mp.workdps(DPS):
        a, b = mp.sqrt(mp.mpf(sp.Rational(p).p)/sp.Rational(p).q), mp.sqrt(mp.mpf(sp.Rational(L2).p)/sp.Rational(L2).q)
        def num(x_):
            x_ = sp.Rational(x_)
            return mp.mpf(int(x_.p))/mp.mpf(int(x_.q))
        co = [num(r4[0]) + num(r4[1])*a + num(r4[2])*b + num(r4[3])*a*b for r4 in coeffs]
        rt = mp.polyroots(co, maxsteps=500, extraprec=300)
        pairs = [(i, j) for i in range(6) for j in range(i)]
        im_abs = max(abs(mp.im(z)) for z in rt)
        sep_abs = min(abs(rt[i] - rt[j]) for i, j in pairs)
        im_rel = max(abs(mp.im(z))/max(1, abs(z)) for z in rt)
        sep_rel = min(abs(rt[i] - rt[j])/max(1, abs(rt[i]), abs(rt[j])) for i, j in pairs)
        # exact: c = 0 is a root iff the constant term of Q_new is exactly 0 in Q(a, b) (minimal-polynomial test;
        # the four rational components alone do not decide it when p or Lambda/2 is a rational square)
        zero_root = exact_zero(coeffs, sp.Integer(0), p, L2)
        im_pure = None if zero_root else max(abs(mp.im(z))/abs(z) for z in rt)
        sep_pure = None if zero_root else min(abs(rt[i] - rt[j])/max(abs(rt[i]), abs(rt[j])) for i, j in pairs)

        def cls(mi, se):
            if mi is None:
                return 'NOT-DEFINED (root exactly 0)'
            return 'COMPLEX-ROOT' if mi > TOL else ('SIX-SIMPLE-REAL' if se > TOL else 'REPEATED-REAL-INCONCLUSIVE')
        c_pri, c_abs, c_pure = cls(im_rel, sep_rel), cls(im_abs, sep_abs), cls(im_pure, sep_pure)
        roots = sorted(rt, key=lambda z: (float(mp.re(z)), float(mp.im(z))))
        nd = [min(abs(z - roots[j])/max(1, abs(z), abs(roots[j])) for j in range(6) if j != i) for i, z in enumerate(roots)]
        return dict(roots=[[mp.nstr(mp.re(z), 20), mp.nstr(mp.im(z), 5)] for z in roots],
                    nearest_other_root_rel_distance=[mp.nstr(x_, 6) for x_ in nd],
                    n_real=sum(1 for z in rt if abs(mp.im(z))/max(1, abs(z)) <= TOL),
                    n_complex=sum(1 for z in rt if abs(mp.im(z))/max(1, abs(z)) > TOL),
                    max_rel_imag=mp.nstr(im_rel, 6), min_rel_root_distance=mp.nstr(sep_rel, 6),
                    max_abs_imag=mp.nstr(im_abs, 6), min_abs_root_distance=mp.nstr(sep_abs, 6),
                    max_pure_rel_imag=(mp.nstr(im_pure, 6) if im_pure is not None else None),
                    min_pure_rel_root_distance=(mp.nstr(sep_pure, 6) if sep_pure is not None else None),
                    class_primary_spec_relative=c_pri,
                    diag_class_registered_code_absolute=c_abs, diag_class_pure_relative=c_pure,
                    diag_absolute_agrees=bool(c_abs == c_pri), diag_pure_relative_agrees=bool(c_pure == c_pri)), rt, co


# ---------------- certified cross-check ----------------
def ivq(x_):
    x_ = sp.Rational(x_)
    return mp.iv.mpf(int(x_.p))/mp.iv.mpf(int(x_.q))


def ivsign(X):
    if X.a > 0:
        return 1
    if X.b < 0:
        return -1
    if X.a == 0 and X.b == 0:
        return 0
    return None


def iv_poly(coeffs, X, ia, ib):
    acc = mp.iv.mpf(0)
    for r4 in coeffs:
        acc = acc*X + (ivq(r4[0]) + ivq(r4[1])*ia + ivq(r4[2])*ib + ivq(r4[3])*ia*ib)
    return acc


def deriv(coeffs):
    n = len(coeffs) - 1
    return [[x_*(n - i) for x_ in r4] for i, r4 in enumerate(coeffs[:-1])]


def exact_zero(coeffs, x_, p, L2):
    val = sum(to_sym(r4, p, L2)*x_**(len(coeffs) - 1 - i) for i, r4 in enumerate(coeffs))
    X = sp.Symbol('X')
    return sp.minimal_polynomial(val, X) == X


def certified(coeffs, p, L2, c):
    N = norm_poly(coeffs, p, L2, c)
    Ns = sp.Poly(sp.quo(N, sp.gcd(N, N.diff(c))), c)
    dco = deriv(coeffs)
    ivs = sorted(Ns.intervals(), key=lambda t_: (t_[0][0], t_[0][1]))
    n_simple, undecided, levels_used = 0, 0, []
    for (a0, b0), _ in ivs:
        a0, b0 = sp.Rational(a0), sp.Rational(b0)
        if a0 == b0:
            if exact_zero(coeffs, a0, p, L2):
                if exact_zero(dco, a0, p, L2):
                    undecided += 1            # exactly repeated rational root: not a simple root
                else:
                    n_simple += 1
            levels_used.append(-1)
            continue
        decided = False
        for lvl, (ivdps, nexp) in enumerate(LEVELS):
            lo, hi = Ns.refine_root(a0, b0, eps=sp.Rational(1, 10**nexp))
            lo, hi = sp.Rational(lo), sp.Rational(hi)
            if Ns.count_roots(lo, hi) != 1 or not (a0 <= lo <= hi <= b0):
                raise SystemExit("ABORT: refinement lost the root")
            if lo != hi and (Ns.eval(lo) == 0 or Ns.eval(hi) == 0):
                continue                      # an endpoint is itself a root of N: refine further at the next level
            old = mp.iv.dps
            mp.iv.dps = ivdps
            try:
                ia, ib = mp.iv.sqrt(ivq(p)), mp.iv.sqrt(ivq(L2))
                X = mp.iv.mpf([ivq(lo).a, ivq(hi).b])
                slo = ivsign(iv_poly(coeffs, mp.iv.mpf([ivq(lo).a, ivq(lo).b]), ia, ib)) if lo != hi else None
                shi = ivsign(iv_poly(coeffs, mp.iv.mpf([ivq(hi).a, ivq(hi).b]), ia, ib)) if lo != hi else None
                if lo == hi:                  # refinement hit the rational root exactly
                    if exact_zero(coeffs, lo, p, L2):
                        if not exact_zero(dco, lo, p, L2):
                            n_simple += 1
                            decided = True
                    else:
                        decided = True
                elif slo not in (None, 0) and shi not in (None, 0):
                    if slo != shi:
                        if ivsign(iv_poly(dco, X, ia, ib)) not in (None, 0):
                            n_simple += 1
                            decided = True
                    elif ivsign(iv_poly(coeffs, X, ia, ib)) not in (None, 0):
                        decided = True
            finally:
                mp.iv.dps = old
            if decided:
                levels_used.append(lvl)
                break
        if not decided:
            undecided += 1
            levels_used.append(None)
    if undecided:
        cls = 'UNDETERMINED'
    elif n_simple == 6:
        cls = 'SIX-SIMPLE-REAL'
    elif n_simple < 6:
        cls = 'COMPLEX-ROOT'
    else:
        raise SystemExit("ABORT: more than six certified roots")
    return dict(class_certified=cls, n_simple_real_certified=n_simple, n_nonreal_certified=(6 - n_simple if not undecided else None),
                n_intervals_norm_poly=len(ivs), n_undecided_intervals=undecided,
                max_level_used=max((x_ for x_ in levels_used if x_ is not None), default=-1))


# ---------------- control: eigenvalues of P_inf ----------------
def control_pinf(ns, rec, sg, rt):
    with mp.workdps(DPS):
        def num(x_):
            x_ = sp.Rational(x_)
            return mp.mpf(int(x_.p))/mp.mpf(int(x_.q))
        P, KAP, L = num(rec['p']), num(rec['kappa']), num(rec['Lambda'])
        gp = mp.sqrt(P)
        vals = {ns['M']: 1, ns['gp']: gp, ns['gs']: -KAP/gp, ns['gu']: sg*mp.sqrt(L/2),
                ns['w']: num(rec['w']), ns['q']: num(rec['q']), ns['v']: num(rec['v']), ns['r']: num(rec['r'])}
        syms = list(vals)
        f = sp.lambdify(syms, ns['Pinf'], 'mpmath')
        Pm = mp.matrix(f(*[vals[s_] for s_ in syms]))
        ev = mp.eig(Pm, left=False, right=False)
        cs = [1j*z for z in ev]                     # eigenvalue -i c  ->  c = i * eigenvalue
        return match_dev(cs, list(rt))


def match_dev(cs, rt):
    """Optimal (minimax) matching of eigenvalue-derived speeds to roots over all permutations (6! = 720):
    min over permutations of max_j |cs[perm j] - rt[j]| / max(1, |rt[j]|)."""
    import itertools
    best = None
    for perm in itertools.permutations(range(len(rt))):
        d_ = max(abs(cs[perm[j]] - rt[j])/max(1, abs(rt[j])) for j in range(len(rt)))
        if best is None or d_ < best:
            best = d_
    return best


def main(outdir):
    T0 = time.time()
    check_hashes()
    st2 = load_stage2()
    ns = load_checker()
    cases = st2.load()                               # hash-checked inputs, rationalised, rank-one ABORT rule
    c = ns['c']
    rows = []
    worst = mp.mpf(0)
    for tag, rec in cases:
        row = dict(tag=tag, cls=rec['cls'], mode=rec['mode'], delta_rank1=str(sp.N(rec['delta_rank1'], 10)))
        for sg in SIGNS:
            coeffs, p, L2 = qnew_exact(ns, rec, sg)
            reg, rt, _ = registered_rule(coeffs, p, L2)
            cert = certified(coeffs, p, L2, c)
            dev = control_pinf(ns, rec, sg, rt)
            worst = max(worst, dev)
            reg.update(cert)
            reg['primary_and_certified_agree'] = (None if cert['class_certified'] == 'UNDETERMINED'
                                                  else bool(cert['class_certified'] == reg['class_primary_spec_relative']))
            reg['control_eig_Pinf_vs_roots_max_rel'] = mp.nstr(dev, 5)
            reg['control_ok'] = bool(dev <= TOL)
            row[str(sg)] = reg
        rows.append(row)
    summ = {}
    for sg in SIGNS:
        k_ = str(sg)
        summ[k_] = dict(
            primary_spec_relative={x_: sum(1 for r_ in rows if r_[k_]['class_primary_spec_relative'] == x_)
                                   for x_ in ('SIX-SIMPLE-REAL', 'COMPLEX-ROOT', 'REPEATED-REAL-INCONCLUSIVE')},
            diag_registered_code_absolute={x_: sum(1 for r_ in rows if r_[k_]['diag_class_registered_code_absolute'] == x_)
                                           for x_ in ('SIX-SIMPLE-REAL', 'COMPLEX-ROOT', 'REPEATED-REAL-INCONCLUSIVE')},
            diag_pure_relative={x_: sum(1 for r_ in rows if r_[k_]['diag_class_pure_relative'] == x_)
                                for x_ in ('SIX-SIMPLE-REAL', 'COMPLEX-ROOT', 'REPEATED-REAL-INCONCLUSIVE',
                                           'NOT-DEFINED (root exactly 0)')},
            certified={x_: sum(1 for r_ in rows if r_[k_]['class_certified'] == x_)
                       for x_ in ('SIX-SIMPLE-REAL', 'COMPLEX-ROOT', 'UNDETERMINED')},
            diag_absolute_disagreements=[r_['tag'] for r_ in rows if not r_[k_]['diag_absolute_agrees']],
            diag_pure_relative_disagreements=[r_['tag'] for r_ in rows if not r_[k_]['diag_pure_relative_agrees']],
            primary_certified_disagreements=[r_['tag'] for r_ in rows if r_[k_]['primary_and_certified_agree'] is False],
            control_failures=[r_['tag'] for r_ in rows if not r_[k_]['control_ok']])
    result = dict(
        spec="docs/PROCA_TIMELIKE_CURRENT_ANALYTIC_SPEC_v0.5.md (87fc31b)", base="53e2227",
        stage="Q_new(c;1) root arithmetic (V3 input only)",
        not_included="V3 PASS/FAIL, condition number of V_inf, Stage 3: not authorised / not computed",
        inputs={k_: sha(k_) for k_ in EXPECTED_SHA256}, script_sha256=sha(__file__),
        environment=dict(python=platform.python_version(), sympy=sp.__version__, mpmath=mp.__version__),
        conventions="M = 1; g_phi = +sqrt(p); g_psi = -kappa/sqrt(p); gamma_u = sigma sqrt(Lambda/2), sigma = +1, -1",
        rule=("primary: relative 1e-9 (|x - y| / max(1, |y|), as in the registered T3-fund relative tests), 50 digits; "
              "diagnostics: absolute 1e-9 (registered code line for this rule) and pure relative"),
        control_eig_Pinf_vs_roots_max_rel=mp.nstr(worst, 5), summary=summ,
        class_note=("SIX-SIMPLE-REAL makes the D5 addendum v0.3 section 6 argument applicable for that background and sign; "
                    "COMPLEX-ROOT meets NOT-VIABLE (d) for that background and sign. Model class not changed here."),
        backgrounds=rows, seconds=round(time.time() - T0, 1))
    od = Path(outdir)
    od.mkdir(parents=True, exist_ok=True)
    (od / "proca_timelike_qnew_roots.json").write_text(json.dumps(result, indent=1, default=str) + "\n", encoding='utf-8', newline='\n')
    print(json.dumps(dict(summary=summ, control_max_rel=result['control_eig_Pinf_vs_roots_max_rel']), indent=1))
    bad = any(summ[k_]['control_failures'] for k_ in summ)
    return 1 if bad else 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(64)
    sys.exit(main(sys.argv[1]))
