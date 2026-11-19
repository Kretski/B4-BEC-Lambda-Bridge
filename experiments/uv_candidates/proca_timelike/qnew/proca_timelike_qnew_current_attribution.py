"""Narrow attribution check for NOT-VIABLE (d) (reviewer decision 2026-11-19 05:10): does the complex pair of
Q_new(c;1) appear when the timelike current is switched on, in ONE AND THE SAME formulation of the principal symbol?
Only the four backgrounds rank1_minus_stable_5713, rank1_minus_boundary_2020, rank1_plus_boundary_329,
rank1_plus_boundary_3437, for sigma = +1 and -1. No population run, no V3 PASS/FAIL, no Stage 3.

Decomposition (exact, from the registered D5 checker): Q_new(c;1) = Q0(c) - 4 gamma_u q B(c)/M^4 with
  Q0(c) = F_eff(c) - c^2 d1 d2,  F_eff and d1 built with Lambda_eff = 2 gamma_u^2/M^2 = Lambda (registered).
Q0 is the principal-symbol polynomial of the same reduced system WITHOUT the exchange term of the current; with
Lambda_eff = Lambda it is, coefficient by coefficient, the registered T3-fund Q(c;1) (bare Lambda, no current), which
the check verifies against the registered code itself (t3fund_proca_stage1.Qcoeffs(P, eps = 1)).

Per background and sign:
 (A) identity: max relative coefficient deviation |Q0 - Q_reg| / max(1, |Q_reg|) over the 7 coefficients, Q_reg from
     the registered code; must be <= 1e-9 (expected ~1e-16: the registered code uses the stored s, Q0 uses kappa^2/p).
 (B) class of Q0 (t = 0) and of Q_new (t = 1): primary rule (relative 1e-9, 50 digits) and the certified cross-check,
     both from the reviewed Q_new script (91ea9152), plus the registered T3-fund strictly_hyperbolic flag for the case.
 (C) homotopy Q_t = Q0 - t 4 gamma_u q B(c), t in [0, 1]: number of real roots (primary rule) on 401 equally spaced t;
     each change of the count between neighbouring grid points is located by bisection in t to 1e-15 and reported
     with the roots involved. These are the OBSERVED transitions along the discretely tracked homotopy: two changes
     between neighbouring grid points can be missed, and polyroots may lose accuracy next to a multiple root (grid
     points where it does not converge are listed; a change of the count whose bisection hits a non-converging point
     is listed separately with its last bisection interval and is not reported as localised). Completeness would need a separate analysis of the discriminant
     in t. Diagnostic only; the endpoint classes (B) are certified.
 (D) diagnostic only: Q_new at gamma_u -> 0 literally (Lambda_eff -> 0 as well), primary and certified class.
Attribution per background and sign (reviewer rule, 2026-11-19 05:16):
  (A) fails                                                                      -> NOT ESTABLISHED
  (A) holds, registered flag consistent with the class of Q0 (both strictly hyperbolic), Q0 SIX/SIX (primary and
      certified), Q_new COMPLEX/COMPLEX                                          -> CURRENT-INDUCED
  (A) holds, Q_new SIX/SIX                                                       -> NOT-INDUCED: hyperbolic with the current
  every other combination (including a Q0 class inconsistent with the registered T3-fund flag) -> NOT ESTABLISHED,
  with the reason recorded. No model class is assigned here.
Usage (from the repository root): python3 -I proca_timelike_qnew_current_attribution.py <qnew script> <outdir>
"""
import hashlib
import importlib.util
import json
import platform
import sys
import time
from pathlib import Path

import mpmath as mp
import sympy as sp

TAGS = ("rank1_minus_stable_5713", "rank1_minus_boundary_2020", "rank1_plus_boundary_329", "rank1_plus_boundary_3437")
SIGNS = (+1, -1)
REG_DIR = "t3fund_proca_stage1_impl_bundle"
REG_RESULT = "t3fund_proca_stage1_result_bundle/proca_stage1_result.json"
EXPECTED_SHA256 = {
    "data/B3_max_speed_7_315096.json": "99a6b04bc34821ddd5c5d58a78b3474d369aa9fc0ee264163096398b0028ccbd",
    "data/e3_manifest.json": "1bae892a08c3de1f928e4057ee0ad801294a90462478aaedf5b9c724e23d87e0",
    REG_DIR + "/t3fund_proca_stage1.py": "7be1486a4fde2fe5ae8480816d5d46ba76f54e3e3b98d77bca708a70142482ff",
    REG_DIR + "/smeta_ref.py": "6f915e9b93c29e8e863b8e26a5967e687a2ef643a26b67c28805d9bbeb9cd7a3",
    REG_RESULT: "27c9af3267fee535f368543dc64c0d65e05edc3c8f3a20546176958952726614",
}
QNEW_SHA = "91ea9152b1aef0277369044d9e5d04bbd68b5f6d29d83eaec3db01b68b795a30"
DPS = 50
TOL = mp.mpf('1e-9')
NT = 401


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def coeffs_of(qn, ns, expr, rec, sg, gu_override=None):
    """Exact coefficients (highest first) over (1, a, b, ab) of expr(c) at the background, M = 1."""
    p, kap = rec['p'], rec['kappa']
    L2 = rec['Lambda']/2
    gu = sg*qn.B_ if gu_override is None else gu_override
    subs = {ns['M']: 1, ns['gp']: qn.A_, ns['gs']: -kap*qn.A_/p, ns['gu']: gu,
            ns['w']: rec['w'], ns['q']: rec['q'], ns['v']: rec['v'], ns['r']: rec['r']}
    P = sp.Poly(sp.expand(expr.subs(subs)), ns['c'])
    return [qn.reduce_ab(co, p, L2) for co in P.all_coeffs()], p, L2


def num_coeffs(cf, p, L2):
    a, b = mp.sqrt(mp.mpf(sp.Rational(p).p)/sp.Rational(p).q), mp.sqrt(mp.mpf(sp.Rational(L2).p)/sp.Rational(L2).q)
    def nm(x_):
        x_ = sp.Rational(x_)
        return mp.mpf(int(x_.p))/mp.mpf(int(x_.q))
    return [nm(r[0]) + nm(r[1])*a + nm(r[2])*b + nm(r[3])*a*b for r in cf]


def n_real(co):
    """Number of real roots by the primary rule; (None, None) if polyroots does not converge."""
    try:
        rt = mp.polyroots(co, maxsteps=500, extraprec=300)
    except mp.libmp.NoConvergence:
        return None, None
    return sum(1 for z in rt if abs(mp.im(z))/max(1, abs(z)) <= TOL), rt


def attribute(qn, ns, reg, rec, sg, registered_strictly_hyperbolic):
    """Checks (A)-(D) and the attribution for one background and one sign."""
    c = ns['c']
    Bterm = 4*ns['gu']*ns['q']*ns['Bc']/ns['M']**4
    Q0expr = sp.expand(ns['Qnew1'] + Bterm)
    out = {}
    raw = rec['raw']
    cf1, p, L2 = coeffs_of(qn, ns, ns['Qnew1'], rec, sg)
    cf0, _, _ = coeffs_of(qn, ns, Q0expr, rec, sg)
    cfB, _, _ = coeffs_of(qn, ns, Bterm, rec, sg)
    with mp.workdps(DPS):
        # (A) identity with the registered code
        case = dict(p=raw['p'], s=raw['s'], kappa=raw['kappa'], Lambda=raw['Lambda'], w=raw['w'], q=raw['q'],
                    v=raw['v'], r=raw['r'], gphi=mp.sqrt(mp.mpf(str(raw['p']))),
                    gpsi=-mp.mpf(str(raw['kappa']))/mp.sqrt(mp.mpf(str(raw['p']))))
        P = reg.pars(case, DPS)
        qreg = [mp.re(z) for z in reg.Qcoeffs(P, mp.mpf(1))][::-1]        # highest first
        q0n = num_coeffs(cf0, p, L2)
        if len(qreg) != len(q0n):
            raise SystemExit("ABORT: degree mismatch with the registered Q(c;1)")
        devA = max(abs(x_ - y_)/max(1, abs(y_)) for x_, y_ in zip(q0n, qreg))
        out['A_identity_Q0_vs_registered_Q_max_rel'] = mp.nstr(devA, 5)
        out['A_identity_ok'] = bool(devA <= TOL)
        # (B) classes at t = 0 and t = 1
        for lab, cf in (('t0_Q0_no_exchange_term', cf0), ('t1_Q_new', cf1)):
            rr, _, _ = qn.registered_rule(cf, p, L2)
            ce = qn.certified(cf, p, L2, c)
            out['B_' + lab] = dict(class_primary=rr['class_primary_spec_relative'], class_certified=ce['class_certified'],
                                  n_real_primary=rr['n_real'], min_rel_root_distance=rr['min_rel_root_distance'],
                                  max_rel_imag=rr['max_rel_imag'], roots=rr['roots'])
        # (C) homotopy in t
        bn = num_coeffs(cfB, p, L2)
        if len(bn) > len(q0n):
            raise SystemExit("ABORT: exchange term of higher degree than Q0")
        bn = [mp.mpf(0)]*(len(q0n) - len(bn)) + bn                 # align degrees (highest first): B term has degree 4
        if max(abs(x_ - y_ - z_)/max(1, abs(z_)) for x_, y_, z_ in zip(q0n, bn, num_coeffs(cf1, p, L2))) > TOL:
            raise SystemExit("ABORT: Q0 - exchange term does not reproduce Q_new")
        def co_t(t_):
            return [x_ - t_*y_ for x_, y_ in zip(q0n, bn)]
        grid = [mp.mpf(i)/(NT - 1) for i in range(NT)]
        counts = [n_real(co_t(t_))[0] for t_ in grid]
        nonconv = [mp.nstr(grid[i], 6) for i in range(NT) if counts[i] is None]
        trans = []
        trans_nc = []
        for i in range(NT - 1):
            if counts[i] is None or counts[i + 1] is None:
                continue
            if counts[i] != counts[i + 1]:
                lo, hi = grid[i], grid[i + 1]
                nlo = counts[i]
                converged = True
                while hi - lo > mp.mpf('1e-15'):
                    mid = (lo + hi)/2
                    nm_ = n_real(co_t(mid))[0]
                    if nm_ is None:
                        converged = False
                        break
                    if nm_ == nlo:
                        lo = mid
                    else:
                        hi = mid
                if not converged:
                    # the change of the count between the grid points stays recorded, but it is NOT a localised transition
                    trans_nc.append(dict(grid_interval=[mp.nstr(grid[i], 15), mp.nstr(grid[i + 1], 15)],
                                         n_real_at_grid_points=[counts[i], counts[i + 1]],
                                         last_bisection_interval=[mp.nstr(lo, 17), mp.nstr(hi, 17)],
                                         reason='polyroots did not converge during the bisection'))
                    continue
                _, rt = n_real(co_t(hi))
                _, rtl = n_real(co_t(lo))
                reals_lo = sorted(mp.re(z) for z in (rtl or []) if abs(mp.im(z))/max(1, abs(z)) <= TOL)
                cplx_hi = [z for z in (rt or []) if abs(mp.im(z))/max(1, abs(z)) > TOL]
                trans.append(dict(t_star=mp.nstr((lo + hi)/2, 15), n_real_before=counts[i], n_real_after=counts[i + 1],
                                  complex_roots_after=[[mp.nstr(mp.re(z), 12), mp.nstr(mp.im(z), 5)] for z in cplx_hi],
                                  real_roots_before=[mp.nstr(x_, 12) for x_ in reals_lo]))
        out['C_homotopy'] = dict(n_grid=NT, n_real_at_t0=counts[0], n_real_at_t1=counts[-1],
                                 observed_transitions_discrete_tracking=trans,
                                 transition_bisection_without_convergence=trans_nc,
                                 grid_points_without_convergence=nonconv,
                                 note=("observed transitions along the discretely tracked homotopy; not a proof that all "
                                       "transitions were found (that needs the discriminant in t); diagnostic only"))
        # (D) literal gamma_u -> 0 (Lambda_eff -> 0 too), diagnostic only
    cfD, _, _ = coeffs_of(qn, ns, ns['Qnew1'], rec, sg, gu_override=0)
    with mp.workdps(DPS):
        rr, _, _ = qn.registered_rule(cfD, p, L2)
    ce = qn.certified(cfD, p, L2, c)
    out['D_diag_literal_gamma_u_zero'] = dict(class_primary=rr['class_primary_spec_relative'], class_certified=ce['class_certified'])
    b0, b1 = out['B_t0_Q0_no_exchange_term'], out['B_t1_Q_new']
    agree_reg = bool((b0['class_primary'] == 'SIX-SIMPLE-REAL') == bool(registered_strictly_hyperbolic))
    out['t0_class_agrees_with_registered_t3fund'] = agree_reg
    q0_six = b0['class_primary'] == b0['class_certified'] == 'SIX-SIMPLE-REAL'
    q1_six = b1['class_primary'] == b1['class_certified'] == 'SIX-SIMPLE-REAL'
    q1_cx = b1['class_primary'] == b1['class_certified'] == 'COMPLEX-ROOT'
    if not out['A_identity_ok']:
        att, why = 'NOT ESTABLISHED', 'identity (A) with the registered Q(c;1) failed'
    elif q0_six and q1_cx and agree_reg and registered_strictly_hyperbolic:
        att, why = 'CURRENT-INDUCED', ''
    elif q1_six:
        att, why = 'NOT-INDUCED: hyperbolic with the current', ''
    else:
        reasons = []
        if not agree_reg:
            reasons.append('class of Q0 inconsistent with the registered T3-fund strictly_hyperbolic flag')
        if not q0_six:
            reasons.append('Q0 not SIX-SIMPLE-REAL by both the primary rule and the certified check')
        if not q1_cx:
            reasons.append('Q_new not COMPLEX-ROOT (or not SIX-SIMPLE-REAL) by both the primary rule and the certified check')
        att, why = 'NOT ESTABLISHED', '; '.join(reasons) or 'other combination'
    out['attribution'] = att
    out['attribution_reason'] = why
    return out


def main(qnew_path, outdir):
    T0 = time.time()
    for pth, h in EXPECTED_SHA256.items():
        if sha(pth) != h:
            raise SystemExit("ABORT: sha256 mismatch for %s" % pth)
    if sha(qnew_path) != QNEW_SHA:
        raise SystemExit("ABORT: Q_new script is not the reviewed version 91ea9152")
    qn = load_module(qnew_path, "qnew")
    qn.check_hashes()
    st2 = qn.load_stage2()
    ns = qn.load_checker()
    sys.path.insert(0, str(Path(REG_DIR).resolve()))
    reg = load_module(REG_DIR + "/t3fund_proca_stage1.py", "t3fund_reg")
    regres = {(x['case']['mode'], x['case']['label']): x for x in json.load(open(REG_RESULT, encoding='utf-8'))
              if x['case']['pop'] == 'PE'}
    cases = dict(st2.load())
    rows = []
    for tag in TAGS:
        rec = cases[tag]
        parts = tag.split('_')
        rkey = (parts[0] + '_' + parts[1], 'E3_' + '_'.join(parts[2:]))
        row = dict(tag=tag, registered_t3fund_outcome=regres[rkey]['outcome'],
                   registered_t3fund_strictly_hyperbolic=regres[rkey]['uv']['strictly_hyperbolic'])
        for sg in SIGNS:
            row[str(sg)] = attribute(qn, ns, reg, rec, sg, row['registered_t3fund_strictly_hyperbolic'])
        rows.append(row)
    summary = {str(sg): {r_['tag']: r_[str(sg)]['attribution'] for r_ in rows} for sg in SIGNS}
    result = dict(stage="narrow attribution check for NOT-VIABLE (d); four backgrounds x sigma = +-1",
                  base="53e2227 (+ reviewed Q_new script 91ea9152)", inputs={k_: sha(k_) for k_ in EXPECTED_SHA256},
                  qnew_script_sha256=sha(qnew_path), script_sha256=sha(__file__),
                  environment=dict(python=platform.python_version(), sympy=sp.__version__, mpmath=mp.__version__),
                  summary=summary, all_identity_ok=all(r_[str(sg)]['A_identity_ok'] for r_ in rows for sg in SIGNS),
                  backgrounds=rows, seconds=round(time.time() - T0, 1),
                  note="No model class is assigned here; the decision on (d) is the reviewer's.")
    od = Path(outdir)
    od.mkdir(parents=True, exist_ok=True)
    (od / "proca_timelike_qnew_current_attribution.json").write_text(json.dumps(result, indent=1, default=str) + "\n",
                                                                       encoding='utf-8', newline='\n')
    print(json.dumps(dict(summary=summary, all_identity_ok=result['all_identity_ok']), indent=1))
    return 0 if result['all_identity_ok'] else 1


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(64)
    sys.exit(main(sys.argv[1], sys.argv[2]))
