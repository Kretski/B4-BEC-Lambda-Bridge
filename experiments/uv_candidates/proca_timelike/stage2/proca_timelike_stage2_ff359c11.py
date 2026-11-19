"""Stage 2 (arithmetic only) for docs/PROCA_TIMELIKE_CURRENT_ANALYTIC_SPEC_v0.5.md (registered 87fc31b).
Evaluates, on the 41 registered E3 backgrounds and for both signs of gamma_u, the Stage 1 formulas:
  (i)  D4: the extended common-M condition (c'_max = max(sqrt p, |kappa|/sqrt p, sqrt(Lambda/2)), G' = G),
       sets A, B, C and the BORDERLINE rule of the registered common-M specification (3cba04e), unchanged;
  (ii) D3: k_der per real leading-order branch, with R-UNDEFINED, no-crossing branches, poles and H3 reported separately.
No new derivation and no new criterion. NOT included (not authorised): roots of Q_new(c;1) / V3 / NOT-VIABLE (d).

Fixed conventions (from the Stage 1 report):
  M = 1 (all quantities are ratios to M); g_phi = +sqrt(p) (T7a matching), g_psi = -kappa/sqrt(p), so g_phi g_psi = -kappa;
  gamma_u = sigma sqrt(Lambda/2), sigma = +1 and -1 evaluated separately (Lambda_eff = registered Lambda);
  only the relative sign of gamma_u and g_phi is physical; it is fixed by sigma with g_phi > 0.
  K_lead = B4 kernel with (p, s_eff = kappa^2/p, kappa, Lambda); F_eff(c) = det Khat(c); real branches = real roots of F_eff.
  N_i = |uhat(c_i)|^2 / |tr Khat(c_i)|, uhat = (2 gamma_u q - c g_phi (w - c q), -c g_psi (v - c r)).
  Crossing iff N_i + c_i^2 > 1; then k_der,i/M = (N_i + c_i^2 - 1)^(-1/2). Pole on the branch: k/M = (c_i^2 - 1)^(-1/2) for |c_i| > 1.
  V2 on a background iff max over real branches of (N_i + c_i^2) < 2 (branches without crossing satisfy it automatically).
  R-UNDEFINED iff tr Khat and F_eff share a real root: decided EXACTLY (rationalised inputs, g = gcd(F_eff, tr Khat);
  per distinct real root of F_eff a rational isolating interval; R-UNDEFINED iff g has a root in it, decided by an
  exact Sturm count on the open interval or exact evaluation at a degenerate interval; no tolerance matching).
  V2 / H3 are reported as None (undetermined) when a background has no real branch (no vacuous all([])).
  k_der status per background and sign: 'defined' | 'no crossing on all real branches' |
  'undetermined (R-UNDEFINED)' | 'undetermined (no real branches)'.
  H3 on branch i iff gamma_u q [gamma_u q - c_i g_phi (w - c_i q)] <= 0.
Usage (from the repository root): python3 -I proca_timelike_stage2.py <outdir>
Synthetic-only smoke test (does not load registered backgrounds): python3 smoke_stage2_fix.py ./proca_timelike_stage2.py
"""
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import mpmath as mp
import sympy as sp

DATA_B3 = "data/B3_max_speed_7_315096.json"
DATA_MAN = "data/e3_manifest.json"
EXPECTED_SHA256 = {DATA_B3: "99a6b04bc34821ddd5c5d58a78b3474d369aa9fc0ee264163096398b0028ccbd",
                   DATA_MAN: "1bae892a08c3de1f928e4057ee0ad801294a90462478aaedf5b9c724e23d87e0"}
COMMON_M_CSV = "t3fund_proca_b412_t7_drafts/results_b412_common_m/b412_common_m.csv"   # registered result, read for control only
SETS = {"A_g1_rho0.3": (1.0, 0.3), "B_g1_rho0.1": (1.0, 0.1), "C_gStrong_rho0.3": (math.sqrt(4*math.pi), 0.3)}
BORDER_REL_TOL = 1e-12          # registered common-M rule, unchanged
DPS = 50
SIGNS = (+1, -1)

c = sp.symbols('c')

# ---- control kernel: exact 4x4 Lagrangian kernel (Stage 1 construction, rebuilt here independently of the N formula) ----
_w, _q, _v, _r, _gp, _gs, _gu, _om, _k = sp.symbols('w q v r g_phi g_psi gamma_u omega k', real=True)
_J = sp.symbols('f f_t f_x g g_t g_x A0 A0_t A0_x A1 A1_t A1_x', real=True)


def _control_kernels():
    f_, ft, fx, g_, gt, gx, a0, a0t, a0x, a1, a1t, a1x = _J
    F = a1t - a0x
    Lag = (ft**2/2 - fx**2/2 + gt**2/2 - gx**2/2 + F**2/2 + (a0**2 - a1**2)/2
           + a1*(_gp*ft*fx + _gs*gt*gx) + a0*_gu*fx**2)                       # M = 1
    bg = {f_: 0, ft: _w, fx: _q, g_: 0, gt: _v, gx: _r, a0: -_gu*_q**2, a0t: 0, a0x: 0, a1: _gp*_w*_q + _gs*_v*_r,
          a1t: 0, a1x: 0}
    H = sp.hessian(Lag, _J).subs(bg)
    T = sp.zeros(12, 4)
    for i in range(4):
        T[3*i, i] = 1
        T[3*i + 1, i] = -sp.I*_om
        T[3*i + 2, i] = sp.I*_k
    K = (T.H*H*T).applyfunc(sp.expand)
    KLL, KLA, KAL, KAA = K[:2, :2], K[:2, 2:], K[2:, :2], K[2:, 2:]
    Kex = KLL - KLA*KAA.inv()*KAL
    Kld = KLL - KLA*KAA.subs({_om: 0, _k: 0}).inv()*KAL
    args = [_w, _q, _v, _r, _gp, _gs, _gu, _om, _k]
    return sp.lambdify(args, Kex, 'mpmath'), sp.lambdify(args, Kld, 'mpmath')


KEX, KLD = _control_kernels()
X_CONTROL = mp.mpf('0.5')        # one fixed test point x = k^2/M^2 for the control only (not a criterion)


def spec_norm(Mx):
    Mx = (Mx + Mx.H)/2
    return max(abs(e_) for e_ in mp.eighe(Mx, eigvals_only=True))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load():
    for pth, h in EXPECTED_SHA256.items():
        if sha(pth) != h:
            raise SystemExit("ABORT: sha256 mismatch for %s" % pth)
    b3 = json.load(open(DATA_B3, encoding='utf-8'))
    man = json.load(open(DATA_MAN, encoding='utf-8'))['selected']
    rows = [("ANCHOR_B3", b3)] + [("%s_%s_%s" % (r['mode'], r['cls'], r['index']), r) for r in man]
    if len(rows) != 41:
        raise SystemExit("ABORT: expected 41 backgrounds")
    out = []
    for tag, d in rows:
        rec = {key: sp.Rational(str(d[key])) for key in ('p', 's', 'kappa', 'Lambda', 'w', 'q', 'v', 'r')}
        rec['raw'] = {key: d[key] for key in ('p', 's', 'kappa', 'Lambda', 'w', 'q', 'v', 'r')}
        rec['cls'] = d.get('cls', 'anchor')
        rec['mode'] = d.get('mode')
        if rec['p'] <= 0:
            raise SystemExit("ABORT: p <= 0 for %s" % tag)
        out.append((tag, rec))
    return out


def lead_kernel_polys(rec):
    """Khat(c) for the B4 kernel with (p, kappa^2/p, kappa, Lambda), exact rational: returns (F_eff, tr Khat)."""
    p, kap, L, w, q, v, r = (rec[k_] for k_ in ('p', 'kappa', 'Lambda', 'w', 'q', 'v', 'r'))
    s = kap**2/p
    om, k = sp.symbols('omega k', real=True)
    J = sp.symbols('f f_t f_x g g_t g_x', real=True)
    Lag = (J[1]**2/2 - J[2]**2/2 - L/4*J[2]**4 + J[4]**2/2 - J[5]**2/2 + p/2*J[1]**2*J[2]**2
           - kap*J[1]*J[2]*J[4]*J[5] + s/2*J[4]**2*J[5]**2)
    H = sp.hessian(Lag, J).subs({J[0]: 0, J[1]: w, J[2]: q, J[3]: 0, J[4]: v, J[5]: r})
    T = sp.zeros(6, 2)
    for i in range(2):
        T[3*i, i] = 1
        T[3*i + 1, i] = -sp.I*om
        T[3*i + 2, i] = sp.I*k
    K = (T.H*H*T).applyfunc(sp.expand)
    Kh = (K.subs(om, c*k)/k**2).applyfunc(sp.expand)
    F = sp.Poly(sp.expand(Kh.det()), c)
    tr = sp.Poly(sp.expand(Kh.trace()), c)
    if not all(co.is_rational for co in F.all_coeffs() + tr.all_coeffs()):
        raise SystemExit("ABORT: non-rational lead-kernel coefficients")
    return F, tr


def mpq(xr):
    xr = sp.Rational(xr)
    return mp.mpf(int(xr.p))/mp.mpf(int(xr.q))


def exact_branches(F, common):
    """Distinct real roots of F_eff with multiplicity, rational isolating interval, and EXACT R-UNDEFINED flag.
    R-UNDEFINED <=> the root is a root of g = gcd(F_eff, tr Khat). Since g | F_eff, inside an isolating interval of
    F_eff g can only vanish at that root; this is decided with exact rational arithmetic (Sturm count on the open
    interval, exact evaluation for a degenerate interval [r, r]). No tolerance is used."""
    Fs = sp.Poly(sp.quo(F, sp.gcd(F, F.diff(c))), c)        # square-free part, same distinct roots
    ivs = sorted(F.intervals(), key=lambda t_: (t_[0][0], t_[0][1]))
    gdeg = sp.degree(common, c)

    def open_count(P, a, b):                                 # number of distinct roots of P in the open interval (a, b)
        n_ = P.count_roots(a, b)
        return n_ - (1 if P.eval(a) == 0 else 0) - (1 if P.eval(b) == 0 else 0)

    out_ = []
    for (a, b), mult in ivs:
        a, b = sp.Rational(a), sp.Rational(b)
        if a == b:
            if F.eval(a) != 0:
                raise SystemExit("ABORT: degenerate isolating interval is not a root")
            und_ = bool(gdeg > 0 and common.eval(a) == 0)
            with mp.workdps(DPS + 20):
                z_ = mpq(a)
        else:
            if open_count(Fs, a, b) != 1:
                raise SystemExit("ABORT: isolating interval does not contain exactly one root in its interior")
            und_ = bool(gdeg > 0 and open_count(common, a, b) == 1)
            lo, hi = F.refine_root(a, b, eps=sp.Rational(1, 10**(DPS + 15)))
            lo, hi = sp.Rational(lo), sp.Rational(hi)
            with mp.workdps(DPS + 20):
                z_ = (mpq(lo) + mpq(hi))/2
        out_.append((z_, int(mult), (a, b), und_))
    if len(out_) != Fs.count_roots():
        raise SystemExit("ABORT: number of isolating intervals differs from the number of distinct real roots")
    if gdeg > 0 and sum(1 for *_, u_ in out_ if u_) != sp.Poly(common, c).count_roots():
        raise SystemExit("ABORT: R-UNDEFINED count differs from the number of distinct real roots of gcd(F_eff, tr Khat)")
    return out_


def evaluate(tag, rec):
    F, tr = lead_kernel_polys(rec)
    res = sp.resultant(F.as_expr(), tr.as_expr(), c)
    common = sp.gcd(F, tr)
    sqf = sp.degree(sp.gcd(F, F.diff(c)), c) == 0
    common = sp.Poly(common, c)
    branch_list = exact_branches(F, common)
    nreal = sum(m_ for _, m_, _, _ in branch_list)          # with multiplicity
    roots = [z_ for z_, _, _, _ in branch_list]
    p, kap, L, w, q, v, r = (rec[k_] for k_ in ('p', 'kappa', 'Lambda', 'w', 'q', 'v', 'r'))
    out = dict(tag=tag, cls=rec['cls'], mode=rec['mode'], n_real_branches=nreal, n_distinct_real_branches=len(branch_list),
               F_squarefree=bool(sqf), gcd_F_tr_degree=int(sp.degree(common, c)),
               resultant_F_tr_is_zero=bool(res == 0), s_minus_kappa2_over_p=str(sp.N(rec['s'] - kap**2/p, 20)))
    # ---- D4 (sign independent) ----
    pf, kf, Lf = float(p), float(kap), float(L)
    G = max(abs(float(q)), abs(float(w)), abs(float(r)), abs(float(v)))
    cmax = max(math.sqrt(pf), abs(kf)/math.sqrt(pf))
    cL = math.sqrt(Lf/2)
    cmax2 = max(cmax, cL)
    d4 = {}
    for sname, (gmax, rmax) in SETS.items():
        scale = gmax*rmax
        def cls_of(m):
            return "PASS" if m > BORDER_REL_TOL*scale else ("FAIL" if m < -BORDER_REL_TOL*scale else "BORDERLINE")
        m0, m1 = scale - cmax*G, scale - cmax2*G
        d4[sname] = dict(T7a_margin=m0, T7a_class=cls_of(m0), extended_margin=m1, extended_class=cls_of(m1))
    out['D4'] = dict(G=G, c_max=cmax, c_Lambda=cL, c_max_extended=cmax2, c_Lambda_dominates=bool(cL > cmax), sets=d4)
    # ---- D3 per sign ----
    out['D3'] = {}
    with mp.workdps(DPS):
        P, KAP, LL, W, Q, V, R = (mpq(z_) for z_ in (p, kap, L, w, q, v, r))
        gphi = mp.sqrt(P)
        gpsi = -KAP/gphi
        Fc = [mpq(co) for co in F.all_coeffs()]
        trc = [mpq(co) for co in tr.all_coeffs()]
        for sg in SIGNS:
            gu = sg*mp.sqrt(LL/2)
            branches = []
            for ci, mult, (ia, ib), rundef in branch_list:
                ci = +ci
                trv = mp.polyval(trc, ci)
                b = dict(c=mp.nstr(ci, 20), multiplicity=int(mult), isolating_interval=[str(ia), str(ib)],
                         tr_Khat=mp.nstr(trv, 20), R_UNDEFINED=bool(rundef))
                if rundef:
                    b['R_UNDEFINED_reason'] = ("common real root of F_eff and tr Khat: gcd(F_eff, tr Khat) has a root in the "
                                               "rational isolating interval of this root of F_eff (exact Sturm count / exact "
                                               "evaluation, no tolerance)")
                ux = 2*gu*Q - ci*gphi*(W - ci*Q)
                ug = -ci*gpsi*(V - ci*R)
                h3 = gu*Q*(gu*Q - ci*gphi*(W - ci*Q))
                b['H3_holds'] = bool(h3 <= 0)
                b['H3_value'] = mp.nstr(h3, 15)
                if abs(ci) > 1:
                    b['pole_k_over_M'] = mp.nstr(1/mp.sqrt(ci**2 - 1), 20)
                if not rundef:
                    N = (ux**2 + ug**2)/abs(trv)
                    # control: direct ratio of spectral norms from the exact kernels at x = X_CONTROL vs the formula
                    kk = mp.sqrt(X_CONTROL)
                    argsv = [W, Q, V, R, gphi, gpsi, gu, ci*kk, kk]
                    Kx, Kl = mp.matrix(KEX(*argsv)), mp.matrix(KLD(*argsv))
                    direct = spec_norm(Kx - Kl)/spec_norm(Kl)
                    formula = N*X_CONTROL/abs(1 + (1 - ci**2)*X_CONTROL)
                    b['control_R_direct_vs_formula_rel'] = mp.nstr(abs(direct - formula)/max(abs(formula), mp.mpf(10)**-40), 5)
                    S = N + ci**2
                    b['N'] = mp.nstr(N, 20)
                    b['N_plus_c2'] = mp.nstr(S, 20)
                    if S > 1:
                        b['crossing'] = True
                        b['k_der_over_M'] = mp.nstr(1/mp.sqrt(S - 1), 20)
                    else:
                        b['crossing'] = False
                    b['V2_branch_margin_2_minus_N_plus_c2'] = mp.nstr(2 - S, 20)
                branches.append(b)
            und_any = any(b_['R_UNDEFINED'] for b_ in branches)
            crossing = [b_ for b_ in branches if b_.get('crossing')]
            kd = min((mp.mpf(b_['k_der_over_M']) for b_ in crossing), default=None)
            # V2: None (undetermined) if there is no real branch (no vacuous all([])) or any branch is R-UNDEFINED
            if not branches or und_any:
                v2 = None
            else:
                v2 = all(mp.mpf(b_['V2_branch_margin_2_minus_N_plus_c2']) > 0 for b_ in branches)
            if not branches:
                kstat = 'undetermined (no real branches)'
            elif und_any:
                kstat = 'undetermined (R-UNDEFINED)'
            elif crossing:
                kstat = 'defined'
            else:
                kstat = 'no crossing on all real branches'
            out['D3'][str(sg)] = dict(
                branches=branches,
                n_crossing=len(crossing), n_no_crossing=sum(1 for b_ in branches if b_.get('crossing') is False),
                n_R_UNDEFINED=sum(1 for b_ in branches if b_['R_UNDEFINED']),
                k_der_status=kstat,
                k_der_over_M=(mp.nstr(kd, 20) if kstat == 'defined' else None),
                k_der_over_M_min_over_defined_branches_only=(mp.nstr(kd, 20) if (kstat != 'defined' and kd is not None) else None),
                V2_holds=v2,
                H3_all_branches=(all(b_['H3_holds'] for b_ in branches) if branches else None))
    return out


def main(outdir):
    T0 = time.time()
    od = Path(outdir)
    od.mkdir(parents=True, exist_ok=True)
    cases = load()
    rows = [evaluate(tag, rec) for tag, rec in cases]
    # control: T7a-side counts reproduce the registered common-M result (CSV, read only here)
    import csv
    reg = {r_['tag']: r_ for r_ in csv.DictReader(open(COMMON_M_CSV, encoding='utf-8'))}
    ctrl = []
    for row in rows:
        rr = reg.get(row['tag'])
        if rr is None:
            ctrl.append((row['tag'], 'missing in registered CSV'))
            continue
        for key in rr:
            if key.lower().startswith('c_max'):
                if abs(float(rr[key]) - row['D4']['c_max']) > 1e-12*max(1, row['D4']['c_max']):
                    ctrl.append((row['tag'], 'c_max differs'))
    summ = {}
    for sname in SETS:
        summ[sname] = dict(
            T7a={k_: sum(1 for r_ in rows if r_['D4']['sets'][sname]['T7a_class'] == k_) for k_ in ('PASS', 'FAIL', 'BORDERLINE')},
            extended={k_: sum(1 for r_ in rows if r_['D4']['sets'][sname]['extended_class'] == k_) for k_ in ('PASS', 'FAIL', 'BORDERLINE')})
    worst = max((mp.mpf(b_['control_R_direct_vs_formula_rel']) for r_ in rows for sg in SIGNS
                 for b_ in r_['D3'][str(sg)]['branches'] if 'control_R_direct_vs_formula_rel' in b_), default=mp.mpf(0))
    d3s = {}
    for sg in SIGNS:
        key = str(sg)
        d3s[key] = dict(
            V2_true=sum(1 for r_ in rows if r_['D3'][key]['V2_holds'] is True),
            V2_false=sum(1 for r_ in rows if r_['D3'][key]['V2_holds'] is False),
            V2_undetermined=sum(1 for r_ in rows if r_['D3'][key]['V2_holds'] is None),
            backgrounds_with_R_UNDEFINED=sum(1 for r_ in rows if r_['D3'][key]['n_R_UNDEFINED'] > 0),
            backgrounds_without_any_crossing=sum(1 for r_ in rows if r_['D3'][key]['n_crossing'] == 0),
            H3_all_branches_true=sum(1 for r_ in rows if r_['D3'][key]['H3_all_branches']))
    result = dict(
        spec="docs/PROCA_TIMELIKE_CURRENT_ANALYTIC_SPEC_v0.5.md (87fc31b)", stage="2 (arithmetic only)",
        not_included="roots of Q_new(c;1), V3, NOT-VIABLE (d): not authorised",
        inputs={k_: sha(k_) for k_ in EXPECTED_SHA256}, script_sha256=sha(__file__),
        environment=dict(python=platform.python_version(), sympy=sp.__version__, mpmath=mp.__version__),
        conventions="M = 1; g_phi = +sqrt(p); g_psi = -kappa/sqrt(p); gamma_u = sigma sqrt(Lambda/2), sigma = +1, -1",
        control_c_max_vs_registered_common_M=("PASS" if not ctrl else ctrl),
        control_R_direct_vs_formula_max_rel=mp.nstr(worst, 5),
        summary_D4=summ, summary_D3=d3s,
        class_note=("Class after Stage 2 cannot be VIABLE-FOR-NUMERICS without V3 (not evaluated). If V2 fails on any "
                    "background for a sign, that sign is OPEN regardless of V3."),
        backgrounds=rows, seconds=round(time.time() - T0, 1))
    txt = json.dumps(result, indent=1, default=str) + "\n"
    (od / "proca_timelike_stage2.json").write_text(txt, encoding='utf-8', newline='\n')
    print(json.dumps(dict(summary_D4=summ, summary_D3=d3s, control=result['control_c_max_vs_registered_common_M']), indent=1))
    return 0 if not ctrl else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(64)
    sys.exit(main(sys.argv[1]))
