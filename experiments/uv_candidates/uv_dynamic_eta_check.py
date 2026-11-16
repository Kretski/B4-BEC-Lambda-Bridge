"""uv_dynamic_eta_check.py -- script A of the dynamic-eta UV candidate (revision 1.1, FOR STATIC REVIEW; NOT RUN).

Revision 1.1 (review of 2026-11-16): P2a on the full determinant (multiplicities kept); T6 needs the full-system
confirmation; C1 analytic only with a block-diagonal principal symbol; fixed input hashes; singular F_eta_eta handled
locally; DC5 union matching; exact root multiplicities (square-free decomposition), unresolved roots INCONCLUSIVE.

Implements docs/UV_DYNAMIC_ETA_SCRIPT_SPEC_v0.3.2.md (commit aece78a, sha256 8bbd57f8...) for the candidate
docs/UV_DYNAMIC_ETA_SPEC_v0.4.1.md (commit 5fdd4f6, sha256 1379e599...) and the comparison document
docs/UV_DYNAMIC_ETA_VS_DYNAMIC_U_COMPARISON_v0.3.1.md (commit 5fdd4f6, sha256 024aaae1...).
Registered protocol: docs/UV_ACCEPTANCE_PROTOCOL.md, docs/UV_ACCEPTANCE_PROTOCOL_ADDENDUM_T3_T6.md (220bf76, d70c1e5).

Evidence level: every output of this script is level 1 (theoretical consistency), the internal checks DC0-DC6 level 2
(numerical consistency). Nothing here is experimental validation of B4.  ALL_DC_PASS covers DC0-DC6 only; it is not a
physics verdict and not an independent proof.  C4 (time domain, script B) is NOT EVALUATED here.

Usage:  python3 -I uv_dynamic_eta_check.py <outdir>
        (run from the repository root; inputs data/B3_max_speed_7_315096.json and data/e3_manifest.json, key 'selected')

------------------------------------------------------------------------------------------------------------------------
FIXED BEFORE ANY RUN (all values from the registered specifications; none is changed after results are seen)
  signature (+,-); 1+1D; M = 1; f in {0.1, 1, 10}; c14 = c123 = f (Einstein-aether variant implemented, NOT run, S8)
  couplings g1 = sqrt(p), g2 = -kappa/sqrt(p) (rank one); stop if p <= 0
  action B1:  L = 1/2(dphi)^2 + 1/2(dpsi)^2 + 1/2(dsigma)^2 - 1/2 M^2 sigma^2 + 1/2 f (deta)^2
                  + sigma sum_i g_i (u.dphi_i)(n.dphi_i) - (Lam/4)(n.dphi)^4,
              u.dF = cosh(eta) F_t + sinh(eta) F_x,  n.dF = sinh(eta) F_t + cosh(eta) F_x  (as in U1-U4)
  S1.2  exact enumeration of eta_bar via z = e^eta, rationalised inputs sympy.Rational(str(x)), Sturm count on (0, inf),
        roots located at 130 digits; double roots / relative separation < 1e-20 flagged
  S1.3  C0 residual <= 1e-10 at 50 digits (sigma and eta equations from the action; phi, psi identically zero)
  S2.3  route R2: mpmath 100 digits, central differences (3-point pure, 4-point mixed), h1 = 1e-20, h2 = 1e-25 relative to
        max(1, |jet|); comparison |a-b| <= atol + rtol max(|a|,|b|), rtol = 1e-30, atol = 1e-30 S, S = max |K_R1| at the
        point; ZERO (agree) needs both <= atol AND the comparison rule; points (omega,k) = (0.3,0.1), (1.1,0.7), (2.0,3.0),
        (0.5,10.0), (7.0,2.0)
  S3    FUNDAMENTAL regime, evaluated in the frame of u (preferred time; see FRAME below); P1 symbolic factorisation of
        det K; P2 grid 61 log points k/M in [1e-3,1e3] + 31 in (1e3,1e6]; P3 bounds 1e-7 / 1e-5, confirmation at 50 digits
        r > 1e-12; P2a at k = 0 with r = |Im omega|/max(M,|omega|); P4 fit on [1e4,1e6]; P5 cond threshold 1e8, coincidence
        1e-6 relative; P6 u(k) = |gamma50 - gamma100|, margin 10, refinement factor 4 (3 extra log points per interval of the
        band and one neighbouring interval on each side), FAIL-UNSTABLE needs gamma > 10 u on >= 3 consecutive refined points
  S4    T5 tolerance -1e-12; k in {0.01, 0.1, 1, 10} for H2 and Krein; u frame and lab frame (lab only where eta_bar != 0)
  S5    |eta_bar| > 10 flagged as unbounded-rapidity edge
  S6    classes TACHYONIC (m2 < 0), LIGHT (0 <= m2 < 1e-2 M^2), HEAVY (m2 >= 1e-2 M^2); BOUNDARY within 1e-40 of 0 or 1e-2
  DC5   J' = 0 roots identical to 1e-30
------------------------------------------------------------------------------------------------------------------------
IMPLEMENTATION CHOICES (stated here so that the reviewer can check them before the run)
  FRAME   C1-C3 (S3) are evaluated in the frame of u, i.e. with respect to the preferred time, as T3-fund requires.  The
          u-frame system is the action with eta -> eta - eta_bar and boosted derivatives (exact: the action is invariant
          under a boost of the coordinates together with eta -> eta + beta); background gradients become
          (w', q') = (u.dphi, n.dphi) and (v', r') = (u.dpsi, n.dpsi).  Check X1 confirms det K_lab(omega,k) =
          det K_u(omega',k') with omega' = cosh(eb) omega - sinh(eb) k, k' = cosh(eb) k - sinh(eb) omega.  DC1 is done in
          the lab frame (it tests the eta_bar dependence of the derivation).
  C1      PASS [ANALYTIC] when the principal symbol is real with a complete eigenvector set (strong hyperbolicity: then
          Im omega(k) stays bounded as k -> infinity, the lower-order terms being bounded perturbations); otherwise the
          numerical result is "numerically passed on [1e-3,1e6]" = INCONCLUSIVE, or FAIL if the principal roots are complex.
  C2      first-order form U = (F, F_t, F_x) per field; principal matrix block-diagonal (derived: every off-diagonal entry of
          K has no degree-2 part, checked symbolically); per field the 2x2 companion of A c^2 + B c + C; the F component has
          no principal part (eigenvalue 0, own eigenvector).  Roots coinciding between DIFFERENT blocks: complete by the
          direct-sum argument (block-diagonal matrix with each block diagonalisable).  A double root inside one block (or a
          vanishing time coefficient A): INCONCLUSIVE.
  P2a     ZERO MODE when the multiplicity of omega = 0 in det K(omega, 0) is exactly 4 (phi and psi massless, one double
          root each) and every other root is real.  Any other multiplicity of omega = 0: INCONCLUSIVE.
  P6      the 'largest change between a coarse point and its refined neighbours' part of u(k) is implemented, as in the
          approved T3-fund stage 1 script, only through the sign-change / three-consecutive-points rule (the registered text
          is ambiguous; reported in the JSON).
  T5      E2(k) = blockdiag(A, C(k)) with A = velocity block, C(k) = potential block (the velocity-coordinate cross terms
          drop out of the energy).  T5-static by route 2 of Part D: grid over P2 at 50 digits + analytic large-k statement
          (C(k)/k^2 -> C2) + analytic k -> 0 limit (Schur complement of the light block).  PASS needs A > 0, C2 > 0, the
          k -> 0 Schur complement >= -1e-12 and every grid point >= -1e-12 max(1, |C(k)|); FAIL if any of these is violated
          beyond the tolerance; INCONCLUSIVE if A or C2 is singular within the tolerance.
  T6      per (background, root, f): STATED-DYNAMICAL if eta_bar exists, m2 >= 0 (not BOUNDARY), confirmed in the full
          system (P2a passes) and C1-C3 pass; FAIL if C1-C3 fail; INCONCLUSIVE otherwise.  The eta sector cannot be
          separated from the full coupled system, so C1-C3 of the whole system are used.  Set level: PARTIAL when eta_bar
          exists on a strict subset (subset and complement listed).
  T7      per (background, root, f), after F2: LIGHT -> NOT-REPRODUCED (F2(ii)); TACHYONIC -> NOT-REPRODUCED (F2(ii'),
          elimination invalid at a maximum of V_eff); BOUNDARY -> INCONCLUSIVE; HEAVY -> eta eliminated at tree level by its
          own (algebraic at leading order) equation; REPRODUCED-RANK-ONE only if, in addition to the exact frozen-eta
          operator identity, the eta-exchange correction to the quadratic operator vanishes (relative <= 1e-30), otherwise
          NOT-REPRODUCED (operator of the same order modified, coefficient functions of the background).  Set level:
          CONDITIONAL when classes are mixed.  The frozen-eta comparison is reported separately as a diagnostic.
  extra checks X1 (frame consistency) and X2 (low-k roots of the full system = roots of the eta- and sigma-eliminated
          effective operator) are reported but are NOT part of ALL_DC_PASS.
"""
import hashlib
import json
import math
import platform
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

REVISION = "1.1"
SPEC = "docs/UV_DYNAMIC_ETA_SCRIPT_SPEC_v0.3.2.md"
SPEC_COMMIT = "aece78a"
CANDIDATE_SPEC = "docs/UV_DYNAMIC_ETA_SPEC_v0.4.1.md (5fdd4f6)"
DATA_B3 = "data/B3_max_speed_7_315096.json"
DATA_MAN = "data/e3_manifest.json"
EXPECTED_SHA256 = {DATA_B3: "99a6b04bc34821ddd5c5d58a78b3474d369aa9fc0ee264163096398b0028ccbd",
                   DATA_MAN: "1bae892a08c3de1f928e4057ee0ad801294a90462478aaedf5b9c724e23d87e0"}   # repo fa19637

F_VALUES = (sp.Rational(1, 10), sp.Integer(1), sp.Integer(10))
M_VAL = sp.Integer(1)
TOL_REAL, TOL_UNRES, TOL_CONF = 1e-7, 1e-5, 1e-12
KGRID = np.concatenate([np.geomspace(1e-3, 1e3, 61), np.geomspace(1e3, 1e6, 32)[1:]])
R2_POINTS = ((0.3, 0.1), (1.1, 0.7), (2.0, 3.0), (0.5, 10.0), (7.0, 2.0))
R2_STEPS = ("1e-20", "1e-25")
R2_DPS = 100
RTOL, ATOL_FAC = "1e-30", "1e-30"
C0_TOL = 1e-10
SEP_FLAG = "1e-20"
ROOT_DPS = 130
ETA_EDGE = 10
COND_THR, COINC = 1e8, 1e-6
KS_H = (0.01, 0.1, 1.0, 10.0)
T5_TOL = 1e-12
LIGHT_THR = sp.Rational(1, 100)
BOUNDARY = "1e-40"
DC5_TOL = "1e-30"
ETA_EXCH_ZERO = "1e-30"

# ======================================================================================================================
# symbols and the action (route R1)
# ======================================================================================================================
w, q, v, r, g1, g2, M, Lam, f, sb, eb = sp.symbols('w q v r g1 g2 M Lam f sb eb', real=True)
c14, c123 = sp.symbols('c14 c123', real=True)
om, k, c = sp.symbols('omega k c', real=True)
JETS = sp.symbols('P Pt Px Q Qt Qx S St Sx E Et Ex', real=True)
P_, Pt, Px, Q_, Qt, Qx, S_, St, Sx, E_, Et, Ex = JETS
FIELDS = ('phi', 'psi', 'sigma', 'eta')


def lagrangian(kin):
    """Lagrangian density in jet variables.  kin = 'eta' (B1) or 'aether' (c14, c123 kinetic term of the comparison doc)."""
    ch, sh = sp.cosh(E_), sp.sinh(E_)
    uP, nP = ch*Pt + sh*Px, sh*Pt + ch*Px
    uQ, nQ = ch*Qt + sh*Qx, sh*Qt + ch*Qx
    L = (sp.Rational(1, 2)*(Pt**2 - Px**2) + sp.Rational(1, 2)*(Qt**2 - Qx**2) + sp.Rational(1, 2)*(St**2 - Sx**2)
         - sp.Rational(1, 2)*M**2*S_**2 + S_*(g1*uP*nP + g2*uQ*nQ) - Lam/4*nP**4)
    if kin == 'eta':
        return L + f/2*(Et**2 - Ex**2)
    Du, Dn = ch*Et + sh*Ex, sh*Et + ch*Ex
    return L + c14/2*Du**2 - c123/2*Dn**2


BG = {P_: 0, Pt: w, Px: q, Q_: 0, Qt: v, Qx: r, S_: sb, St: 0, Sx: 0, E_: eb, Et: 0, Ex: 0}


def tmatrix(omega, kk, I=sp.I):
    """12x4 map from plane-wave amplitudes to jets (value, t, x): exp(i(k x - omega t))."""
    T = sp.zeros(12, 4)
    for i in range(4):
        T[3*i, i] = 1
        T[3*i + 1, i] = -I*omega
        T[3*i + 2, i] = I*kk
    return T


def build_symbolic():
    t0 = time.time()
    Leta, Lae = lagrangian('eta'), lagrangian('aether')
    H = sp.hessian(Leta, JETS).subs(BG)
    Hae = sp.hessian(Lae, JETS).subs(BG)
    grad = [sp.diff(Leta, j).subs(BG) for j in JETS]
    K = (tmatrix(om, k).H*H*tmatrix(om, k)).applyfunc(sp.expand)
    Kae = (tmatrix(om, k).H*Hae*tmatrix(om, k)).applyfunc(sp.expand)
    Ku = K.subs(eb, 0).applyfunc(sp.expand)
    det_u = sp.expand(Ku.det(method='berkowitz'))
    return dict(Leta=Leta, Lae=Lae, H=H, Hae=Hae, grad=grad, K=K, Kae=Kae, Ku=Ku, det_u=det_u, seconds=time.time() - t0)


# ======================================================================================================================
# S0  symbolic checks of the comparison document
# ======================================================================================================================
def canon_pow(e):
    """Rewrite in exp and canonicalise the base of every non-integer power (e.g. sqrt(cosh^2 - sinh^2) -> 1)."""
    e = sp.sympify(e).rewrite(sp.exp)
    return e.replace(lambda x_: x_.is_Pow and not x_.exp.is_integer,
                     lambda x_: sp.cancel(sp.expand(x_.base))**x_.exp)


def zero_expr(e):
    """Symbolic zero test: rewrite in exp, together, expand the numerator."""
    e = sp.together(sp.expand(sp.sympify(e).rewrite(sp.exp)))
    return sp.expand(sp.numer(e)) == 0


def s0_checks():
    t, x = sp.symbols('t x', real=True)
    c1, c2, c3, c4, lam_ = sp.symbols('c1 c2 c3 c4 lam', real=True)
    X = (t, x)
    g = sp.diag(1, -1)
    eta = sp.Function('eta')(t, x)
    u = [sp.cosh(eta), sp.sinh(eta)]
    n = [sp.sinh(eta), sp.cosh(eta)]
    du = [[sp.diff(u[mu], X[al]) for mu in range(2)] for al in range(2)]          # du[al][mu] = d_al u^mu

    def terms(uu, d):
        T1 = sum(g[al, al]*g[mu, mu]*d[al][mu]**2 for al in range(2) for mu in range(2))
        T2 = (d[0][0] + d[1][1])**2
        T3 = sum(d[al][mu]*d[mu][al] for al in range(2) for mu in range(2))
        T4 = sum(uu[al]*uu[be]*g[mu, mu]*d[al][mu]*d[be][mu] for al in range(2) for be in range(2) for mu in range(2))
        return T1, T2, T3, T4

    T1, T2, T3, T4 = terms(u, du)
    Lkin = -sp.Rational(1, 2)*(c1*T1 + c2*T2 + c3*T3 + c4*T4)
    Du = u[0]*sp.diff(eta, t) + u[1]*sp.diff(eta, x)
    Dn = n[0]*sp.diff(eta, t) + n[1]*sp.diff(eta, x)
    red = sp.Rational(1, 2)*(c1 + c4)*Du**2 - sp.Rational(1, 2)*(c1 + c2 + c3)*Dn**2
    out = {}
    out['reduction_to_c14_c123'] = zero_expr(Lkin - red)
    out['c2_equals_c3_pointwise_unit_field'] = zero_expr(T2 - T3)
    # general (non-unit) field: c3-term - c2-term = d_al(u^mu d_mu u^al - u^al d_mu u^mu)
    U0, U1 = sp.Function('U0')(t, x), sp.Function('U1')(t, x)
    U = [U0, U1]
    dU = [[sp.diff(U[mu], X[al]) for mu in range(2)] for al in range(2)]
    G1, G2, G3, G4 = terms(U, dU)
    divv = sum(sp.diff(sum(U[mu]*sp.diff(U[al], X[mu]) for mu in range(2))
                       - U[al]*sum(sp.diff(U[mu], X[mu]) for mu in range(2)), X[al]) for al in range(2))
    out['c3_minus_c2_total_derivative_general_field'] = sp.expand(G3 - G2 - divv) == 0
    out['c14_eq_c123_eq_f_gives_half_f_deta2'] = zero_expr(
        red.subs({c1: f, c2: 0, c3: 0, c4: 0}) - f/2*(sp.diff(eta, t)**2 - sp.diff(eta, x)**2))
    dc3 = all(out.values())

    # ---------------- S0.2  V1 and V2 (jet variables; total derivatives built explicitly) ----------------
    def jet(name):
        return dict(zip(('0', 't', 'x', 'tt', 'tx', 'xx'), sp.symbols(' '.join(name + s_ for s_ in ('', '_t', '_x', '_tt', '_tx', '_xx')), real=True)))
    jf = {nm_: jet(nm_) for nm_ in ('phi', 'psi', 'sig', 'eta', 'U0', 'U1')}
    nxt = {('0', 't'): 't', ('0', 'x'): 'x', ('t', 't'): 'tt', ('t', 'x'): 'tx', ('x', 't'): 'tx', ('x', 'x'): 'xx'}

    def D(g_, a_):
        """total derivative d/da of an expression in field values and first derivatives."""
        out_ = 0
        for jj in jf.values():
            for o_ in ('0', 't', 'x'):
                if g_.has(jj[o_]):
                    out_ += sp.diff(g_, jj[o_])*jj[nxt[(o_, a_)]]
        return out_

    def euler(L, nm_):
        jj = jf[nm_]
        return sp.diff(L, jj['0']) - D(sp.diff(L, jj['t']), 't') - D(sp.diff(L, jj['x']), 'x')

    def matter(uu, nn):
        def ud(F): return uu[0]*F['t'] + uu[1]*F['x']
        def nd(F): return nn[0]*F['t'] + nn[1]*F['x']
        P_, Q_, S_ = jf['phi'], jf['psi'], jf['sig']
        return (sp.Rational(1, 2)*(P_['t']**2 - P_['x']**2) + sp.Rational(1, 2)*(Q_['t']**2 - Q_['x']**2)
                + sp.Rational(1, 2)*(S_['t']**2 - S_['x']**2) - sp.Rational(1, 2)*M**2*S_['0']**2
                + S_['0']*(g1*ud(P_)*nd(P_) + g2*ud(Q_)*nd(Q_)) - Lam/4*nd(P_)**4)

    def kin(uu, du_):
        """-1/2 K du du with du_[al][mu] = d_al u^mu."""
        T1_, T2_, T3_, T4_ = terms(uu, du_)
        return -sp.Rational(1, 2)*(c1*T1_ + c2*T2_ + c3*T3_ + c4*T4_)

    e_ = jf['eta']
    ue = [sp.cosh(e_['0']), sp.sinh(e_['0'])]
    ne = [sp.sinh(e_['0']), sp.cosh(e_['0'])]
    due = [[ne[mu]*e_[a_] for mu in range(2)] for a_ in ('t', 'x')]
    L1 = matter(ue, ne) + kin(ue, due)                                    # V1: substitution u(eta), n(eta)
    EL1 = euler(L1, 'eta')
    Uj = [jf['U0'], jf['U1']]
    Uv = [Uj[0]['0'], Uj[1]['0']]
    nrm = sp.sqrt(Uv[0]**2 - Uv[1]**2)
    nU = [Uv[1]/nrm, Uv[0]/nrm]                                           # V2: off-shell extension of n(u)
    dU_ = [[Uj[mu][a_] for mu in range(2)] for a_ in ('t', 'x')]
    L2 = matter(Uv, nU) + kin(Uv, dU_) + lam_*(Uv[0]**2 - Uv[1]**2 - 1)
    E0, E1 = euler(L2, 'U0'), euler(L2, 'U1')                             # lower-index E_mu = dL/dU^mu
    ch_, sh_ = sp.cosh(e_['0']), sp.sinh(e_['0'])
    on = {Uj[0]['0']: ch_, Uj[1]['0']: sh_,
          Uj[0]['t']: sh_*e_['t'], Uj[0]['x']: sh_*e_['x'], Uj[1]['t']: ch_*e_['t'], Uj[1]['x']: ch_*e_['x'],
          Uj[0]['tt']: ch_*e_['t']**2 + sh_*e_['tt'], Uj[0]['tx']: ch_*e_['t']*e_['x'] + sh_*e_['tx'], Uj[0]['xx']: ch_*e_['x']**2 + sh_*e_['xx'],
          Uj[1]['tt']: sh_*e_['t']**2 + ch_*e_['tt'], Uj[1]['tx']: sh_*e_['t']*e_['x'] + ch_*e_['tx'], Uj[1]['xx']: sh_*e_['x']**2 + ch_*e_['xx']}
    yv = sp.symbols('y', positive=True)
    toy = {sp.cosh(e_['0']): (yv + 1/yv)/2, sp.sinh(e_['0']): (yv - 1/yv)/2}

    def ycanon(ex_):
        ex_ = ex_.xreplace(on).xreplace(toy)
        return ex_.replace(lambda x_: x_.is_Pow and not x_.exp.is_integer, lambda x_: sp.cancel(x_.base)**x_.exp)
    E0s, E1s = ycanon(E0), ycanon(E1)
    ny, uy = [(yv - 1/yv)/2, (yv + 1/yv)/2], [(yv + 1/yv)/2, (yv - 1/yv)/2]
    tang = ny[0]*E0s + ny[1]*E1s                                          # projection on n^mu
    norm_proj = uy[0]*E0s + uy[1]*E1s                                     # projection on u^mu
    dfr = sp.together(tang - EL1.xreplace(toy))
    dc0 = sp.expand(sp.numer(dfr)) == 0
    lam_sol = sp.solve(sp.Eq(sp.numer(sp.together(norm_proj)), 0), lam_)
    lam_const = None
    if len(lam_sol) == 1:
        consts = {}
        for nm_ in ('phi', 'psi', 'sig', 'eta'):
            for o_ in ('tt', 'tx', 'xx'):
                consts[jf[nm_][o_]] = 0
        for nm_ in ('sig', 'eta'):
            for o_ in ('t', 'x'):
                consts[jf[nm_][o_]] = 0
        lam_const = sp.factor(sp.simplify(lam_sol[0].xreplace(consts)))
    return dict(DC3=dc3, S0_1=out, DC0=dc0,
                S0_2=dict(V1="Euler-Lagrange of L(u(eta), n(eta)) w.r.t. eta, general c1..c4",
                          V2="n^mu(u) = (u^1, u^0)/sqrt(u.u); free variation in u^mu and lam; projections on n and u",
                          tangential_projection_equals_V1=dc0,
                          lambda_from_V2=[str(s_) for s_ in lam_sol],
                          lambda_on_constant_background=str(lam_const),
                          variables="jets: field, first and second derivatives; cosh(eta), sinh(eta) -> (y +- 1/y)/2, y = e^eta",
                          note="lambda is reported from V2 only; V1 does not determine it"))


# ======================================================================================================================
# backgrounds: inputs and exact enumeration of eta_bar (S1)
# ======================================================================================================================
def sha_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def verify_inputs(paths=None):
    """ABORT unless both registered input files have their fixed SHA-256."""
    paths = paths or {DATA_B3: DATA_B3, DATA_MAN: DATA_MAN}
    for key, path in paths.items():
        if not Path(path).is_file():
            raise SystemExit("ABORT: input file missing: %s" % path)
        h = sha_file(path)
        if h != EXPECTED_SHA256[key]:
            raise SystemExit("ABORT: sha256 mismatch for %s: %s (expected %s)" % (path, h, EXPECTED_SHA256[key]))
    return True


def load_backgrounds(paths=None):
    paths = paths or {DATA_B3: DATA_B3, DATA_MAN: DATA_MAN}
    verify_inputs(paths)
    b3 = json.load(open(paths[DATA_B3], encoding='utf-8'))
    man = json.load(open(paths[DATA_MAN], encoding='utf-8'))['selected']
    rows = [('B3_anchor', b3)] + [('manifest_%02d_index_%s' % (i, d.get('index')), d) for i, d in enumerate(man)]
    if len(rows) != 41:
        raise SystemExit("ABORT: expected 41 backgrounds, found %d" % len(rows))
    out = []
    for name, d in rows:
        rec = {}
        for key in ('p', 's', 'kappa', 'Lambda', 'w', 'q', 'v', 'r'):
            val = d[key]
            if isinstance(val, bool) or not isinstance(val, (int, float)) or not math.isfinite(val):
                raise SystemExit("ABORT: invalid %s in %s" % (key, name))
            rec[key] = sp.Rational(str(val))
            rec[key + '_str'] = str(val)
        if rec['p'] <= 0:
            raise SystemExit("ABORT: p <= 0 in %s" % name)
        out.append((name, rec))
    return out


zz = sp.symbols('z', positive=True)


def background_functions(rec, Lv=None):
    """Exact a = sqrt(p) A, b = sqrt(p) B and the polynomial P(z) of the eta equation (S1.2)."""
    p_, kap, w_, q_, v_, r_ = rec['p'], rec['kappa'], rec['w'], rec['q'], rec['v'], rec['r']
    La = rec['Lambda'] if Lv is None else Lv
    a = p_*w_*q_ - kap*v_*r_
    b = sp.Rational(1, 2)*(p_*(w_**2 + q_**2) - kap*(v_**2 + r_**2))
    J2 = ((a + b)*zz**2 + (a - b)/zz**2)**2/(4*p_)              # J^2, rational (1/p carries the 1/sqrt(p)^2)
    N = ((w_ + q_)*zz + (q_ - w_)/zz)/2
    V = J2/(2*M_VAL**2) - La/4*N**4
    expr = sp.together(zz*sp.diff(V, zz))                        # dV/deta = z dV/dz
    num = sp.Poly(sp.expand(sp.numer(expr)), zz)
    nz = 0
    while num.eval(0) == 0 and not num.is_zero:
        num = sp.Poly(sp.cancel(num.as_expr()/zz), zz)
        nz += 1
    return dict(a=a, b=b, V=V, P=num, stripped_z_power=nz)


def roots_with_multiplicity(P):
    """distinct real roots z > 0 of P with exact multiplicities (square-free decomposition), at ROOT_DPS digits."""
    out = []
    for fac, mult in P.sqf_list()[1]:
        for rt in sp.Poly(fac, zz).real_roots():          # fac is square-free: every root listed once
            val = sp.N(rt, ROOT_DPS)
            if val > 0:
                out.append((val, int(mult)))
    return sorted(out, key=lambda t_: t_[0])


def enumerate_roots(rec):
    bf = background_functions(rec)
    P = bf['P']
    sturm = int(P.count_roots(0, None))                    # distinct real roots in [0, inf) (z = 0 stripped)
    sqf = sp.degree(sp.gcd(P, P.diff(zz)), zz) == 0
    located = roots_with_multiplicity(P)
    roots = []
    with mp.workdps(ROOT_DPS):
        zs = [mp.mpf(str(z_)) for z_, _ in located]
        close = [False]*len(zs)
        for i in range(len(zs)):
            for j in range(i + 1, len(zs)):
                if abs(zs[i] - zs[j]) <= mp.mpf(SEP_FLAG)*max(abs(zs[i]), abs(zs[j])):
                    close[i] = close[j] = True
        for (z_, mult), zm, cl in zip(located, zs, close):
            roots.append(dict(z=mp.nstr(zm, ROOT_DPS - 5), eta_bar=mp.nstr(mp.log(zm), ROOT_DPS - 5),
                              multiplicity=mult, simple=bool(mult == 1), close_to_other_root=bool(cl),
                              resolved=bool(mult == 1 and not cl)))
    return dict(P_degree=int(P.degree()), stripped_z_power=int(bf['stripped_z_power']), sturm_count=sturm,
                n_located=len(roots), squarefree=bool(sqf), separation_flag=bool(any(r_['close_to_other_root'] for r_ in roots)),
                multiplicity_method="exact square-free decomposition of P (sympy sqf_list)", roots=roots,
                a=str(bf['a']), b=str(bf['b'])), bf


def mpq(x):
    """exact sympy Rational -> mpf at the current precision."""
    x = sp.Rational(x)
    return mp.mpf(int(x.p))/mp.mpf(int(x.q))


def numeric_bg(rec, eta_str, fval, frame):
    """Background values at the current mp precision.  frame 'lab': (w,q,v,r,eb); frame 'u': boosted gradients, eb = 0."""
    p_ = mpq(rec['p'])
    G1 = mp.sqrt(p_)
    G2 = -mpq(rec['kappa'])/G1
    W, Qv, Vv, R = (mpq(rec[k_]) for k_ in ('w', 'q', 'v', 'r'))
    e = mp.mpf(eta_str)
    ch, sh = mp.cosh(e), mp.sinh(e)
    A = G1*W*Qv + G2*Vv*R
    B = (G1*(W**2 + Qv**2) + G2*(Vv**2 + R**2))/2
    J = A*mp.cosh(2*e) + B*mp.sinh(2*e)
    sig = J/mpq(M_VAL)**2
    vals = dict(g1=G1, g2=G2, M=mpq(M_VAL), Lam=mpq(rec['Lambda']), f=mpq(fval), sb=sig, c14=mpq(fval), c123=mpq(fval))
    if frame == 'lab':
        vals.update(w=W, q=Qv, v=Vv, r=R, eb=e)
    else:
        vals.update(w=ch*W + sh*Qv, q=sh*W + ch*Qv, v=ch*Vv + sh*R, r=sh*Vv + ch*R, eb=mp.mpf(0))
    return vals


SYMS = dict(w=w, q=q, v=v, r=r, g1=g1, g2=g2, M=M, Lam=Lam, f=f, sb=sb, eb=eb, c14=c14, c123=c123)
ARGS = ('w', 'q', 'v', 'r', 'g1', 'g2', 'M', 'Lam', 'f', 'sb', 'eb', 'c14', 'c123')
ARGSYM = [SYMS[a_] for a_ in ARGS]


def argv(vals):
    return [vals[a_] for a_ in ARGS]


# ======================================================================================================================
# route R2: independent Lagrangian in mpmath and finite differences (S2.3)
# ======================================================================================================================
def lagr_r2(jet, par):
    """Independent mpmath implementation of B1 (not generated from the SymPy expression)."""
    p_t, p_x, q_t, q_x = jet[1], jet[2], jet[4], jet[5]
    s0, s_t, s_x = jet[6], jet[7], jet[8]
    e0, e_t, e_x = jet[9], jet[10], jet[11]
    C, S = mp.cosh(e0), mp.sinh(e0)
    up, np_ = C*p_t + S*p_x, S*p_t + C*p_x
    uq, nq = C*q_t + S*q_x, S*q_t + C*q_x
    val = (p_t*p_t - p_x*p_x)/2 + (q_t*q_t - q_x*q_x)/2 + (s_t*s_t - s_x*s_x)/2 - par['M']**2*s0*s0/2
    val += s0*(par['g1']*up*np_ + par['g2']*uq*nq) - par['Lam']*np_**4/4
    val += par['f']*(e_t*e_t - e_x*e_x)/2
    return val


def hessian_r2(par, h):
    base = [mp.mpf(0), par['w'], par['q'], mp.mpf(0), par['v'], par['r'], par['sb'], mp.mpf(0), mp.mpf(0),
            par['eb'], mp.mpf(0), mp.mpf(0)]
    hs = [h*max(mp.mpf(1), abs(b_)) for b_ in base]
    Hn = mp.matrix(12, 12)

    def Lat(shift):
        jet = list(base)
        for i_, d_ in shift:
            jet[i_] += d_
        return lagr_r2(jet, par)
    L0 = Lat([])
    for i in range(12):
        Hn[i, i] = (Lat([(i, hs[i])]) - 2*L0 + Lat([(i, -hs[i])]))/hs[i]**2
        for j in range(i + 1, 12):
            val = (Lat([(i, hs[i]), (j, hs[j])]) - Lat([(i, hs[i]), (j, -hs[j])])
                   - Lat([(i, -hs[i]), (j, hs[j])]) + Lat([(i, -hs[i]), (j, -hs[j])]))/(4*hs[i]*hs[j])
            Hn[i, j] = Hn[j, i] = val
    return Hn


def kmat_num(Hn, omega, kk):
    T = mp.matrix(12, 4)
    for i in range(4):
        T[3*i, i] = 1
        T[3*i + 1, i] = -1j*omega
        T[3*i + 2, i] = 1j*kk
    return T.H*Hn*T


def compare(a, b, atol, rtol):
    return abs(a - b) <= atol + rtol*max(abs(a), abs(b))


def dc1_background(Hfun, par):
    rows = []
    with mp.workdps(R2_DPS):
        rtol, atolf = mp.mpf(RTOL), mp.mpf(ATOL_FAC)
        H1 = mp.matrix(Hfun(*argv(par)))
        HR = [hessian_r2(par, mp.mpf(h_)) for h_ in R2_STEPS]
        for (o_, k_) in R2_POINTS:
            o_, k_ = mp.mpf(str(o_)), mp.mpf(str(k_))
            K1 = kmat_num(H1, o_, k_)
            Ka, Kb = (kmat_num(Hh, o_, k_) for Hh in HR)
            S = max(abs(K1[i, j]) for i in range(4) for j in range(4))
            atol = atolf*S
            for i in range(4):
                for j in range(4):
                    a_, b_, ba = K1[i, j], Kb[i, j], Ka[i, j]
                    if not compare(ba, b_, atol, rtol):
                        cls = 'R2 NOT CONVERGED'
                    elif abs(a_) <= atol and abs(b_) <= atol:
                        cls = 'ZERO (agree)' if compare(a_, b_, atol, rtol) else 'disagree'
                    else:
                        cls = 'agree' if compare(a_, b_, atol, rtol) else 'disagree'
                    rel = mp.nstr(abs(a_ - b_)/max(abs(a_), abs(b_)), 5) if max(abs(a_), abs(b_)) > 0 else None
                    rows.append(dict(point=[float(o_), float(k_)], entry=[i, j], cls=cls, abs_diff=mp.nstr(abs(a_ - b_), 5),
                                     rel_diff=rel, S=mp.nstr(S, 8), abs_R1=mp.nstr(abs(a_), 8), abs_R2=mp.nstr(abs(b_), 8)))
    cnt = {}
    for r_ in rows:
        cnt[r_['cls']] = cnt.get(r_['cls'], 0) + 1
    verdict = 'FAIL' if cnt.get('disagree') else ('INCONCLUSIVE' if cnt.get('R2 NOT CONVERGED') else 'PASS')
    return verdict, cnt, rows


# ======================================================================================================================
# characteristic polynomial, P1, P2a, P3, P4, P6 (S3)
# ======================================================================================================================
class DetPoly:
    """det K_u(omega, k) = sum d_ab omega^a k^b with coefficient lambdas (mpmath)."""

    def __init__(self, det_u):
        Pd = sp.Poly(det_u, om, k)
        self.terms = [(int(a_), int(b_), sp.lambdify(ARGSYM, cf, 'mpmath')) for (a_, b_), cf in Pd.terms()]
        self.deg_om = int(Pd.degree(om))

    def coeffs(self, vals):
        return [(a_, b_, fn(*argv(vals))) for a_, b_, fn in self.terms]

    @staticmethod
    def in_c(cf, kv):
        """coefficients of det(c k, k) in c, ascending."""
        deg = max(a_ for a_, _, _ in cf)
        out = [mp.mpf(0)]*(deg + 1)
        for a_, b_, d_ in cf:
            out[a_] += d_*kv**(a_ + b_)
        return out

    @staticmethod
    def at_k0(cf):
        deg = max(a_ for a_, _, _ in cf)
        out = [mp.mpf(0)]*(deg + 1)
        for a_, b_, d_ in cf:
            if b_ == 0:
                out[a_] += d_
        return out


def roots_from(co, dps):
    """roots of ascending coefficient list co at dps digits (None on non-convergence)."""
    while len(co) > 1 and co[-1] == 0:
        co = co[:-1]
    try:
        return mp.polyroots(co[::-1], maxsteps=600, extraprec=4*dps)
    except mp.libmp.NoConvergence:
        return None


def rel_imag(z, scale=1):
    return float(abs(mp.im(z))/max(mp.mpf(scale), abs(z)))


def classify_grid(DP, rec, eta_str, fval):
    """P3 at every grid point (double), with 50-digit confirmation of candidates and unresolved roots."""
    with mp.workdps(60):
        cf60 = DP.coeffs(numeric_bg(rec, eta_str, fval, 'u'))
    G = []
    for kv in KGRID:
        with mp.workdps(60):
            co = DP.in_c(cf60, mp.mpf(kv))
            sc = max(abs(x_) for x_ in co)
            cd = [complex(x_/sc) for x_ in co]
        rd = np.roots(cd[::-1])
        r_d = max(abs(z_.imag)/max(1.0, abs(z_)) for z_ in rd) if len(rd) else 0.0
        cls = 'real'
        g50 = None
        if r_d > TOL_REAL:
            with mp.workdps(50):
                cf = DP.coeffs(numeric_bg(rec, eta_str, fval, 'u'))
                rm = roots_from(DP.in_c(cf, mp.mpf(kv)), 50)
                if rm is None:
                    cls = 'nodata'
                else:
                    r_m = max(rel_imag(z_) for z_ in rm)
                    g50 = float(max(mp.mpf(0), max(mp.im(z_) for z_ in rm))*mp.mpf(kv))
                    if r_d <= TOL_UNRES:
                        cls = 'real' if r_m <= TOL_REAL else ('unresolved_stuck' if r_m <= TOL_UNRES else
                                                            ('complex' if r_m > TOL_CONF else 'real'))
                    else:
                        cls = 'complex' if r_m > TOL_CONF else 'real'
        G.append(dict(k=float(kv), r_double=float(r_d), cls=cls, gamma50=g50))
    return G


def gamma_at(DP, rec, eta_str, fval, kv, dps):
    with mp.workdps(dps):
        cf = DP.coeffs(numeric_bg(rec, eta_str, fval, 'u'))
        rt = roots_from(DP.in_c(cf, mp.mpf(kv)), dps)
        if rt is None:
            return None
        return +(max(mp.mpf(0), max(mp.im(z_) for z_ in rt))*mp.mpf(kv))


def p6(DP, rec, eta_str, fval, G):
    idx = [i for i, g_ in enumerate(G) if g_['cls'] == 'complex']
    if not idx:
        return None, None
    ks = [g_['k'] for g_ in G]
    lo, hi = max(0, min(idx) - 1), min(len(ks) - 1, max(idx) + 1)
    fine = []
    for i in range(lo, hi):
        fine += list(np.geomspace(ks[i], ks[i + 1], 5)[:-1])
    fine.append(ks[hi])
    rows = []
    for kv in fine:
        g50 = gamma_at(DP, rec, eta_str, fval, kv, 50)
        g100 = gamma_at(DP, rec, eta_str, fval, kv, 100)
        if g50 is None or g100 is None:
            rows.append((kv, None, None, None))
            continue
        with mp.workdps(120):
            rows.append((kv, g50, g100, abs(g50 - g100)))
    ok = [r_ for r_ in rows if r_[1] is not None]
    flags = [0 if g50 <= u_ else (1 if g50 <= 10*u_ else 2) for _, g50, _, u_ in ok]
    run = best = 0
    for fl in flags:
        run = run + 1 if fl == 2 else 0
        best = max(best, run)
    signchg = any((a_[1] > 0) != (b_[1] > 0) for a_, b_ in zip(ok, ok[1:])) and any(fl == 1 for fl in flags)
    if best >= 3:
        verdict = 'FAIL-UNSTABLE'
    elif any(fl == 1 for fl in flags) or signchg or best > 0 or len(ok) < len(rows):
        verdict = 'INCONCLUSIVE'
    else:
        verdict = 'PASS'
    with mp.workdps(120):
        gmax = max((r_[1] for r_ in ok), default=mp.mpf(0))
        info = dict(n_refined=len(fine), n_failed=len(rows) - len(ok), gamma_hat_max=mp.nstr(gmax, 10),
                    max_u=mp.nstr(max((r_[3] for r_ in ok), default=mp.mpf(0)), 8),
                    k_band=[G[min(idx)]['k'], G[max(idx)]['k']], max_consecutive_gt10u=best,
                    e_fold_time=(mp.nstr(1/gmax, 8) if gmax > 0 else None),
                    u_definition="|gamma50 - gamma100|; refinement part via sign-change / 3-consecutive rule (as stage 1)")
    return verdict, info


def p4_fit(G):
    pts = [(g_['k'], g_['gamma50']) for g_ in G if g_['cls'] == 'complex' and g_['k'] >= 1e4 and g_['gamma50']]
    if len(pts) < 2:
        return None
    s_, _ = np.polyfit(np.log([p_[0] for p_ in pts]), np.log([p_[1] for p_ in pts]), 1)
    return dict(slope=float(s_), n_points=len(pts), indicates_illposed=bool(s_ > 0.05))


def p2a(DP, rec, eta_str, fval):
    with mp.workdps(50):
        cf = DP.coeffs(numeric_bg(rec, eta_str, fval, 'u'))
        co = DP.at_k0(cf)
        nz = 0
        while nz < len(co) and co[nz] == 0:
            nz += 1
        rest = co[nz:]
        rt = roots_from(rest, 50) if len(rest) > 1 else []
        if rt is None:
            return dict(cls='INCONCLUSIVE', reason='no convergence')
        rmax = max((rel_imag(z_, 1) for z_ in rt), default=0.0)
        om2 = sorted({mp.nstr(mp.re(z_**2), 20) for z_ in rt}, key=lambda s_: float(s_))
        growth = [mp.nstr(mp.im(z_), 12) for z_ in rt if mp.im(z_) > 0 and rel_imag(z_, 1) > TOL_CONF]
    if rmax > TOL_CONF:
        cls = 'FAIL-UNSTABLE'
    elif nz == 4:
        cls = 'ZERO MODE'
    else:
        cls = 'INCONCLUSIVE'
    return dict(cls=cls, zero_root_multiplicity=nz, max_r_nonzero=rmax, nonzero_omega2=om2, growing_roots_im=growth)


# ======================================================================================================================
# principal symbol, C1, C2 (S3.2, S3.3)
# ======================================================================================================================
def c1_from_c2(c2):
    """C1 [ANALYTIC] only from a C2 PASS (which itself requires the block-diagonal principal symbol)."""
    if c2['status'] == 'PASS':
        return 'PASS [ANALYTIC]'
    return 'FAIL' if c2['status'] == 'FAIL' else 'INCONCLUSIVE'


def principal_symbolic(Ku):
    """degree-2 part of each entry of K_u; off-diagonal degree-2 parts must vanish (derived, not assumed)."""
    def deg2(e):
        Pe = sp.Poly(sp.expand(e), om, k)
        return sum(cf*om**a_*k**b_ for (a_, b_), cf in Pe.terms() if a_ + b_ == 2)
    offdiag_zero = all(sp.expand(deg2(Ku[i, j])) == 0 for i in range(4) for j in range(4) if i != j)
    blocks = []
    for i in range(4):
        d2 = sp.Poly(deg2(Ku[i, i]), om, k)
        blocks.append(tuple(sp.lambdify(ARGSYM, d2.coeff_monomial(mon), 'mpmath') for mon in (om**2, om*k, k**2)))
    diag_str = [str(sp.factor(deg2(Ku[i, i]))) for i in range(4)]
    return offdiag_zero, blocks, diag_str


def c2_numeric(blocks, vals, offdiag_zero):
    """principal roots c of A c^2 + B c + C per field (u frame) and the first-order principal matrix.
    The per-block analysis and the direct-sum argument are valid only if the principal symbol is block-diagonal."""
    if not offdiag_zero:
        return dict(status='INCONCLUSIVE', reason='principal symbol not block-diagonal: per-block analysis not valid')
    out = []
    Mfull = np.zeros((12, 12))
    complete = True
    real = True
    with mp.workdps(50):
        for i, (fa, fb, fc) in enumerate(blocks):
            A_, B_, C_ = fa(*argv(vals)), fb(*argv(vals)), fc(*argv(vals))
            if A_ == 0:
                return dict(status='INCONCLUSIVE', reason='vanishing time coefficient in block %s' % FIELDS[i])
            disc = B_**2 - 4*A_*C_
            if disc < 0:
                real = False
            rts = sorted([(-B_ + mp.sqrt(disc))/(2*A_), (-B_ - mp.sqrt(disc))/(2*A_)], key=lambda z_: float(mp.re(z_)))
            double = abs(rts[0] - rts[1]) <= COINC*max(1, abs(rts[0]), abs(rts[1]))
            if double:
                complete = False
            # first-order variables (F, F_t, F_x): d_t F = F_t (no principal part); companion block for (F_t, F_x)
            Mfull[3*i + 1, 3*i + 1] = float(-B_/A_)
            Mfull[3*i + 1, 3*i + 2] = float(-C_/A_)
            Mfull[3*i + 2, 3*i + 1] = 1.0
            out.append(dict(field=FIELDS[i], A=mp.nstr(A_, 15), B=mp.nstr(B_, 15), C=mp.nstr(C_, 15),
                            roots=[mp.nstr(z_, 15) for z_ in rts], discriminant=mp.nstr(disc, 15), double_in_block=bool(double)))
    allr = sorted(float(mp.mpf(z_)) for x_ in out for z_ in x_['roots']) if real else []
    across = bool(real and any(abs(a_ - b_) <= COINC*max(1, abs(a_), abs(b_)) for a_, b_ in zip(allr, allr[1:])))
    ev, V = np.linalg.eig(Mfull)
    cond = float(np.linalg.cond(V))
    if not real:
        status = 'FAIL'
    elif not complete:
        status = 'INCONCLUSIVE'
    else:
        status = 'PASS'
    return dict(status=status, blocks=out, principal_roots_real=real, complete_eigenvectors=bool(complete and real),
                coincident_across_blocks=across,
                coincidence_argument=("block-diagonal principal matrix; each block has distinct real roots and is diagonalisable;"
                                      " the direct sum is diagonalisable, so roots shared by different blocks keep a complete"
                                      " eigenvector set") if (across and complete and real) else None,
                cond_eigvec=cond, cond_above_1e8=bool(cond > COND_THR),
                first_order_form="U = (F, F_t, F_x) per field, u frame; dU/dt + A_p dU/dx = lower order",
                normalisation="unit coefficient of 1/2 F_t^2 except eta (f) -- the companion block divides by A")


# ======================================================================================================================
# T5: energy form, Hamiltonian, Krein (S4)
# ======================================================================================================================
def t5_frame(Hfun, vals):
    with mp.workdps(50):
        Hn = mp.matrix(Hfun(*argv(vals)))
        tol = mp.mpf(T5_TOL)

        def blocks(kk):
            Tt = mp.matrix(12, 4)
            Ty = mp.matrix(12, 4)
            for i in range(4):
                Tt[3*i + 1, i] = 1
                Ty[3*i, i] = 1
                Ty[3*i + 2, i] = 1j*kk
            A = Tt.H*Hn*Tt
            B = Tt.H*Hn*Ty
            C = -(Ty.H*Hn*Ty)
            return A, B, C

        def herm(X):
            return (X + X.H)/2

        def mineig(X):
            return min(mp.re(e_) for e_ in mp.eighe(herm(X), eigvals_only=True))

        A, _, C0 = blocks(mp.mpf(0))
        _, _, Cp = blocks(mp.mpf(1))
        _, _, Cm = blocks(mp.mpf(-1))
        C1 = (Cp - Cm)/2
        C2 = (Cp + Cm)/2 - C0
        a_min = mineig(A)
        c2_min = mineig(C2)
        # k -> 0: light (phi, psi) and heavy (sigma, eta) blocks of the potential form
        li, hv = [0, 1], [2, 3]
        sub = lambda X, I_, J_: mp.matrix([[X[i, j] for j in J_] for i in I_])
        light0 = max(abs(C0[i, j]) for i in li for j in range(4)) + max(abs(C1[i, j]) for i in li for j in li)
        mh = sub(C0, hv, hv)
        mh_min = mineig(mh)
        schur = None
        if mh_min > tol:
            schur = mineig(sub(C2, li, li) - sub(C1, li, hv)*mp.inverse(mh)*sub(C1, hv, li))
        grid_min, worst_k = None, None
        for kv in KGRID:
            _, _, Ck = blocks(mp.mpf(kv))
            e_ = mineig(Ck)/max(1, mp.mnorm(Ck, 1))
            if grid_min is None or e_ < grid_min:
                grid_min, worst_k = e_, float(kv)
        if a_min < -tol or c2_min < -tol or (schur is not None and schur < -tol) or grid_min < -tol or mh_min < -tol:
            static = 'FAIL'
        elif abs(a_min) <= tol or abs(c2_min) <= tol or schur is None or light0 != 0:
            static = 'INCONCLUSIVE'
        else:
            static = 'PASS'
        ham, krein = {}, {}
        if a_min > tol:
            Jm = mp.matrix(8, 8)
            for i in range(4):
                Jm[i, 4 + i] = 1
                Jm[4 + i, i] = -1
            for kv in KS_H:
                A_, B_, C_ = blocks(mp.mpf(kv))
                Ai = mp.inverse(A_)
                Hm = mp.matrix(8, 8)
                X11 = C_ + B_.H*Ai*B_
                X12 = -(B_.H*Ai)
                X21 = -(Ai*B_)
                for i in range(4):
                    for j in range(4):
                        Hm[i, j], Hm[i, 4 + j], Hm[4 + i, j], Hm[4 + i, 4 + j] = X11[i, j], X12[i, j], X21[i, j], Ai[i, j]
                Hm = herm(Hm)
                ev = mp.eighe(Hm, eigvals_only=True)
                ham[str(kv)] = dict(n_negative=int(sum(1 for e_ in ev if mp.re(e_) < -tol)), min_eig=mp.nstr(min(mp.re(e_) for e_ in ev), 10))
                omg, Vv = mp.eig(1j*Jm*Hm)
                neg = pos = cpx = 0
                for j in range(8):
                    if rel_imag(omg[j], 1) > TOL_CONF:
                        cpx += 1
                        continue
                    vec = Vv[:, j]
                    en = mp.re((vec.H*Hm*vec)[0, 0])
                    if en < 0:
                        neg += 1
                    else:
                        pos += 1
                krein[str(kv)] = dict(real_modes_positive_energy=pos, real_modes_negative_energy=neg, complex_modes=cpx)
        return dict(T5_static=static, A_min_eig=mp.nstr(a_min, 12), C2_min_eig=mp.nstr(c2_min, 12),
                    heavy_block_k0_min_eig=mp.nstr(mh_min, 12), schur_k0_min_eig=(mp.nstr(schur, 12) if schur is not None else None),
                    light_rows_at_k0_zero=bool(light0 == 0), grid_min_normalised=mp.nstr(grid_min, 12), grid_worst_k=worst_k,
                    T5_hamiltonian=ham if ham else 'NOT EVALUATED (velocity block A not positive definite)',
                    T5_krein=krein if krein else 'NOT EVALUATED (velocity block A not positive definite)')


# ======================================================================================================================
# S6, S7: eta mass, eliminations, comparison with B4
# ======================================================================================================================
def dedupe(xs, tol):
    out = []
    for x_ in xs:
        if not any(abs(x_ - y_) <= tol for y_ in out):
            out.append(x_)
    return out


def dc5_match(found, jprime, jzero, tol):
    """DC5: the distinct found roots must equal the union of the J' = 0 and J = 0 classes (no double counting)."""
    F_ = dedupe(found, tol)
    U_ = dedupe(list(jprime) + list(jzero), tol)
    every_theory_found = all(any(abs(f_ - u_) <= tol for f_ in F_) for u_ in U_)
    every_found_explained = all(any(abs(f_ - u_) <= tol for u_ in U_) for f_ in F_)
    jp_ok = all(any(abs(f_ - u_) <= tol for f_ in F_) for u_ in jprime)
    return dict(n_found_distinct=len(F_), n_union_distinct=len(U_), Jprime_reproduced=bool(jp_ok),
                all_union_found=bool(every_theory_found), no_extra_roots=bool(every_found_explained),
                pass_=bool(every_theory_found and every_found_explained and len(F_) == len(U_)))


def t6_class(ok, outc, ecls, k0cls, conf):
    """T6 per (background, root, f).  STATED-DYNAMICAL needs the full-system confirmation conf is True."""
    if not ok:
        return 'INCONCLUSIVE'
    if outc in ('FAIL-ILLPOSED', 'FAIL-UNSTABLE'):
        return 'FAIL'
    if ecls in ('LIGHT', 'HEAVY') and k0cls == 'ZERO MODE' and outc.startswith('PASS') and conf is True:
        return 'STATED-DYNAMICAL'
    return 'INCONCLUSIVE'


def t7_class(ok, ecls, s7n, s73):
    """T7 per (background, root, f).  TACHYONIC: the elimination procedure is not admissible at a maximum of V_eff;
    this is not a proof that no UV origin exists."""
    if not ok or ecls == 'BOUNDARY' or s7n.get('singular') or s7n.get('error'):
        return 'INCONCLUSIVE'
    if ecls in ('LIGHT', 'TACHYONIC'):
        return 'NOT-REPRODUCED'
    if s7n['eta_exchange_zero'] and s73['operator_identity_with_B4_rank_one'] and s73['characteristic_polynomial_identical']:
        return 'REPRODUCED-RANK-ONE'
    return 'NOT-REPRODUCED'


def veff_derivs(rec, eta_str):
    """V(eta) = J^2/(2M^2) - Lam/4 N^4 and its first two eta-derivatives, 60 digits."""
    with mp.workdps(60):
        p_ = mpq(rec['p'])
        G1 = mp.sqrt(p_)
        G2 = -mpq(rec['kappa'])/G1
        W, Qv, Vv, R = (mpq(rec[k_]) for k_ in ('w', 'q', 'v', 'r'))
        La, Mv = mpq(rec['Lambda']), mpq(M_VAL)
        A = G1*W*Qv + G2*Vv*R
        B = (G1*(W**2 + Qv**2) + G2*(Vv**2 + R**2))/2
        Vf = lambda e: (A*mp.cosh(2*e) + B*mp.sinh(2*e))**2/(2*Mv**2) - La/4*(mp.sinh(e)*W + mp.cosh(e)*Qv)**4
        e = mp.mpf(eta_str)
        return Vf(e), mp.diff(Vf, e, 1), mp.diff(Vf, e, 2)


def eta_class(m2):
    with mp.workdps(60):
        b_ = mp.mpf(BOUNDARY)
        thr = mpq(LIGHT_THR)*mpq(M_VAL)**2
        if abs(m2) <= b_ or abs(m2 - thr) <= b_:
            return 'BOUNDARY'
        if m2 < 0:
            return 'TACHYONIC'
        return 'LIGHT' if m2 < thr else 'HEAVY'


X4 = sp.symbols('X1 X2 X3 X4', real=True)      # (u.dphi, n.dphi, u.dpsi, n.dpsi) = (phi_t', phi_x', psi_t', psi_x')
ETA_S = sp.symbols('eta_s', real=True)


def s7_symbolic():
    """S7.1 elimination of sigma; S7.3(a),(b) frozen-eta comparison with B4 (u frame), exact."""
    t, x = sp.symbols('t x', real=True)
    Jf = sp.Function('J')(t, x)
    boxJ = sp.diff(Jf, t, 2) - sp.diff(Jf, x, 2)
    sig2 = Jf/M**2 - boxJ/M**4                                     # classical solution to O(M^-4)
    Lsig = sp.Rational(1, 2)*(sp.diff(sig2, t)**2 - sp.diff(sig2, x)**2) - M**2/2*sig2**2 + sig2*Jf
    target = Jf**2/(2*M**2) + sp.Rational(1, 2)*(sp.diff(Jf, t)**2 - sp.diff(Jf, x)**2)/M**4
    eps = sp.symbols('eps', positive=True)
    diff = sp.expand((Lsig - target).subs(M, 1/eps))
    rem = sp.expand(sum(diff.coeff(eps, n_) * eps**n_ for n_ in range(0, 5)))   # terms up to O(M^-4)
    from sympy.calculus.euler import euler_equations
    td = euler_equations(rem, Jf, [t, x])[0].lhs if rem != 0 else 0
    s71 = dict(classical_leading="sigma = J_tot/M^2 + O(d^2/M^4)",
               gaussian_nonlocal="+1/2 J_tot (box + M^2)^-1 J_tot  (signature (+,-), sign +)",
               expansion="(box + M^2)^-1 = M^-2 - box M^-4 + O(M^-6)",
               leading_local_term="J_tot^2/(2 M^2)",
               remainder_O_M4="+1/2 (dJ_tot)^2 / M^4  (= -1/2 J_tot box J_tot / M^4 up to a total derivative)",
               check_classical_equals_expansion_to_M4=bool(sp.simplify(td) == 0))
    # frozen eta, u frame (eta_bar' = 0): leading local effective Lagrangian
    Jt = g1*X4[0]*X4[1] + g2*X4[2]*X4[3]
    Leff = sp.Rational(1, 2)*(X4[0]**2 - X4[1]**2 + X4[2]**2 - X4[3]**2) + Jt**2/(2*M**2) - Lam/4*X4[1]**4
    pS, sS, kS = sp.symbols('p_B4 s_B4 kappa_B4', real=True)
    LB4 = (sp.Rational(1, 2)*(X4[0]**2 - X4[1]**2) - Lam/4*X4[1]**4 + sp.Rational(1, 2)*(X4[2]**2 - X4[3]**2)
           + pS/2*X4[0]**2*X4[1]**2 - kS*X4[0]*X4[1]*X4[2]*X4[3] + sS/2*X4[2]**2*X4[3]**2)
    match = {pS: g1**2/M**2, sS: g2**2/M**2, kS: -g1*g2/M**2}
    op_identity = sp.expand(Leff - LB4.subs(match)) == 0
    pp, kk_ = sp.symbols('p kappa', positive=True)
    rank_one = sp.simplify((g2**2/M**2).subs({g1: sp.sqrt(pp), g2: -kk_/sp.sqrt(pp), M: 1}) - kk_**2/pp) == 0
    coeffs = {str(mon): str(sp.factor(cf)) for mon, cf in sp.Poly(sp.expand(Leff - sp.Rational(1, 2)*(X4[0]**2 - X4[1]**2 + X4[2]**2 - X4[3]**2)), *X4).terms()}
    # characteristic polynomial of the frozen-eta effective operator vs F of B4 (lambda_eff.py), exact
    Wb, Qb, Vb, Rb = sp.symbols('wp qp vp rp', real=True)
    bgX = {X4[0]: Wb, X4[1]: Qb, X4[2]: Vb, X4[3]: Rb}
    Hx = sp.hessian(Leff, X4).subs(bgX)
    tv = sp.Matrix([-c, 1])
    Km = sp.zeros(2, 2)
    for i in range(2):
        for j in range(2):
            Km[i, j] = (tv.T*Hx[2*i:2*i + 2, 2*j:2*j + 2]*tv)[0, 0]
    det_eff = sp.expand(Km.det())
    pB, sB, kB = g1**2/M**2, g2**2/M**2, -g1*g2/M**2
    aphi, apsi = pB*Wb*Qb - kB*Vb*Rb, -kB*Wb*Qb + sB*Vb*Rb
    d1 = c**2 - 2*aphi*c - (1 + 3*Lam*Qb**2)
    d2 = c**2 - 2*apsi*c - 1
    Fc = d1*d2 + pB*(Wb - Qb*c)**2*d2 + sB*(Vb - Rb*c)**2*d1
    ratio = sp.simplify(det_eff/sp.expand(Fc))
    phi_x2 = sp.expand(sp.diff(Leff, X4[1], 2).subs(bgX).subs({X4[0]: Wb}))
    lam_norm = sp.simplify(-(sp.diff(sp.Rational(1, 2)*(-X4[1]**2) - Lam/4*X4[1]**4, X4[1], 2)).subs(X4[1], Qb) - (1 + 3*Lam*Qb**2)) == 0
    s73 = dict(frame="u frame, eta_bar fixed (frozen), signature (+,-); B4 form docs/notation_and_conventions.md with m = 0",
               operator_identity_with_B4_rank_one=bool(op_identity),
               matching="p = g1^2/M^2, s_eff = g2^2/M^2 = kappa^2/p, kappa = -g1 g2/M^2, Lambda unchanged (bare)",
               rank_one_identity_s_eff_eq_kappa2_over_p=bool(rank_one),
               effective_quartic_coefficients=coeffs,
               characteristic_polynomial_ratio_to_F=str(ratio),
               characteristic_polynomial_identical=bool(ratio == 1),
               lambda_term_gives_1_plus_3Lam_q2_in_d1=bool(lam_norm),
               psi_mass_term="not generated; B4 records carry no m; compared at m = 0",
               d2Leff_dphix2=str(phi_x2))
    return dict(S7_1=s71, S7_3=s73)


def s7_numeric(rec, eta_str, DP, fval):
    """eta-exchange correction (Schur) and low-k roots: full system vs eliminated vs frozen-eta B4 polynomial."""
    with mp.workdps(60):
        vals = numeric_bg(rec, eta_str, fval, 'u')
        G1, G2, La, Mv = vals['g1'], vals['g2'], vals['Lam'], vals['M']
        Xb = [vals['w'], vals['q'], vals['v'], vals['r']]

        def Ffun(e, X):
            C, S = mp.cosh(e), mp.sinh(e)
            u1, n1 = C*X[0] + S*X[1], S*X[0] + C*X[1]
            u2, n2 = C*X[2] + S*X[3], S*X[2] + C*X[3]
            Jv = G1*u1*n1 + G2*u2*n2
            return Jv**2/(2*Mv**2) - La/4*n1**4
        e0 = mp.mpf(0)
        Fee = mp.diff(lambda e: Ffun(e, Xb), e0, 2)
        if abs(Fee) <= mp.mpf(BOUNDARY):
            return dict(singular=True, F_eta_eta=mp.nstr(Fee, 15),
                        reason='F_eta_eta = 0 within 1e-40: eta exchange undefined; T7 INCONCLUSIVE for this root and f')
        HX = mp.matrix(4, 4)
        bx = mp.matrix(4, 1)
        for i in range(4):
            def Fi(e, xi, i=i):
                X = list(Xb)
                X[i] = xi
                return Ffun(e, X)
            bx[i] = mp.diff(lambda e, xi: Fi(e, xi), (e0, Xb[i]), (1, 1))
            for j in range(4):
                def Fij(a_, b_, i=i, j=j):
                    X = list(Xb)
                    X[i] = X[i] + a_
                    X[j] = X[j] + b_
                    return Ffun(e0, X)
                HX[i, j] = mp.diff(Fij, (0, 0), (1, 1)) if i != j else mp.diff(lambda a_: Fij(a_, 0), 0, 2)
        corr = bx*bx.T/Fee
        rel = mp.mnorm(corr, 1)/max(mp.mnorm(HX, 1), mp.mpf('1e-300'))
        kin = mp.diag([1, -1, 1, -1])

        def char_roots(Hm):
            # 2x2 operator for plane waves, (dF_t, dF_x) = i k (-c, 1) F
            cs = [mp.mpf(i_) for i_ in range(-2, 3)]
            vals_ = []
            for cc in cs:
                tvec = [-cc, 1]
                Km = mp.matrix(2, 2)
                for a_ in range(2):
                    for b_ in range(2):
                        Km[a_, b_] = sum(tvec[i_]*Hm[2*a_ + i_, 2*b_ + j_]*tvec[j_] for i_ in range(2) for j_ in range(2))
                vals_.append(mp.det(Km))
            # the determinant is a polynomial of degree 4 in c: exact interpolation through 5 nodes
            V = mp.matrix([[cc**p_ for p_ in range(5)] for cc in cs])
            sol = mp.lu_solve(V, mp.matrix(vals_))
            return sorted(mp.polyroots([sol[i_] for i_ in range(4, -1, -1)], maxsteps=400, extraprec=200), key=lambda z_: float(mp.re(z_)))
        frozen = char_roots(kin + HX)
        elim = char_roots(kin + HX - corr)
        cf = DP.coeffs(vals)
        nmin = min(a_ + b_ for a_, b_, d_ in cf if d_ != 0)
        q0 = [mp.mpf(0)]*(DP.deg_om + 1)
        for a_, b_, d_ in cf:
            if a_ + b_ == nmin:
                q0[a_] += d_
        full = sorted(roots_from(q0, 60) or [], key=lambda z_: float(mp.re(z_)))
        dif = lambda A_, B_: (mp.nstr(max(abs(a_ - b_) for a_, b_ in zip(A_, B_)), 6) if len(A_) == len(B_) and A_ else 'length mismatch')
        return dict(eta_exchange_relative_norm=mp.nstr(rel, 8), eta_exchange_zero=bool(rel <= mp.mpf(ETA_EXCH_ZERO)),
                    F_eta_eta=mp.nstr(Fee, 15),
                    lowk_full_system_roots=[mp.nstr(z_, 15) for z_ in full],
                    roots_eta_eliminated=[mp.nstr(z_, 15) for z_ in elim],
                    roots_eta_frozen_B4=[mp.nstr(z_, 15) for z_ in frozen],
                    X2_full_vs_eliminated_maxdiff=dif(full, elim),
                    full_vs_frozen_B4_maxdiff=dif(full, frozen),
                    lowk_order_in_k=int(nmin))


# ======================================================================================================================
# main
# ======================================================================================================================
def main(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    T0 = time.time()
    log = []

    def note(msg):
        line = "[%8.1fs] %s" % (time.time() - T0, msg)
        log.append(line)
        print(line, flush=True)

    res = dict(evidence_level="1 (theoretical consistency); DC0-DC6 level 2; not experimental validation of B4",
               spec=SPEC, spec_commit=SPEC_COMMIT, candidate=CANDIDATE_SPEC, revision=REVISION,
               environment=dict(python=platform.python_version(), sympy=sp.__version__, mpmath=mp.__version__, numpy=np.__version__),
               inputs={DATA_B3: sha_file(DATA_B3), DATA_MAN: sha_file(DATA_MAN)},
               fixed=dict(f=[str(x_) for x_ in F_VALUES], M=1, P2=[1e-3, 1e3, 61, 1e6, 31], P3=[TOL_REAL, TOL_UNRES, TOL_CONF],
                          R2=dict(dps=R2_DPS, steps=R2_STEPS, rtol=RTOL, atol_factor=ATOL_FAC, points=R2_POINTS),
                          C0_tol=C0_TOL, light_threshold="1e-2 M^2", boundary=BOUNDARY, T5_tol=T5_TOL, k_H=KS_H,
                          aether="implemented, NOT run; c14 = c123 = f"),
               C4="NOT EVALUATED (script B)")
    checks = {}
    verify_inputs()                                     # ABORT before any computation if an input file changed
    note("spec %s (%s); revision %s; sympy %s; mpmath %s" % (SPEC, SPEC_COMMIT, REVISION, sp.__version__, mp.__version__))
    note("inputs verified against fixed SHA-256: %s" % EXPECTED_SHA256)

    # ---------------- S0 ----------------
    s0 = s0_checks()
    checks['DC0'] = bool(s0['DC0'])
    checks['DC3'] = bool(s0['DC3'])
    res['S0'] = s0
    note("S0: DC3 %s, DC0 %s" % (checks['DC3'], checks['DC0']))

    # ---------------- symbolic system ----------------
    SYM_ = build_symbolic()
    note("symbolic K and det K_u built in %.1fs" % SYM_['seconds'])
    # S1.1 background equation from the action
    J_ = g1*(sp.cosh(eb)*w + sp.sinh(eb)*q)*(sp.sinh(eb)*w + sp.cosh(eb)*q) + g2*(sp.cosh(eb)*v + sp.sinh(eb)*r)*(sp.sinh(eb)*v + sp.cosh(eb)*r)
    A_s, B_s = g1*w*q + g2*v*r, sp.Rational(1, 2)*(g1*(w**2 + q**2) + g2*(v**2 + r**2))
    N_ = sp.sinh(eb)*w + sp.cosh(eb)*q
    Vsym = J_**2/(2*M**2) - Lam/4*N_**4
    eta_eq_action = SYM_['grad'][JETS.index(E_)]
    sig_eq_action = SYM_['grad'][JETS.index(S_)]
    s11 = dict(J_form=zero_expr(J_ - (A_s*sp.cosh(2*eb) + B_s*sp.sinh(2*eb))),
               eta_eq_equals_dV_after_sigma=zero_expr(eta_eq_action.subs(sb, J_/M**2) - sp.diff(Vsym, eb)),
               eta_eq_form=zero_expr(eta_eq_action - (sb*sp.diff(J_, eb) - Lam*N_**3*sp.diff(N_, eb))),
               sigma_eq=zero_expr(sig_eq_action - (J_ - M**2*sb)))
    res['S1_1'] = s11
    note("S1.1: %s" % s11)
    # DC2: aether path at c14 = c123 = f
    dc2_H = all(zero_expr(SYM_['Hae'][i, j].subs({c14: f, c123: f}) - SYM_['H'][i, j]) for i in range(12) for j in range(12))
    Kae_u = SYM_['Kae'].subs(eb, 0).subs({c14: f, c123: f}).applyfunc(sp.expand)
    dc2_det = sp.expand(Kae_u.det(method='berkowitz') - SYM_['det_u']) == 0
    m2_ae = -(1/c14)*sp.diff(Vsym, eb, 2)
    dc2_m2 = sp.simplify(m2_ae.subs(c14, f) - (-(1/f)*sp.diff(Vsym, eb, 2))) == 0
    checks['DC2'] = bool(dc2_H and dc2_det and dc2_m2)
    res['DC2'] = dict(hessian_identical=dc2_H, det_K_u_identical=dc2_det, m2_formula_identical=dc2_m2,
                      aether_eta_block_u_frame=str(sp.factor(sp.expand(Kae_u[3, 3].subs({w: 0, q: 0, v: 0, r: 0})))))
    note("DC2: %s" % res['DC2'])
    # S2.1 mixing blocks (lab frame, symbolic)
    Hs = SYM_['H']
    nm = [str(j_) for j_ in JETS]
    mix = {}
    for (fa, ia), (fb, ib) in ((('sigma', 6), ('phi', 0)), (('sigma', 6), ('psi', 3)), (('eta', 9), ('phi', 0)), (('eta', 9), ('psi', 3)), (('sigma', 6), ('eta', 9))):
        mix['%s-%s' % (fa, fb)] = {'%s,%s' % (nm[i], nm[j]): str(sp.simplify(Hs[i, j])) for i in range(ia, ia + 3) for j in range(ib, ib + 3) if Hs[i, j] != 0}
    res['S2_1_mixing_blocks_lab'] = mix
    # P1
    t1 = time.time()
    fl = sp.factor_list(SYM_['det_u'], om, k)
    res['P1'] = dict(n_factors=len(fl[1]), factor_degrees=[[int(sp.Poly(a_, om, k).total_degree()), int(m_)] for a_, m_ in fl[1]],
                     constant=str(fl[0]), seconds=round(time.time() - t1, 1),
                     note="factorisation over Q(background symbols); a factor that exists only for special numeric values is not divided out")
    det_used = sp.expand(sp.Mul(*[a_ for a_, m_ in fl[1] if sp.Poly(a_, om, k).total_degree() > 0]))   # square-free product:
    # used ONLY for distinct roots (P3, P4, P6, low-k roots); P2a and X1 use the full determinant DP_full
    if len(fl[1]) > 1:
        note("P1: det K_u factorises: %s" % res['P1']['factor_degrees'])
    note("P1: %s" % res['P1'])
    DP = DetPoly(det_used)
    DP_full = DetPoly(SYM_['det_u'])
    off0, pblocks, pdiag = principal_symbolic(SYM_['Ku'])
    res['S3_2_principal_symbolic'] = dict(offdiagonal_degree2_parts_zero=bool(off0), diagonal_principal_parts_u_frame=dict(zip(FIELDS, pdiag)))
    note("S3.2: off-diagonal principal parts zero: %s; %s" % (off0, pdiag))
    Hfun = sp.lambdify(ARGSYM, SYM_['H'], 'mpmath')
    gradfun = sp.lambdify(ARGSYM, sp.Matrix(SYM_['grad']), 'mpmath')
    detlab = sp.lambdify(ARGSYM + [om, k], sp.Matrix(SYM_['K']).det(method='berkowitz'), 'mpmath')
    res.update(s7_symbolic())
    note("S7 symbolic: %s" % {k_: v_ for k_, v_ in res['S7_3'].items() if isinstance(v_, bool)})

    # ---------------- backgrounds ----------------
    bgs = load_backgrounds()
    rows = []
    dc4_ok = True
    dc5_ok = True
    dc1_all = []
    dc6_ok = True
    for bi, (name, rec) in enumerate(bgs):
        en, bf = enumerate_roots(rec)
        row = dict(name=name, inputs={k_: rec[k_ + '_str'] for k_ in ('p', 's', 'kappa', 'Lambda', 'w', 'q', 'v', 'r')},
                   rationalisation="sympy.Rational(str(x))", S1_2=en)
        s_minus = rec['s'] - rec['kappa']**2/rec['p']
        row['s_minus_kappa2_over_p'] = str(sp.N(s_minus, 20))
        # S1.4 lab-frame condition
        F4 = 2*bf['a']*bf['b']/(rec['p']*M_VAL**2) - rec['Lambda']*rec['q']**3*rec['w']
        row['S1_4'] = dict(value=str(sp.N(F4, 30)), eta_zero_is_root=bool(F4 == 0))
        # DC5 Lambda = 0
        bf0 = background_functions(rec, Lv=sp.Integer(0))
        a_, b_ = bf['a'], bf['b']
        z0 = [sp.N(z_, 60) for z_ in bf0['P'].real_roots() if sp.N(z_, 60) > 0]
        with mp.workdps(60):
            etas0 = [mp.log(mp.mpf(str(z_))) for z_ in z0]
            jprime = [] if abs(b_) >= abs(a_) else [mp.atanh(-mpq(b_)/mpq(a_))/2]
            jzero = [] if abs(a_) >= abs(b_) else [mp.atanh(-mpq(a_)/mpq(b_))/2]
            m5 = dc5_match(etas0, jprime, jzero, mp.mpf(DC5_TOL))
        dc5_bg = m5['pass_']
        dc5_ok = dc5_ok and dc5_bg
        row['DC5'] = dict(m5, n_roots_Lambda0=len(etas0), Jprime_roots=[mp.nstr(x_, 40) for x_ in jprime],
                          J_zero_roots=[mp.nstr(x_, 40) for x_ in jzero], historical_U1_U4_has_stationary=bool(len(jprime) == 1))
        # DC4 part 1
        dc4_bg = en['sturm_count'] == en['n_located']
        row['per_root'] = []
        if en['n_located'] == 0:
            row['T3fund'] = {str(fv): 'INCONCLUSIVE (C0: no real eta_bar)' for fv in F_VALUES}
            row['T6'] = 'INCONCLUSIVE (item 4: no stationary frame)'
        for ri, rt in enumerate(en['roots']):
            eta_str = rt['eta_bar']
            resolved = rt['resolved']
            rr = dict(root_index=ri, eta_bar=mp.nstr(mp.mpf(eta_str), 30), multiplicity=rt['multiplicity'],
                      resolved=resolved, edge_flag=bool(abs(mp.mpf(eta_str)) > ETA_EDGE),
                      note=(None if resolved else 'multiple or unresolved root: every class INCONCLUSIVE'))
            # C0 residual (S1.3) from the action
            with mp.workdps(50):
                vals = numeric_bg(rec, eta_str, F_VALUES[1], 'lab')
                gv = gradfun(*argv(vals))
                res_sig, res_eta = abs(gv[JETS.index(S_)]), abs(gv[JETS.index(E_)])
                c0_ok = bool(res_sig <= C0_TOL and res_eta <= C0_TOL)
                rr['C0'] = dict(residual_sigma=mp.nstr(res_sig, 5), residual_eta=mp.nstr(res_eta, 5),
                                phi_psi_equations="identically zero for constant gradients (analytic)", pass_=c0_ok)
                sig_bar = vals['sb']
                rr['sigma_bar'] = mp.nstr(sig_bar, 20)
                rr['u_frame_gradients'] = {k_: mp.nstr(v_, 20) for k_, v_ in numeric_bg(rec, eta_str, 1, 'u').items() if k_ in ('w', 'q', 'v', 'r')}
            dc4_bg = dc4_bg and c0_ok
            V0, V1_, V2_ = veff_derivs(rec, eta_str)
            rr['S5'] = dict(Veff_second_derivative=mp.nstr(-V2_, 20), extremum=('minimum of V_eff' if -V2_ > 0 else ('maximum of V_eff' if -V2_ < 0 else 'degenerate')),
                            dV_at_root=mp.nstr(V1_, 5))
            rr['per_f'] = {}
            for fv in F_VALUES:
                tag = str(fv)
                pf = {}
                with mp.workdps(60):
                    m2 = -V2_/mpq(fv)
                pf['m2_eta_eff'] = mp.nstr(m2, 30)
                pf['eta_class'] = eta_class(m2)
                # DC1 (lab frame, R1 vs R2)
                with mp.workdps(R2_DPS):
                    parv = numeric_bg(rec, eta_str, fv, 'lab')
                v_, cnt, drows = dc1_background(Hfun, parv)
                dc1_all.append(v_)
                pf['DC1'] = dict(verdict=v_, counts=cnt, entries=drows)
                # X1 frame consistency
                with mp.workdps(50):
                    pl = numeric_bg(rec, eta_str, fv, 'lab')
                    pu = numeric_bg(rec, eta_str, fv, 'u')
                    e_ = pl['eb']
                    x1 = []
                    for (o_, k_) in R2_POINTS:
                        o_, k_ = mp.mpf(str(o_)), mp.mpf(str(k_))
                        dl = detlab(*argv(pl), o_, k_)
                        op_, kp_ = mp.cosh(e_)*o_ - mp.sinh(e_)*k_, mp.cosh(e_)*k_ - mp.sinh(e_)*o_
                        du = sum(d_*op_**a_*kp_**b_ for a_, b_, d_ in DP_full.coeffs(pu))
                        x1.append(float(abs(dl - du)/max(abs(dl), mp.mpf('1e-300'))))
                pf['X1_frame_consistency_max_rel'] = max(x1)
                # C2, C1
                with mp.workdps(50):
                    c2 = c2_numeric(pblocks, numeric_bg(rec, eta_str, fv, 'u'), off0)
                pf['C2'] = c2
                c1 = c1_from_c2(c2)
                # P2a, P3, P4, P6
                k0 = p2a(DP_full, rec, eta_str, fv)      # full determinant: multiplicity of omega = 0 kept
                pf['P2a'] = k0
                G = classify_grid(DP, rec, eta_str, fv)
                cnts = {}
                for g_ in G:
                    cnts[g_['cls']] = cnts.get(g_['cls'], 0) + 1
                pf['P3_counts'] = cnts
                pf['P4'] = p4_fit(G)
                if c1 != 'PASS [ANALYTIC]' and not any(g_['cls'] in ('complex', 'unresolved_stuck', 'nodata') for g_ in G):
                    c1 = c1 if c1 == 'FAIL' else 'INCONCLUSIVE (numerically passed on [1e-3, 1e6])'
                pf['C1'] = c1
                v6, info6 = p6(DP, rec, eta_str, fv, G)
                pf['P6'] = info6
                if v6 is None:
                    c3 = 'INCONCLUSIVE' if any(g_['cls'] in ('unresolved_stuck', 'nodata') for g_ in G) else 'PASS'
                else:
                    c3 = v6
                if k0['cls'] == 'FAIL-UNSTABLE':
                    c3 = 'FAIL-UNSTABLE'
                elif k0['cls'] == 'INCONCLUSIVE' and c3 == 'PASS':
                    c3 = 'INCONCLUSIVE'
                pf['C3'] = c3
                if not (c0_ok and resolved):
                    outc = 'INCONCLUSIVE'
                elif c1 == 'FAIL' or c2['status'] == 'FAIL':
                    outc = 'FAIL-ILLPOSED'
                elif c3 == 'FAIL-UNSTABLE':
                    outc = 'FAIL-UNSTABLE'
                elif c1 == 'PASS [ANALYTIC]' and c2['status'] == 'PASS' and c3 == 'PASS':
                    outc = 'PASS (C0-C3; C4 not evaluated)'
                else:
                    outc = 'INCONCLUSIVE'
                pf['T3fund_outcome'] = outc
                # T5
                pf['T5_u_frame'] = t5_frame(Hfun, numeric_bg(rec, eta_str, fv, 'u'))
                pf['T5_lab_frame'] = (t5_frame(Hfun, numeric_bg(rec, eta_str, fv, 'lab')) if mp.mpf(eta_str) != 0
                                      else 'identical to u frame (eta_bar = 0)')
                # S6.3 confirmation from the full system at k = 0
                om2 = [mp.mpf(x_) for x_ in k0.get('nonzero_omega2', [])]
                conf = None
                if om2:
                    if pf['eta_class'] == 'TACHYONIC':
                        conf = any(x_ < 0 for x_ in om2)
                    elif pf['eta_class'] == 'LIGHT':
                        conf = any(0 <= x_ < 1e-2 for x_ in om2)
                    elif pf['eta_class'] == 'HEAVY':
                        conf = all(x_ >= 1e-2 for x_ in om2)
                pf['S6_3_full_system_confirms_class'] = conf
                t6 = t6_class(c0_ok and resolved, outc, pf['eta_class'], k0['cls'], conf)
                pf['T6_class'] = t6
                # S7.2, S7.4
                try:
                    s7n = s7_numeric(rec, eta_str, DP, fv)
                except (ZeroDivisionError, mp.libmp.NoConvergence, ValueError) as exc_:
                    s7n = dict(error='%s: %s' % (type(exc_).__name__, exc_), reason='S7 numeric step failed locally; T7 INCONCLUSIVE')
                pf['S7_numeric'] = s7n
                f2 = {'LIGHT': 'F2(ii) light: not eliminable', 'TACHYONIC': "F2(ii') maximum of V_eff: elimination procedure not admissible in this regime (not a proof that no UV origin exists)",
                      'HEAVY': 'heavy: eliminable at tree level', 'BOUNDARY': 'BOUNDARY: not classified'}[pf['eta_class']]
                pf['F2'] = f2
                t7 = t7_class(c0_ok and resolved, pf['eta_class'], s7n, res['S7_3'])
                pf['T7_class'] = t7
                rr['per_f'][tag] = pf
                note("%s root %d f=%s: eta %s, m2 %s %s, C0 %s, C1 %s, C2 %s, C3 %s -> %s; T6 %s; T7 %s; DC1 %s; X2 maxdiff %s"
                     % (name, ri, tag, rr['eta_bar'][:12], pf['m2_eta_eff'][:12], pf['eta_class'], c0_ok, c1, c2['status'],
                        c3, outc, t6, t7, v_, s7n.get('X2_full_vs_eliminated_maxdiff', 'n/a (' + str(s7n.get('reason')) + ')')))
            rr['F2_iii_f_to_0'] = dict(implicit_function_condition_Veff2_nonzero=bool(V2_ != 0), regular=bool(V2_ != 0),
                                       effective_operator=("W(dphi, dpsi) = stationary value over eta of J_tot^2/(2M^2) - Lam/4 (n.dphi)^4:"
                                                           " homogeneous of degree 4 in the gradients, coefficients are functions of the"
                                                           " gradient direction (coefficient functions, not constants)"))
            row['per_root'].append(rr)
        row['DC4'] = bool(dc4_bg)
        dc4_ok = dc4_ok and dc4_bg
        if en['n_located'] > 0:
            row['T6'] = 'eta_bar exists (%d root(s))' % en['n_located']
        rows.append(row)
    checks['DC1'] = 'PASS' if all(v_ == 'PASS' for v_ in dc1_all) else ('FAIL' if 'FAIL' in dc1_all else 'INCONCLUSIVE')
    checks['DC4'] = bool(dc4_ok)
    checks['DC5'] = bool(dc5_ok)

    # ---------------- summaries ----------------
    summ = {}
    n_with = sum(1 for r_ in rows if r_['per_root'])
    for fv in F_VALUES:
        tag = str(fv)
        cnt = {}
        for key in ('T3fund_outcome', 'T6_class', 'T7_class', 'eta_class'):
            cc = {}
            for r_ in rows:
                if not r_['per_root']:
                    lab = 'INCONCLUSIVE (no eta_bar)'
                    cc[lab] = cc.get(lab, 0) + 1
                for rr in r_['per_root']:
                    lab = rr['per_f'][tag][key]
                    cc[lab] = cc.get(lab, 0) + 1
            cnt[key] = cc
        t7v = [rr['per_f'][tag]['T7_class'] for r_ in rows for rr in r_['per_root']]
        cnt['T7_set'] = ('REPRODUCED-RANK-ONE' if t7v and all(x_ == 'REPRODUCED-RANK-ONE' for x_ in t7v) and n_with == len(rows)
                         else ('NOT-REPRODUCED' if t7v and all(x_ == 'NOT-REPRODUCED' for x_ in t7v) else 'CONDITIONAL / mixed (see counts)'))
        summ[tag] = cnt
    res['T6_set'] = dict(backgrounds_with_eta_bar=n_with, backgrounds_without=len(rows) - n_with,
                         set_class=('PARTIAL' if 0 < n_with < len(rows) else ('see per-background classes' if n_with == len(rows) else 'INCONCLUSIVE')),
                         subset=[r_['name'] for r_ in rows if r_['per_root']], complement=[r_['name'] for r_ in rows if not r_['per_root']],
                         item3_foliation=("1+1D: local integrability automatic (Frobenius condition vanishes identically); constant eta_bar:"
                                          " leaves are the straight lines orthogonal to u_bar; non-constant eta out of scope"),
                         item1_item2=("eta_bar is fixed by the matter gradients and Lambda through dV/deta = 0 (S1); property of the"
                                      " construction, not an interpretation of nature"))
    res['historical_note_DC5'] = dict(backgrounds_with_Jprime_root=sum(1 for r_ in rows if r_['DC5']['historical_U1_U4_has_stationary']),
                                      recorded_U1_U4="17/41", note="J = 0 roots also solve sigma_bar J' = 0; listed, not a correction")
    checks['DC6'] = bool(dc6_ok)
    res['backgrounds'] = rows
    res['summary_per_f'] = summ
    DC = ['DC0', 'DC1', 'DC2', 'DC3', 'DC4', 'DC5', 'DC6']
    res['checks'] = {k_: (checks[k_] if isinstance(checks[k_], str) else ('PASS' if checks[k_] else 'FAIL')) for k_ in DC}
    res['ALL_DC_PASS'] = all(res['checks'][k_] == 'PASS' for k_ in DC)
    res['scope'] = ("ALL_DC_PASS covers DC0-DC6 only; not a physics verdict, not an independent proof. T3-fund, T5, T6, T7 are reported"
                    " separately per background set and per f; no single verdict. C4 not evaluated: 'passes T3-fund' cannot be stated.")
    res['seconds'] = round(time.time() - T0, 1)
    txt = json.dumps(res, indent=1, default=str) + "\n"
    import re
    if re.search(r'NaN|Infinity|\bnan\b|\binf\b', txt):
        res['checks']['DC6'] = 'FAIL'
        res['ALL_DC_PASS'] = False
        txt = json.dumps(res, indent=1, default=str) + "\n"
    (out / "uv_dynamic_eta_results.json").write_text(txt, encoding="utf-8", newline="\n")
    for k_ in DC:
        note("%s %s" % (k_, res['checks'][k_]))
    note("summary per f: %s" % json.dumps(summ))
    note("ALL_DC_PASS" if res['ALL_DC_PASS'] else "DC_FAILED_OR_INCONCLUSIVE")
    (out / "uv_dynamic_eta_output.txt").write_text("\n".join(log) + "\n", encoding="utf-8", newline="\n")
    return 0 if res['ALL_DC_PASS'] else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(64)
    sys.exit(main(sys.argv[1]))
