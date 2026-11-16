"""Unit tests for uv_dynamic_eta_check.py revision 1.1 (review of 2026-11-16, items 4-8).
Synthetic inputs only; no registered background is evaluated (the registered data files are only hashed in T8).
Usage: python3 -I test_script_a_rev1_1.py <script_dir> <repo_root>"""
import inspect
import os
import shutil
import sys
import tempfile

sys.path.insert(0, sys.argv[1])
REPO = sys.argv[2]
import mpmath as mp
import sympy as sp
import uv_dynamic_eta_check as U

results = []


def check(name, cond, detail=''):
    results.append((name, bool(cond)))
    print('%-4s %s %s' % ('OK' if cond else 'FAIL', name, detail))


REC = {}
for k_, v_ in dict(p=0.3, s=0.2**2/0.3, kappa=-0.2, Lambda=0.5, w=0.4, q=0.9, v=-0.3, r=0.7).items():
    REC[k_] = sp.Rational(str(v_))
    REC[k_ + '_str'] = str(v_)
ETA = '0.39'

# ---- 4. P2a keeps the multiplicity of omega = 0 --------------------------------------------------------------------
om, k = U.om, U.k
full = sp.expand(om**4*(om**2 - 1)*(om**2 - 4)*(om**2 - k**2 - 9))
fl = sp.factor_list(full, om, k)
sqf = sp.expand(sp.Mul(*[a_ for a_, m_ in fl[1] if sp.Poly(a_, om, k).total_degree() > 0]))
r_full = U.p2a(U.DetPoly(full), REC, ETA, 1)
r_sqf = U.p2a(U.DetPoly(sqf), REC, ETA, 1)
check('T4a P2a on full determinant: multiplicity 4 -> ZERO MODE', r_full['zero_root_multiplicity'] == 4 and r_full['cls'] == 'ZERO MODE', r_full['zero_root_multiplicity'])
check('T4b square-free product would lose it (multiplicity 1 -> INCONCLUSIVE)', r_sqf['zero_root_multiplicity'] == 1 and r_sqf['cls'] == 'INCONCLUSIVE')
src = inspect.getsource(U.main)
check('T4c main calls P2a with the full determinant', 'p2a(DP_full, rec, eta_str, fv)' in src and 'p2a(DP, ' not in src)

# ---- 5. T6 needs conf is True ----------------------------------------------------------------------------------------
P = 'PASS (C0-C3; C4 not evaluated)'
check('T5a conf True -> STATED-DYNAMICAL', U.t6_class(True, P, 'HEAVY', 'ZERO MODE', True) == 'STATED-DYNAMICAL')
check('T5b conf False -> not STATED-DYNAMICAL', U.t6_class(True, P, 'HEAVY', 'ZERO MODE', False) == 'INCONCLUSIVE')
check('T5c conf None -> not STATED-DYNAMICAL', U.t6_class(True, P, 'LIGHT', 'ZERO MODE', None) == 'INCONCLUSIVE')
check('T5d unresolved root -> INCONCLUSIVE', U.t6_class(False, P, 'HEAVY', 'ZERO MODE', True) == 'INCONCLUSIVE')
check('T5e main uses t6_class with conf', 't6_class(c0_ok and resolved, outc, pf[\'eta_class\'], k0[\'cls\'], conf)' in src)

# ---- 6. C1 analytic only with block-diagonal principal symbol ----------------------------------------------------------
SYM = U.build_symbolic()
off0, blocks, _ = U.principal_symbolic(SYM['Ku'])
with mp.workdps(50):
    vals = U.numeric_bg(REC, ETA, 1, 'u')
    c2_ok = U.c2_numeric(blocks, vals, True)
    c2_bad = U.c2_numeric(blocks, vals, False)
check('T6a valid block diagonality -> C1 PASS [ANALYTIC]', c2_ok['status'] == 'PASS' and U.c1_from_c2(c2_ok) == 'PASS [ANALYTIC]')
check('T6b off-diagonal principal parts present -> no analytic PASS', c2_bad['status'] == 'INCONCLUSIVE' and U.c1_from_c2(c2_bad) == 'INCONCLUSIVE', c2_bad.get('reason'))
check('T6c main passes off0 to c2_numeric', "c2_numeric(pblocks, numeric_bg(rec, eta_str, fv, 'u'), off0)" in src)

# ---- 7a. singular F_eta_eta ------------------------------------------------------------------------------------------
DP = U.DetPoly(SYM['det_u'])
orig = mp.diff
state = {'n': 0}


def fake_diff(*a_, **kw_):
    state['n'] += 1
    if state['n'] == 1:            # the first call in s7_numeric is F_eta_eta
        return mp.mpf(0)
    return orig(*a_, **kw_)


mp.diff = fake_diff
try:
    s7n = U.s7_numeric(REC, ETA, DP, 1)
finally:
    mp.diff = orig
S73 = dict(operator_identity_with_B4_rank_one=True, characteristic_polynomial_identical=True)
check('T7a singular F_eta_eta handled locally', s7n.get('singular') is True)
check('T7b singular -> T7 INCONCLUSIVE', U.t7_class(True, 'HEAVY', s7n, S73) == 'INCONCLUSIVE')
check('T7c local failure record -> T7 INCONCLUSIVE', U.t7_class(True, 'HEAVY', dict(error='NoConvergence'), S73) == 'INCONCLUSIVE')
check('T7d main catches local S7 failures', 'except (ZeroDivisionError, mp.libmp.NoConvergence, ValueError)' in src)

# ---- 7b. DC5 with overlapping classes ----------------------------------------------------------------------------------
tol = mp.mpf('1e-30')
a, b = mp.mpf('0.1'), mp.mpf('0.2')
check('T7e overlap J\'=J root, one root found -> PASS (no double counting)', U.dc5_match([a], [a], [a], tol)['pass_'])
check('T7f overlap, extra root found -> FAIL', not U.dc5_match([a, b], [a], [a], tol)['pass_'])
check('T7g disjoint classes both found -> PASS', U.dc5_match([a, b], [a], [b], tol)['pass_'])
check('T7h missing J\' root -> FAIL', not U.dc5_match([b], [a], [b], tol)['pass_'])

# ---- 7c. multiple and close roots --------------------------------------------------------------------------------------
z = U.zz
for label, Pz, expect in (('double root z=2', sp.Poly((z - 2)**2*(z - 3)*(z + 1), z), {2: False, 3: True}),
                          ('close roots', sp.Poly((z - 2)*(z - 2 - sp.Rational(1, 10**25))*(z - 5), z), {2: False, 5: True})):
    saved = U.background_functions
    U.background_functions = lambda rec, Lv=None, P_=Pz: dict(P=P_, stripped_z_power=0, a=0, b=0)
    try:
        en, _ = U.enumerate_roots(REC)
    finally:
        U.background_functions = saved
    got = {int(round(float(r_['z']))): r_['resolved'] for r_ in en['roots']}
    ok = all(any(got_k == k_ and got_v == v_ for got_k, got_v in got.items()) for k_, v_ in expect.items())
    if label == 'double root z=2':
        ok = ok and [r_['multiplicity'] for r_ in en['roots']] == [2, 1] and en['sturm_count'] == en['n_located'] == 2
    check('T7i ' + label + ' -> unresolved, INCONCLUSIVE path', ok, [(r_['z'][:8], r_['multiplicity'], r_['resolved']) for r_ in en['roots']])
check('T7j main makes unresolved roots INCONCLUSIVE', "if not (c0_ok and resolved):" in src)

# ---- 8. changed input files -> ABORT -------------------------------------------------------------------------------------
tmp = tempfile.mkdtemp()
b3 = os.path.join(tmp, 'b3.json')
mn = os.path.join(tmp, 'man.json')
shutil.copy(os.path.join(REPO, U.DATA_B3), b3)
shutil.copy(os.path.join(REPO, U.DATA_MAN), mn)
check('T8a unchanged copies pass the hash check', U.verify_inputs({U.DATA_B3: b3, U.DATA_MAN: mn}))
for victim in (b3, mn):
    shutil.copy(os.path.join(REPO, U.DATA_B3), b3)
    shutil.copy(os.path.join(REPO, U.DATA_MAN), mn)
    with open(victim, 'ab') as fh:
        fh.write(b' ')
    try:
        U.load_backgrounds({U.DATA_B3: b3, U.DATA_MAN: mn})
        check('T8 modified %s -> ABORT' % os.path.basename(victim), False)
    except SystemExit as exc:
        check('T8 modified %s -> ABORT' % os.path.basename(victim), 'sha256 mismatch' in str(exc), str(exc)[:60])
check('T8d main verifies inputs before any computation', src.index('verify_inputs()') < src.index('s0_checks()'))
shutil.rmtree(tmp)

n_fail = sum(1 for _, ok in results if not ok)
print('TOTAL %d tests, %d failed' % (len(results), n_fail))
sys.exit(1 if n_fail else 0)
