"""Self-test on random backgrounds with complex bands (not manifest cases): the FAIL path (P6, verdicts) and the Pred-4 contrapositive; plus the ghost control path."""
import sys, json, numpy as np, mpmath as mp
sys.path.insert(0, '.')
import t3fund_proca_stage1 as S
rng = np.random.default_rng(20261116); found = []; tried = 0
while len(found) < 4 and tried < 400:
    tried += 1
    gp, gs = [float(z) for z in rng.uniform(-1.5, 1.5, 2)]; Lam = float(rng.uniform(0, 1)); w, q, v, r = [float(z) for z in rng.uniform(-3, 3, 4)]
    case = dict(pop='TEST', scope=True, label='f%d' % tried, p=gp**2, s=gs**2, kappa=-gp*gs, Lambda=Lam, w=w, q=q, v=v, r=r, gphi=gp, gpsi=gs)
    o = S.evaluate(case)
    if o['outcome'] != 'PASS': found.append((case, o))
for case, o in found:
    print(case['label'], o['outcome'], 'n_complex', o['n_complex'], 'kmin', o['complex_at_kmin'], 'uvband', o['complex_in_uv_band'], 'P6', o.get('P6'), 'E', o['pred4']['E'], 'F_max_imag', o['pred6']['F_max_imag'], 'hyper', o['uv']['strictly_hyperbolic'], 'mu', o['mu'].get('defined'), 'hidden', o['hidden_by_double'], 'tiny', o['real_with_tiny_imag'], 'P6complete', o['P6']['complete'], flush=True)
    assert not (o['pred4']['E'] and o['outcome'] != 'PASS'), 'Pred-4 counterexample'
print('tried', tried, 'non-PASS found', len(found))
g = dict(found[0][0] if found else case, pop='H3', kinetic_sign_psi=-1.0, scope=False); print('ghost control:', S.evaluate(g)['outcome'])

# ---- unit tests of the corrected logic (monkeypatching; random background, not a manifest case) ----
case, o0 = found[0]
orig_mp = S.roots_Q_mp; orig_dbl = S.roots_Q_dbl
# (1) hidden-by-double: double precision reports purely real roots while the 50-digit roots are complex -> the point must NOT be classified 'real'
S.roots_Q_dbl = lambda P, x: np.array([float(i) for i in range(6)], dtype=complex)
o1 = S.evaluate(dict(case, label='hid')); S.roots_Q_dbl = orig_dbl
print('hidden-by-double test: hidden_by_double =', o1['hidden_by_double'], ' n_complex =', o1['n_complex'], ' unresolved =', o1['unresolved_stuck'], ' outcome =', o1['outcome'])
assert o1['n_grid'] == 92 and o1['hidden_by_double'] == 92 and o1['n_real'] == 0 and o1['n_complex'] == 92, 'hidden-by-double must hit all 92 points and leave no real point'
# (2) P6 completeness: a refinement point that cannot be evaluated at 100 digits -> INCONCLUSIVE (never PASS / never unqualified FAIL)
def mp_fail_100(P, x, dps):
    return None if dps == 100 else orig_mp(P, x, dps)
S.roots_Q_mp = mp_fail_100; o2 = S.evaluate(dict(case, label='p6')); S.roots_Q_mp = orig_mp
print('P6 completeness test: complete =', o2['P6']['complete'], ' n_rows =', o2['P6']['n_rows'], ' n_refined =', o2['P6']['n_refined'], ' C3 =', o2['C3_numeric'], ' outcome =', o2['outcome'])
assert (not o2['P6']['complete']) and o2['C3_numeric'] == 'INCONCLUSIVE'
# (3) unified classification by r_m: double screening says complex (r_d large) while the 50-digit value is <=1e-7 -> 'real' (no 'complex' from the r_d branch)
def dbl_complex(P, x): return np.array([1+1j, 1-1j, 2, 3, 4, 5], dtype=complex)
S.roots_Q_dbl = dbl_complex; o3 = S.evaluate(dict(found[0][0], label='cls')); S.roots_Q_dbl = orig_dbl
# on a REAL background (50-digit r_m<=1e-7) with double forced complex: all points must be real
case_real = dict(pop='TEST', scope=False, label='real', p=0.09, s=0.04, kappa=-0.06, Lambda=0.3, w=0.2, q=0.1, v=-0.1, r=0.15, gphi=0.3, gpsi=0.2)
S.roots_Q_dbl = dbl_complex; o4 = S.evaluate(case_real); S.roots_Q_dbl = orig_dbl
print('classification test on a real background with double forced complex: n_real =', o4['n_real'], ' n_complex =', o4['n_complex'], ' outcome =', o4['outcome'])
assert o4['n_real'] == 92 and o4['n_complex'] == 0 and o4['outcome'] == 'PASS'
# (4) nodata cannot be masked by a P6 PASS: complex window with dps-dependent imaginary part (P6 -> PASS) + one unevaluated point -> INCONCLUSIVE
KG = S.KGRID; win = (KG[9]*0.999, KG[11]*1.001); nod = KG[60]
def mock(P, x, dps):
    if abs(x-nod) < 1e-12*nod and dps == 50: return None
    rr = orig_mp(P, x, dps)
    if rr is None: return None
    if win[0] <= x <= win[1] and dps == 50:
        rr = list(rr); a = mp.re(rr[0]); rr[0] = mp.mpc(a, mp.mpf('1e-3')); rr[1] = mp.mpc(a, -mp.mpf('1e-3'))
    return rr
S.roots_Q_mp = mock; o5 = S.evaluate(dict(case_real, label='nod')); S.roots_Q_mp = orig_mp
print('nodata-masking test: n_nodata =', o5['n_nodata'], ' n_complex =', o5['n_complex'], ' P6 =', o5.get('P6'), ' C3 =', o5['C3_numeric'], ' outcome =', o5['outcome'])
assert o5['n_nodata'] == 1 and o5['n_complex'] > 0 and o5['P6']['max_consecutive_gt10u'] == 0 and o5['P6']['complete'] and o5['C3_numeric'] == 'INCONCLUSIVE'
print('unit tests passed')
