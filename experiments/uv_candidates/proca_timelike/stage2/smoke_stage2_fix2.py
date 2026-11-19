"""Smoke test of the Stage 2 revision (certified decisions, rank-one ABORT) on SYNTHETIC inputs only.
No registered background is read. Usage: python3 -I smoke_stage2_fix2.py <stage2 script> <previous stage2 script>"""
import sys, importlib.util
import sympy as sp, mpmath as mp


def load_mod(path, name):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m


st2 = load_mod(sys.argv[1], "st2"); prev = load_mod(sys.argv[2], "prev")
c = st2.c
ok = []
def chk(name, cond): ok.append(bool(cond)); print(("OK  " if cond else "FAIL"), name)
R = sp.Rational
rec0 = dict(p=R(1, 2), kappa=R(1, 4), Lambda=R(1, 10), w=R(3, 10), q=R(1, 5), v=R(1, 10), r=R(2, 5), cls='synthetic', mode='smoke')
rec0['s'] = rec0['kappa']**2/rec0['p']
# ---- exact identification (unchanged logic, re-checked)
F = sp.Poly((c - 1)**2*(c + 2)*(c**2 - 3), c)
tr = sp.Poly((c**2 - 3)*(c - 1)*(c + 5), c)
br = st2.exact_branches(F, sp.Poly(sp.gcd(F, tr), c))
fl = {round(float(z), 6): u for z, _, _, u in br}
chk("E1 exact R-UNDEFINED at -sqrt3, 1, sqrt3; not at -2", fl[1.0] and fl[round(3**.5, 6)] and fl[round(-3**.5, 6)] and not fl[-2.0])
# ---- rank-one ABORT
r_ok = dict(rec0); d = st2.rank1_check('SYN', r_ok)
chk("E2 exact rank-one input: delta_rank1 = 0, no ABORT", d == 0)
r_f = dict(rec0); r_f['s'] = rec0['s']*(1 + R(1, 10**16))
chk("E3 float-level deviation 1e-16 (as in registered K4, max 2.2e-16): no ABORT", st2.rank1_check('SYN', r_f) < R(1, 10**9))
r_bad = dict(rec0); r_bad['s'] = rec0['s']*(1 + R(1, 10**6))
try:
    st2.rank1_check('SYN_BAD', r_bad); ab = False
except SystemExit as e_:
    ab = 'rank-one' in str(e_)
chk("E4 deviation 1e-6 > 1e-9: ABORT", ab)
def aborts(rec_):
    try:
        st2.rank1_check('SYN_EDGE', rec_); return False
    except SystemExit as e_:
        return 'rank-one' in str(e_)
chk("E4b s < 0 with kappa = 0: ABORT", aborts(dict(rec0, s=R(-1, 10), kappa=R(0))))
chk("E4c s = 0 with kappa != 0: ABORT", aborts(dict(rec0, s=R(0))))
chk("E4d s < 0 with kappa != 0: ABORT", aborts(dict(rec0, s=R(-1, 10))))
r00 = dict(rec0, s=R(0), kappa=R(0))
chk("E4e s = 0 and kappa = 0: no ABORT, delta_rank1 = 0", (not aborts(dict(r00))) and st2.rank1_check('SYN00', r00) == 0)
# ---- certified sign helper
mp.iv.dps = 30
chk("E5 ivsign: [-1e-40, 1e-40] -> None, [0,0] -> 0, [1e-40, 1] -> +1",
    st2.ivsign(mp.iv.mpf(['-1e-40', '1e-40'])) is None and st2.ivsign(mp.iv.mpf(0)) == 0 and st2.ivsign(mp.iv.mpf(['1e-40', '1'])) == 1)
# ---- pole decision, exact
Fp = sp.Poly((c - 1)*(c**2 - 2)*(c + sp.Rational(1, 2)), c)
dec = {}
for (a, b), m in Fp.intervals():
    lo, hi = st2.enclose_root(Fp, sp.Rational(a), sp.Rational(b), 30)
    dec[round(float((lo + hi)/2), 6)] = st2.pole_decision(Fp, lo, hi)
chk("E6 pole: c = 1 exactly -> no pole; |c| = sqrt2 -> pole; c = -1/2 -> no pole",
    dec[1.0] is False and dec[round(2**.5, 6)] is True and dec[round(-2**.5, 6)] is True and dec[-0.5] is False)
# ---- exact zero H3 (q = 0): certified True via exact [0, 0]
rq = dict(rec0); rq['q'] = R(0); rq['s'] = rq['kappa']**2/rq['p']; st2.rank1_check('SYNQ', rq)
outq = st2.evaluate('SYN_Q0', rq)
chk("E7 q = 0: h3 is exactly 0 -> H3 certified True on every branch, both signs",
    all(b['H3_holds'] is True and b['H3_enclosure'] == ['0.0', '0.0'] for sg in ('1', '-1') for b in outq['D3'][sg]['branches']))
# ---- full synthetic background: agreement with the previous (ff359c11) decisions, enclosures, widths
st2.rank1_check('SYN', rec0)
new = st2.evaluate('SYN_REAL', rec0); old = prev.evaluate('SYN_REAL', rec0)
for sg in ('1', '-1'):
    dn, do = new['D3'][sg], old['D3'][sg]
    same = len(dn['branches']) == len(do['branches'])
    for bn, bo in zip(dn['branches'], do['branches']):
        same &= bn['certified'] and bn['crossing'] == bo['crossing'] and bn['H3_holds'] == bo['H3_holds']
        same &= bn['V2_branch'] == (mp.mpf(bo['V2_branch_margin_2_minus_N_plus_c2']) > 0)
        same &= bn['pole'] == ('pole_k_over_M' in bo)
        same &= sp.N(sp.Rational(bn['root_enclosure'][1]) - sp.Rational(bn['root_enclosure'][0])) <= 10.0**-65
    same &= dn['V2_holds'] == do['V2_holds'] and dn['k_der_status'] == do['k_der_status'] and dn['H3_all_branches'] == do['H3_all_branches']
    if dn['k_der_status'] == 'defined':
        kd = mp.mpf(do['k_der_over_M']); enc = [mp.mpf(x) for x in dn['k_der_over_M_enclosure']]
        same &= enc[0] - mp.mpf(10)**-18 <= kd <= enc[1] + mp.mpf(10)**-18
    chk("E8[%s] certified decisions = previous decisions on all branches; refined enclosure width <= 1e-65; k_der inside enclosure "
        "(status %s, V2 %s)" % (sg, dn['k_der_status'], dn['V2_holds']), same)
    worst = max(mp.mpf(b['control_R_direct_vs_formula_rel']) for b in dn['branches'] if 'control_R_direct_vs_formula_rel' in b)
    chk("E9[%s] control direct vs formula rel = %s < 1e-40" % (sg, mp.nstr(worst, 3)), worst < mp.mpf(10)**-40)
# ---- propagation of a non-certified sign (forced): status and V2 must be undetermined
# force: make the crossing sign undecidable on every branch (logic test of the propagation only)
cb_orig = st2.certified_branch
def cb_forced(*a, **k):
    r_ = cb_orig(*a, **k); r_['crossing'] = None; r_['certified'] = False; r_.pop('_kd', None); return r_
st2.certified_branch = cb_forced
outf = st2.evaluate('SYN_FORCED', rec0)
chk("E10 forced uncertified crossing: k_der_status 'undetermined (sign not certified)', V2 None, no k_der value",
    all(outf['D3'][sg]['k_der_status'] == 'undetermined (sign not certified)' and outf['D3'][sg]['V2_holds'] is None
        and outf['D3'][sg]['k_der_over_M_enclosure'] is None and outf['D3'][sg]['n_not_certified'] > 0 for sg in ('1', '-1')))
st2.certified_branch = cb_orig
# ---- no real branches
st2.lead_kernel_polys = lambda rec: (sp.Poly(c**2 + 1, c), sp.Poly(c**2 + 2, c))
outn = st2.evaluate('SYN_NOREAL', rec0)
chk("E11 no real branches: V2 None, H3 None, status 'undetermined (no real branches)'",
    all(outn['D3'][sg]['V2_holds'] is None and outn['D3'][sg]['H3_all_branches'] is None
        and outn['D3'][sg]['k_der_status'] == 'undetermined (no real branches)' for sg in ('1', '-1')))
chk("E12 delta_rank1 and D4 margins present in output", 'delta_rank1' in new and all('T7a_margin' in v and 'extended_margin' in v for v in new['D4']['sets'].values()))
print("TOTAL %d checks, %d failed" % (len(ok), ok.count(False)))
sys.exit(0 if all(ok) else 1)
