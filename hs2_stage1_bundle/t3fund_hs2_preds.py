# Registered predictions Pred-1..4 of T3_FUND_HS2_SPEC.md (section 4), evaluated on hs2_stage1_result.json.
import json, sys, numpy as np
sys.path.insert(0, '.')
import t3fund_hs2_stage1 as S
sys.path.insert(0, '/tmp/hs2/t3fund_hs2_spec_bundle')
man = json.load(open('t3fund_hs2_manifest.json')); R = json.load(open(sys.argv[1]))
# F_disc from the manifest (same order as build())
disc = {}
for d in man['H1']+man['H3_control']: disc[(d['cs'], d['idx'], 'w')] = d
def key(tag, th): return (tag, tuple(round(x, 12) for x in th))
Fd = {}
for tagp, lst in (('H1_', man['H1']), ('H3_', man['H3_control'])):
    for d in lst: Fd[(tagp+d['cs'], tuple(round(x, 12) for x in (d['a'], d['b'], d['e'], d['w'], d['q'], d['v'], d['r'])))] = d['F_disc']
out = {}
# Pred-1/2 on H1
h1 = [r for r in R if r['set'].startswith('H1_')]
tab = {}; ex1 = []; ex2 = 0; uvband = 0
for r in h1:
    fd = Fd[(r['set'], tuple(round(x, 12) for x in r['theta']))]; unst = r['outcome'] == 'FAIL-UNSTABLE'; ir_unst = fd < 0
    t = tab.setdefault(r['set'], dict(n=0, PASS=0, FAIL=0, INCONCLUSIVE=0, other=0, match=0, exceptions=0))
    t['n'] += 1
    t['PASS' if r['outcome'] == 'PASS' else 'FAIL' if unst else 'INCONCLUSIVE' if r['outcome'] == 'INCONCLUSIVE' else 'other'] += 1
    if unst == ir_unst: t['match'] += 1
    else: t['exceptions'] += 1; ex1.append((r['set'], r['theta'], r['outcome'], fd))
    if unst and not r['complex_at_kmin']: ex2 += 1
    if r['complex_in_uv_band']: uvband += 1
out['Pred1_H1'] = tab; out['Pred1_exceptions'] = ex1[:10]
out['Pred2'] = dict(unstable_not_complex_at_kmin=ex2, H1_with_complex_in_uv_band=uvband)
h3 = [r for r in R if r['set'].startswith('H3_')]; h2 = [r for r in R if r['set'].startswith('H2_')]; an = [r for r in R if r['set'].startswith('ANCHOR_')]
out['Pred3'] = dict(H3=dict(n=len(h3), fail_unstable=sum(r['outcome'] == 'FAIL-UNSTABLE' for r in h3)), H2=dict(n=len(h2), fail_ghost=sum(r['outcome'] == 'FAIL-GHOST' for r in h2)),
                    anchors={r['set']: r['outcome'] for r in an})
# Pred-4: K1 vs one-field HS (g=2a, theta=(phi+psi)/sqrt2): Q roots = one-field quartic roots
worst = 0; n4 = 0
for r in [x for x in R if x['set'].startswith('H1_K1') or x['set'] == 'ANCHOR_B2h'][:60]:
    a, b, e, w, q, v, rr = r['theta']; g = 2*a; wt = (w+v)/np.sqrt(2); qt = (q+rr)/np.sqrt(2)
    for kv in (0.05, 0.5, 2.0, 30.0):
        eps = 1/kv**2; Z = 1+2*g*g*(wt**2-qt**2)
        d = Z*np.polymul([-1, 0, 1], [-1, 0, 1+eps]) - 4*g*g*eps*np.concatenate([[0, 0], np.polymul([wt, qt], [wt, qt])])
        r1 = np.roots(d); r2 = np.roots([complex(x) for x in S.Qpoly([float(t) for t in r['theta']], kv)][::-1])
        worst = max(worst, max(min(abs(x-y) for y in r1) for x in r2), max(min(abs(x-y) for y in r2) for x in r1)); n4 += 1
out['Pred4'] = dict(n_checked=n4, max_root_distance=float(worst))
# diagnostics
out['diag'] = dict(max_poly_vs_K=max(r.get('poly_vs_K_maxrelerr', 0) for r in R), max_fo_dist=max(r['first_order']['max_dist_poly_root_to_fo_eig'] for r in R if 'first_order' in r),
                   max_C0_residual=max(r['C0']['residual_chi_eq'] for r in R), hidden_by_double=sum(r.get('hidden_by_double', 0) for r in R), unresolved_stuck=sum(r.get('unresolved_stuck', 0) for r in R),
                   max_u=max((float(r['P6']['max_u_str']) for r in R if 'P6' in r), default=0), min_ratio=min((float(r['P6']['min_ratio_gamma_over_u']) for r in R if 'P6' in r and r['P6']['min_ratio_gamma_over_u']), default=None),
                   max_consec_ge3_all=all(r['P6']['max_consecutive_gt10u'] >= 3 for r in R if r.get('outcome') == 'FAIL-UNSTABLE' and 'P6' in r),
                   outcome_counts={k: sum(1 for r in R if r['outcome'] == k) for k in set(r['outcome'] for r in R)})
print(json.dumps(out, indent=1, default=str)); json.dump(out, open('hs2_preds.json', 'w'), indent=1, default=str)
