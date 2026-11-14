"""Aggregate the registered Proca Stage 1 evaluation from the raw JSONL, with a UNIQUE case key (pop, cut, mode, label).
Reason: the driver t3fund_proca_stage1.py (commit 0649531) keys cases by label only; the label 'E3_stable_6052' occurs twice in the manifest (rank1_minus and rank1_plus), so the aggregated JSON written by the driver contains a duplicate and misses one case.
The raw JSONL (one line per evaluated case, both collided cases present) is the source of truth. Usage: python aggregate_from_jsonl.py MANIFEST JSONL OUT_JSON OUT_CELL_TABLE OUT_STATS"""
import sys, json, collections as C, numpy as np
man = json.load(open(sys.argv[1])); rows = [json.loads(l) for l in open(sys.argv[2])]
key = lambda c: (c['pop'], c.get('cut'), c.get('mode'), c['label'])
ids = [key(e) for k in ('PE', 'PN', 'H3') for e in man[k]]
assert len(ids) == len(set(ids)) == 10081
d = {}
for r in rows:
    k = key(r['case']); assert k not in d, k; d[k] = r
assert set(d) == set(ids) and len(rows) == 10081
R = [d[i] for i in ids]                      # manifest order
json.dump(R, open(sys.argv[3], 'w'), indent=1)
G = C.OrderedDict()
for r in R: G.setdefault((r['case']['pop'], r['case'].get('cut'), r['case'].get('mode')), []).append(r)
tab = '| population | cut | mode | N | PASS | FAIL-UNSTABLE | FAIL-GHOST | INCONCLUSIVE |\n|---|---|---|---|---|---|---|---|\n'
for (p, cu, m), g in G.items():
    o = C.Counter(r['outcome'] for r in g)
    tab += f"| {p} | {cu if cu is not None else '-'} | {m} | {len(g)} | {o['PASS']} | {o['FAIL-UNSTABLE']} | {o['FAIL-GHOST']} | {o['INCONCLUSIVE']} |\n"
open(sys.argv[4], 'w').write(tab)
nh = [r for r in R if r['case']['pop'] != 'H3']; S = {}
S['identities'] = len(d); S['outcomes'] = dict(C.Counter(r['outcome'] for r in R))
S['PE'] = dict(C.Counter(r['outcome'] for r in R if r['case']['pop'] == 'PE'))
S['PE_by_mode'] = {m: dict(C.Counter(r['outcome'] for r in R if r['case']['pop'] == 'PE' and r['case']['mode'] == m)) for m in ('rank1_minus', 'rank1_plus')}
S['PN'] = dict(C.Counter(r['outcome'] for r in R if r['case']['pop'] == 'PN'))
S['C0_max'] = max(max(abs(v) for v in r['C0'].values()) for r in R)
k0 = [r for r in R if r.get('k0')]; S['k0'] = [len(k0), sum(r['k0']['agree'] for r in k0), max(r['k0']['max_rel_dev'] for r in k0)]
c1 = [r for r in R if r.get('C1b')]; S['C1b'] = [len(c1), sum(r['C1b']['ok'] for r in c1), max(r['C1b']['max_rel_dev_M_vs_Q'] for r in c1)]
h = [r for r in R if r.get('hessian_route')]; S['hessian'] = [len(h), sum(r['hessian_route']['ok'] for r in h), max(r['hessian_route']['max_rel_ratio_dev'] for r in h)]
m = [r for r in R if r.get('mu')]; S['mu'] = dict(n=len(m), defined=sum(bool(r['mu'].get('defined')) for r in m), violations=sum(r['mu'].get('violations', 0) for r in m if r['mu'].get('defined')))
p4 = [r for r in R if r.get('pred4')]; S['pred4'] = {f"{r_o}|E={e}": n for (r_o, e), n in C.Counter((r['outcome'], r['pred4']['E']) for r in p4).items()}
p3 = [r['uv']['pred3_max_err'] for r in nh if 'pred3_max_err' in r['uv']]; S['pred3'] = [len(p3), max(p3)]
S['hyper'] = dict(C.Counter(r['uv']['strictly_hyperbolic'] for r in nh)); S['kappa_kin'] = dict(C.Counter(str(r['kappa_kin']) for r in nh))
hi = [r for r in nh if r['pred6']['F_max_imag'] > 1e-3]; S['pred6'] = [len(hi), sum(r['pred6']['first_point_complex'] and r['outcome'] == 'FAIL-UNSTABLE' for r in hi)]
fl = [r for r in nh if r['outcome'] == 'FAIL-UNSTABLE']
S['fail'] = dict(n=len(fl), first_complex=sum(r['pred6']['first_point_complex'] for r in fl), nonhyper=sum(not r['uv']['strictly_hyperbolic'] for r in fl), ncomplex=[min(r['n_complex'] for r in fl), max(r['n_complex'] for r in fl)], maximag=[min(r['uv']['max_imag'] for r in fl), max(r['uv']['max_imag'] for r in fl)])
ps = [r for r in nh if r['outcome'] == 'PASS']; S['pass_n_real'] = dict(C.Counter(r['n_real'] for r in ps)); S['pass_hyper'] = dict(C.Counter(r['uv']['strictly_hyperbolic'] for r in ps))
S['diag'] = {k: sum(r.get(k, 0) for r in nh) for k in ('hidden_by_double', 'real_with_tiny_imag', 'n_nodata', 'unresolved_stuck', 'n_unresolved')}
pn = [r for r in nh if r['case']['pop'] == 'PN']; f = [r['uv']['front'] for r in pn if r['outcome'] == 'PASS']
S['PN_PASS_front'] = [min(f), float(np.median(f)), max(f), sum(x <= 2.0 for x in f), len(f)]
pe = [r for r in nh if r['case']['pop'] == 'PE' and r['outcome'] == 'PASS']; S['PE_PASS_front_le2'] = sum(r['uv']['front'] <= 2.0 for r in pe)
S['C4_PN'] = [r['case']['label'] for r in pn if r['outcome'] == 'PASS' and r['uv']['front'] <= 2.0][:2]
S['C4_PE'] = [r['case']['label'] + '/' + r['case']['mode'] for r in nh if r['case']['pop'] == 'PE' and r['outcome'] == 'PASS' and r['uv']['front'] <= 2.0][:1]
S['PN95'] = {mm: [sum(r['outcome'] != 'PASS' for r in pn if r['case']['cut'] == 0.95 and r['case']['mode'] == mm), 1000] for mm in ('rank1_minus', 'rank1_plus')}
json.dump(S, open(sys.argv[5], 'w'), indent=1); print(json.dumps(S, indent=1))
