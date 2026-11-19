"""Post-hoc, read only: compares the Q_new(c;1) classes of this run with the REGISTERED T3-fund Proca Stage 1 result
for the same 41 PE backgrounds (vector mediator without the timelike current, gamma_u = 0; registered commit 34953d0,
field uv.strictly_hyperbolic). No recomputation, no new criterion. Note: the registered result is the registered
T3-fund evaluation of its own Q(c;1); this is a comparison of two registered/authorised outcomes, not a polynomial identity.
Usage (from the repository root): python3 -I qnew_crosstab_registered_t3fund.py <qnew json> <out csv>"""
import collections, csv, json, sys
reg = json.load(open('t3fund_proca_stage1_result_bundle/proca_stage1_result.json', encoding='utf-8'))
pe = {(x['case']['mode'], x['case']['label']): x for x in reg if x['case']['pop'] == 'PE'}
q = json.load(open(sys.argv[1], encoding='utf-8'))
def key(tag):
    if tag == 'ANCHOR_B3':
        return ('rank1_minus', 'B3_anchor')
    parts = tag.split('_')
    return (parts[0] + '_' + parts[1], 'E3_' + '_'.join(parts[2:]))
rows, tab = [], {s: collections.Counter() for s in ('1', '-1')}
for r in q['backgrounds']:
    x = pe[key(r['tag'])]
    g0 = x['uv']['strictly_hyperbolic']
    row = dict(tag=r['tag'], registered_label=key(r['tag'])[1], registered_t3fund_outcome=x['outcome'],
               registered_t3fund_strictly_hyperbolic=g0)
    for s in ('1', '-1'):
        row['qnew_class_sigma%s' % ('p' if s == '1' else 'm')] = r[s]['class_primary_spec_relative']
        tab[s][(g0, r[s]['class_primary_spec_relative'])] += 1
    rows.append(row)
assert len(rows) == 41 and len({rw['registered_label'] + rw['tag'][:11] for rw in rows}) == 41
with open(sys.argv[2], 'w', newline='\n', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n'); w.writeheader(); w.writerows(rows)
for s in tab:
    print('sigma %s: (registered strictly hyperbolic at gamma_u = 0, Q_new class) ->' % s, dict(tab[s]))
print('COMPLEX only with the current (registered strictly hyperbolic at gamma_u = 0):')
for rw in rows:
    if rw['registered_t3fund_strictly_hyperbolic'] and 'COMPLEX' in (rw['qnew_class_sigmap'] + rw['qnew_class_sigmam']):
        print('  %s  sigma=+1: %s  sigma=-1: %s' % (rw['tag'], rw['qnew_class_sigmap'], rw['qnew_class_sigmam']))
print('SIX-SIMPLE-REAL only with the current (registered not strictly hyperbolic at gamma_u = 0):')
for rw in rows:
    if not rw['registered_t3fund_strictly_hyperbolic'] and 'SIX' in (rw['qnew_class_sigmap'] + rw['qnew_class_sigmam']):
        print('  %s' % rw['tag'])
print('registered outcomes of the 13 not strictly hyperbolic:', collections.Counter(rw['registered_t3fund_outcome'] for rw in rows if not rw['registered_t3fund_strictly_hyperbolic']))
print('classes (cls) of the 13:', collections.Counter('_'.join(rw['tag'].split('_')[2:-1]) for rw in rows if not rw['registered_t3fund_strictly_hyperbolic']))
