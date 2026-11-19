"""Post-hoc tables from proca_timelike_stage2.json (read only; no recomputation, no new criterion).
Usage: python3 -I stage2_summary_tables.py <proca_timelike_stage2.json> <out.csv>"""
import csv, json, sys
import mpmath as mp
mp.mp.dps = 30
d = json.load(open(sys.argv[1], encoding='utf-8'))
rows = d['backgrounds']
def lo(x): return mp.mpf(x[0])
def hi(x): return mp.mpf(x[1])
def minabs(encs):  # smallest distance of an enclosure from 0 (enclosures here exclude 0)
    return min((min(abs(lo(e)), abs(hi(e))) for e in encs), default=None)
out = []
for r in rows:
    row = dict(tag=r['tag'], n_real=r['n_real_branches'], squarefree=r['F_squarefree'], gcd_deg=r['gcd_F_tr_degree'],
               delta_rank1=r['delta_rank1'])
    for sname, v in r['D4']['sets'].items():
        row['D4_%s_T7a_margin' % sname] = '%.6g' % v['T7a_margin']
        row['D4_%s_T7a' % sname] = v['T7a_class']
        row['D4_%s_ext_margin' % sname] = '%.6g' % v['extended_margin']
        row['D4_%s_ext' % sname] = v['extended_class']
    row['c_Lambda_dominates'] = r['D4']['c_Lambda_dominates']
    for sg in ('1', '-1'):
        x = r['D3'][sg]; br = x['branches']; t = 'sp' if sg == '1' else 'sm'
        row[t + '_crossing'] = x['n_crossing']; row[t + '_nocross'] = x['n_no_crossing']
        row[t + '_RUNDEF'] = x['n_R_UNDEFINED']; row[t + '_notcert'] = x['n_not_certified']
        row[t + '_poles'] = sum(1 for b in br if b['pole'])
        row[t + '_kder_status'] = x['k_der_status']
        row[t + '_kder_over_M'] = mp.nstr(lo(x['k_der_over_M_enclosure']), 12) if x['k_der_over_M_enclosure'] else ''
        row[t + '_V2'] = x['V2_holds']
        row[t + '_V2_max_N_plus_c2'] = mp.nstr(max(hi(b['N_plus_c2_enclosure']) for b in br), 12)
        row[t + '_H3_all'] = x['H3_all_branches']
        row[t + '_min_dist0_S_minus_1'] = mp.nstr(min(abs(lo(b['N_plus_c2_enclosure']) - 1) for b in br), 6)
        row[t + '_min_dist0_V2_margin'] = mp.nstr(minabs([b['V2_branch_margin_enclosure'] for b in br]), 6)
        row[t + '_min_dist0_H3'] = mp.nstr(minabs([b['H3_enclosure'] for b in br]), 6)
    out.append(row)
with open(sys.argv[2], 'w', newline='\n', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]), lineterminator='\n'); w.writeheader(); w.writerows(out)
def g(k, f=float): return [f(r_[k]) for r_ in out]
print('max delta_rank1:', d['max_delta_rank1'])
for s in ('A_g1_rho0.3', 'B_g1_rho0.1', 'C_gStrong_rho0.3'):
    for kind in ('T7a', 'ext'):
        m = g('D4_%s_%s_margin' % (s, kind))
        print('D4 %s %s: min |margin| = %.4g (scale %s)' % (s, kind, min(abs(x) for x in m), s))
flips = [r_['tag'] for r_ in out if r_['D4_C_gStrong_rho0.3_T7a'] != r_['D4_C_gStrong_rho0.3_ext']]
print('set C T7a -> extended class changes:', flips)
print('c_Lambda dominates on', sum(1 for r_ in out if r_['c_Lambda_dominates']), 'backgrounds')
for t in ('sp', 'sm'):
    print(t, 'smallest separations from 0: S-1 %s, V2 margin %s, H3 %s' % (
        min(g(t + '_min_dist0_S_minus_1')), min(g(t + '_min_dist0_V2_margin')), min(g(t + '_min_dist0_H3'))))
    print(t, 'k_der/M range', min(g(t + '_kder_over_M')), max(g(t + '_kder_over_M')),
          'min max(N+c2)', min(g(t + '_V2_max_N_plus_c2')))
