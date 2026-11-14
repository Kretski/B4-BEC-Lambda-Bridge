"""C4 outcome per designated background, rule of T3_FUND_PROCA_SPEC.md section 6 incl. Clarification 1 (as registered):
 PASS: norm ratio <= 3 at both resolutions, total-energy drift <= 1e-9 (coupled run, both resolutions), free control reaches the same target (both resolutions), 1% agreement of the norm ratio between resolutions.
 FAIL-UNSTABLE: norm ratio > 3 with agreement between resolutions (never FAIL on energy alone).
 INCONCLUSIVE: control misses target / agreement fails / drift missed with control passing / non-finite values.
Usage: python aggregate_c4.py DIR TAG   (reads DIR/c4proca_TAG_{coupled,control}_N{16384,32768}.json)"""
import json, sys
TARGET = 1e-9; RMAX = 3.0; AGREE = 0.01
def outcome(c1, c2, k1, k2):
    reasons = []
    if not all(r['finite'] for r in (c1, c2, k1, k2)): return 'INCONCLUSIVE', ['non-finite values']
    agree = abs(c1['norm_ratio_max'] - c2['norm_ratio_max']) <= AGREE*abs(c2['norm_ratio_max'])
    ratio_hi = c1['norm_ratio_max'] > RMAX and c2['norm_ratio_max'] > RMAX
    ratio_ok = c1['norm_ratio_max'] <= RMAX and c2['norm_ratio_max'] <= RMAX
    ctrl = k1['drift_total_max'] <= TARGET and k2['drift_total_max'] <= TARGET
    drift = c1['drift_total_max'] <= TARGET and c2['drift_total_max'] <= TARGET
    if not ctrl: return 'INCONCLUSIVE', ['control misses energy target (P7/Clarification: INCONCLUSIVE, not a physical FAIL)']
    if ratio_hi and agree: return 'FAIL-UNSTABLE', ['ratio>3 at both resolutions, agreement']
    if not agree: reasons.append('1% agreement of norm ratio fails')
    if not ratio_ok and not ratio_hi: reasons.append('ratio >3 at one resolution only')
    if not ctrl: reasons.append('control misses energy target')
    if not drift: reasons.append('coupled drift misses target' + (' (control passes)' if ctrl else ''))
    if not reasons and ratio_ok and agree and ctrl and drift: return 'PASS', []
    return 'INCONCLUSIVE', reasons
if __name__ == '__main__':
    d, tag = sys.argv[1], sys.argv[2]; L = lambda k, n: json.load(open('%s/c4proca_%s_%s_N%d.json' % (d, tag, k, n)))
    c1, c2, k1, k2 = L('coupled', 16384), L('coupled', 32768), L('control', 16384), L('control', 32768)
    o, why = outcome(c1, c2, k1, k2)
    print(json.dumps(dict(tag=tag, outcome=o, reasons=why, ratio=[c1['norm_ratio_max'], c2['norm_ratio_max']], drift_total=[c1['drift_total_max'], c2['drift_total_max']],
        drift_fl=[c1['drift_fl_max'], c2['drift_fl_max']], control_drift=[k1['drift_total_max'], k2['drift_total_max']], control_ratio=[k1['norm_ratio_max'], k2['norm_ratio_max']],
        ratio_agreement_rel=abs(c1['norm_ratio_max']-c2['norm_ratio_max'])/c2['norm_ratio_max']), indent=1))
