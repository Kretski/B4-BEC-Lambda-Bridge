"""Reproducible tests: (1) aggregator rule on synthetic inputs, (2) constant background is a fixed point of the C4 rhs, (3) short coupled run (T=25) is finite with energy drift small."""
import sys, json, numpy as np; sys.path.insert(0, '.')
from aggregate_c4 import outcome
import t3fund_proca_c4 as C
b = lambda r, d, f=True: dict(norm_ratio_max=r, drift_total_max=d, finite=f)
T = [((b(1.3, 1e-16), b(1.3, 1e-16), b(1, 1e-16), b(1, 1e-16)), 'PASS'),
     ((b(5, 1e-16), b(5.02, 1e-16), b(1, 1e-16), b(1, 1e-16)), 'FAIL-UNSTABLE'),
     ((b(5, 1e-7), b(5.02, 1e-7), b(1, 1e-16), b(1, 1e-16)), 'FAIL-UNSTABLE'),            # never FAIL on energy alone: ratio decides, coupled drift irrelevant
     ((b(5, 1e-16), b(5.02, 1e-16), b(1, 1e-6), b(1, 1e-6)), 'INCONCLUSIVE'),             # control misses target -> INCONCLUSIVE even with ratio>3
     ((b(1.3, 1e-16), b(1.4, 1e-16), b(1, 1e-16), b(1, 1e-16)), 'INCONCLUSIVE'),          # 1% agreement fails
     ((b(1.3, 1e-7), b(1.3, 1e-7), b(1, 1e-16), b(1, 1e-16)), 'INCONCLUSIVE'),            # drift missed with control passing
     ((b(1.3, 1e-16), b(1.3, 1e-16), b(1, 1e-6), b(1, 1e-6)), 'INCONCLUSIVE'),            # control misses target
     ((b(1.3, 1e-16), b(1.3, 1e-16, False), b(1, 1e-16), b(1, 1e-16)), 'INCONCLUSIVE'),   # non-finite
     ((b(2.9, 1e-16), b(3.2, 1e-16), b(1, 1e-16), b(1, 1e-16)), 'INCONCLUSIVE')]         # ratio > 3 at one resolution only
res = [outcome(*a)[0] == e for a, e in T]; print('aggregator tests', sum(res), '/', len(res)); assert all(res)
man = json.load(open(sys.argv[1])); N = 16384
for lab in ('PN_cut0.01_rank1_minus_0', 'PN_cut0.01_rank1_minus_1'):
    c = [e for e in man['PN'] if e['label'] == lab][0]; th = C.theta_from_case(c, True); S, x, dx, kk, eps = C.initial_state(th, N)
    B = np.array([np.full(N, th['w']), np.full(N, th['q']), np.full(N, th['v']), np.full(N, th['r']), np.full(N, th['Abar']), np.zeros(N), np.zeros(N), np.zeros(N)])
    r = np.abs(C.rhs(B, th, kk)[:6]).max(); print(lab, 'background fixed-point residual', r); assert r <= 1e-14
    res_, ts, nr, Es = C.run(c, True, N, 0.025, T=25.0); print(lab, 'T=25: finite', res_['finite'], 'ratio', res_['norm_ratio_max'], 'drift', res_['drift_total_max']); assert res_['finite'] and res_['drift_total_max'] < 1e-9
print('ALL TESTS PASS')
