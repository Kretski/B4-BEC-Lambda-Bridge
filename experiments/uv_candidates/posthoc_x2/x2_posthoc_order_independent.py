"""POST-HOC diagnostic (not part of the Script A run): recompute X2 order-independently from the saved roots.
Usage: python3 -I x2_posthoc_order_independent.py <uv_dynamic_eta_results.json>   (writes to stdout)"""
import json
import sys

d = json.load(open(sys.argv[1], encoding='utf-8'))
C = lambda s: complex(s.replace(' ', ''))
print("POST-HOC (not part of the run): X2 recomputed order-independently from the saved roots (15 significant digits as stored)")
worst = 0
for b in d['backgrounds']:
    for rr in b['per_root']:
        for f, pf in rr['per_f'].items():
            s7 = pf['S7_numeric']; v = s7.get('X2_full_vs_eliminated_maxdiff')
            fr = [C(x) for x in s7['lowk_full_system_roots']]; el = [C(x) for x in s7['roots_eta_eliminated']]
            md = max(min(abs(a - x) for x in el) for a in fr)
            md2 = max(min(abs(a - x) for x in fr) for a in el)
            worst = max(worst, md, md2)
            if float(v) > 1e-30:
                print(f"{b['name']} root {rr['root_index']} f={f}: stored X2 {v}; order-independent {max(md, md2):.2e}; complex conjugate pair present: {any(abs(z.imag) > 0 for z in fr)}")
print("worst order-independent difference over all 117 (root, f):", worst)
