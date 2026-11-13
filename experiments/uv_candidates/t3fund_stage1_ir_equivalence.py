# Post-run analysis (not part of the registered procedure): does FAIL-UNSTABLE coincide with complex roots of the k -> 0 limit of det (the local B4-type polynomial F(c))?
# usage: python t3fund_stage1_ir_equivalence.py t3fund_stage1_result.json
import json, sys, numpy as np
R = json.load(open(sys.argv[1])); p = np.polynomial.polynomial
def Floc(model, th):
    if model == 'B':
        g, w, q = th; Z = 1+2*g*g*(w*w-q*q)
        return np.roots(p.polyadd([Z, 0, -Z], [-4*g*g*q*q, -8*g*g*w*q, -4*g*g*w*w])[::-1])
    gp, gs, w, q, v, r = th; sb = gp*w*q+gs*v*r
    Aa = [1, 2*gp*sb, -1]; Ab = [1, 2*gs*sb, -1]
    F = p.polysub(p.polysub(p.polymul(Aa, Ab), gp*gp*p.polymul([w*w, -2*w*q, q*q], Ab)), gs*gs*p.polymul([v*v, -2*v*r, r*r], Aa))
    return np.roots(F[::-1])
for name, model in (('A1_J0eq0_seed11', 'A'), ('A2_manifest40', 'A'), ('A3_general_seed41', 'A'), ('B1_Zpos_seed31', 'B')):
    rows = R[name]['rows']; agree = 0
    for r in rows:
        rt = Floc(model, r['theta']); cplx = np.abs(rt.imag).max() > 1e-9*max(1, np.abs(rt).max())
        agree += (cplx == (r['outcome'] != 'PASS'))
    print(name, 'agree', agree, 'of', len(rows))
