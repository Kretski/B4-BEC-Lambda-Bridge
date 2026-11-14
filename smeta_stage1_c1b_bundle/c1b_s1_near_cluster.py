# Diagnostic for the 5 S1 cases whose registered double-precision check (b) failed but whose det K roots are not exactly multiple.
# Reports: min root separation of det K at k=0.7 (60 digits), conditioning of the double-precision first-order eigenproblem, and the roots.
import json,sys,numpy as np,mpmath as mp
sys.path.insert(0,'.')
import t3fund_smeta_stage1 as S
man=json.load(open('t3fund_smeta_manifest.json')); P={json.dumps([p['pop'],p['cs'],p.get('idx'),p['f'],p.get('label')]):p for p in S.build(man)}
ids=['["S1", "G1", 124, 10.0, null]','["S1", "G1", 146, 1.0, null]','["S1", "G1", 146, 10.0, null]','["S1", "G1", 165, 10.0, null]','["S1", "G3", 188, 10.0, null]']
for cid in ids:
    mp.mp.dps=60; g,f,bg,H,one=S.to_mp(P[cid],60); Hd=[[float(x) for x in r] for r in H]
    d=S.poly_c(H,mp.mpf('0.7'),False); co=[mp.re(x) for x in d]; r=mp.polyroots(co[::-1],maxsteps=2000,extraprec=500)
    r=[x*mp.mpf('0.7') for x in r]
    sep=min(abs(r[i]-r[j]) for i in range(len(r)) for j in range(i))
    A0,A1,B=S.first_order(Hd); P0=-np.linalg.inv(A0)@A1; Cs=-np.linalg.inv(A0)@B; Mx=1j*0.7*P0+Cs
    w,V=np.linalg.eig(Mx)
    print(cid,'min_sep',float(sep),'cond_V_FO_k0.7',np.linalg.cond(V),'norm_M',np.linalg.norm(Mx),'roots',[round(float(mp.re(x)),6) for x in r])
