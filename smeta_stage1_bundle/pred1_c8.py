import json,sys,numpy as np,mpmath as mp
sys.path.insert(0,'.')
import t3fund_smeta_stage1 as S
R=json.load(open('smeta_stage1_result.json')); man=json.load(open('t3fund_smeta_manifest.json'))
P={json.dumps([p['pop'],p['cs'],p.get('idx'),p['f'],p.get('label')]):p for p in S.build(man)}
n=0; bad=0; worst=1e9
for r in R:
    k=r['kappa']
    if min(abs(x) for x in k)<=S.TAU: continue
    p=P[r['_cid']]; g,f,bg,Hmp,one=S.to_mp(p,50)
    d=S.poly_c(Hmp,mp.mpf(1000),False); lead=mp.re(d[-1]); kap=[mp.mpf(x) for x in k]
    pred=kap[0]*kap[1]*kap[3]  # Z_t,phi Z_t,psi f (kappa=(Zphi,Zpsi,1,f))
    n+=1
    if (lead>0)!=(pred>0): bad+=1
    ratio=abs(lead)/abs(pred); 
print('non-degenerate cases',n,'sign mismatches',bad)
