import json, sys, numpy as np, mpmath as mp
sys.path.insert(0,'.')
import t3fund_smeta_stage1 as S
R=json.load(open('smeta_stage1_result.json')); man=json.load(open('t3fund_smeta_manifest.json'))
P={json.dumps([p['pop'],p['cs'],p.get('idx'),p['f'],p.get('label')]):p for p in S.build(man)}
bad=[r for r in R if 'C1' in r and not r['C1']['gamma_within_bound']]
import random; random.seed(1); sample=random.sample([r for r in bad if r['case']['pop']=='S1'],4)+random.sample([r for r in bad if r['case']['pop']=='Z0'],3)
for r in sample:
    p=P[r['_cid']]; g,f,bg,Hmp,one=S.to_mp(p,50); Hd=[[float(x) for x in row] for row in Hmp]
    gb=r['C1']['gamma_max']; worst=None
    for kv in S.KGRID:
        rd,ir=S.roots_dbl(Hd,kv); gd=max(0.0,rd.imag.max())*kv
        if gd>gb*(1+1e-9)+1e-9 and (worst is None or gd>worst[1]): worst=(kv,gd)
    kv,gd=worst; rm=S.roots_mp(Hmp,kv,50); gm=max(0.0,float(max(mp.im(c) for c in rm)))*kv if rm is not None else None
    print(r['case']['pop'],r['case']['cs'],r['case']['idx'],r['case']['f'],'gamma_max',round(gb,4),'worst k',kv,'dbl gamma',round(gd,4),'mp50 gamma',gm)
