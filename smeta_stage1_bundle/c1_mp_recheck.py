# Supplementary (post-hoc, diagnostic only; does not alter any registered outcome):
# recompute max_k gamma(k) in 50-digit arithmetic at the grid points where the double-precision value exceeded the C1 bound.
import json, sys, os, numpy as np, mpmath as mp
from multiprocessing import Pool
sys.path.insert(0,'.')
import t3fund_smeta_stage1 as S
R=json.load(open('smeta_stage1_result.json')); man=json.load(open('t3fund_smeta_manifest.json'))
P={json.dumps([p['pop'],p['cs'],p.get('idx'),p['f'],p.get('label')]):p for p in S.build(man)}
bad=[r for r in R if 'C1' in r and not r['C1']['gamma_within_bound']]
def work(r):
    p=P[r['_cid']]; g,f,bg,Hmp,one=S.to_mp(p,50); Hd=[[float(x) for x in row] for row in Hmp]
    gb=r['C1']['gamma_max']; worst=0.0; npts=0; nodata=0
    for kv in S.KGRID:
        rd,ir=S.roots_dbl(Hd,kv); gd=max(0.0,rd.imag.max())*kv
        if gd>gb*(1+1e-9)+1e-9:
            npts+=1; rm=S.roots_mp(Hmp,kv,50)
            if rm is None: nodata+=1; continue
            worst=max(worst,max(0.0,float(max(mp.im(c) for c in rm)))*kv)
    return dict(_cid=r['_cid'],pop=r['case']['pop'],gamma_max=gb,n_points=npts,n_nodata=nodata,max_gamma_mp=worst,within=bool(worst<=gb*(1+1e-9)+1e-9 and nodata==0))
if __name__=='__main__':
    out='c1_mp_recheck.jsonl'; done=set()
    if os.path.exists(out): done={json.loads(l)['_cid'] for l in open(out)}
    todo=[r for r in bad if r['_cid'] not in done]; print('todo',len(todo),'of',len(bad),flush=True)
    with Pool(2) as pool, open(out,'a') as fh:
        for x in pool.imap_unordered(work,todo):
            fh.write(json.dumps(x)+'\n'); fh.flush()
