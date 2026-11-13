# T3-fund HS2: T7-b numerical operator matching, per registered spec section 6 (SPEC_HASH 9d52c927a9acd3bd).
# HS2 k->0 local polynomial: P(c) = D(c) * F(c), D=1-c^2, F(c)=M^2 D detZ - 4 s^T adj(Z) s (monic-normalised); B4 polynomial as in heavy_scalar_two_field_matching.py.
# Residual per background: ||(P_B4(x)-P_target)[:4]|| / ||(P_target-free)[:4]||, free=(c^2-1)^2. Joint least-squares over 15 backgrounds, 60 random starts.
import numpy as np, json, sys
from scipy.optimize import least_squares
P_=np.polynomial.polynomial
FREE=np.array([1.,0.,-2.,0.,1.])
def hs2_poly(a,b,e,w,q,v,r,M=1.0):
    Cp=np.array([[a,e/2],[e/2,b]]); wv=np.array([w,v]); qv=np.array([q,r])
    chi=(a*(w*w-q*q)+b*(v*v-r*r)+e*(w*v-q*r))/M**2
    Z=np.eye(2)+2*chi*Cp; dZ=np.linalg.det(Z); adj=np.array([[Z[1,1],-Z[0,1]],[-Z[1,0],Z[0,0]]])
    s0=Cp@qv; s1=Cp@wv; ss=np.zeros(3)
    for i in range(2):
        for j in range(2): ss=P_.polyadd(ss,adj[i,j]*P_.polymul([s0[i],s1[i]],[s0[j],s1[j]]))
    D=np.array([1.,0.,-1.]); F=P_.polyadd(M**2*dZ*D,-4*ss); Pq=P_.polymul(D,F)
    return Pq/Pq[-1]
def b4_poly(bg,x):
    w,q,v,r=bg; p,s,kap,Lam=x
    ap=p*w*q-kap*v*r; aps=-kap*w*q+s*v*r
    d1=np.array([-(1+3*Lam*q*q),-2*ap,1.]); d2=np.array([-1.,-2*aps,1.])
    F=P_.polyadd(P_.polyadd(P_.polymul(d1,d2),p*P_.polymul(P_.polymul([w,-q],[w,-q]),d2)),s*P_.polymul(P_.polymul([v,-r],[v,-r]),d1))
    return F/F[-1]
def fit(bgs,T,rng,nstart=60):
    dev=[np.linalg.norm((t-FREE)[:4]) for t in T]
    def f(x):
        out=[]
        for bg,t,d in zip(bgs,T,dev):
            try: pb=b4_poly(bg,x)
            except Exception: return np.full(4*len(bgs),1e6)
            out.append((pb[:4]-t[:4])/d)
        o=np.concatenate(out); return np.where(np.isfinite(o),o,1e6)
    best=None
    for _ in range(nstart):
        x0=rng.normal(0,0.5,4)
        try: sol=least_squares(f,x0,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
        except Exception: continue
        if best is None or sol.cost<best.cost: best=sol
    per=[np.linalg.norm(c) for c in best.fun.reshape(len(bgs),4)]
    return best.x,max(per),float(np.sqrt(np.mean(np.square(per)))),per
if __name__=='__main__':
    man=json.load(open('t3fund_hs2_manifest.json')); rng=np.random.default_rng(71); res={}
    for cs in ['K1','K3','K4']:
        sel=[d for d in man['H1'] if d['cs']==cs and d['F_real']][:15]
        bgs=[(d['w'],d['q'],d['v'],d['r']) for d in sel]; a,b,e=sel[0]['a'],sel[0]['b'],sel[0]['e']
        T=[hs2_poly(a,b,e,*bg) for bg in bgs]
        x,mx,rms,per=fit(bgs,T,rng)
        Tc=[b4_poly(bg,(0.05,0.04,0.0447,0.3)) for bg in bgs]; xc,mxc,rmsc,_=fit(bgs,Tc,rng)
        dev=[float(np.linalg.norm((t-FREE)[:4])) for t in T]
        res[cs]=dict(idx=[d['idx'] for d in sel],fit_x=list(map(float,x)),max_rel_resid=float(mx),rms_rel_resid=rms,per_bg=list(map(float,per)),
                     dev_from_free_min=min(dev),dev_from_free_max=max(dev),control_max=float(mxc),control_x=list(map(float,xc)),
                     outcome='MATCH' if (mx<=1e-6) else 'NO-MATCH',control_ok=bool(mxc<=1e-12))
        print(cs,res[cs]['outcome'],'max %.3e rms %.3e control %.2e'%(mx,rms,mxc),'x',np.round(x,4),flush=True)
    json.dump(res,open('hs2_t7b_result.json','w'),indent=1)
