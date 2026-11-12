# B4 local characteristic cone (principal symbol of the second-order B4 EFT) vs the unit (fundamental) cone.
# F(c) = d1 d2 + p (w - q c)^2 d2 + s (v - r c)^2 d1,  d1 = c^2 - 2 a_phi c - (1 + 3 Lam q^2),  d2 = c^2 - 2 a_psi c - 1,  kappa^2 = p s.
import numpy as np, json
P=np.polynomial.polynomial
def roots(p,s,kap,Lam,w,q,v,r):
    ap=p*w*q-kap*v*r; aps=-kap*w*q+s*v*r
    d1=np.array([-(1+3*Lam*q*q),-2*ap,1.]); d2=np.array([-1.,-2*aps,1.])
    F=P.polyadd(P.polyadd(P.polymul(d1,d2),p*P.polymul(P.polymul([w,-q],[w,-q]),d2)),s*P.polymul(P.polymul([v,-r],[v,-r]),d1))
    return P.polyroots(F)
rng=np.random.default_rng(5); out={}
for Lam in (0.0,0.01,0.1,1.0):
    for cut in (0.09,0.95):
        n=0; real=0; sub=0; mx=[]
        while n<4000:
            p,s=rng.uniform(0.02,1,2); kap=rng.choice([-1,1])*np.sqrt(p*s); w,q,v,r=rng.uniform(-1.5,1.5,4)
            if max(p*q*q,s*r*r,p*w*w,s*v*v)>cut: continue
            n+=1; rt=roots(p,s,kap,Lam,w,q,v,r)
            if np.abs(rt.imag).max()>1e-9: continue
            real+=1; m=np.abs(rt.real).max(); mx.append(m); sub+= m<=1+1e-9
        out[f"Lam={Lam},cut={cut}"]=dict(n=n,hyperbolic=real,max_speed_le_1=int(sub),min_max_speed=float(min(mx)),median_max_speed=float(np.median(mx)))
print(json.dumps(out,indent=1))
# q = r = 0 (purely time-dependent background): analytic roots +-1, +-sqrt(1 - p w^2 - s v^2)
print("q=r=0 check:",np.sort(roots(0.3,0.2,np.sqrt(0.06),0.7,0.8,0.0,0.5,0.0).real), "expected +-1, +-",np.sqrt(1-0.3*0.64-0.2*0.25))
