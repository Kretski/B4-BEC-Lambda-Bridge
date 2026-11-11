# Scalar-mediator (Lambda=0) evaluation, T1-T6 as pre-registered in B4_CLOSURE_and_B5_PLAN.md.
# Thresholds fixed BEFORE running: T4 sqrtX_max<=0.3 (g<=1,rho<=0.3); T3 front<=1+1e-9; T2 real roots on x grid [1e-3,1e3] (121 pts), rel tol 1e-7;
# T5 static: min eigenvalue of quadratic lab-frame Hamiltonian H2(k) >= -1e-12 on k/M in [1e-3,1e3].
import numpy as np, json
P=np.polynomial.polynomial
def setup(p,s,kap,w,q,v,r):
    ap=p*w*q-kap*v*r; aps=-kap*w*q+s*v*r
    d1=np.array([-1.,-2*ap,1.]); d2=np.array([-1.,-2*aps,1.])
    return ap,aps,d1,d2
def Qt(x2,p,s,kap,w,q,v,r):
    ap,aps,d1,d2=setup(p,s,kap,w,q,v,r)
    s1=np.array([w,-q]); s2=np.array([v,-r])
    t=P.polymul(np.array([1+x2,0,-x2]),P.polymul(d1,d2))
    return P.polyadd(P.polyadd(t,p*P.polymul(P.polymul(s1,s1),d2)),s*P.polymul(P.polymul(s2,s2),d1))
def T2(*a,xs=np.geomspace(1e-3,1e3,121)):
    for x in xs:
        rt=P.polyroots(Qt(x*x,*a))
        if np.abs(rt.imag).max()>1e-7*max(1,np.abs(rt).max()): return False
    return True
def T3front(*a):
    ap,aps,d1,d2=setup(*a)
    return max(1.0,np.abs(P.polyroots(d1)).max(),np.abs(P.polyroots(d2)).max())
def T5min(p,s,kap,w,q,v,r,M=1.0,ks=np.geomspace(1e-3,1e3,61)):
    gp=M*np.sqrt(p); gs=M*np.sqrt(s)*np.sign(kap)*-1  # kappa=-g_phi g_psi/M^2 -> g_psi=-kap*M^2/g_phi
    gs=-kap*M**2/gp
    sb=(gp*w*q+gs*v*r)/M**2
    mn=np.inf
    for k in ks*M:
        # z=(pi_phi,pi_psi,pi_sig,phi,psi,sig); H=1/2 z^dag Hm z
        def row(coefs): return np.array(coefs,dtype=complex)
        Apm=row([0,0,0, sb*gp*1j*k,0, gp*q]); Aps=row([0,0,0,0, sb*gs*1j*k, gs*r])
        u1=row([1,0,0,0,0,0])-Apm; u2=row([0,1,0,0,0,0])-Aps
        us=[u1,u2,row([0,0,1,0,0,0]),row([0,0,0,k,0,0]),row([0,0,0,0,k,0]),row([0,0,0,0,0,k]),row([0,0,0,0,0,M])]
        H=sum(np.outer(u.conj(),u) for u in us)
        # cross terms -g w sig phi_x - g v sig psi_x (Hermitian)
        for idx,c in ((3,gp*w),(4,gs*v)):
            H[5,idx]+=-c*1j*k; H[idx,5]+=np.conj(-c*1j*k)
        mn=min(mn,np.linalg.eigvalsh(H).min())
    return mn
def beta(p,s,kap,w,q,v,r):
    """O(X^2) shift of local roots: c=c0+A X^2, A=-(1-c0^2) d1 d2/F'(c0); beta=A/M^2 -> returns A for each root"""
    ap,aps,d1,d2=setup(p,s,kap,w,q,v,r)
    s1=np.array([w,-q]); s2=np.array([v,-r])
    F=P.polyadd(P.polyadd(P.polymul(d1,d2),p*P.polymul(P.polymul(s1,s1),d2)),s*P.polymul(P.polymul(s2,s2),d1))
    out=[]
    for c0 in P.polyroots(F):
        c0=c0.real if abs(c0.imag)<1e-9 else c0
        Fp=P.polyval(c0,P.polyder(F)); out.append((c0, -(1-c0**2)*P.polyval(c0,d1)*P.polyval(c0,d2)/Fp if abs(Fp)>1e-9 else np.nan))
    return out
rng=np.random.default_rng(11)
rows=[]
while len(rows)<400:
    p,s=rng.uniform(0.02,0.09,2); kap=rng.choice([-1,1])*np.sqrt(p*s)
    w,q,v=rng.uniform(-1.5,1.5,3)
    if abs(kap*v)<1e-3: continue
    r=p*w*q/(kap*v)
    if max(np.sqrt(p)*abs(q),np.sqrt(s)*abs(r),np.sqrt(p)*abs(w),np.sqrt(s)*abs(v))>0.3: continue
    rows.append((p,s,kap,w,q,v,r))
res=dict(n=len(rows)); 
res['T2']=sum(T2(*a) for a in rows)
res['T3_front1']=sum(T3front(*a)<=1+1e-9 for a in rows)
m5=[T5min(*a) for a in rows]; res['T5_pos']=int(sum(m>=-1e-12 for m in m5)); res['T5_min']=float(min(m5))
loc=[np.abs(P.polyroots(Qt(0.0,*a))).max() for a in rows]; res['local_speed_max_range']=(float(min(loc)),float(max(loc)))
# beta for non-luminal roots
bl=[]; ll=[]
for a in rows:
    for c0,A in beta(*a):
        if np.isfinite(A) and abs(abs(c0)-1)>1e-6: bl.append(A); ll.append(2*A/c0 if abs(c0)>1e-6 else np.nan)
bl=np.array(bl); ll=np.array(ll)
res['beta_nonluminal']=dict(n=len(bl),pos=int((bl>0).sum()),neg=int((bl<0).sum()),absmed=float(np.median(np.abs(bl))))
res['Lam_eff_pos']=int((ll[np.isfinite(ll)]>0).sum()); res['Lam_eff_neg']=int((ll[np.isfinite(ll)]<0).sum())
lum=[A for a in rows for c0,A in beta(*a) if abs(abs(c0)-1)<1e-6]; res['beta_luminal_maxabs']=float(np.max(np.abs(lum))) if lum else None
print(json.dumps(res,indent=1))
# near-J0=0 deformation: scale r by (1+d)
for d in (1e-3,1e-2,1e-1):
    fr=[]; 
    for (p,s,kap,w,q,v,r) in rows[:100]:
        fr.append(T3front(p,s,kap,w,q,v,r*(1+d)))
    print('r*(1+%g): front median %.4f max %.4f'%(d,np.median(fr),np.max(fr)))
