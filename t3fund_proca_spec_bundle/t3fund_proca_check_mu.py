import numpy as np
rng=np.random.default_rng(5)
def mats(gp,gs,Lam,w,q,v,r):
    Ab=(gp*w*q+gs*v*r)
    Ps=np.zeros((6,6),complex); Pk=np.zeros((6,6),complex)
    Ps[0,1]=1j*(1+3*Lam*q**2); Ps[0,0]=-2j*gp*Ab; Ps[0,4]=-gp*1j*w
    Ps[1,0]=1j; Ps[2,3]=1j; Ps[2,2]=-2j*gs*Ab; Ps[2,4]=-gs*1j*v; Ps[3,2]=1j
    Pk[0,5]=-gp*q; Pk[2,5]=-gs*r; Pk[4,5]=1
    Pk[5,0]=gp*q; Pk[5,1]=gp*w; Pk[5,2]=gs*r; Pk[5,3]=gs*v; Pk[5,4]=-1
    return Ps,Pk
viol=0; n=0; stats=[]
for t in range(400):
    gp,gs=rng.uniform(-1,1,2); Lam=rng.uniform(0,1); w,q,v,r=rng.uniform(-1.2,1.2,4)
    Ps,Pk=mats(gp,gs,Lam,w,q,v,r); Pinf=Ps+Pk
    lam,V=np.linalg.eig(Pinf)
    if np.abs((1j*lam).imag).max()>1e-9 or np.linalg.cond(V)>1e6: continue     # need real simple speeds
    sep=min(abs(lam[i]-lam[j]) for i in range(6) for j in range(i)); 
    if sep<1e-6: continue
    Vi=np.linalg.inv(V); n+=1
    for k in np.geomspace(1e-3,1e4,40):
        kh=np.sqrt(1+k*k); Mx=k*Ps+kh*Pk
        gam=np.linalg.eigvals(Mx).real.max()
        Rr=Vi@((kh-k)*Pk)@V; mu=np.linalg.eigvalsh((Rr+Rr.conj().T)/2).max()
        if gam>mu+1e-9: viol+=1
        stats.append((gam,mu))
print('backgrounds with strictly hyperbolic P_inf:',n,' violations of gamma<=mu:',viol,' max gamma',max(s[0] for s in stats))
