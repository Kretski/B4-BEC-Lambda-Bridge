# T3-fund SM+dynamical-eta spec checks (sympy/numpy). Verifies: (1) background equations and closed forms; (2) linearised operator K(omega,k) built from the jet Hessian equals
# the linearisation of the full Euler-Lagrange equations computed independently with explicit plane waves; (3) principal-symbol closed forms; (4) k=0 Schur closed form.
import numpy as np, sympy as sp
rng=np.random.default_rng(5)
# ---- Lagrangian (M=1, signature +-): L=1/2(dphi)^2+1/2(dpsi)^2+1/2(dsig)^2-1/2 M^2 sig^2+1/2 f (deta)^2+sig*(gp*ab(phi)+gs*ab(psi)),  ab=(ch ft+sh fx)(sh ft+ch fx)
names=['P','Q','S','E']; J={n:sp.symbols(f'{n} {n}t {n}x') for n in names}; gp,gs,f,M=sp.symbols('gp gs f M',real=True)
def ab(ft,fx,eta): return (sp.cosh(eta)*ft+sp.sinh(eta)*fx)*(sp.sinh(eta)*ft+sp.cosh(eta)*fx)
Pt,Px=J['P'][1],J['P'][2]; Qt,Qx=J['Q'][1],J['Q'][2]; sg=J['S'][0]; et=J['E'][0]
L=sp.Rational(1,2)*(Pt**2-Px**2)+sp.Rational(1,2)*(Qt**2-Qx**2)+sp.Rational(1,2)*(J['S'][1]**2-J['S'][2]**2)-sp.Rational(1,2)*M**2*sg**2+sp.Rational(1,2)*f*(J['E'][1]**2-J['E'][2]**2)+sg*(gp*ab(Pt,Px,et)+gs*ab(Qt,Qx,et))
allj=[x for n in names for x in J[n]]; Hm=sp.hessian(L,allj); Hf=sp.lambdify([x for x in allj]+[gp,gs,f,M],Hm,'numpy')
Lf=sp.lambdify([x for x in allj]+[gp,gs,f,M],L,'numpy'); dL=[sp.lambdify([x for x in allj]+[gp,gs,f,M],sp.diff(L,x),'numpy') for x in allj]
def bgS(g1,g2,w,q,v,r):
    A=g1*w*q+g2*v*r; B=0.5*(g1*(w*w+q*q)+g2*(v*v+r*r))
    if abs(B)>=abs(A): return None
    eta=0.5*np.arctanh(-B/A); sig=np.sign(A)*np.sqrt(A*A-B*B); return A,B,eta,sig
def jets(w,q,v,r,sig,eta): return [w,w,q][0:0]+[0,w,q, 0,v,r, sig,0,0, eta,0,0]
def K_from_H(vals,om,k,g1,g2,fv):
    H=np.array(Hf(*vals,g1,g2,fv,1.0),dtype=complex); T=np.zeros((12,4),dtype=complex)
    for i in range(4): T[3*i,i]=1;T[3*i+1,i]=-1j*om;T[3*i+2,i]=1j*k
    return T.conj().T@H@T
# ---- (2) independent linearisation with explicit plane waves, numeric floats
t,x,eps,om,k=sp.symbols('t x eps om k',real=True); a=sp.symbols('a0:4')
def lin_matrix(vals_bg,g1,g2,fv):
    w,q,v,r,sig,eta=vals_bg; Eexp=sp.exp(sp.I*(k*x-om*t))
    fld=[w*t+q*x+eps*a[0]*Eexp, v*t+r*x+eps*a[1]*Eexp, sig+eps*a[2]*Eexp, eta+eps*a[3]*Eexp]
    sub={}
    for n,F in zip(names,fld):
        sub[J[n][0]]=F; sub[J[n][1]]=sp.diff(F,t); sub[J[n][2]]=sp.diff(F,x)
    par={gp:g1,gs:g2,f:fv,M:1.0}
    rows=[]
    for n,F in zip(names,fld):
        d0=sp.diff(L,J[n][0]).subs(par).subs(sub); dt=sp.diff(L,J[n][1]).subs(par).subs(sub); dx=sp.diff(L,J[n][2]).subs(par).subs(sub)
        EL=d0-sp.diff(dt,t)-sp.diff(dx,x)       # Euler-Lagrange expression
        rows.append(EL)
    # background residual (eps=0) and first order
    bgres=[float(sp.N(r_.subs(eps,0).subs({t:0.3,x:0.7}))) for r_ in rows]
    lin=[sp.simplify(sp.diff(r_,eps).subs(eps,0)/Eexp) for r_ in rows]
    Mat=sp.Matrix([[sp.diff(l,a[j]) for j in range(4)] for l in lin]); return bgres,sp.lambdify((om,k),Mat,'numpy')
worst_bg=0; worst_det=0; worst_cf=0; n=0
while n<4:
    g1,g2=rng.uniform(-0.5,0.5,2); fv=rng.choice([0.1,1.0,10.0]); w,q,v,r=rng.uniform(-1.5,1.5,4); b=bgS(g1,g2,w,q,v,r)
    if b is None: continue
    A,B,eta,sig=b; n+=1
    bgres,Mf=lin_matrix((w,q,v,r,sig,eta),g1,g2,fv); worst_bg=max(worst_bg,max(abs(z) for z in bgres))
    vals=[0,w,q,0,v,r,sig,0,0,eta,0,0]
    for (c,kk) in ((0.4,1.3),(-1.7,0.6),(2.2,3.1)):
        Kx=K_from_H(vals,c*kk,kk,g1,g2,fv); d1=np.linalg.det(Kx); d2=np.linalg.det(np.array(Mf(c*kk,kk),dtype=complex))
        worst_det=max(worst_det,min(abs(d1-d2),abs(d1+d2))/max(abs(d1),1e-30))
    # (1) closed forms
    H=np.array(Hf(*vals,g1,g2,fv,1.0),dtype=float)
    for gi,(al,be),idx in ((g1,(w,q),0),(g2,(v,r),3)):
        Zt=H[idx+1,idx+1]; Xg=-H[idx+2,idx+2]; Y=H[idx+1,idx+2]
        worst_cf=max(worst_cf,abs(Zt-(1-gi*B)),abs(Xg-(1+gi*B)),abs(Y-gi*A),abs(Y**2+Zt*Xg-(1+gi**2*sig**2)))
        Nt=H[idx+1,9]; worst_cf=max(worst_cf,abs(Nt-2*gi*(A*al-B*be)))
    worst_cf=max(worst_cf,abs(H[9,9]-4*sig**2),abs(H[6,9]))     # eta-eta value part 4 sigma^2 ; sigma-eta value coupling vanishes on stationary branch
print("max |EL residual on background| (phi,psi,sigma,eta):",worst_bg)
print("max rel diff det K (jet Hessian) vs det(linearised EL):",worst_det)
print("max abs error of closed forms (Zt=1-gB, X=1+gB, Y=gA, Y^2+ZtX=1+g^2 sig^2, N=2g(A al-B be), H_etaeta=4 sig^2, H_sigma-eta=0):",worst_cf)
# k=0 identity: det K(omega,0) = omega^4 Zt1 Zt2 * [ (omega^2-M^2-sum m^2/Zt)(f omega^2+4 sig^2-sum N^2/Zt) - (sum m N/Zt)^2 ],  m_i=g_i (A be_i - B al_i)/sig, N_i=2 g_i (A al_i - B be_i)
worst=0; m_=0
while m_<30:
    g1,g2=rng.uniform(-0.5,0.5,2); fv=rng.choice([0.1,1.0,10.0]); w,q,v,r=rng.uniform(-1.5,1.5,4); b=bgS(g1,g2,w,q,v,r)
    if b is None: continue
    A,B,eta,sig=b; Zts=[1-g1*B,1-g2*B]; m_+=1; vals=[0,w,q,0,v,r,sig,0,0,eta,0,0]
    al=[w,v]; be=[q,r]; gg=[g1,g2]
    Ns=[2*gg[i]*(A*al[i]-B*be[i]) for i in range(2)]; ms=[gg[i]*(A*be[i]-B*al[i])/sig for i in range(2)]
    for omv in (0.37,1.9,-0.8):
        d=np.linalg.det(K_from_H(vals,omv,0.0,g1,g2,fv)).real
        S1=sum(mi**2/Zi for mi,Zi in zip(ms,Zts)); S2=sum(Ni**2/Zi for Ni,Zi in zip(Ns,Zts)); S3=sum(mi*Ni/Zi for mi,Ni,Zi in zip(ms,Ns,Zts))
        pred=omv**4*Zts[0]*Zts[1]*((omv**2-1.0-S1)*(fv*omv**2+4*sig**2-S2)-S3**2)
        worst=max(worst,abs(d-pred)/max(abs(pred),1e-30))
print("k=0 identity (quadratic in omega^2 after Schur elimination of phi,psi): max rel err",worst)
# principal part (k->infinity) closed forms: K(ck,k)/k^2 -> diag(Zt_i c^2-2Y_i c-X_i (i=phi,psi), c^2-1, f(c^2-1)); leading coefficient of det K in c^8 is Zt_phi Zt_psi f (up to sign (-1)^0)
worst_pr=0; worst_lead=0; m_=0
while m_<10:
    g1,g2=rng.uniform(-0.5,0.5,2); fv=rng.choice([0.1,1.0,10.0]); w,q,v,r=rng.uniform(-1.5,1.5,4); b=bgS(g1,g2,w,q,v,r)
    if b is None: continue
    A,B,eta,sig=b; m_+=1; vals=[0,w,q,0,v,r,sig,0,0,eta,0,0]; Zt=[1-g1*B,1-g2*B]; Xg=[1+g1*B,1+g2*B]; Y=[g1*A,g2*A]
    for c in (0.3,-1.4,2.0):
        kk=1e5; Kx=K_from_H(vals,c*kk,kk,g1,g2,fv)/kk**2
        pred=np.diag([Zt[0]*c*c-2*Y[0]*c-Xg[0],Zt[1]*c*c-2*Y[1]*c-Xg[1],c*c-1,fv*(c*c-1)])
        worst_pr=max(worst_pr,np.abs(Kx-pred).max()/max(np.abs(pred).max(),1e-12))
    cs_=np.linspace(-3,3,17); dv=np.array([np.linalg.det(K_from_H(vals,c*7.0,7.0,g1,g2,fv)).real for c in cs_]); co=np.polyfit(cs_,dv,8)
    worst_lead=max(worst_lead,abs(co[0]/7.0**8/(Zt[0]*Zt[1]*fv)-1))
print("principal part K(ck,k)/k^2 vs closed form (k=1e5): max rel diff (O(1/k) expected):",worst_pr)
print("leading c^8 coefficient of det K / k^8 vs Zt_phi Zt_psi f: max rel err",worst_lead)
