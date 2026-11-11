# Probe on a stationary background (Hypothesis A in dynamics). Lorentz-invariant heavy-scalar two-field model, M=1, lab-frame periodic box.
# Background: phi0=w t+q x, psi0=v t+r x, chi_bar=O/M^2 (exact solution; checked below). Linear probe: spectral time series at fixed k, compared with roots of det K.
import numpy as np, sympy as sp
a,b,e,M=0.2,0.2,0.4,1.0
w,q,v,r=1.0,0.3,0.8,-0.2
# ---------- analytic: det K for fields (phi,psi,chi), from jets Hessian
ft,fx,gt,gx,Ct,Cx,Cc=sp.symbols('ft fx gt gx Ct Cx Cc',real=True)
Xp=ft**2-fx**2; Xs=gt**2-gx**2; Yy=ft*gt-fx*gx
L=sp.Rational(1,2)*Xp+sp.Rational(1,2)*Xs+sp.Rational(1,2)*(Ct**2-Cx**2)-sp.Rational(1,2)*M**2*Cc**2+Cc*(a*Xp+b*Xs+e*Yy)
jets=[ft,fx,gt,gx,Ct,Cx,Cc]
Hf=sp.lambdify(jets,sp.hessian(L,jets),'numpy')
Obar=a*(w*w-q*q)+b*(v*v-r*r)+e*(w*v-q*r); chibar=Obar/M**2
H=np.array(Hf(w,q,v,r,0.0,0.0,chibar),dtype=float)
def K(om,k):
    # fields order: phi(ft,fx), psi(gt,gx), chi(Ct,Cx,Cc): map plane-wave
    T=np.zeros((7,3),dtype=complex)
    T[0,0]=-1j*om;T[1,0]=1j*k;T[2,1]=-1j*om;T[3,1]=1j*k;T[4,2]=-1j*om;T[5,2]=1j*k;T[6,2]=1
    return T.conj().T@H@T
def roots(k):
    cs=np.linspace(-4,4,41); d=np.array([np.linalg.det(K(c*k,k)).real for c in cs])
    # det is polynomial of degree 6 in c (after dividing nothing): fit
    co=np.polyfit(cs,d,6); rt=np.roots(co); rt=rt[np.abs(rt.imag)<1e-6].real
    return np.sort(rt)
# self-consistency of background: chi EOM  -M^2 chibar + O = 0 by construction; phi,psi EOM trivially satisfied for constants.
print("chibar = %.4f ; Z=[[1+2a chi, e chi],[e chi,1+2b chi]] eigenvalues:"%chibar,np.linalg.eigvalsh(np.array([[1+2*a*chibar,e*chibar],[e*chibar,1+2*b*chibar]])))
ks_probe=[0.1,0.2,0.3,0.5]
for k in ks_probe: print("k=%.2f analytic c roots:"%k,np.round(roots(k),4))
# low-k fit of light modes: c(k)=c0+A k^2
xs=np.array([0.02,0.04,0.06,0.08]); rs=np.array([roots(x) for x in xs])
print("light roots at small k:",[np.round(r_,4) for r_ in rs[:1]])

# ---------- time-domain simulation of a small probe on the background (full nonlinear equations, fluctuation variables)
def simulate(eps,N=2048,Lbox=1024.0,T=1600.0,dt=0.1,s=60.0,seed=0):
    x=np.linspace(-Lbox/2,Lbox/2,N,endpoint=False); dx=x[1]-x[0]; kk=2*np.pi*np.fft.fftfreq(N,dx)
    D=lambda f:np.real(np.fft.ifft(1j*kk*np.fft.fft(f)))
    rng=np.random.default_rng(seed)
    env=np.exp(-x**2/(2*s*s))
    dphi=eps*env*np.cos(0.3*x)+eps*env*np.cos(0.1*x)+eps*env*np.cos(0.2*x)+eps*env*np.cos(0.5*x)
    dpsi=0.5*dphi; zero=np.zeros(N)
    y=np.array([dphi,zero,dpsi,zero,zero,zero])   # (dphi, dphi_t, dpsi, dpsi_t, dchi, dchi_t); dchi=0 initially
    def rhs(y):
        dp,dpt,ds_,dst,dc,dct=y
        px=q+D(dp);pt=w+dpt;sx=r+D(ds_);st=v+dst
        chi=chibar+dc; cx=D(dc); ct=dct
        pxx=D(D(dp)); sxx=D(D(ds_)); cxx=D(cx)
        Xp_=pt**2-px**2; Xs_=st**2-sx**2; Y_=pt*st-px*sx
        O=a*Xp_+b*Xs_+e*Y_
        cc=cxx-M*M*chi+O
        Rp=2*a*(ct*pt-cx*px)+e*(ct*st-cx*sx); Rs=2*b*(ct*st-cx*sx)+e*(ct*pt-cx*px)
        Z11=1+2*a*chi;Z12=e*chi;Z22=1+2*b*chi;det=Z11*Z22-Z12**2
        bp=(-Rp*Z22+Rs*Z12)/det; bs=(-Rs*Z11+Rp*Z12)/det
        return np.array([dpt,pxx+bp,dst,sxx+bs,dct,cc])
    # sanity: background alone (zero fluctuation) must stay stationary
    n=int(round(T/dt)); ik={kj:int(round(kj/(kk[1]))) for kj in ks_probe}
    series={kj:[] for kj in ks_probe}
    for i in range(n):
        k1=rhs(y);k2=rhs(y+.5*dt*k1);k3=rhs(y+.5*dt*k2);k4=rhs(y+dt*k3); y=y+dt/6*(k1+2*k2+2*k3+k4)
        Fk=np.fft.fft(y[0])
        for kj in ks_probe: series[kj].append(Fk[ik[kj]])
    return np.array(dt*np.arange(1,n+1)),{kj:np.array(v_) for kj,v_ in series.items()},kk,ik
def peaks(t,sig,dt,kj,ik_k,kgrid):
    kx=kgrid[ik_k]
    S=np.fft.fft(sig*np.hanning(len(sig)),8*len(sig)); om=2*np.pi*np.fft.fftfreq(8*len(sig),dt); P=np.abs(S)
    # local maxima
    idx=[i for i in range(1,len(P)-1) if P[i]>P[i-1] and P[i]>=P[i+1] and P[i]>0.02*P.max()]
    return kx,sorted([(om[i]/kx) for i in idx])   # phase speeds omega/k (sign from e^{-i omega t} convention)
for eps in (1e-6,1e-4):
    t,ser,kk,ik=simulate(eps)
    print("eps=%g"%eps)
    for kj in ks_probe:
        kx,cs=peaks(t,ser[kj],0.1,kj,ik[kj],kk)
        print("  k=%.4f measured c (omega/k) peaks:"%kx,np.round(cs,4),"| analytic light:",np.round([x for x in roots(kx) if abs(x)<0.95],4))
