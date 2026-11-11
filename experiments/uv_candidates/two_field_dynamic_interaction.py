# TwoField_DynamicInteraction (Hypothesis B, first pass). Lorentz-invariant heavy-scalar model, lab-frame periodic grid.
# Pre-registered measurements: (i) centroid shift of each packet vs g=0 control, (ii) change of central wavenumber and bandwidth after separation,
# (iii) residual energy in chi (heavy radiation), (iv) energy exchange phi<->psi, (v) amplitude scaling of the shift.  No thresholds tuned after the run.
import numpy as np
def run(A_phi=1.0,A_psi=1.0,a=0.2,b=0.2,e=0.4,M=1.0,k0=0.3,s=15.0,Lbox=1024.0,N=8192,T=200.0,dt=0.05,x1=-100.0,x2=100.0):
    x=np.linspace(-Lbox/2,Lbox/2,N,endpoint=False); dx=x[1]-x[0]; kk=2*np.pi*np.fft.fftfreq(N,dx)
    D=lambda f:np.real(np.fft.ifft(1j*kk*np.fft.fft(f)))
    g1=lambda x0,A:(A*np.exp(-(x-x0)**2/(2*s*s))*np.cos(k0*(x-x0)))
    phi=g1(x1,A_phi); psi=g1(x2,A_psi); chi=np.zeros(N)
    # right-moving phi: phi_t=-phi_x ; left-moving psi: psi_t=+psi_x
    pt=-D(phi); st=D(psi); ct=np.zeros(N)
    def rhs(y):
        phi,pt,psi,st,chi,ct=y
        px,sx,cx=D(phi),D(psi),D(chi)
        pxx,sxx,cxx=D(px),D(sx),D(cx)
        Xp=pt**2-px**2; Xs=st**2-sx**2; Y=pt*st-px*sx
        O=a*Xp+b*Xs+e*Y
        cc=cxx-M*M*chi+O          # chi_tt
        Rp=2*a*(ct*pt-cx*px)+e*(ct*st-cx*sx); Rs=2*b*(ct*st-cx*sx)+e*(ct*pt-cx*px)
        # Z [boxphi, boxpsi] = -[Rp, Rs] with box = f_tt - f_xx ; Z=[[1+2a chi, e chi],[e chi, 1+2b chi]]
        Z11=1+2*a*chi; Z12=e*chi; Z22=1+2*b*chi; det=Z11*Z22-Z12**2
        bp=(-Rp*Z22+Rs*Z12)/det; bs=(-Rs*Z11+Rp*Z12)/det
        return np.array([pt,pxx+bp,st,sxx+bs,ct,cc])
    y=np.array([phi,pt,psi,st,chi,ct]); n=int(round(T/dt)); snaps={}
    for i in range(n):
        k1=rhs(y);k2=rhs(y+.5*dt*k1);k3=rhs(y+.5*dt*k2);k4=rhs(y+dt*k3); y=y+dt/6*(k1+2*k2+2*k3+k4)
    phi,pt,psi,st,chi,ct=y; px,sx,cx=D(phi),D(psi),D(chi)
    return x,dx,kk,dict(phi=phi,pt=pt,px=px,psi=psi,st=st,sx=sx,chi=chi,ct=ct,cx=cx)
def meas(x,dx,kk,S,which):
    f,ft,fx=(S['phi'],S['pt'],S['px']) if which=='phi' else (S['psi'],S['st'],S['sx'])
    en=0.5*(ft**2+fx**2); E=en.sum()*dx; cen=(x*en).sum()*dx/E
    F=np.fft.fft(f); P=np.abs(F)**2; m=kk>0; kc=(kk[m]*P[m]).sum()/P[m].sum(); kw=np.sqrt(((kk[m]-kc)**2*P[m]).sum()/P[m].sum())
    return E,cen,kc,kw
def Echi(S,dx,M=1.0): return (0.5*(S['ct']**2+S['cx']**2+M*M*S['chi']**2)).sum()*dx
if __name__=="__main__":
    out={}
    for tag,(Ap,As,a,b,e) in {'ctrl':(1,1,0,0,0),'g':(1,1,.2,.2,.4),'g_Apsi2':(1,2,.2,.2,.4),'ctrl_Apsi2':(1,2,0,0,0)}.items():
        x,dx,kk,S=run(Ap,As,a,b,e)
        out[tag]=dict(phi=meas(x,dx,kk,S,'phi'),psi=meas(x,dx,kk,S,'psi'),Echi=Echi(S,dx),maxchi=float(np.abs(S['chi']).max()))
        print(tag,"phi (E,x_c,k_c,k_w)=%s psi=%s Echi=%.3e"%(np.round(out[tag]['phi'],5),np.round(out[tag]['psi'],5),out[tag]['Echi']),flush=True)
    for tag,ctl in (('g','ctrl'),('g_Apsi2','ctrl_Apsi2')):
        dxphi=out[tag]['phi'][1]-out[ctl]['phi'][1]; dxpsi=out[tag]['psi'][1]-out[ctl]['psi'][1]
        dk=out[tag]['phi'][2]-out[ctl]['phi'][2]; dw=out[tag]['phi'][3]-out[ctl]['phi'][3]
        dE=(out[tag]['phi'][0]-out[ctl]['phi'][0])/out[ctl]['phi'][0]
        print(tag,"shift phi %.4e psi %.4e ; d k_c %.3e ; d k_w %.3e ; rel dE_phi %.3e"%(dxphi,dxpsi,dk,dw,dE))
