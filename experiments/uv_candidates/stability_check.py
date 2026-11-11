import numpy as np
a,b,e,M=0.2,0.2,0.4,1.0; w,q,v,r=1.0,0.3,0.8,-0.2
chibar=(a*(w*w-q*q)+b*(v*v-r*r)+e*(w*v-q*r))/M**2
N=2048;Lbox=1024.0;x=np.linspace(-Lbox/2,Lbox/2,N,endpoint=False);dx=x[1]-x[0];kk=2*np.pi*np.fft.fftfreq(N,dx)
D=lambda f:np.real(np.fft.ifft(1j*kk*np.fft.fft(f)))
def rhs(y):
    dp,dpt,ds,dst,dc,dct=y
    px=q+D(dp);pt=w+dpt;sx=r+D(ds);st=v+dst;chi=chibar+dc;cx=D(dc);ct=dct
    pxx=D(D(dp));sxx=D(D(ds));cxx=D(cx)
    O=a*(pt**2-px**2)+b*(st**2-sx**2)+e*(pt*st-px*sx)
    cc=cxx-M*M*chi+O
    Rp=2*a*(ct*pt-cx*px)+e*(ct*st-cx*sx);Rs=2*b*(ct*st-cx*sx)+e*(ct*pt-cx*px)
    Z11=1+2*a*chi;Z12=e*chi;Z22=1+2*b*chi;det=Z11*Z22-Z12**2
    return np.array([dpt,pxx+(-Rp*Z22+Rs*Z12)/det,dst,sxx+(-Rs*Z11+Rp*Z12)/det,dct,cc])
def energy(y):
    dp,dpt,ds,dst,dc,dct=y
    px=q+D(dp);pt=w+dpt;sx=r+D(ds);st=v+dst;chi=chibar+dc;cx=D(dc);ct=dct
    L=0.5*(pt**2-px**2)+0.5*(st**2-sx**2)+0.5*(ct**2-cx**2)-0.5*M*M*chi**2+chi*(a*(pt**2-px**2)+b*(st**2-sx**2)+e*(pt*st-px*sx))
    pphi=pt*(1+2*a*chi)+e*chi*st; ppsi=st*(1+2*b*chi)+e*chi*pt
    return ((pphi*pt+ppsi*st+ct*ct)-L).sum()*dx
eps=1e-3; env=np.exp(-x**2/(2*60.0**2)); dphi=eps*env*(np.cos(.3*x)+np.cos(.1*x)+np.cos(.2*x)+np.cos(.5*x)); z=np.zeros(N)
y=np.array([dphi,z,0.5*dphi,z,z,z]); dt=0.1;T=1600.0;n=int(T/dt)
E0=energy(y); nrm=lambda y:np.sqrt(((y[0]**2+y[2]**2+y[4]**2)).sum()*dx)
rec=[(0,nrm(y),0.0)]
for i in range(1,n+1):
    k1=rhs(y);k2=rhs(y+.5*dt*k1);k3=rhs(y+.5*dt*k2);k4=rhs(y+dt*k3);y=y+dt/6*(k1+2*k2+2*k3+k4)
    if i%2000==0: rec.append((i*dt,nrm(y),(energy(y)-E0)/abs(E0)))
Eb=energy(np.array([z,z,z,z,z,z]))
for t,nv,de in rec: print("t=%6.0f  ||fluct||=%.4e  (E-E0)/E0=%.2e"%(t,nv,de))
print("background energy density x box = %.4f ; E0 = %.4f"%(Eb,E0))
