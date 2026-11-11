import sys,json; sys.path.insert(0,'.')
import numpy as np, sympy as sp
from b410.fixed_eft import load_b3
# jets
names=['P','Q','S','E']  # delta phi, psi, sigma, eta
J={n:sp.symbols(f'{n} {n}t {n}x') for n in names}
w,q,v,r,gp,gs,M,f,sb,eb,eps=sp.symbols('w q v r gp gs M f sb eb eps')
def ab(ft,fx,eta):
    return (sp.cosh(eta)*ft+sp.sinh(eta)*fx)*(sp.sinh(eta)*ft+sp.cosh(eta)*fx)
Pt=w+eps*J['P'][1];Px=q+eps*J['P'][2];Qt=v+eps*J['Q'][1];Qx=r+eps*J['Q'][2]
sig=sb+eps*J['S'][0]; eta=eb+eps*J['E'][0]
L=sp.Rational(1,2)*(Pt**2-Px**2)+sp.Rational(1,2)*(Qt**2-Qx**2)+sp.Rational(1,2)*((eps*J['S'][1])**2-(eps*J['S'][2])**2)-sp.Rational(1,2)*M**2*sig**2+sp.Rational(1,2)*f*((eps*J['E'][1])**2-(eps*J['E'][2])**2)+sig*(gp*ab(Pt,Px,eta)+gs*ab(Qt,Qx,eta))
L2=sp.diff(L,eps,2).subs(eps,0)/2
z=[x for n in names for x in J[n]]
Hm=sp.hessian(L2,z)  # L2 = 1/2 z^T Hm z (L2 quadratic)
Hf=sp.lambdify((w,q,v,r,gp,gs,M,f,sb,eb),Hm,'numpy')
def K(vals,om,k):
    H=np.array(Hf(*vals),dtype=complex)
    T=np.zeros((12,4),dtype=complex)
    for i in range(4):
        T[3*i,i]=1;T[3*i+1,i]=-1j*om;T[3*i+2,i]=1j*k
    return T.conj().T@H@T
def roots(vals,X,Mv=1.0):
    k=X*Mv
    # det K as polynomial in c: sample & fit degree-8 poly (omega=c k)
    cs=np.linspace(-3,3,25); d=np.array([np.linalg.det(K(vals,c*k,k)).real for c in cs])
    co=np.polyfit(cs,d,8); rt=np.roots(co); return rt
b=load_b3(); man=json.load(open('data/e3_manifest.json'))['selected']
rows=[dict(p=b.p,s=b.s,kappa=b.kappa,w=b.w,q=b.q,v=b.v,r=b.r)]+man
out=[]
for d in rows:
    p,s,k,w_,q_,v_,r_=[d[x] for x in('p','s','kappa','w','q','v','r')]
    g1=np.sqrt(p); g2=-k/g1
    A=g1*w_*q_+g2*v_*r_; B=.5*(g1*(w_**2+q_**2)+g2*(v_**2+r_**2))
    t=-B/A
    if abs(t)>=1: continue
    e=0.5*np.arctanh(t); Jv=A*np.cosh(2*e)+B*np.sinh(2*e)
    out.append((d,e,Jv,g1,g2))
print("n backgrounds with stationary eta:",len(out))
res=[]
for fv in (0.1,1.0,10.0):
    fronts=[];cplx=0
    for d,e,Jv,g1,g2 in out:
        vals=(d['w'],d['q'],d['v'],d['r'],g1,g2,1.0,fv,Jv,e)
        rt=roots(vals,X=200.0)
        # imaginary part => hyperbolicity failure
        if np.abs(rt.imag).max()>1e-3*max(1,np.abs(rt).max()): cplx+=1
        fronts.append(np.abs(rt.real).max())
    print("f=%g: front (large X) min %.3f median %.3f max %.3f; front<=1+1e-6: %d/%d; non-real: %d"%(fv,min(fronts),np.median(fronts),max(fronts),sum(x<=1+1e-6 for x in fronts),len(fronts),cplx))
# analytic prediction: |sigma_bar|=|J(eb)| nonzero?
print("sigma_bar min |J| over stationary backgrounds: %.4g"%min(abs(o[2]) for o in out))
