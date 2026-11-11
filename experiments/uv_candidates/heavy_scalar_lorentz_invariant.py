import sympy as sp, numpy as np
w,q,g,M,om,k,c=sp.symbols('w q g M om k c',real=True)
# fields: Phi,chi ; jets (f, f_t, f_x); L=1/2(dphi)^2+1/2(dchi)^2-1/2 M^2 chi^2 + g chi (dphi)^2   (signature +-)
P=sp.symbols('P Pt Px'); C=sp.symbols('C Ct Cx'); e=sp.symbols('e')
Xp=(w+e*P[1])**2-(q+e*P[2])**2
Chi=g*(w**2-q**2)/M**2+e*C[0]
L=sp.Rational(1,2)*Xp+sp.Rational(1,2)*((e*C[1])**2-(e*C[2])**2)-sp.Rational(1,2)*M**2*Chi**2+g*Chi*Xp
L2=sp.diff(L,e,2).subs(e,0)/2
z=list(P)+list(C); H=sp.hessian(L2,z)
T=sp.zeros(6,2); T[0,0]=1;T[1,0]=-sp.I*om;T[2,0]=sp.I*k;T[3,1]=1;T[4,1]=-sp.I*om;T[5,1]=sp.I*k
K=sp.simplify((T.H*H*T)); D=sp.factor(sp.simplify(K.det().subs(om,c*k)))
print(D)
import numpy as np
def roots(wv,qv,gv,Mv,X):
    d=sp.Poly(sp.expand(D.subs({w:wv,q:qv,g:gv,M:Mv,k:X*Mv})),c); return sorted(np.roots([complex(a) for a in d.all_coeffs()]),key=lambda z:z.real)
for (wv,qv,gv) in [(1,0.3,0.4),(1,0.3,-0.4),(0.3,1,0.4),(0.3,1,-0.4)]:
    print('w,q,g=',wv,qv,gv)
    for X in (1e-3,0.1,1,10,1e3):
        r=roots(wv,qv,gv,1.0,X); print('  X=%g'%X,[np.round(z,4) for z in r])
