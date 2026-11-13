# Algebra identity check (no population, no verdict): linearise the Euler-Lagrange equations of the HS2 Lagrangian about the exact background with sympy,
# build the 3x3 symbol K(omega,k) and compare det K with the closed form  k^6 * D * [ (D+eps M^2) D detZ - 4 eps s^T adj(Z) s ].
import sympy as sp, random
t,x=sp.symbols('t x'); a,b,e,M,w,q,v,r,om,k,ep=sp.symbols('a b e M w q v r omega k epsilon')
phi,psi,chi=[sp.Function(n)(t,x) for n in ('phi','psi','chi')]
d=lambda f,g:sp.diff(f,t)*sp.diff(g,t)-sp.diff(f,x)*sp.diff(g,x)
L=sp.Rational(1,2)*(d(phi,phi)+d(psi,psi)+d(chi,chi))-M**2*chi**2/2+chi*(a*d(phi,phi)+b*d(psi,psi)+e*d(phi,psi))
def EL(f): return sp.diff(sp.diff(L,sp.diff(f,t)),t)+sp.diff(sp.diff(L,sp.diff(f,x)),x)-sp.diff(L,f)
eqs=[EL(phi),EL(psi),EL(chi)]
chib=(a*(w*w-q*q)+b*(v*v-r*r)+e*(w*v-q*r))/M**2
al,be,xi=[sp.Function(n)(t,x) for n in ('al','be','xi')]
ee=sp.symbols('ee'); sub={phi:w*t+q*x+ee*al,psi:v*t+r*x+ee*be,chi:chib+ee*xi}
lin=[sp.diff(E.subs(sub).doit(),ee).subs(ee,0) for E in eqs]
A0,B0,X0=sp.symbols('A0 B0 X0'); pw=sp.exp(sp.I*(k*x-om*t))
K=sp.zeros(3,3)
for i,E in enumerate(lin):
    ex=E.subs({al:A0*pw,be:B0*pw,xi:X0*pw}).doit()/pw
    ex=sp.expand(sp.simplify(ex))
    for j,s in enumerate((A0,B0,X0)): K[i,j]=ex.coeff(s)
Cp=sp.Matrix([[a,e/2],[e/2,b]]); Zb=sp.eye(2)+2*chib*Cp
c=sp.symbols('c'); D=1-c**2
sv=Cp*sp.Matrix([c*w+q,c*v+r]); adj=Zb.adjugate()
closed=D*((D+ep*M**2)*D*Zb.det()-4*ep*(sv.T*adj*sv)[0])
detK=K.det()
vals={a:.21,b:-.13,e:.37,M:1.3,w:.9,q:-.4,v:.55,r:.7,k:1.7}
cv=0.83; omv=cv*vals[k]
lhs=complex(detK.subs(vals).subs(om,omv).evalf())
rhs=complex((vals[k]**6*closed).subs(vals).subs({c:cv,ep:1/vals[k]**2}).evalf())
print('det K',lhs,' closed',rhs,' relerr',abs(lhs-rhs)/abs(rhs))
random.seed(1); mx=0
for _ in range(20):
    vv={a:random.uniform(-.5,.5),b:random.uniform(-.5,.5),e:random.uniform(-.5,.5),M:random.uniform(.5,2),w:random.uniform(-1,1),q:random.uniform(-1,1),v:random.uniform(-1,1),r:random.uniform(-1,1),k:random.uniform(.1,5)}
    cv=random.uniform(-2,2); l=complex(detK.subs(vv).subs(om,cv*vv[k]).evalf()); rr=complex((vv[k]**6*closed).subs(vv).subs({c:cv,ep:1/vv[k]**2}).evalf())
    mx=max(mx,abs(l-rr)/abs(rr))
print('max relerr over 20 random',mx)
# k=0 polynomial
K0=K.subs(k,0); P0=sp.factor(sp.simplify(K0.det())); print('k=0 det:',P0)
