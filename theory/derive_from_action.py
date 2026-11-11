"""Derive A_i from the Proca action coefficients (SymPy), then check against b410.fixed_eft.
Action (1+1D, signature -+):
 L = 1/2 f_t^2 - 1/2 f_x^2 - (Lam/4) f_x^4 + 1/2 g_t^2 - 1/2 g_x^2
     + 1/2 (A1_t - A0_x)^2 + 1/2 M^2 (A0^2 - A1^2) + A1 (gphi f_t f_x + gpsi g_t g_x)
Background: f = w t + q x, g = v t + r x, A0 = 0, A1 = Abar = (gphi w q + gpsi v r)/M^2.
"""
import sympy as sp
t,x=sp.symbols('t x',real=True)
M,gphi,gpsi,w,q,v,r,Lam,e=sp.symbols('M gphi gpsi w q v r Lam e')
k,om,c,eps_=sp.symbols('k omega c eps')
f=sp.Function('f')(t,x); g=sp.Function('g')(t,x); a0=sp.Function('a0')(t,x); a1=sp.Function('a1')(t,x)
Abar=(gphi*w*q+gpsi*v*r)/M**2
F_=w*t+q*x+e*f; G_=v*t+r*x+e*g; A1=Abar+e*a1; A0=e*a0
L=(sp.Rational(1,2)*sp.diff(F_,t)**2-sp.Rational(1,2)*sp.diff(F_,x)**2-Lam/4*sp.diff(F_,x)**4
   +sp.Rational(1,2)*sp.diff(G_,t)**2-sp.Rational(1,2)*sp.diff(G_,x)**2
   +sp.Rational(1,2)*(sp.diff(A1,t)-sp.diff(A0,x))**2+sp.Rational(1,2)*M**2*(A0**2-A1**2)
   +A1*(gphi*sp.diff(F_,t)*sp.diff(F_,x)+gpsi*sp.diff(G_,t)*sp.diff(G_,x)))
L2=sp.expand(sp.diff(L,e,2).subs(e,0)/2)
# linear terms must vanish on the background (check)
L1=sp.expand(sp.diff(L,e,1).subs(e,0))
print("L2 =",sp.simplify(L2))
eqs=sp.calculus.euler.euler_equations(L2,[f,g,a0,a1],[t,x])
amp=sp.symbols('Uf Ug Ua0 Ua1'); ph=sp.exp(sp.I*(k*x-om*t))
sub={f:amp[0]*ph,g:amp[1]*ph,a0:amp[2]*ph,a1:amp[3]*ph}
K=sp.Matrix([[sp.simplify(sp.expand(eq.lhs.subs(sub).doit()/ph).coeff(a)) for a in amp] for eq in eqs])
D=sp.factor(sp.simplify(K.det().subs(om,c*k).subs(k,M*sp.Symbol('x_'))))
print("det K(omega=c k, k=M x) =",D)
p,s,kap=gphi**2/M**2,gpsi**2/M**2,-gphi*gpsi/M**2
aphi=p*w*q-kap*v*r; apsi=-kap*w*q+s*v*r
d1=c**2-2*aphi*c-(1+3*Lam*q**2); d2=c**2-2*apsi*c-1
Fc=d1*d2+p*(w-q*c)**2*d2+s*(v-r*c)**2*d1
X=sp.Symbol('x_'); epsx=X**2/(1+X**2)
Q=Fc-epsx*c**2*d1*d2
ratio=sp.simplify(D/Q)
print("det K / Q =",ratio)
# A_i: implicit differentiation  Q(c,eps)=0  ->  dc/deps|0 = c^2 d1 d2 / F'(c);  eps = x^2 + O(x^4)
Ai=c**2*d1*d2/sp.diff(Fc,c)
print("A_i(c) =",Ai)
import pickle; pickle.dump((D,Q,Ai),open('/tmp/val/sym.pkl','wb'))
