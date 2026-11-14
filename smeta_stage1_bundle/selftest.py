import sys,json,time,numpy as np
sys.path.insert(0,'.')
import t3fund_smeta_stage1 as S
import sympy as sp
# sympy Hessian (independent of closed form)
names=['P','Q','Sg','E']; J={n:sp.symbols(f'{n} {n}t {n}x') for n in names}; gp,gs,f,Mm=sp.symbols('gp gs f M',real=True)
ab=lambda ft,fx,eta:(sp.cosh(eta)*ft+sp.sinh(eta)*fx)*(sp.sinh(eta)*ft+sp.cosh(eta)*fx)
L=sp.Rational(1,2)*(J['P'][1]**2-J['P'][2]**2)+sp.Rational(1,2)*(J['Q'][1]**2-J['Q'][2]**2)+sp.Rational(1,2)*(J['Sg'][1]**2-J['Sg'][2]**2)-sp.Rational(1,2)*Mm**2*J['Sg'][0]**2+sp.Rational(1,2)*f*(J['E'][1]**2-J['E'][2]**2)+J['Sg'][0]*(gp*ab(J['P'][1],J['P'][2],J['E'][0])+gs*ab(J['Q'][1],J['Q'][2],J['E'][0]))
allj=[x for n in names for x in J[n]]; Hf=sp.lambdify(allj+[gp,gs,f,Mm],sp.hessian(L,allj),'numpy')
man=json.load(open('t3fund_smeta_manifest.json')); P=S.build(man)
import math
for pop in ('S1','Z0'):
    c=[p for p in P if p['pop']==pop][3]; g=[c['gphi'],c['gpsi']]
    bg=S.background(g,c['w'],c['q'],c['v'],c['r'],math,1.0); H=np.array(S.hessian(g,c['f'],bg,math,1.0))
    vals=[0,c['w'],c['q'],0,c['v'],c['r'],bg['sig'],0,0,bg['eta'],0,0]
    Hs=np.array(Hf(*vals,g[0],g[1],c['f'],1.0),dtype=float); print(pop,'Hessian closed form vs sympy:',abs(H-Hs).max())
    # poly vs det K
    worst=0
    for k,cc in ((0.3,0.4),(1.1,-0.7),(2.5,1.6),(7.0,0.2)):
        om=cc*k; T=np.zeros((12,4),complex)
        for i in range(4): T[3*i,i]=1;T[3*i+1,i]=-1j*om;T[3*i+2,i]=1j*k
        Kx=T.conj().T@Hs@T; d=np.linalg.det(Kx)
        co=S.poly_c(H.tolist(),k,True); pv=sum(co[i]*cc**i for i in range(len(co)))*k**8
        worst=max(worst,abs(d-pv)/max(abs(pv),1e-300))
    print(pop,'det K vs poly:',worst)
t=time.time(); c=[p for p in P if p['pop']=='S1' and p['cs']=='G3'][0]; r=S.evaluate(c); print(time.time()-t)
print({k:v for k,v in r.items() if k not in ('P6',)})
