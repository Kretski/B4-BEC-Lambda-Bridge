import sys, numpy as np
sys.path.insert(0,"/mnt/user-data/outputs")
import t3fund_stage1 as T
rng=np.random.default_rng(5)
# (A) explicit matrix K(omega,k) from linearised equations
for _ in range(3):
    gp,gs=rng.uniform(-1,1,2); w,q,v,r=rng.uniform(-1.5,1.5,4); k=rng.uniform(0.3,3); c=complex(rng.uniform(-2,2),rng.uniform(-1,1)); om=c*k; M=1.
    sb=gp*w*q+gs*v*r
    K=np.array([[k*k-om*om+2*gp*sb*om*k,0,1j*gp*(k*w-om*q)],
                [0,k*k-om*om+2*gs*sb*om*k,1j*gs*(k*v-om*r)],
                [-1j*gp*(w*k-q*om),-1j*gs*(v*k-r*om),k*k-om*om+M]],dtype=complex)
    d=np.linalg.det(K)/k**6
    p=T.poly('A',[gp,gs,w,q,v,r],k); pv=sum(co*c**i for i,co in enumerate(p))
    print('A',abs(d-pv)/abs(pv))
# k=0 polynomial
gp,gs=0.4,-0.7; w,q,v,r=0.3,1.1,-0.5,0.8
for om in (0.7,1.9+0.3j):
    K=np.array([[-om*om,0,1j*gp*(-om*q)],[0,-om*om,1j*gs*(-om*r)],[-1j*gp*(-q*om),-1j*gs*(-r*om),-om*om+1]],dtype=complex)
    co=T.k0_check('A',[gp,gs,w,q,v,r]); pv=sum(float(x)*om**i for i,x in enumerate(co)); print('A k0',abs(np.linalg.det(K)-pv)/abs(pv))
# (B) explicit
for _ in range(3):
    g=rng.uniform(-1,1); w,q=rng.uniform(-1.5,1.5,2); k=rng.uniform(.3,3); c=complex(rng.uniform(-2,2),rng.uniform(-1,1)); om=c*k
    Z=1+2*g*g*(w*w-q*q); s=w*om+q*k
    K=np.array([[Z*(k*k-om*om),-2j*g*s],[2j*g*s,k*k-om*om+1]],dtype=complex)
    d=np.linalg.det(K)/k**4; p=T.poly('B',[g,w,q],k); pv=sum(co*c**i for i,co in enumerate(p)); print('B',abs(d-pv)/abs(pv))
g,w,q=0.4,0.9,0.5; Z=1+2*g*g*(w*w-q*q)
for om in (0.7,1.9+0.3j):
    K=np.array([[Z*(-om*om),-2j*g*w*om],[2j*g*w*om,-om*om+1]],dtype=complex); co=T.k0_check('B',[g,w,q]); pv=sum(float(x)*om**i for i,x in enumerate(co)); print('B k0',abs(np.linalg.det(K)-pv)/abs(pv))
