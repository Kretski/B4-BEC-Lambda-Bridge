# Manifest generator for HeavyScalar_TwoLight (HS2). Draws only; no C3 evaluation. Z>0 and the IR quadratic F are analytic pre-classifications (Z positive definite is a C2 precondition).
import numpy as np, json, hashlib
M=1.0
CPL={'K1':(0.2,0.2,0.4),'K3':(0.3,-0.2,0.1),'K4':(0.15,0.25,0.2)}
def zbar(a,b,e,w,q,v,r):
    O=a*(w*w-q*q)+b*(v*v-r*r)+e*(w*v-q*r); ch=O/M**2
    Z=np.array([[1+2*a*ch,e*ch],[e*ch,1+2*b*ch]]); return ch,Z
def Fcoef(a,b,e,w,q,v,r):
    ch,Z=zbar(a,b,e,w,q,v,r); Cp=np.array([[a,e/2],[e/2,b]]); dZ=np.linalg.det(Z); adj=np.array([[Z[1,1],-Z[0,1]],[-Z[1,0],Z[0,0]]])
    # s = Cp (c*(w,v)+(q,r)) = c*s1 + s0 ; F(c) = M^2 (1-c^2) dZ - 4 (c s1+s0)^T adj (c s1+s0)
    s1=Cp@np.array([w,v]); s0=Cp@np.array([q,r])
    A2=-M**2*dZ-4*s1@adj@s1; A1=-8*s1@adj@s0; A0=M**2*dZ-4*s0@adj@s0
    return A2,A1,A0
rows=[]
for ci,(name,(a,b,e)) in enumerate(CPL.items()):
    rng=np.random.default_rng(31+ci)
    for n in range(300):
        w,q,v,r=rng.uniform(-1.5,1.5,4)
        ch,Z=zbar(a,b,e,w,q,v,r); ev=np.linalg.eigvalsh(Z)
        A2,A1,A0=Fcoef(a,b,e,w,q,v,r); disc=A1*A1-4*A2*A0
        rows.append(dict(cs=name,a=a,b=b,e=e,idx=n,w=float(w),q=float(q),v=float(v),r=float(r),M=M,chibar=float(ch),zmin=float(ev[0]),zmax=float(ev[1]),F_disc=float(disc),
                         zpos=bool(ev[0]>0),F_real=bool(disc>=0 and abs(A2)>0)))
H1=[x for x in rows if x['zpos']]                   # population for C3 (Z positive definite)
H2=[x for x in rows if not x['zpos']]               # ghost-type population (C2 expected FAIL: Z not positive definite)
H3=[x for x in rows if x['zpos'] and not x['F_real']][:50]   # IR-unstable control (analytic F complex roots), first 50 in order
c4={}
for cs in ('K3','K4'):
    c4[cs]=[x for x in H1 if x['cs']==cs and x['F_real']][0]
anchors=[dict(label='B2h',a=0.2,b=0.2,e=0.4,w=1.0,q=0.3,v=0.8,r=-0.2,M=1.0),dict(label='B2k',a=0.15,b=0.25,e=0.2,w=1.0,q=0.3,v=0.8,r=-0.2,M=1.0),
         dict(label='B2s',a=0.3,b=-0.2,e=0.1,w=1.0,q=0.3,v=0.8,r=-0.2,M=1.0),dict(label='B2z',a=0.2,b=0.2,e=0.0,w=1.0,q=0.3,v=0.8,r=-0.2,M=1.0)]
man=dict(H1=H1,H2=H2,H3_control=H3,anchors=anchors,C4_designated=c4)
s=json.dumps(man,sort_keys=True); open('t3fund_hs2_manifest.json','w').write(json.dumps(man,indent=1))
print('MANIFEST_HASH',hashlib.sha256(s.encode()).hexdigest()[:16])
for cs in CPL: print(cs,'raw 300; H1',sum(1 for x in H1 if x['cs']==cs),'H2',sum(1 for x in H2 if x['cs']==cs),'F_real in H1',sum(1 for x in H1 if x['cs']==cs and x['F_real']))
print('H3 control',len(H3)); print('C4',json.dumps(c4,indent=0))
