# T3-fund SM + dynamical eta: population manifest (M=1). Coupling sets, background draws and analytic pre-classification (closed forms of the spec, section 2).
import numpy as np, json, hashlib
SETS={'G1':(0.3,-0.2,41),'G2':(0.25,-0.35,42),'G3':(0.8,-0.6,43)}; NRAW=200
FMAIN=[0.1,1.0,10.0]
def k0_class(g,f,al,be,A,B,sig,Zt):
    """closed-form k=0 class: det K(omega,0)=omega^4 Zt1 Zt2 [(x-1-S1)(f x+4 sig^2-S2)-S3^2], x=omega^2. Returns ('real'|'complex'), roots x."""
    N=[2*g[i]*(A*al[i]-B*be[i]) for i in range(2)]; m=[g[i]*(A*be[i]-B*al[i])/sig for i in range(2)]
    S1=sum(m[i]**2/Zt[i] for i in range(2)); S2=sum(N[i]**2/Zt[i] for i in range(2)); S3=sum(m[i]*N[i]/Zt[i] for i in range(2))
    a2=f; a1=(4*sig**2-S2)-f*(1+S1); a0=-(1+S1)*(4*sig**2-S2)-S3**2
    disc=a1*a1-4*a2*a0
    if disc<0: return 'complex',None
    x=sorted([(-a1-np.sqrt(disc))/(2*a2),(-a1+np.sqrt(disc))/(2*a2)])
    return ('real' if x[0]>=0 else 'complex'),[float(t) for t in x]
def build():
    man={'sets':{k:dict(gphi=v[0],gpsi=v[1],seed=v[2]) for k,v in SETS.items()},'nraw':NRAW,'fmain':FMAIN,'draws':[]}
    for cs,(g1,g2,seed) in SETS.items():
        rng=np.random.default_rng(seed); X=rng.uniform(-1.5,1.5,(NRAW,4))
        for idx,(w,q,v,r) in enumerate(X):
            g=[g1,g2]; al=[w,v]; be=[q,r]
            A=g1*w*q+g2*v*r; B=0.5*(g1*(w*w+q*q)+g2*(v*v+r*r))
            d=dict(cs=cs,idx=idx,gphi=g1,gpsi=g2,w=float(w),q=float(q),v=float(v),r=float(r),A=float(A),B=float(B))
            if abs(B)<abs(A):
                sig=float(np.sign(A)*np.sqrt(A*A-B*B)); eta=float(0.5*np.arctanh(-B/A)); d['branch']='S'
                Zt=[1-g[i]*B for i in range(2)]; Xg=[1+g[i]*B for i in range(2)]
                d.update(sigma=sig,eta=eta,Zt=Zt,Xg=Xg,Y=[g[i]*A for i in range(2)])
                if min(Zt)<=0: d['pop']='S2g'
                elif min(Xg)<=0: d['pop']='S2x'
                else: d['pop']='S1'
                if min(Zt)>0:
                    d['front']=float(max(abs(g[i]*A+s*np.sqrt(1+(g[i]*sig)**2))/Zt[i] for i in range(2) for s in (1,-1)))
                    d['k0class_f1']=k0_class(g,1.0,al,be,A,B,sig,Zt)[0]
            elif abs(A)<abs(B):
                eta=float(0.5*np.arctanh(-A/B)); d.update(branch='Z',sigma=0.0,eta=eta,pop='Z0')
            else: d['branch']='boundary'; d['pop']='excluded'
            man['draws'].append(d)
    # anchors: fixed background, G1 and G2, f=1; must be S1
    man['anchors']=[]
    for cs in ('G1','G2','G3'):
        g1,g2,_=SETS[cs]; w,q,v,r=1.0,0.3,0.8,-0.2; A=g1*w*q+g2*v*r; B=0.5*(g1*(w*w+q*q)+g2*(v*v+r*r))
        man['anchors'].append(dict(label='A_'+cs,cs=cs,gphi=g1,gpsi=g2,w=w,q=q,v=v,r=r,f=1.0,S_branch=bool(abs(B)<abs(A))))
    S1=[d for d in man['draws'] if d.get('pop')=='S1']
    man['H3_control_fneg']=[dict(cs=d['cs'],idx=d['idx'],f=-1.0) for d in S1[:40]]
    man['H3_boundary_f0']=[dict(cs=d['cs'],idx=d['idx'],f=0.0) for d in S1[:10]]
    # C4 designated: the first two S1 draws in manifest order (across sets G1,G2,G3) whose closed-form k=0 class at f=1 is real
    man['C4_designated']=[]
    for d in S1:
        if d.get('k0class_f1')=='real' and len(man['C4_designated'])<2:
            man['C4_designated'].append(dict(cs=d['cs'],idx=d['idx'],f=1.0,w=d['w'],q=d['q'],v=d['v'],r=d['r'],gphi=d['gphi'],gpsi=d['gpsi'],sigma=d['sigma'],eta=d['eta']))
    return man
if __name__=='__main__':
    man=build(); s=json.dumps(man,sort_keys=True,separators=(',',':')); h=hashlib.sha256(s.encode()).hexdigest()[:16]
    json.dump(man,open('t3fund_smeta_manifest.json','w'),indent=1,sort_keys=True)
    print('MANIFEST_HASH',h)
    from collections import Counter
    print(Counter((d['cs'],d.get('pop')) for d in man['draws']))
    print('anchors',man['anchors']); print('C4',man['C4_designated'])
    print('S1 k0class f=1:',Counter((d['cs'],d.get('k0class_f1')) for d in man['draws'] if d.get('pop')=='S1'))
    print('S1 front>1+1e-6:',Counter((d['cs'],d['front']>1+1e-6) for d in man['draws'] if d.get('pop')=='S1'))
