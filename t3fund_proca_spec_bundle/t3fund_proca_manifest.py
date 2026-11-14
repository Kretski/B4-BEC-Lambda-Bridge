"""Population manifest for the T3-fund evaluation of the Proca mediator (see T3_FUND_PROCA_SPEC.md).
Prescription (fixed-EFT, rank-one): M=1, g_phi = sqrt(p), g_psi = -kappa/sqrt(p), kappa = -/+ sqrt(p s) (rank1_minus / rank1_plus).
PE : 40 backgrounds of data/e3_manifest.json (sha256 recorded) + the B3 anchor (data/B3_max_speed_7_315096.json).
PN : new draws, generator identical in distribution to vendor_t173/robustness_scan.generate_samples (p,s log-uniform on [0.05,1], Lambda U[0,1],
     X_q,X_r,X_w,X_v U[0,cut], random signs); cuts (0.01,0.04,0.09,0.25,0.95) x modes (rank1_minus, rank1_plus), N=1000 per cell,
     seed = 2026111400 + 10*cut_index + mode_index (numpy default_rng).
H3 : ghost control, first 40 PN draws of cell (cut 0.09, rank1_minus), evaluated with the psi kinetic term sign-flipped (time-derivative Hessian entry -1).
"""
import json, hashlib, numpy as np
CUTS = (0.01, 0.04, 0.09, 0.25, 0.95); MODES = (('rank1_minus', -1.0), ('rank1_plus', 1.0)); N = 1000
def sha(path): return hashlib.sha256(open(path, 'rb').read()).hexdigest()
def gen(n, cut, seed):
    rng = np.random.default_rng(seed)
    for i in range(n):
        p = float(np.exp(rng.uniform(np.log(0.05), np.log(1.0)))); s = float(np.exp(rng.uniform(np.log(0.05), np.log(1.0)))); Lam = float(rng.uniform(0.0, 1.0))
        Xq, Xr, Xw, Xv = rng.uniform(0, cut, 4); sg = rng.choice([-1.0, 1.0], 4)
        yield dict(index=i, p=p, s=s, Lambda=Lam, q=float(sg[0]*np.sqrt(Xq/p)), r=float(sg[1]*np.sqrt(Xr/s)), w=float(sg[2]*np.sqrt(Xw/p)), v=float(sg[3]*np.sqrt(Xv/s)), Xq=float(Xq), Xr=float(Xr), Xw=float(Xw), Xv=float(Xv))
def couplings(p, s, kappa):
    return float(np.sqrt(p)), float(-kappa/np.sqrt(p))
def build():
    man = dict(version=1, M=1.0, prescription='g_phi=sqrt(p), g_psi=-kappa/sqrt(p), kappa=-/+sqrt(p s) (rank-one); spectrum depends on x=k/M only',
               cuts=list(CUTS), modes=[m for m, _ in MODES], N_per_cell=N, e3_source_sha256=sha('e3_manifest_source.json'), anchor_source_sha256=sha('B3_anchor_source.json'))
    pe = []
    e3 = json.load(open('e3_manifest_source.json'))
    for e in e3['selected']:
        gp, gs = couplings(e['p'], e['s'], e['kappa'])
        pe.append(dict(pop='PE', label='E3_%s_%d' % (e['cls'], e['index']), mode=e['mode'], cls=e['cls'], p=e['p'], s=e['s'], kappa=e['kappa'], Lambda=e['Lambda'], w=e['w'], q=e['q'], v=e['v'], r=e['r'], gphi=gp, gpsi=gs))
    a = json.load(open('B3_anchor_source.json')); gp, gs = couplings(a['p'], a['s'], a['kappa'])
    pe.append(dict(pop='PE', label='B3_anchor', mode=a['mode'], cls='anchor', p=a['p'], s=a['s'], kappa=a['kappa'], Lambda=a['Lambda'], w=a['w'], q=a['q'], v=a['v'], r=a['r'], gphi=gp, gpsi=gs))
    pn = []
    for ci, cut in enumerate(CUTS):
        for mi, (mode, sg) in enumerate(MODES):
            for d in gen(N, cut, 2026111400 + 10*ci + mi):
                kap = sg*float(np.sqrt(d['p']*d['s'])); gp, gs = couplings(d['p'], d['s'], kap)
                pn.append(dict(pop='PN', label='PN_cut%g_%s_%d' % (cut, mode, d['index']), cut=cut, mode=mode, p=d['p'], s=d['s'], kappa=kap, Lambda=d['Lambda'], w=d['w'], q=d['q'], v=d['v'], r=d['r'], gphi=gp, gpsi=gs,
                               Xmax=max(d['Xq'], d['Xr'], d['Xw'], d['Xv']), omegaP2_k0=1 + d['p']*d['q']**2 + d['s']*d['r']**2))
    h3 = [dict(e, pop='H3', label=e['label'].replace('PN_', 'H3_'), kinetic_sign_psi=-1.0) for e in [x for x in pn if x['cut'] == 0.09 and x['mode'] == 'rank1_minus'][:40]]
    man.update(PE=pe, PN=pn, H3=h3)
    man['counts'] = dict(PE=len(pe), PN=len(pn), H3=len(h3), total=len(pe)+len(pn)+len(h3))
    return man
if __name__ == '__main__':
    man = build(); json.dump(man, open('t3fund_proca_manifest.json', 'w'), indent=1, sort_keys=True)
    h = hashlib.sha256(json.dumps(man, sort_keys=True, separators=(',', ':')).encode()).hexdigest()[:16]
    print(man['counts'], 'MANIFEST_HASH', h)
