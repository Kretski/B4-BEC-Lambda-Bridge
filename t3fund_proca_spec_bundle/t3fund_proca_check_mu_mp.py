"""50-digit verification of the C1(c) bound gamma(k) <= mu(k) of T3_FUND_PROCA_SPEC.md section 3(iv), at the registered C3 grid points.
Matrices (M=1): Mlin(k) = k*Ps + khat*Pk, khat=sqrt(1+k^2), Pinf=Ps+Pk; V = eigenvectors of Pinf (50 digits), Rr = V^-1 [(khat-k) Pk] V,
mu = lambda_max(Herm(Rr)), gamma = max Re eig Mlin(k).  Check gamma <= mu*(1+1e-9)+1e-9 at all 92 grid points k/M of the C3 grid.
Backgrounds: 12 random (seed 20261114) with strictly hyperbolic Pinf (real simple eigenvalues, separation>1e-6 in double pre-screen, confirmed in 50 digits).
This is a verification of the implementation on random backgrounds, NOT an evaluation of the registered populations."""
import numpy as np, mpmath as mp, json, time
mp.mp.dps = 50
KGRID = np.concatenate([np.geomspace(1e-3, 1e3, 61), np.geomspace(1e3, 1e6, 32)[1:]])
def mats(gp, gs, Lam, w, q, v, r):
    I = mp.mpc(0, 1); Ab = gp*w*q + gs*v*r
    Ps = mp.zeros(6, 6); Pk = mp.zeros(6, 6)
    Ps[0, 1] = I*(1 + 3*Lam*q**2); Ps[0, 0] = -2*I*gp*Ab; Ps[0, 4] = -gp*I*w
    Ps[1, 0] = I; Ps[2, 3] = I; Ps[2, 2] = -2*I*gs*Ab; Ps[2, 4] = -gs*I*v; Ps[3, 2] = I
    Pk[0, 5] = -gp*q; Pk[2, 5] = -gs*r; Pk[4, 5] = 1
    Pk[5, 0] = gp*q; Pk[5, 1] = gp*w; Pk[5, 2] = gs*r; Pk[5, 3] = gs*v; Pk[5, 4] = -1
    return Ps, Pk
rng = np.random.default_rng(20261114); out = []; t0 = time.time(); done = 0
while done < 12:
    gp, gs = [mp.mpf(float(z)) for z in rng.uniform(-1, 1, 2)]; Lam = mp.mpf(float(rng.uniform(0, 1))); w, q, v, r = [mp.mpf(float(z)) for z in rng.uniform(-1.2, 1.2, 4)]
    Ps, Pk = mats(gp, gs, Lam, w, q, v, r); Pinf = Ps + Pk
    lam, V = mp.eig(Pinf)
    cspeeds = [(mp.mpc(0, 1)*l) for l in lam]            # eigenvalue of Pinf is -i c, so c = i*lambda
    if max(abs(mp.im(z)) for z in cspeeds) > mp.mpf('1e-30'): continue
    sep = min(abs(lam[i] - lam[j]) for i in range(6) for j in range(i))
    if sep < mp.mpf('1e-6'): continue
    Vi = mp.inverse(V); worst_ratio = -mp.inf; maxgam = mp.mpf(0); nviol = 0
    for k in KGRID:
        kk = mp.mpf(float(k)); kh = mp.sqrt(1 + kk*kk); Mx = kk*Ps + kh*Pk
        ev = mp.eig(Mx, left=False, right=False); gam = max(mp.re(e) for e in ev)
        Rr = Vi*((kh - kk)*Pk)*V; H = (Rr + Rr.H)/2; mu = max(mp.eigh(H, eigvals_only=True))
        if gam > mu*(1 + mp.mpf('1e-9')) + mp.mpf('1e-9'): nviol += 1
        maxgam = max(maxgam, gam)
    cond = mp.mnorm(V, 1)*mp.mnorm(Vi, 1)
    out.append(dict(n_grid=len(KGRID), violations=nviol, max_gamma=mp.nstr(maxgam, 5), cond_V_inf=mp.nstr(cond, 5), sep_inf=mp.nstr(sep, 5)))
    done += 1; print(done, out[-1], round(time.time() - t0), 's', flush=True)
print('SUMMARY: backgrounds', len(out), 'grid points each', len(KGRID), 'total violations', sum(o['violations'] for o in out))
json.dump(out, open('proca_check_mu_mp.json', 'w'), indent=1)
