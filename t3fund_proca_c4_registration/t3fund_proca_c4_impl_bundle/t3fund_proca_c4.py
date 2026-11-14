# T3-fund Proca C4 (time domain). Implements T3_FUND_PROCA_SPEC.md section 6 (registered in commit e355af2; SPEC_HASH 86c360326da7406f).
# Reduced first-order system of section 3(iii) (A0 eliminated by the Gauss constraint; F_xx by Fourier multiplier), periodic box 1024, pseudo-spectral D, RK4, T=1000, records every 5.
# State: (a=f_t, b=f_x, c=g_t, d=g_x, A1, F) plus the fluctuation fields (df, dg) with d(df)/dt = a - w, d(dg)/dt = c - v (they do not feed back; needed for the norm of section 6).
import numpy as np, json, time, argparse, os
Lbox = 1024.0; s = 15.0; k0p = 0.5; EPS_REL = 1e-3; T = 1000.0; REC = 5.0

def make_ops(N, L):
    kk = 2*np.pi*np.fft.rfftfreq(N, L/N)
    return kk
def rhs(S, th, kk, M=1.0):
    """S shape (8,N): a,b,c,d,A1,F,df,dg. th: dict gphi,gpsi,Lam,w,v. Returns dS/dt. Uses D = i k (rfft)."""
    gp, gs, Lam = th['gphi'], th['gpsi'], th['Lambda']; N = S.shape[1]
    a, b, c, d, A, F = S[:6]
    fa, fb, fc, fd, fA, fF = (np.fft.rfft(u) for u in (a, b, c, d, A, F))
    D = lambda fu: np.fft.irfft(1j*kk*fu, N)
    Dxx_F = np.fft.irfft(-kk**2*fF, N)
    A1t = F - Dxx_F/M**2
    Db, Da, Dd, Dc, DA = D(fb), D(fa), D(fd), D(fc), D(fA)
    da = (1 + 3*Lam*b**2)*Db - gp*(A1t*b + DA*a + 2*A*Da)
    dc = Dd - gs*(A1t*d + DA*c + 2*A*Dc)
    dF = gp*a*b + gs*c*d - M**2*A
    return np.array([da, Da, dc, Dc, A1t, dF, a - th['w'], c - th['v']])
def energy(S, th, kk, dx, M=1.0):
    Lam = th['Lambda']; N = S.shape[1]; a, b, c, d, A, F = S[:6]
    DF = np.fft.irfft(1j*kk*np.fft.rfft(F), N)
    return dx*np.sum(0.5*a**2 + 0.5*b**2 + Lam/4*b**4 + 0.5*c**2 + 0.5*d**2 + 0.5*F**2 + DF**2/(2*M**2) + 0.5*M**2*A**2)
def background_energy(th, L, M=1.0):
    Abar = th['Abar']
    return L*(0.5*th['w']**2 + 0.5*th['q']**2 + th['Lambda']/4*th['q']**4 + 0.5*th['v']**2 + 0.5*th['r']**2 + 0.5*M**2*Abar**2)

def theta_from_case(c, coupled):
    th = dict(gphi=c['gphi'], gpsi=c['gpsi'], Lambda=c['Lambda'], w=c['w'], q=c['q'], v=c['v'], r=c['r'])
    if not coupled: th['gphi'] = th['gpsi'] = 0.0
    th['Abar'] = th['gphi']*th['w']*th['q'] + th['gpsi']*th['v']*th['r']        # M=1
    return th
def initial_state(th, N, L=Lbox):
    x = np.linspace(-L/2, L/2, N, endpoint=False); dx = x[1]-x[0]; kk = make_ops(N, L)
    D = lambda f: np.fft.irfft(1j*kk*np.fft.rfft(f), N)
    eps = EPS_REL*max(abs(th['q']), abs(th['r']))
    gau = eps*np.exp(-x**2/(2*s*s))*np.cos(k0p*x)
    dfl, dgl = gau.copy(), gau.copy(); dA = gau.copy()
    S = np.array([th['w'] - D(dfl), th['q'] + D(dfl), th['v'] + D(dgl), th['r'] + D(dgl), th['Abar'] + dA, np.zeros(N), dfl, dgl])
    return S, x, dx, kk, eps
def norm_of(S, th, dx):
    return np.sqrt(dx*np.sum(S[6]**2 + S[7]**2 + (S[4]-th['Abar'])**2))

def run(case, coupled, N, dt, T=T, ckpt=None, ckpt_every=4000):
    th = theta_from_case(case, coupled); S, x, dx, kk, eps = initial_state(th, N)
    n = int(round(T/dt)); every = int(round(REC/dt)); i0 = 0
    ts = [0.]; Es = [energy(S, th, kk, dx)]; nr = [norm_of(S, th, dx)]
    if ckpt and os.path.exists(ckpt):
        z = np.load(ckpt); S = z['S']; i0 = int(z['i']); ts = list(z['ts']); Es = list(z['Es']); nr = list(z['nr'])
    t0 = time.time()
    for i in range(i0+1, n+1):
        k1 = rhs(S, th, kk); k2 = rhs(S+.5*dt*k1, th, kk); k3 = rhs(S+.5*dt*k2, th, kk); k4 = rhs(S+dt*k3, th, kk); S = S + dt/6*(k1+2*k2+2*k3+k4)
        if i % every == 0:
            ts.append(i*dt); Es.append(energy(S, th, kk, dx)); nr.append(norm_of(S, th, dx))
        if ckpt and i % ckpt_every == 0:
            np.savez(ckpt+'.tmp.npz', S=S, i=i, ts=np.array(ts), Es=np.array(Es), nr=np.array(nr)); os.replace(ckpt+'.tmp.npz', ckpt)
    ts, Es, nr = map(np.array, (ts, Es, nr)); Ebg = background_energy(th, Lbox); Efl = Es - Ebg
    finite = bool(np.all(np.isfinite(S)) and np.all(np.isfinite(Es)) and np.all(np.isfinite(nr)))
    res = dict(coupled=bool(coupled), N=N, dt=dt, eps=float(eps), n_records=len(ts), t_end=float(ts[-1]), finite=finite, E0=float(Es[0]), Ebg=float(Ebg), Efl0=float(Efl[0]),
               drift_total_max=float(np.abs(Es-Es[0]).max()/abs(Es[0])), drift_total_final=float((Es[-1]-Es[0])/Es[0]), drift_fl_max=float(np.abs(Efl-Efl[0]).max()/abs(Efl[0])),
               norm0=float(nr[0]), norm_ratio_max=float(nr.max()/nr[0]), norm_ratio_final=float(nr[-1]/nr[0]))
    return res, ts, nr, Es
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--case'); ap.add_argument('--manifest'); ap.add_argument('--coupled', type=int); ap.add_argument('--N', type=int); ap.add_argument('--dt', type=float)
    ap.add_argument('--out'); ap.add_argument('--T', type=float, default=T); a = ap.parse_args()
    man = json.load(open(a.manifest)); allc = man['PE'] + man['PN']
    case = [e for e in allc if e['label'] == a.case.split('/')[0] and (len(a.case.split('/')) < 2 or e['mode'] == a.case.split('/')[1])]; assert len(case) == 1, len(case); case = case[0]
    tag = 'c4proca_%s_%s_N%d' % (case['label'] + ('_' + case['mode'] if case['pop'] == 'PE' else ''), 'coupled' if a.coupled else 'control', a.N); t0 = time.time()
    res, ts, nr, Es = run(case, a.coupled, a.N, a.dt, a.T, ckpt='%s/%s.ckpt.npz' % (a.out, tag))
    res['seconds'] = time.time()-t0; res['case'] = case['label']; res['mode'] = case['mode']
    json.dump(res, open('%s/%s.json' % (a.out, tag), 'w'), indent=1); np.savetxt('%s/%s.csv' % (a.out, tag), np.c_[ts, nr, Es], header='t norm_fluct energy', delimiter=',')
    print(tag, json.dumps(res), flush=True)
