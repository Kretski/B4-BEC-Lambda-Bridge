# T3-fund HS2 C4 (time domain).  Implements T3_FUND_HS2_SPEC.md section 5 (registered, commit 7b0e867, SPEC_HASH 9d52c927a9acd3bd).
# Backgrounds: manifest C4_designated (K3 idx 0, K4 idx 0), both C3 PASS.  Method as registered C4: periodic box 1024, pseudo-spectral D, RK4, T=1000, records every 5.
import numpy as np, json, time, argparse
M = 1.0; Lbox = 1024.0; s = 15.0; k0p = 0.5; EPS_REL = 1e-3; T = 1000.0; REC = 5.0

def run(th, coupled, N, dt, T=T):
    a_, b_, e_, w, q, v, r = th
    if not coupled: a_ = b_ = e_ = 0.0
    x = np.linspace(-Lbox/2, Lbox/2, N, endpoint=False); dx = x[1]-x[0]; kk = 2*np.pi*np.fft.fftfreq(N, dx)
    D = lambda f: np.real(np.fft.ifft(1j*kk*np.fft.fft(f)))
    gau = lambda A: A*np.exp(-x**2/(2*s*s))*np.cos(k0p*x)
    C = np.array([[a_, e_/2], [e_/2, b_]]); wv = np.array([w, v]); qv = np.array([q, r])
    Obg = wv@C@wv-qv@C@qv; chib = Obg/M**2                      # exact background: M^2 chi = O
    eps = EPS_REL*max(abs(q), abs(r))
    A = gau(eps); Bf = gau(eps)
    y = np.array([A, -D(A), Bf, D(Bf), np.zeros(N), np.zeros(N)])   # a (right-moving), b (left-moving), xi, ct
    def pieces(y):
        a, pa, b, pb, xi, ct = y
        Pt = np.array([w+pa, v+pb]); Px = np.array([q+D(a), r+D(b)]); chi = chib+xi
        Z = np.empty((2, 2, N)); 
        for i in range(2):
            for j in range(2): Z[i, j] = (1.0 if i == j else 0.0)+2*chi*C[i, j]
        return Pt, Px, chi, Z, xi, ct
    def rhs(y):
        a, pa, b, pb, xi, ct = y; Pt, Px, chi, Z, xi, ct = pieces(y); sx = D(xi)
        rh = np.empty((2, N))
        for i in range(2): rh[i] = 2*sum(C[i, j]*(ct*Pt[j]-sx*Px[j]) for j in range(2))
        det = Z[0, 0]*Z[1, 1]-Z[0, 1]*Z[1, 0]
        s0 = (Z[1, 1]*rh[0]-Z[0, 1]*rh[1])/det; s1 = (-Z[1, 0]*rh[0]+Z[0, 0]*rh[1])/det       # Z^-1 rh
        Oq = sum(C[i, j]*(Pt[i]*Pt[j]-Px[i]*Px[j]) for i in range(2) for j in range(2))
        return np.array([pa, D(D(a))-s0, pb, D(D(b))-s1, ct, D(D(xi))-M*M*chi+Oq])
    def E(y):
        Pt, Px, chi, Z, xi, ct = pieces(y)
        kin = sum(Z[i, j]*(Pt[i]*Pt[j]+Px[i]*Px[j]) for i in range(2) for j in range(2))
        return dx*(0.5*kin+0.5*ct**2+0.5*D(xi)**2+0.5*M*M*chi**2).sum(), float(min(np.linalg.eigvalsh(np.moveaxis(Z, 2, 0))[:, 0].min(), 1e9))
    Zbg = np.eye(2)+2*chib*C; Ebg = Lbox*(0.5*(wv@Zbg@wv+qv@Zbg@qv)+0.5*M*M*chib**2)
    n = int(round(T/dt)); every = int(round(REC/dt)); ts = [0.]; e0, z0 = E(y); Es = [e0]; Zm = [z0]; nr = [np.sqrt(dx*(y[0]**2+y[2]**2+y[4]**2).sum())]
    for i in range(1, n+1):
        k1 = rhs(y); k2 = rhs(y+.5*dt*k1); k3 = rhs(y+.5*dt*k2); k4 = rhs(y+dt*k3); y = y+dt/6*(k1+2*k2+2*k3+k4)
        if i % every == 0:
            ee, zz = E(y); ts.append(i*dt); Es.append(ee); Zm.append(zz); nr.append(np.sqrt(dx*(y[0]**2+y[2]**2+y[4]**2).sum()))
    ts, Es, nr, Zm = map(np.array, (ts, Es, nr, Zm)); Efl = Es-Ebg
    return dict(coupled=bool(coupled), N=N, dt=dt, eps=float(eps), n_records=len(ts), t_end=float(ts[-1]), E0=float(Es[0]), Ebg=float(Ebg), Efl0=float(Efl[0]),
                drift_total_max=float(np.abs(Es-Es[0]).max()/Es[0]), drift_total_final=float((Es[-1]-Es[0])/Es[0]), drift_fl_max=float(np.abs(Efl-Efl[0]).max()/abs(Efl[0])),
                norm0=float(nr[0]), norm_ratio_max=float(nr.max()/nr[0]), norm_ratio_final=float(nr[-1]/nr[0]), zmin_min=float(Zm.min())), ts, nr, Es

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--bg'); ap.add_argument('--coupled', type=int); ap.add_argument('--N', type=int); ap.add_argument('--dt', type=float)
    ap.add_argument('--out'); ap.add_argument('--T', type=float, default=T); ap.add_argument('--manifest', default='t3fund_hs2_manifest.json'); a = ap.parse_args()
    d = json.load(open(a.manifest))['C4_designated'][a.bg]; th = [d[k] for k in 'abewqvr']; t0 = time.time()
    res, ts, nr, Es = run(th, a.coupled, a.N, a.dt, a.T); res['seconds'] = time.time()-t0; res['bg'] = a.bg; res['theta'] = th
    tag = 'c4hs2_%s_%s_N%d' % (a.bg, 'coupled' if a.coupled else 'control', a.N)
    json.dump(res, open('%s/%s.json' % (a.out, tag), 'w'), indent=1); np.savetxt('%s/%s.csv' % (a.out, tag), np.c_[ts, nr, Es], header='t norm_fluct energy', delimiter=',')
    print(tag, json.dumps(res), flush=True)
