# T3-fund C4 (time-domain confirmation) for the designated backgrounds; registered protocol: UV_ACCEPTANCE_PROTOCOL_ADDENDUM_T3_T6.md, P7; calibration: P7_CALIBRATION_RESULT.md.
# SPECIFICATION (frozen before the first C4 run; SPEC_HASH printed):
#  Designated backgrounds (project lead, 2026-11-13, before any C4 run):
#     SM : Scalar mediator, first background of set A1 (seed 11): g_phi=0.15789059608093756, g_psi=-0.1705408351625817, w=1.3449853598753254, q=0.3656507783891485,
#          v=-0.39302062881062705, r=-1.1585014776652514, M=1.
#     HS : HeavyScalar, one light field, B2 anchor: w=1, q=0.3, g=+0.4, M=1.
#  Method (same as the P7 calibration): periodic box L=1024, pseudo-spectral D = real(ifft(i k fft)), RK4, D(D(f)) for second derivatives, T=1000, records every 5 (201).
#  Resolutions: base (N,dt)=(16384,0.025); convergence (32768,0.0125).
#  The exact background (constant gradients, constant mediator value) is carried analytically; the state holds the periodic fluctuations (a, pa=a_t, [b, pb], xi, ct):
#     Phi_t = w + pa, Phi_x = q + D(a), (Psi_t = v + pb, Psi_x = r + D(b)), mediator = bar + xi.  Full nonlinear equations are integrated:
#     HS : pa_t = D(D(a)) - (2g/Z)(ct Phi_t - D(xi) Phi_x),  ct_t = D(D(xi)) - M^2 (chi_bar + xi) + g (Phi_t^2 - Phi_x^2),  Z = 1 + 2 g (chi_bar + xi)
#     SM : pa_t = D(D(a)) - g_phi (ct Phi_x + D(xi) Phi_t) - 2 g_phi sigma D(pa),  (same for b with g_psi, Psi),  ct_t = D(D(xi)) - M^2 sigma + g_phi Phi_t Phi_x + g_psi Psi_t Psi_x,  sigma = sigma_bar + xi
#     Energy: HS: sum 1/2 Z (Phi_t^2 + Phi_x^2) + 1/2 ct^2 + 1/2 D(xi)^2 + 1/2 M^2 chi^2;  SM: sum 1/2 (Phi_t^2 + Phi_x^2 + Psi_t^2 + Psi_x^2) + 1/2 ct^2 + 1/2 D(xi)^2 + 1/2 M^2 sigma^2 (times dx).
#  Initial perturbation (same shape as P7): Gaussian packets s=15, k0p=0.5, centre x=0, amplitude eps = 1e-3 * max |d_x phi_bg| over the light fields (HS: |q|; SM: max(|q|,|r|)),
#     in the light fields only (a: right-moving pa=-D(a); b: left-moving pb=+D(b)); mediator unperturbed (xi=0, ct=0).
#  Norm: N(t) = sqrt(dx sum(a^2 (+b^2) + xi^2)) of the fluctuation fields. Criterion: max_t N/N(0) <= 3. Energy: max_t |E(t)-E(0)|/E(0) <= 1e-9 required (target); diagnostic E_fl = E - E_bg(analytic).
#  Free control (P7 (a)-(d)): same code with all couplings set to zero, mediator background 0 (solves the free equations), same gradients, same perturbation, same (N,dt), then (2N,dt/2).
#  Outcome for the candidate on the designated background: PASS if at both resolutions N-ratio <= 3 and the energy-drift target is reached AND the free control reaches the target at the same resolution and
#     the max N-ratio agrees between resolutions within 1 %;  INCONCLUSIVE if the free control does not reach the target or the 1 % agreement fails;  FAIL-UNSTABLE if N-ratio > 3 at either resolution
#     with agreement between resolutions (an actual growth); the physical model is never failed on the energy criterion alone.
#  Clarification (project lead, 2026-11-13, before the first C4 run): if the free control meets the energy criterion but the coupled system does not reach the drift target for the total energy, the outcome
#     is INCONCLUSIVE. FAIL of the physical model is never based on the energy drift alone. A separate diagnostic of the integrator, of the energy functional and of the convergence is then required.
#  Order of runs: (1) Scalar mediator, coupled and free control, both resolutions; (2) HeavyScalar, coupled and free control, both resolutions; (3) energy drift, norm, convergence; (4) report and commit.
#  After registration the amplitude 1e-3, T = 1000, the box size and the resolutions are not changed.
import numpy as np, json, hashlib, sys, time, argparse

M = 1.0; Lbox = 1024.0; s = 15.0; k0p = 0.5; EPS_REL = 1e-3; T = 1000.0; REC = 5.0
BG = dict(
    SM=dict(gp=0.15789059608093756, gs=-0.1705408351625817, w=1.3449853598753254, q=0.3656507783891485, v=-0.39302062881062705, r=-1.1585014776652514),
    HS=dict(g=0.4, w=1.0, q=0.3))

def run(model, coupled, N, dt, T=T):
    x = np.linspace(-Lbox/2, Lbox/2, N, endpoint=False); dx = x[1]-x[0]; kk = 2*np.pi*np.fft.fftfreq(N, dx)
    D = lambda f: np.real(np.fft.ifft(1j*kk*np.fft.fft(f)))
    gau = lambda A: A*np.exp(-(x-0.)**2/(2*s*s))*np.cos(k0p*x)
    P = BG[model]
    if model == 'HS':
        g = P['g'] if coupled else 0.0; w, q = P['w'], P['q']; chib = g*(w*w-q*q)/M**2
        eps = EPS_REL*abs(q); a = gau(eps); pa = -D(a); xi = np.zeros(N); ct = np.zeros(N)
        y = np.array([a, pa, xi, ct])
        def rhs(y):
            a, pa, xi, ct = y; Pt = w+pa; Px = q+D(a); chi = chib+xi; Z = 1+2*g*chi
            return np.array([pa, D(D(a))-(2*g/Z)*(ct*Pt-D(xi)*Px), ct, D(D(xi))-M*M*chi+g*(Pt**2-Px**2)])
        def E(y):
            a, pa, xi, ct = y; Pt = w+pa; Px = q+D(a); chi = chib+xi; Z = 1+2*g*chi
            return dx*(0.5*Z*(Pt**2+Px**2)+0.5*ct**2+0.5*D(xi)**2+0.5*M*M*chi**2).sum()
        Ebg = Lbox*(0.5*(1+2*g*chib)*(w*w+q*q)+0.5*M*M*chib**2)
        fl = lambda y: y[[0, 2]]
    else:
        gp = P['gp'] if coupled else 0.0; gs = P['gs'] if coupled else 0.0; w, q, v, r = P['w'], P['q'], P['v'], P['r']
        sb = (gp*w*q+gs*v*r)/M**2
        eps = EPS_REL*max(abs(q), abs(r)); a = gau(eps); pa = -D(a); b = gau(eps); pb = D(b); xi = np.zeros(N); ct = np.zeros(N)
        y = np.array([a, pa, b, pb, xi, ct])
        def rhs(y):
            a, pa, b, pb, xi, ct = y; Pt = w+pa; Px = q+D(a); Qt = v+pb; Qx = r+D(b); sg = sb+xi; sx = D(xi)
            return np.array([pa, D(D(a))-gp*(ct*Px+sx*Pt)-2*gp*sg*D(pa), pb, D(D(b))-gs*(ct*Qx+sx*Qt)-2*gs*sg*D(pb), ct,
                             D(D(xi))-M*M*sg+gp*Pt*Px+gs*Qt*Qx])
        def E(y):
            a, pa, b, pb, xi, ct = y; Pt = w+pa; Px = q+D(a); Qt = v+pb; Qx = r+D(b); sg = sb+xi
            return dx*(0.5*(Pt**2+Px**2+Qt**2+Qx**2)+0.5*ct**2+0.5*D(xi)**2+0.5*M*M*sg**2).sum()
        Ebg = Lbox*(0.5*(w*w+q*q+v*v+r*r)+0.5*M*M*sb**2)
        fl = lambda y: y[[0, 2, 4]]
    n = int(round(T/dt)); every = int(round(REC/dt)); ts, Es, Ns = [0.], [E(y)], [np.sqrt(dx*(fl(y)**2).sum())]
    for i in range(1, n+1):
        k1 = rhs(y); k2 = rhs(y+.5*dt*k1); k3 = rhs(y+.5*dt*k2); k4 = rhs(y+dt*k3); y = y+dt/6*(k1+2*k2+2*k3+k4)
        if i % every == 0:
            ts.append(i*dt); Es.append(E(y)); Ns.append(np.sqrt(dx*(fl(y)**2).sum()))
            if not np.isfinite(Es[-1]): break
    ts, Es, Ns = map(np.array, (ts, Es, Ns))
    Efl = Es-Ebg
    return dict(model=model, coupled=bool(coupled), N=N, dt=dt, eps=float(eps), n_records=len(ts), t_end=float(ts[-1]), E0=float(Es[0]), Ebg=float(Ebg), Efl0=float(Efl[0]),
                drift_total_max=float(np.abs(Es-Es[0]).max()/Es[0]), drift_total_final=float((Es[-1]-Es[0])/Es[0]), drift_fl_max=float(np.abs(Efl-Efl[0]).max()/abs(Efl[0])),
                norm0=float(Ns[0]), norm_ratio_max=float(Ns.max()/Ns[0]), norm_ratio_final=float(Ns[-1]/Ns[0])), ts, Ns, Es

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument('--model', required=True, choices=['SM', 'HS']); ap.add_argument('--coupled', type=int, required=True)
    ap.add_argument('--N', type=int, required=True); ap.add_argument('--dt', type=float, required=True); ap.add_argument('--out', default='.'); ap.add_argument('--T', type=float, default=T)
    a = ap.parse_args()
    spec = open(__file__).read().split("import numpy")[0]; print("SPEC_HASH", hashlib.sha256(spec.encode()).hexdigest()[:16], flush=True)
    t0 = time.time(); r, ts, Ns, Es = run(a.model, a.coupled, a.N, a.dt, a.T); r['seconds'] = time.time()-t0
    tag = "c4_%s_%s_N%d" % (a.model, 'coupled' if a.coupled else 'control', a.N)
    np.savetxt("%s/%s.csv" % (a.out, tag), np.c_[ts, Ns, Es], header="t norm_fluct E_total", delimiter=",")
    json.dump(r, open("%s/%s.json" % (a.out, tag), "w"), indent=1); print(json.dumps(r), flush=True)
