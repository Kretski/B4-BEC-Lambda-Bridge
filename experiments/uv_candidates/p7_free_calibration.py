# P7 calibration on the FREE system (all couplings zero).  Registered protocol: UV_ACCEPTANCE_PROTOCOL_ADDENDUM_T3_T6.md, P7 (a)-(d).
# SPECIFICATION (frozen before the first run; see SPEC_HASH printed at start; nothing below is changed after results):
#   Method  : same as two_field_dynamic_interaction.py: periodic box, pseudo-spectral D = ifft(i k fft), RK4, D(D(f)) for second derivatives.
#   Fields  : phi, psi (massless, light), chi (heavy, M=1), free equations  f_tt = f_xx (light),  chi_tt = chi_xx - M^2 chi.
#   Grid    : Lbox=1024; (N,dt) = (8192,0.05) and (2N,dt/2) = (16384,0.025); T=1000 (=1000/M); records every 5 time units (201 times).
#   Background (solves the free equations): phi_bg = right-moving packet, A=1, s=15, k0=0.3, x=-100 ; psi_bg = left-moving packet, same, x=+100 ;
#            chi_bg = 0 and chi_t = 0 (solves the free heavy equation).
#   Perturbation (same shape in the candidate runs): Gaussian packet, s=15, k0p=0.5, centred x=0, right-moving in phi, left-moving in psi;
#            field amplitude eps = 1e-3 * max|d_x phi_bg|  (background gradient scale); heavy field unperturbed.
#   Fluctuation = (total run) - (background-only run).  Fluctuation norm = sqrt(dx * sum(phi^2+psi^2+chi^2)) of the fluctuation fields.
#   Energy  : E = dx*sum[ 1/2(phi_t^2+phi_x^2) + 1/2(psi_t^2+psi_x^2) + 1/2(chi_t^2+chi_x^2+M^2 chi^2) ] of the total system.
#   Drift   : max_t |E(t)-E(0)|/E(0) (and signed final). Diagnostic E_fl = E_total - E_bg(0), drift relative to E_fl(0) (energy of perturbation + cross term).
#   Extra (NOT in the protocol, reported separately and labelled): same run with a perturbed chi packet, to exercise the heavy-field integrator.
import numpy as np, json, hashlib, sys, time

M = 1.0; Lbox = 1024.0; s = 15.0; k0 = 0.3; k0p = 0.5; EPS_REL = 1e-3; T = 1000.0; REC = 5.0

def run(N, dt, with_pert, chi_pert=False):
    x = np.linspace(-Lbox/2, Lbox/2, N, endpoint=False); dx = x[1]-x[0]; kk = 2*np.pi*np.fft.fftfreq(N, dx)
    D = lambda f: np.real(np.fft.ifft(1j*kk*np.fft.fft(f)))
    g = lambda x0, A, k: A*np.exp(-(x-x0)**2/(2*s*s))*np.cos(k*(x-x0))
    phi = g(-100., 1., k0); psi = g(100., 1., k0)
    eps = EPS_REL*np.abs(D(phi)).max()           # gradient scale of the background, computed on the same grid
    chi = np.zeros(N); ct = np.zeros(N)
    pt = -D(phi); st = D(psi)                     # right-moving phi, left-moving psi
    if with_pert:
        pp = g(0., eps, k0p); ps = g(0., eps, k0p)
        phi = phi+pp; psi = psi+ps; pt = pt-D(pp); st = st+D(ps)
        if chi_pert: chi = g(0., eps, k0p)
    def E(y):
        phi, pt, psi, st, chi, ct = y
        return dx*(0.5*(pt**2+D(phi)**2)+0.5*(st**2+D(psi)**2)+0.5*(ct**2+D(chi)**2+M*M*chi**2)).sum()
    def rhs(y):
        phi, pt, psi, st, chi, ct = y
        return np.array([pt, D(D(phi)), st, D(D(psi)), ct, D(D(chi))-M*M*chi])
    y = np.array([phi, pt, psi, st, chi, ct]); n = int(round(T/dt)); every = int(round(REC/dt))
    ts, Es, Ys = [0.], [E(y)], [y[[0, 2, 4]].copy()]
    for i in range(1, n+1):
        k1 = rhs(y); k2 = rhs(y+.5*dt*k1); k3 = rhs(y+.5*dt*k2); k4 = rhs(y+dt*k3); y = y+dt/6*(k1+2*k2+2*k3+k4)
        if i % every == 0: ts.append(i*dt); Es.append(E(y)); Ys.append(y[[0, 2, 4]].copy())
    return np.array(ts), np.array(Es), np.array(Ys), dx, eps

def analyse(N, dt, chi_pert=False):
    t0 = time.time()
    ts, Et, Yt, dx, eps = run(N, dt, True, chi_pert)
    _, Eb, Yb, _, _ = run(N, dt, False)
    fl = Yt-Yb                                    # fluctuation fields (phi,psi,chi)
    nrm = np.sqrt(dx*(fl**2).sum(axis=(1, 2)))
    drift_t = np.abs(Et-Et[0])/Et[0]
    Efl = Et-Eb[0]                                # E_bg(0) analytic = energy of the background at t=0 (conserved)
    dflt = np.abs(Efl-Efl[0])/abs(Efl[0])
    return dict(N=N, dt=dt, eps=float(eps), n_records=len(ts), drift_total_max=float(drift_t.max()), drift_total_final=float((Et[-1]-Et[0])/Et[0]),
                drift_fl_max=float(dflt.max()), norm0=float(nrm[0]), norm_ratio_max=float(nrm.max()/nrm[0]), norm_final_ratio=float(nrm[-1]/nrm[0]),
                E0=float(Et[0]), Efl0=float(Efl[0]), seconds=time.time()-t0), ts, nrm, drift_t

if __name__ == "__main__":
    spec = open(__file__).read().split("import numpy")[0]
    print("SPEC_HASH", hashlib.sha256(spec.encode()).hexdigest()[:16], flush=True)
    out = {}
    for tag, cp in (("protocol", False), ("extra_chi_perturbed_NOT_PROTOCOL", True)):
        res = {}
        for (N, dt) in ((8192, 0.05), (16384, 0.025)):
            r, ts, nrm, dr = analyse(N, dt, cp); res["N%d_dt%g" % (N, dt)] = r
            np.savetxt("p7_%s_N%d.csv" % (tag, N), np.c_[ts, nrm, dr], header="t norm_fluct drift_total", delimiter=",")
            print(tag, json.dumps(r), flush=True)
        a, b = res["N8192_dt0.05"], res["N16384_dt0.025"]
        res["refinement"] = dict(drift_total_ratio=a["drift_total_max"]/b["drift_total_max"] if b["drift_total_max"] > 0 else None,
                                 norm_ratio_max_rel_diff=abs(a["norm_ratio_max"]-b["norm_ratio_max"])/a["norm_ratio_max"],
                                 drift_decreases=bool(b["drift_total_max"] < a["drift_total_max"]),
                                 target_reached_coarse=bool(a["drift_total_max"] <= 1e-9), target_reached_fine=bool(b["drift_total_max"] <= 1e-9))
        print(tag, "refinement", json.dumps(res["refinement"]), flush=True)
        out[tag] = res
    json.dump(out, open("p7_free_calibration_result.json", "w"), indent=1)
