"""
Dipole-oscillation period vs amplitude for a Raman SOC BEC, single lower band.
Model: Zhang et al., PRL 109, 115301 (2012), arXiv:1201.6018, Eq. (1)-(2):
    H0 = [[(k+kr)^2/2m - d/2, W/2], [W/2, (k-kr)^2/2m + d/2]]
    kdot = -m w^2 x ,  xdot = dE_-(k)/dk
Units: E_R = kr^2/2m, kappa = k/kr, eps = E/E_R, tau = w t.
    eps_-(kappa) = kappa^2 + 1 - sqrt((2 kappa - d/2)^2 + (W/2)^2)
    dkappa/dtau = -xi ,  dxi/dtau = eps'(kappa)/2      (H = xi^2/2 + eps(kappa)/2)
Free case eps = kappa^2 + 1 gives period 2*pi in tau (= trap period T0).
Two-level only: the perturbative m_F=+1 correction used by the authors is NOT included.
Question: does the two-level E_- alone give a period that grows with amplitude
at (W, d) = (4.8, 1.0) E_R, as in their Fig. 2(b)?
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

K = sp.symbols('kappa', real=True)


def band(W, d):
    e = K**2 + 1 - sp.sqrt((2*K - sp.Rational(1, 2)*d)**2 + (W/2)**2)
    f = [sp.lambdify(K, sp.diff(e, K, n), 'numpy') for n in range(5)]
    return f


def analyse(W, d, amps, label):
    f = band(sp.nsimplify(W), sp.nsimplify(d))
    r = minimize_scalar(f[0], bounds=(-3, 3), method='bounded',
                        options={'xatol': 1e-13})
    k0 = r.x
    # check for a second local minimum on a grid
    g = np.linspace(-3, 3, 60001)
    e = f[0](g)
    nmin = int(np.sum((e[1:-1] < e[:-2]) & (e[1:-1] < e[2:])))
    e2, e3, e4 = f[2](k0), f[3](k0), f[4](k0)
    w0 = np.sqrt(e2/2)
    a3, b4 = e3/12, e4/48              # U = eps/2 = w0^2 u^2/2 + a3 u^3 + b4 u^4
    Kq = 3*b4/(2*w0)                   # Landau-Lifshitz: w = w0 + (Kq + Kc) A^2
    Kc = -15*a3**2/(4*w0**3)
    print(f"\n=== {label}: Omega={W} E_R, delta={d} E_R ===")
    print(f"local minima on grid: {nmin};  kappa0 = {k0:.6f}")
    print(f"eps''={e2:.6f}  eps'''={e3:.6f}  eps''''={e4:.6f}")
    print(f"small-amplitude T/T0 = 1/w0 = {1/w0:.5f}")
    print(f"LL shift dw/A^2: quartic {Kq:+.5f}, cubic {Kc:+.5f}, total {Kq+Kc:+.5f}"
          f"  -> period {'grows' if Kq+Kc < 0 else 'shrinks'} with A (2nd order)")
    rhs = lambda t, y: [-y[1], 0.5*f[1](y[0])]
    ev = lambda t, y: y[1]
    ev.direction = 1
    rows = []
    for side in (+1, -1):
        for a in amps:
            sol = solve_ivp(rhs, (0, 40*np.pi/w0), [k0 + side*a, 0.0],
                            events=ev, rtol=1e-11, atol=1e-13, max_step=0.05)
            te = sol.t_events[0]
            T = np.mean(np.diff(te[:6]))
            m = (sol.t >= te[0]) & (sol.t <= te[1])
            A = 0.5*(sol.y[0][m].max() - sol.y[0][m].min())
            rows.append((side, a, A, T/(2*np.pi)))
            print(f"start side {side:+d}  init offset {a:.2f}  A={A:.4f} k_r"
                  f"  T/T0={T/(2*np.pi):.5f}")
    return rows


if __name__ == "__main__":
    amps = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    main = analyse(4.8, 1.0, amps, "Fig.2(b) point")
    analyse(4.8, 0.0, [0.1, 0.3, 0.6, 0.9], "same Omega, delta=0 (no cubic)")
    analyse(4.0, 0.0, [0.1, 0.3, 0.6, 0.9], "transition point, delta=0")
    print("\n=== shape test vs Fig.2(b) (read-off: ~38 ms at A~0.3, ~42 ms at A~0.9) ===")
    for side in (+1, -1):
        pts = sorted((A, t) for s, a, A, t in main if s == side)
        As, ts = np.array(pts).T
        t3, t9 = np.interp(0.3, As, ts), np.interp(0.9, As, ts)
        print(f"side {side:+d}: T(0.9)/T(0.3) = {t9/t3:.4f}  (figure ~1.10);"
              f"  T(0.3)/T0 = {t3:.4f}, T(0.9)/T0 = {t9:.4f}")
        for T0 in (19.5, 22.2):
            print(f"    with T0={T0} ms: T(0.3)={t3*T0:.1f} ms, T(0.9)={t9*T0:.1f} ms")
