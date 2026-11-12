"""T6.3: K(z) for the cosmological (comoving) foliation, and a numerical check of the first-order formula.

Local dispersion measured by the comoving observer (u_mu = -d_mu T, T = cosmic time):  omega = c k sqrt(1 + Lam k^2),  k = kappa/a (physical).
Group velocity: v_g = c (1 + 2 Lam k^2)/sqrt(1 + Lam k^2) = c (1 + 3/2 Lam k^2 + O(Lam^2 k^4)).
Null-like ray: d chi = v_g dt / a.  Same comoving distance for all frequencies  =>  arrival-time shift at observed frequency f (k0 = 2 pi f / c, a0 = 1):
    delta t(f) = -(3/2) Lam k0^2 * Int_0^z (1+z')^2 / H(z') dz' = -6 pi^2 Lam f^2 K(z) / c^3,     K(z) = c * Int_0^z (1+z')^2 / H(z') dz'   [length]
Stationary phase with h~(f) = A exp(i Psi), Psi = 2 pi f t_c - ... :  dPsi/df = 2 pi t(f)  =>  delta Psi(f) = -4 pi^3 Lam K(z) f^3 / c^3 .
(The overall sign flips if the opposite Fourier convention is used.)
"""
import numpy as np
from scipy.integrate import quad, solve_ivp

c = 299792.458            # km/s
H0, Om = 67.7, 0.31       # km/s/Mpc, flat LCDM (illustrative values)
Hz = lambda z: H0*np.sqrt(Om*(1+z)**3 + 1-Om)
K = lambda z: c*quad(lambda x: (1+x)**2/Hz(x), 0, z)[0]          # Mpc
DC = lambda z: c*quad(lambda x: 1/Hz(x), 0, z)[0]                # comoving distance, Mpc
print(" z     K(z) [Mpc]   D_C [Mpc]   D_L [Mpc]   K/D_L")
for z in (0.01, 0.1, 0.2, 0.5, 1.0, 2.0):
    print("%4.2f  %10.1f  %10.1f  %10.1f  %6.3f" % (z, K(z), DC(z), (1+z)*DC(z), K(z)/((1+z)*DC(z))))

# ---- numerical check of delta t: integrate the exact group velocity through the expansion (dimensionless: c = H0 = 1)
E = lambda z: np.sqrt(Om*(1+z)**3 + 1-Om)
def arrival_shift(eps, z_s):
    """eps = Lam k0^2 (dimensionless, observed). Returns (t_arrival(eps) - t_arrival(0)) * H0, emitting at the same time from the same comoving distance."""
    chi_s = quad(lambda x: 1/E(x), 0, z_s)[0]                       # comoving distance travelled by the eps = 0 ray
    # integrate d chi / dz = v_g/(c E) backwards from the source; find the redshift offset where chi = chi_s, convert to time
    vg = lambda z: (1+2*eps*(1+z)**2)/np.sqrt(1+eps*(1+z)**2)
    # chi(z_s -> z_end) with exact v_g; solve for z_end such that chi = chi_s
    from scipy.optimize import brentq
    f = lambda zend: quad(lambda x: vg(x)/E(x), zend, z_s)[0] - chi_s
    zend = brentq(f, -0.5, z_s)                                     # zend > 0 means it arrives before "today"
    # time between zend and z = 0:  dt = dz / ((1+z) E)
    return -quad(lambda x: 1/((1+x)*E(x)), 0, zend)[0]
for z_s in (0.1, 1.0):
    I = quad(lambda x: (1+x)**2/E(x), 0, z_s)[0]
    for eps in (1e-6, 1e-4):
        num = arrival_shift(eps, z_s); ana = -1.5*eps*I
        print("z=%.1f eps=%.0e  delta t numeric %.6e  formula %.6e  ratio %.5f" % (z_s, eps, num, ana, num/ana))

# ---- scale: Lambda that gives |delta Psi| = 1 rad at f = 100 Hz
cm = 299792458.0; Mpc = 3.0856775814913673e22
for z in (0.1, 1.0):
    Lam = cm**3/(4*np.pi**3*K(z)*Mpc*100.0**3)
    print("z=%.1f: |dPsi| = 1 rad at 100 Hz  <=>  Lambda = %.3e m^2  (sqrt = %.3e m)" % (z, Lam, np.sqrt(Lam)))
