# Scalar mediator (Lambda = 0), evaluation against T1-T6

Construction: L_sigma = 1/2 sigma_t^2 - 1/2 sigma_x^2 - 1/2 M^2 sigma^2 + sigma (g_phi phi_t phi_x + g_psi psi_t psi_x); no bare Lambda.
Script: scalar_mediator_T1_T6.py. Verdicts hold for this construction and the sampled backgrounds only.

| Test | B3 anchor + 41 manifest backgrounds (J0 != 0) | J0 = 0 population (400 backgrounds, sqrtX<=0.3) |
|---|---|---|
| T1 | yes: same local polynomial F as Proca (det K / Q_sigma = -M^6 X^4) | yes |
| T2 | not scanned | 399/400 real on the x grid |
| T3 | 0/41 pass; front min 1.008, median 1.728, max 6.044 | 400/400 pass: front = 1 |
| T4 | same couplings as Proca, depends on the background only | passes by construction of the sample (sqrtX<=0.3) |
| T5 | not run | 400/400 positive, min eigenvalue 8.2e-7 |
| T6 | the coupling to phi_t phi_x singles out a frame; origin external, not dynamical; J0 depends on the frame | same |

Further results on the J0 = 0 population:
- Local (x->0) maximal phase speed is 1 (two luminal modes, beta = 0 to 1e-29); the other two modes are subluminal.
- Non-luminal modes: beta nonzero, median |beta| M^2 = 2.6e-3; Lambda_eff = 2 beta/c0 > 0 for all 800 roots (frame dependent; not identified with the bare Lambda).
- Moving off J0 = 0 (r -> r(1+d)): front median 1.0001 (d=0.01), 1.0011 (d=0.1); maximum 1.0056 at d=0.1.
Not established: T7 (operator matching), T8 (observable), behaviour outside the sampled range.

## ScalarMediator_DynamicU (U1-U4, first pass)
L = 1/2(d phi)^2 + 1/2(d psi)^2 + 1/2(d sigma)^2 - 1/2 M^2 sigma^2 + 1/2 f (d eta)^2 + sigma * sum_i g_i (u.d phi_i)(n.d phi_i),
u=(cosh eta, sinh eta), n=(sinh eta, cosh eta). Script: scalar_mediator_dynamicU_check.py.
- Equilibrium: sigma_bar = J(eta_bar)/M^2, J(eta)=A cosh 2eta + B sinh 2eta; stationarity J'(eta_bar)=0 gives tanh 2eta_bar = -B/A (needs |B|<|A|).
  Stationary eta_bar exists on 17 of 41 manifest backgrounds.
- sigma_bar = 0 together with J' = 0 requires A^2 = B^2 (degenerate, lightlike eta_bar): no generic solution.
- The principal symbol of the phi, psi blocks in the u frame equals the scalar-mediator one with a_i = g_i sigma_bar; f and the eta-sector enter only at lower order.
- Numerical front at X = k/M = 200 (degree-8 fit of det K, approximate): 0/17 with front <= 1; min 1.007, median 3.07, max 5.85; identical for f = 0.1, 1, 10.
- All 17 spectra real (no Im omega ~ k).
Verdict (this construction, sampled backgrounds): a dynamical rapidity eta does not repair T3.

## HeavyScalar_LorentzInvariant (one light field, first pass)
L = 1/2 (d phi)^2 + 1/2 (d chi)^2 - 1/2 M^2 chi^2 + g chi (d phi)^2, background phi0 = w t + q x, chi_bar = g (w^2-q^2)/M^2.
Script: heavy_scalar_lorentz_invariant.py (sympy det K, roots c = omega/k as functions of X = k/M).
- Principal symbol of both fields is proportional to (omega^2 - k^2): front = 1 for any background and any g (Lorentz-invariant principal part).
- Roots real at all scanned X (1e-3 .. 1e3) for (w,q) = (1,0.3), (0.3,1) (timelike and spacelike gradient), g = +-0.4. Spectrum even in g.
- Light mode: |c| < 1 at low X (0.706 / -0.905 for (1,0.3); 0.14 / -0.641 for (0.3,1)), increasing monotonically toward 1 from below. Heavy mode: |c| > 1 at low X (c ~ M/k, massive branch), -> 1 from above at high X.
- Sign: the effective operator is +g^2/(2M^2) X^2 > 0 for any real g. This is the subluminal sign (cf. positivity bounds, recalled from literature, not re-derived here) and the same sign as the GP quartic term; B4 with Lambda > 0 has local speed sqrt(1+3 Lambda q^2) > 1.
Scope: one light field, tree level, this coupling; not yet two fields, kappa cross term or T4/T5.

## HeavyScalar_TwoField (T1/T7 operator-level test, first pass)
L_int = chi O, O = a (d phi)^2 + b (d psi)^2 + e d phi.d psi; low-energy L_eff = 1/2 X_phi + 1/2 X_psi + O^2/(2 M^2). Script: heavy_scalar_two_field_matching.py.
- Light-cone weights: with a = d_+ phi, b = d_- phi, phi_x^4 and (phi_t phi_x)^2 contain a^4, b^4 (boost weight +-4); every Lorentz-invariant quartic has weight 0
  (X^2 = 16 a^2 b^2). So the B4 operators Lambda phi_x^4 and p (phi_t phi_x)^2 are not Lorentz invariant as operators.
- Numerical: 15 random backgrounds (|w,q,v,r| <= 0.6), fit constant (p, s, kappa, Lambda) of the B4 local polynomial to the Lorentz-invariant one
  (60 random starts). Control (B4 -> B4): residual 1e-17. Lorentz-invariant couplings (a,b,e) = (0.2,0.2,0), (0.3,-0.2,0.1), (0.2,0.1,0.4), (0.5,0.5,0):
  residual / deviation-from-free = 0.18, 0.24, 0.60, 0.47; fitted kappa^2/(ps) = 9.7, 13.0, 9.7, 29 (not rank one).
- Not done: principal-symbol positivity (ghost freedom) of the full phi, psi, chi system, T4/T5, couplings beyond the four samples.
Scope: this class (one heavy scalar, tree level, Lorentz-invariant coupling) and the sampled backgrounds/couplings.

## B4 Lorentz classification (analytic + manifest check)
Q1 (covariant form with a frame field): yes. With u = (cosh eta, sinh eta), n = eps u, the B4 terms become Lambda (n.d phi)^4, p (u.d phi)^2 (n.d phi)^2 etc.;
the free part is Lorentz invariant. With a frame defined by the background ( u = d psi / sqrt(X_psi) ) the operator is (n.d phi)^2 = G / X_psi,
G = X_phi X_psi - Y^2, i.e. Lambda/4 * G^2 / X_psi^2: Lorentz invariant, non-analytic at X_psi -> 0, defined only for timelike d psi.
Q2 (frame from the background itself, no external structure): the covariant model equals lab-frame B4 only when the chosen u equals the lab frame.
Manifest check (anchor + 40, n=41): X_phi>0 on 27, X_psi>0 on 24; |q/w| in 0.085..25, |r/v| in 0.12..6.0;
the constant lambda with q + lambda r = 0 (u from d phi + lambda d psi) ranges -8.6..4.0 (median -0.45): no fixed lambda works.
Q3: constant p, s, kappa, Lambda across backgrounds with a background-defined frame is not obtained for the manifest backgrounds.

## TwoField_DynamicInteraction (Hypothesis B, first pass)
Model: Lorentz-invariant heavy scalar, L = 1/2 X_phi + 1/2 X_psi + 1/2 (d chi)^2 - 1/2 M^2 chi^2 + chi (a X_phi + b X_psi + e Y), (a,b,e) = (0.2,0.2,0.4), M = 1.
Setup: counter-propagating Gaussian packets (k0 = 0.3, width 15, starting at -+100), periodic box 1024, RK4 + FFT derivatives, T = 200, compared with the g = 0 control.
Script: two_field_dynamic_interaction.py.
- Equal amplitudes (A = 1): centroid shift -1.42 for phi and +1.42 for psi (time delay, no advance); change of central wavenumber 1.2e-6; bandwidth 2.7e-6;
  relative energy change of phi -1.0e-5; residual chi energy 2.7e-5. Doubling the amplitude of psi: phi shift -5.51 (ratio 3.9, expected ~4 for A^2 scaling), but
  larger spectral and energy changes (dk_c -1.0e-3, dE_phi -6.2%, E_chi 9.4e-4) and unequal shifts (-5.51 vs +4.32): the nonlinear/EFT regime is left.
- Convergence: (N, dt) = (8192, 0.05) and (16384, 0.025) give the same shift -1.42063 and the same E_phi.
- Conclusion for this model: the interaction produces a scattering time delay ~ A_other^2 and small heavy-field radiation; in vacuum the asymptotic fields are free (omega = |k|),
  so no persistent change of dispersion is seen. A Lambda_eff-like effect exists only in the overlap region, as a background-induced effect (Hypothesis A localised).
Not done: probe-packet dispersion measurement after the collision, other couplings and k0/M, energy-conservation diagnostics, B4 lab-frame models.

## Probe on a stationary background (Hypothesis A in dynamics, first pass)
Model/couplings as above, (a,b,e)=(0.2,0.2,0.4), M=1; background (w,q,v,r)=(1,0.3,0.8,-0.2), chi_bar = O/M^2 = 0.646 (exact constant solution; Z = [[1+2a chi, e chi],[e chi, 1+2b chi]] has eigenvalues 1 and 1.517, positive definite).
Script: probe_on_background.py. Analytic roots c = omega/k of det K (3 fields) vs time-domain evolution of a weak probe (full nonlinear equations in fluctuation variables, T=1600, dt=0.1, N=2048, box 1024).
- Analytic: light modes c = 0.7484 / -0.7933 (k=0.1), 0.7642 / -0.8047 (k=0.5), |c| increasing with k; one exactly luminal pair (+-1) at all k; heavy branch |c| >> 1 at small k, decreasing toward 1.
- Time domain (probe amplitude 1e-6 and 1e-4, identical spectra): peaks agree with the analytic light roots to about 1e-3 (k=0.3: 0.7543 vs 0.7541, 0.7967 vs 0.7974; k=0.5: 0.7862?-0.7644 and 0.8049 vs 0.8045; limited by frequency resolution 2 pi/T); sidelobes of the Hann window are also listed.
- beta from two analytic points: right mover ~ +0.066 (c0 ~ 0.748), left mover ~ -0.048 (c0 ~ -0.793); Lambda_eff = 2 beta / c0 ~ +0.18 and +0.12 (frame dependent, not the bare B4 Lambda).
Not done: explicit growth-rate / energy-conservation check; amplitudes where the probe changes the background; comparison to the B4 operator structure (T1/T7 mismatch from HeavyScalar_TwoField stays).

## Stability / energy check of the heavy scalar probe (same background)
Fluctuation amplitude 1e-3, T = 1600: total Hamiltonian drift (E - E0)/E0 = -5.2e-13 (dominated by the background energy 1549.88), fluctuation norm bounded between 1.0e-2 and 1.2e-2 (no growth).
Principal symbol of (phi, psi, chi) = diag(Z, 1)(omega^2 - k^2), Z eigenvalues 1 and 1.517 > 0 on this background: ghost free there. Script: experiments/uv_candidates/stability_check.py (same equations as probe_on_background.py).
