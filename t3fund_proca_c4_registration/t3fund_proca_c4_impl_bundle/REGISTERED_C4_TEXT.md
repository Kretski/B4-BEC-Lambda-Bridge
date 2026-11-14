# Registered C4 text (verbatim extracts; nothing here is new)
Source 1: `T3_FUND_PROCA_SPEC.md` as committed in e355af2; recomputed SPEC_HASH of the bytes before "## 9. Hashes" = 86c360326da7406f (matches the registered value).

## 6. C4 (time domain)
Designated backgrounds (rule fixed before any evaluation, independent of outcomes): (a) the first two PN draws in manifest order (cell order cut 0.01 → 0.95, rank1_minus before rank1_plus) that are C3 PASS and have UV front speed max|c| ≤ 2.0; (b) the first PE case in manifest order that is C3 PASS with max|c| ≤ 2.0. The bound max|c|≤2.0 is a numerical-stability rule for the registered resolutions: with pseudo-spectral k_max=π/dx=50.3 at (16384, dt=0.025) the RK4 condition ω_max·dt<2.83 gives max|c|<2.25. A designated background that is not C3 PASS is NOT-RUN with no replacement. If a rule selects nothing, C4 is NOT-AVAILABLE for that rule (reported; e.g. the B3 anchor, whose front exceeds 2, is outside the C4 domain).
Method as registered C4: periodic box 1024, pseudo-spectral D, RK4, T=1000, records every 5, base resolution (16384, dt 0.025), convergence (32768, 0.0125). State (a,b,c̃,d,A1,F) evolved with the reduced system of §3(iii) (Fourier multipliers for F_xx); the implementation is checked against the independent sympy Euler–Lagrange operator at random states (≤1e-12) before any run. Energy E of §3(iii) (full nonlinear, conserved), E_fl = E − E_bg.
Perturbation: δf right-moving (δf_t=−Dδf), δg left-moving (δg_t=+Dδg), δA1 Gaussian with δF=0; Gaussian s=15, k0p=0.5, centre 0, amplitude ε=1e-3·max(|q|,|r|). Norm 𝒩=√(dx Σ(δf²+δg²+δA1²)); criterion max 𝒩/𝒩(0) ≤ 3. Free control: g_φ=g_ψ=0 (Ā=0), same gradients, Λ, perturbation, resolutions. Outcome rule as registered C4 incl. Clarification 1 (PASS: ratio ≤3 at both resolutions, total energy drift ≤1e-9 reached, control reaches target, 1% agreement; INCONCLUSIVE: control misses target / 1% agreement fails / drift missed with control passing / non-finite values; FAIL-UNSTABLE: ratio >3 with agreement; never FAIL on energy alone). Stated limitation: total-energy drift is dominated by the background; the norm criterion is substantive. Passing C4 on small-X backgrounds is weak evidence for the strong-X ones.

Source 2: `UV_ACCEPTANCE_PROTOCOL_ADDENDUM_T3_T6.md` (registered protocol; P7 and the outcome section), verbatim:

```
P7 C4 (clarification D).
   - Initial data: Gaussian wave packet(s) in the light fields, amplitude 1e-3 of the background gradient scale, heavy field unperturbed; fully specified in the script.
   - Norm: L2 norm over the box of all fluctuation fields (light and heavy), recorded at >= 100 times over a duration >= 1000/M; the time series is plotted.
     Test criterion: max over time of norm / initial norm <= 3. This is a criterion for this test, not a definition of stability.
   - Energy: the total conserved energy of the full system (all fields). For Lagrangian models without dissipation its conservation is a requirement on
     the integrator, not on the model. For a model with dissipation or sources the stated balance law replaces conservation.
     Target: relative drift of the total energy <= 1e-9 over the run. The target is not assumed to be attainable: it is calibrated on the free system
     (all couplings set to zero) with the same numerical method, time step and spatial grid, repeated at (2N, dt/2); the drift must decrease under refinement.
     If the free system does not reach the target, C4 is INCONCLUSIVE for the candidates, not a FAIL of the physical model.
     Calibration run (with a perturbation; a free system without perturbation conserves energy trivially and does not test the integrator):
       (a) the background gradients of the light fields are kept, to the extent that they also solve the free equations;
       (b) the background values of heavy and auxiliary fields are chosen so that they satisfy the FREE equations (they are not copied from the coupled
           background: with the coupling switched off, the coupled value of a sourced heavy field is in general no longer a solution);
       (c) the same profile, shape and relative amplitude of the initial perturbation as in the candidate run;
       (d) the same numerical method, time step and spatial grid, and then (2N, dt/2).
     Additional diagnostic, reported but not a pass criterion: drift of the fluctuation energy E_fl = E_total - E_background with the background part
     subtracted analytically. (The total energy is dominated by the background, so the total-energy drift alone is a weak test of the fluctuations.)
   - Convergence: all of the above repeated at (N, dt) and (2N, dt/2); results agree within 1 %.


### Outcomes per background
- PASS: C0-C3 hold with the required evidence level (C4 on designated backgrounds).
- FAIL-ILLPOSED: C1 or C2 fails.
- FAIL-UNSTABLE: C1 and C2 hold, and a growing mode is confirmed in a finite band. Counts as a T3-fund failure.
- INCONCLUSIVE: C0 not satisfied; unresolved roots; no analytic argument for C1 (FUNDAMENTAL) or for a degenerate C2; C4 not converged.
Per candidate: the counts of the four outcomes on each background set, together with the declared regime. "Passes T3-fund on a set" is used only when
there is no FAIL and no INCONCLUSIVE on that set.
```

Clarification 1 of that addendum (verbatim first lines): `Clarification 1 (project lead, 2026-11-13, recorded before any evaluation): the three-band rule for gamma against u(k) in P6 and the explicit ZERO MODE`. It concerns P6 (three-band rule) and P2a (zero mode), not C4; the parenthesis "incl. Clarification 1" in Proca spec section 6 is followed by the explicit rule list, which is what the aggregator implements.

## Two points of the review, resolved by the registered text
A. Initial excitation: Proca spec section 6 states explicitly "delta A1 Gaussian with delta F=0" (and delta f, delta g). The general protocol P7 says "heavy field unperturbed" ("Gaussian wave packet(s) in the light fields ... heavy field unperturbed"). The two texts differ; the Proca spec (registered later, model-specific, written before any evaluation) is what the code implements. If the reviewer rules that the general P7 governs, delta A1 must be set to zero in `initial_state` (one line) and the registered Proca text would need an explicit clarification before runs; no change has been made.
B. Control failure: the Proca spec lists "control misses target" among the INCONCLUSIVE conditions, and P7 says: "If the free system does not reach the target, C4 is INCONCLUSIVE ... not a FAIL of the physical model." The aggregator (v1) returned FAIL-UNSTABLE before checking the control; this was wrong against the registered text and is corrected in v2: control misses target -> INCONCLUSIVE with priority over FAIL-UNSTABLE.
