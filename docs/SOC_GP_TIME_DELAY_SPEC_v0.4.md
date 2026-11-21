# SOC BEC: Gross–Pitaevskii test of R = 1 + 12Q — pre-registered spec v0.4

**Status:** a specification only. There is no script and nothing has been run. It is registered as a pre-registration document in the commit that adds this file, on top of 94362ec.

**Supersedes:** v0.3 (512a1047…). The third review of 2026-11-21 asked for one last change:
- **C6:** PC4/PC5 run at **both amplitudes**, m_a = 0.005 and 0.01, for all seven configurations and both widths. σ_num takes the largest relative difference over this whole set, and also the difference of a₂ itself (§6, §6.1).

Nothing else changed.

**Earlier:** v0.2 (9171ff68…) was conditionally approved by the second review with two changes:
- **C4:** PC4/PC5 now also run at w = 200, so both widths used in the extrapolation are checked for convergence (§6, §6.1).
- **C5:** T1 is described as a diagnostic only. The interaction with preparation transients is stated as a potential systematic effect that is not bounded (§4.3, §6, §8).

A scope paragraph on three separate lines of work is added to §1 (S1). Nothing else changed.

**Earlier:** v0.1 (08ec719d…) was conditionally approved with three methodological corrections:
- **C1:** PC4/PC5 (time step and grid) are now run at **every** primary SOC point and every Galilean control, not only at S2 (§6).
- **C2:** σ_num has an explicit formula. All terms are relative and summed linearly; numerator and denominator are combined linearly for the ratio (§6.1).
- **C3:** the claim about the pulse-only subtraction is narrowed to what it guarantees (§4.3, §8).

Not changed: the prediction R = 1 + 12Q, the three primary points, and every threshold in §6–§7.

**Base:** 94362ec (`SOC_DETUNING_AND_PULSE_CROSSING_REPORT_v0.2`, §II.4 and §III.4).

**Workflow after approval:**
1. script and tests;
2. static review;
3. permission;
4. **one** run on Lenovo, from a fresh clean LF clone, with the input hashes verified first;
5. results report.

**The criteria in §7 are fixed by this document and are not changed after any result.**

## 1. Question

In the dispersionless mean-field hydrodynamics, a weak right-moving sound probe crosses a localised left-moving simple wave u = εS(x₀/w), with ∫S = 0. Its arrival-time shift is

δT = −(1 + 12Q − 5τ²) · ∫u² dx₀ / (32 c³ M²)  (registered, 94362ec).

At δ = 0 (τ = 0), the ratio to a Galilean BEC with the same M and the same μ is

**R ≡ δT_SOC / δT_Gal(M, μ) = 1 + 12Q,  Q = 8μM²c₄,  Λ > 0 ⇔ Q > 1 ⇔ R > 13.**

**Question:** does R survive in the full 1D two-component Gross–Pitaevskii (GP) dynamics?

What the GP adds compared with the hydrodynamics:
- the quantum pressure (dispersion), i.e. finite pulse length against the healing length;
- the upper band and its virtual corrections, of order μ/Ω;
- the exact non-parabolic band at all orders;
- the real preparation of the pulse and the probe.

**Three separate lines (S1).** The project keeps three questions apart, and this test addresses **only the third link of the third line**:
1. **The method.** Does the matching/EFT method have independent, reproducible predictive power outside the SOC example? This test does **not** answer that.
2. **Λ.** What is Λ mathematically and physically? Can it be extracted or bounded from independent observables, without assuming the B4 interpretation? This test does **not** answer that either.
3. **The SOC chain.** It runs: contact coefficient C′ (sign of Λ in the SOC EFT, registered) → hydrodynamic prediction R = 1 + 12Q (registered, machine-checked) → **survival of R in GP dynamics (this test)**.

A PASS here would **not** count as validation of the method as a general tool, and **not** as evidence for the physical reality or universality of Λ.

**What a PASS would mean.** It would support the transfer of the quantitative relation to mean-field GP dynamics at δ = 0, g_↑↓ = g, in 1D.
- It does **not** test beyond-mean-field physics, because the GP is a mean-field theory.
- It is **not** an experiment.
- It is **not** evidence for B4 and does **not** change Stage 2 (NO MATCH).

## 2. Model and units

**Units:** ħ = 1, E_R = k_R = 1 (so m = 1/2). The GP is quasi-1D along x, with periodic boxes.

**SOC GP** (gauge-transformed Raman form, δ = 0, SU(2)-symmetric g):

i∂ₜψ↑ = [(−i∂ₓ + 1)² + g n] ψ↑ + (Ω/2) ψ↓
i∂ₜψ↓ = [(−i∂ₓ − 1)² + g n] ψ↓ + (Ω/2) ψ↑,  n = |ψ↑|² + |ψ↓|², Ω = 4x.

**Galilean GP:** i∂ₜψ = −∂ₓ²ψ/(2M) + g|ψ|²ψ.

**Background:** uniform, n₀ = 1, g = μ = γ (so g n₀ = γ E_R), at the band minimum κ₀ = 0.
- The SOC case has M = m* = x/(2(x−1)) and c = √(μ/m*).
- The Galilean control uses the **same M and μ**, hence the same c.

## 3. Configurations

| id | system | x | γ | M | role |
|---|---|---|---|---|---|
| G0 | Galilean | – | 0.2 | 0.5 | method validation (absolute) |
| S1 / G1 | SOC / Galilean(m*) | 1.10 | 0.2 | 5.5 | primary, Q = 9.091, R_th = 110.1 |
| S2 / G2 | SOC / Galilean(m*) | 1.20 | 0.2 | 3.0 | primary, Q = 2.083, R_th = 26.00 |
| S3 / G3 | SOC / Galilean(m*) | 1.30 | 0.2 | 2.1667 | primary, Q = 0.8547, R_th = 11.26 (Λ < 0 side) |
| E1 / GE1 | SOC / Galilean(m*) | 1.20 | 0.1 | 3.0 | diagnostic only (μ/Ω halved; Q = 1.042, R_th = 13.50) |

Q(x) = γ/(2x(x−1)²) is the exact value at δ = 0 for this band. The three primary points lie on both sides of R = 13.

## 4. Observable and its measurement

**4.1 Initial states.** All initial states use the local-density form with the exact lower-band spinor.

- **Pulse** (left-moving simple wave), supported on [x_P, x_P + w], with S(ξ) = sin(2πξ) sin²(πξ):
  - u(x) = εS((x − x_P)/w);
  - θ(x) = ∫u dx, which returns to 0 after the pulse because ∫S = 0;
  - n(x) = y(u(x))/g, with y from dy/du = −√(y Ẽ″(u)), y(0) = μ. This is integrated numerically with the exact band, as in the registered pulse-crossing run.
  - For SOC: ψ_σ = √n e^{iθ} χ₋,σ(u(x)), with χ₋(κ) the lower-band eigenvector of the 2×2 kinetic matrix at quasimomentum κ.
- **Probe** (linear right-moving wave), centred at x_pr = x_P − 1.5w, with width w_p = w/2:
  - δn = η n₀ sin²(π(x − x_pr + w_p/2)/w_p) on its support;
  - δu = g δn/c, δθ = ∫δu;
  - η = 1e-4.

**4.2 Runs per (configuration, w).** For each amplitude:
- **probe-only:** probe on the uniform background;
- **pulse-only(ε):** pulse without the probe;
- **pulse+probe(ε):** both.

**4.3 Probe signal and shift.**
- The probe signal is Δn(x, t_f) = n_{pulse+probe} − n_{pulse-only}.
  - **What the subtraction guarantees (C3):** it removes the stand-alone contribution of the pulse and of its preparation transients to the density.
  - **What it does not remove:** any interaction of those transients with the probe. That stays in the measured signal and is **not** removed exactly.
  - The long pulse and the extrapolation in w reduce the sensitivity to finite-size effects. The odd part (§4.4) and T1 (§6) are diagnostic indicators only (C5).
  - **Any residual contribution from this interaction remains a potential systematic effect.** It is not bounded by this spec and is not part of σ_num.
- The reference is Δn_ref = n_{probe-only} − n₀.
- Both are evaluated at the same final time t_f, inside a window of half-width 3w_p around the expected probe position x_pr + c t_f. The left-moving pulse and any reflected wave are outside the window by construction.
- The centroid is X = ∫x Δn dx / ∫Δn dx, with ΔX = X_{with pulse} − X_{ref}, and

  **δT(ε, w) ≡ −ΔX / c.**

- t_f is the time at which the probe centre is 2w beyond the trailing edge of the pulse, computed from c: t_f = 2.25 w/c. The window then stays at least 0.5w away from the pulse.

**4.4 Amplitude extraction.**
- The pulse amplitude is fixed by its peak Mach number m_a = (ε/M)/c, i.e. ε = m_a M c. This keeps the nonlinearity, and so the no-breaking margin, equal across configurations.
- m_a ∈ {±0.005, ±0.01}. S_x and G_x have the same M and c, so they use identical ε.
- The even part δT₊(ε) = [δT(ε) + δT(−ε)]/2 removes all odd orders.
- a₂(w) is obtained by Richardson extrapolation in ε²: with ε₁ (m_a = 0.005) and ε₂ = 2ε₁, a₂ = [4 δT₊(ε₁)/ε₁² − δT₊(ε₂)/ε₂²]/3.
- The odd part is reported as a diagnostic. Its O(ε) term must vanish because ∫S = 0.

**4.5 Dispersionless limit.**
- w ∈ {100, 200}. The healing length is ξ = 1/√(2Mμ) ≤ 2.3 in all configurations. The effective dispersive length ξ√(1+Q) is also ≤ 2.3 (largest: G0, 2.24; S1, 2.14).
- With A(w) ≡ a₂(w)/(w ∫S² dξ), the extrapolation A_∞ = [4A(200) − A(100)]/3 assumes O((ξ/w)²) corrections.
- The extrapolation uncertainty is σ_w ≡ |A_∞ − A(200)|.

**4.6 Primary quantity:**

**R_GP(x) = A_∞(S_x) / A_∞(G_x).**

The hydrodynamic predictions are A_∞(G) = −1/(32c³M²) and R_th = 1 + 12Q(x).

## 5. Numerics

- **Method:** Strang split-step Fourier.
  - The kinetic step is exact in k-space. For SOC this is the analytic exponential of the 2×2 matrix at each k.
  - The nonlinear step is exact: a common phase e^{−i g n dt}, which is exact because the interactions are SU(2)-symmetric.
- **Grid:** box length L = 16w, N = 2^13 for w = 100 and 2^14 for w = 200, so dx ≈ 0.195. Time step dt = 0.05. The layout keeps every wave at least 2w from the box edge until t_f, so nothing wraps around.
- **Runtime:** roughly 10–20 h in total on Lenovo, CPU only, including PC4/PC5 at both widths and both amplitudes. The dx/2 runs at w = 200 use N = 2^15. This is an estimate, not a measurement. The single invocation may run unattended overnight. This is to be measured, not guaranteed. A single invocation processes all configurations in sequence and writes one JSON per configuration plus a summary.
- **Environment:** Lenovo `stage2` conda environment, NumPy and SciPy pinned and hash-verified as in Stage 2.

## 6. Positive controls and gates (all evaluated in the same run)

| id | control | criterion |
|---|---|---|
| PC1 | probe speed in the probe-only runs (centroid drift) | equals c within 1e-3 (relative) |
| PC2 | probe linearity: S2, w = 100, m_a = 0.01, with η and 2η | ΔX equal within 1% |
| PC3 | no wave breaking while the probe is inside the pulse: the hydrodynamic bound max|λ₋′| t_exit over the pulse, for the largest m_a. Estimate ≈ 13 m_a ≈ 0.13 | < 0.2 for every configuration (checked before the runs) |
| PC4 | time step (C1, C4, C6): **each of S1, S2, S3, G0, G1, G2, G3**, at **both w = 100 and w = 200** and **both m_a = ±0.005 and ±0.01**, at dt and dt/2 | δT₊ equal within 1% for every configuration, width and amplitude; a₂ equal within 1% for every configuration and width |
| PC5 | grid (C1, C4, C6): the same configurations, widths and amplitudes, at dx and dx/2 | the same criteria |
| PC6 | lower-band fidelity (SOC): maximum upper-band population over the run | < 1e-4 |
| PC7 | conservation: norm drift / energy drift | < 1e-9 / < 1e-6 (relative) |
| **V0** | **method validation, G0:** A_∞ against −1/(32c³M²) | within 2% |
| **V1** | **each control G1–G3 and GE1:** A_∞ against −1/(32c³M²) | within 2% |

**Diagnostic T1 (reported, no threshold, not a gate).** In every pulse-only run: the largest |n − n₀| inside the probe window at t_f, in units of η n₀. It indicates whether transients are present where the probe is read. **It does not measure their interaction with the probe and gives no bound on it** (C5).

### 6.1 σ_num (C2), exact definition

For each configuration K ∈ {S_x, G_x}, define four **relative** contributions. Every absolute difference is divided by the magnitude of the quantity it perturbs:

- r_w(K) = |A_∞(K) − A(K; 200)| / |A_∞(K)| — the extrapolation in w;
- r_R(K) = max over w ∈ {100, 200} of |a₂(K; w) − δT₊(K; ε₁, w)/ε₁²| / |a₂(K; w)| — the Richardson step;
- r_dt(K) is the largest of these (PC4 of K itself, C4, C6):
  - |δT₊^{dt/2}(K, w, ε) − δT₊^{dt}(K, w, ε)| / |δT₊^{dt}(K, w, ε)|, over w ∈ {100, 200} and ε ∈ {ε₁, ε₂};
  - |a₂^{dt/2}(K, w) − a₂^{dt}(K, w)| / |a₂^{dt}(K, w)|, over w ∈ {100, 200}.

  The a₂ term covers the different weights with which the two amplitudes enter a₂.
- r_dx(K) is defined in the same way for the grid (PC5 of K itself).

Then:
- σ(K) = r_w(K) + r_R(K) + r_dt(K) + r_dx(K), a **linear, conservative** sum (no quadrature);
- **σ_num(x) = σ(S_x) + σ(G_x)**, the linear propagation of relative errors for the ratio R_GP = A_∞(S_x)/A_∞(G_x), to first order;
- σ_num(x) is dimensionless and is compared directly with D(x), which is also relative.

No contribution is taken from a configuration other than K itself (C1).

**What σ_num covers.** It estimates only the numerical contributions listed above. It is not a guarantee that every physical systematic is covered. In particular it excludes the transient–probe interaction (C5) and the upper-band (μ/Ω) effect, which is assessed only by diagnostic E1.

## 7. Decision rules (fixed now)

Let D(x) ≡ |R_GP(x)/R_th(x) − 1|, for x ∈ {1.10, 1.20, 1.30}.

| outcome | condition |
|---|---|
| **INCONCLUSIVE (method)** | any of PC1–PC7, V0 or V1 fails, or σ_num(x) > 5% at any primary x. **No statement about R is made.** |
| **PASS** | the method gates pass and **D(x) ≤ 10% at all three primary points** |
| **FAIL** | the method gates pass and at **at least one** primary point D(x) > 25% **and** D(x) > 3σ_num(x) |
| **INCONCLUSIVE** | every other case (some 10% < D ≤ 25%, or D > 25% but within 3σ_num) |

**Why these tolerances.**
- The EFT omits virtual upper-band corrections of relative order μ/Ω ≈ 4–5% at γ = 0.2. Q multiplies small differences, so the 10% band allows about twice that.
- 25% is far outside any effect identified in §8.
- V0 and V1 use 2% because the Galilean case has no upper band and its extrapolation in w is the only systematic.

**Reported descriptively** (not part of the decision): whether R_GP crosses 13 between x = 1.2 and 1.3; the odd parts; the absolute A_∞(S_x) against the hydrodynamic bracket; and E1/GE1.

## 8. Attribution of deviations (item 4 of the review)

| source | how it is separated |
|---|---|
| finite pulse length and healing length (quantum pressure) | the w dependence (100 vs 200) and its extrapolation; σ_w enters σ_num |
| numerical dispersion and splitting | PC4 (dt) and PC5 (dx); they enter σ_num |
| preparation transients and spurious waves | their stand-alone density is removed by the pulse-only subtraction. Their interaction with the probe is **not** removed (C3) and **not bounded** (C5). The long pulse and the w extrapolation reduce the sensitivity; the odd part and T1 are diagnostics. It remains a potential systematic effect outside σ_num |
| nonlinearity beyond O(ε²) | the ±ε symmetrisation and the Richardson step; the step enters σ_num |
| upper band, i.e. the EFT reduction (μ/Ω) | **diagnostic E1**: at the same x, γ = 0.1 halves μ/Ω. A deviation D that roughly halves from S2 to E1 points to the upper band. This is not used in the decision. |
| beyond mean field | **outside the scope**: the GP cannot test it, and it is stated as untested |

## 9. Outputs

- One JSON per configuration: all δT(ε, w), a₂, A, A_∞ and the PC values.
- `summary.json`: R_GP, R_th, D, σ_num, the decision and the descriptive items.
- stdout, stderr, RUN_CONTEXT (machine, versions, input hashes, commit), SHA256SUMS.
- Registration afterwards uses the usual .ps1, selecting files by SHA-256.

## 10. Not in scope

- δ ≠ 0;
- g_↑↓ ≠ g;
- 2D/3D or trap geometry;
- finite temperature;
- beyond-mean-field corrections;
- any B4 statement.
