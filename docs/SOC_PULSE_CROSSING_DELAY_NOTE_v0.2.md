# Time delay of a sound probe crossing a counter-propagating pulse: consistent treatment — note v0.2 (pre-registration)

**Status:** hand derivation and the plan of **one** check run. Nothing has been run. Not registered.

**Supersedes:** v0.1 (36a3f84a…), which the review of 2026-11-21 marked REVISION REQUIRED.

**Changes in v0.2 (review items):**
- **R1 (g₂).** λ₋ is now kept to second order, λ₋ = −c + g₁u + g₂u². §4 shows the g₂ terms explicitly.
  - The review proposed an additional term +g₂/(2c²). By the derivation below and by a structural argument, the g₂ contributions **cancel**.
  - The disagreement is not settled by argument. A new synthetic test, **B0**, measures the g₂ coefficient directly, and decision row D4 (§6) covers both possible answers.
- **R2 (A4).** A4 now derives f₁, f₂, g₁ **and g₂** explicitly.
- **R3 (B1).** The numeric part compares with the formula including any g₂ term. It builds its solution only from the exact λ±(u, y(u)) and the numerically integrated simple wave; it never uses the symbolic coefficients.
- **R4 (PC3).** The no-crossing condition is checked over the **whole pulse** and the whole time interval up to probe exit, not only along the probe path.
- **R5.** D1–D3 are unchanged. D4 is added for the g₂ question only.

**Base:** 0f9752c. Also refers to the unregistered `SOC_DETUNING_TOLERANCE_NOTE_v0.1` (49508fdb) and its run (stdout 1ae64edd).

**Purpose.** To replace the amplitude-level prescription of `SOC_BEC_QUARTIC_NOTE_v0.3` §7 by a well-defined observable and a consistent calculation, and to check whether its conclusions survive.

**Warning, stated before any check.** The hand derivation below (§4) suggests that, with a consistent treatment, even an **ordinary Galilean BEC** shows a time **advance** in this observable. If that is confirmed, the time-delay sign is **not** a proxy for the sign of Λ in systems with a cubic vertex. The amplitude-level reading of v0.3 §7 would then have to be withdrawn. The decision table (§6) is fixed now, before the run.

## 1. What is wrong with the v0.3 §7 prescription

v0.3 §7 computed the speed c₊ of a right-mover on a **constant-gradient** background with θₜ = cε and θₓ = ε. It then read the delay from the O(ε²) coefficient of 1/c₊, as for a **static** inhomogeneity: Δt = ∫(1/c₊ − 1/c) dx.

For a localised left-moving pulse in a system with a cubic vertex, three effects enter at the same order O(ε²) and are missing:

- **(a) The pulse relation.** A localised pulse is a *simple wave*. The relation between its density and its gradient is fixed by the Riemann invariant. It is not θₜ = cε. The difference is O(ε²) in the density, so O(ε²) in c₊.
- **(b) The crossing kinematics.** The probe and the pulse move relative to each other at λ₊ − λ₋, and that speed depends on the local u = θₓ at O(ε). The probe spends more time where the relative speed is lower. The O(ε) speed shift times this O(ε) dwell-time bias gives O(ε²).
- **(c) The pulse deformation.** Each level of the simple wave moves at its own speed λ₋(u).

In a relativistic P(X) theory with an exactly null background, all three are absent at this order. That is why the prescription works there. **The B4 Lagrangian has no cubic term either. The SOC phonon EFT has one, w q²/(2gm*), and the Galilean BEC has the same structure.**

## 2. The observable

**Setting:** the dispersionless 1D hydrodynamics of v0.3 §1–§2, for pulses much longer than the healing length:

∂ₜu + ∂ₓ(y + Ẽ(u)) = 0,  ∂ₜy + ∂ₓ(y Ẽ′(u)) = 0,  u = θₓ, y = g n.

**Characteristic speeds:** λ± = Ẽ′(u) ± √(y Ẽ″(u)). The unperturbed state is u = 0, y = μ, c² = μ/M.

**Pulse.** A localised **left-moving simple wave** of amplitude ε, with ∫u dx = 0, i.e. no net phase jump. On it, the Riemann invariant of the λ₊ family is constant. This gives

dy/du = −√(y Ẽ″(u)),  y(0) = μ.

**Probe.** A weak right-moving disturbance, linear in its own amplitude. To linear order it travels exactly along the λ₊ characteristic of the pulse field. This is exact because the λ₊ invariant of the pulse is uniform.

**Observable.** The shift δT of the probe's arrival time at a far detector, compared with no pulse, to O(ε²). A positive δT is a delay; a negative one is an advance.

## 3. Method

**Pulse.** Write the pulse in Lagrangian labels x₀ (its profile at preparation time t = 0). Each level u = S(x₀) moves on the straight line x = x₀ + λ₋(S(x₀)) t. This is exact until the characteristics cross.

**Probe.** Its position obeys dX/dt = λ₊(S(x₀(t))), where x₀(t) is the label it is crossing at time t. This gives

dt/dx₀ = (1 + λ₋′(x₀) t) / (λ₊ − λ₋).

The arrival shift is δT = −ΔX/c, with ΔX = w − 2c (t_exit − t_entry) for a pulse supported on labels [0, w].

## 4. Hand derivation (candidate, to be checked)

**Notation.** On the simple wave, write
- λ₊ = c + f(u), with f = f₁u + f₂u² + …;
- λ₋ = −c + g(u), with g = g₁u + g₂u² + … (g₂ kept, R1).

**Expansion.** Write h ≡ f − g and G ≡ g(S(x₀)), so that λ₊ − λ₋ = 2c + h and λ₋′ = G′. Then

dt/dx₀ = (1 + G′t)/(2c + h),  with t = t₁ + x₀/(2c) + τ₁(x₀) + …,

where t₁ is the entry time and τ₁ = O(ε).

**Where g₂ enters, and why it cancels (R1).** g₂ appears exactly twice at O(ε²):
- in the first-order term −h/(2c), through h ⊃ −g₂u²;
- in the deformation term G′·x₀/(2c), through G ⊃ g₂u².

After integration by parts (G vanishes at the pulse edges), these first-order pieces combine to

(1/(2c)) ∫ [−h/(2c) + G′x₀/(2c)] dx₀ = −(1/(4c²)) ∫ (h + G) dx₀ = −(1/(4c²)) ∫ f dx₀.

This identity is exact in g, because h + G = f. **g₁ and g₂ drop out of this term entirely.** g₁ survives only through the products of two O(ε) factors (h², G′τ₁, hG′x₀). These are quadratic in u, so g₂ would contribute there only at O(ε³).

**Structural argument.** If λ₊ ≡ c identically (f ≡ 0), the probe moves at c everywhere and δT = 0 exactly, however the pulse moves or deforms. So **every g-dependent term must carry a factor of f**. A term +g₂/(2c²) with no factor f would violate this.

**Results.** With ∫u dx₀ = 0 and the pulse vanishing at the edges:
- The O(ε) term is −(f₁/(2c²))∫u dx₀ = 0.
- The terms that depend on the preparation time t_entry cancel at O(ε²), because ∫u u′ dx₀ = 0.
- What remains is

  **δT = [ (2f₁² − f₁g₁)/(8c³) − f₂/(2c²) ] ∫u² dx₀ + O(ε³).**

**Galilean BEC** (Ẽ = u²/(2M), so c₃ = c₄ = 0). The Riemann invariant gives y = M(c − u/(2M))² exactly. Then
- λ₊ = c + u/(2M), so f₁ = 1/(2M) and f₂ = 0;
- λ₋ = −c + 3u/(2M), so g₁ = 3/(2M).

Hence 2f₁² − f₁g₁ = −1/(4M²) and

**δT = −∫u² dx₀ / (32 c³ M²) < 0, i.e. an ADVANCE.**

**Kinematic picture.** Where u > 0, the probe is faster (f > 0). There the relative speed 2c − 3u/(2M) + u/(2M) = 2c − u/M is **lower**, so the probe dwells longer in the regions where it is faster. The net effect is an advance, and the sign is fixed by kinematics, not by any contact coefficient.

**If confirmed, this contradicts the v0.3 §7 sanity limit.** That limit said "c₄ = 0 gives a delay, consistent with K3a" (5/8 > 0). Under the consistent treatment the ordinary BEC gives an advance. The v0.3 prescription would then have produced the "healthy" Galilean sign only because its background was not a pulse.

**SOC.** f₁, f₂ and g₁ follow from Ẽ = u²/(2M) + c₃u³ + c₄u⁴ and the simple-wave relation. By hand, with y₂ ≡ y″(0)/2:

y = μ − c u + y₂ u² + O(u³),  y₂ = 1/(4M) − (3/2) c₃ √(μM).

The full f₂, which depends on c₄, is left to the symbolic check. The SOC-specific information sits in the dependence of δT on c₄ (and on c₃ when δ ≠ 0), measured **against the Galilean baseline**.

## 5. Planned check (one run, after review and permission)

**Script:** `soc_pulse_crossing_check.py`. SymPy and NumPy/SciPy, `python3 -I`, one run, local workspace (as before).

**Part A, symbolic:**
- **A1:** the eigenvalues λ± of the system in §2, and the simple-wave relation dy/du = −√(yẼ″) (right eigenvector of λ₋).
- **A2:** the Galilean reduction. y = M(c − u/(2M))² solves A1 exactly, and f₁ = 1/(2M), f₂ = 0, g₁ = 3/(2M).
- **A3:** the series y(u) to O(u³) for general M, c₃, c₄. This includes y₂ above.
- **A4:** f₁, f₂, g₁ **and g₂** for general M, c₃, c₄, μ (R2), and the bracket of the §4 formula. The v0.3 quantity (s₁² − c s₂)/c³ is printed next to it for comparison.
- **A5:** a symbolic re-derivation of the §4 expansion, done independently of the hand derivation:
  - expand dt/dx₀ to O(ε²) with **generic** symbols f₁, f₂, g₁, g₂;
  - use a polynomial test profile u = ε·P(x₀) on [0, 1] with P(0) = P(1) = 0 and ∫P = 0, for example P = x₀(1−x₀)(1−2x₀);
  - integrate exactly;
  - report the coefficient of g₂ in δT, which §4 predicts is 0, and compare the rest with the §4 bracket.

**Part B, numeric** (independent of the §4 expansion: exact simple wave plus probe ODE). Per R3, the solution is built only from the exact λ±(u, y(u)), with y(u) integrated numerically from dy/du = −√(yẼ″). The symbolic coefficients are used only in the final comparison.

- **B0 (synthetic kinematics, new; decides R1).** Replace the physical λ± by prescribed polynomials λ₊ = c + f₁u + f₂u² and λ₋ = −c + g₁u + g₂u², with chosen numbers. Integrate the same probe ODE and extract δT/ε² by Richardson extrapolation. Do this for:
  - **(i)** f ≡ 0 with (g₁, g₂) = (1.5, 0) and (1.5, 2.0). Prediction: δT = O(ε³), so the extrapolated δT/ε² is 0 to 1e-6.
  - **(ii)** (f₁, f₂, g₁) = (0.5, 0.3, 1.5) with g₂ ∈ {−2, 0, 2}. Prediction: δT/ε² does not depend on g₂, and equals the §4 bracket to 1e-3 relative.

  The fitted slope of δT/ε² against g₂ is printed together with the review's predicted slope, (1/(2c²))∫P² dx₀ over the pulse support (pulse u = εP(x₀)).
- **Pulse:** smooth and compact, u(x₀) = ε·sin(2πx₀/w)·sin²(πx₀/w) on [0, w], so ∫u dx₀ = 0.
- **B1:** solve dt/dx₀ = (1 + λ₋′t)/(λ₊ − λ₋) with the exact λ±(u) (with y(u) from the exact ODE of A1). Do this for ε ∈ {1e-2, 5e-3, 2.5e-3}, extract δT/ε² by Richardson extrapolation, and compare with the §4 formula. Criterion: relative error < 1e-3.
- **Cases:**
  - Galilean (M = 1, μ = 1);
  - SOC at δ = 0, (γ, x) = (0.2, 1.2) and (0.2, 1.10);
  - SOC at δ = 0.05 E_R, (0.2, 1.2), both directions, the left-moving probe being obtained by mirror.
- **B2:** the same with two preparation times, t_entry = 0 and t_entry = w/(4c). The difference must be O(ε³).
- **B3:** the scan x = 1.02…1.50 for γ ∈ {0.1, 0.2, 0.5} at δ = 0. Print:
  - the §4 bracket;
  - the Galilean bracket with the same M and μ;
  - their difference, which is the SOC-specific part;
  - the v0.3 prescription value.

**Positive controls (if any fails, the run is FAIL and is not interpreted):**
- **PC1:** with no pulse (ε = 0), δT = 0 to machine precision.
- **PC2:** a pulse with ∫u ≠ 0, u = ε sin²(πx₀/w), gives δT linear in ε with coefficient −(f₁/(2c²))∫(u/ε)dx₀, to 1e-3.
- **PC3 (R4):** no crossing of the pulse characteristics, 1 + λ₋′(x₀)t > 0, checked for **every** label x₀ of the pulse (a grid of 4001 labels) and **every** t from 0 to the probe's exit time. The minimum over this whole region is printed for every case; a non-positive minimum is a FAIL.

## 6. Decision table (fixed before the run)

| outcome | consequence |
|---|---|
| **D1:** the Galilean bracket is negative (advance), confirmed by A and B | The time-delay sign of this observable is **not** a proxy for Λ in systems with a cubic vertex. The amplitude-level reading of v0.3 §7, and the amplitude columns of the detuning run, are **withdrawn as an interpretation**. They remain as computed numbers of a prescription. **The contact-level result (Λ > 0, C′ sign, detuning tolerance 0.16 E_R at γ = 0.2, x = 1.2) is unaffected.** The GP spec then targets a **quantitative** test: δT(x) against this formula, including its c₄-dependent part relative to the Galilean baseline. It does not target a sign. |
| **D2:** the Galilean bracket is positive (delay), i.e. §4 contains an error | The corrected formula replaces both v0.3 §7 and §4. The SOC sign question is re-read from B3 with the same rule as v0.3 (advance where the bracket is negative). |
| **D3:** A and B disagree with each other | FAIL: nothing is interpreted, and the cause is found before any further step. |
| **D4 (R1):** B0 shows a non-zero g₂ slope | The §4 formula is incomplete. The term measured by B0 is added. The A5 symbolic expansion must reproduce it, otherwise D3 applies. D1/D2 are then read with the corrected formula; the Galilean case is unaffected because g₂ = 0 there. If B0 confirms a zero slope, the §4 formula stands as written. |

**Note on K2.** In case D1, a time-delay-based bound like the provisional K2 bound would not discriminate in non-relativistic media with a cubic vertex, because an ordinary BEC already shows an advance. This is consistent with §9.3 of v0.3: no causality principle forbids it. It would weaken the K2 link of the SOC result. It would not weaken the contact-level coefficient result.

## 7. What this note does not do

- No dispersion: the quantum pressure and the k³ terms of the phonon branch are omitted, so the pulses must be much longer than the healing length.
- No g_↑↓ ≠ g.
- No beyond-mean-field corrections.
- No Feynman 2→2 amplitude.
- The GP simulation remains the independent physical test.
