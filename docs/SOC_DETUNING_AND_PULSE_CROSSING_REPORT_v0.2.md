# SOC BEC: detuning tolerance and pulse-crossing delay — results report v0.2

**Status:** report on two pre-registered, reviewed, single runs (2026-11-20 and 2026-11-21), plus one requested symbolic check of a hand step (§II.4). It is registered in the same commit as the two notes and the three run folders, on top of 0f9752c.

**Changes from v0.1 (1f9fd1a9…, not registered).**
- The hand rewrite in §II.4 is now machine-checked (run `bracket_rewrite_check_v0.1`, R1–R4 PASS).
- Two points from the 2026-11-21 review are stated explicitly: (i) PASS here is not independent reproduction; (ii) the GP prediction is a model prediction, not an experimental confirmation.

**Scope.** This report records the results, which earlier readings are withdrawn, and what is left open. It contains no new derivation. The one rewrite of a printed symbolic result (§II.4) is machine-checked.

| item | file | sha256 |
|---|---|---|
| detuning note (pre-registration) | `docs/SOC_DETUNING_TOLERANCE_NOTE_v0.1.md` | 49508fdb… |
| detuning run | `experiments/soc_bec/detuning_tolerance_v0.1/` | script 3bb7eff7…, stdout 1ae64edd… |
| pulse-crossing note (pre-registration, reviewed: v0.1 REVISION REQUIRED, v0.2 APPROVED) | `docs/SOC_PULSE_CROSSING_DELAY_NOTE_v0.2.md` | c2349f9d… (supersedes v0.1, 36a3f84a…, not registered) |
| pulse-crossing run | `experiments/soc_bec/pulse_crossing_v0.2/` | script 0566b8b3…, stdout d00758a4… |
| check of the §II.4 rewrite | `experiments/soc_bec/bracket_rewrite_check_v0.1/` | script c67f23cb…, stdout 3a9d4c4e… |

**All three runs:**
- one attempt each, exit code 0, stderr empty;
- local workspace, not a clean-clone run;
- Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1, SymPy 1.14.0;
- **ALL CHECKS: PASS.**

**What PASS means here.** ALL_CHECKS_PASS is a statement about the pre-registered checks of each script. It is not independent proof. All three runs were made by the same author in the same local workspace, not from a clean clone, and none has been reproduced by anyone else. Their scripts, outputs and hashes are registered so that this can be done.

## Part I — Detuning tolerance (contact level)

**I.1 Exact reduction.**
- At δ ≠ 0 the single-particle Hamiltonian stays real and the interactions stay SU(2)-symmetric. So P(Y) = Y²/(2g) remains exact at leading derivative order.
- Only the band changes. Around its shifted minimum p₀(δ) it acquires a cubic term: Ẽ = u²/(2M) + c₃u³ + c₄′u⁴.
- Symbolic checks A1–A3: PASS.

**I.2 Exact local sign statement (A1).** The quartic coefficient of the band at its minimum is negative if and only if |D₀| > Ω/4. This explains the negative local quartic seen in the registered dipole check (`dipole_check_v0.1`) at Ω = 4.8 E_R, δ = 1 E_R.

**I.3 Small parameter.** d ≈ −δ/(4E_R(x − 1)) (A2). The tolerance therefore shrinks towards the transition x → 1, which is where the Λ > 0 window lies.

**I.4 Contact-level result** (Part B, exact in δ). Define Q′ ≡ 8gnM²c₄′; then Λ > 0 ⇔ Q′ > 1.

The critical detuning at which Q′ = 1, at γ = 0.2:

| x | 1.02 | 1.10 | 1.16 | 1.20 | 1.26 |
|---|---|---|---|---|---|
| δ_c / E_R | 0.108 | 0.161 | 0.174 | **0.164** | 0.098 |

For γ = 0.5 the contact tolerance reaches 0.29 E_R. The full tables are in the run's stdout.

**Positive controls (all passed):**
- PC1: the δ = 0 thresholds of v0.3 are reproduced to 1e-3;
- PC2: the registered dipole check is reproduced to 6 digits;
- PC3: agreement with the O(d²) forms within 2·10⁻⁴.

**I.5 Withdrawn (see Part II).** The run also printed amplitude-level columns: the "advance" thresholds for the two directions, and the pre-registered category "sensitive" (δ_tol = 0.094 E_R). These were computed with the v0.3 §7 prescription. Part II shows that this prescription is not a consistent treatment of a pulse crossing.
- **The amplitude columns and the category are withdrawn as an interpretation.** They remain in the registered stdout as numbers of that prescription only.
- The same applies to the direction asymmetry (σ = ±1) in §3 of the detuning note. Under the consistent treatment the two directions give the same result (§II.3).
- **The contact-level results I.1–I.4 are unaffected.**

## Part II — Pulse-crossing delay (consistent treatment)

**II.1 Problem with v0.3 §7.** For a localised pulse in a system with a cubic vertex, three effects enter at O(ε²) and were missing:
- the simple-wave relation of the pulse;
- the crossing kinematics, i.e. the dwell-time bias;
- the pulse deformation.

They are absent in a relativistic theory with a null background, and in B4, which has no cubic term. They are present in the SOC phonon EFT and also in an ordinary Galilean BEC.

**II.2 Formula.** In dispersionless 1D hydrodynamics (pulses much longer than the healing length), consider a weak right-moving probe crossing a localised left-moving simple wave with ∫u dx₀ = 0. With λ₊ = c + f₁u + f₂u² and λ₋ = −c + g₁u + g₂u² on the simple wave:

**δT = [ (2f₁² − f₁g₁)/(8c³) − f₂/(2c²) ] ∫u² dx₀ + O(ε³).**

**Evidence:**
- **A5:** an independent symbolic expansion with generic coefficients and a polynomial profile reproduces the bracket. The coefficient of g₂ is exactly 0, and there is no dependence on the preparation time.
- **B0, synthetic kinematics:**
  - with f ≡ 0, δT/ε² = O(10⁻¹¹);
  - the fitted slope against g₂ is 2.95·10⁻¹², against 0.078 for the additional term +g₂/(2c²) proposed in the first review.
  - **The g₂ question is settled: no g₂ term. The reviewer has withdrawn the proposal.**
- **B1, exact simple wave plus probe ODE**, built from the exact λ±(u, y(u)) only: agreement with the formula within 2·10⁻⁸ relative, in all five cases.
- **B2:** no dependence on the preparation time.
- **PC1–PC3:** all passed. The minimum of 1 + λ₋′t over the whole pulse up to probe exit is 0.89, so the characteristics never cross.

**II.3 Galilean BEC: decision D1.** For Ẽ = u²/(2M) the run gives f₁ = 1/(2M), f₂ = 0, g₁ = 3/(2M), g₂ = 0 (A2), hence

**δT_Gal = −∫u² dx₀ / (32 c³ M²) < 0, i.e. an advance.**

The symbolic and the numeric values agree (−0.03125 per unit ∫u² at M = μ = 1). The v0.3 prescription gives the opposite sign for the same case: +0.0977 against −0.00488 for the profile used.

**Decision D1 (pre-registered):**
- The sign of this time shift is **not** a proxy for the sign of Λ in systems with a cubic vertex.
- **The amplitude-level reading of `SOC_BEC_QUARTIC_NOTE_v0.3` §7 is withdrawn.** The v0.3 sanity limit ("c₄ = 0 gives a delay, consistent with K3a") holds only under that prescription.

**II.4 SOC: a quantitative handle on Q.**
- **What the run printed (A4).** The bracket is

  (180 M³ c₃² μ − 96 M² c₄ μ − 1) / (32 √M μ^{3/2}).

- **Rewrite, machine-checked.** With c³M² = μ^{3/2}M^{1/2}, Q = 8μM²c₄ (μ = gn) and τ = 6c₃M√(μM), this becomes

  **δT = −(1 + 12Q − 5τ²) · ∫u² dx₀ / (32 c³ M²).**

- **At δ = 0** (c₃ = 0, τ = 0), the ratio to the Galilean baseline with the same M and c is

  **R ≡ δT_SOC / δT_Gal = 1 + 12Q,**

  where Q is exactly the quantity with Λ > 0 ⇔ Q > 1. **Λ > 0 corresponds to R > 13.** The observable does not discriminate by sign, but it measures Q quantitatively.

- **Check of the rewrite** (`bracket_rewrite_check_v0.1`, requested before registration):
  - **R1:** the printed A4 bracket, taken as text from stdout d00758a4, is identically equal to the rewrite;
  - **R2:** the bracket recomputed from scratch (simple wave, f₁, f₂, g₁) is identically equal to the rewrite;
  - **R3:** at c₃ = 0, the ratio to the Galilean bracket with the same M and μ is exactly 1 + 12Q;
  - **R4:** the rewrite reproduces the four registered B1 "formula" values within 3·10⁻¹⁰ relative.
  - All four: PASS.

- **Direction.** For δ ≠ 0 the bracket depends on c₃ only through c₃². The two directions are therefore equal. This is confirmed numerically at δ = 0.05 E_R, where both give −0.714043.

- **Limits of the statement:**
  - It is derived for δ = 0 and for this effective theory: dispersionless, mean field, g_↑↓ = g, quasi-1D.
  - For δ ≠ 0 the sign of δT is set by 1 + 12Q − 5τ². Q′ can fall, and c₄′ can turn negative for |D₀| > Ω/4. So the advance is **not** claimed for all SOC parameters. It holds in all printed cases (B1, and B3 at δ = 0).

## Part III — Consequences and limits

**III.1 What changes:**

| statement | status |
|---|---|
| Galilean BEC shows an advance in this observable | **confirmed** (D1) |
| a g₂ term in the coefficient | **rejected** (B0, A5) |
| amplitude-level reading of v0.3 §7 | **withdrawn** |
| amplitude columns and the "sensitive" category of the detuning run | **withdrawn as an interpretation** |
| SOC allows Λ > 0 at contact level | **unaffected** |
| contact-level detuning tolerance ≈ 0.16 E_R at γ = 0.2, x = 1.2 | **unaffected** |
| quantitative relation R = 1 + 12Q at δ = 0 | **supported within the model** |
| full realisation of B4 | **not achieved** (p = 0, cubic term present, second field not derived) |
| Stage 2 | **NO MATCH** (unchanged) |
| link of the SOC result to the provisional K2 bound | **weakened** |

**III.2 K2.**
- A time-delay-based bound such as the provisional K2 bound does not discriminate in non-relativistic media with a cubic vertex, because an ordinary BEC already shows an advance.
- This is consistent with v0.3 §9.3: no causality principle forbids an advance relative to c_s.
- The contact coefficient sign remains a statement about the EFT coefficient. It is not about this observable.

**III.3 B4.** Nothing here is evidence for B4. B4 has no cubic vertex. The SOC phonon EFT has one, and lacks the w²q² term (p = 0) and a derived second field.

**III.4 Next independent test (not authorised yet).** A Gross–Pitaevskii simulation, kept as a separate post-registration step with its own pre-registered spec, a Galilean control and a quantitative criterion for R_GP. **R = 1 + 12Q is a prediction of the dispersionless mean-field model. It is not experimentally or numerically confirmed beyond that model.** The target is:

**R_GP(x) ≟ 1 + 12Q(x) at δ = 0,** with the Galilean BEC as the control.

The open conditions are:
- dispersion: pulse length against healing length;
- g_↑↓ ≠ g;
- beyond-mean-field corrections near x ≈ 1.

## Sources of the numbers

All numbers are quoted from the three registered stdout files. The rewrite in §II.4 is checked by `bracket_rewrite_check_v0.1`.
