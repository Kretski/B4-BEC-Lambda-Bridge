# Raman SOC BEC: tolerance of the Λ > 0 result to the detuning δ — note v0.1 (pre-registration of a symbolic/numeric check)

**Status:** a hand derivation, plus the plan of **one** check run. Nothing has been run. Not registered.

**Base:** 0f9752c (`SOC_BEC_QUARTIC_NOTE_v0.3`, open item 2, validity row V6).

**What this note does:**
1. It shows that the reduction of v0.3 §1 stays exact for δ ≠ 0.
2. It derives the δ ≠ 0 generalisation of the contact condition (v0.3 §2) and of the amplitude condition (v0.3 §7). This is by hand.
3. It fixes, **before** any computation, the following:
   - the check script's content;
   - the pass criteria;
   - the positive controls;
   - the categories used to interpret the result.

**What it does not do:**
- no Gross–Pitaevskii simulation;
- no g_↑↓ ≠ g;
- no beyond-mean-field corrections;
- no Feynman-diagram amplitude.

The amplitude level uses the same eikonal/background method as v0.3 §7, with the same caveat.

## 1. The reduction stays exact at δ ≠ 0

**Model** (quasi-1D, ħ = 1):

H₀(p) = [[ (p + k_R)²/2m − δ/2 , Ω/2 ], [ Ω/2 , (p − k_R)²/2m + δ/2 ]].

**Lower band.** With a ≡ Ω/2, b ≡ k_R/m and D(p) ≡ b p − δ/2:

E₋(p) = p²/2m + E_R − R(p),  R ≡ √(D² + a²).

This is the same convention as the registered dipole check (`experiments/soc_bec/dipole_check_v0.1`, κ = p/k_R, energies in E_R).

**Why the reduction of v0.3 §1 still holds:**
- **(i)** H₀(p) is real and symmetric for every δ. So the lower-band spinor χ₋(p) can be chosen real, and no Berry-phase term arises.
- **(ii)** The interactions are SU(2)-symmetric (g_↑↑ = g_↓↓ = g_↑↓ = g). The interaction energy then depends only on the total density, not on the polarisation, which now depends on both p and δ.
- **(iii)** The upper band stays gapped by at least 2a = Ω.

Hence, as in v0.3 §1–§2, the result is still exact at leading derivative order:

**P(Y) = Y²/(2g),  Y = μ − ∂ₜθ − Ẽ(∂ₓθ),  Ẽ(u) ≡ E₋(p₀ + u) − E₋(p₀).**

Here p₀ is the (single) minimum of E₋. For x ≡ Ω/(4E_R) > 1 the band has one minimum for every δ.

**What changes.** E₋ is no longer even about its minimum:

Ẽ(u) = u²/(2M) + c₃ u³ + c₄′ u⁴ + O(u⁵).

The three coefficients are evaluated at the shifted minimum p₀(δ). At the minimum the group velocity E₋′(p₀) is 0, so the ground state carries no current.

## 2. Local band coefficients at the shifted minimum (exact)

**Derivatives of R:**
- dR/dp = bD/R;
- d²R/dp² = b²a²/R³;
- d³R/dp³ = −3b³a²D/R⁵;
- d⁴R/dp⁴ = −3b⁴a²(a² − 4D²)/R⁷.

The last form uses R² − 5D² = a² − 4D².

**Exact coefficients, with D₀ ≡ D(p₀) and R₀ ≡ R(p₀):**
- **1/M** = 1/m − b²a²/R₀³;
- **c₃** = b³a²D₀/(2R₀⁵);
- **c₄′** = b⁴a²(a² − 4D₀²)/(8R₀⁷).

**Minimum condition:** p₀/m = b D₀/R₀, with D₀ = b p₀ − δ/2.

**Exact sign statement.** The local quartic coefficient of the band **at its minimum is negative if and only if |D₀| > a/2 = Ω/4.**
- This explains the registered §10 observation of v0.3.
- At Ω = 4.8 E_R and δ = 1 E_R: κ₀ = −0.559 gives D₀ = 2κ₀ − δ/2 ≈ −1.62 E_R, while a/2 = 1.2 E_R.

**Leading order in δ:**
- From the minimum condition, D₀(1 − m b²/R₀) = −δ/2.
- Also m b²/a = 4E_R/Ω = 1/x.
- Hence

  **d ≡ D₀/a ≈ −δ/(4E_R (x − 1)).**

The effective small parameter is therefore δ/(x − 1), not δ. **The tolerance shrinks towards the transition x → 1**, which is exactly where the Λ > 0 window of v0.3 lies.

**To O(d²):**
- 1/M ≈ (1/m*)(1 + 3d²/(2(x−1)));
- c₃ ≈ b³d/(2a²);
- c₄′ ≈ c₄(1 − 15d²/2), with c₄ = b⁴/(8a³) as in v0.3.

## 3. The conditions with the cubic term

**Contact level.** Y²/(2g) gives:
- the pure quartic C′ = 1/(8gM²) − n c₄′ (c₃ does not enter it);
- a **new pure-spatial cubic vertex** −n c₃ u³, which B4 does not have.

Define **Q ≡ 8 g n M² c₄′.** Then Λ > 0 ⇔ **Q > 1**. At δ = 0 this is v0.3 §2: Q = (γ/2)/(x(x−1)²).

**Amplitude level** (eikonal, as in v0.3 §7).
- The exact speed of a right-mover on a left-moving background (θₜ = cε, θₓ = ε) is unchanged in form:

  c₊ = Ẽ′(ε) + √(g n₀(ε) Ẽ″(ε)),  g n₀ = μ − cε − Ẽ(ε),  c² = μ/M.

- Write c₊ = c + s₁ε + s₂ε². Then:
  - s₁ = 1/(2M) + t/2, with t ≡ 6c₃μ/c;
  - c s₂ = 6c₄′μ − 1/(4M²) − (t − 1/M)²/8.
- The O(ε²) coefficient of 1/c₊ is K₂ = (s₁² − c s₂)/c³. This gives

  **K₂ c³ M² = 5/8 + τ/4 + 3τ²/8 − (3/4) Q,  τ ≡ tM = 6 c₃ M √(μM).**

**Direction dependence.** The mirror x → −x maps Ẽ(u) → Ẽ(−u), i.e. c₃ → −c₃. So a left-mover crossing a right-moving pulse has τ → −τ.

**Advance condition, direction σ = ±1:**

**(3/4) Q > 5/8 + σ τ/4 + 3τ²/8.**

**Checks and reading:**
- At τ = 0 this gives Q > 5/6, which is v0.3 B5: (48/5) g n c₄ m*² > 1. ✓
- The **τ² term always adds delay**: cubic exchange is "healthy".
- The **linear term has opposite signs in the two directions.** At small δ it is the dominant correction, because τ = O(d) while the Q corrections are O(d²).
- The tolerance is set by the **worse** direction.

**Leading order of τ** (with μ = γE_R and m* = m x/(x−1)):

τ ≈ (3/√2) d √(γ/x) (x − 1)^{−3/2} ≈ −(3/(4√2)) (δ/E_R) √(γ/x) (x − 1)^{−5/2}.

## 4. Rough hand estimate (to be replaced by the exact computation)

**Inputs:** γ = 0.2, x = 1.2.
- δ = 0: Q = 2.083, so (3/4)Q = 1.56 against 5/8. The margin is 0.94.
- τ ≈ −12 δ/E_R.

**Two estimates:**
- **Linear term only:** the worse direction loses the sign at about |τ| ≈ 1.28, i.e. **|δ| ≈ 0.1 E_R**.
- **Including the O(d²) drop of Q:** at |δ| = 0.1 E_R, d ≈ 0.125, so Q falls by roughly 35%. This moves the estimate to **|δ| ≈ 0.05–0.1 E_R**.

**Reading.** At this d the O(d²) expansion is already marginal. **The estimate is not a result.** The exact, non-perturbative computation of §5 Part B decides.

**Scale.** In the Zhang et al. experiment E_R = 2π × 2.21 kHz, so 0.1 E_R ≈ 2π × 220 Hz.

## 5. The planned check (one run, after review and permission)

**Script:** `soc_detuning_tolerance_check.py`. SymPy plus NumPy/SciPy, `python3 -I`, one run.

**Part A, symbolic (exact identities):**
- **A1:** the derivatives of R in §2, including the form (a² − 4D²).
- **A2:** the minimum condition, and the leading-order d ≈ −δ/(4E_R(x−1)), checked by series in δ.
- **A3:** the O(d²) forms of 1/M, c₃ and c₄′.
- **A4:** from c₊ = Ẽ′(ε) + √(g n₀ Ẽ″(ε)) with symbolic M, c₃, c₄′, the O(ε²) coefficient of 1/c₊ equals the §3 formula for K₂.
- **A5:** at c₃ = 0, A4 reduces to the registered B3 of v0.3: (1/c)(5/(8Mμ) − 6c₄′M).
- **A6:** the mirror u → −u maps the right-mover formula to the left-mover formula with τ → −τ.
  - Implementation: compute the left-mover speed on a right-moving background directly, with no symmetry assumption.

**Part B, numeric (exact in δ, no expansion in δ):**
- **Grid:**
  - γ ∈ {0.1, 0.2, 0.5};
  - x from 1.02 to 1.50 in steps of 0.02;
  - δ/E_R ∈ [0, 0.5].
- **For each point:**
  - find p₀ as the root of E₋′ (bracketed; the uniqueness of the minimum is checked on a grid);
  - compute M, c₃ and c₄′ from the exact §2 formulas;
  - compute Q and τ.
- **Outputs:**
  - **B1:** for each (γ, x), the critical δ at which Q = 1 (contact level), found by bisection. If there is none in [0, 0.5], that is stated.
  - **B2:** the critical δ for each direction σ = ±1 (amplitude level), and the minimum of the two.
  - **B3:** a printed table of B1 and B2.
  - **B4:** the local-quartic sign switch, i.e. the δ at which |D₀| = a/2, for comparison.

**Positive controls (must pass, otherwise the run is reported as FAIL and not interpreted):**
- **PC1:** at δ = 0 the code reproduces the registered δ = 0 thresholds of v0.3 to 1e-3:
  - contact x < 1.204, 1.280, 1.420;
  - amplitude x < 1.222, 1.303, 1.454;
  - for γ = 0.1, 0.2, 0.5.
- **PC2:** at Ω = 4.8 E_R, δ = 1 E_R the code reproduces the registered dipole check (stdout 3584bd7f…) to 1e-4. The quantities are κ₀ = −0.559…, ε″, ε‴ and ε⁗.
- **PC3:** the exact Part B values agree with the §2–§3 O(d²) forms at small δ (δ = 1e-3 E_R, relative error < 1e-3).

**Interpretation categories, fixed now, before any number exists.**

Let δ_tol be the smallest critical δ of B2 (worse direction, amplitude level) at γ = 0.2, x = 1.2, the centre of the v0.3 window. Then:
- **robust:** δ_tol ≥ 0.1 E_R;
- **sensitive:** 0.01 E_R ≤ δ_tol < 0.1 E_R;
- **fragile:** δ_tol < 0.01 E_R. That is about 2π × 22 Hz at the E_R of Zhang et al.

**The categories describe experimental fragility only.** None of them changes the mean-field sign result at δ = 0. No category stops the project. "Fragile" would move the planned GP spec to include δ as a primary scan variable, not a nuisance parameter.

**Run conditions.**
- As for the v0.3 symbolic check: a local workspace, `python3 -I`, and a RUN_CONTEXT with the script hash, versions, exit code and the stdout/stderr hashes.
- Registration afterwards uses the usual .ps1 by SHA-256.
- If you want this run on Lenovo from a clean clone instead, say so in the review.

## 6. Open after this note

- g_↑↓ ≠ g (v0.3 §5, item 3).
- The GP time-delay spec, which will use the δ_tol from this check.
- The two-field structure against B4.
