# Raman spin–orbit-coupled BEC: sign of the one-derivative quartic (Λ) — note v0.3

**Status:** analytic note with one registered symbolic check (§6) and one local numerical cross-check (§10). v0.2 is registered at d3100f0. v0.3 and the §10 artefacts (`experiments/soc_bec/dipole_check_v0.1/`) are registered in the commit that adds this file, on top of d3100f0. It does not change Stage 2, K2, K3b or the CONTACT note.

**Base:** d3100f0. **Supersedes:** v0.2 (18671d71…).

**Changes in v0.3.** No new derivation. Every formula of §1, §2, §6 and §7 is unchanged. The changes are:
- R1: the status line is corrected (v0.2 claimed "no code and no run, not registered").
- R2: prior art is stated. The quartic curvature of the lower band is known; see §9.
- R3: the meaning of "wrong sign" is made precise (§4, §7, new §9.3). It is a statement relative to the provisional K2 bound. It is not a causality violation.
- R4–R6: three validity rows are added to §3: detuning δ, phonon damping, and beyond-mean-field corrections.
- R7: a possible experimental anchor for c₄ is added (§10).
- R8: the open item g_↑↓ ≠ g is made concrete (§5).
- R9: references are added; the literature check is recorded in §9.

**Motivation.**
- Lemma L2 (`K3B_MEAN_FIELD_LEMMA_v0.2`) gives Λ < 0 for every stable Galilean mean-field mixture. The premise that carries it is the **parabolic** single-particle kinetic energy.
- The nonlocal photon-fluid note showed that a mediator which only renormalises the local couplings cannot change this.
- A Raman SOC condensate breaks the parabolic premise already at mean field. This note asks whether that is enough for Λ > 0.

**Stopping rule, fixed before the derivation (v0.1).** Stop if Λ > 0 requires a regime in which any of these fails:
- the projection onto the lower band, gn ≪ Ω;
- the quartic expansion;
- stability of the uniform state.

The validity rows added in v0.3 (V6–V8) are robustness conditions. They are **not** part of the stopping rule fixed in v0.1.

## 0. Model

- **Coordinates and units:** quasi-1D along x, ħ = 1. Raman coupling Ω, detuning δ = 0, recoil k_R, E_R = k_R²/(2m).
- **Single-particle bands** (standard form; e.g. Li, Martone and Stringari, arXiv:1410.5526, Eq. 3; the Raman-SOC reviews arXiv:2608.08121 and arXiv:2108.10159):

  E_±(p) = p²/2m + E_R ± √((Ω/2)² + (k_R p/m)²).

- **Phase:** single-minimum phase, Ω > 4E_R. Write x ≡ Ω/(4E_R) > 1.
- **Interactions:** spin-symmetric, g_↑↑ = g_↓↓ = g_↑↓ = g > 0 (SU(2)-symmetric). This is approximately the case for ⁸⁷Rb. The case g_↑↓ ≠ g is open (§5).
- **Geometry:** the claim is restricted to quasi-1D propagation along the SOC axis (see V7).

## 1. Elimination of the heavy band (the "mediator" step)

**Exact uniform-flow energy.**
- Take a uniform flowing state ψ = √n e^{iqx} χ₋(q), with χ₋ the lower-band spinor.
- With SU(2)-symmetric interactions, the interaction energy does not depend on the spin polarisation.
- The mean-field energy density is therefore **exactly** e(n, q) = n·Ẽ(q) + (g/2) n², with Ẽ(q) ≡ E₋(q) − E₋(0).

**Slowly varying fields.**
- The upper band is gapped by Ω at q = 0. Virtual excitation of it gives corrections suppressed by ω/Ω and gn/Ω.
- At δ = 0 the spinor χ₋(q) is real, so no Berry-phase term arises.
- The spatial variation of χ₋(q(x)) produces (∂ₓq)²-type terms. These are higher-derivative.

**Effective Lagrangian at leading derivative order:**

L = n·Y − (g/2) n²,  with Y ≡ μ − ∂ₜθ − Ẽ(∂ₓθ).

This is the T7a step, with the upper band in place of the Proca field. The only memory the heavy mode leaves at this order is the **non-parabolic** Ẽ(q).

**Relation to the literature.** The linear limit of this description is the standard SOC hydrodynamics of Li, Martone and Stringari (arXiv:1410.5526, §4.3, Eqs. 57–59). It uses the same two assumptions: the relative phase is locked because phonon frequencies are well below the gap ~Ω, and g_↑↓ = g. That hydrodynamics is **linear**. It keeps only m/m* and gives c² = gn/m*. The present note keeps the full non-parabolic Ẽ(q) in Y, which is what produces the nonlinear terms below.

## 2. Elimination of the density (as in K3b)

n is algebraic: n = Y/g. Hence **P(Y) = Y²/(2g)**.

**Expansion of the band.** With a = Ω/2 and b = k_R/m:

Ẽ(q) = q²/(2m*) + c₄ q⁴ + O(q⁶),
- m/m* = 1 − 4E_R/Ω = 1 − 1/x (agrees with arXiv:1410.5526, Eq. 5, where 2k₀² = 4E_R);
- c₄ = b⁴/(8a³) = 4E_R²/(m²Ω³) **> 0**.

**Prior art (R2).** The quartic curvature of the lower band is **not** new. Li, Martone and Stringari state that at Ω = 4E_R the spectrum near the minimum has a p_x⁴ dispersion (arXiv:1410.5526, below Eq. 5). The anharmonicity beyond the effective-mass approximation has been observed in dipole oscillations (Zhang et al., PRL 109, 115301 (2012), arXiv:1201.6018). The general-x coefficient c₄ is an elementary Taylor coefficient. What this note adds is **where c₄ enters**: the competition in C below, the threshold, and the sign of the sound time delay (§7). That combination was not found in the literature checked (§9).

**Pure spatial quartic** (one derivative per field, coefficient of (∂ₓθ)⁴ in L). It comes from Y²/(2g), with Y ⊃ μ − q²/(2m*) − c₄q⁴:

**C = 1/(8 g m*²) − n c₄.**

- The first term is the compressibility term of L2. It is positive.
- The second term comes from the non-parabolic band. It is **negative**.

**Sign of Λ.** Λ ∝ −C (as in K3a and K3b). So **Λ > 0 ⇔ 8 g n m*² c₄ > 1.** With γ ≡ gn/E_R this becomes

**(Λ > 0)  γ/2 > x (x − 1)².**

**Examples:**

| γ | Λ > 0 for x below about |
|---|---|
| 0.1 | 1.2 |
| 0.2 | 1.28 (i.e. Ω ≲ 5.1 E_R) |
| 0.5 | 1.42 |

## 3. Validity checks

V1–V5 are the stopping-rule checks of v0.1/v0.2. V6–V8 are added in v0.3.

| check | condition | status at γ = 0.2, x ≈ 1.1–1.28 |
|---|---|---|
| V1 band projection | gn ≪ Ω, i.e. γ ≪ 4x | 0.2 ≪ 4.4–5.1: **satisfied** |
| V2 radius of the q-expansion | b\|q\| < a, i.e. \|q\| < x·k_R | the coefficient C is a Taylor coefficient and is well defined. The backgrounds must stay at \|q\| well below k_R; this bounds the B4 backgrounds that can be represented |
| V3 uniform state | single-minimum phase for x > 1 with SU(2)-symmetric g: the uniform p = 0 state is the ground state | **satisfied** (x > 1). The phonon speed c² = gn/m* softens as x → 1, which agrees with the measured softening (Ji et al., arXiv:1408.1755) |
| V4 thermodynamic stability | g > 0, i.e. P_YY = 1/g > 0 | **satisfied** |
| V5 dilute limit, mean field | √(na³) ≪ 1 | **assumed.** Caution (R6): quantum fluctuations are enhanced by SOC and modify the sound velocity and the superfluid density (Liang and Törmä, PRA 100, 023619 (2019), arXiv:1903.08182). Near x ≈ 1–1.3 the plain √(na³) criterion may not be sufficient. Not quantified here |
| V6 detuning (R4) | δ = 0 is assumed throughout | **open, and sensitive.** For δ ≠ 0 the band is asymmetric: the minimum moves and odd terms appear. In the local check of §10, at Ω = 4.8 E_R and δ = 1 E_R the minimum sits at κ₀ = −0.56 k_R and the **local quartic coefficient there is negative**. The sign result of §2 and §7 therefore needs a tolerance on δ. Not derived |
| V7 phonon damping (R5) | the phonon must be long-lived on the scale of the delay measurement | **favourable along x.** For G₂ = 0 and propagation along the SOC axis (θ = 0), the zero-temperature Beliaev rate carries the factor (m/m*)² relative to the SOC-free rate, so it is **suppressed** towards the transition (Wu and Liang, PRL 121, 180401 (2018), arXiv:1802.10295, Eqs. 5 and 7). For θ ≠ 0 the rate diverges at the phase boundary. The claim is therefore restricted to quasi-1D propagation along x. Finite-temperature Landau damping is not assessed (Liang and Törmä find it dominant) |
| V8 band projection with g_↑↓ ≠ g | see §5 item 3 | open |

**Result of the stopping rule.** It is **not triggered**. Λ > 0 holds in a controlled, stable window 1 < x ≲ 1.3 for γ ≈ 0.2, and the window widens with γ. It does not require the edge x → 1, where the phonon softens completely. The added rows V5–V8 do not trigger it either, but V5 and V6 are open and V6 is sensitive.

## 4. What this does and does not give for B4

**Gives:** a physical, stable mean-field system with **Λ > 0**, confirmed symbolically (§6). It is the first such system in the project.

**It is outside every exclusion so far:**
- C1 assumes relativistic scalar mediators;
- K1 assumes a Lorentz-invariant UV;
- K3a, K3b (L1, L2) and the photon-fluid note assume a parabolic kinetic term or a Galilean UV.

**Does not give the rest of B4.** For the single phonon field in P(Y) = Y²/(2g):
- **p = 0.** P is exactly quadratic in Y, so there is no w²q² term. The registered B4 backgrounds have p > 0.
- **A cubic term w q²/(2gm*) is present.** B4 has none.

The physically meaningful combination is therefore the **2→2 amplitude including exchange**, not the contact coefficients alone (K3b §7, `matching_protocol.md` §4).

**Meaning of "wrong sign" (R3).** In this note "wrong sign" means: the sign that the **provisional K2 bound** 2p + 3Λ ≤ 0 excludes. It does **not** mean a violation of causality. See §9.3: in a non-relativistic medium, microcausality refers to the light cone. A sound time advance relative to the sound speed is not excluded by it.

Consequence for K2: see §7.

**Second field.** The gapped spin mode, with gap of order Ω at q = 0, is a natural ψ with a mass term, as in L_B4 ⊃ −½m²ψ². Its couplings and its speed relative to the phonon are not derived here.

## 5. Open items, in order of importance

1. **Independent cross-check of the amplitude sign (§7).** Either a Feynman-diagram computation of the 1+1D 2→2 amplitude, or a numerical Gross–Pitaevskii time-delay experiment (planned spec `SOC_GP_TIME_DELAY_SPEC_v0.1`). The spec must include the δ tolerance of item 2.
2. **Detuning tolerance (V6).** How large may |δ| be before the sign of C, and of the §7 coefficient, changes?
3. **g_↑↓ ≠ g.** The interaction then couples to the spin polarisation s_z(q) ∝ q. This adds a q-dependent coupling and shifts the transition point. Starting point from the literature: with the spin-density interaction G₂ (∝ g − g_↑↓), m/m* = 1 − 2k₀²/(Ω + 4G₂) and the phase boundary is 4G₂ + Ω = 2k₀² (ħ = m = 1; Wu and Liang, Eqs. 5 and 8; Li, Martone and Stringari). The G₂-dependence of c₄ and of C is not available and has to be derived.
4. **The two-field structure** (phonon plus gapped spin mode) against the rest of B4: p, s, κ, rank one, speeds and β signs.
5. **Beyond-mean-field corrections near the transition (V5).**

## References (checked for this note)

- Raman-SOC single-minimum phase for ħΩ_R ≥ 4E_R: arXiv:2608.08121 (review) and arXiv:1605.02108.
- Y. Li, G. I. Martone and S. Stringari, *Spin-orbit-coupled Bose–Einstein condensates*, review, arXiv:1410.5526 (bands, m/m*, p⁴ at the transition, linear hydrodynamics).
- J.-Y. Zhang et al., *Collective dipole oscillation of a spin-orbit coupled BEC*, PRL 109, 115301 (2012), arXiv:1201.6018 (anharmonicity beyond the effective mass, observed).
- Phonon softening near the transition, measured: S.-C. Ji et al., arXiv:1408.1755.
- R. Wu and Z. Liang, *Beliaev damping of a spin-orbit-coupled BEC*, PRL 121, 180401 (2018), arXiv:1802.10295.
- L. Liang and P. Törmä, *Quantum corrections to a spin-orbit coupled BEC*, PRA 100, 023619 (2019), arXiv:1903.08182.
- Q. Zhu, C. Zhang and B. Wu, *Exotic superfluidity in spin-orbit coupled BECs*, EPL 100, 50003 (2012), arXiv:1109.5811 (Rashba SOC; frame-dependent critical velocity).
- Negative-mass hydrodynamics from the engineered SOC dispersion: M. A. Khamehchi et al., arXiv:1612.04055.
- M. E. Mossman et al., *Shock waves in a superfluid with higher-order dispersion*, arXiv:2004.00832.
- P. Creminelli, O. Janssen and L. Senatore, *Positivity bounds on EFTs with spontaneously broken Lorentz invariance*, arXiv:2207.14224.
- L. Hui, A. Nicolis, A. Podo and S. Zhou, *Microcausality without Lorentz invariance*, arXiv:2502.04215.
- L. Hui et al., *S-matrix positivity without Lorentz invariance: a case study*, arXiv:2312.08440.
- Coherently coupled mixtures, single-minimum dispersion: arXiv:2108.10159.
- `K3B_MEAN_FIELD_LEMMA_v0.2.md` (L1, L2); `NONLOCAL_PHOTON_FLUID_VS_B4_NOTE_v0.1.md`; T7a.

## 6. Symbolic check of §1–§2 (authorised 2026-11-20; registered d3100f0)

**Script and run.**
- Script: `soc_quartic_and_amplitude_check.py`, sha256 8d7abd61….
- One run with `python3 -I` (SymPy 1.14.0): rc 0. stdout 4f218cd4…, RUN_CONTEXT 3666eca9….
- **13/13 PASS.** Part A has 7 checks; Part B has 6 and belongs to §7.

**Part A confirms:**
- m/m* = 1 − 4E_R/Ω;
- c₄ = 4E_R²/(m²Ω³) > 0;
- no odd powers of q at δ = 0;
- C = 1/(8g m*²) − n c₄;
- **p = 0**: there is no w²q² term in Y²/(2g);
- the cubic term wq²/(2g m*) is present;
- 8g n m*² c₄ = (γ/2)/(x(x−1)²).

**§2 is therefore confirmed as written.**

## 7. Amplitude level, including exchange through the cubic vertex

**Method.** The background (eikonal) method of `K2_BOOST_BREAKING_2D_NOTE_v0.5` §1 is used:
- Take a constant-gradient **left-moving** background, θₜ = cε and θₓ = ε (a left-moving wave of the linear theory).
- Find the **exact** speed c₊ of a right-moving fluctuation from the quadratic operator of the full P(Y). This includes the cubic vertex wq² exactly at the level of the background.
- The time delay of a right-mover crossing a localised left-moving pulse comes from the expansion of 1/c₊.
  - The **O(ε)** term is odd in the background. Its integral across a localised pulse is ∫ε = Δθ = 0.
  - The delay is therefore set by the **O(ε²)** coefficient of 1/c₊. It includes the square of the O(ε) shift, which is the exchange effect at this level.

**Exact speed.** c₊ = Ẽ′(ε) + √(g n₀ Ẽ″(ε)), with n₀ = (μ − cε − Ẽ(ε))/g and c² = μ/m*.

**Results (Part B, all PASS):**
- **B2:** the O(ε) shift of c₊ is ε/(2m*).
- **B3:** the O(ε²) coefficient of 1/c₊ is **(1/c)·(5/(8 m* μ) − 6 c₄ m*)**.
  - **c₄ = 0** (parabolic band, an ordinary superfluid): the coefficient is positive. That means a time delay, consistent with K3a. This sanity limit is reproduced. Check B4 in the script is only a placeholder, so this limit rests on reading B3 at c₄ = 0, not on B4.
  - **c₄ > 0** large enough: the coefficient is negative. That means a **time advance** relative to propagation at the unperturbed sound speed.
- **B5:** the advance condition is **(48/5) g n c₄ m*² > 1 ⇔ (3/5)γ > x(x − 1)².**

**Thresholds, compared:**

| γ | contact level: x below (γ/2 > x(x−1)²) | amplitude level: x below (3γ/5 > x(x−1)²) |
|---|---|---|
| 0.1 | 1.204 | 1.222 |
| 0.2 | 1.280 | 1.303 |
| 0.5 | 1.420 | 1.454 |

**Reading.**
- Including exchange **does not remove** the sign. The amplitude-level window is slightly **wider** than the contact-level one.
- The stable single-minimum Raman-SOC condensate therefore has, at the amplitude level, the sign that the provisional K2 bound (2p + 3Λ ≤ 0) excludes.

**Status:** this is a symbolic result **under the eikonal/background method**. It is **not** a Feynman-diagram computation of the 1+1D 2→2 amplitude. In 1+1D, the t- and u-channel exchanges of massless quanta are collinear and need separate treatment.

**Consequence for K2, conditional on that method:**
- A physical, stable boost-breaking system shows the sign that K2 excludes. So the premises of the provisional K2 exclusion **do not hold** for this class: analyticity, crossing or boundedness fail.
- This is consistent with `K2_BOOST_BREAKING_2D_NOTE_v0.5` §3. That note found none of those premises established for a non-relativistic UV.
- It is **not** a causality violation (§9.3).

## 8. Updated status

| item | status |
|---|---|
| Λ > 0 at contact level, stable, controlled window | **symbolically confirmed** (Part A) |
| time advance at amplitude level, including exchange (eikonal) | **symbolically confirmed** (Part B) |
| quartic band curvature itself | **prior art** (arXiv:1410.5526; observed indirectly, arXiv:1201.6018) |
| competition in C, threshold, sign of the sound delay | **not found in the literature checked** (§9) |
| full B4 structure | **not reproduced**: p = 0 (B4 needs p > 0); a cubic term is present; the second field is not derived |
| δ tolerance | **open, sensitive** (V6) |
| g_↑↓ ≠ g | open (starting point in §5) |
| Feynman-diagram or GP cross-check of §7 | open |

## 9. Literature check (2026-11-20)

**9.1 Scope.** The question was whether the effect is already described: the quartic curvature of the lower band, its competition with compressibility in the flow-quartic coefficient, and the sign of the sound time delay.

Read in full (sha256 of the PDFs used):
- arXiv:1109.5811 (f1998cb0);
- arXiv:1410.5526 (7896ce16);
- arXiv:2207.14224 (eeddf031);
- arXiv:2502.04215 (a6d5b885);
- arXiv:1201.6018 (0563e431);
- arXiv:1802.10295 (6918c128).

arXiv:2004.00832 was read through a web summary only, not the full PDF.

**9.2 Findings on the SOC side.**
- **Known:**
  - the bands and m/m*;
  - p⁴ at the transition;
  - the linear hydrodynamics with m*;
  - the anisotropic sound speed and its softening;
  - Galilean breaking and frame-dependent critical velocities (Rashba case);
  - observed dipole anharmonicity beyond the effective mass;
  - Beliaev damping with SOC;
  - shock waves with higher-order dispersion (at Ω = 1.5 E_R, δ = 0.54 E_R, i.e. in the double-minimum regime);
  - quantum corrections.
- **Not found:**
  - the flow-quartic coefficient C of the phonon EFT;
  - its sign change and the threshold γ/2 > x(x−1)²;
  - the sign of the time delay of a phonon crossing a counter-propagating pulse.
- **Caveat:** absence in a finite search is not proof of absence.

**9.3 Findings on the positivity/causality side (R3).**
- **Creminelli, Janssen and Senatore (arXiv:2207.14224).**
  - The bounds come from two-point functions of conserved currents and assume a conformal UV.
  - The authors state that they cannot constrain operators that start cubic or higher in fluctuations.
  - The term C·Y⁴ around the ground state is quartic in fluctuations. It lies outside these bounds, and a BEC does not have a conformal UV.
- **Hui, Nicolis, Podo and Zhou (arXiv:2502.04215).**
  - Microcausality is formulated with respect to the light cone.
  - For a superfluid phonon it bounds the k² dispersion coefficient. The bound becomes loose in the non-relativistic limit c_s ≪ 1 (their Eq. 6.39).
  - Quartic interactions do not enter their two-point analysis.
- **Hui et al. (arXiv:2312.08440).** For S-matrix positivity without Lorentz invariance the question is stated to be open.
- **Consequence.** No result found forbids a sound time advance in a non-relativistic BEC. The "wrong sign" of §7 is wrong only relative to the provisional K2 bound.

## 10. Local numerical cross-check: dipole period vs amplitude (authorised 2026-11-20; artefacts in `experiments/soc_bec/dipole_check_v0.1/`)

**Purpose.** Test whether the two-level lower band E₋ alone reproduces the amplitude-dependent dipole period of Zhang et al. (arXiv:1201.6018, Fig. 2b; Ω = 4.8 E_R, δ = 1 E_R). Model: their Eq. 2, k̇ = −mω²x, ẋ = ∂E₋/∂k. The perturbative m_F = +1 correction used by the authors is **not** included.

**Runs.**
- Attempt 1 aborted. Script a3e66ac5…; empty stdout and stderr. The integration span 40π/ω₀ was infinite at Ω = 4 E_R, δ = 0, where ω₀ = 0.
- Attempt 2: script a5eb205e… (fixed: chunked span). rc 0, stdout 3584bd7f…, stderr empty.
- Environment: local workspace, Python 3.13, NumPy 2.5.3, SciPy 1.18.1. Not a clean-clone run.

**Results.**

| case | κ₀ | cubic part of dω/A² | quartic part | period vs amplitude A |
|---|---|---|---|---|
| Ω = 4.8, δ = 1 (Fig. 2b) | −0.559 k_R | −0.083 | −0.033 | grows; T(0.9)/T(0.3) ≈ 1.07–1.08 |
| Ω = 4.8, δ = 0 | 0 | 0 | +0.266 | falls; T/T₀ = 2.43 → 1.84 for A = 0.1 → 0.9 |
| Ω = 4.0, δ = 0 | 0 | 0 | pure k⁴ | finite, T ∝ 1/A approximately |

**Reading.**
- The growing period of Fig. 2b is reproduced in shape. The figure gives a ratio of about 1.10, read off by eye.
- The absolute period is about 20% shorter than measured: about 31 ms against about 38 ms at A ≈ 0.3. Likely reasons are the missing m_F = +1 level and the unknown trap period at Ω = 4.8.
- **Fig. 2b does not test the c₄ of this note.** At δ = 1 the oscillation is around κ₀ = −0.56 k_R, where the local quartic coefficient is negative and the cubic term dominates.
- At δ = 0 the derivatives reproduce the §2 formulas exactly: ε⁗/24 = 16(E_R/Ω)³ = 0.14468 and ε″/2 = 1/6 = 1 − 4E_R/Ω.
- **Possible experimental anchor (R7).** A dipole amplitude scan at δ = 0, Ω > 4E_R should show a period that **falls** with amplitude. This would test the single-particle c₄ only. The interaction term 1/(8gm*²) of C does not enter dipole motion.
