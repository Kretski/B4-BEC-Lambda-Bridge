# Raman spin–orbit-coupled BEC: sign of the one-derivative quartic (Λ) — note v0.2

**Status:** analytic note for review. Hand derivations; no code and no run. Not registered. It does not change Stage 2, K2, K3b or the CONTACT note.

**Base:** 30f882a. **Supersedes:** v0.1 (26a572ab…).

**Added in v0.2:**
- §6: the authorised symbolic check of §1–§2.
- §7: the amplitude-level combination, including exchange through the cubic vertex.

§1–§5 are unchanged.

**Motivation.**
- Lemma L2 (`K3B_MEAN_FIELD_LEMMA_v0.2`) gives Λ < 0 for every stable Galilean mean-field mixture. The premise that carries it is the **parabolic** single-particle kinetic energy.
- The nonlocal photon-fluid note showed that a mediator which only renormalises the local couplings cannot change this.
- A Raman SOC condensate breaks the parabolic premise already at mean field. This note asks whether that is enough for Λ > 0.

**Stopping rule, fixed before the derivation.** Stop if Λ > 0 requires a regime in which any of these fails:
- the projection onto the lower band, gn ≪ Ω;
- the quartic expansion;
- stability of the uniform state.

## 0. Model

- **Coordinates and units:** quasi-1D along x, ħ = 1. Raman coupling Ω, detuning δ = 0, recoil k_R, E_R = k_R²/(2m).
- **Single-particle bands** (standard form, e.g. the Raman-SOC reviews arXiv:2608.08121 and arXiv:2108.10159):

  E_±(p) = p²/2m + E_R ± √((Ω/2)² + (k_R p/m)²).

- **Phase:** single-minimum phase, Ω > 4E_R. Write x ≡ Ω/(4E_R) > 1.
- **Interactions:** spin-symmetric, g_↑↑ = g_↓↓ = g_↑↓ = g > 0 (SU(2)-symmetric). This is approximately the case for ⁸⁷Rb. The case g_↑↓ ≠ g is an open item (§5).

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

## 2. Elimination of the density (as in K3b)

n is algebraic: n = Y/g. Hence **P(Y) = Y²/(2g)**.

**Expansion of the band.** With a = Ω/2 and b = k_R/m:

Ẽ(q) = q²/(2m*) + c₄ q⁴ + O(q⁶),
- m/m* = 1 − 4E_R/Ω = 1 − 1/x;
- c₄ = b⁴/(8a³) = 4E_R²/(m²Ω³) **> 0**.

**Pure spatial quartic** (one derivative per field, coefficient of (∂ₓθ)⁴ in L). It comes from Y²/(2g), with Y ⊃ μ − q²/(2m*) − c₄q⁴:

**C = 1/(8 g m*²) − n c₄.**

- The first term is the compressibility term of L2. It is positive.
- The second term is new. It is **negative**, because the band is non-parabolic.

**Sign of Λ.** Λ ∝ −C (as in K3a and K3b). So **Λ > 0 ⇔ 8 g n m*² c₄ > 1.** With γ ≡ gn/E_R this becomes

**(Λ > 0)  γ/2 > x (x − 1)².**

**Examples:**

| γ | Λ > 0 for x below about |
|---|---|
| 0.1 | 1.2 |
| 0.2 | 1.28 (i.e. Ω ≲ 5.1 E_R) |
| 0.5 | 1.42 |

## 3. The stopping-rule checks

| check | condition | status at γ = 0.2, x ≈ 1.1–1.28 |
|---|---|---|
| V1 band projection | gn ≪ Ω, i.e. γ ≪ 4x | 0.2 ≪ 4.4–5.1: **satisfied** |
| V2 radius of the q-expansion | b\|q\| < a, i.e. \|q\| < x·k_R | the coefficient C is a Taylor coefficient and is well defined. The backgrounds must stay at \|q\| well below k_R; this bounds the B4 backgrounds that can be represented |
| V3 uniform state | single-minimum phase for x > 1 with SU(2)-symmetric g: the uniform p = 0 state is the ground state | **satisfied** (x > 1). The phonon speed c² = gn/m* softens as x → 1, which agrees with the measured softening (Ji et al., arXiv:1408.1755) |
| V4 thermodynamic stability | g > 0, i.e. P_YY = 1/g > 0 | **satisfied** |
| V5 dilute limit, mean field | √(na³) ≪ 1 | assumed |

**Result of the stopping rule.** It is **not triggered**. Λ > 0 holds in a controlled, stable window 1 < x ≲ 1.3 for γ ≈ 0.2, and the window widens with γ. It does not require the edge x → 1, where the phonon softens completely.

## 4. What this does and does not give for B4

**Gives (hand derivation):** a physical, stable mean-field system with **Λ > 0**. This is the first one in the project.

**It is outside every exclusion so far:**
- C1 assumes relativistic scalar mediators;
- K1 assumes a Lorentz-invariant UV;
- K3a, K3b (L1, L2) and the photon-fluid note assume a parabolic kinetic term or a Galilean UV.

**Does not give the rest of B4.** For the single phonon field in P(Y) = Y²/(2g):
- **p = 0.** P is exactly quadratic in Y, so there is no w²q² term. The registered B4 backgrounds have p > 0.
- **A cubic term w q²/(2gm*) is present.** B4 has none.

The physically meaningful combination is therefore the **2→2 amplitude including exchange**, not the contact coefficients alone (K3b §7, `matching_protocol.md` §4).

Consequence for K2:
- If the amplitude-level combination has the sign of 3Λ > 0, this system has the "wrong" sign for the provisional K2 bound 2p + 3Λ ≤ 0. It would then be a **physical counterexample** to K2's premises for this boost-breaking class.
- This is **not established** here, because exchange has not been computed.

**Second field.** The gapped spin mode, with gap of order Ω at q = 0, is a natural ψ with a mass term, as in L_B4 ⊃ −½m²ψ². Its couplings and its speed relative to the phonon are not derived here.

## 5. Open items, in order of importance

1. **Symbolic check** of §1–§2: the band expansion, m*, c₄, C and the condition γ/2 > x(x − 1)².
2. **The 2→2 amplitude including the exchange from the cubic term w q².** It decides whether Λ > 0 survives at the physical level, and settles the K2 question.
3. **g_↑↓ ≠ g.** The interaction then couples to the spin polarisation s_z(q) ∝ q. That adds a q-dependent coupling and shifts the transition point. The effect on C has to be derived.
4. **The two-field structure** (phonon plus gapped spin mode) against the rest of B4: p, s, κ, rank one, speeds and β signs.

## References (checked for this note)

- Raman-SOC single-minimum phase for ħΩ_R ≥ 4E_R: arXiv:2608.08121 (review) and arXiv:1605.02108.
- Phonon softening near the transition, measured: S.-C. Ji et al., arXiv:1408.1755.
- Negative-mass hydrodynamics from the engineered SOC dispersion: M. A. Khamehchi et al. (Engels group), arXiv:1612.04055.
- Coherently coupled mixtures, single-minimum dispersion: arXiv:2108.10159.
- `K3B_MEAN_FIELD_LEMMA_v0.2.md` (L1, L2); `NONLOCAL_PHOTON_FLUID_VS_B4_NOTE_v0.1.md`; T7a.


## 6. Symbolic check of §1–§2 (authorised 2026-11-20)

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
  - **c₄ = 0** (parabolic band, an ordinary superfluid): the coefficient is positive. That means a time delay, the healthy sign, consistent with K3a. This sanity limit is reproduced. Check B4 in the script is only a placeholder, so this limit rests on reading B3 at c₄ = 0, not on B4.
  - **c₄ > 0** large enough: the coefficient is negative. That means a **time advance**, i.e. the "wrong" sign at amplitude level.
- **B5:** the advance condition is **(48/5) g n c₄ m*² > 1 ⇔ (3/5)γ > x(x − 1)².**

**Thresholds, compared:**

| γ | contact level: x below (γ/2 > x(x−1)²) | amplitude level: x below (3γ/5 > x(x−1)²) |
|---|---|---|
| 0.1 | 1.204 | 1.222 |
| 0.2 | 1.280 | 1.303 |
| 0.5 | 1.420 | 1.454 |

**Reading.**
- Including exchange **does not remove** the wrong sign. The amplitude-level window is slightly **wider** than the contact-level one.
- The stable single-minimum Raman-SOC condensate therefore has, at the amplitude level, the sign that the provisional K2 bound (2p + 3Λ ≤ 0) excludes.

**Status:** this is a symbolic result **under the eikonal/background method**. It is **not** a Feynman-diagram computation of the 1+1D 2→2 amplitude. In 1+1D, the t- and u-channel exchanges of massless quanta are collinear and need separate treatment.

**Consequence for K2, conditional on that method:**
- A physical, stable boost-breaking system shows the wrong sign. So the premises of the provisional K2 exclusion **do not hold** for this class: analyticity, crossing or boundedness fail.
- This is consistent with `K2_BOOST_BREAKING_2D_NOTE_v0.5` §3. That note found none of those premises established for a non-relativistic UV.

## 8. Updated status

| item | status |
|---|---|
| Λ > 0 at contact level, stable, controlled window | **symbolically confirmed** (Part A) |
| wrong sign at amplitude level, including exchange (eikonal) | **symbolically confirmed** (Part B) |
| full B4 structure | **not reproduced**: p = 0 (B4 needs p > 0); a cubic term is present; the second field is not derived |
| g_↑↓ ≠ g | open |
| Feynman-diagram 2→2 amplitude in 1+1D | open (cross-check of §7) |
