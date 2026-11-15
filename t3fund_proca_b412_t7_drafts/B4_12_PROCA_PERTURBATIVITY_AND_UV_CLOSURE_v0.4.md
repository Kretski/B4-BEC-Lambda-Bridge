# B4.12 — Proca: perturbativity and UV closure (PROPOSAL v0.4)

**Status: PROPOSAL, NOT REGISTERED.** Nothing is executed before this document is reviewed and registered. It does not change:
- B4.9, B4.11 or the historical B4.12 (v0.5);
- T3-fund Stage 1 or C4, which remain as registered and approved;
- the committed B4.12 v0.3 (`B4_12_PROCA_PERTURBATIVITY_AND_UV_CLOSURE.md`), which is kept unchanged as the previous version;
- the approved T7a drafts and the T7a manuscript.

**Changes from v0.3**
1. Written in English, the language of the repository. v0.3 was in Bulgarian.
2. **Source of ρ.** The definition is now taken from the primary source, the evidence package committed at 7ae9c5c/511f6fe (§0).
3. **Coupling of ψ (W4).** The coupling is stated as the code defines it, \|g_ψ\| = M\|κ\|/√p. This is equal to M√s only under the rank-one condition (§1, item 1).
4. **Reading of the historical counts (W4).** The historical set counts are read as a necessary condition. The common-M result (commits 893d005, 0a7c7ab, 6901214) is added (§1, items 3–4).
5. **Q1** is answered as a table, without new computation, as the approved procedure requires (§2).
6. **Q2 and Q3** carry the decisions of the review of the working draft: v0.2/v0.2.1, review of 2026-11-15, decisions 1–7 (§3, §4, §5).
7. **Derivation (§0).** The derivation and conventions of v0.3 §3 Q2 item 2a are replaced by a reference to the approved T7a step 3, metric (−,+).
8. **Review of 2026-11-15, 15:15.**
   - Q3 item 2: the contact operator is written as it follows from the Proca matching, with coefficient κ²/(2p), and κ²/p = s is stated only for the rank-one case.
   - §0 Notation: the corresponding consistency edit, so that s is not identified with g_ψ²/M² unconditionally.

---

## 0. Conventions and definitions

**Derivation.** T7a step 3, v0.5.2 (commit 61fefff, approved as a draft), metric (−,+), D = ∂_t² − ∂_x²:
- (D + M²)A1 = (1 − ∂_x²/M²) J, with J = g_φ f_t f_x + g_ψ g_t g_x;
- exact tree-level kernel K(ω,k) = (1 + k²/M²)/(M² + k² − ω²);
- A1 = J/M² − ∂_t²J/M⁴ + O(M⁻⁶);
- L_eff = J²/(2M²) + (∂_tJ)²/(2M⁴) + O(M⁻⁶), after integration by parts.

v0.3 was written in (+,−). The results agree, and from v0.4 on only (−,+) is used.

**ρ (primary source).** `evidence/B4_10_proca_mass_sweep_v0_5/b410/uv.py`, lines 3–4 and 16–36 (SHA-256 `db32b8ec…b7b0c3be`):
- ρ_j = \|u_j\|/M for u_j ∈ {q, w, r, v}, where w = f_t, q = f_x, v = g_t, r = g_x;
- M is a single Proca mass, common to both fields;
- the couplings are g_φ = M√p (for q, w) and \|g_ψ\| = M\|κ\|/√p (for r, v);
- identity, component by component: g_j ρ_j = √X_j.

**Name clash.** `docs/STAGE2_PREREGISTRATION.md` uses ρ = a2/a1, a BEC parameter. That is a different quantity.

**Notation.** s denotes only the B4 EFT coefficient. In the Proca construction g_ψ²/M² = κ²/p, which equals s only under the rank-one condition. The kinematic (Mandelstam) variable is s_cm.

## 1. Established

1. **Fixed-EFT prescription** (`fixed_eft.couplings`):
   - g_φ = M√p, g_ψ = −Mκ/√p, so \|g_ψ\| = M\|κ\|/√p;
   - this equals M√s only under the rank-one condition κ² = ps;
   - on the 41 B4.12 backgrounds the deviation is δ_rank1 ≤ 2.2·10⁻¹⁶;
   - the identity g_j ρ_j = √X_j is verified to 2·10⁻¹⁶;
   - consequently g ∝ M and ρ ∝ 1/M, and the EFT data fix only the products.
2. **Historical control sets.** The thresholds come from `uv.SETS`:
   - A: g ≤ 1, ρ ≤ 0.3;
   - B: g ≤ 1, ρ ≤ 0.1;
   - C: g ≤ √(4π) (α ≤ 1), ρ ≤ 0.3.

   The historical code tests √X_max ≤ g_max·ρ_max. This is a **necessary** condition for membership.
3. **Historical counts, read as the necessary condition** on the 41 E3 backgrounds (strong-X by construction, from the cut-0.95 population): A 0/41, B 0/41, and the historical set-C count of 41/41 (29/29 without the negative controls). This statement concerns these 41 backgrounds only.
4. **Common-M result** (spec v0.1.2, 3cba04e; script 893d005; outputs 0a7c7ab; result document 6901214; all controls PASS). For each background the test asks whether a single M exists with g_φ, \|g_ψ\| ≤ g_max and all ρ_j ≤ ρ_max, i.e. whether c_max·G ≤ g_max·ρ_max with c_max = max(√p, \|κ\|/√p) and G = max(\|q\|,\|w\|,\|r\|,\|v\|).

   | Set | Common-M PASS / FAIL / BORDERLINE (all 41) | Without negative controls (29) |
   |---|---|---|
   | A | 0 / 41 / 0 | 0 / 29 / 0 |
   | B | 0 / 41 / 0 | 0 / 29 / 0 |
   | C | 21 / 20 / 0 | 15 / 14 / 0 |

   The B3 anchor fails C. Its admissible M interval [14.12, 6.205] is empty. As a diagnostic example only, its minimum common-M coupling at ρ ≤ 0.3 is g_min ≈ 8.07, with α_min ≈ 5.18.
5. **Population** (10 000 per cell, seed 172026). The principal-symbol front velocity is c_front > 1 in every sample of every cell. In the set-A region (√X ≤ 0.3) this holds in 20 000/20 000 samples, with an empirical minimum excess ≈ 0.5√X. The excess shrinks only together with the EFT effect (c − 1 ∝ gρ). This result is historical and unchanged.
6. **Recorded limitations.**
   - Λ_* ~ 1/√p is a dimensional estimate, not a computed bound.
   - The backgrounds are pointwise.
   - The setting is 1+1D.
   - c_front > 1 is not proof of acausality under preferred-foliation causality (decision of 2026-11-12).
   - The common-M result concerns compatibility of the stated constraints on these 41 backgrounds only. It does not show that the Proca construction is physically impossible.

## 2. Q1 — interpretation of the sets A/B/C and of the counts

Explanation and table only, with no new computation, as approved in v0.3.

| Set | Definition | Historical count (necessary condition) | Common-M count | Scope |
|---|---|---|---|---|
| A | g ≤ 1, ρ ≤ 0.3 | 0/41 | 0/41 | the 41 E3 backgrounds |
| B | g ≤ 1, ρ ≤ 0.1 | 0/41 | 0/41 | the 41 E3 backgrounds |
| C | g ≤ √(4π), ρ ≤ 0.3 | 41/41 (29/29) | 21/41 (15/29) | the 41 E3 backgrounds |

Reading:
- A and B fail even the necessary condition on every background. The failure is robust.
- For C, the historical count establishes only the necessary condition. Membership with a common M holds for about half of the backgrounds.
- None of these is a property of the model as a whole. The backgrounds are strong-X by construction. The thresholds are operational (see Q2), and α = 1 is at the border of perturbativity, not comfortably inside it.

**Proposed status: Q1 ANSWERED** (subject to review).

## 3. Q2 — validity of the M⁻² approximation: ε as a function of (ρ, k/M)

### 3.1 Quantities (not to be mixed)

| Symbol | Meaning | Type |
|---|---|---|
| ρ_j | background gradient scale, per component (§0) | input (from the background) |
| θ | direction/structure of the background gradient (w : q, v : r, φ/ψ mixing) | input (from the background) |
| k/M | wavenumber of the FLUCTUATION | input (grid) |
| ω/M | fluctuation frequency, from the dispersion relation | derived |
| ν = 1/(Mℓ) | scale of background variation; ν = 0 for constant-gradient backgrounds (C4) | input; inhomogeneous case only |
| ε_EFT, ε_tot, ε_abs, ε_op | approximation errors (§3.3) | derived |

### 3.2 Hypothesis H-Q2 (to be derived; not a result)

Expand in fluctuations, J = J̄ + J₁ + J₂ (orders 0, 1, 2). For ν = 0, at tree level, the exact quadratic part is expected to be ½·J₁ K(ω,k) J₁ + J̄·J₂/M². This would mean that all corrections beyond M⁻² act only on the J₁² part, with the factor

  M²K − 1 = (ω²/M²) / (1 + (k² − ω²)/M²)   (an exact algebraic identity from K).

The approximation M²K − 1 ≃ (k/M)² is **not** a general equality. It requires an on-shell condition, (ω² − k²)/M² ≪ 1 with ω ≈ k, and a suitable fluctuation regime. ρ and θ then enter through ω² − k², which is the B4 effect, at subleading order.

K points towards this structure, so the choice of error metric is not blind. H-Q2 is a preliminary hypothesis. Its derivation and check are separate from the choice of metric.

### 3.3 Error metrics (decided)

- **ε_EFT (primary):** the error of the M⁻² approximation to the dispersion correction (ω² − k²), i.e. to the B4 effect, relative to the exact K.
- **ε_tot (secondary):** the same error relative to the full ω².
- **ε_abs (secondary):** the absolute error.
- **ε_op (diagnostic only):** the norm of v0.3 §3 Q2 item 2c. It does not replace the primary definition.

### 3.4 Regime, reference, grids

- Tree level, classical. Only A is integrated out, using the A equation and the Gauss constraint.
- Fluctuations are linearised about the background.
- The dispersion relation is an on-shell observable. Q2 uses L_eff or the full K as they stand, with no reduction by the light-field equations and no use of the T7a off-shell basis. Option A and option B are not mixed.
- **Reference:** the exact K for ν = 0 only. For ν > 0 the reference is not carried over; the inhomogeneous background is a separate case, with the additional term ∂_tJ̄·∂_tJ₂/M⁴.
- **Grids and thresholds: not fixed in this document.** First the following must be established: the exact formula for ε_EFT; the range in which the fluctuation dispersion is physically admissible; and the relation to the C4 and E3 backgrounds. k/M ∈ 10⁻³…10⁰ is only a proposal for a preliminary grid. τ and ε_max are not proposed. No threshold is tuned after a result is seen.
- **Possible outcomes** (from v0.3):
  - (i) ε is a function of ρ alone over the admissible k/M range;
  - (ii) ε also depends on k/M, so a ρ threshold is not a standalone criterion;
  - (iii) no unique answer, so the threshold 0.3 stays OPERATIONAL.

**Next step for Q2:** a separate specification of the analytic derivation (H-Q2, the formula for ε_EFT), before any numbers.

## 4. Q3 — unitarity in 1+1D (option (a) only)

**Notation:** s_cm = 4E_L E_R for massless left- and right-moving particles.

| # | Condition | Content (decided) | Open |
|---|---|---|---|
| 1 | Normalisation | L = −½(∂f)² − ½(∂g)² in (−,+); [f] = [g] = 0 and [A] = 0 in 1+1D; g_φ, g_ψ dimensionless; p, s, κ of dimension M⁻² | confirmation of dimensions |
| 2 | Operator and coefficient | **primary:** the contact 4-point vertex from J²/(2M²), written as it follows from the Proca matching (g_φ = M√p, g_ψ = −Mκ/√p, so g_φg_ψ/M² = −κ and g_ψ²/M² = κ²/p): **L_contact = (p/2)(f_tf_x)² − κ(f_tf_x)(g_tg_x) + (κ²/2p)(g_tg_x)²**. Only under the rank-one condition is κ²/p = s, so that the last coefficient becomes s/2. The M⁻⁴ terms (R_S from T7a step 4) are reported separately | — |
| 3 | Channels | full coupled-channel matrix (ff, gg, fg, including ff↔gg via κ), then its eigenvalues | — |
| 4 | Kinematics | **primary: vacuum (φ̄ = 0).** Massless external lines. The nontrivial two-particle collision is between a left-moving and a right-moving particle, with s_cm = 4E_L E_R > 0. States moving in the same direction (LL or RR) have s_cm = 0 and do not describe the same collision. Scattering about the B4 background is a **separate sub-case** (cubic vertices J₁J₂/M², broken Lorentz invariance) | — |
| 5 | Amplitude conventions | to be derived: normalisation of two-particle states, 1+1D phase space, the S ↔ 𝒯 relation | derivation |
| 6 | Unitarity bound | derived from S†S = 1 in the two-particle sector: the dimensionless a(s_cm), the exact bound, and the rules for identical particles and coupled channels. No constant is carried over from 3+1D or from the standard partial-wave expansion | derivation |
| 7 | What is computed | primary quantity (ii): violation of the condition of item 6 → Λ_unit. Quantities (i) and (iii) are reported separately, without identification | — |

**Full tree-level Proca.** Computed separately, as a check of the tree-level UV dynamics. It is not identified with the EFT bound, and it does not mask loss of validity of the EFT expansion. A1·J is not a covariant coupling A_μJ^μ, so the high-energy behaviour is derived from this action, not assumed.

**Comparison with Λ_*.** Λ_* ~ 1/√p remains a diagnostic scale, and Λ_unit/Λ_* is a diagnostic, not a condition. Hypothesis H-Q3: with κ² = ps the channel matrix may have one nonzero eigenvalue that depends on a combination of p, s and κ. If this is confirmed, the comparison is reformulated, without a forced return to 1/√p.

**Further limitation.** Massless scalars in 2D have IR problems (Coleman). This is no obstacle for derivative couplings at tree level, but it is recorded.

**Next step for Q3:** a separate specification of the derivation (items 5 and 6), before any numbers.

## 5. Joint Q2/Q3 check

- **Common variable.** Q2 concerns fluctuations about a background, Q3 concerns vacuum scattering. A common variable is a device for comparing scales, not an identification of the two set-ups. The relation k = √s_cm/2 is admissible only for symmetric kinematics with a fixed frame and a stated identification.
- **Which scale binds first.** The dimensional estimate (Λ_unit ~ M/g versus k ~ M) is motivation only. The conclusion must use the derived formulas.
- **Uniqueness.** Λ_unit is unique only when the following are all fixed: the eigenvalue (item 3), the sector and constant (item 6), the vacuum (item 4) and the EFT contact (item 2).

## 6. Q4 — unchecked assumptions (approved in v0.3; unchanged)

Three levels are kept apart:
- (L) local validity about a given background;
- (D) dynamical attainment of that background;
- (S) robustness for time-dependent and spatially inhomogeneous backgrounds.

C4 covers only a limited part of (D) and (S), and only for the three registered backgrounds. Nothing of this is generalised. Every claim is marked "unchecked" or "partially checked in C4".

## 7. Possible end states (document level; unchanged from v0.3)

- **B4.12-CLOSED-NEGATIVE:** Proca is not a perturbative model for B4 in the controlled region, and it remains a negative control (as already accepted on 2026-11-10).
- **B4.12-CLOSED-POSITIVE:** requires BOTH (a) a computed unitarity that is acceptable in some regime, AND (b) that this regime satisfies the pre-defined validity criterion of the approximation (from Q2). Finding a regime with acceptable unitarity alone is NOT sufficient.
- **B4.12-OPEN:** not resolved, for lack of justification or computation.

The choice of state is the reviewer's, after review. The executor does not prejudge it.

**Current state: OPEN.**
- Q1 is answered (subject to review).
- Q2 and Q3 have fixed definitions but no derivation yet.
- Q4 is unchanged.

## 8. Rules

- No new code, simulations or changes to the criteria of C0–C4 until this document is registered.
- Every numerical or analytic step (Q2, Q3) requires a separate specification before it is executed.
- Conclusions are limited to the realisation considered (Proca, 1+1D).
- Registration of T7a still requires a registered decision on B4.12 (T7 protocol draft, §7).
- No Zenodo.

## 9. Requested from the reviewer

1. Approval of v0.4 as the B4.12 proposal, for registration.
2. Confirmation that Q1 is ANSWERED.
3. The order of the next two specifications: Q2 (analytic derivation) and Q3 (1+1D unitarity derivation).
