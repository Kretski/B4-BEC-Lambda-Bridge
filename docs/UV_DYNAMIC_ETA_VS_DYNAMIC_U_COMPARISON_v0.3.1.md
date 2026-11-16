# Dynamic η versus dynamic u^μ (Einstein-aether type) — comparison of the actions (v0.3.1)

**Changes from v0.3 (final review of 2026-11-16, one correction).** Only §2 is changed. The restricted variation V1 is tangential and does not determine λ. Only V2 gives the normal projection that determines λ. No other change.

**Changes from v0.2 (review of 2026-11-16, REVISION REQUIRED).**
1. §2: the dependence n(u) off the constraint surface is specified. Two admissible procedures are stated, and one is to be fixed in the script specification.
2. §5: c_η is stated as the characteristic speed of the isolated principal kinetic block of δη only, not of the full coupled system.
3. §10: the alternative of a separate candidate is removed, and the closing sentence on the UV order is aligned with the reviewer's decision.

**Changes from v0.1 (review of 2026-11-16, REVISION REQUIRED).** The parametrisation in §6, §9 and §10 is corrected. The general Einstein-aether case has **two** independent kinetic combinations (c₁₄, c₁₂₃). It is a one-parameter extension of dynamic η only after fixing c₁₄ = f, because c₁₄ cannot be absorbed into η: η enters u(η) and n(η) nonlinearly. §2 now requires the explicit check of the variation before and after the parametrisation, including the projections on u and n. The reduction c₁₄ = c₁₂₃ = f is unchanged.

**Status: PROPOSAL FOR REVIEW.** Comparison of two actions, required by the review of 2026-11-16. It contains hand derivations, each step stated. It contains no code, no run and no numerical result. The derivations are to be verified symbolically in the script stage, before any result relies on them. No registered criterion is changed. The Einstein-aether action is a comparison candidate, not a physical claim (reviewer's decision, 2026-11-16). B4.12 remains OPEN. T7a remains blocked. T7b is not started.

## 0. Setting and conventions

- Flat 1+1D Minkowski space, signature (+,−), as in `UV_DYNAMIC_ETA_SPEC_v0.4.1.md` (B1). No gravity: the metric is not dynamical. A coupling of u^μ to curvature is out of scope.
- A unit timelike vector in 1+1D is parametrised by its rapidity:

      u = (cosh η, sinh η),   n = (sinh η, cosh η),   u·u = 1,   n·n = −1,   u·n = 0,   du/dη = n,   dn/dη = u.

- Completeness: g^{αβ} = u^α u^β − n^α n^β.
- Directional derivatives: D_u = u^α ∂_α, D_n = n^α ∂_α, so that (∂η)² = (D_u η)² − (D_n η)².
- The dot products of the matter coupling follow `UV_DYNAMIC_ETA_SPEC_v0.4.1.md` B1 and the U1–U4 script.

## 1. The two actions

**Dynamic η** (`UV_DYNAMIC_ETA_SPEC_v0.4.1.md`, B1):

    L_η = L_matter(φ, ψ, σ; u(η), n(η)) + ½ f (∂η)²,   f > 0,

where L_matter = ½(∂φ)² + ½(∂ψ)² + ½(∂σ)² − ½M²σ² + σ Σ_i g_i (u·∂φ_i)(n·∂φ_i) − (Λ/4)(n·∂φ)⁴.

**Dynamic u^μ, Einstein-aether type** (reviewer's choice):

    L_u = L_matter(φ, ψ, σ; u, n) − ½ K^{αβ}_{μν} ∂_α u^μ ∂_β u^ν + λ (u·u − 1),
    K^{αβ}_{μν} = c₁ g^{αβ} g_{μν} + c₂ δ^α_μ δ^β_ν + c₃ δ^α_ν δ^β_μ + c₄ u^α u^β g_{μν}.

The overall sign −½K follows the reviewer's formula for signature (+,−). Sign conventions for Einstein-aether differ in the literature. With this one, c₁₄ > 0 gives a positive time-kinetic term (§5). In flat space ∇ = ∂. n is the unit spacelike vector orthogonal to u with a fixed orientation. In 1+1D it is determined by u.

The matter sector L_matter is **identical** in both. The comparison is about the dynamics of the frame field only.

## 2. Constraint and degrees of freedom

- u^μ has two components, and the constraint u·u = 1 (imposed by λ) removes one. With the future-directed branch, the general solution of the constraint is u = (cosh η, sinh η) with real η.
- **Equivalence of the variational problems: to be verified symbolically in the script stage, not accepted on the holonomic-constraint argument.** Three derivations are compared:
  - the variation of L_u with respect to u^μ and λ;
  - the variation with respect to η after the substitution u = (cosh η, sinh η);
  - the projections of the vector equation on u^μ and on n^μ.

  The expected outcome is that the u projection determines λ and the n projection is the η equation. The holonomic argument alone is not sufficient here, because the matter sector depends also on n^μ, whose orientation is derived from u^μ. The dependence of n on u must be included in the variation with respect to u.
- **How n(u) is defined off the constraint surface.** On the set of unit timelike vectors, n is fixed uniquely by u: n^μ = (u¹, u⁰), with the orientation of §0. For an arbitrary variation of both components of u, the extension of n(u) off that set must be stated, otherwise the variation with respect to u^μ is not defined. The script specification fixes one of the following two procedures and records which:
  - **(V1) Restricted variation.** Vary only along the constraint surface: δu^μ = n^μ δη and δn^μ = u^μ δη. The multiplier λ then drops out, and the result is the η equation (see below for what it is compared with).
  - **(V2) Explicit off-shell extension.** Define n^μ(u) off the surface, e.g. n^μ(u) = (u¹, u⁰)/√(u·u), so that n is unit for any timelike u and orthogonal to u. Then vary with respect to u^μ and λ freely.

  The two procedures do not have the same capabilities:
  - **V1** is a tangential variation along the surface u·u = 1, so λ drops out. It gives the η equation but **not** the normal equation. V1 by itself does not determine λ. Its η equation is compared with the tangential projection (along n) of the vector equation obtained in V2.
  - **V2** gives the full vector equation. Both projections are checked: the tangential projection (along n) must coincide with the η equation, and the normal projection (along u) determines λ.

  A check of components without a stated n(u) is not sufficient. If both procedures are used, any difference between them must be reported, not resolved by choosing the convenient result.
- **Degrees of freedom.** One field component (η) in both constructions. Whether it propagates is decided by the signs of its kinetic coefficients (§3).

## 3. Reduction of the Einstein-aether kinetic term in 1+1D

With u = u(η): ∂_α u^μ = (∂_α η) n^μ. Term by term:

| Term | Definition | Value in 1+1D |
|---|---|---|
| c₁ | g^{αβ} g_{μν} ∂_αu^μ ∂_βu^ν | (∂η)² (n·n) = −(∂η)² = −(D_uη)² + (D_nη)² |
| c₂ | (∂_μ u^μ)² | (n^μ ∂_μ η)² = (D_nη)² |
| c₃ | ∂_α u^μ ∂_μ u^α | (∂_α η n^μ)(∂_μ η n^α) = (D_nη)² |
| c₄ | u^α u^β g_{μν} ∂_αu^μ ∂_βu^ν | (D_uη)² (n·n) = −(D_uη)² |

Inserting into −½K∂u∂u:

    L_u,kin = ½ c₁₄ (D_u η)² − ½ c₁₂₃ (D_n η)²,   c₁₄ ≡ c₁ + c₄,   c₁₂₃ ≡ c₁ + c₂ + c₃.

**Answers to the reviewer's items 1 and 2:**
1. **Independent combinations.** Only c₁₄ and c₁₂₃ appear. Of the four parameters, two combinations remain in flat 1+1D.
2. **Equivalent and boundary terms.**
   - The c₂ and c₃ terms are equal **pointwise** in 1+1D for a unit field: both equal (D_nη)².
   - In flat space of any dimension, c₃-term − c₂-term = ∂_α(u^μ ∂_μ u^α − u^α ∂_μ u^μ), a total derivative.
   - The c₁ term splits between the two surviving combinations.
   - No twist (vorticity) term exists in 1+1D.

   These statements are to be verified symbolically.

## 4. Field equations and backgrounds

**η equation.**
- **Dynamic η:** f □η = ∂L_matter/∂η.
- **Einstein-aether:** with L_u,kin, the kinetic coefficients depend on η through u and n. The η equation therefore contains terms quadratic in ∂η in addition to the second-derivative terms.

**Constant backgrounds.** For constant η̄ all kinetic terms and their η derivatives vanish in both constructions. The background equations are therefore **identical**: the C0 equation of `UV_DYNAMIC_ETA_SPEC_v0.4.1.md`,

    d/dη [ J(η)²/(2M²) − (Λ/4) N(η)⁴ ] = 0 at η = η̄,   σ̄ = J(η̄)/M².

The set of stationary backgrounds, their counts and the domain of definition (T6 item 4) are the same for both actions. In Einstein-aether, λ on the background is fixed by the u component of the u equation. It does not affect η̄.

## 5. Quadratic η sector and characteristic speeds

Expand about a constant background, η = η̄ + δη, with ū, n̄ the background frame.

**Kinetic part of δη:**

| Construction | Kinetic part of δη |
|---|---|
| dynamic η | ½ f [(D_ū δη)² − (D_n̄ δη)²] |
| Einstein-aether | ½ c₁₄ (D_ū δη)² − ½ c₁₂₃ (D_n̄ δη)² |

The η-dependence of u and n in the aether kinetic coefficients contributes only at cubic order, because the kinetic terms are quadratic in ∂η.

**Characteristic speed of the isolated δη block.** c_η, below, is the characteristic speed of the isolated principal kinetic block of δη, linearised about a constant background. The characteristic speeds of the full coupled system (δφ, δψ, δσ, δη) are determined by the principal symbol of the whole system. They are **not** assumed to equal c_η, or 1, because δη mixes with the other fields (see "Mixing" below).

In the rest frame of ū (D_ū → ∂_t′, D_n̄ → ∂_x′), from the principal part of the δη block alone:

| Construction | Speed |
|---|---|
| dynamic η | c_η² = 1, the B4 unit cone |
| Einstein-aether | c_η² = c₁₂₃ / c₁₄ |

**Signs.**
- T5-type positivity of the δη kinetic energy and gradient energy: c₁₄ > 0 and c₁₂₃ > 0.
- Hyperbolicity of the δη principal part (T3-fund C2): c₁₂₃/c₁₄ > 0.
- Under the registered preferred-foliation convention, c_η > 1 is not by itself a failure. T3-fund tests well-posedness, not a speed bound.

**Mixing.** The mixing of δη with δφ, δψ and δσ comes from L_matter and is identical in the two constructions. The full fluctuation system (δφ, δψ, δσ, δη) of `UV_DYNAMIC_ETA_SPEC_v0.4.1.md` E5 applies to both, with the δη kinetic block replaced as above.

## 6. Reduction to the dynamic-η parametrisation

**Answer to the reviewer's item 5.**
- The Einstein-aether construction is **exactly** the dynamic-η construction, at the level of the action, if and only if

      c₁₄ = c₁₂₃ = f.

  Then L_u,kin = ½ f [(D_uη)² − (D_nη)²] = ½ f (∂η)².
- **General case.** The Einstein-aether kinetic sector has **two** independent combinations: c₁₄ > 0 (time normalisation of the frame mode) and c₁₂₃ > 0 (space normalisation). c₁₄ **cannot** be absorbed by a rescaling of η, because η enters the matter sector nonlinearly through u(η) = (cosh η, sinh η) and n(η). A rescaled field would change those couplings. Relative to dynamic η, which has one kinetic parameter f, the general aether case therefore has **one more kinetic parameter**, not one more ratio.
- **Fixed time normalisation.** Only after identifying c₁₄ = f is the aether case a one-parameter extension of dynamic η:

      c₁₄ = f,   c₁₂₃ = f c_η²,

  with c_η² = c₁₂₃/c₁₄ the speed of the frame mode relative to the B4 unit cone.
- Both constructions have one degree of freedom and the same matter sector.
- **Consequence for the candidate list.** In flat 1+1D, dynamic u^μ of Einstein-aether type is **not an independent UV candidate**. It is a generalisation of dynamic η with different time and space kinetic normalisations (c₁₄, c₁₂₃). Dynamic η is its special case c₁₄ = c₁₂₃ = f. This holds in flat 1+1D only:
  - in higher dimensions the twist term and the separate c₂ and c₃ structures do not reduce in this way;
  - with a dynamical metric, the aether couples to gravity.

## 7. Symmetries

- **Shared.** Both actions are invariant under global Lorentz boosts acting on the coordinates together with u → Λu, i.e. η → η + β, since u is a dynamical field in both. A background with constant η̄ ≠ 0, or with nonzero gradients, breaks the boost invariance **spontaneously**, not explicitly. Spacetime translations are unbroken by constant-gradient backgrounds, up to the shift symmetries of φ and ψ.
- **Dynamic η only.** The kinetic term ½f(∂η)² is also invariant under η → η + const **alone**, without a coordinate transformation. The matter coupling breaks this through u(η) and n(η).
- **Einstein-aether.** For c₁₄ ≠ c₁₂₃ the kinetic term depends on η through u and n. It is not invariant under the shift alone. Only the combined boost symmetry remains.

## 8. Elimination and EFT operator structure

**σ.** Its elimination is identical in both, by `UV_DYNAMIC_ETA_SPEC_v0.4.1.md` F1:
- the exact Gaussian integral is nonlocal;
- the leading local term is J_tot²/(2M²);
- the remainders are O(M⁻⁴);
- the rank-one structure s = κ²/p holds by F3.

**η.** The condition of F2 applies with the time kinetic coefficient f replaced by c₁₄:

    m²_η,eff = −(1/c₁₄) d²/dη² [ J²/(2M²) − (Λ/4) N⁴ ] at η = η̄   (Einstein-aether),

with f in place of c₁₄ for dynamic η. The limit in which η becomes auxiliary is c₁₄ → 0 together with c₁₂₃ → 0 for aether, and f → 0 for dynamic η. Whether these limits are regular is part of F2(iii) and is not presumed.

**EFT structure at tree level.** Identical matter operators from σ exchange in both. Any difference in the low-energy theory comes from the η sector:
- if η is light, it remains a light field with speed c_η;
- if it is eliminable, the elimination rule of F2 applies with c₁₄ in place of f.

## 9. Summary of the comparison

| Item | Dynamic η | Dynamic u^μ (Einstein-aether, flat 1+1D) |
|---|---|---|
| Kinetic term | ½f(∂η)² | ½c₁₄(D_uη)² − ½c₁₂₃(D_nη)² |
| Constraint | none (u unit by parametrisation) | u·u = 1 via λ; solved by u = (cosh η, sinh η) |
| Parameters of the frame sector | f | c₁₄, c₁₂₃: two independent combinations (c₁…c₄ reduce to two; c₂ and c₃ coincide pointwise). With c₁₄ = f fixed: one extra parameter, c_η² = c₁₂₃/f |
| Degrees of freedom | 1 | 1 |
| Symmetries | boosts with η → η + β; η-shift of the kinetic term alone | boosts with η → η + β only (for c₁₄ ≠ c₁₂₃) |
| Background equations (constant η̄) | C0 equation | identical |
| Characteristic speed of the isolated δη block (not of the full coupled system) | 1 | √(c₁₂₃/c₁₄) |
| Conditions (positivity, hyperbolicity) | f > 0 | c₁₄ > 0, c₁₂₃ > 0 |
| Elimination of σ | F1 | identical |
| Elimination of η | F2 with f | F2 with c₁₄ |
| EFT operator structure | rank-one B4 from σ exchange, plus η sector | identical matter part; η sector differs by c_η |
| Reduction | — | equals dynamic η iff c₁₄ = c₁₂₃ = f |

## 10. Proposed decision (for the reviewer)

1. In flat 1+1D, treat dynamic u^μ of Einstein-aether type **not as a separate fifth candidate**, but as a generalisation of dynamic η with different time and space kinetic normalisations (c₁₄, c₁₂₃). Dynamic η is the special case c₁₄ = c₁₂₃ = f. This is the reviewer's decision of 2026-11-16.
2. If it is tested within the dynamic-η programme, the two parameters are treated as such: either the general pair (c₁₄, c₁₂₃), or c₁₄ = f with c₁₂₃ = f c_η². The choice and the values must be fixed in a separate registered document, together with the criteria, before any run. No value of c₁₄ or c₁₂₃ is proposed here.

**Consequence for the UV order (reviewer's decision, 2026-11-16).** Dynamic u^μ of Einstein-aether type in flat 1+1D is a comparison parametrisation and generalisation within the dynamic-η programme, not a separate fifth candidate. The UV order of separate candidates is therefore:
- dynamic η, including the Einstein-aether generalisation (c₁₄, c₁₂₃) as a variant within that programme;
- HeavyScalar, one light field;
- HeavyScalar, two light fields;
- followed by the comparative analysis and the elimination of the mediators.

A dynamic u^μ with a different action (e.g. a free vector with a potential, or another construction) would be a new candidate and would require its own specification.

## Out of scope

- Gravity: a dynamical metric, curvature couplings of u, and the aether–graviton sector.
- Dimensions above 1+1, where the twist term and the separate c₂ and c₃ structures matter.
- Other dynamic-u^μ actions (a free vector with a potential, or other constructions). They are not compared here.
- Any computation, code or numerical evaluation.
- T7b, T7a and B4.12.

**B4.12: OPEN. T7a: blocked. T7b: not started.**
