# Dynamic u^μ as a UV route to T7b — admissibility note v0.7

**Status:** proposal for review. Analytic only: hand derivations, no code, no run, no numerical result. Each derived formula is to be checked symbolically before registration, in the way the registered structural note was checked (§7).

**Base:** 73347c3.

**Changes from v0.6 (review of 2026-11-19 11:07).** Exactly three changes, all in §4.3a. The classification, condition (R), G1-restricted, the G1/G2 boundary and §5–§7 are unchanged.

1. **Step 6.1.** Condition 6.1′ is added: x − X* < 2X*·m/M, with m = min_{[X*,x]}|V″| and M = max_{[X*,x]}|V″|. It gives D < 0 on the A₁ = 0 part of the outer path.
2. **Step 6.3.** The exclusion of individual amplitudes is removed, with the reason: by 6.1, D < 0 for every A² ≥ X*. The two places in step 7 that referred to the excluded values now cite 6.1 instead.
3. **Step 8.** The O(1) remainder is justified by reference to `DYNAMIC_U_G1_LAYER_LEMMA_v0.3.md` (4ae480c5…), with its scope (A ≥ 4, A > √x) and its premises.

**Changes from v0.5 (reviews of 2026-11-19 06:42 and 10:27).** Three text corrections only. The classification, condition (R), the G1/G2 boundary and the scope of the results are unchanged.

1. **Scope of the unboundedness statement.** It now requires all conditions of §4.3a, step 6.1, not only V′(X*) > 0. This is changed in §4.3a, step 9, and in the G1-standard summary of §4.3, which carried the same wording.
2. **Vacuum normalisation.** The energy is measured relative to the vacuum density V(x) (§4.3a, step 8).
3. **Energy formula.** The full expression E_rel = 2ℓ[V(X*) − V(x)] + A²vℓ(−4 + (4/3)vℓ²) is used in step 8 and in §7, item 5.

**Symbolic checks.** All §7 identities and the two scope checks were verified by `dynamic_u_admissibility_check.py` (SHA-256 152a7e39…): one authorised run from a clean clone of 73347c3, 27 checks, 0 failed (stdout 61f2d0e8…). The checks confirm symbolic identities only. They do not prove (R), the existence of the layer profiles, or the remainder estimates.

**Changes from v0.4 (review of 2026-11-19 06:40).** Two corrections in §4.3a. Nothing else changes.

1. **(R) is quantified over amplitudes that grow without bound.** It now requires a sequence A_n → ∞ for which the two plateau points lie in the same connected component of {D < 0}. Connectedness for a single A is not enough.
2. **The layer width is ε(A) ≤ 1/(A·C(A)).** With this choice the layers contribute O(1/A) and the cross term π_plateau·δπ in ½π² is O(1), so the stated E(A) = … + O(1) holds. With the v0.4 choice that cross term was only O(A).

**Changes from v0.3 (review of 2026-11-19 06:34).** Only §4.3a, steps 6–9, the G1-standard statement and §7 item 5 change. Three things are made explicit:

1. the order of choices and limits;
2. an exact zero of ∫2A₀V′ dx, from an odd profile;
3. the regularity of the second-class constraints.

The third point **weakens** the v0.3 claim:

- Unboundedness below holds on the full constraint set (both constraints satisfied pointwise) for every smooth V with V′(X*) > 0 at some timelike X*.
- On the regular constraint surface, where the second-class matrix is invertible everywhere, it is proven only under a condition (R) on V away from the vacuum.

On the line A₁ = 0, (R) can never be met: there is always a degenerate point (Rolle). Whether (R) holds is open; if it fails, the question belongs to G1-restricted. The classification of U2 is stated accordingly.

**Changes from v0.2 (review of 2026-11-19 06:30).** Only §4.3 (G1-standard), §7 items 5–6 and the U2 row of the brief result change. G1-standard is now derived directly in 1+1D, with our variables (§4.3a). The 3+1D result is no longer carried over; it serves only as agreement. The configurations used for the no-go are shown to be admissible, i.e. on the constraint surface. They are an explicit finite-energy family, replacing the phrase "arbitrarily small perturbations".

**Changes from v0.1 (review of 2026-11-19).** Only §4.3, §5, §6, §7 and the references change. §1–§4.2 and §4.4 are unchanged.

1. **G1 is split.**
   - **G1-standard:** the Hamiltonian is bounded below on the standard physical phase space. For a smooth V this is now a no-go for U2, from the necessary condition V′(X) ≤ 0 of Bailey, Potting, Walter-Cardona.
   - **G1-restricted:** a restricted phase space that keeps the δA₀ channel. This remains open; the known restriction X = x removes the channel.
2. **The three conditions are kept apart:** boundedness of V(X), boundedness of the Hamiltonian after the constraints, and linear stability of the vacuum.
3. **Citations corrected.**
   - arXiv:0812.1049 is now attributed to Carroll, Dulaney, Gresham, Tam, with its actual scope (fixed-norm aether).
   - The instability of the quadratic bumblebee potential is now cited to Bailey et al.
   - The cubic example is described as a degenerate minimum on which the dynamics is Maxwell in a nonlinear gauge.

**Unchanged:** the registered constraints A, B and C (`docs/UV_STRUCTURAL_CONSTRAINTS_NOTE_v0.3.md`, c7b92ba), B4.12, T7a, Script B/C4, EXT-1.

**Question.** Can a dynamical frame field u^μ produce Λ > 0 at tree level, leading local order, with no bare Λ, while respecting A, B and C and the lesson of the closed Proca line (73347c3, NOT-VIABLE (d))?

**Result in brief:**

| Realisation | Class |
|---|---|
| U1 — unit-norm frame (Einstein-aether type) | **EXCLUDED** for T7b |
| U2 — vector condensate (Maxwell kinetic term, smooth potential, spontaneous timelike vacuum value) | **EXCLUDED on the full constraint set**, and on the regular constraint surface if (R) holds (G1-standard, §4.3a). **CONDITIONAL** only on an unidentified restricted phase space that keeps the δA₀ channel (G1-restricted) |

**Recommendation:** no full specification, no G2 and no numerical work. U2 stays open only through G1-restricted (§6).

## 0. Conventions

- Flat 1+1D, signature (+,−), as in the registered structural note.
- X₁ = u·∂φ, X₂ = n·∂φ, X₃ = u·∂ψ, X₄ = n·∂ψ.
- B4 contains −(Λ/4)X₂⁴; we need Λ > 0.
- "Tree level, leading local order": classical elimination of the mediator, gradient terms of the mediator dropped (as in C1 and C2).
- Labels: **Proven** (follows from the stated premises), **Conditional**, **Open**.

## 1. What is already registered

- **A (Proven):** the B4 operators need a preferred frame. Its origin is open item T6.
- **C1 (Proven):** exchange of healthy scalar modes gives Λ_eff ≤ 0.
- **C2 (Conditional):** the time component of a vector mediator gives the opposite sign. For Proca, Λ_eff = 2(γ_u² − γ_n²)/M².
- **The mechanism behind C2.** Λ > 0 at tree level comes only from a **non-dynamical (auxiliary) field component whose quadratic term has the sign of the Proca A₀ term**. Eliminating such a component gives −J₀²/(2M²) rather than the +J²/(2M²) of a healthy mode.
- **`UV_DYNAMIC_ETA_VS_DYNAMIC_U_COMPARISON_v0.3.1` §6.** In flat 1+1D, the unit-norm (Einstein-aether) frame is not an independent candidate. It is the dynamic-η construction with two kinetic constants (c₁₄, c₁₂₃). Dynamic η is the special case c₁₄ = c₁₂₃ = f.
- **Dynamic η, Script A (636d278).** Closed: rank one, as B1 predicts. Its Λ was inserted by hand.
- **Proca with a timelike current (73347c3).** NOT-VIABLE (d). The exchange term 4γ_u qB(c) puts a pair of characteristic speeds off the real axis on backgrounds that are strictly hyperbolic without it.

## 2. Origin of the frame: the two realisations

A dynamical u^μ in 1+1D has two field components. The frame comes from its timelike background value in either case. The two realisations differ in what happens to the component **along** the background.

- **U1, unit norm.** A multiplier imposes u·u = 1. The component along u is removed. One mode remains: the rapidity (frame) mode η, with u = (cosh η, sinh η).
- **U2, condensate.** L = −¼F_{μν}F^{μν} − V(X) + A^μJ_μ, with X = A_μA^μ and a timelike vacuum value Ā = (a, 0), X̄ = a² ≡ x. "Spontaneous" means V′(x) = 0, with no source needed to hold the vacuum value. Both components of the fluctuation are present:
  - δA₀ along the vacuum value. It is auxiliary, because F contains only ∂₁A₀.
  - δA₁ across it.

  The registered Proca model is not of this kind: its background Ā₀ = −J̄_u/M² is held by the source (V(X) = −½M²X, V′ ≠ 0).

## 3. U1 — unit-norm frame: EXCLUDED for T7b

**Statement U1-a (Proven, given §1).** The unit-norm constraint removes the only kind of component that can give Λ > 0 at tree level: an auxiliary component along u. The remaining mediator is the frame mode η. Its kinetic term is ½c₁₄(D_uη)² − ½c₁₂₃(D_nη)², and T5 positivity requires c₁₄ > 0 and c₁₂₃ > 0 (comparison note §5). So η is a healthy scalar mode.

**Statement U1-b (Proven).** The exchange of η falls under C1. The η mass² comes from the background (the C0 stationarity condition of the comparison note §4). Three cases:

- **m_η² > 0:** C1 applies. The η source contributes Λ_eff = −2h²/m_η² ≤ 0, where h is its X₂² coefficient.
- **m_η² < 0:** the background is tachyonically unstable (T5 fails).
- **m_η² = 0:** the exchange is non-local, so it gives no B4 contact operator.

In addition, any healthy massive scalar in the matter sector (σ in the dynamic-η action) contributes ≤ 0 by C1.

**Consequence (Proven).** In U1, Λ > 0 cannot be generated at tree level, leading local order. This is the dynamic-η result: there, Λ had to be inserted by hand. U1 cannot solve T7b. Whatever further role U1 has as the origin of the frame (T6) does not touch T7b.

## 4. U2 — vector condensate: CONDITIONAL

### 4.1 The Λ channel (Conditional: hand derivation, leading local order)

Expand about Ā = (a, 0), with V′(x) = 0:

X = x + 2aδA₀ + δA₀² − δA₁²,  so  −V ⊃ −½V″(x)(2a δA₀)² = −2a²V″(x) δA₀².

The vector couples directly to the current, A^μJ_μ ⊃ δA₀ J₀, with J₀ = u·J = γ_u X₂² for the T7b channel. Eliminating the auxiliary δA₀ (gradient term ∂₁δA₀ in F dropped):

δA₀ = J₀/(4a²V″),  L_eff = J₀²/(8a²V″),

so that

**Λ_eff = −γ_u² / (2a² V″(x)).**

**Cross-check against C2.** Proca is V = −½M²X; its A₀ term is +½M²A₀², which gives L_eff = −J₀²/(2M²) and Λ = 2γ_u²/M² (registered C4). In U2 the same sign needs the A₀ quadratic coefficient −2a²V″ to be positive.

**Necessary condition (Conditional).** Λ > 0 if and only if **V″(x) < 0**. Together with V′(x) = 0, this means the vacuum value sits at a **local maximum of V along X**.

Two further points:

- If the matter couples only to the **direction** u^μ = A^μ/√X, and not to A^μ itself, then δA₀ does not enter u at linear order. The channel is absent, and U2 reduces to U1 for T7b.
- A quadratic potential V = (λ/2)(X − x)² (λ > 0) gives Λ_eff < 0, the wrong sign.

### 4.2 The transverse channel (Conditional)

With V′(x) = 0 the δA₁ mass term vanishes. δA₁ is the massless Goldstone mode of the broken boosts, the U2 counterpart of η.

- A spatial current J₁ = n·J ≠ 0 would therefore produce a non-local term J₁ (1/∂²) J₁, not a B4 contact operator.
- A local B4 match then requires J₁ = 0. In that case p, s, κ cannot come from this field and need a separate mediator (B1, B1-V). That is a two-field construction, and the registered rank and sign constraints apply to it.

### 4.3 Stability: three separate conditions

These three conditions are not equivalent and are kept apart:

- **(i) Boundedness of V(X)** as a function. This is not sufficient.
- **(ii) Boundedness of the Hamiltonian after the constraints.** This is the condition that matters for the physical phase space.
- **(iii) Linear stability of the chosen timelike vacuum.**

**G1-standard (§4.3a, smooth V).**

- **On the full constraint set (Proven):** if there is a timelike X* satisfying **all** conditions of §4.3a, step 6.1, the Hamiltonian is unbounded below. The conditions are V′(X*) = v > 0, V″ < 0 on [X*, x], and 2v + 4X*V″(X*) < 0; V′(X*) > 0 alone is not what is proven. For U2, V′(x) = 0 and V″(x) < 0 give V′(X) = V″(x)(X − x) + O((X − x)²) > 0 for X slightly below x, and all three conditions hold for X* close enough to x. So such an X* exists. **On the full constraint set, U2 cannot combine V′(x) = 0, V″(x) < 0 and a Hamiltonian bounded below.**
- **On the regular constraint surface (D ≠ 0 everywhere): Proven under condition (R)** of §4.3a, step 7. (R) requires a sequence of amplitudes A_n → ∞ and is open. A crossing of A₀ = 0 along A₁ = 0 is never regular. Where D = 0 the rank of the constraint matrix changes, so the set D = 0 is not counted as part of the regular physical phase space.

This agrees with the necessary condition V′(X) ≤ 0 of Bailey, Potting, Walter-Cardona (arXiv:2610.21066, §II, Eqs. 13–14), derived there in 3+1D. Only the 1+1D derivation below is used as the argument.

### 4.3a Derivation of G1-standard in 1+1D (hand derivation; symbolic check in §7, item 5)

**Setting.** No source. Signature (+,−), A_μ = (A₀, A₁), X = A₀² − A₁², F₀₁ = ∂₀A₁ − ∂₁A₀, so

L = −¼F_{μν}F^{μν} − V(X) = ½F₀₁² − V(X).

**1. Momenta.**
- π ≡ ∂L/∂(∂₀A₁) = F₀₁.
- π⁰ ≡ ∂L/∂(∂₀A₀) = 0: the primary constraint.
- ∂₀A₁ = π + ∂₁A₀.

**2. Canonical Hamiltonian.** After an integration by parts with fields of compact support:

H = π∂₀A₁ − L = ½π² + π∂₁A₀ + V(X) ≃ ½π² − A₀∂₁π + V(X).

**3. Secondary constraint.** Preservation of π⁰ = 0 requires ∂₀π⁰ = −∂H/∂A₀ = 0:

Ω₂:  ∂₁π = 2A₀V′(X).

The pair (π⁰, Ω₂) is second class wherever ∂Ω₂/∂A₀ = 2V′ + 4A₀²V″ ≠ 0. On the constraint surface the canonical Hamiltonian is the conserved energy.

**4. Reduced Hamiltonian.** Substituting Ω₂ into the term −A₀∂₁π:

H ≈ ½π² + V(X) − 2A₀²V′(X).

This is the same structure as in 3+1D, without the magnetic term, which does not exist in 1+1D.

**5. What is different in 1+1D.** π has a single component and **no transverse part**. Ω₂ therefore fixes π completely, up to a constant: π(x) = ∫^x 2A₀V′(X) dx′. Finite energy needs π → 0 at both ends, i.e. ∫2A₀V′ dx = 0. So the kinetic term ½π² is not free and could in principle compensate the negative term. The 3+1D argument cannot simply be carried over, and the compensation has to be checked.

**6. Choices, in this order.**

1. Fix X* < x close enough to x that:
   - V′(X*) = v > 0;
   - V″ < 0 on [X*, x];
   - 2v + 4X*V″(X*) < 0;
   - **(6.1′)** x − X* < 2X*·m/M, with m = min_{[X*,x]}|V″| and M = max_{[X*,x]}|V″|.

   All four are possible because V′(x) = 0, V″(x) < 0 and V″ is continuous: as X* → x, v → 0, m/M → 1, and x − X* → 0. Condition 6.1′ is used for the A₁ = 0 part of the outer path (step 7). For X ∈ [X*, x], V′(X) ≤ M(x − X) and 4XV″(X) ≤ −4Xm, so D ≤ 2M(x − X*) − 4X*m < 0 (Layer Lemma v0.3, L2).
2. Fix ℓ > 0 with ℓ² < 3/v.
3. Take any amplitude A > √x. No individual values of A need to be excluded. Since V″(X*) < 0, for every A² ≥ X*

   2v + 4A²V″(X*) ≤ 2v + 4X*V″(X*) < 0  (by 6.1).

   So D < 0 on the plateaus and on the hyperbola part of the outer path, for every A; A > √x is stronger than A² ≥ X*. (v0.6 excluded at most two values of A here; that exclusion was redundant.) This does not concern the central layer, whose regularity is condition (R).
4. For each A, choose a layer width ε(A) ≤ 1/(A·C(A)) (step 8).
5. Only then let A → ∞. On the regular surface, take A along the sequence A_n of (R).

**7. The odd profile.** The boundary values are the two vacua, A₀(−∞) = −a and A₀(+∞) = +a, with A₁ → 0. This is a finite-energy sector with vacuum ends.

- **Region R, x ∈ [ε, ε + ℓ]:** A₀ = +A and A₁ = √(A² − X*), so X = X*.
- **Region L:** the mirror image, A₀(−x) = −A₀(x) and A₁(−x) = A₁(x).
- **Outer layer, of width ε:** joins region R to the vacuum (a, 0). Its path in the (A₀, A₁) plane runs along A₁ = 0 from (a, 0) to (√X*, 0), with X between X* and x, and then along the hyperbola X = X* out to (A, √(A² − X*)).
- **Central layer, of width 2ε:** joins L and R through A₀ = 0.

**Exact zero of the integral.** A₀ is odd and A₁ is even, so X and V′(X) are even and A₀V′(X) is odd. Hence ∫2A₀V′ dx = 0 **exactly**. π, defined by Ω₂ with π(−∞) = 0, is even and returns to 0 at +∞ exactly. Both constraints hold pointwise.

**Regularity of the second-class pair.** The matrix is invertible where D ≡ 2V′(X) + 4A₀²V″(X) ≠ 0:

- **On the plateaus:** D = 2v + 4A²V″(X*) < 0 for every A > √x, by step 6.1 (see step 6.3).
- **On the hyperbola part of the outer layer:** D = 2v + 4A₀²V″(X*) ≤ 2v + 4X*V″(X*) < 0 by step 6.1.
- **On the A₁ = 0 part of the outer layer:** D = 4√X · d/dX[√X V′(X)] = 2V′ + 4XV″. This is < 0 for X ∈ [X*, x] by condition 6.1′ (step 6.1). At X = x it equals 4xV″(x) < 0.
- **At the vacuum:** D = 4xV″(x) < 0.

So the whole outer layer, the plateaus and the vacuum ends are regular, with D < 0.

**The central layer is the exception.** It must cross A₀ = 0, where D = 2V′(X) with X = −A₁² ≤ 0. Along the line A₁ = 0, with g(X) ≡ √X V′(X), D = 4√X g′(X). Since g(0) = g(x) = 0, Rolle's theorem gives g′ = 0, and therefore **D = 0, somewhere in (0, x)**. A crossing along A₁ = 0 is therefore never regular. A regular crossing exists if and only if:

**(R)** there is a sequence A_n → ∞ (A_n > √x) such that for every n the two plateau points (A_n, √(A_n² − X*)) and (−A_n, √(A_n² − X*)) lie in the same connected component of {D < 0} in the (A₀, A₁) plane. The connecting path necessarily passes through the spacelike or null part (X ≤ 0), since A₀ changes sign.

Connectedness for one fixed A is not enough: the construction needs amplitudes without bound.

The Rolle argument shows only that a crossing **along A₁ = 0** meets D = 0. It does not show that every crossing is singular. Whether a detour through the spacelike region exists is exactly the content of (R).

(R) depends on V away from the vacuum and is not implied by the U2 assumptions.

**8. Energy, with the order of limits.**

**Normalisation.** The energy is measured relative to the vacuum. The vacuum (π = 0, A₀ = a, A₁ = 0) has reduced density V(x) − 2a²V′(x) = V(x). The relative density is H − V(x). Without this subtraction the integral over the infinite line diverges when V(x) ≠ 0. Equivalently, one may set V(x) = 0. The coefficient of A² does not depend on this choice.

- **Plateaus:** kinetic (4/3)A²v²ℓ³; potential, relative to the vacuum, 2ℓ[V(X*) − V(x)] − 4A²vℓ. On the plateaus π_plateau = O(A). Outside the configuration the relative density is 0.
- **Layers:** H contains no gradient of A₀ or A₁. For fixed A the layer paths are bounded, so with C(A) finite:
  - the layer energy is at most ε·C(A);
  - the change of π across a layer is δπ = O(ε·C(A)).
- **Choice ε(A) ≤ 1/(A·C(A)).** Then:
  - the layer energy is O(1/A);
  - the cross term ∫π_plateau δπ dx over the plateaus is O(A·ε·C(A)) = O(1);
  - the term ∫(δπ)²/2 dx is O(ε²C(A)²) = O(1/A²).

  (With only ε(A) ≤ 1/C(A), the cross term would be O(A). That is still o(A²) and would suffice for the proof, but not for the O(1) written below.)

**Justification of the O(1) remainder.** The bullets above are made precise in `DYNAMIC_U_G1_LAYER_LEMMA_v0.3.md` (SHA-256 4ae480c5…; reviewed and accepted on 2026-11-19 within its scope, not registered).

- **Premises:** V ∈ C²; X*, v, ℓ fixed as in steps 6.1 (including 6.1′) and 6.2; the odd profile of step 7, with π defined by Ω₂. π is C¹. The profile is C^∞ under the lemma's conditions on the path pieces and the reparametrisation.
- **Construction:** C(A) = max{C_out(A), C_c(A), ½(2Avℓ + 1)², 1}, where C_out and C_c are suprema of F = |V − V(x)| + 2A₀²|V′| + 2|A₀V′| along the outer and central paths, and ε(A) = 1/(A·C(A)).
- **Bound:** for **A ≥ 4 and A > √x**,

  |E_rel(A) − E_main(A)| ≤ 8/A + 16vℓ² + 16ℓ/A²,  with E_main(A) = 2ℓ[V(X*) − V(x)] + A²vℓ(−4 + (4/3)vℓ²).
- **Scope:** the lemma covers the full constraint set and the regularity of the outer layers and plateaus. It does **not** cover the regularity of the central layer, which is condition (R).

Hence

E_rel(A) = 2ℓ[V(X*) − V(x)] + A²vℓ(−4 + (4/3)vℓ²) + O(1),  and E_rel(A) → −∞ as A → ∞, because ℓ² < 3/v. The A² coefficient equals (4/3)vℓ(vℓ² − 3), which has the sign of vℓ² − 3.

**9. Conclusion and scope.**
- **On the full constraint set (Proven, smooth V).** The Hamiltonian, measured relative to the vacuum, is unbounded below whenever there is a timelike X* satisfying all conditions of step 6.1: V′(X*) = v > 0, V″ < 0 on [X*, x], and 2v + 4X*V″(X*) < 0. The construction uses all three, through the regularity of the outer layers and the plateaus. It is not proven for an arbitrary potential with V′(X*) > 0 at some point. For U2 such an X* exists just below x.
- **On the regular constraint surface (Proven under (R)).** If (R) holds, the family taken along A_n has D ≠ 0 everywhere and E(A_n) → −∞. Whether (R) holds is **open**. If it fails, every finite-energy configuration whose A₀ changes sign passes through the degenerate set. In that case the "regular phase space" is itself a restriction of the phase space, and the question moves to G1-restricted.
- **What is not claimed.** Nothing about (iii), linear stability of the vacuum, and nothing about nonlinear evolution. The family uses large field values. V′ ≤ 0 is necessary only.

**G1-restricted (Open).**

Bluhm, Gagne, Potting, Vrublevskis (arXiv:0802.4071) find, for certain bumblebee models with a timelike vector, restrictions of the initial data that give a ghost-free model with a positive Hamiltonian. The restricted phase space then matches electromagnetism in a nonlinear gauge.

That restriction is X = x. Linearised about Ā = (a, 0), δX = 2a δA₀ = 0, so **δA₀ = 0. This removes exactly the channel that generates Λ_eff in §4.1.** Bailey et al. reach the same structure:

- a field that starts at the minimum stays there (their §IV A);
- on that set the solutions are those of Maxwell theory in the nonlinear gauge B_μB^μ = x (their §V A);
- at linear order about the minimum the equations are those of unperturbed Maxwell theory (their §IV B).

Open question: is there an admissible restriction of the phase space that keeps a non-zero δA₀ exchange? None is known, and none is proposed here.

**Context, not used as the argument:**

- **Cubic potential.** The stable cubic potential V = −(λ/3)(X − x)³ of arXiv:2610.21066 has V′(x) = V″(x) = 0, so it gives no local X₂⁴ term from §4.1. Its stability rests on a degenerate minimum, on which the dynamics reduces to Maxwell in a nonlinear gauge. The radial channel is absent for that reason, not by accident.
- **Quadratic potential.** It is unstable (Bailey et al., §II–III) and gives the wrong sign here (§4.1).
- **Fixed-norm æther.** Carroll, Dulaney, Gresham, Tam (arXiv:0812.1049) treat fixed-norm æther fields. Generic kinetic terms give ghosts or tachyons. For the Maxwell and scalar kinetic terms the Hamiltonian is unbounded below. The timelike sigma model is stable when the norm is fixed by a constraint, but not when it is set by minimising a potential. This concerns U1-type models and the general picture. It is not the source of the U2 no-go.

### 4.4 The lesson of (d) (Open)

In Proca the failure came from the exchange term in the principal symbol. With J_n = 0 the registered B(c) reduces to M²γ_u q(1 − c²), so in Proca this term survives even without a spatial current:

Q_new = F_eff − c²d1d2 − 2Λq²(1 − c²).

**Open for U2:** whether the corresponding term arises when the transverse mode is massless and δA₀ is auxiliary with coefficient −2a²V″. The principal symbol comes from k → ∞, where the Proca propagator term k_μk_ν/M² produced the γ_u²/M² ∝ Λ contribution. The U2 propagator has a different mass structure (0 across, −4a²V″ along), so the term **must be derived, not carried over from Proca**. It is possible that U2 gives the same mechanism, since Λ is generated by the same auxiliary channel.

## 5. Classification

| Item | U1 | U2 |
|---|---|---|
| A: origin of the frame | dynamical, spontaneous | dynamical, spontaneous |
| C: Λ > 0 at tree level, local | **impossible** (U1-a, U1-b) | possible only if V″(x) < 0 and the coupling is to A^μ itself (§4.1) |
| B: p, s, κ | rank one at most, as dynamic η | not from this field (J₁ = 0 is forced); needs a second mediator (§4.2) |
| Stability | as dynamic η (Script A: 35/41 FAIL-UNSTABLE with Λ inserted by hand) | **G1-standard fails** for smooth V on the full constraint set, and on the regular surface under (R); G1-restricted is open (including the case that (R) fails), and the known restriction X = x removes the Λ channel (§4.3, §4.3a) |
| Lesson of (d) | not relevant: no Λ channel | **open** (§4.4) |
| **Class** | **EXCLUDED** for T7b | **EXCLUDED on the full constraint set**, and on the regular surface if (R) holds; **CONDITIONAL** only through G1-restricted (which includes the case that (R) fails) |

## 6. Recommendation

No full specification of dynamic u^μ, no G2 and no numerical work.

- **U1:** closed for T7b.
- **U2:** excluded for smooth V on the full constraint set, and on the regular constraint surface if (R) holds (G1-standard). If (R) fails, the regular surface counts as a restricted phase space under G1-restricted.
- **The only open path is G1-restricted:** an admissible restriction of the phase space, possibly the regular surface itself if (R) fails, that is stable and keeps a non-zero δA₀ exchange. Until one is shown, U2 is treated as excluded in practice. It is not recorded as NOT-VIABLE, because G1-restricted is open.
- **G2** (the principal symbol, §4.4) is not considered unless G1-restricted is answered positively.

**Consequence for T7b.** The tree-level mediator routes analysed so far give no admissible Λ > 0:

| Route | Result |
|---|---|
| Healthy scalars | C1 |
| Proca with a timelike current | NOT-VIABLE (d) |
| Unit-norm frame | U1 |
| Vector condensate on the full constraint set (and on the regular surface under (R)) | U2 |

What remains are the two routes the structural note already names, a justified contact term and quantum (loop) corrections, plus U2 through G1-restricted.

## 7. To be checked symbolically at registration

1. U2 expansion: the δA₀ coefficient −2a²V″(x) and the δA₁ mass term V′(x) = 0.
2. L_eff = J₀²/(8a²V″) and Λ_eff = −γ_u²/(2a²V″), including the cross-check with C4 for V = −½M²X.
3. The reduction of the registered B(c) to M²γ_u q(1 − c²) at g_φ = g_ψ = 0, and Q_new = F_eff − c²d1d2 − 2Λq²(1 − c²) in that case.
4. Coupling through u = A/√X: δA₀ drops out of u at linear order.
5. G1-standard in 1+1D (§4.3a): the momenta, H ≃ ½π² − A₀∂₁π + V, the secondary constraint Ω₂: ∂₁π = 2A₀V′, the reduced H ≈ ½π² + V − 2A₀²V′, the second-class determinant D = 2V′ + 4A₀²V″, its form D = 4√X g′(X) on A₁ = 0 together with the Rolle argument, the exact zero of ∫2A₀V′ for the odd profile, the vacuum density V(x), and the plateau energy relative to the vacuum, E_rel = 2ℓ[V(X*) − V(x)] + A²vℓ(−4 + (4/3)vℓ²) + O(1). Verified symbolically (27/27, script 152a7e39…).
6. G1-restricted: on X = x, δA₀ = 0 at linear order about Ā = (a, 0).

## References

- Registered: `docs/UV_STRUCTURAL_CONSTRAINTS_NOTE_v0.3.md` (c7b92ba); `docs/UV_DYNAMIC_ETA_VS_DYNAMIC_U_COMPARISON_v0.3.1.md`; `docs/UV_DYNAMIC_ETA_SCRIPT_A_RESULTS_NOTE_v0.2.md` (636d278); `docs/PROCA_TIMELIKE_QNEW_DECISION_NOTE_v0.1.md` (73347c3).
- R. Bluhm, N. L. Gagne, R. Potting, A. Vrublevskis, *Constraints and Stability in Vector Theories with Spontaneous Lorentz Violation*, arXiv:0802.4071.
- S. M. Carroll, T. R. Dulaney, M. I. Gresham, H. Tam, *Instabilities in the Æther*, arXiv:0812.1049 (fixed-norm æther).
- Q. G. Bailey, R. Potting, D. T. Walter-Cardona, *Stability and dynamics of bumblebee models*, arXiv:2610.21066 (2026): reduced Hamiltonian H_V = V − 2B₀²V′ and the necessary condition V′ ≤ 0 (§II, Eqs. 13–14); restricted set X = x and Maxwell in a nonlinear gauge (§IV A, §IV B, §V A).
