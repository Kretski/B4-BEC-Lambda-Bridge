# G1-standard: lemma on the layers and the O(1) remainder — v0.3 (separate supplement, for independent review)

**Status:** proposal for review. It is **not** accepted as closing the O(1) question. It does not modify `DYNAMIC_U_UV_ADMISSIBILITY_NOTE_v0.6.md` (899bc4ef…), which stays unchanged and unregistered.

**Base:** 73347c3.

**Changes from v0.2 (review of 2026-11-19 10:58).** One precision only. The amplitude range of L3 is stated as A ≥ 4 **and** A > √x. The construction of v0.6 needs A > √x, and the estimate of item 1 needs A ≥ 4; the condition was previously written as A ≥ 4 alone. The asymptotic result is unchanged.

**Changes from v0.1 (review of 2026-11-19 10:46).** L1 and L2 were accepted; L3 was accepted in substance. Two corrections:

1. **L3:** the two sup bounds C_out(A) and C_c(A) are now defined explicitly, as the supremum of one and the same function F along the respective path. This function contains both the potential terms and the π-source term, so one C(A) controls every estimate used.
2. **L1:** the C^∞ claim is stated as an explicit condition. Each piece of the path must be C^∞, and the reparametrisation must be flat to all orders at every junction. Continuity of γ alone is not enough.

Nothing else changes.

**Scope:**
- the admissibility of the family of §4.3a on the full constraint set;
- the regularity of the **outer** layers;
- the energy of all layers and the O(1) remainder.

**Not covered:** the regularity of the **central** layer, which remains condition (R).

**Notation** as in v0.6 §4.3a: (+,−); X = A₀² − A₁²; vacuum Ā = (±a, 0), a² = x; V′(x) = 0, V″(x) < 0; reduced density H = ½π² + V(X) − 2A₀²V′(X); Ω₂: ∂₁π = 2A₀V′(X); D = 2V′ + 4A₀²V″. Energy is measured relative to the vacuum: H_rel = H − V(x). X*, v and ℓ are fixed as in step 6 (6.1 and 6.2), and V is C².

## L0. What "admissible initial data" means here

A configuration is the triple (A₀(s), A₁(s), π(s)) on the line, together with π⁰ = 0. It is admissible on the **full constraint set** if both constraints hold pointwise (π⁰ = 0 and Ω₂) and the relative energy is finite. It is admissible on the **regular surface** if, in addition, D(s) ≠ 0 for every s.

Where D ≠ 0, the implicit function theorem lets Ω₂ be solved locally for A₀ in terms of (A₁, ∂₁π). This is the usual second-class reduction, and the data (A₁, π) are then the free Cauchy data.

## L1. Construction from a path (for review: is every such profile admissible?)

Let γ: [0, 1] → ℝ², τ ↦ (α₀(τ), α₁(τ)), be a continuous path made of finitely many pieces, **each of them C^∞ up to its end points**. The paths used here are real-analytic pieces: segments of A₁ = 0, branches of the hyperbola X = X*, and straight segments.

Let φ: [s₀, s₀ + w] → [0, 1] be a C^∞ monotone reparametrisation that is **flat to all orders** (all derivatives zero) at the end points and at every junction of the pieces. Set (A₀, A₁)(s) = γ(φ(s)).

**Condition for C^∞:** both requirements together, C^∞ pieces and flatness at every junction. Under them the profile is C^∞ in s, even where γ has a corner: the "stop-and-go" parametrisation. Continuity of γ alone does not give C^∞. Without these requirements only a continuous profile results; the energy estimates of L3 still hold for it, but not the smoothness claim.

**π is not free.** It is **defined** by Ω₂:

π(s) = π(s₀) + ∫_{s₀}^{s} 2A₀V′(X) ds′.

The integrand is continuous, so π is C¹ and Ω₂ holds pointwise. π⁰ = 0 holds identically. So every piece built this way, with π defined by Ω₂, lies on the full constraint set.

**The global condition** is π → 0 at ±∞ for finite energy. For the odd profile it holds exactly, because A₀ is odd and A₁ even, so A₀V′(X) is odd and ∫2A₀V′ = 0 (v0.6, step 7; check 5i).

So the only conditions on the profile are continuity, the C^∞ reparametrisation and the global oddness. No gradient energy constrains it, because H contains no ∂A₀ or ∂A₁.

*Point for the reviewer:* this rests on the reduced H. The term π∂₁A₀ was removed by integration by parts and replaced through Ω₂ (v0.6 §4.3a, steps 2–4). On the full constraint set this is an identity, so no hidden gradient term remains. It assumes compact support of π·A₀ − (vacuum) at the ends. That holds: π = 0 outside the configuration, and the boundary term π·A₀ vanishes there.

## L2. Regularity of the outer layers (D < 0 along the whole outer path)

**Outer path** γ_out, from the vacuum (a, 0) to the plateau point P_A = (A, √(A² − X*)):

- **Segment S:** A₁ = 0, A₀ from √x down to √X*. Here X ∈ [X*, x].
- **Hyperbola branch H_A:** X = X*, A₀ from √X* up to A.

The left outer path is the mirror image (A₀ → −A₀).

**On S.** D = 2V′(X) + 4XV″(X) with X ∈ [X*, x]. On [X*, x], V″ < 0 by 6.1. Also V′(X) = −∫_X^x V″ ≥ 0, with V′(X) ≤ |V″|_max (x − X). Hence

D ≤ 2|V″|_max(x − X) − 4X·min|V″| < 0,

provided x − X* < 2X*·min|V″|/|V″|_max on [X*, x]. This holds for X* close enough to x, since min|V″| and |V″|_max both tend to |V″(x)| > 0 as X* → x. **Requirement 6.1′:** this condition is added explicitly to step 6.1. It is implied by "X* close enough to x", but v0.6 did not state it. At X = x, D = 4xV″(x) < 0, which is the vacuum.

**On H_A.** X = X* is fixed, so D = 2v + 4A₀²V″(X*). This is decreasing in A₀² because V″(X*) < 0. At the vertex A₀² = X* it equals 2v + 4X*V″(X*) < 0 by 6.1. Hence D < 0 on all of H_A, for every A.

**At the plateau point.** D = 2v + 4A²V″(X*) < 0, by the same bound. Step 6.3 of v0.6, which excludes at most two values of A, is therefore **not needed** once 6.1 holds: D < 0 for every A ≥ √X*.

**Conclusion of L2.** For every A, the outer layers and the plateaus lie in {D < 0}. The vacuum ends have D = 4xV″(x) < 0. The central layer is not covered here; its regularity is condition (R).

## L3. Energy bound with an explicit ε(A)

**Plateau quantities.** On the plateaus |A₀| = A, X = X*, |∂₁π| = 2Av. The plateau excursion of π has magnitude at most 2Avℓ.

**The function to bound.** Along any path in the (A₀, A₁) plane define

F(A₀, A₁) ≡ |V(X) − V(x)| + 2A₀²|V′(X)| + 2|A₀V′(X)|,  X = A₀² − A₁².

The first two terms bound the potential part of the relative energy density. The third bounds the π source, since |∂₁π| = |2A₀V′|. A single bound on F therefore controls both the layer energy and the change of π.

**Explicit definitions.**

C_out(A) ≡ sup over γ_out(A) ∪ γ̃_out(A) of F,  C_c(A) ≡ sup over γ_c(A) of F.

Here γ_out(A) is the right outer path of L2 (segment S, then the hyperbola branch H_A up to P_A), γ̃_out(A) is its mirror image (A₀ → −A₀), and γ_c(A) is the central path. Both suprema are finite for each A, because the paths are compact and V, V′ are continuous.

**Outer path, explicit growth.**
- On S, A₀ ≤ √x and X ∈ [X*, x]. So F is bounded there by a constant K_S that does not depend on A.
- On H_A, X = X*, so V and V′ are constant there and |A₀| ≤ A. Hence F ≤ |V(X*) − V(x)| + 2A²v + 2Av.

Therefore C_out(A) ≤ K_S + |V(X*) − V(x)| + 2A²v + 2Av = O(A²).

**Central path.** On the full constraint set, γ_c(A) can be taken as the straight segment from (A, √(A² − X*)) to (−A, √(A² − X*)), at fixed A₁. C_c(A) is finite for each A. Its growth in A depends on V away from the vacuum, including the spacelike region crossed at A₀ = 0, and is **not** controlled. It does not need to be: ε is chosen after A.

**π inside the layers.** |π| ≤ 2Avℓ + Σ|δπ| along the configuration.

**Choice of the layer width.** Define

C(A) = max{ C_out(A), C_c(A), ½(2Avℓ + 1)², 1 },  ε(A) = 1/(A·C(A)).

So on every layer F ≤ C(A), and, once the π bound of item 1 holds, ½π² ≤ C(A). The same C(A) is used in every estimate below.

Then:

1. **π change across the layers:** |δπ| ≤ width·sup F ≤ width·C(A). The three layers have widths ε, 2ε and ε, total 4ε. So the total change is at most 4εC(A) = 4/A ≤ 1 for A ≥ 4 (A > √x is assumed throughout, as in the construction of v0.6). Hence |π| ≤ 2Avℓ + 1 everywhere, and the offset entering any plateau is |δ| ≤ 4/A.
2. **Relative energy of the layers:** pointwise ½π² ≤ ½(2Avℓ + 1)² ≤ C(A), and |V − V(x)| + 2A₀²|V′| ≤ F ≤ C(A). So ∫_layers [½π² + |V − V(x)| + 2A₀²|V′|] ds ≤ (total width 4ε)·2C(A) = 8/A.
3. **Plateau correction.** On a plateau π = π₀(s) + δ, where π₀ is the v0.6 plateau profile and the constant offset satisfies |δ| ≤ 4/A. The change of ∫½π² over the two plateaus is at most 2ℓ(|π₀|_max|δ| + ½δ²) ≤ 2ℓ(2Avℓ·4/A + 8/A²) = 16vℓ² + 16ℓ/A².
4. **Vacuum tails:** relative energy 0.

**Result (L3).** For every A with A ≥ 4 and A > √x,

|E_rel(A) − {2ℓ[V(X*) − V(x)] + A²vℓ(−4 + (4/3)vℓ²)}| ≤ 8/A + 16vℓ² + 16ℓ/A².

The remainder is bounded by explicit constants, uniformly for A ≥ 4, A > √x: **O(1)**. Since ℓ² < 3/v, E_rel(A) → −∞.

## L4. What this supplement claims, and what it does not

**Claimed, for independent checking:**
- (i) Every profile of L1 lies on the full constraint set.
- (ii) The outer layers and the plateaus are regular (D < 0) for every A, under 6.1 and 6.1′.
- (iii) The remainder bound of L3 has explicit constants.

**Not claimed:**
- The regularity of the central layer (condition (R)).
- Anything about linear stability or nonlinear evolution.

**Changes to v0.6 suggested if accepted** (not made here):
- add requirement 6.1′;
- note that step 6.3 is redundant;
- refer to this supplement for step 8.

**Points where the reviewer's check matters most:**
- L1: no hidden gradient term, and the boundary terms of the integration by parts vanish.
- L3, item 1: the bound on π inside the layers, which is circular-looking (it uses ε, chosen with the plateau value of π). It closes because C(A) ≥ ½(2Avℓ + 1)² is fixed before ε.
- L3, item 3: the cross term.
