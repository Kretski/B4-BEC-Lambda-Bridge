# UV structural constraints on B4: preferred frame, single-mediator rank, tree-level sign of Λ — v0.3

**Status:** v0.3, for registration. v0.2 was accepted by the independent review of 2026-11-17 subject to two editorial corrections (parity-even scope in A; mediator-exchange scope in E), both applied here. The mathematics and the checks are unchanged from v0.2. Not yet registered. No numerical run; the claims are algebraic.

**Evidence level:** 1 (theoretical). The companion script checks the algebra symbolically (17 checks, 0 failed, see §6). It has no data, no thresholds and no tolerances.

**Purpose:** fix what follows already, from structure alone, for any UV candidate of B4. The aim is to choose the next candidate on this basis instead of by another large numerical run.

## 0. Conventions

- Signature (+,−), 1+1D.
- In the frame of u, write X₁ = u·∂φ, X₂ = n·∂φ, X₃ = u·∂ψ, X₄ = n·∂ψ.
- B4, from `docs/notation_and_conventions.md` with m = 0:

  L_B4 = ½(X₁² − X₂²) − (Λ/4)X₂⁴ + ½(X₃² − X₄²) + (p/2)X₁²X₂² − κ X₁X₂X₃X₄ + (s/2)X₃²X₄².

- "Healthy mediator" means positive kinetic term and M² > 0.
- "Tree level, leading local order" means the classical elimination of the mediator with its gradient terms dropped. Those terms are O(∂²/M⁴); Script A, S7.1, derived the remainder.

Each statement below is labelled **Proven** (follows from the stated premises), **Conditional** (true only under named extra conditions) or **Open** (needs further work).

## A. A preferred frame is required

**Statement A (Proven).** Take parity-even local quartic polynomials in the first derivatives of two scalars, invariant under the full Lorentz group O(1,1). Every such polynomial is a combination of X_φ², X_ψ², Y², X_φX_ψ, X_φY and X_ψY, where X_φ = (∂φ)², X_ψ = (∂ψ)², Y = ∂φ·∂ψ.

The B4 quartic part lies in this span only for p = s = κ = Λ = 0 (check A1). The same holds for each B4 quartic operator taken on its own (checks A2: X₂⁴, X₁²X₂², X₁X₂X₃X₄, X₃²X₄²).

Short reason: X_φ² contains X₁⁴ with the same coefficient as X₂⁴, and B4 has no X₁⁴ term. The mixed sector works the same way.

**Proper group only.** Under the proper group SO⁺(1,1) alone there is also the pseudoscalar E = ε^{μν}∂_μφ ∂_νψ = X₁X₄ − X₂X₃. Products such as E·X_φ are invariant but parity-odd. They cannot contribute to the B4 operators, which are all parity-even (under x → −x, X₂ and X₄ change sign; every B4 term has an even number of them). The parity-even product E² = Y² − X_φX_ψ (an identity in 1+1D) is already in the span above. So the conclusion of A also holds for SO⁺(1,1).

**Scope of Statement A.** It classifies local quartic polynomials in the first derivatives of two scalars in 1+1D. It is not a statement about every Lorentz-invariant theory: non-local actions, additional fields and spontaneous symmetry breaking are outside it (see Conditional below).

**Consequence.** Within that class, no Lorentz-invariant action reproduces any non-trivial B4 operator. This agrees with the registered result for HeavyScalar_TwoField (T1/T7: X, A in `docs/UV_CANDIDATE_ASSESSMENT.md`).

**Conditional.** Statement A concerns the action. A preferred direction can still enter in three ways:

1. an externally given field u^μ;
2. a dynamical field with a non-zero background vector (dynamic η / u^μ);
3. spontaneous breaking, i.e. a Lorentz-invariant theory expanded about a background with a timelike gradient, which then defines the frame.

**Requirement.** Every successful UV realisation must explain where its preferred frame comes from. This is the open item T6.

## B. One mediator gives rank one

**Statement B1 (Proven — scalar mediator).** Take one healthy scalar mediator χ, linearly coupled through a single scalar source J = Σᵢ gᵢOᵢ, at tree level and leading local order:

L_eff = J²/(2M²)   (check B1).

The coefficient of OᵢOⱼ is gᵢgⱼ/M², an outer product, so of rank one. With O₁ = X₁X₂ and O₂ = X₃X₄, matching to B4 gives

p = g₁²/M²,   s = g₂²/M²,   κ = −g₁g₂/M²   (check B2),

and therefore **κ² = ps exactly** (check B3).

This is a property of one scalar mediator with a single linear source. Dynamic η (Script A, S7.3) reaches exactly this rank-one subfamily, as B1 predicts, because σ is a scalar.

**Statement B1-V (Conditional — vector mediator).** For a vector mediator the source is J^μ = Σᵢ gᵢ^μ Oᵢ, and the coefficient matrix is proportional to the Gram matrix gᵢ^μ g_{jμ}. Its rank is **not** fixed by the number of mediators: in 1+1D it can be up to two.

- With a common polarisation (all Oᵢ along one direction), one Proca mediator again gives κ² = ps (check B8).
- With currents along independent directions it does not. Example: J^μ = g₁u^μO₁ + g₂n^μO₂ gives p = −g₁²/M², s = g₂²/M², κ = 0, so κ² ≠ ps (check B7). Here p < 0, i.e. the time component enters with the opposite sign.

The rank-one result of T7a therefore comes from its specific current, J^μ = (0, J): a single spatial component, a common polarisation (T7a manuscript draft v0.3, definition of the coupling A₁J = A_μJ^μ with J^μ = (0, J)). It does not follow from the number of mediators. The frozen T7a manuscript draft (v0.3, discussion of the rank-one relation) states that "a single Proca mediator at tree level reaches only the rank-one subfamily". That sentence is broader than B1-V allows and should be narrowed to "with this current" in any later revision of the manuscript.

**Statement B2 (Proven).** With N healthy scalar mediators,

ps − κ² = Σ_{a<b} (g₁ₐg₂_b − g₁_b g₂ₐ)²/(Mₐ²M_b²) ≥ 0   (check B4, written out for N = 3; Cauchy–Binet).

So healthy scalar exchange at tree level reaches only **κ² ≤ ps**. Equality holds for one mediator, or for parallel coupling vectors. Reaching general (p, s, κ) with κ² < ps needs at least two independent mediators. **κ² > ps cannot be reached** by healthy scalar exchange at tree level.

**Conditional.**

- The positivity uses healthy mediators. A wrong-sign mediator flips the sign of its contribution (check B5).
- For vector mediators the contraction has an indefinite metric, so B2 does not carry over (see B1-V and check B7: p < 0 is reachable).
- A direct contact term is not restricted by B1 or B2.

**Statement B3 (Proven).** If the same mediator also couples to (n·∂φ)², it generates X₁X₂³, which is not a B4 operator (check B6). Couplings outside the B4 operator set create extra operators of the same order.

**Scope.** B1–B3 concern the leading local quartic structure. They make no claim about higher EFT terms or loops.

## C. Tree-level sign of Λ

**Statement C1 (Proven).** Let a healthy scalar mediator couple to any bilinear O = αX₁² + βX₁X₂ + hX₂². Then the coefficient of X₂⁴ in L_eff is h²/(2M²) ≥ 0. Comparing with −(Λ/4)X₂⁴:

**Λ_eff = −2h²/M² ≤ 0**   (check C1).

For N healthy scalars, Λ_eff = −2 Σₐ hₐ²/Mₐ² ≤ 0 (check C2).

**Consequence.** Within tree-level exchange of healthy scalar mediators, **Λ > 0 cannot be generated.**

**Statement C2 (Conditional — vector mediator).** For a Proca mediator at leading local order:

L_eff = −J·J/(2M²)   (check C3).

With J^μ = γ_u u^μ (n·∂φ)² + γ_n n^μ (n·∂φ)²:

**Λ_eff = 2(γ_u² − γ_n²)/M²**   (check C4).

So the sign is positive for a current along u and negative for a current along n. **Λ > 0 at tree level is therefore not excluded in general.** It is possible with a vector mediator whose current has a component along the preferred time direction.

This is a mechanism for the sign of **one** operator. A positive sign does **not** establish:

- stability;
- perturbativity;
- a full operator-by-operator match with B4.

It is not a working UV realisation of B4.

This is the static exchange of the time component, the same sign mechanism as Yukawa repulsion between like charges.

**Open, for C2.**

- (i) A non-conserved current couples to the longitudinal Proca mode. This lowers the strong-coupling scale and is tied directly to B4.12 (perturbativity).
- (ii) Stability (T3-fund, T5) of such a coupling has not been checked.
- (iii) Loop corrections are not considered.

**Open, generally.**

- The full sign conventions must be compared operator by operator for any concrete candidate.
- The X₂⁴ coefficient is fixed only in the on-shell operator basis used here; mixing with other same-order operators and field redefinitions have not been analysed.
- Three statements must be kept apart:
  1. the sign of the Gross–Pitaevskii quartic term in its own normalisation (`docs/UV_CANDIDATE_ASSESSMENT.md` records Λ_BEC < 0);
  2. the sign of Λ_eff obtained from a UV mediator (this note);
  3. the measured BEC dispersion.

  They are not equivalent. Statement 1 stays a separate hypothesis until its exact matching to the B4 normalisation is checked. Measured BEC dispersion data (e.g. Shammass et al. 2012) are not evidence for C1 or C2. A dispersion correction of either sign would first have to be separated from the standard 2D GPE / LDA control, and even then it would not identify B4 or a specific mediator.

## D. Summary by status

| Statement | Status | Premises |
|---|---|---|
| A: B4 operators are not Lorentz invariant; a frame is required | **Proven** | local quartic polynomials in first derivatives of two scalars, 1+1D (not non-local actions, extra fields or spontaneous breaking) |
| A: origin of the frame (external / dynamical / spontaneous) | **Open** (T6) | — |
| B1: one scalar mediator gives κ² = ps | **Proven** | one healthy scalar mediator, single linear source, tree level, leading local order |
| B1-V: one vector mediator: rank = rank of the Gram matrix of the current directions (≤ 2 in 1+1D); rank one only with a common polarisation | **Conditional** | Proca at leading local order |
| B2: N healthy scalars give κ² ≤ ps | **Proven** | healthy scalar mediators |
| B2 for vectors, contact terms, wrong-sign mediators | **Conditional** / not covered | — |
| B3: extra couplings give non-B4 operators of the same order | **Proven** | as B1 |
| C1: healthy scalars give Λ_eff ≤ 0 | **Proven** | as B1, scalar mediators |
| C2: Proca gives Λ_eff = 2(γ_u² − γ_n²)/M², so Λ > 0 is possible | **Conditional** | leading local order; says nothing about stability, perturbativity or the full operator match |
| Loops, operator mixing, field redefinitions, Λ_BEC | **Open** | — |

## E. Consequences for UV candidates

Every new candidate must state in advance:

1. **The origin of the preferred frame** (A, T6).
2. **Scalar or vector mediators, and how many.**
   - **Scalar:** one healthy scalar mediator with a single linear source J = Σᵢ gᵢOᵢ produces a rank-one coefficient matrix, and therefore κ² = ps for the specified B4 matching (B1). With N healthy scalars, κ² ≤ ps (B2). Scalars always give Λ ≤ 0 (C1).
   - **Vector:** the rank depends on the independent spacetime directions of the currents and must be evaluated explicitly (B1-V). The sign of each operator depends on the direction of the current (C2, B7). A current along u can give Λ > 0, but also flips the sign of the operators it multiplies. A timelike component therefore cannot simply be added to the T7a current on the assumption that the rest of the match survives; this must be derived.
3. **Its mechanism for the sign of Λ.** Scalar mediators only: Λ ≤ 0 (C1). Within the mediator-exchange mechanisms analysed here, a tree-level route to Λ > 0 is provided by a vector mediator with a suitable timelike current (C2). Independently justified contact interactions and non-tree-level mechanisms remain separate possibilities.
4. **Which extra operators of the same order it generates** (B3), and how they are treated.
5. **Any direct contact term** must be justified independently, not added to fit.

Status of the existing candidates:

| Candidate | Status |
|---|---|
| Dynamic η, Script A | Closed (636d278); rank-one, as B1 predicts |
| HeavyScalar, Lorentz-invariant (one or two light fields) | Not repeated: excluded by A, already registered T1/T7 negative |
| Single-mediator scalar with frame couplings | Restricted to κ² = ps and Λ ≤ 0 |
| Two or more independent scalar mediators with frame couplings | κ² < ps possible, Λ ≤ 0 still; does **not** address T7b with Λ > 0 |
| Proca with a current along u | The most concrete next candidate for **testing** T7b, not a successful candidate. Perturbativity (longitudinal mode, B4.12), stability and the full operator match are open; to be analysed analytically against B4.12 before any new code |
| Proca / T7a | Kept as control. Rank one because of its single spatial current (B1-V), not because it has one mediator. Its Λ is inserted by hand; its spatial current would give Λ_eff ≤ 0 by C2 |
| Dynamic u^μ | Open; frame origin plus stability |

**Implication for T7b.** Healthy scalar mediators, however many, cannot give Λ > 0 in this mechanism (C1). A Proca mediator with a current along u can give the right sign (C2). Two further routes stay open: a direct contact term (must be justified independently) and quantum corrections (need a separate loop analysis).

The next step is therefore analytic: a Proca–B4.12 review, with no new code, answering three questions:

1. the exact relation between the timelike current and the longitudinal Proca mode;
2. the perturbativity bound;
3. whether Λ > 0 can be obtained without violating the other stability and operator-matching conditions.

## 6. Symbolic checks

Script: `uv_structural_constraints_check.py` (SymPy only, no data), SHA-256 `956fa676…`, unchanged since v0.2. Its docstring names note v0.2; the checks are the same for v0.3.

```
python3 -I uv_structural_constraints_check.py
```

Result in the development check: 17 checks, 0 failed (A1, A2 ×4, B1–B8, C1–C4). Checks B7 and B8 were added in v0.2 for the vector-mediator rank.

One test-logic error was fixed during development: A2 initially treated coefficients absent from an operator as unsolved. The mathematics was not changed.
