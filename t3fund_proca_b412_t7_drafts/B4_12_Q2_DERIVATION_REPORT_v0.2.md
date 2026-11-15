# B4.12 Q2 — Validity of the M⁻² approximation: derivation report (draft v0.2)

**Status: DRAFT FOR REVIEW. Not committed.** This report states the results of the registered Q2 derivation. It adds no computation, no numerical value and no physical conclusion beyond `results_b412_q2/b412_q2_results.json` and the registered specification `B4_12_Q2_DERIVATION_SPEC_v0.2.md`. **B4.12 remains OPEN.**

*Changes from v0.1 (review of 2026-11-15, 17:36): wording only. The three requested places are §2.6 (twice) and §2.8. In addition, "independent methods" in §5 is changed to "distinct methods", for consistency with §2.8. No result changed.*

## 0. Artefacts

| Item | Reference |
|---|---|
| Specification | `t3fund_proca_b412_t7_drafts/B4_12_Q2_DERIVATION_SPEC_v0.2.md`, registered at `641121c` |
| Script | `t3fund_proca_b412_t7_drafts/b412_q2_derivation.py`, revision 7, commit `a7aec47`, SHA-256 `68b47e47…346cf` |
| Run outputs | `t3fund_proca_b412_t7_drafts/results_b412_q2/`, commit `e8852dd`; `b412_q2_results.json` SHA-256 `24dcc0e6…9dea` |
| Registration record | `t3fund_proca_b412_t7_drafts/B4_12_Q2_RUN_REGISTRATION.md`, commit `35465c2` |
| Environment | Python 3.13.18, SymPy 1.14.0; 11.4 s; exit code 0; script verdict `ALL_CHECKS_PASS` (12 mandatory checks) |

Below, "JSON key" refers to `results` or `checks` in `b412_q2_results.json`. Expressions too long to reproduce here are cited by key.

## 1. Setting (from spec v0.2, unchanged)

- **Action and conventions.** The registered T3-fund Proca action in 1+1D, metric (−,+), with generic couplings g_φ, g_ψ.
- **Proca matching.** g_φ = M√p and g_ψ = −Mκ/√p, so g_ψ²/M² = κ²/p.
- **Background.** Constant gradients: f = wt + qx, g = vt + rx (ν = 0).
- **Light relation.** Written in c = ω/k and u = (k/M)², at fixed p, κ, Λ and fixed background. Therefore g ∝ M.
- **Method.** Symbolic only. No registered background was evaluated.

## 2. Derived results

All items in this section come from the script's symbolic computation and its checks (§5 says what that does and does not establish).

### 2.1 D1 — background and quadratic system
- The background Proca field is **solved** from the O(ε) equations: Ā0 = 0 and Ā1 = (g_φ q w + g_ψ r v)/M² (key `D1_background_A`).
- The light O(ε) Euler–Lagrange expressions are identically zero (`D1_light_background_EOM_zero`). They are the only identically zero equations recorded: `df in L1_background` and `dg in L1_background` (`identically_zero_EL_equations`).
- The O(ε²) fluctuation system is homogeneous (`D1_homogeneous_fluctuation_system`).

### 2.2 Q2-C1 — kernel
The kernel derived from the Proca block with an external source equals the T7a kernel identically (`Q2_C1_kernel`):

  K(ω,k) = (M² + k²) / (M⁴ + M²k² − M²ω²) = (1 + k²/M²)/(M² + k² − ω²)

(key `Q2_C1_K_derived`).

### 2.3 Q2-C2 — hypothesis H-Q2
For ν = 0, the reduced 2×2 operator obtained by tree-level elimination of a_μ (Schur complement) equals M₀ + K(ω,k)·N. Here M₀ + c_loc·N is the operator of the local Lagrangian L₂^light + (c_loc/2)·J₁², and L₂^light includes Ā1·J₂ = J̄J₂/M². The remainder is identically zero in all four entries (`Q2_C2_HQ2`, `Q2_C2_remainder_srepr`).

**Hence H-Q2 holds:** at tree level and for constant-gradient backgrounds, the exact quadratic action is S₂^light + ½·J₁KJ₁ + J̄·J₂/M².

### 2.4 Q2-C3 — relation to B4
1. The Proca M⁻² truncation equals the B4 EFT fluctuation operator of the **rank-one subfamily**, s = κ²/p, identically (`Q2_C3_rank_one_match`). The comparison is with B4 at m = 0, because the registered Proca action has no ψ mass term.
2. For general B4 with independent s, the difference is non-zero and proportional to (s − κ²/p). It is confined to the g–g entry and equals (s − κ²/p)·(k²v² − 4kω r v + ω²r²) in the script's convention (`Q2_C3_general_B4_differs_by_(s-kappa^2/p)`, key `Q2_C3_general_B4_difference_over_(s-kappa^2/p)`).

### 2.5 Q2-C4 — decoupling of the full system
For g_φ = g_ψ = 0, the full 4×4 system (δf, δg, a0, a1) decouples (`Q2_C4_decoupling`) into four signed light branches and the free Proca branch (key `Q2_C4_branches`):

| Block | Roots |
|---|---|
| light f | ω = ±k√(1 + 3Λq²) |
| light g | ω = ±k |
| Proca | ω = ±√(M² + k²) |

### 2.6 D2 — dispersion relations and the heavy branch
- **Full system:**
  - P_full has degree 6 in ω (key `D2b_full_polynomial_degree_in_omega`);
  - det M_aa = −M²(M² + k² − ω²) (key `D2b_Daa`);
  - the Schur identity P_full = det M_aa · det M_red holds (`D2b_schur_identity`).
- **Poles of K:** P_full does not vanish at ω = +√(M² + k²) or at ω = −√(M² + k²) (key `D2b_Pfull_at_pole_is_zero`: both false), and the numerator of the reduced relation is not divisible by M² + k² − ω² (`D2b_reduced_numerator_divisible_by_pole`: false). Within this symbolic check, i.e. as identities in the symbols of the computation (generic parameters), no root of the reduced relation lies at a pole of K.
- **Light relation:** F(c, u) has degree 6 in c and 1 in u; F0(c) = F(c, 0) has degree 4 (keys `light_relation_F_c_u`, `F0_c_(M^-2 light polynomial)`, `F0_degree_in_c`). Four light branches are therefore given by the roots c0 of F0.
- **Heavy branch:** defined by continuation from the decoupled Proca root, with coupling scale η (key `D2b_heavy_branch`).
  - For each sign, the decoupled root ±√(M² + k²) is a root at η = 0.
  - The derivative of the characteristic polynomial with respect to ω, ∂P_full/∂ω, evaluated at the heavy root and η = 0, is ∓2M⁴√(M² + k²)·(3Λk²q² − M²). This is not identically zero, so the root is simple for generic parameters (`D2b_heavy_continuation_defined_generic`).
  - The leading correction is of order η². The coefficient a_n is recorded in the JSON and is not identically zero.

### 2.7 Q2-C6 — the factor M²K − 1
In the light variables, M²K − 1 = u c²/(1 + u(1 − c²)) identically (`Q2_C6_M2K_minus_1`). This is the identity (ω²/M²)/(1 + (k² − ω²)/M²) of B4.12 v0.4 §3.2.

### 2.8 D3 — expansion of the light branches and error coefficients
- **Expansion:** for a simple root c0 of F0, the light branch is c(u) = c0 + c1·u + c2·u² + O(u³). c1 and c2 are given as rational functions of c0 and the background parameters (keys `D3_c1`, `D3_c2`).
- **Q2-C5:** two distinct analytic derivation routes give the same c1 and c2 modulo F0(c0) = 0 (`Q2_C5_order_u1`, `Q2_C5_order_u2`). Route A uses implicit differentiation; route B uses perturbation theory on the dispersion polynomial.
- **Leading error coefficients**, for c0 ≠ 0 and c0² ≠ 1 (keys `D3_eps_EFT_leading`, `D3_eps_tot_leading`, `D3_eps_abs_leading`):
  - ε_EFT = |2c0c1/(c0² − 1)|·u + O(u²);
  - ε_tot = |2c1/c0|·u + O(u²);
  - ε_abs = |2c0c1|·u² + O(u³).
- **Special cases** (key `D3_special_cases`):
  - **c0² = 1:** ε_EFT is undefined, and only ε_abs is reported. F0(±1) is given symbolically.
  - **c0 = 0:** this branch exists iff F0(0) = 3Λκ²q²v² − 3Λpq² + κ²v² + p²w² − p = 0. There ε_tot = 1 exactly whenever c(u) ≠ 0, and it is undefined if c(u) ≡ 0.
    - Conditional on c1 ≠ 0 (c1(c0 = 0) given symbolically, valid where F_c(0,0) ≠ 0): ε_EFT = c1²u² + O(u³) and ε_abs = c1²u³ + O(u⁴).
    - If c1 = 0, the case is **unresolved**.

### 2.9 D3 (iii) — dependence on ρ at fixed X
Under the scaling (w, q, v, r) → λ(w, q, v, r), (p, κ, Λ) → (p, κ, Λ)/λ², every X_j = (g_jρ_j)² and Λq² stays fixed while ρ changes by λ at fixed M. Under this scaling the light relation changes only by the overall factor λ⁻² (`D3_iii_F_invariant_at_fixed_X`, key `D3_iii_overall_scaling_factor`). The light branches, and hence the error coefficients, are invariant under this transformation. The statement is limited to this stated scaling transformation.

## 3. Conditions of validity

The four limitations, as registered:
1. **The heavy branch is defined only for generic, non-resonant parameters.** The exceptional set is M² = 3Λq²k², where the heavy root ω² = k² + M² meets the light f branch ω² = (1 + 3Λq²)k². The continuation is not established there.
2. **D4(a) and D4(c) are not derived** (§4).
3. **D5 holds only for generic simple light branches** (c0 ≠ 0, c0² ≠ 1). The special branches and degenerate roots are excluded.
4. **A single Proca mediator reproduces the rank-one subfamily of B4, not general B4** (§2.4).

Further conditions of the derivation:
- tree level;
- 1+1D;
- constant-gradient backgrounds (ν = 0);
- the fixed-EFT prescription with u varied at fixed p, κ, Λ and fixed background;
- simple roots (F_c(c0, 0) ≠ 0) for D3.

## 4. Not derived

- **D4(a) — reality of the light branches at finite k/M: NOT DERIVED.** Only the limit k/M → 0 is covered: real light branches there iff all roots of F0(c) are real. T3-fund C3 is cited only for its registered grid, not as a proof here.
- **D4(c) — separation of light and heavy branches: NOT DERIVED.** The pole condition is not claimed to be sufficient for it.
- **D4(b)** is a **necessary** condition only: 1 + u(1 − c_b²) > 0 along a light branch, i.e. (k/M)² < 1/(c_b² − 1) for c_b² > 1, evaluated at leading order c_b → c0. It is not a sufficient condition for spectral admissibility.

## 5. Limits of the verification

- `ALL_CHECKS_PASS` is the verdict of the script's own checks. It means that the 12 mandatory checks evaluated as stated. It is **not** an independent proof of every physical statement.
- **What the checks establish:** exact symbolic identities (kernel, H-Q2, rank-one match, decoupling, Schur identity, the M²K − 1 identity) and the agreement of two analytic routes for c1, c2.
- **What they do not establish:** an independent re-derivation with a second code. The reviewer checked the integrity and internal consistency of the run package, not a re-execution.
- **The two routes in Q2-C5** are distinct methods, but they share the same light relation F(c, u). An error in building F would affect both.
- **Zero tests:** they act on numerators modulo F0, and their interpretation as equality of rational functions assumes non-vanishing denominators. In particular they assume F_c ≠ 0 and the non-degenerate cases stated above.

## 6. Outcome for Q2 (definitions of B4.12 v0.4 §3.4)

**Classification:** outcome **(ii)**, on generic simple light branches (key `D5_classification`, scope `D5_scope`). The leading ε_EFT is of order (k/M)², and its coefficient is not identically zero. At fixed X-invariants it does not depend on ρ separately. By the definition in v0.4 §3.4, outcome (ii) means that a ρ threshold is not a standalone criterion.

Not concluded here:
- no value of ε_max, τ or a grid;
- no evaluation on the C4, E3 or population backgrounds;
- no statement on unitarity (Q3) or on the joint Q2/Q3 check;
- no change of the B4.12 status, which remains **OPEN**.
