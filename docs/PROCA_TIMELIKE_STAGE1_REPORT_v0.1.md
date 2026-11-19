# Proca with a timelike current — Stage 1 report v0.1

**Specification:** `docs/PROCA_TIMELIKE_CURRENT_ANALYTIC_SPEC_v0.5.md`, registered in 87fc31b. Followed without change.

**Status:** Stage 1 result, for review. Not registered.

**Evidence level:** 1 (theoretical).

**Symbolic checks:** `proca_timelike_stage1_check.py` (SHA-256 `b59393d2…`). Output `stage1_output.txt` (`138550c5…`): **33 checks, 0 failed**. Pure algebra. None of the 41 backgrounds is used, and there is no numerical threshold.

During development, one check (D1d) first failed. That was an error in my check: it expected the whole k² coefficient of K_lead[f,f] to be −(1 + 3Λ_eff q²), but that coefficient also contains the contact term p w². The check now tests only the γ_u-dependent part, as D1 requires. The full kernel equality D1c had passed all along. No mathematics was changed.

**Regression.** With γ_u = 0 and a bare Λ, the jet-Hessian route reproduces the registered T3-fund identity det K(ck, Mx) = −M⁸x⁴(1 + x²)Q(c; ε) exactly (check R1). So the conventions and the action are those of the registered specification.

## Summary

| Item | Result | Status |
|---|---|---|
| D0 | Λ_eff = 2γ_u²/M² > 0 for every real γ_u ≠ 0 | **Proven** |
| D1 / H1 | no J_nJ_u cross term; p, s, κ exactly as T7a (rank one); no non-B4 operator at O(M⁻²); full leading kernel = B4 kernel with Λ_eff; D1 consistency check passed | **Proven** |
| D2 | T7a term reproduced; two new terms; no undecided terms | **Proven** |
| D3 | closed-form k_der per branch; R-UNDEFINED test; H3 decided per branch by an explicit inequality, not universally | formula **proven**; values need Stage 2 |
| D4 / H4 | G′ = G, c′_max ≥ c_max: the B4.12 counts cannot improve | **Proven** (from the registered definitions) |
| D5 | the timelike current changes the principal symbol; principal speeds = roots of an explicit degree-6 polynomial Q_new(c; 1) | formula **proven**; hyperbolicity per background undetermined |
| **Class (§3)** | **OPEN** | V2 and V3 undetermined analytically; NOT-VIABLE excluded except possibly through (d) |

## D0 / D1 — leading order (checks C0a–c, D0a–c, D1a–e)

- **Background.** The background Ā0 = −J̄_u/M², Ā1 = J̄_n/M², F̄ = 0 solves all four Euler–Lagrange equations. The Gauss constraint is M²A0 + J_u + F_x = 0.
- **Leading effective Lagrangian.** L_eff^(M⁻²) = (J_n² − J_u²)/(2M²). It has no cross term. This is H1.
- **f_x⁴ coefficient.** −γ_u²/(2M²), so Λ_eff = 2γ_u²/M².
- **Leading light kernel.** K_lead (free part + local exchange + background term) equals the B4 kernel with Λ = Λ_eff exactly (D1c). Its determinant is k⁴F_eff(c) (D1e).
- **D1 consistency.** The γ_u part of K_lead[f,f] is −3Λ_eff q²k². It is assembled from the background term Ā0·δ²J_u (−2γ_u²q²k²/M²) and the δA0 exchange (−4γ_u²q²k²/M²) (D1d).

## D2 — O(∂²/M⁴) (checks D2a–e)

The O(M⁻⁴) exchange is

L⁽⁴⁾ = [(∂_tJ_n)² + (∂_xJ_u)² + 2 ∂_tJ_u ∂_xJ_n]/(2M⁴)   (symbol identity, D2b).

| Term | Class | Reason |
|---|---|---|
| (∂_tJ_n)²/(2M⁴) | **1: T7a-equivalent** | identical to the T7a step-1 term; grading (t,x) = (4,2), the grading of the T7a canonical basis |
| (∂_xJ_u)²/(2M⁴) | **2: new** | grading (0,6); not a total derivative |
| ∂_tJ_u ∂_xJ_n/M⁴ | **2: new** | grading (2,4); not a total derivative; odd in γ_u |
| — | **3: undecided** | none |

IBP relations D_t(·) and D_x(·) are homogeneous in the (t,x) grading. So a term of grading (0,6) or (2,4) cannot be IBP-equivalent to a (4,2) term unless it is itself a total derivative. The Euler-operator check (D2e) shows it is not.

The effect of the two new terms on the dispersion is contained exactly in the kernel K_rem of D3, so their effect on B4.12 is decided there. **V4 holds.**

## D3 — k_der (checks D3a–g, H3a)

**Background term.** It sits in the light–light block K_LL, which is common to K_exact and K_lead. So it contributes only to the leading part, as expected in the specification.

**Remainder, exact (D3b).**

K_rem = u uᵀ / (M²(M² + k² − ω²)),

with u = (2γ_u k² q − g_φ ω(k w − ω q), g_ψ ω(ω r − k v)).

K_rem is **rank one**. The reason is that k_μk_ν − η_μν k² = k² e eᵀ with e = (1, −c), which is rank one in 1+1D (D3c).

**On a branch ω = c k (D3d, D3e).** With x = k²/M²:

- K_rem = k² · x/(1 + (1 − c²)x) · û ûᵀ/M², where û = (2γ_u q − c g_φ(w − c q), −c g_ψ(v − c r));
- K_lead = k² K̂(c), with det K̂(c_i) = F_eff(c_i) = 0.

Hence ‖K̂‖₂ = |tr K̂(c_i)|.

**Definitions.**

- N_i = |û(c_i)|² / (M² |tr K̂(c_i)|), dimensionless.
- R = N_i x / |1 + (1 − c_i²)x|.

**Result (D3f, D3g).** For each real branch c_i:

- **R-UNDEFINED** iff K̂(c_i) = 0, i.e. tr K̂(c_i) = 0 (together with det = 0). This is decided exactly.
- **No crossing (below)** if N_i + c_i² ≤ 1.
- **Crossing** if N_i + c_i² > 1, at k_der,i = M (N_i + c_i² − 1)^(−1/2).
- For |c_i| > 1 the propagator pole on the branch is at x_s = 1/(c_i² − 1). It always lies beyond the crossing when N_i > 0, and it is reported.

So k_der(Φ̄) = min over real branches of k_der,i. It is never 0 where it is defined. Therefore NOT-VIABLE (c) can only arise if every background is R-UNDEFINED.

**V2 in closed form.** k_der > M on Φ̄ ⟺ max over real branches of (N_i + c_i²) < 2.

A necessary condition is that every real leading-order speed satisfies |c_i| < √2. This is a statement about the formula, not a prediction. The values need the Stage 2 arithmetic.

**H3 (comparison with "T7a alone").** Read here as γ_u = 0 with bare Λ = Λ_eff, so K_lead and the branches are identical. Then

|û|²_with − |û|²_without = 4γ_u q [γ_u q − c g_φ(w − c q)]   (H3a).

So H3 holds on branch i iff γ_u q [γ_u q − c_i g_φ(w − c_i q)] ≤ 0. This is **not universal**. It depends on the branch, the background and the **sign of γ_u**. It is decided per background in Stage 2.

## D4 — perturbativity (checks D4a–b)

- |γ_u| = M√(Λ_eff/2). With ρ_q = |q|/M, this gives |γ_u|ρ_q = √(Λ_eff/2)|q|, the analogue of g_jρ_j = √X_j.
- J_u = γ_u f_x² contains only q, which is already in G, so **G′ = G**.
- c′_max = max(√p, |κ|/√p, √(Λ_eff/2)) ≥ c_max. Hence c′_max·G′ ≥ c_max·G: **H4 holds**.

The B4.12 counts for sets A, B and C cannot improve through this mechanism. At best, Proca with a timelike current inherits **B4.12-OPEN**. The per-background values are Stage 2.

## D5 — principal symbol (checks D5a–e)

The full 4×4 system (δf, δg, δA0, δA1) gives

det K(ck, Mx) = −M⁸x⁴(1 + x²) [F_eff(c) − ε(c² d1 d2 + 4γ_u q B(c)/M⁴)],   ε = x²/(1 + x²).

This is the registered form, plus one O(ε) term with the explicit polynomial

B(c) = M²[−c⁴g_φq + c³g_φw + c²g_φq − c²γ_u q − c g_φ w + γ_u q] + 2c³g_φ²g_ψq²w + 2c³g_φg_ψ²qrv − 2c²g_φ²g_ψqw² − 2c²g_φg_ψ²rvw − c²g_ψ²γ_u q r² + 2c g_φg_ψγ_u q² w + 4c g_ψ²γ_u q r v − g_ψ²γ_u q v².

- The leading-order dispersion (ε → 0) is exactly F_eff.
- **The timelike current changes the principal symbol.** The principal speeds are the six roots of Q_new(c; 1) = F_eff − c²d1d2 − 4γ_u q B(c)/M⁴.
- B(c) has terms odd in γ_u, so **the sign of γ_u enters the principal symbol**.
- A0 is eliminated exactly through the Gauss row, so this formulation has no constraint-violating directions.
- **Not derived:** whether the registered 6×6 first-order system keeps the pure k P_s + k̂ P_k structure. The determinant shows that the dispersion relation keeps the registered form with Q → Q_new.
- **Strong hyperbolicity (six simple real roots of Q_new(c; 1)) cannot be decided analytically in general. It depends on the background.**

## Findings not anticipated by the specification

1. **The sign of γ_u matters beyond leading order.** §1 fixed only |γ_u| through Λ_eff. The sign is irrelevant at O(M⁻²) (D0, D1), but it enters the O(M⁻⁴) mixed term, K_rem (hence N_i, k_der and H3), and the principal symbol (B(c)). The relative sign of γ_u and g_φ is therefore a model parameter. Both signs must be evaluated, and the result reported for each. No sign is chosen here.
2. **K_rem is rank one in 1+1D** (e eᵀ structure). This is what makes k_der available in closed form.

## Classification (§3 of the specification)

1. **NOT-VIABLE:**
   - (a) does not apply: Λ_eff > 0 (proven).
   - (b) does not apply: H1 proven.
   - (c) does not apply unless every background is R-UNDEFINED; k_der > 0 wherever defined.
   - (d) is undetermined: hyperbolicity of Q_new(c; 1) is background-dependent.
2. **VIABLE-FOR-NUMERICS:**
   - V1 holds;
   - V4 holds;
   - V2 is undetermined: it needs the Stage 2 arithmetic of max_i(N_i + c_i²) < 2;
   - V3 is undetermined: it needs the roots of Q_new(c; 1) per background.
3. **Outcome at Stage 1: OPEN.**

## What is needed to go further (not done; needs authorisation)

- **Stage 2, as registered:** on the 41 backgrounds, for both signs of γ_u:
  - (i) the D4 condition;
  - (ii) k_der from the D3 formula, with R-UNDEFINED, no-crossing cases, poles and H3 reported separately.
- **Not covered by the registered Stage 2 text:** the roots of Q_new(c; 1) on the 41 backgrounds, i.e. V3 and NOT-VIABLE (d). This is arithmetic of a Stage 1 formula. It needs an explicit decision: either an extension of the Stage 2 authorisation, or a separate Stage 3 specification.

## Not changed

B4.12, T7a, the thresholds A/B/C and the registered results are unchanged. Stage 2 and Stage 3 were not started. Nothing here is evidence for B4.
