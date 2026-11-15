# B4.12 Q2 — Specification of the analytic derivation of H-Q2 and ε_EFT (proposal v0.2)

**Status: PROPOSAL FOR REVIEW. Nothing has been derived or computed.** This specification fixes, before any work, what is derived, by which method, which checks must pass, and what is out of scope. It implements Q2 of the registered B4.12 v0.4 (`B4_12_PROCA_PERTURBATIVITY_AND_UV_CLOSURE_v0.4.md`, commit 74c2a77; registration record 1e0ad68).

It changes nothing registered:
- the definitions of v0.4 §3 (quantities, error metrics, reference, regime) are used unchanged;
- no grid, threshold, τ or ε_max is fixed here.

Q3 is a later, separate specification.

**Changes from v0.1** (review of 2026-11-15, 15:28, REVISION REQUIRED):
1. D2 and Q2-C4: the light branches come from the reduced action. The heavy Proca branch is checked in the full quadratic system before a_μ is eliminated. Roots of the reduced relation at the pole of K are spurious until confirmed in the full system.
2. Q2-C3: the comparison is with the rank-one subfamily of B4. The non-coincidence with general B4 (independent p, s, κ) is stated explicitly.
3. D3: ε_tot is undefined when ω_exact,b² = 0. ε_abs is reported instead, and no regulator is introduced.

---

## 1. Setting (from v0.4 and T7a step 3; unchanged)

- **Action.** The registered T3-fund Proca action (`T3_FUND_PROCA_SPEC.md` §1, e355af2), in 1+1D, metric (−,+), with D = ∂_t² − ∂_x².
- **Couplings** (Proca matching): g_φ = M√p and g_ψ = −Mκ/√p. Hence g_ψ²/M² = κ²/p, which equals s only under rank one.
- **Exact tree-level kernel:** K(ω,k) = (1 + k²/M²)/(M² + k² − ω²) (T7a step 3, 61fefff).
- **Background:** constant gradients (w̄, q̄, v̄, r̄), i.e. **ν = 0**. This is the primary and only case of this specification. ν > 0 is a separate later case.
- **Quantities and error metrics:** as in v0.4 §3.1 and §3.3:
  - ε_EFT is primary;
  - ε_tot and ε_abs are secondary;
  - ε_op is diagnostic only.

## 2. What is derived

**D1 — exact quadratic fluctuation action (test of H-Q2).**
- Expand the fields, f = f̄ + δf, g = ḡ + δg and A = Ā + a, and derive the quadratic action in (δf, δg) after integrating out a exactly at tree level (A equation and Gauss constraint, as in step 3).
- Write J = J̄ + J₁ + J₂ (order 0, 1, 2 in the fluctuations).
- Establish whether, for ν = 0, the result is exactly

  S₂ = S₂^light + ½·J₁ K J₁ + J̄·J₂/M²,

  where S₂^light is the quadratic part of the light action, including the Λ term.
- The derivation must state explicitly: the background solution (Ā0, Ā1) for constant gradients; why J̄·K·J₂ reduces to J̄·J₂·K(0,0) = J̄·J₂/M² (zero mode, self-adjointness of K); and every assumption on boundary terms.
- If H-Q2 fails, the derived form replaces it and the failure is reported. It is not repaired.

**D2 — dispersion relations.** In Fourier space (δf, δg ∝ e^{i(kx − ωt)}).

**D2(a) — light branches, from the reduced action of D1.** The 2×2 operator in (δf, δg) gives three dispersion relations:
- (i) **exact**, with the full K;
- (ii) **M⁻² truncation**, K → 1/M², the B4 EFT dispersion of the rank-one subfamily (see Q2-C3);
- (iii) **M⁻⁴ truncation**, K → (1 + ω²/M²)/M², reported for reference.

There are four light branches (ω = ±… for each of two modes). Branch identification follows a rule fixed here: a light branch is the root that tends continuously to the corresponding M⁻² root as M → ∞ at fixed k. Tracking at finite k/M follows that continuation analytically. No numerical nearest-root assignment is used in this step.

**D2(b) — heavy Proca branch, from the full quadratic system.** The heavy branch is **not** obtained from the reduced relation. When A_μ is eliminated, its own dynamics are no longer contained in the effective action for (δf, δg). The heavy branch is therefore derived from the full quadratic system in (δf, δg, a0, a1) **before** the elimination of a_μ.

**Spurious roots.** Multiplying the reduced relation of D2(a)(i) by the denominator of K can introduce roots at the pole of K, M² + k² − ω² = 0. Such roots are not treated as physical. Each is confirmed or rejected in the full system of D2(b), and the result is reported.

**D3 — ε_EFT and the secondary metrics, per light branch b.** With Δ_X(k) = ω_X,b(k)² − k² for X ∈ {exact, 2}:
- ε_EFT,b = |Δ_2 − Δ_exact| / |Δ_exact|, the relative error of the B4 effect itself (v0.4 §3.3);
- ε_tot,b = |ω_2,b² − ω_exact,b²| / ω_exact,b²; where ω_exact,b² = 0 (e.g. at k = 0 on a massless branch), ε_tot,b is reported as **undefined**, ε_abs,b is reported instead, and no regulator or small-denominator substitute is introduced;
- ε_abs,b = |ω_2,b² − ω_exact,b²| / M².

Deliverables:
- (i) the closed form where one exists;
- (ii) the leading terms of the expansion in k/M, with explicit coefficients in terms of the background (p, κ, Λ, w̄, q̄, v̄, r̄) and of the structure θ;
- (iii) the **explicit dependence on ρ_j** (v0.4 §0), stated as derived: whether ε_EFT depends on ρ at leading order, and if so how.

For any branch where Δ_exact → 0 (no B4 effect, e.g. J₁ ≡ 0 for that mode), ε_EFT is undefined. The derivation states these cases explicitly; only ε_abs is reported for them (v0.3 §3 Q2 item 2c). τ is not fixed here.

**D4 — physically admissible regime.** Derive analytic conditions, in terms of k/M and the background, for:
- (a) real ω for the light branches;
- (b) no crossing of the K pole, i.e. 1 + (k² − ω²)/M² > 0 along each light branch;
- (c) separation of the light branches from the heavy ones (no avoided crossing below the stated k/M);
- (d) consistency with the T3-fund C3 result for the same action (reality and hyperbolicity), cited, not recomputed.

Where a condition cannot be put in closed form, this is stated, and the condition is left for a later numerical step with its own specification.

**D5 — outcome classification of v0.4 §3.4.** From D3 and D4, state which of the following holds over the admissible regime, with the derivation that supports it:
- (i) ε_EFT depends on ρ alone;
- (ii) it also depends on k/M, so a ρ threshold is not a standalone criterion;
- (iii) not determinable analytically.

No ε_max is used; this is a structural classification only.

## 3. Method

- **Tool:** symbolic derivation with SymPy, in a new standalone script `b412_q2_derivation.py`. Like T7a step 3, it derives the equations from the full action, not from a pre-reduced system. The arbiter is the registered action and its conventions; SymPy is a tool.
- **No floating-point evaluation** on any registered background (the three C4 backgrounds, the 41 E3 backgrounds or the populations). Numerical evaluation of the derived formulas on those backgrounds is a later step that needs its own specification.
- Symbolic checks may use exact rational test values (e.g. p = 1/4) only to verify identities. This is not an evaluation.

## 4. Checks (all mandatory; a derivation without them is INCONCLUSIVE)

| # | Check | PASS condition |
|---|---|---|
| Q2-C1 | Kernel | the kernel derived from the action equals K(ω,k) of T7a step 3, identically |
| Q2-C2 | H-Q2 identity | S₂(derived) − [S₂^light + ½J₁KJ₁ + J̄J₂/M²] ≡ 0 symbolically, or the explicit non-zero remainder is reported (then H-Q2 FAILS, as a reported result) |
| Q2-C3 | EFT limit | (1) the Proca M⁻² truncated dispersion equals the B4 EFT dispersion of the **rank-one subfamily** (κ² = ps). The Proca coefficient of (g_tg_x)² is κ²/(2p) (from g_ψ²/M² = κ²/p), which equals s/2 only under κ² = ps. (2) For general B4 with independent p, s, κ, the report states explicitly that the coincidence does **not** hold: it is not a general property of the single-mediator Proca construction |
| Q2-C4 | Decoupling (full system) | in the **full quadratic system** of D2(b), for g_φ = g_ψ = 0, the system decouples into the two free light branches (the free dispersion of S₂^light, including the Λ term) and the free Proca branch |
| Q2-C5 | Two routes for the series | the k/M expansion of ε_EFT from (a) the series of the exact root and (b) perturbation theory in 1/M² on the dispersion polynomial agree term by term to the stated order |
| Q2-C6 | Consistency with M²K − 1 | on the light shell, the leading correction factor equals (ω²/M²)/(1 + (k² − ω²)/M²) (v0.4 §3.2), identically |

The verdict `ALL_CHECKS_PASS` covers only these checks.

## 5. Outputs

- `b412_q2_derivation.py` (standalone; SymPy).
- `b412_q2_output.txt`: each check with ZERO / non-zero, and the final verdict.
- `b412_q2_results.json`: the derived expressions (`sympy.srepr` and LaTeX) for D1–D4, the series coefficients of D3, and the classification of D5.
- `B4_12_Q2_DERIVATION_REPORT.md`: the derivation in readable form, the D5 classification and its scope.

All outputs are in English.

## 6. Out of scope

- ν > 0 (inhomogeneous backgrounds).
- Any numerical grid, τ, ε_max or threshold.
- Evaluation on registered backgrounds.
- Q3 (unitarity) and the joint Q2/Q3 check.
- Changes to registered documents, to the evidence, or to the T7a results.
- Any statement on whether Proca is a valid UV completion of B4.
- No Zenodo.

## 7. Order

1. Review of this specification.
2. Commit.
3. Script, with a review of the script before any run, as for the common-M check.
4. Run, then the report, review and commit.
5. Then the Q3 specification.
