# Proca with a timelike current — Stage 1 D5 addendum v0.3

**Specification:** `docs/PROCA_TIMELIKE_CURRENT_ANALYTIC_SPEC_v0.5.md` (87fc31b), unchanged.

**Status:** addendum to the Stage 1 report v0.1, for review. Not registered.

**Evidence level:** 1.

**What it covers:** only D5, plus the four checks of the Stage 1 script that were recorded as an unconditional `True`. D0–D4 are not changed.

**Symbolic checks:** `proca_timelike_d5_addendum_check.py` (SHA-256 `2064fa84…`). Output `d5_addendum_output.txt` (`80d8201b…`): **13 checks, 0 failed**. Pure algebra; no backgrounds and no thresholds.

**Changes from v0.2:** only the checker. W2 now requires every limit of M̃(k) − kP_∞ to be exactly 0 (`rem == zeros(6, 6)`), matching the claim B(k) → 0; before, it required only finite limits. The text and all results are unchanged.

**Changes from v0.1:** §6 rewritten in answer to the review (pseudodifferential reduction, weighted norm, uniform control of the transformation, conditioning of the diagonaliser). Two computed checks added: W3 (sup over k ≥ 0 of the lower-order part is finite) and W4 (reality: k < 0 reduces to k > 0). §1–§5, §7 and the results of the eleven v0.1 checks are unchanged.

During development, G2 first failed. The check had tried to solve the Gauss equation for A0 while A0 still appeared inside F. It now substitutes A0 from the constraint and A1_t = F + A0_x directly. No mathematics was changed.

## 1. Reduced first-order system (registered variables)

The variables are U = (a, b, c̃, d, A1, F) = (f_t, f_x, g_t, g_x, A1, F), as in T3_FUND_PROCA_SPEC §3(iii). They obey:

- **Gauss constraint:** A0 = −(F_x + J_u)/M².
- **Definition of F:** A1_t = F + A0_x = F − (F_xx + ∂_xJ_u)/M².
- **A1 equation:** F_t = J_n − M²A1.
- **Compatibility:** b_t = a_x and d_t = c̃_x.
- **f and g equations:** a_t and c̃_t come from the f and g Euler–Lagrange equations. The f equation contains the background Proca field, A1_t, and ∂_x(2γ_u A0 b).

**G1.** The Gauss constraint holds identically for every state, because A0 is defined by it and is not an evolved variable.

**G2, G3.** The reduced f, g and F equations equal the Euler–Lagrange equations of the full action. This is proven symbolically for arbitrary fields, not only for the linearisation.

## 2. Spectrum = Lagrangian spectrum, for all k (S1)

The linearisation about the background gives dU/dt = M_lin(k)U. Then

det(M_lin(k) + iω I) = −det K(ω, k)/M²   (exact).

Here K is the 4×4 Lagrangian kernel of Stage 1. So the six eigenvalues −iω of M_lin(k) are exactly the six roots of the dispersion relation of the full constrained Lagrangian system, at every k. The reduction neither loses nor adds modes.

## 3. Principal symbol under the registered scaling (W1, W2)

**Scaling.** With the registered scaling Ã = M δA1, F̃ = (k̂/M) δF, k̂ = √(M² + k²), every entry of the scaled matrix M̃(k) is O(k). The limit

P_∞ = lim_{k→∞} M̃(k)/k

exists, and M̃(k) − k P_∞ → 0 as k → ∞ (W1, W2). **No Douglis–Nirenberg weights beyond the registered scaling are needed.**

The timelike current fits the registered structure. It enters P_∞ through three entries:

- (a, b): (M² + 2g_φγ_u q² + 6γ_u²q²)/M²;
- (c̃, b): 2g_ψγ_u q r/M²;
- (Ã, b): −2γ_u q/M, which comes from ∂_xJ_u in the A1 equation, as anticipated in the specification.

It also changes the (a, F̃) entry to −q(g_φ + 2γ_u)/M. The full P_∞ is printed in the output file.

## 4. Spectrum of P_∞ = roots of Q_new(c; 1) (P1)

det(λI − P_∞) evaluated at λ = −ic equals Q_new(c; 1) exactly: the proportionality constant is 1.

So **the eigenvalues of P_∞ are −i c_j, where c_j are the six roots of Q_new(c; 1) = F_eff − c²d1d2 − 4γ_u q B(c)/M⁴.**

## 5. Physical and non-physical directions (D5e, now computed)

- The six reduced variables are free Cauchy data:
  - f and g are recovered from (a, b) and (c̃, d) up to constants, and in 1+1D there is no integrability constraint;
  - (A1, F) are the single Proca degree of freedom.
- A0 is not evolved, and the Gauss constraint holds identically (G1).
- The spectrum equals the Lagrangian spectrum (S1).

Together these show: **there are no constraint-violating directions, and all six directions are physical.** In this formulation, hyperbolicity of the reduced system and physical hyperbolicity are the same thing.

## 6. Strong hyperbolicity of the linearised system: norm, transformation, diagonaliser

**Scope.** This section concerns only the linearisation about one homogeneous background, in 1+1D. Its coefficients are constant, so the Fourier transform in x decouples it exactly into one 6×6 ODE system per wavenumber k: dÛ/dt = M_lin(k)Û. No symbol calculus and no commutator remainders occur. Nonlinear or variable-coefficient well-posedness is **not** claimed.

**(1) The reduction is pseudodifferential, not local.** The scaling uses k̂ = √(M² + k²), i.e. the Fourier multiplier (M² − ∂_x²)^{1/2}. In the scaled variables Ũ = S U the system is dŨ/dt = M̃(k)Ũ with M̃(k) = S(k)M_lin(k)S(k)⁻¹, S(k) = diag(1, 1, 1, 1, M, k̂/M). M̃(k) is a Fourier multiplier of order one, not the symbol of a local first-order differential operator. Every statement below is therefore made in Fourier space, for each k separately, and then transferred to x-space by Parseval. This is exact for constant coefficients.

**(2) The weighted norm.** For U = (a, b, c̃, d, A1, F) define

‖U‖²_E = ‖a‖² + ‖b‖² + ‖c̃‖² + ‖d‖² + M²‖A1‖² + M⁻²‖(M² − ∂_x²)^{1/2} F‖²,

with ‖·‖ the L² norm in x. By Parseval, ‖U‖_E = ‖S(k)Û(k)‖_{L²(dk)}. So the derivatives f_t, f_x, g_t, g_x and the Proca potential A1 are in L², and F is in H¹. This is the natural energy space of the Proca field: with the Gauss constraint A0 = −F_x/M² (free part), the free Proca energy density ½F² + ½M²A0² + ½M²A1² equals ½(M²A1² + M⁻²(M²F² + F_x²)), i.e. exactly the last two terms. (Elementary algebra; stated, not a new check.)

**(3) The transformation is uniformly controlled.** By the definition in (2), S is an isometry from E onto L²(ℝ)⁶, and S⁻¹ is its inverse isometry. Both have operator norm 1, for all k at once, with no k-dependent constant. The matrix S(k) itself grows like k; that growth is absorbed exactly by choosing E as above, which is the whole purpose of the scaling. E is the usual space L²⁵ × H¹ with the M-dependent equivalent norm.

**(4) The diagonaliser has bounded conditioning.** Write M̃(k) = kP_∞ + B(k). Here:

- P_∞ = lim_{k→∞} M̃(k)/k does not depend on k (W1). Its eigenvalues are −ic_j, where c_j are the six roots of Q_new(c; 1) (P1).
- B(k) = M̃(k) − kP_∞ is continuous on k ≥ 0 (its denominators are only powers of M and of k̂ ≥ M, W3) and tends to 0 as k → ∞ (W2). Hence β := sup_{k≥0} ‖B(k)‖ < ∞.
- M̃(−k) = conj M̃(k) (W4), so ‖e^{tM̃(−k)}‖ = ‖e^{tM̃(k)}‖ and k ≥ 0 is enough.

**If Q_new(c; 1) has six simple real roots,** then P_∞ = V D V⁻¹ with D = −i diag(c_1, …, c_6) purely imaginary and V a **fixed** matrix, independent of k, with finite condition number κ(V). Only the principal part kP_∞ is diagonalised; the eigenvectors of the full M̃(k) are not used, so their possible degeneration at finite k does not matter. Then

V⁻¹M̃(k)V = kD + V⁻¹B(k)V.

kD is skew-Hermitian, so the logarithmic norm satisfies μ(V⁻¹M̃V) ≤ ‖V⁻¹BV‖ ≤ κ(V)β. Hence, for all k and all t ≥ 0,

‖e^{tM̃(k)}‖ ≤ κ(V) e^{κ(V)β t},

and by Parseval ‖U(t)‖_E ≤ κ(V) e^{κ(V)β t} ‖U(0)‖_E. This is well-posedness of the linearised problem in E, i.e. strong hyperbolicity in the sense of a uniformly bounded symmetriser (V⁻¹)*V⁻¹. It is the registered T3_FUND_PROCA_SPEC §3(iv) argument γ(k) ≤ μ(k) with the fixed V_∞, with Q replaced by Q_new and the scaling and norm stated explicitly.

**Other cases (unchanged from v0.1 and from the registered rule):**

- **Any complex root:** P_∞ has an eigenvalue with non-zero real part, the growth rate is unbounded in k, and the problem is ill posed. This is NOT-VIABLE (d), FAIL-UNSTABLE in the registered sense.
- **A repeated real root:** INCONCLUSIVE. Diagonalisability with bounded κ(V) is not implied by the root count. The rank of (P_∞ + icI) could decide it, but that is not part of the registered rule and is not proposed.
- The registered tolerance for "six simple real roots" is 1e-9 relative at 50 digits. It is used unchanged.

**Note for the root arithmetic.** The registered §3(iv) records κ(V_∞) per background. The authorised Q_new arithmetic lists roots, counts, simplicity, minimum distance and class. Recording κ(V_∞) as well would complete the bound above. It is proposed, not added without approval.

## 7. The four Stage 1 checks recorded as unconditional `True` — now computed

| Check | Before | Now |
|---|---|---|
| D2c | `True` with a comment | L4 with J^0 = 0 reduces to (∂_tJ_n)²/2 (·M⁻⁴): computed |
| D3a | `True` (structural) | with the background Proca values as free symbols, K_LA and K_AA do not depend on them; only K_LL does. So the background term cancels in K_rem: computed |
| D5b | `True` | Q_new(c; ε = 0) = F_eff(c): computed |
| D5e | `True` | G1 + S1 + the six free Cauchy variables: computed |

## 8. Status

**D5 is complete analytically.** The principal symbol, the scaling, the norm, the Gauss constraint, the physical directions, and the relation between the spectrum and Q_new(c; 1) are all derived. Six simple real roots of Q_new(c; 1) imply the uniform bound of §6 in the norm E. Whether Q_new(c; 1) has six simple real roots on the 41 backgrounds, for each sign of γ_u, is arithmetic of this formula. It needs the authorisation the review reserved.

The class remains **OPEN**.
