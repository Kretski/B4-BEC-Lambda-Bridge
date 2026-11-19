# Symbolic check of the Stage 2 procedural annex, items A1 and A3 — report v0.1

**Status:** report for review.

**Authorised scope** (review of 2026-11-19): a symbolic check of A1 and A3 of `STAGE2_PREREGISTRATION_PROCEDURAL_ANNEX_v0.1.md` (ad965017…) only. Not authorised: a numerical Stage 2 run, a background search, a fit, any change to B4, registration or commit.

**Script:** `stage2_annex_A1_A3_symbolic_check.py`, sha256 8fdb07e864e523c4db44b025eac6703fc2ca330907de10e80500b107c4443141. It uses SymPy only, reads no input files and uses no randomness.

**Run:** one run with `python3 -I`, Python 3.13.18, SymPy 1.14.0, rc 0.
- stdout sha256: 8ba52ccc06d42540744d2afbcae8140403b1a224010369c3796418f4966a5faf
- RUN_CONTEXT sha256: 199b1660…

**Result:** **15/15 PASS.**

**Note on procedure.** The reviewer authorised the check directly. The script was not sent for static review before the run, so the run is reported together with the script.

## 1. What is proven (symbolically, by the script)

| ID | Statement | Result |
|---|---|---|
| A1.0 | For a k ≠ 0 Fourier mode, the first-order term of E − Σμ_iN_i averages to zero over one period | PASS |
| A1.0′ | Uniform stationarity gives μ_i = g_iin_i + g₁₂n_j + m_iV_i²/2, for i = 1, 2 | PASS |
| A1.1 | The second-order period average of the hydrodynamic GP functional (quantum pressure ħ²(n′)²/(8mn), kinetic (m/2)nu², g_iin²/2, g₁₂n₁n₂, −μn), for a real mode with cos and sin parts, equals ½[E₂(cos) + E₂(sin)], with E₂ exactly as in annex A1 | PASS |
| A1.2–A1.3 | The matrix M of E₂ in z = (δn₁, δu₁, δn₂, δu₂) is as stated, and symmetric | PASS |
| A1.4 | D·M·D, with D = diag(g₁₁^{−½}, (m₁n₁)^{−½}, g₂₂^{−½}, (m₂n₂)^{−½}), equals the annex matrix Ẽ(k) written in a_i = g_iin_i/m_i, b_i = ħ²/(4m_i²), f = g₁₂/√(g₁₁g₂₂) and V_i | PASS |
| A1.5 | D is invertible for g_ii, m_i, n_i > 0 | PASS |
| A1.6 | a_i > 0 and b_i > 0 for g_ii, n_i, m_i, ħ > 0 | PASS |
| A3.1 | Ẽ(k) − Ẽ(0) = diag(b₁k²/a₁, 0, b₂k²/a₂, 0) | PASS |
| A3.2 | This difference is positive semidefinite | PASS |
| A3.3 | For 0 ≤ k₁ < k₂: Ẽ(k₂) − Ẽ(k₁) = diag(b_i(k₂² − k₁²)/a_i on the δn entries, 0 elsewhere), which is positive semidefinite | PASS |
| A3.4 | One component at k → 0: the 2×2 block has det = 1 − V²/a and trace = 2 | PASS |
| A3.5 | One component: positive definite ⇔ V² < a, i.e. −√a < V < √a | PASS |
| A3.6 | The component-1 block of Ẽ(0) at f = 0 equals that one-component matrix | PASS |

## 2. What is cited, not proven by the script

- **Monotonicity of λ_min.** If P ⪰ 0, then λ_min(A + P) ≥ λ_min(A). This is the Courant–Fischer theorem (Weyl's inequality), cited, not re-derived. With A3.3 it gives: **λ_min(Ẽ(k)) is non-decreasing in k²**, and over any k grid that contains k = 0, the minimum is attained at k = 0.
- **2×2 positive definiteness.** A symmetric 2×2 matrix is positive definite iff its trace and determinant are positive. This standard fact is used in A3.5.

## 3. Assumptions and scope limits

- **Physical positivity.** g_ii > 0, n_i > 0, m_i > 0, ħ > 0. The registered `bec_dispersion.BECParams.validate` enforces g_ii > 0 and m_i, n_i > 0. g₁₂ and V_i are arbitrary real numbers.
- **Mean field, 1D, lab frame.** No LHY term and no trap. The energy is E − Σμ_iN_i in the lab frame, with uniform background flows V_i. Fluctuations are treated as **independent** δn_i and δu_i = (ħ/m_i)∂_xδθ_i. This is legitimate for a k ≠ 0 mode, where any δu amplitude can be reached.
- **Meaning of "k = 0".** Ẽ(0) is the limit k → 0⁺ of the k ≠ 0 form. It does not describe a uniform (k = 0) change of the flow, which would change the winding or current and is not part of this diagnostic. The grid point "k = 0" in A3 of the annex is to be read in this sense.
- **Frame dependence.** Positive definiteness of Ẽ depends on V₁ and V₂ separately, i.e. on the lab frame (annex A2). The check does not decide the choice of frame.
- **Sylvester's law of inertia.** The congruence preserves the **inertia**: the numbers of positive, negative and zero eigenvalues. So positive definiteness of M and of Ẽ coincide. **It does not preserve eigenvalue magnitudes.** The v2 thresholds −10⁻⁹ and −10⁻⁶ therefore do not carry over between M and Ẽ. As the reviewer requires, λ_min(Ẽ) must be written in the annex as an **explicitly chosen operational quantity** (a choice of normalisation), not as a consequence of Sylvester's law. Annex v0.1 §A1 currently suggests the latter. This is to be corrected in v0.2.

## 4. Consequence for the annex

- **A1 is correct** as a derivation of the normalised 4×4 form from the stated functional, under the assumptions of §3.
- **A3 is correct.** The k grid cannot lower λ_min below its k = 0 value. Energetic stability in this diagnostic reduces to positive definiteness of Ẽ(0), a two-component Landau-type condition. The one-component limit is the Landau criterion V² < a.
- **Grid in an implementation.** It stays only as an implementation check (ABORT if λ_min(k) < λ_min(0) − 10⁻¹²), as already written in the annex.

## 5. Not done

- No numerical Stage 2 work, no search, no fit, no change to B4, no registration, no commit.
- The other two clarifications required before registration go into annex v0.2, after this report is reviewed:
  - **B2:** how starts are generated and accepted, and how the solver stays in the admissible region under hard rejection;
  - **C1:** rejection of perturbed targets with ps < 0, and the rule for κ = 0.
- Also for annex v0.2: the A1 wording on thresholds (§3 above).

## 6. Files

- `check/stage2_annex_A1_A3_symbolic_check.py` (8fdb07e8…)
- `run_out/stdout.txt` (8ba52ccc…)
- `run_out/RUN_CONTEXT.txt` (199b1660…)
- `run_out/SHA256SUMS.txt`
