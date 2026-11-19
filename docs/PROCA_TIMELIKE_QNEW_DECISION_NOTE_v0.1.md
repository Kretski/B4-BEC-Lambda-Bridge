# Proca with a timelike current — NOT-VIABLE (d): decision note v0.1

**Specification:** `docs/PROCA_TIMELIKE_CURRENT_ANALYTIC_SPEC_v0.5.md` (87fc31b), unchanged.

**Base:** 53e2227 (Stage 1 report, D5 addendum v0.3, Stage 2 results).

**Decision (reviewer, 2026-11-19 05:48 EET):** **NOT-VIABLE by criterion (d), for both signs of γ_u.** It applies within the registered linear model, the wording of the criterion and the checked backgrounds.

**Evidence level:** 1 (exact and certified arithmetic of the registered principal symbol; no new physics input). Not evidence for or against B4.

## 1. What was run (each run authorised once, after static review)

| Run | Script (SHA-256) | Scope | Result files |
|---|---|---|---|
| Q_new(c; 1) roots | `proca_timelike_qnew_roots.py` `91ea9152…` | 41 registered E3 backgrounds × σ = ±1 | `results_qnew/` (JSON `ff5f724c…`) |
| Current attribution | `proca_timelike_qnew_current_attribution.py` `21c0934f…` | 4 backgrounds × σ = ±1 | `results_attribution/` (JSON `b044acaf…`) |

Both were run from a fresh, clean LF clone. The Q_new run used 53e2227. The attribution run used 53e2227 together with the reviewed Q_new script. All input hashes were verified before each start. Conventions as in Stage 2: M = 1, g_φ = +√p, g_ψ = −κ/√p, γ_u = σ√(Λ/2).

**Revisions before the runs (all statically reviewed):**

- Q_new script:
  - `2c475aba` → `f34ea02d`: exact check of the leading coefficient; relative threshold made primary (reviewer decision); optimal 720-permutation matching in the P_∞ control;
  - `f34ea02d` → `91ea9152`: exact test for a zero root via the minimal polynomial.
- Attribution script:
  - `de4aa952` → `c599a80e`: a failed identity check (A), or a disagreement with the registered flag, gives NOT ESTABLISHED; the homotopy is labelled as observed transitions;
  - `c599a80e` → `21c0934f`: a bisection that does not converge is reported separately, not as a localised transition.

Synthetic smoke tests: 14/14 for the Q_new script and 9/9 for the attribution script.

## 2. Q_new(c; 1) on the 41 backgrounds

Q_new is taken verbatim from the registered D5 checker (13/13). Its coefficients are computed exactly in Q(√p, √(Λ/2)).

**Primary rule.** T3_FUND_PROCA_SPEC §3(iv), "tolerance 1e-9 relative", read as |x − y| / max(1, |y|) at 50 digits. This is the reading used by every other relative 1e-9 test of the registered T3-fund code.

**Diagnostics.** Absolute 1e-9 (the registered code line for this rule) and pure relative.

**Certified cross-check.** Exact norm polynomial of degree 24, interval arithmetic.

| | σ = +1 | σ = −1 |
|---|---|---|
| SIX-SIMPLE-REAL | 25 | 26 |
| COMPLEX-ROOT (always exactly one complex pair) | 16 | 15 |
| REPEATED / UNDETERMINED | 0 / 0 | 0 / 0 |

- The three readings of the threshold and the certified check agree on every background.
- Smallest margins: |Im c| / max(1, |c|) ≥ 0.034 on COMPLEX; relative root distance ≥ 0.040 on SIX.
- Control: the eigenvalues of P_∞ equal −i c_j to 7.5·10⁻⁴⁹.

**Read-only comparison with the registered T3-fund Proca result** (34953d0; vector mediator, bare Λ, no current):

- 13 of the 41 are already not strictly hyperbolic there (12 negative controls and 1 boundary). All 13 are COMPLEX here for both signs.
- Strictly hyperbolic there but COMPLEX here: 3 backgrounds for σ = +1 and 2 for σ = −1.
- No background changes the other way.

## 3. Attribution: is the complex pair caused by the current?

**Decomposition** (exact, from the registered checker):

Q_new(c; 1) = Q0(c) − 4γ_u q B(c),  with Q0 = F_eff − c²d1d2 and Λ_eff = 2γ_u² = Λ.

Q0 is the principal-symbol polynomial of the same reduced system without the exchange term of the current.

The literal limit γ_u → 0 is not used for the attribution, because Λ_eff would vanish too: that is a different theory. It is reported only as diagnostic (D).

| Background | σ | (A) Q0 vs registered `Qcoeffs` | Q0 primary/certified | Q_new primary/certified | Observed transition 6 → 4 at t* | Attribution |
|---|---|---|---|---|---|---|
| rank1_minus_stable_5713 | +1 | 3.7e-16 | SIX/SIX | COMPLEX/COMPLEX | 0.0034 | **CURRENT-INDUCED** |
| rank1_minus_stable_5713 | −1 | 3.7e-16 | SIX/SIX | SIX/SIX | — | NOT-INDUCED |
| rank1_minus_boundary_2020 | +1 | 3.2e-16 | SIX/SIX | SIX/SIX | — | NOT-INDUCED |
| rank1_minus_boundary_2020 | −1 | 3.2e-16 | SIX/SIX | COMPLEX/COMPLEX | 0.9066 | **CURRENT-INDUCED** |
| rank1_plus_boundary_329 | +1 | 5.1e-16 | SIX/SIX | COMPLEX/COMPLEX | 0.3862 | **CURRENT-INDUCED** |
| rank1_plus_boundary_329 | −1 | 5.1e-16 | SIX/SIX | COMPLEX/COMPLEX | 0.7228 | **CURRENT-INDUCED** |
| rank1_plus_boundary_3437 | +1 | 6.3e-16 | SIX/SIX | COMPLEX/COMPLEX | 0.5239 | **CURRENT-INDUCED** |
| rank1_plus_boundary_3437 | −1 | 6.3e-16 | SIX/SIX | SIX/SIX | — | NOT-INDUCED |

- On all four backgrounds the registered T3-fund flag `strictly_hyperbolic` is True, and it agrees with the class of Q0.
- Smallest relative root distance of Q0: 0.050. Smallest |Im c| / max(1, |c|) of Q_new where COMPLEX: 0.083.
- **Homotopy Q_t = Q0 − t·4γ_u qB.** These are the observed transitions along the discretely tracked homotopy (401 points of t, bisection to 1e-15). Each COMPLEX case shows exactly one transition in which one pair of real roots leaves the real axis. No grid point or bisection failed to converge. This is a diagnostic; a proof that all transitions were found would need the discriminant in t.
- For 5713 (σ = +1) the transition comes early (t* = 0.0034), next to the closest pair of roots of Q0. The endpoint is clearly COMPLEX (0.342).

## 4. Decision

**Interpretation of (d).** The reviewer's decision of 2026-11-19 05:10 was taken after the Q_new run and before the attribution run.

- (d) is existential. Unlike (c), it carries no "on all applicable backgrounds". One applicable background on which **the current** makes the physical principal symbol non-hyperbolic is sufficient.
- Backgrounds that are not hyperbolic already without the current do not count.

**Physical symbol.** D5 addendum (D5e, G1, S1) shows that the Gauss constraint holds identically and that the six reduced variables are free Cauchy data. So no eigenvalue belongs to constraint-violating directions, and the reduced principal symbol is the physical one.

**Result.** There are 5 CURRENT-INDUCED cases: 3 for σ = +1 (5713, 329, 3437) and 2 for σ = −1 (2020, 329). Hence:

- **NOT-VIABLE (d) for σ = +1 and for σ = −1.**
- By the decision rule of v0.5 §3, applied in order, NOT-VIABLE is the outcome. The Stage 2 V2 result (σ = +1: 0/41; σ = −1: 2/41) is not needed for this class, and it is unchanged.

**Scope.** The decision applies to the registered linear model about homogeneous backgrounds in 1+1D, the registered wording of (d) and the checked backgrounds. It does not cover other configurations of the theory, and it does not claim a nonlinear loss of well-posedness.

**Not done, and not to be done to rescue the result:** Stage 3, population runs, V3 PASS/FAIL, κ(V_∞).

**What remains open:**

- T7b: a UV-derived Λ coefficient, not inserted by hand. This Proca configuration with a timelike current is closed as a negative result.
- EXT-1, T7a, Script B/C4 and B4.12 are unchanged.
- The specification and all thresholds are unchanged.

**Scientific content.** A concrete mechanism was identified: the exchange term 4γ_u q B(c) of the timelike current moves one pair of characteristic speeds off the real axis on backgrounds that are strictly hyperbolic without it. This is a structural loss of hyperbolicity, not a numerical failure.

## 5. Files

`experiments/uv_candidates/proca_timelike/qnew/`:

- `proca_timelike_qnew_roots.py` `91ea9152…`; smoke test `smoke_qnew.py` `dcbebce6…` with its output `47f1c2ef…` (14/14); `revision_diffs/`;
- `proca_timelike_qnew_current_attribution.py` `21c0934f…`; smoke test `smoke_attribution.py` `19865f06…` with its output `85aeb652…` (9/9); `revision_diffs/`;
- `results_qnew/`: JSON, `run_log`, `RUN_CONTEXT`, `SHA256SUMS` (`575ca1e7…`), and `posthoc/` (read-only comparison with the registered T3-fund result);
- `results_attribution/`: JSON, `run_log`, `RUN_CONTEXT`, `SHA256SUMS` (`d6f88742…`).
