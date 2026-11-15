# B4.12 — Common-M test: result and interpretation (draft v0.1)

**Status: RESULT DOCUMENT (draft v0.1, including the edits of the review of 2026-11-15, 15:07).** The computation was accepted by the reviewer on 2026-11-15, 14:58. The interpretation (§3) and the decisions on W1–W4 (§4) were approved on 2026-11-15, 15:07. This document changes no approved or registered file. B4.12 remains OPEN until the documentary integration and the review of the B4.12 proposal v0.4.

## 1. What was run

| Item | Value |
|---|---|
| Specification | `t3fund_proca_b412_t7_drafts/B4_12_COMMON_M_TEST_PROPOSAL_v0.1.2.md` (commit 3cba04e) |
| Script | `t3fund_proca_b412_t7_drafts/b412_common_m_check.py`, SHA-256 `5a8ad4466b4b6734aa26491be020dc72151799cd10a0c90c5b6c793317bf777e` |
| Evidence | `evidence/B4_10_proca_mass_sweep_v0_5/`, extracted with `git archive 511f6fe`; `SHA256SUMS.txt` 83/83 before the run |
| Environment | Python 3.13.18, standard library only |
| Outputs | `t3fund_proca_b412_t7_drafts/results_b412_common_m/`: `b412_common_m.json`, `b412_common_m.csv`, `terminal_output.txt`, `exit_code.txt` (0), `RUN_CONTEXT.txt`. These are the files of the reviewed run archive `b412_common_m_run_20261115.zip` (SHA-256 `042cf43899fedf1b45b8adf997e8e07fd0110e5ef38b724abef8ecca402205a1`), unchanged |

Controls: K1 PASS (denominators 40/40/41/12/29, anchor among the 29; historical 0, 0, 41 and 29 reproduced), K2 PASS, K3 PASS (max deviation 2.2·10⁻¹⁶), K4 max δ_rank1 = 2.2·10⁻¹⁶ (|κ|/√p agrees numerically with √s), K5 PASS. Verdict `ALL_CHECKS_PASS`. The reviewer independently recomputed the set-C margins, M intervals, classes and c_max for all 41 backgrounds from the CSV and found no discrepancy.

## 2. Results

| Set | (g_max, ρ_max) | Original, all 41 | Original, non-control 29 | Common-M, all 41 (PASS / FAIL / BORDERLINE) | Common-M, non-control 29 |
|---|---|---|---|---|---|
| A | (1, 0.3) | 0 | 0 | 0 / 41 / 0 | 0 / 29 / 0 |
| B | (1, 0.1) | 0 | 0 | 0 / 41 / 0 | 0 / 29 / 0 |
| C | (√(4π), 0.3) | 41 | 29 | **21 / 20 / 0** | **15 / 14 / 0** |

Set C, Common-M, by class:

| Class | PASS | FAIL |
|---|---|---|
| stable (16) | 8 | 8 |
| boundary (12) | 7 | 5 |
| negative_control (12) | 6 | 6 |
| anchor (1) | 0 | 1 |

By mode: rank1_minus 12 PASS / 9 FAIL; rank1_plus 9 PASS / 11 FAIL. The case closest to the border is `rank1_minus_boundary_4487`, with m = −0.0110 against g_max·ρ_max = 1.0635. That is far outside the BORDERLINE tolerance of 1.06·10⁻¹², so its FAIL is stable.

**ANCHOR_B3 (K6, recomputed by hand and by the reviewer).** c_max = max(√p, |κ|/√p) = max(0.2252, 0.5713) = 0.5713 and G = |q| = 4.2372. So c_max·G = 2.421 > 1.0635, and m = −1.357. The admissible interval for M is [G/ρ_max, g_max/c_max] = [14.12, 6.205], which is empty.

## 3. Interpretation

**Mechanism.** The set-C failure is a conflict between two constraints on the single mass M:
- the largest background gradient (for the anchor, q) sets a lower bound, M ≥ G/ρ_max;
- the largest coupling coefficient (for the anchor, that of ψ) sets an upper bound, M ≤ g_max/c_max.

When the gradient and the coefficient belong to different fields, the two bounds can be incompatible. The original criterion compares √X_j = (coefficient of field)·(gradient of the same field) component by component, so it cannot see this conflict. This is **not** an additional physical term in the Lagrangian; it is incompatibility of constraints on one shared parameter.

**What the result establishes**
1. The historical check is reproduced: Original set C is 41/41, and 29/29 without the negative controls.
2. The stricter condition is different: a common M with g ≤ √(4π) and ρ ≤ 0.3 for both fields exists for 21 of 41 backgrounds, and for 15 of 29 without the negative controls.
3. The difference is not caused by a rank-one mismatch. δ_rank1 ≤ 2.2·10⁻¹⁶, so |κ|/√p and √s are numerically identical on these backgrounds.

**What the result does NOT establish**
- It does not show that the Proca construction is physically impossible, or that it is UV-incompatible in general.
- It is a statement about compatibility of the stated constraints (g ≤ g_max, ρ ≤ ρ_max, one common M) on these 41 backgrounds only. The backgrounds are strong-X by construction, selected from the cut-0.95 population.
- It does not change the thresholds, the definition of ρ, or the historical results. It does not touch the front-velocity population result of B4_12_REPORT v0.5.
- It says nothing about Q2/Q3.

**Precise reading of the historical 41/41.** The historical set-C count of 41/41 means that the necessary condition √X_max ≤ √(4π)·0.3 holds on all 41 backgrounds. It does not mean that each of these backgrounds lies in control set C, i.e. that a common Proca mass exists with g ≤ √(4π) and ρ ≤ 0.3 simultaneously. That holds for 21 of 41.

## 4. Audit of existing wording (no file is changed by this draft)

| # | Location | Status of file | Wording | Issue | Proposed handling |
|---|---|---|---|---|---|
| W1 | `evidence/B4_10_proca_mass_sweep_v0_5/B4_12_REPORT.md`, lines 12 and 15 | historical evidence, must not be edited | "C: g≤√(4π) (α≤1), ρ≤0.3 (√X≤1.06)" … "set C: 41/41 (29/29 non-controls)" | Presents the set as (g, ρ) bounds, but the count is of the necessary condition only | No edit. The present result document is the correction of record; `PROVENANCE.md` already states that the set test is a necessary condition |
| W2 | same file, line 16 | historical evidence | "Minimum g required for ρ≤0.3: 1.85 … 3.25 (anchor 3.18, α≈0.81)" | Computed as √X_max/0.3 (per component). With one common M, the minimum is g_min = c_max·G/ρ_max. For the anchor that is 2.421/0.3 ≈ 8.07, with α_min = g_min²/4π ≈ 5.18, not 3.18 (α ≈ 0.81) | No edit. **Decision:** only the anchor value is reported, as a diagnostic example. g_min is a derived quantity not defined in spec v0.1.2. No table for the 41 backgrounds is given. If one is wanted later, it needs a separately defined derived diagnostic, without changing v0.1.2 |
| W3 | `t3fund_proca_spec_bundle/T3_FUND_PROCA_SPEC.md`, line 53 | **registered** (e355af2), must not be edited | "B4.12: control sets A,B 0/41, set C only at α≈1" | "Only at α≈1" follows W2. Under the common-M condition, for 20 of 41 backgrounds even α = 1 (g = √(4π)) is not enough at ρ ≤ 0.3. The sentence is labelled "for orientation only" and is not a criterion of the protocol | No edit to the registered file. **Decision:** a dated clarification note (§4a), no formal erratum. The note does not change the registered file and is not a correction of it |
| W4 | `t3fund_proca_b412_t7_drafts/B4_12_PROCA_PERTURBATIVITY_AND_UV_CLOSURE.md` v0.3, §1 items 1–3 | unregistered draft | under "Established": "\|g_ψ\| = M√s" and "C 41/41 (29/29 without the controls)" (translated from the Bulgarian original) | The code coefficient is \|κ\|/√p (equal to √s only at rank one); 41/41 is listed as established without the necessary-condition qualifier | **Decision:** corrected in the B4.12 proposal v0.4, a separate later task. v0.4 must state \|g_ψ\| = M\|κ\|/√p (equivalent to M√s only under the rank-one condition), and give the necessary-condition reading of the historical set-C count of 41/41 together with the common-M result. Q1 can cite this document |

### 4a. Clarification (2026-11-15)

The historical B4.12 set-C count of 41/41 refers to the necessary component-wise condition only. The subsequent common-M test finds simultaneous compatibility with g ≤ √(4π) and ρ ≤ 0.3 for 21/41 backgrounds (15/29 non-controls). The registered orientation statement in `t3fund_proca_spec_bundle/T3_FUND_PROCA_SPEC.md`, line 53 ("set C only at α≈1"), is not a common-M result.

*This note does not change the registered file and is not a correction of it.*

### 4b. Not affected

Checked and **not** affected: `docs/UV_CANDIDATE_ASSESSMENT.md` (refers to `B4_12_REPORT.md` without repeating 41/41); the C4 and Stage 1 result bundles (they cite only "0/41 for g≤1", which the common-M result confirms); the T7a manuscript v0.3 (does not cite B4.12 numbers).

## 5. Next steps

1. Commits, separately and in this order (approved 2026-11-15, 15:07): (1) the script; (2) the run outputs in `t3fund_proca_b412_t7_drafts/results_b412_common_m/`, unchanged; (3) this document.
2. Not part of these commits: any edit to the registered `T3_FUND_PROCA_SPEC.md`, to the historical evidence, or to spec v0.1.2.
3. Separate later task: the B4.12 proposal v0.4 (W4, Q1), then the working draft.
