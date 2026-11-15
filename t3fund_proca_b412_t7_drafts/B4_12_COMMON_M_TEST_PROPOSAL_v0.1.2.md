# B4.12 — Strict common-mass (common-M) test and the B4.10 v0.5 evidence package (proposal v0.1.2)

**Status: PROPOSAL. Nothing has been computed.** This document implements steps 2–4 of the reviewer decision of 2026-11-15, 14:16. It does not change:
- the definition of ρ;
- the thresholds of sets A, B and C;
- the historical result of B4.12 v0.5;
- B4.12 v0.3, the working draft v0.2.1, T7a or the T7 protocol.

The computation is run only after this proposal is approved and committed, and after the committed files have been checked on GitHub.

**Changes from v0.1.1** (review of 2026-11-15, 14:24, APPROVED AFTER K4 CLARIFICATION):
1. The document is in English, the language of the repository. v0.1 and v0.1.1 were Bulgarian drafts and were not committed.
2. K4 is made unambiguous: the code coefficient |κ|/√p is **always** used in the common-M criterion; K4 only measures and reports the deviation from the rank-one identity (§5).

**Changes from v0.1** (review of 2026-11-15, 14:20, APPROVED WITH MINOR CLARIFICATIONS):
1. K1: the report shows the denominator, the number passing, the number of negative controls, and how 29/29 is obtained (§5, K1).
2. BORDERLINE is a third category, separate from PASS and FAIL (§3).
3. The full evidence package was chosen, with protection of the bytes against line-ending conversion (§6).
4. The two limits stated in the review were added (§7).

---

## 1. Primary source (item 1, closed at document level)

Archive `B4_10_proca_mass_sweep_v0_5.zip`. SHA-256 of the archive as uploaded: `714cf23e22df30ee319d36a8de849f53a523e04bd0786e527e468786bf34663c`.

| File | Lines | SHA-256 |
|---|---|---|
| `B4_12_REPORT.md` | 5–16 | `31e0f43bf369303c9596273ef320cf9218ad739d30449139fb0c21aa9f04cd96` |
| `b410/uv.py` | 3–4, 16–36 | `db32b8ec76e3ccf876f49d2303e4d23ef8cf1b50f26e6157d6f2fb0b7db0c3be` |
| `run_b412.py` | 13–39 | `e0b83bbb4ed52e12364a14d5a30919ad28852c123f913fbd86626ed3f15aecbd` |
| `results_b412/b412_summary.json` | — | `995eac0251da75d0f98b2d51183ad88e0b694e034f27466b97095a7b7b414ff9` |

**Definition (accepted by the reviewer).** ρ_j = |u_j|/M for u_j ∈ {q, w, r, v}, where M is a single mass common to both fields. The couplings are g_φ = M√p for q and w, and |g_ψ| for r and v. The identity g_j ρ_j = √X_j holds for each component separately.

**Detail from the code.** `fixed_eft.couplings` sets g_ψ = −Mκ/√p, not M√s. The two agree only when κ² = ps. All 41 backgrounds are rank1_minus or rank1_plus. The code coefficient is used throughout, and the deviation from κ² = ps is reported by control K4 (§5).

## 2. Inputs (unchanged)

The 41 backgrounds are exactly those of `run_b412.py`, lines 15–17:
- `ANCHOR_B3`, from `data/B3_max_speed_7_315096.json`;
- 40 backgrounds from `data/e3_manifest.json`, key `selected`: 16 stable, 12 boundary and 12 negative_control; 20 rank1_minus and 20 rank1_plus.

**Checked on the inputs only, before any result:** p > 0 and s > 0 for all 41 backgrounds (min p = 0.0507, min s = 0.0524). The question of gradients of a field with zero coupling therefore does not arise.

The population part of B4.12 (2 × 10 000 samples, seed 172026, front velocity) is out of scope.

## 3. The two criteria

Let G = max(|q|, |w|, |r|, |v|) and c_max = max(√p, |κ|/√p), where |κ|/√p is the coefficient of |g_ψ|/M in the code.

- **Original criterion (historical, reproduced only):** √X_max ≤ g_max·ρ_max, with √X_max = max_j √X_j (`uv.controlled`, lines 31–34). This is a necessary condition.
- **Common-M criterion (new, strict):** there exists M > 0 such that g_φ ≤ g_max, |g_ψ| ≤ g_max and ρ_j ≤ ρ_max for all four j, simultaneously. This is equivalent to

  **c_max · G ≤ g_max · ρ_max**,

  and the admissible interval for M is [G/ρ_max, g_max/c_max].

Both criteria are applied to A, B and C with the thresholds of `uv.SETS`, unchanged:
- A: (1, 0.3);
- B: (1, 0.1);
- C: (√(4π), 0.3).

**Common-M classification (fixed in advance).** The margin is m = g_max·ρ_max − c_max·G. For each background and set:
- **PASS:** m > 10⁻¹²·g_max·ρ_max;
- **FAIL:** m < −10⁻¹²·g_max·ρ_max;
- **BORDERLINE:** |m| ≤ 10⁻¹²·g_max·ρ_max.

BORDERLINE is a third, separate category. It is counted neither as PASS nor as FAIL, and it is not reclassified. For each set, PASS + FAIL + BORDERLINE = 41. The original criterion is reproduced as a boolean flag exactly as in `uv.controlled`, with no BORDERLINE category.

## 4. Implementation

A new standalone script `b412_common_m_check.py` in `t3fund_proca_b412_t7_drafts/`:
- **Standard library only:** `json`, `math` and `hashlib`. It does not import `b410`, numpy or scipy, so it is independent of the code it checks.
- **Reads only two input files** from the evidence package (§6): `data/B3_max_speed_7_315096.json` and `data/e3_manifest.json`.
- **Checks the inputs first:** the SHA-256 of both files against the package's `SHA256SUMS.txt`. On a mismatch it stops with FAIL.
- **Outputs:**
  - `b412_common_m.json`: counts per set under both criteria, BORDERLINE cases, and the controls;
  - `b412_common_m.csv`: per background tag, cls, mode, p, s, √X_max, G, c_max, the Original and Common-M flags for A/B/C, and for set C the margin m and the M interval.

## 5. Controls (all mandatory)

| # | Control | PASS condition |
|---|---|---|
| K1 | Reproduction of Original | exact agreement with `results_b412/b412_summary.json` (keys `n`, `controlled_counts`, `controlled_counts_noncontrol`), reported as the table below |
| K2 | Implication | for every background and set: Common-M PASS ⇒ Original true. The converse is not required. A violation indicates an error in the script |
| K3 | Identity | \|g_j ρ_j − √X_j\| ≤ 10⁻¹², computed at M = 7 as in `run_b412.py` |
| K4 | Rank-one deviation (report only) | δ_rank1 = \|κ²/(ps) − 1\| is computed for every background and its maximum is reported. The common-M criterion **always** uses the code coefficient \|κ\|/√p, whatever the value of δ_rank1. If max δ_rank1 ≤ 10⁻⁹, the report states that \|κ\|/√p agrees numerically with √s. If max δ_rank1 > 10⁻⁹, the report gives the deviation and the backgrounds concerned. The formula is never switched during the test |
| K5 | Inputs | the SHA-256 of both input files match `SHA256SUMS.txt` |
| K6 | Hand check | for `ANCHOR_B3`, the values of c_max, G, m and the M interval are printed in the report so they can be checked by hand |

**Format of K1 (for each of the sets A, B, C).** The report shows explicitly:

| Column | Value |
|---|---|
| all backgrounds: denominator | 41 (`ANCHOR_B3` plus 40 from the manifest) |
| all backgrounds: passing | number with Original flag true |
| negative controls | number with `cls == "negative_control"` (from the inputs: 12) |
| without negative controls: denominator | 41 − 12 = 29 (`ANCHOR_B3` has `cls = "anchor"` and is counted here) |
| without negative controls: passing | number with Original flag true and `cls != "negative_control"` |
| historical values | `controlled_counts[set]` and `controlled_counts_noncontrol[set]` from `b412_summary.json` |
| agreement | yes/no |

The historical 41/41 means 41 passing out of a denominator of 41. The historical 29/29 means 29 passing out of a denominator of 29, i.e. all backgrounds except the 12 negative controls. The same table is given for Common-M, with three columns PASS / FAIL / BORDERLINE.

The verdict `ALL_CHECKS_PASS` covers only these controls; K4 is a report, not a pass/fail control.

## 6. Evidence package B4.10 v0.5

**Why the minimal set is not enough.** A minimal set (`B4_12_REPORT.md`, `b410/uv.py`, `run_b412.py`, `results_b412/`, `SHA256SUMS.txt`) cannot run `run_b412.py`:
- `uv.py` imports `fixed_eft` and `experiment2`;
- `run_b412.py` imports `_b49` (which checks the vendored B4.9 against its manifest), `causality`, `experiment3` and `vendor_t173/robustness_scan.py`;
- the data are in `data/`.

**Decision (review of 14:20): full package.** The package goes in whole and unmodified to `evidence/B4_10_proca_mass_sweep_v0_5/`. The archive holds 83 files, about 5.1 MB, including the vendored B4.9 and all reports (B4.11, B4.12, E3, audits). The archive contains no `__pycache__` or `.pytest_cache` entries, so nothing is excluded.

**Protection of the bytes.** 19 of the files are CSV files with CRLF line endings. They include `results_b412/b412_control_parameters.csv` and the five CSV files of the vendored B4.9. Without protection, Git on Windows (`core.autocrlf`) could store them with LF. Their content would then no longer be the original, and `_b49.verify_manifest()` in `run_b412.py` would reject the vendored B4.9. An `evidence/.gitattributes` file with the line `* -text` is therefore added. It disables all line-ending conversion under `evidence/`.

The following are added to the package:
- `SHA256SUMS.txt`, generated from the extracted files, with paths relative to the package root;
- `PROVENANCE.md`: the SHA-256 of the original archive, the date it was located (2026-11-15), and a note that the content is unedited historical material, not a registered result.

`SHA256SUMS.txt` lists only the 83 original files. The three added files (`SHA256SUMS.txt`, `PROVENANCE.md`, `evidence/.gitattributes`) are identified as added in `PROVENANCE.md`.

## 7. What this step does NOT do

- It does not change the definition of ρ, the thresholds, or the historical 0/41, 0/41, 41/41.
- It does not reinterpret the front-velocity conclusion of B4_12_REPORT v0.5.
- It does not touch Q2/Q3. Updating the working draft and the decision table is step 5, after the comparison.
- If set C does not pass Common-M for all 41 backgrounds, the thresholds are not changed and no way is sought to preserve 41/41. The actual result is reported.
- The check concerns only the compatibility of the constraints on the common mass M. On its own it does not establish causality, a UV completion, or the physical admissibility of the model.
- No Zenodo.

## 8. Order

1. Review of this proposal.
2. Commit of the proposal and of the evidence package (two separate commits).
3. Check of the committed files on GitHub: all files byte for byte, the 19 CRLF CSV files, and `SHA256SUMS.txt`.
4. Script and run.
5. Report comparing Original and Common-M; review; commit.
6. Update of the working draft.
