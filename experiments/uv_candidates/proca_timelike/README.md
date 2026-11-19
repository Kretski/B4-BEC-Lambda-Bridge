# Proca with a timelike current: Stage 1, D5 addendum and Stage 2

Specification: `docs/PROCA_TIMELIKE_CURRENT_ANALYTIC_SPEC_v0.5.md` (registered in 87fc31b, unchanged).

Documents (in `docs/`):

- `PROCA_TIMELIKE_STAGE1_REPORT_v0.1.md` — the Stage 1 derivation. Its checks D2c, D3a, D5b and D5e are recorded there as unconditional; they are computed in the D5 addendum;
- `PROCA_TIMELIKE_STAGE1_D5_ADDENDUM_v0.3.md` — D5, accepted by the reviewer on 2026-11-19 within its stated scope (linear system about one homogeneous background, 1+1D);
- `PROCA_TIMELIKE_STAGE2_RESULTS_NOTE_v0.1.md` — Stage 2 results. Stage 2 accepted; model class **OPEN**.

The headers of the two Stage 1 documents still say "for review"; the documents were accepted by hash and are registered unchanged.

Not run: roots of Q_new(c; 1), V3, NOT-VIABLE (d), Stage 3.

## Files

| Path | SHA-256 | Role |
|---|---|---|
| `stage1/proca_timelike_stage1_check.py` | `b59393d2…` | Stage 1 symbolic checks (33/33) |
| `stage1/stage1_output.txt` | `138550c5…` | its output |
| `stage1/proca_timelike_d5_addendum_check.py` | `2064fa84…` | D5 addendum symbolic checks (13/13); loads the Stage 1 check from the same folder |
| `stage1/d5_addendum_output.txt` | `80d8201b…` | its output |
| `stage2/proca_timelike_stage2.py` | `585143fa50d390b27e40b6aa1d8a2114a32e2b881600d6ee21d0cc22e27a5ff6` | **the approved Stage 2 script; the only script used for the registered results** |
| `stage2/proca_timelike_stage2_ff359c11.py` | `ff359c11131fefafc90d43436a6add37c7f79a4fe33811eae6a3fcf8c7df8048` | **archived earlier revision, NOT the approved Stage 2 script.** Kept only as the comparison input of the smoke test. Not a current alternative implementation; never used to produce results |
| `stage2/smoke_stage2_fix2.py` | `d2ea9272…` | synthetic-only smoke test (no registered background is read) |
| `stage2/smoke2_output.txt` | `e2fd161a…` | its output (18/18) |
| `stage2/revision_diffs/` | — | the reviewed revisions 16fa3506 → 9c26bce4 → ff359c11 → bd6eb9fc → 585143fa |
| `results_stage2/` | see `SHA256SUMS` | the single authorised run (immutable); `posthoc/` holds presentation tables derived read-only from the JSON |

## Reproduction

From the repository root (LF checkout):

    python3 -I experiments/uv_candidates/proca_timelike/stage1/proca_timelike_stage1_check.py
    python3 -I experiments/uv_candidates/proca_timelike/stage1/proca_timelike_d5_addendum_check.py
    python3 -I experiments/uv_candidates/proca_timelike/stage2/smoke_stage2_fix2.py \
        experiments/uv_candidates/proca_timelike/stage2/proca_timelike_stage2.py \
        experiments/uv_candidates/proca_timelike/stage2/proca_timelike_stage2_ff359c11.py

The Stage 2 run itself (`python3 -I experiments/uv_candidates/proca_timelike/stage2/proca_timelike_stage2.py <outdir>`) was authorised once and has been done. Its result is in `results_stage2/`. A rerun is not part of the registration. The run was made with the script outside the checkout, from the root of a clean LF clone of 87fc31b (`results_stage2/RUN_CONTEXT.txt`).
