# Proca with a timelike current: Q_new(c; 1) roots and current attribution

Decision note: `docs/PROCA_TIMELIKE_QNEW_DECISION_NOTE_v0.1.md`. Result: **NOT-VIABLE by criterion (d), for both signs of γ_u** (reviewer, 2026-11-19), within the registered linear model and the checked backgrounds. The specification v0.5 (87fc31b) is unchanged.

| Path | SHA-256 | Role |
|---|---|---|
| `proca_timelike_qnew_roots.py` | `91ea9152b1aef0277369044d9e5d04bbd68b5f6d29d83eaec3db01b68b795a30` | Q_new(c; 1) root arithmetic, 41 backgrounds × σ = ±1 (one authorised run) |
| `smoke_qnew.py`, `smoke_qnew_output.txt` | `dcbebce6…`, `47f1c2ef…` | synthetic-only smoke test (14/14) |
| `proca_timelike_qnew_current_attribution.py` | `21c0934f676ce8e2de0a94bf47d8201ddcc5c518eb549824477a69d711c7a0b4` | attribution check, 4 backgrounds × σ = ±1 (one authorised run) |
| `smoke_attribution.py`, `smoke_attribution_output.txt` | `19865f06…`, `85aeb652…` | synthetic-only smoke test (9/9) |
| `revision_diffs/` | — | the reviewed revisions of both scripts and both smoke tests |
| `results_qnew/` | see `SHA256SUMS` | result of the Q_new run (immutable); `posthoc/` is a read-only comparison with the registered T3-fund Proca result |
| `results_attribution/` | see `SHA256SUMS` | result of the attribution run (immutable) |

## Reproduction of the smoke tests (synthetic inputs only)

From the repository root, on an LF checkout:

    python3 -I experiments/uv_candidates/proca_timelike/qnew/smoke_qnew.py \
        experiments/uv_candidates/proca_timelike/qnew/proca_timelike_qnew_roots.py
    python3 -I experiments/uv_candidates/proca_timelike/qnew/smoke_attribution.py \
        experiments/uv_candidates/proca_timelike/qnew/proca_timelike_qnew_current_attribution.py \
        experiments/uv_candidates/proca_timelike/qnew/proca_timelike_qnew_roots.py

The two authorised runs have been done; reruns are not part of the registration. Both runs used the scripts outside the checkout, from the root of a fresh clean LF clone (see each `RUN_CONTEXT.txt`).
