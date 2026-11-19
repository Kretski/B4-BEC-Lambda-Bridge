# BEC Stage 2: authorised test runs 1 and 2 — record v0.1

**Scope:** the test suite of the Stage 2 script, which uses synthetic inputs only. **No Stage 2 run.** No registered target was evaluated, searched or sampled. The Stage 2 run itself still needs a separate authorisation.

## Test run 1: registered script v0.4.1 (commit 2410aef)

- **Authorised:** 2026-11-20.
- **Setup:** fresh GitHub clone, LF line endings, clean before and after the run. Input hashes and the SciPy pin (1.18.1, trf.py 0c7c29d5…) were verified first.
- **Result:** rc = 1. **15 passed, 13 errors at fixture setup.**
- **Cause:** `load_module` did not register the module in `sys.modules` before `exec_module`. `@dataclass` in `theory/lambda_eff.py` and `theory/bec_dispersion.py` then failed under Python 3.13 with `AttributeError`.
- **Consequence:** a real run of v0.4.1 would have stopped at module loading, before any computation.
- **Records:** `test_runs/run1_v041/` (stdout 6e151d5e…, RUN_CONTEXT a085bbdf…).

## Fix: script and tests v0.4.2

- `load_module` registers the module in `sys.modules` before execution and removes it if execution fails. It ABORTs if the name is already in use.
- The test helper `_load` gets the same fix.
- New test T27 covers the fix.
- No other change to logic, criteria or constants.
- Script d315eda3…, tests c9b7555e….

## Test run 2: candidate v0.4.2, before registration

- **Authorised:** 2026-11-20.
- **Setup:** fresh clone at 2410aef (clean before and after). The candidate files were placed outside the clone. Inputs and the SciPy pin were verified.
- **Result:** rc = 0. **29 passed**: 28 test functions, with T24 parametrised over regions A and B.
- **T15:** the expectation `n_calls == nfev + 4·njev` under scipy 1.18.1 is **confirmed**. It was not adjusted.
- **Records:** `test_runs/run2_v042/` (stdout 58d52bd8…, RUN_CONTEXT ad48fa8d…).

## What this does and does not establish

- **Established:** the implementation behaves as specified on the synthetic cases the tests cover. This includes the energy form, classification, criteria, the residual wrapper and its call counting, R10, C1c, the perturbation rules, the bookkeeping of the positive control and module loading.
- **Not established:** anything about B4 or the BEC family. No physical result exists yet.
- **Not executed:**
  - the positive control at full size (29 targets × 80 starts per region);
  - the search on the 29 registered targets;
  - the sign-tuple sampling.
