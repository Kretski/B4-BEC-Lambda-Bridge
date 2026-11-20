# BEC Stage 2: test runs 1–3 and the first run attempt — record v0.2

**Supersedes:** v0.1 (afabc20c…). v0.2 adds the first Stage 2 run attempt, the v0.4.3 fix and test run 3. Runs 1 and 2 are summarised, and their records stay in `test_runs/run1_v041` and `test_runs/run2_v042`.

**Scope:** test suite and one aborted run attempt. **No Stage 2 computation has taken place.** No registered target was evaluated, searched or sampled.

## Test runs 1 and 2 (summary; details in v0.1)

| run | code | result |
|---|---|---|
| run 1 | registered v0.4.1 (2410aef) | rc = 1: 15 passed, 13 setup errors. `load_module` did not register the module in `sys.modules`, so `@dataclass` failed under Python 3.13 |
| run 2 | candidate v0.4.2 | rc = 0: 29 passed. T15 expectation confirmed |

v0.4.2 was registered in 0eeecaa.

## First Stage 2 run attempt (Lenovo, registered v0.4.2, 0eeecaa) — ABORT before any computation

- **Setup:** fresh LF clone of 0eeecaa; conda environment with scipy 1.18.1.
- **Outcome:**
  - The input-hash check passed, and so did the SciPy **version** check.
  - The **trf.py hash check ABORTed**: raw sha256 a174a6e7… against the pinned 0c7c29d5….
  - `out/` was created and is empty. `run_stdout.txt` is empty.
  - The stderr is recorded as `run_attempts/attempt1_lenovo_v042/run_stderr.txt` (as uploaded, sha256 3eaae0ab…).
- **Diagnosis:**
  - The Windows trf.py of scipy 1.18.1 has CRLF line endings. Converting the pinned LF file to CRLF gives exactly a174a6e7…, so the source text is identical.
  - A fresh pip-wheel environment on the same machine gave the same CRLF hash.
- **Reading:** the guard worked as intended. The ABORT came before any computation, and no result was produced.

## Fix: script and tests v0.4.3

- `verify_scipy` compares the sha256 of trf.py **with CRLF normalised to LF** against the pin. It records the raw sha256 and the line-ending kind (LF, CRLF or mixed).
- Only CRLF pairs are rewritten. A lone CR is kept.
- New tests:
  - **T28:** the LF and CRLF forms give equal normalised hashes; different text and a lone CR do not.
  - **T29:** the installed SciPy matches the pin, whatever its line endings.
- No other change.
- Script 3f5d5ea1…, tests 13099c65….

## Test run 3: candidate v0.4.3, before registration

- **Setup:** fresh clone at 0eeecaa, clean before and after. The candidate files were placed outside the clone. Input hashes and the SciPy pin were verified.
- **Result:** rc = 0. **31 passed**: 30 test functions, with T24 parametrised over regions A and B.
- **Records:** `test_runs/run3_v043/` (stdout d73c655f…, RUN_CONTEXT f5c9b2ae…).
- **Limitation:** the test machine has an LF trf.py. The CRLF path of the pin is covered there only by the synthetic T28. On Lenovo, T29, or the start of the run itself, exercises it on the real file.

## What this does and does not establish

**Established:** the implementation behaves as specified on the synthetic cases the tests cover.

**Not established:** any physical result about B4 or the BEC family.

**Not executed:**
- the full-size positive control;
- the 29-target search;
- the sign-tuple sampling.
