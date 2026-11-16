# Script A — revision 1.1: corrections from the review of 2026-11-16

**Status: FOR FINAL STATIC REVIEW. Not run on the 41 registered backgrounds.**

| File | SHA-256 |
|---|---|
| `uv_dynamic_eta_check.py` rev 1.1 (1373 lines) | `fa46712ba0136f63e207ea95f6f21e3de8fa9cccce2514cdc48aedc271bb874e` |
| `test_script_a_rev1_1.py` | `bc5a2467e9de5552d0c62a6533a277d9b68516034d953b2844af9805b6cb35ef` |
| `rev1_0_to_rev1_1.diff` | unified diff against rev 1.0 (`f920ecea…`) |

`py_compile`: OK. Unit tests: 26 run, 0 failed. Command:

```
python3 -I test_script_a_rev1_1.py <script_dir> <repo_root>
```

## What was not changed

The following are unchanged from rev 1.0:

- the Lagrangian;
- the 41 backgrounds;
- f ∈ {0.1, 1, 10};
- the P2 grid;
- all thresholds and tolerances;
- the R2 points and steps;
- the P6 rule.

The diff contains only changes to check logic, multiplicity handling, input protection and boundary-case handling.

## Changes

| Review item | Change in rev 1.1 | Test |
|---|---|---|
| 1.1 P2a multiplicity | `p2a(DP_full, …)` now uses the full determinant with all multiplicities. The square-free product is used only for distinct roots (P3, P4, P6, low-k roots). | T4a–c: the full determinant gives multiplicity 4 → ZERO MODE; the square-free product would give 1 → INCONCLUSIVE; `main` calls P2a with `DP_full` |
| 1.2 T6 confirmation | New function `t6_class(...)`. STATED-DYNAMICAL is assigned only when `conf is True`. False or None → INCONCLUSIVE. | T5a–e |
| 1.3 C1 and block diagonality | `c2_numeric(blocks, vals, off0)` returns INCONCLUSIVE ("principal symbol not block-diagonal") when `off0` is False. `c1_from_c2` gives an analytic PASS only from a C2 PASS. | T6a–c |
| 1.4 input hashes | `EXPECTED_SHA256` is fixed for both files (B3 `99a6b04b…`, manifest `1bae892a…`, from repo fa19637). `verify_inputs()` runs first in `main` and again in `load_backgrounds`; any mismatch → ABORT. | T8a–d: unchanged copies pass; one byte added to B3 or to the manifest → ABORT; the check runs before any computation |
| 2.1 singular F_ηη | If \|F_ηη\| ≤ 1e-40 (the registered BOUNDARY threshold), `s7_numeric` returns `singular`. Local exceptions in S7 are caught and recorded. In both cases T7 is INCONCLUSIVE for that root and f, and the remaining backgrounds continue. | T7a–d |
| 2.2 DC5 overlap | `dc5_match`: the distinct roots found must equal the de-duplicated union of the J′ = 0 and J = 0 classes. A shared root is counted once. | T7e–h |
| 2.3 multiple roots | `roots_with_multiplicity` takes exact multiplicities from the square-free decomposition (`sqf_list`). A root with multiplicity > 1, or with relative separation from another root below 1e-20, gets `resolved = False`, and every class for it is INCONCLUSIVE (T3-fund, T6, T7). | T7i (double root z = 2; close roots), T7j |
| choice 6 | The F2(ii′) text now says the elimination procedure is not admissible in this regime. This is not a proof that no UV origin exists. The same wording is in the `t7_class` docstring. | — |
| X2 visibility | X2 (maximum difference between the low-k roots of the full system and those of the eliminated operator) is now printed on every per-root log line, in addition to the JSON. It stays outside ALL_DC_PASS. | log line |

## Other checks

The end-to-end smoke test on the synthetic background (not a registered record) was repeated: rc 0, behaviour as in rev 1.0. Its physics outputs are not results.

## Status

Running on the 41 registered backgrounds is NOT authorised.
