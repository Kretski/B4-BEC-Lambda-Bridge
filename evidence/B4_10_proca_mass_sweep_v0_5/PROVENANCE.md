# Provenance — B4.10 Proca mass sweep, package v0.5

**Status: HISTORICAL EVIDENCE. NOT A REGISTERED RESULT. CONTENT NOT EDITED.**

This folder is the unmodified content of the archive `B4_10_proca_mass_sweep_v0_5.zip`. The package predates the repository: `docs/UV_CANDIDATE_ASSESSMENT.md` already referred to `B4_12_REPORT.md` when it was added at 05b9a6b (2026-11-11). It is the primary source for:
- the B4.12 definition of ρ (`b410/uv.py`, lines 3–4 and 16–36);
- the historical B4.12 results (sets A/B/C: 0/41, 0/41, 41/41; 29/29 for set C without negative controls), produced by `run_b412.py` into `results_b412/`.

## Source archive

| Item | Value |
|---|---|
| File name | `B4_10_proca_mass_sweep_v0_5.zip` |
| SHA-256 | `714cf23e22df30ee319d36a8de849f53a523e04bd0786e527e468786bf34663c` |
| Located | 2026-11-15, on the author's original laptop, during the search for the B4.12 source of ρ |
| Files | 83; no `__pycache__` or `.pytest_cache` entries, so nothing was excluded |

## What was done to the content

The files were extracted byte for byte and nothing else was changed. The archive's top-level folder `B4_10_proca_mass_sweep_v0_5/` is this folder. `SHA256SUMS.txt` lists all 83 original files with paths relative to this folder; check them with `sha256sum -c SHA256SUMS.txt`.

## Files added (not part of the original archive)

| File | Purpose |
|---|---|
| `SHA256SUMS.txt` | SHA-256 of the 83 original files |
| `PROVENANCE.md` | this file |
| `../.gitattributes` | `* -text`: prevents line-ending conversion in `evidence/`. 19 CSV files in the package use CRLF line endings. Without this rule Git on Windows could store them with LF, which would change their content and make the vendored B4.9 manifest check in `run_b412.py` fail. |

## Scope

- Historical material. The results, criteria and conclusions in these files are those of the original B4.10–B4.12 work. They are not re-evaluated here, and they are not registered results under the project's later protocols.
- The B4.12 set test in `b410/uv.py` (`√X_max ≤ g_max·ρ_max`) is a necessary condition for the existence of a common Proca mass M. A separate, stricter common-M test is specified in `t3fund_proca_b412_t7_drafts/B4_12_COMMON_M_TEST_PROPOSAL_v0.1.1.md`. Its results are reported outside this folder and do not modify it.
