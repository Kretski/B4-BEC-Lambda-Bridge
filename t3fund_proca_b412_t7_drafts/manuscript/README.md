# T7a manuscript (draft)

**Status: DRAFT MANUSCRIPT. T7a IS NOT REGISTERED.** This folder holds the technical manuscript that presents the approved T7a draft results. It is not a registered result. It does not evaluate any criterion of the draft T7 protocol (PASS-FORMAL, PASS-VERIFIED, FAIL, INCONCLUSIVE). It is not evidence that the EFT or its UV completion is consistent. B4.12 is OPEN, and nothing from B4.12 is used here. No Zenodo version is made from this folder.

## Contents

| File | What it is |
|---|---|
| `T7A_MANUSCRIPT_DRAFT_v0.2.md` | *Tree-Level Proca Integration and Independent IBP Classification of Higher-Derivative Operators in a Two-Scalar EFT*, draft v0.2, English, Markdown |
| `README.md` | this file |

## What the manuscript contains

- The tree-level integration of the Proca field from the full action, with the Gauss constraint retained:
  - the exact kernel K(ω,k) = (1+k²/M²)/(M²+k²−ω²);
  - p = g_φ²/M², s = g_ψ²/M², κ = −g_φg_ψ/M², κ² = ps (rank-one subfamily of B4 only).
- The first correction ΔL^(M⁻⁴) = (∂_tJ)²/(2M⁴) and its off-shell, IBP-only classification (option A): dim(V_S/~) = 3 / 12 / 3 in the f⁴, f²g², g⁴ sectors.
- The agreement of two algorithmically independent implementations (Euler operator / explicit IBP matrix): identical bases and exact rational coordinates.

It does not contain: numerical verification under a protocol, order M⁻⁶, T7b (emergent Λ), or any statement on validity, unitarity or causality.

## Source commits

Only these commits are used. The files they introduced are not modified by this folder.

| Commit | Content |
|---|---|
| 61fefff | step 3: derivation and SymPy check (v0.5.2) |
| a72c323 | option A equivalence rules (protocol draft v0.6, §2) |
| e7ab8d4 | canonical-representative rule v0.7.1; classification algorithm v0.8.1 |
| 98780c9 | step 4: Euler classification and outputs |
| 4074199 | step 4b specification v0.9 |
| 2fc4d16 | step 4b: IBP-matrix confirmation and outputs |

## Reproducing the results

Requirements: Python 3.13.18; SymPy 1.14.0 for steps 3 and 4 (step 4b uses the standard library only). Run from `t3fund_proca_b412_t7_drafts/` (PowerShell), writing into separate output folders that are created first:

```
New-Item -ItemType Directory -Force out3, out4, out4b
python -I t7a_step3_symbolic_check.py > out3/t7a_step3_output.txt
python -I t7a_step4_classify.py out4
python -I t7a_step4b_ibp_matrix.py t7a_step4_classification.json out4b
```

Compare every output byte by byte with the committed file of the same name. Expected last lines: `ALL_ZERO` for step 3 and `ALL_CHECKS_PASS` for steps 4 and 4b. A re-run from a fresh clone at 2fc4d16 on 2026-11-15 gave byte-identical outputs for all five files. The SHA-256 values of the scripts and JSON outputs are listed in §8 of the manuscript and refer to the LF repository blobs.

## Open items before any publication

1. Physical motivation of B4 in the introduction (marked `AUTHOR TO COMPLETE`).
2. The scope of the README statement "No candidate reproduces the B4 operator structure" relative to this result (marked `AUTHOR TO CONFIRM` in §1.1).
3. Licence, `CITATION.cff` and the untracked files must be resolved before any future Release.
