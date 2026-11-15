# Registration record — B4.12 Q2 derivation run (revision 7)

This is a registration record, not a report. It records the committed Q2 derivation artefacts and their verification. The scientific Q2 report is a separate, later document.

| Item | Value |
|---|---|
| Specification | `t3fund_proca_b412_t7_drafts/B4_12_Q2_DERIVATION_SPEC_v0.2.md`, registered at `641121c`, SHA-256 `0c4b3f70e400f3cb6e410e74f6bde0ee40662453207e5a6e3a939c61dbfe9f49`, unchanged |
| Script commit | `a7aec474c8390056f52021493efcc7e3f63582f9` (parent `641121c`) |
| Script | `t3fund_proca_b412_t7_drafts/b412_q2_derivation.py`, revision 7, SHA-256 `68b47e477e776b2cb10803e467f6fda2a9b59006fe95a72ac223c5ac473346cf` |
| Run-output commit | `e8852dd635cddb9895e75beb6d4b1b50361e8d45` (parent `a7aec47`) |
| Run outputs | `t3fund_proca_b412_t7_drafts/results_b412_q2/` (9 files, including `SHA256SUMS.txt` and a byte-identical copy of the script) |
| `b412_q2_results.json` | SHA-256 `24dcc0e657c7dff1fb83c53a5ba277513f0bd7aef3fa8f3ed59b7367a9b30dea` |
| `b412_q2_output.txt` | SHA-256 `a92f7c02c1e2e4651e0c843a5a4027b5bfc16bc66901b4441a9312a864b80896` (identical to the console output) |
| Run | 2026-11-15, 17:22:39–17:22:51 (Europe/Sofia), 11.4 s, exit code 0, Python 3.13.18, SymPy 1.14.0 |
| Script verdict | `ALL_CHECKS_PASS` (12 mandatory checks) |

## Verification at registration

- All 10 committed files are byte-identical to the approved bundle (SHA-256 `c1fec73365ff4ba4facf4fba4184fc42a1f51598140499d771a5069ae6436ec6`).
- The run's `SHA256SUMS.txt`, checked on the files as stored in GitHub, gives 8/8 OK.
- From `641121c` to `e8852dd` only these 10 files were added. `evidence/`, the registered specifications, Q2 spec v0.2 and B4.12 v0.4 are unchanged.
- The reviewer independently checked the integrity and internal consistency of the run package (2026-11-15, 17:25). This is not an independent re-execution of the derivation by a second code.

## Revision history (for traceability)

Each revision was reviewed before execution.

| Revision | Outcome |
|---|---|
| 1, 2 | static review, revisions required; not run |
| 3 | run stopped in D1: `euler_equations` drops identically zero equations, so the guard misfired; no result |
| 4 | fixed the guard; stalled in D2(b) (Bareiss determinant); no result |
| 5 | Berkowitz determinant; stalled in D3 (`sp.cancel`); no result |
| 6 | `sp.together`; derivation completed but the JSON write failed (`sympy.Integer`); no result |
| 7 | `int()` on two degrees; complete run, `ALL_CHECKS_PASS` |

Revisions 4–7 changed only error handling, performance and serialisation. The formulas and checks are those of revision 3, which was approved after review.

## Scope and status

- **Derived** (by the script's checks): D1, Q2-C1 to Q2-C6, D2(b) with the heavy branch continued for **generic, non-resonant parameters**. The exceptional set is M² = 3Λq²k², where the heavy root meets the light f branch. D3, including D3 (iii) for the stated scaling transformation.
- **Not derived:** D4(a), reality at finite k/M; D4(c), light–heavy separation. D4(b) is a necessary condition only.
- **D5** classification (ii) holds for generic simple light branches only (c0 ≠ 0, c0² ≠ 1).
- **Matching scope:** a single Proca mediator reproduces the rank-one subfamily of B4, not general B4.
- **Not authorised** by this record:
  - numerical evaluation on registered backgrounds;
  - Q3;
  - closing B4.12;
  - T7a registration;
  - Zenodo.
- **B4.12: OPEN.**
