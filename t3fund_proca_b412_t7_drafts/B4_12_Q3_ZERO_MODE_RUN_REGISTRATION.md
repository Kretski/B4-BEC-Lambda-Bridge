# Registration record — B4.12 Q3 zero-mode run (revision 1.1)

This is a registration record, not a report. It records the committed zero-mode artefacts and how they were verified. Their interpretation is the separate memorandum `B4_12_Q3_ZERO_MODE_MEMORANDUM_v0.1.md` (`1b4fdce`).

| Item | Value |
|---|---|
| Specification | `t3fund_proca_b412_t7_drafts/B4_12_Q3_ZERO_MODE_SPEC_v0.2.md`, registered at `fcd7ace`, SHA-256 `2a198f7f1e74360ca909f2e1b4d0a0db826acf6e0c19494491d6cb06f962657a`, unchanged |
| Script commit | `58504143b6a681d6a101d27b434e764ae24dd911` (parent `1b4fdce`) |
| Script | `t3fund_proca_b412_t7_drafts/b412_q3_zero_modes.py`, revision 1.1, SHA-256 `1e2f51ade19aecbca331727bda72ce48f36b9f4044e305c47663f3465212e509` (66 936 B) |
| Run-output commit | `3facc9b8942ba6dfa8336755a0045236512e67d5` (parent `5850414`) |
| Run outputs | `t3fund_proca_b412_t7_drafts/results_b412_q3_zero_modes_rev1_1/`: 9 files, including `SHA256SUMS.txt` and a byte-identical copy of the script |
| `b412_q3_zero_modes_results.json` | SHA-256 `5057c0cdd4fd36eec117d5a121f89d234c6f8503f014b50a7c0750870fa66683` (111 184 B) |
| `b412_q3_zero_modes_output.txt` | SHA-256 `c876e412e5340fd43d4073a175bf40972bfcfbeba8d701e6146f8e281f2fdcc7` (3 441 B), identical to the console output |
| Run bundle (not committed) | `results_b412_q3_zero_modes_rev1_1.tar.gz`, SHA-256 `efce73f3aafd98272fc97941ac72903e1f7df449f50f3cf10b2de49e67103f63` (31 094 B) |
| Run | 2026-11-15, 19:56:17–19:59:19 (Europe/Sofia), 181.4 s, exit code 0, Python 3.13.18, SymPy 1.14.0; one run |
| Command | `python3 -I b412_q3_zero_modes.py results_b412_q3_zero_modes_rev1_1` |
| Inputs | none (symbolic only; registered expressions of run `df4d8b9` embedded in the script) |
| Script verdict | `ALL_CHECKS_PASS`, 7 mandatory checks: ZC1 mode algebra, ZC2 derivative fields, ZC3 continuum consistency, ZC4 zero-mode momenta, ZC5 lattice commutators by two routes, ZC6 box-to-continuum limit, ZC7 Wick remainders |
| Memorandum | `B4_12_Q3_ZERO_MODE_MEMORANDUM_v0.1.md`, `1b4fdce`, SHA-256 `cd0a444c15cf80ea12ad21a01e1f4d4fcc5c2d6fa8e1dee3bfe801fd682a3864` |

**Output folder name.** Spec §5 names the output folder `results_b412_q3_zero_modes/`. The run was executed with the folder name given in the reviewer's authorisation, `results_b412_q3_zero_modes_rev1_1/`. That name is kept so that the path matches `RUN_CONTEXT.txt` and the executed command. The contents are exactly the outputs required by §5, and the registered `results_b412_q3/` is untouched.

## Verification at registration

- The 9 files committed in the run-output directory are byte-identical to the corresponding files in the run bundle (SHA-256 `efce73f3…`). The standalone script committed at `5850414` is byte-identical to the script copy in that bundle.
- The run's `SHA256SUMS.txt`, checked on the files as stored in GitHub, gives 8/8 OK.
- From `1b4fdce` to `3facc9b` only these 10 files were added: the standalone script and the 9 files in the run-output directory. No registered file was changed. The two unrelated untracked files (`tests/test_dimensions.py`, `theory/operator_matching.py`) were not committed.
- The reviewer independently checked the integrity and internal consistency of the run package (2026-11-15, 20:02):
  - the SHA-256 of the JSON, the output and the bundle;
  - the executed script inside the bundle (66 936 B, approved SHA-256);
  - the internal `SHA256SUMS.txt`;
  - exit code 0.

  This is not an independent re-execution or re-derivation.

## Revision history (for traceability)

Each revision was statically reviewed before any execution.

| Revision / document | Outcome |
|---|---|
| specification v0.1 | revision required (items A–H, conditions on D-z1–D-z3) |
| specification v0.2 | approved and registered, `fcd7ace` |
| script revision 1 (SHA-256 `3280f9b0…`) | static review: one mandatory fix. `scal_ok` (Z5 L-scaling versus the registered table) was computed but not part of the ZC6 PASS condition. Not run |
| script revision 1.1 | the fix only (ZC6 requires `scal_ok`); approved for execution; one complete run, `ALL_CHECKS_PASS` |
| run package | accepted for scientific analysis; no rerun, no revision 1.2 |
| memorandum v0.1 | approved and committed, `1b4fdce` |

## Scope and status

**`ALL_CHECKS_PASS`** covers ZC1–ZC7 only. It is not an independent proof and does not select O1, O2 or O3.

**Coverage limits:**
1. Cases with equal final mode labels have not been checked and are not part of the coverage. They need separate checks of the Fock-state normalisations and the symmetry factors.
2. Route 2 (explicit finite matrices) was evaluated only for N_max ∈ {1, 2, 3} and initial sectors 0 ≤ N_R, N_L ≤ 2. Wider statements rest on route 1.
3. Cutoff independence holds only for the comparable blocks: 192 blocks in 11 nonzero commutators; the comparison is trivial (0 blocks) where the commutator vanishes identically. The image of a sector depends on N_max for the nonzero commutators.
4. ZC7 PASS does not mean that every Wick remainder commutes with H₀. The Λ-term remainders ∫:L_f R_f: do not (memorandum §4).
5. Z4 is a classification of first-order matrix elements, not a statement about instability, decay, unitarity or rates.

**Not decided and not authorised by this record:**
- no choice among O1, O2 and O3;
- no change of the candidate status of the zero-mode groups;
- no change of the registered specification, of the registered F2 table or of any registered Q3 result;
- no Λ_unit for the full L_contact;
- no closing of Q3 or of B4.12;
- no T7a registration;
- no Zenodo.

**Next step, to be reviewed separately:** the future Q3 status report. It will carry the new commit identifiers, the two editorial clarifications to the memorandum, and the deferred clarification on E3.

**Q3: OPEN. B4.12: OPEN.**
