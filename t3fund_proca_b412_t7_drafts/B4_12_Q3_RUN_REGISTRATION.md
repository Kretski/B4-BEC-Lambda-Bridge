# Registration record — B4.12 Q3 unitarity run (revision 2.1)

This is a registration record, not a report. It records the committed Q3 artefacts and how they were verified. The analytic memorandum on F3/F4 is a separate, later document.

| Item | Value |
|---|---|
| Specification | `t3fund_proca_b412_t7_drafts/B4_12_Q3_UNITARITY_SPEC_v0.3.1.md`, registered at `2a411d0`, SHA-256 `0c02465d77cfe9c5baa603dd13fcab38da17b7ee5e2923d7dcdb0563602f1bd6`, unchanged |
| Conditions of approval | `B4_12_Q3_SPEC_REGISTRATION.md`, `5169551`, unchanged |
| Scope of revision 2 | `B4_12_Q3_CHIRAL_SECTOR_NOTE_v0.2.md`, `ee1e0c7`, SHA-256 `5b8ff30d455e7a5ce8c775ad5a11dd540e4ad1a33efc6eab109618521d7b11c6`, unchanged |
| Script commit | `59dd10d6d16f68b8efd0155991d45e3df7feffb5` (parent `ee1e0c7`) |
| Script | `t3fund_proca_b412_t7_drafts/b412_q3_unitarity.py`, revision 2.1, SHA-256 `8dd9d47566a52ab79465c8861d2927162c2bd261031a141d40a729931936ed7e` |
| Run-output commit | `df4d8b94877505129825949b43e9c89af1154d2a` (parent `59dd10d`) |
| Run outputs | `t3fund_proca_b412_t7_drafts/results_b412_q3/`: 9 files, including `SHA256SUMS.txt` and a byte-identical copy of the script |
| `b412_q3_results.json` | SHA-256 `8bbdd1d00dc7c9859e2f99efd249add5f81fc4565ae8a30365ae8edca093d928` |
| `b412_q3_output.txt` | SHA-256 `97137b7829f760065b5eed6d0931c4cde99f4b8ed8ec17f5e097ee30626a04cd`, identical to the console output |
| Run | 2026-11-15, 18:55:53–18:55:59 (Europe/Sofia), 4.6 s, exit code 0, Python 3.13.18, SymPy 1.14.0; one run |
| Script verdict | `ALL_CHECKS_PASS`, 14 mandatory checks: Q3-C1 to Q3-C9, F1, F2 (symbolic commutator test), F3 transitions, F3 spectator matrix elements, F4 secular term at first order |

## Verification at registration

- The 9 files committed in the run-output directory are byte-identical to the corresponding files in the approved run bundle (SHA-256 `4edb2c3653bba7fdc3493e7e8bef885fac8edf161993e6f1a16e674146f48f17`). The standalone script committed at `59dd10d` is byte-identical to the script copy included in that bundle.
- The run's `SHA256SUMS.txt`, checked on the files as stored in GitHub, gives 8/8 OK.
- From `ee1e0c7` to `df4d8b9` only these 10 files were added: the standalone script and the 9 files in the run-output directory. The specification, the registration record and the chiral-sector note are unchanged.
- The reviewer independently checked the integrity and internal consistency of the run package (2026-11-15, 19:03). This included the SHA-256 of the executed script inside the archive. It is not an independent re-execution of the derivation.

## Revision history (for traceability)

Each revision was statically reviewed before any execution.

| Revision / document | Outcome |
|---|---|
| script revision 1 | static review: not approved. It reported the finding `CHIRAL_COLLINEAR_SPLITTING`, and R_inel⁽²⁾ could not be set to zero. Not run |
| chiral-sector note v0.1 | revision required (F2 operator statements, F3 table, F4 scope, O3 wording) |
| chiral-sector note v0.2 | approved as the scope of revision 2 with seven conditions; committed `ee1e0c7` |
| script revision 2 | static review: not approved. Four targeted fixes were required: F2 naming and scope, the full spectator matrix element, F4 completeness checks, and the spectator relation in F4. Not run |
| script revision 2.1 | approved for execution; one complete run, `ALL_CHECKS_PASS` |

## Scope and status

**Derived by the script's checks:**
- C1–C9:
  - kinematics;
  - the elastic fibre, with 2s_cm norm and phase space 1/(2s_cm) by two methods;
  - the amplitude matrix by two independent routes;
  - channel normalisation;
  - the normalisation convention c = 1, N = 1/(2s_cm), shown to be convention-independent.
- F1: the exact split L_contact = L_RL + L_ch. L_RL has boost weight 0 and L_ch has weights ±4. Parity invariance holds.
- F3: the first-order transition table. The full spectator matrix elements R L → R R R L and R L → R L L L factorise for all 10 nonzero species combinations, with the spectator factor 2|q|·2π·δ(q′ − q).
- F4: for the stated normalised product packet, at first order, all five channels are ESTABLISHED with probabilities of the form C·t², C finite and nonzero.

**Reading of ESTABLISHED:** nonzero is established for generic values of the parameters. The channels whose coefficients are ∝ κ² (or κ⁴/p²) vanish at κ = 0, and the Λ channel vanishes at Λ = 0. ESTABLISHED is not a statement that every channel is nonzero at every admissible specialisation of the parameters.

**What F4 does not establish:**
- the behaviour under the exact interacting evolution;
- higher orders;
- other wave packets;
- non-existence of any asymptotic construction.

**F2:** this is a symbolic momentum-space commutator test, not a full operator construction in the periodic box. The discrete lattice k_n = 2πn/L_box is not imposed. The zero-mode groups that pass or fail this test are **candidates** for an additional chiral sector. It is not shown that they form one in the full theory.

**Full L_contact:** `NO_UNIQUE_CRITERION`, for the present free-Fock setting. Λ_unit is **not computed**. This is not a proof that the theory is non-unitary, nor that no scattering description exists.

**L_RL (conditional only):**
- The elastic matrix is diagonal, and the leading-order estimate is Λ_RL² = min(8/p, 8p/κ²) = 8/max(p, κ²/p) for κ ≠ 0. The value 8/|κ| is never strictly the smallest term.
- The two remaining terms are equal at p = κ²/p, i.e. κ² = p². This is a special point. It is **not** the Q2 rank-one condition κ² = ps. The Q3 contact operator is already built with s = κ²/p, so the equality needs in addition s = p.
- This result is valid only if a prescription of type O1 or O2 is later shown to leave the first-order RL amplitude unchanged. It is not a result for the full contact operator or for the full theory.

**Not decided and not authorised by this record:**
- no choice among O1, O2 and O3;
- no one-loop computation;
- no change of the registered specification;
- no closing of Q3 or of B4.12;
- no T7a registration;
- no Zenodo.

**Next steps, each to be reviewed separately:**
1. an analytic memorandum on F3/F4;
2. a separate task for the zero modes: a full discrete operator construction in the periodic box;
3. only then, a decision on the asymptotic setting.

**Q3: OPEN. B4.12: OPEN.**
