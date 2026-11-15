# Registration record — B4.12 Q3 unitarity specification v0.3.1

This is a registration record, not a specification. It records the registration of the Q3 specification v0.3.1 and the conditions attached to its approval, without changing the registered file. The header of the registered file reads "PROPOSAL FOR REVIEW", its status before registration; this record supersedes that status line. The file itself is unchanged, so its approved hash is preserved.

| Item | Value |
|---|---|
| Registered document | `t3fund_proca_b412_t7_drafts/B4_12_Q3_UNITARITY_SPEC_v0.3.1.md` |
| Registration commit | `2a411d0` |
| Parent commit | `36f17ce` |
| SHA-256 | `0c02465d77cfe9c5baa603dd13fcab38da17b7ee5e2923d7dcdb0563602f1bd6` |
| Size | 22 447 bytes (LF, ending with a newline) |
| Approval and registration date | 2026-11-15 |
| Review decision | SCIENTIFIC CONTENT APPROVED, with the conditions below |
| Earlier versions (not committed) | v0.1 (`fccc578c…`), v0.2 (`7ee5acff…`), v0.3 (`25df9c77…`); each returned for revision or conditionally approved |

**Verification at registration**
- The file in GitHub at `2a411d0` is byte-identical to the approved file, after pull and `cmp`.
- The copy checked during review had SHA-256 `478b66db4374dc107c0fc9c65f6504e8643114c6e34ad079673303dd489fce88` (22 669 bytes). It is the registered file with CRLF line endings and without the final newline. This was reproduced exactly from the original, so the difference is in line endings only, not in content. The registered file is the original.

## Conditions of approval

These conditions bind the implementation and the later report. They do not change the registered specification.

1. **Registered file.** Only the original `.md` file with the SHA-256 above is the Q3 specification.
2. **§4 step 3.2.**
   - Both relations, −i(A − A†) = c·A†A + R_inel/c and −i(A⁽²⁾ − A⁽²⁾†) = c·a², must be derived for the actual operator L_contact. Reproducing the expected formulas is not enough.
   - The order of R_inel must be shown in the coupling counting: two four-point vertices for 2 → 4, an amplitude of second order, and a contribution to R_inel of fourth order.
   - Any further admissible channels or interactions found must be included and reported.
3. **§4 step 2, formal definition of the elastic sector.**
   - The states are δ-normalised, so the implementation must define the projector P and the elastic sector rigorously. This is done either with momentum wave packets or by decomposing the state space into kinematic sectors.
   - The finite-dimensional channel matrix (over species) must not be confused with the full space of continuous momentum states.
   - This definition must be carried into the relation between matrix elements, phase space and R_inel.
4. **§4 step 3.3 and §5.**
   - The exact result is the bound on the Hermitian part of the exact elastic amplitude, ‖H(A)‖_op ≤ 1/|c|.
   - Its application to the tree-level eigenvalues, |λ_i(a)| ≲ 1/|c|, is a leading-order estimate valid up to perturbative accuracy. It is not an exact constraint on the full amplitude.
   - Λ_unit is reported with this qualifier. NO_UNIQUE_CRITERION remains a valid outcome.

## Status at registration

- **Q2:** completed and registered (`641121c`–`36f17ce`), unchanged.
- **Q3:** specification registered. Execution is **not** yet authorised.
- **Next steps:**
  1. write `b412_q3_unitarity.py` according to the registered specification and the conditions above;
  2. static review of the script;
  3. execution only after approval of the script.
- **B4.12: OPEN.** No conclusion about the full Proca theory, or about the validity of B4 outside the defined contact EFT regime.
- **T7a:** registration remains blocked until a final decision on B4.12.
- **No Zenodo release.**
