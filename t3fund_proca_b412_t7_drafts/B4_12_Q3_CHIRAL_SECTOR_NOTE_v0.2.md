# B4.12 Q3 — Technical note on the chiral sector (`CHIRAL_COLLINEAR_SPLITTING`), proposal v0.2

**Status: PROPOSAL FOR REVIEW. No computation has been run.** This note responds to the static review of `b412_q3_unitarity.py` revision 1 (sha256 d07f6b94…, "NOT APPROVED FOR EXECUTION").

It does **not** change the registered Q3 specification (`B4_12_Q3_UNITARITY_SPEC_v0.3.1.md`, 2a411d0) or its registration record (5169551). Registration condition 2 already requires additional channels to be reported. This note:
- states the finding precisely;
- separates what is algebra from what still has to be derived;
- sets out the options for a decision.

B4.12 remains **OPEN**.

**Changes from v0.1** (review of v0.1, sha256 87714fae…, accepted as the basis for revision 2 after four clarifications):
1. F2 now separates the classical invariance of the integral from the two operator commutation relations. H0 + Q_R + Q_L is stated as a proposal for an asymptotic Hamiltonian, not as the exact interacting Hamiltonian.
2. F3 now has a table: initial state, final state, active part of the operator, and the reason it is allowed or forbidden.
3. F4 now states secularity and the failure of the standard free-Fock perturbative construction. The categorical claim that no S-matrix exists is withdrawn, and so are the claims of unbounded t² growth for every wave packet and of an infinite R_inel⁽²⁾.
4. O3 is stated as the result of the present free asymptotic setting, not as a proof that no unitary theory is possible. §3, §4 and §5 are reworded accordingly.

---

## 1. The finding

On shell in the interaction picture, f = f_R(x − t) + f_L(x + t). Define R_f = ∂_t f_R and L_f = ∂_t f_L, so that ∂_x f_R = −R_f and ∂_x f_L = +L_f, and the same for g. Then, exactly:

  f_t f_x = L_f² − R_f²,  g_t g_x = L_g² − R_g².

There are no mixed R_f L_f terms. Substituting into L_contact splits it into two parts:

**Non-chiral part** (degree (2, 2) in R/L):

  L_RL = −p R_f²L_f² + κ (R_f²L_g² + R_g²L_f²) − (κ²/p) R_g²L_g².

**Chiral part** (degrees (4, 0) and (0, 4)):

  L_ch = (p/2)(R_f⁴ + L_f⁴) − κ (R_f²R_g² + L_f²L_g²) + (κ²/2p)(R_g⁴ + L_g⁴).

This is elementary algebra. Revision 1 already computes it (`chiral_decomposition`), and revision 2 will print both parts explicitly.

## 2. Structural facts

Each item is marked as either algebra already in revision 1, or a derivation to be added in revision 2.

**(F1) Boost behaviour (algebra; to be printed).**
- Under a boost with rapidity η, R² and L² scale with opposite factors e^{±2η}, and d²x is invariant.
- So L_RL is boost-invariant. This agrees with the revision-1 design result that the RL tree amplitude depends on s_cm only.
- L_ch has boost weight ±4 and is therefore the **Lorentz-non-invariant** part of the contact operator.
- This follows from J = g f_t f_x being a single component (T_{01}-like): (T_{01})² = T_{++}² + T_{−−}² − 2T_{++}T_{−−}.

**(F2) Time independence of the chiral integrals. These are three separate statements.**
- **(F2a) Classical; algebra already in revision 1.** Every term of L_ch depends on x and t only through x − t (or only through x + t). For a function G with sufficient decay, ∂_t G(x − t) = −∂_x G(x − t), so ∫dx G(x − t) is independent of t under the free evolution.
- **(F2b) Operator statement, to be derived: [Q_R, H0] = 0** for Q_R := −∫dx :L_ch^R:. The derivation must fix:
  - the operator definition, with normal ordering of the R-sector modes;
  - the boundary conditions, either infinite line with decaying wave packets or a box with stated boundary conditions;
  - the treatment of the zero mode k = 0.

  It must then show the commutator from these definitions, not by appeal to (F2a).
- **(F2c) Operator statement, to be derived: [Q_R, Q_L] = 0.** It is not assumed on the grounds that the classical fields depend on different light-cone coordinates. It must be shown from the mode algebra, with the same boundary and zero-mode conventions, including whether R and L zero modes are shared.
- **Status of H0 + Q_R + Q_L.** This is only a *proposal* for the asymptotic Hamiltonian of the relevant approximation (option O1), not the exact interacting Hamiltonian. The Legendre transformation of a Lagrangian with a nonlinear dependence on time derivatives produces further terms at higher order. Revision 1 (B0) verifies H_int = −L_int only to first order.

**(F3) Complete first-order transition structure from |a_R b_L⟩ (derivation to be added).**
- For any set of free massless particles, E − P = 2 Σ_{L particles} E_i. A single R particle has E − P = 0, so any final state of a transition starting from one R particle consists of R particles only.
- With one four-point vertex (V = 1, six external legs impossible for a connected diagram), the first-order transitions from |a_R b_L⟩ are therefore:
  1. elastic R L → R L, through L_RL only, as computed in C4;
  2. R → R R R with an L spectator, through the (4, 0) part of L_ch, and L → L L L with an R spectator, through the (0, 4) part;
  3. nothing else. R L → R R and R L → L L are excluded by C1. Transitions from L_RL that create extra particles need zero-momentum R or L quanta and are excluded except at the IR end point (Coleman limitation, spec §6).
- The enumeration covers every possible matrix element of a single normal-ordered four-field operator between the initial state |a_R b_L⟩ and any final state. With one R and one L particle initially, a vertex can annihilate n_a = 0, 1 or 2 of them (at most one per direction) and create 4 − n_a particles. The created set must carry the four-momentum of the annihilated set.

| Initial | Annihilated (n_a) | Active part of the operator | Created set | Final state | Verdict and reason |
|---|---|---|---|---|---|
| a_R b_L | none (0) | any | 4 particles with E = P = 0 | — | forbidden: massless particles with E = 0 (zero-mode end point only) |
| a_R b_L | R (1) | (4,0) of L_ch | R R R | R R R + L spectator | **allowed**: E − P = 0 requires all created particles to be R |
| a_R b_L | R (1) | (2,2), (3,1), (1,3) | contains at least one L | — | forbidden: E − P of the created set must be 0, but each L adds 2E_L > 0 |
| a_R b_L | L (1) | (0,4) of L_ch | L L L | R spectator + L L L | **allowed**: mirror of the R case |
| a_R b_L | L (1) | (2,2), (3,1), (1,3) | contains at least one R | — | forbidden: mirror argument (E + P of the created set must be 0) |
| a_R b_L | R and L (2) | (2,2) of L_RL | one R + one L | c_R d_L, forward momenta (C1) | **allowed**: elastic channel, computed in C4 |
| a_R b_L | R and L (2) | (3,1) | two R | R R | forbidden by C1 (needs E_L = 0); present only in the Λ term |
| a_R b_L | R and L (2) | (1,3) | two L | L L | forbidden by C1 (needs E_R = 0); present only in the Λ term |
| a_R b_L | R and L (2) | (4,0) or (0,4) | — | — | impossible: these parts cannot annihilate an L (resp. R) particle |

- Revision 2 must generate this table programmatically:
  - by enumerating, for every monomial of L_contact (and separately of the Λ term), every assignment of its fields to annihilation and creation;
  - by checking each against energy–momentum conservation.

  It must also compute the elements of the allowed splitting rows by route A and route B. This answers the reviewer's point that C4 covers only the four two-particle channels, not the full first-order transition structure.

**(F4) Collinear degeneracy and secularity (derivation to be added).**
- All R states with the same total momentum k have the same free energy, E = k. A one-particle state |k⟩ and the three-particle states |k₁k₂k₃⟩ with k₁ + k₂ + k₃ = k are therefore **exactly degenerate** under H0.
- Q_R couples them with a finite fixed-time matrix element ⟨k₁k₂k₃|Q_R|k⟩. In revision 1 this is `T_13`, which is ∝ p.
- In first-order time-dependent perturbation theory around free Fock states, the transition amplitude therefore contains a **secular term proportional to t**, with no oscillating phase. This is the origin of the δ(0) in the S-matrix element.
- **What this establishes:** the standard perturbative construction of the S-matrix around free massless Fock states has a secular problem in the chiral sectors, so that construction cannot be completed for L_contact as it stands.
- **What it does not establish:**
  - that no asymptotic construction exists;
  - that the probability for every normalised wave packet grows without bound, for example like t², under the exact interacting evolution.

  The latter would need control of the integration over the continuum of final momenta, of the wave-packet normalisation, and of higher orders. Revision 2 derives the secular term and its coefficient for an explicit normalised wave packet at first order only, and states this limitation.
- Classical counterpart, as a remark only: the chiral quartic gives an amplitude-dependent propagation speed for R waves, with no dispersion to balance it.

## 3. Consequences for the registered derivation

- **§4.3.2.**
  - If the splitting states belong to the open final states, the second-order relation becomes −i(A⁽²⁾ − A⁽²⁾†) = c·a² + R_inel⁽²⁾/c, as the reviewer stated. Setting R_inel⁽²⁾ = 0 is not justified by the connected counting.
  - By (F4), within the free-Fock perturbative construction R_inel⁽²⁾ inherits the secular term, so it is not a well-defined second-order quantity in that construction.
- **§4.3.3.**
  - The exact bound ‖H(A)‖_op ≤ 1/|c| presupposes a well-defined elastic block S_el.
  - With the free massless Fock construction blocked in the chiral sectors (F4), S_el is not yet defined for the full L_contact within that construction. The bound, and with it the leading-order criterion, then holds only conditionally on a prescription.
- **What is not affected.**
  - The first-order elastic R L → R L amplitude between bare Fock states. It receives no first-order chiral contribution (F3).
  - The tree-level matrix a computed from L_RL.

## 4. Options

### O1 — Chiral interactions in the asymptotic Hamiltonian

Propose H0′ = H0 + Q_R + Q_L as the asymptotic Hamiltonian (only once (F2b) and (F2c) are established; not the exact Hamiltonian, see F2), and define in/out states as eigenstates of H0′.

What must be derived:
- the stationary states of P_R + Q_R in each degenerate collinear subspace;
- that Møller operators with respect to H0′ exist;
- that the remaining interaction −∫L_RL produces no open channel at second order.

Known risk:
- Because of the exact degeneracy (F4), degenerate perturbation theory makes the zeroth-order asymptotic states eigenvectors of Q_R within each subspace. In general these are **O(1) superpositions of different particle numbers**, not small corrections to Fock states.
- It is then not guaranteed that the tree-level RL amplitude between bare states is the leading amplitude between the new asymptotic states. That would have to be derived.
- Also, H_int = −L_int contains −(p/2)R⁴ for p > 0, so H_R + Q_R is unbounded below. The construction can only hold within the EFT range.

### O2 — Separate collinear sectors and a reduced S-matrix

Define a reduced S-matrix for the non-chiral interaction and treat the chiral sectors separately.

What must be derived:
- factorisation of the evolution into a chiral part and an RL part, at least to second order. This requires control of [Q_R, L_RL]-type terms, which do not vanish in general.
- that the reduced S-matrix is unitary on a well-defined space.

Without this, O2 amounts to discarding L_ch by hand.

### O3 — Limitation of the present free asymptotic setting

Establish (F4) within its stated scope: for L_contact with p ≠ 0 or κ ≠ 0, the standard perturbative construction around free massless Fock states is secular in the chiral sectors and cannot be completed. Then:
- **Full L_contact:** the standard perturbative unitarity argument cannot be closed in the free-Fock setting until the chiral collinear sector is resolved. Within Q3 the result for the full operator is **NO_UNIQUE_CRITERION**, as foreseen in spec §4 step 3.3. This is the result **of the present setting**. It is not a proof that a unitary theory or a scattering description is impossible.
- **Non-chiral part L_RL:** a well-defined, boost-invariant leading-order result exists and can be reported. It is reported as a **conditional** result, valid only if a prescription of type O1 or O2 is later shown to leave the first-order RL amplitude unchanged.
- **§4.3.2:** recorded as not closable for the full operator.

## 5. Proposal for the decision

1. **Revision 2 of the script derives facts, not choices.** It adds (F1)–(F4) as derivations:
   - the explicit L_RL / L_ch split and their boost weights;
   - (F2a) as algebra, and (F2b) and (F2c) as operator statements with stated boundary, zero-mode and normal-ordering conventions;
   - the programmatic first-order transition table of F3, including the R → R R R and L → L L L elements by both routes;
   - the exact degeneracy and the secular (t-linear) first-order amplitude for an explicit normalised wave packet, at first order only.
2. **No option is selected in the script.** On the basis of these derivations the reviewer decides between O1, O2 and O3. O1 and O2 are each separate tasks larger than Q3 and would need their own specification. O3 stays within the registered specification, since NO_UNIQUE_CRITERION is a foreseen outcome, provided (F4) is established in its stated scope. The reviewer's recommendation is not to select O1 or O2 within Q3.
3. **Normalisation constant, as decided by the reviewer.** Revision 2 fixes c = 1 and N = 1/(2 s_cm), so a = 𝒯/(2 s_cm) and the leading-order criterion is |λ_i(a)| ≲ 1. Revision 2 shows explicitly that the criterion is independent of this convention: with c general, only c·N = 1/(2 s_cm) is fixed, and |λ_i(a)| ≤ 1/c ⇔ |λ_i(𝒯)| ≤ 2 s_cm.
4. **No one-loop computation**, as decided by the reviewer.
5. **Until the decision:**
   - §4.3.2 is reported only as the general relation including R_inel⁽²⁾/c. The form c·a² is not claimed for the actual operator.
   - Λ_unit is **not computed** for the full L_contact.
   - For L_RL, revision 2 may derive the elastic matrix and its leading-order criterion as a **conditional** result, explicitly labelled as not a result for the full contact operator.

## 6. Not included

- No run of any script.
- No change of the registered specification or registration record.
- No conclusion on the full Proca theory.
- No change to Q2, T7a, or Zenodo.
