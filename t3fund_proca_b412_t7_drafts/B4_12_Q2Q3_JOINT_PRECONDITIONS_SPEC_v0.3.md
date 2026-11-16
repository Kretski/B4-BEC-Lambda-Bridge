# B4.12 — Joint Q2/Q3 check: specification of the preconditions (proposal v0.3)

**Status: PROPOSAL FOR REVIEW.** This document assesses whether the joint Q2/Q3 check defined in the registered B4.12 v0.4 §5 and §7 can be carried out now. It introduces no calculation, no code, no numerical threshold and no change to v0.4. It selects none of O1, O2, O3. **Q3 and B4.12 remain OPEN.** T7a remains blocked. No Zenodo action.

Scope correction of the review of 2026-11-16 (approved before drafting): the joint check is foreseen by the registered v0.4. It is **not** declared impossible in principle, and no amendment of v0.4 is required merely because Q2 and Q3 use different set-ups. The vacuum case (Part D) is a statement about the registered coverage of Q2, not a proof of incompatibility.

Changes from v0.1 (review of 2026-11-16, REVISION REQUIRED): F1, F3, G1, G3 and H now separate a restricted diagnostic comparison from the fulfilment of the closure criteria of v0.4 §7. Parts A–E are unchanged.

Changes from v0.2 (review of 2026-11-16, REVISION REQUIRED): F4 states that ε_abs at c₀² = 1 is not interchangeable with ε_EFT for v0.4 §7(b); G1 states the O1/O2 analyses as a route to investigate whether a full-operator unitarity scale can be established. No other change.

**Principle used in Parts F–H.** Two kinds of comparison are kept apart:
- a **restricted diagnostic comparison** of scales, with a conditional or partial input stated as such;
- the **closure comparison** of v0.4 §7, which needs both (a) and (b) as defined there.

A restricted diagnostic comparison can be registered and reported. It never fulfils §7(a) or §7(b), and it never permits `B4.12-CLOSED-POSITIVE`.

## Sources

| Source | Identifier |
|---|---|
| B4.12 proposal v0.4 (§3.4, §4, §5, §7, §8) | `74c2a77` |
| Q2 derivation report v0.2 | `36f17ce` |
| Q2 run outputs (`results_b412_q2/b412_q2_results.json`) | `e8852dd` |
| Q3 F3/F4 memorandum v0.2 | `f7962fe` |
| Q3 run registration | `17f53fe` |
| Consolidated Q3 status report v0.2 | `95babd8` |

## Part A — The registered definition (v0.4 §5 and §7)

**v0.4 §5, joint Q2/Q3 check:**
- **Common variable.** Q2 concerns fluctuations about a background; Q3 concerns vacuum scattering. A common variable is a device for comparing scales, not an identification of the two set-ups. The relation k = √s_cm/2 is admissible only for symmetric kinematics with a fixed frame and a stated identification.
- **Which scale binds first.** The dimensional estimate (Λ_unit ~ M/g versus k ~ M) is motivation only. The conclusion must use the derived formulas.
- **Uniqueness.** Λ_unit is unique only when the following are all fixed:
  - the eigenvalue (§4 item 3);
  - the sector and the constant (§4 item 6);
  - the vacuum (§4 item 4);
  - the EFT contact (§4 item 2).

**v0.4 §7, end state B4.12-CLOSED-POSITIVE.** It requires both of the following:
- (a) a computed unitarity that is acceptable in some regime;
- (b) that this regime satisfies the pre-defined validity criterion of the approximation, from Q2.

Finding a regime with acceptable unitarity alone is not sufficient. The choice of end state is the reviewer's.

**v0.4 §3.4, constraint on Q2 thresholds.** Grids and thresholds are not fixed in v0.4. First the following must be established:
- the exact formula for ε_EFT;
- the range in which the fluctuation dispersion is physically admissible;
- the relation to the C4 and E3 backgrounds.

τ and ε_max are not proposed, and no threshold is tuned after a result is seen.

The joint check is therefore a comparison of the Q3 unitarity scale with the Q2 validity regime, under a stated common variable. It is not a search for a common physical mechanism.

## Part B — Status of the Q3 input

- The registered outcome for the full L_contact in the free-Fock vacuum setting is **NO_UNIQUE_CRITERION** (`17f53fe`, `f7962fe`, `95babd8` §4.5). Λ_unit for the full L_contact is **not computed**.
- The only scale available is the conditional estimate Λ_RL² = 8/max(p, κ²/p), for κ ≠ 0, of the elastic RL sector. It is valid only if a prescription of type O1 or O2 is later shown to leave the first-order RL amplitude unchanged. It is not a Λ_unit for the full contact operator.
- Of the four uniqueness conditions of v0.4 §5, three are fixed in the registered Q3 setting: the vacuum, the EFT contact, and the sector with its constant. The fourth is open. With the chiral collinear channels present, no unique eigenvalue criterion exists for the full operator.

**Consequence for §7(a):** no computed unitarity statement for the full L_contact exists that could serve as input (a).

## Part C — Status of the Q2 input

- **Formula.** Q2 derived the leading error for generic simple light branches: ε_EFT = |2c₀c₁/(c₀² − 1)|·u + O(u²), with u = (k/M)². It is of order (k/M)² with a coefficient that is not identically zero (`36f17ce` §2.8 and §6).
- **Outcome (ii).** ε_EFT depends on k/M, so a ρ threshold is not a standalone criterion (`36f17ce` §6).
- **Coverage.** D5 holds only for generic simple light branches, c₀ ≠ 0 and c₀² ≠ 1. The special branches and degenerate roots are excluded (`36f17ce` §3 item 3).
  - At c₀² = 1, ε_EFT is undefined, and only ε_abs is reported (`36f17ce` §2.8, D3 special cases).
- **Admissibility.** D4(a), the reality of the light branches at finite k/M, is not derived. D4(c), the light–heavy separation, is not derived. D4(b) is a necessary condition only (`36f17ce` §4).
- **Thresholds.** No value of ε_max, τ or a grid was given (`36f17ce` §6). By v0.4 §3.4 none may be tuned after a result is seen.
- **Backgrounds.** No evaluation on the C4, E3 or population backgrounds (`36f17ce` §6).

**Consequence for §7(b):** Q2 provides a derived error formula, but no operational PASS/FAIL validity criterion. Two of the three prerequisites that v0.4 §3.4 lists before thresholds can be fixed are not established: the admissible range (D4(a), D4(c)) and the relation to the C4 and E3 backgrounds.

## Part D — The vacuum case c₀² = 1

To be stated formally in this document, by direct substitution into the registered expression. No script and no new computation are needed.

- The registered M⁻² light polynomial F0(c) (`e8852dd`, key `F0_c_(M^-2 light polynomial)`), evaluated at zero background w = q = v = r = 0, reduces to

      F0(c) = −p (c² − 1)².

  This was confirmed by the reviewer on 2026-11-16.
- For p > 0 the light roots are the double roots c₀ = ±1.
- These roots lie outside the coverage of D5 (degenerate roots, c₀² = 1). There ε_EFT is undefined and only ε_abs is defined (`36f17ce` §2.8, D3 special cases).

**Reading.** At the vacuum point, the registered Q2 coverage does not supply ε_EFT. This is a statement about the registered Q2 coverage. It does **not** show that a joint comparison is impossible. By v0.4 §5 the common variable compares scales and does not identify the two set-ups.

## Part E — Three set-ups that must not be mixed

| Set-up | Content | Status |
|---|---|---|
| E1. Q3 vacuum scattering | free massless Fock states about φ̄ = 0, L_contact and the separate Λ term | done and registered; full operator NO_UNIQUE_CRITERION |
| E2. Q2 fluctuations about a background | constant-gradient backgrounds (ν = 0), tree level, generic simple light branches | done and registered within its scope |
| E3. Q3 scattering about the B4 background | v0.4 §4 item 4: a separate sub-case, with cubic vertices J₁J₂/M² and broken Lorentz invariance | foreseen by v0.4; not specified, not performed |

A joint statement must name which set-ups it compares and must state the common variable and its identification (v0.4 §5). Results from E1 must not be carried over to E3, and results from E2 must not be carried over to the vacuum point.

## Part F — Conditions for executability of the joint comparison

The joint comparison of v0.4 §5 and §7 can be carried out only when all of the following hold:

- **F1 — Q3 quantity.** For the closure comparison: a unitarity scale for the full L_contact, with all four uniqueness conditions of v0.4 §5 fixed, for the set-up being compared. Only this can serve as input to §7(a).
  - A restricted quantity, such as the conditional Λ_RL of the elastic RL sector, can be the object of a restricted diagnostic comparison only, after a registered decision that states its conditional status.
  - Λ_RL is not a unitarity bound for the full L_contact. A comparison based on it does not fulfil §7(a).
- **F2 — Q2 criterion.** An operational validity criterion: ε_max and τ fixed in a registered document **before** any comparison is evaluated (v0.4 §3.4).
- **F3 — Q2 admissibility.** For the closure comparison: the prerequisites of v0.4 §3.4 for fixing thresholds are established. These are the physically admissible range (D4(a), D4(c)) and the relation to the C4 and E3 backgrounds.
  - A reviewer decision to proceed without them, with a stated restriction, permits only a restricted diagnostic comparison.
  - Without an established admissible range it cannot be concluded that a regime of acceptable unitarity satisfies the validity criterion of the approximation. Such a decision therefore does not fulfil §7(b) and does not permit `B4.12-CLOSED-POSITIVE`.
- **F4 — Common set-up.** A stated pairing of set-ups from Part E, with the reason why their scales are comparable. Two possibilities are:
  - E1 with Q2 at the vacuum point, which requires a treatment of c₀² = 1, where only ε_abs exists;
  - E2 with E3.

  For the E1–vacuum pairing, ε_abs is the only registered Q2 error measure at c₀² = 1. It is not automatically interchangeable with ε_EFT for the purpose of v0.4 §7(b). Any use of ε_abs in a closure comparison requires a separately justified and registered validity criterion. Without that, the vacuum pairing is diagnostic only.
- **F5 — Common variable.** A stated common variable and its identification, e.g. k = √s_cm/2 for symmetric kinematics in a fixed frame (v0.4 §5). The dimensional estimate is not a substitute.

## Part G — Blocking gaps and possible remedies

| Gap | Blocks | Possible remedies (none selected here) |
|---|---|---|
| G1. No Λ_unit for the full L_contact (Part B) | F1, §7(a) | for §7(a): a specification of the analyses that O1 or O2 would require (e.g. [Q_R, L_RL] for O2), as a route to investigate whether a unitarity scale for the full L_contact can be established. Diagnostic only, not a remedy for §7(a): a registered decision to compare the conditional Λ_RL, with its conditional status carried into every conclusion |
| G2. No operational Q2 criterion (Part C) | F2, §7(b) | a separate registered decision fixing ε_max and τ before any comparison, under v0.4 §3.4 |
| G3. Q2 admissibility not derived (D4(a), D4(c)); no relation to C4/E3 | F3, §7(b) | for §7(b): a separate specification for D4(a) and D4(c), and for the relation to the C4 and E3 backgrounds. Diagnostic only, not a remedy for §7(b): a reviewer decision to proceed without them, with a stated restriction |
| G4. Set-ups not paired (Parts D and E) | F4 | a specification of the Q3 sub-case about the background (E3, v0.4 §4 item 4); or a specification of the Q2 treatment at c₀² = 1 based on ε_abs |
| G5. Common variable not stated | F5 | to be fixed in whichever specification performs the comparison |

The remedies are listed for completeness. Their order and selection are the reviewer's decision. None is authorised by this document.

## Part H — Determination (proposed, subject to review)

**The joint Q2/Q3 check is foreseen in the registered v0.4, but it cannot yet be carried out to completion.** The two reasons are:
- the full Q3 unitarity estimate is missing (G1);
- the operational validity criterion from Q2 is not fixed (G2, with G3).

Pairing the set-ups (G4) and fixing the common variable (G5) are further prerequisites.

The check is not declared impossible in principle. No amendment of v0.4 is required by the difference between the Q2 and Q3 set-ups. The end state of B4.12 cannot be evaluated against v0.4 §7 at present. B4.12 therefore stays OPEN by v0.4 §7, for lack of computation.

A restricted diagnostic comparison (with the conditional Λ_RL, or without the Q2 admissibility prerequisites) may be authorised later by a separate decision. Its result would be reported as diagnostic. It would not fulfil §7(a) or §7(b), and it could not lead to `B4.12-CLOSED-POSITIVE`. Registration of this document records only the finding that the joint check is not yet executable to completion. It authorises no calculation and does not select which of G1–G5 is addressed first.

## Out of scope

- any calculation, code or numerical evaluation;
- any numerical threshold, grid or value of ε_max or τ;
- any change to v0.4 or to a registered specification, result or report;
- any choice among O1, O2 and O3;
- any statement on a common physical mechanism of Q2 and Q3;
- Λ_unit, T7a, Zenodo;
- closing Q3 or B4.12.

## Requested from the reviewer

1. Approval of Parts A–H as the precondition assessment.
2. Confirmation of the formal statement in Part D.
3. A decision on whether this document is registered as it stands, as the record that the joint check is not yet executable.
4. A decision, separately and later, on which gap (G1–G5) is addressed first.

**Q3: OPEN. B4.12: OPEN.**
