# B4.12 — Consolidated Q3 Status Report (proposal v0.2)

**Status: PROPOSAL FOR REVIEW.** This report consolidates the registered Q2 and Q3 derivations and the accepted zero-mode analysis. It introduces no new calculations, changes no registered result, and selects no asymptotic prescription. **Q3 and B4.12 remain OPEN.** O1, O2 and O3 are not selected. T7a remains blocked. No Zenodo action is proposed.

Changes from v0.1, after technical review:
- §3.1: the Proca matching is stated in the registered notation of B4.12 v0.4 and of the Q2 report. This was the blocking item.
- §2: the Q2 report and run-registration identifiers were added, and the chiral-sector note file name was corrected.
- §4.2: the κ² = p² point is clarified.
- §4.5: the registered reading of ESTABLISHED was added.
- §3.2: one added interpretation is marked as such.
- §5.2: the mode-level comparison targets were corrected.
- §5.3 and §5.4: the two editorial clarifications to the zero-mode memorandum were recorded.
- §8.3: the equal-label sentence now uses the approved wording.
- Formulas are written in Unicode.

---

## 1. Purpose and scope

This report consolidates the current mathematical and computational status of the B4.12 programme, with emphasis on Q3.

The work has progressed through three related but distinct investigations:

1. **Q2 — Proca matching and perturbative EFT structure:** whether the tree-level Proca construction reproduces the relevant B4 operator structure, and what can be established about the light and heavy branches.
2. **Q3 — Scattering and perturbative unitarity around free massless Fock states:** the elastic block, the first-order inelastic channels, and the secular behaviour of collinear splitting.
3. **Zero modes — Discrete operator construction in a periodic box:** the operator algebra, the zero-mode contributions, equal-free-energy matrix elements, cutoff comparisons and Wick remainders.

These investigations are related through the B4 operator structure, but they do not establish one another's conclusions automatically. In particular:

- Q2 does not establish Q3 unitarity.
- The F3/F4 calculation does not establish a unique asymptotic-state prescription.
- The zero-mode calculation does not establish a physical instability or a failure of unitarity.
- Passing the registered symbolic checks is not equivalent to an independent proof of the complete physical theory.

The purpose here is to state what is established, under which assumptions, and what remains unresolved.

## 2. Registered basis

The report is based on the following registered artefacts.

| Workstream | Registered artefact | Identifier |
|---|---|---|
| B4.12 base | `B4_12_PROCA_PERTURBATIVITY_AND_UV_CLOSURE_v0.4.md` and `B4_12_v0.4_REGISTRATION.md` | `74c2a77`; `1e0ad68` |
| Q2 specification | `B4_12_Q2_DERIVATION_SPEC_v0.2.md` | `641121c` |
| Q2 script | `b412_q2_derivation.py`, revision 7 | `a7aec47` |
| Q2 results | `results_b412_q2/` | `e8852dd` |
| Q2 run registration | `B4_12_Q2_RUN_REGISTRATION.md` | `35465c2` |
| Q2 report | `B4_12_Q2_DERIVATION_REPORT_v0.2.md` | `36f17ce` |
| Q3 specification | `B4_12_Q3_UNITARITY_SPEC_v0.3.1.md` | `2a411d0` |
| Q3 specification registration | `B4_12_Q3_SPEC_REGISTRATION.md` | `5169551` |
| Chiral-sector note | `B4_12_Q3_CHIRAL_SECTOR_NOTE_v0.2.md` | `ee1e0c7` |
| Q3 script and run | `b412_q3_unitarity.py`, revision 2.1; `results_b412_q3/` | `59dd10d`; `df4d8b9` |
| Q3 run registration | `B4_12_Q3_RUN_REGISTRATION.md` | `17f53fe` |
| F3/F4 interpretation | `B4_12_Q3_F3_F4_MEMORANDUM_v0.2.md` | `f7962fe` |
| Zero-mode specification | `B4_12_Q3_ZERO_MODE_SPEC_v0.2.md` | `fcd7ace` |
| Zero-mode interpretation | `B4_12_Q3_ZERO_MODE_MEMORANDUM_v0.1.md` | `1b4fdce` |
| Zero-mode script | `b412_q3_zero_modes.py`, revision 1.1 | `5850414` |
| Zero-mode results | `results_b412_q3_zero_modes_rev1_1/` | `3facc9b` |
| Zero-mode run registration | `B4_12_Q3_ZERO_MODE_RUN_REGISTRATION.md` | `05da826` |

The Q2, Q3 and zero-mode run artefacts have separate registrations. Their results are not interchangeable, and their respective limitations remain binding.

## 3. Q2 — Proca matching and EFT perturbativity

### 3.1 What was derived

The Q2 derivation concerns the registered 1+1-dimensional T3-fund Proca action with metric (−,+), constant-gradient backgrounds (ν = 0), and tree-level elimination of the Proca field.

The exact tree-level kernel is

    K(ω, k) = (1 + k²/M²) / (M² + k² − ω²).

**Notation (B4.12 v0.4).** s denotes only the B4 EFT coefficient. The kinematic (Mandelstam) variable is s_cm.

**Matching.** The registered Proca matching is g_φ = M√p and g_ψ = −Mκ/√p, so g_ψ²/M² = κ²/p. The Proca M⁻² truncation equals the B4 EFT fluctuation operator of the **rank-one subfamily** s = κ²/p identically, compared with B4 at m = 0 (Q2-C3). For general B4 with independent s, the difference is nonzero and proportional to (s − κ²/p). It is confined to the g–g entry.

So g_ψ²/M² equals s only under the rank-one condition κ² = ps. A single Proca mediator does not reproduce general B4 with independent p, s, κ.

The registered Q2 derivation established:
- the H-Q2 quadratic action;
- the kernel identity;
- the rank-one match;
- the decoupling and Schur identities;
- the M²K − 1 identity;
- the agreement of two analytic routes for the first two light-branch expansion coefficients c₁ and c₂.

### 3.2 Light and heavy branches

The heavy Proca branch is continued from the decoupled root ω = ±√(M² + k²). Its continuation is established for generic, non-resonant parameters. The exceptional set is M² = 3Λq²k², where the heavy root meets the light f branch ω² = (1 + 3Λq²)k². The continuation is not established there.

For generic simple light branches, the leading EFT error is ε_EFT = O((k/M)²), with a coefficient that is not identically zero.

The registered D5 classification is outcome **(ii)**, limited to generic simple light branches. By the definition in B4.12 v0.4 §3.4, this means that a ρ threshold is not a standalone criterion. (Interpretive note, not part of the registered text: this does not mean that ρ is physically irrelevant.)

### 3.3 Q2 limitations

The Q2 report explicitly leaves the following unresolved:

- **D4(a):** reality of the light branches at finite k/M is not derived.
- **D4(c):** light–heavy separation is not derived.
- **D4(b):** the reported condition is necessary, not sufficient.
- **D5:** it does not cover degenerate or special light branches, including c₀ = 0 and c₀² = 1.
- **Backgrounds:** no evaluation on the registered C4, E3 or population backgrounds.
- **Q3:** no Q3 unitarity conclusion and no joint Q2/Q3 conclusion follow from Q2.

The twelve passing checks establish the identities evaluated by the script. The reviewer checked the integrity and internal consistency of the run package. This was not an independent re-derivation with a second implementation. The two Q2-C5 routes are distinct methods, but they share the same light relation F(c, u).

**Q2 conclusion:** the specified Proca matching and the generic-branch perturbative analysis are established within their registered scope. General B4 matching, complete spectral admissibility and Q3 unitarity remain open.

## 4. Q3 — Scattering around free massless Fock states

### 4.1 Setting

The registered Q3 analysis uses:

- 1+1 dimensions;
- massless fields f and g;
- the vacuum (f̄ = ḡ = 0);
- free massless Fock asymptotic states;
- first order in the contact couplings for the registered transition analysis.

The primary contact interaction is

    L_contact = (p/2)(f_t f_x)² − κ(f_t f_x)(g_t g_x) + (κ²/2p)(g_t g_x)²,   p > 0.

The Λ term −(Λ/4) f_x⁴ is treated separately (decision D-a). It is not silently incorporated into the primary contact matrix.

### 4.2 Operator decomposition and elastic block

The registered F1 decomposition is L_contact = L_RL + L_ch, with f_t f_x = L² − R² exactly, where

    L_RL = −p R_f²L_f² + κ(R_f²L_g² + R_g²L_f²) − (κ²/p) R_g²L_g²,
    L_ch = (p/2)(R_f⁴ + L_f⁴) − κ(R_f²R_g² + L_f²L_g²) + (κ²/2p)(R_g⁴ + L_g⁴).

The non-chiral part has boost weight 0. The chiral parts have boost weights +4 and −4. The operator is parity invariant.

Only L_RL contributes to the registered first-order elastic block. Its matrix is diagonal in the channel basis, with eigenvalues of 𝒯

    −p s_cm²/4,   κ s_cm²/4 (twice),   −κ² s_cm²/(4p).

The conditional estimate Λ_RL² = 8/max(p, κ²/p), for κ ≠ 0, belongs to the elastic RL sector:
- It is conditional on an O1- or O2-type prescription leaving the first-order RL amplitude unchanged.
- It is not a unitarity threshold for the full L_contact.
- The two terms of max(p, κ²/p) are equal at κ² = p². This is a special point. It is **not** the Q2 rank-one condition κ² = ps.

### 4.3 F3 — first-order transition structure

**Within the registered F3 transition analysis**, the only allowed first-order transitions from |a_R b_L⟩ are:

1. the elastic transition;
2. R → R R R, with an L spectator;
3. L → L L L, with an R spectator.

The full spectator matrix elements factorise for all ten nonzero species combinations. Two routes to the collinear splitting amplitudes agree (route A = route B).

With K = k₁ + k₂ + k₃, the registered vertices are:

| Channel | T₁₃ / (K k₁k₂k₃) |
|---|---:|
| f → f f f | −12p |
| f → f g g | 4κ |
| g → f f g | 4κ |
| g → g g g | −12κ²/p |
| Separate Λ channel: f → f f f | 6Λ |

The comparison between elastic and splitting entries is a comparison of corresponding coupling structures, not an equality of physical amplitudes or probabilities. Within each corresponding structure, the same coupling coefficient appears up to fixed combinatorial factors, e.g. 12p = (p/2)·4!. There is no independent coupling parameter within L_contact that suppresses a chiral vertex relative to its corresponding elastic interaction.

### 4.4 F4 — secular packet probabilities

For the registered packet φ(k) = √(8π/σ)·√k·e^(−k/σ), the first-order packet calculation gives probabilities proportional to t²:

| Channel | P(t)/t² |
|---|---:|
| f → f f f | 9p²σ⁶/(256π²) |
| f → f g g | 3κ²σ⁶/(256π²) |
| g → f f g | 3κ²σ⁶/(256π²) |
| g → g g g | 9κ⁴σ⁶/(256π²p²) |
| Separate Λ channel: f → f f f | 9Λ²σ⁶/(1024π²) |

For the normalised product packet ψ_R × ψ_L, the sum over the final spectator equals 1. So P(R L → R R R L) = P(R → R R R) at first order, for this packet.

The resulting time scales mark where first-order perturbation theory around free Fock states ceases to be self-consistent for the specified packet. They are not physical lifetimes and do not describe the exact time evolution.

The plane-wave T_int and the finite-time packet factor t encode the same secular degeneracy in their respective formulations. Their exact identification has not been established without a common regulator prescription.

### 4.5 What F3/F4 establish, and what they do not

The registered results establish first-order collinear transitions with a spectator, and secular t² probabilities for the specified packet. They show a limitation of the perturbative description around free massless Fock states at long times.

**Reading of ESTABLISHED.** It holds for generic values of the parameters:
- the channels whose coefficients are ∝ κ² (or κ⁴/p²) vanish at κ = 0;
- the Λ channel vanishes at Λ = 0.

ESTABLISHED is not a statement that every channel is nonzero at every admissible specialisation of the parameters.

They do **not** establish:

- non-unitarity of the theory;
- non-existence of an S-matrix or of any scattering description;
- failure of every possible asymptotic construction;
- the exact interacting evolution or the higher-order behaviour;
- the same probability law for arbitrary initial packets;
- a universal EFT energy threshold;
- a value of Λ_unit for the full L_contact.

The registered outcome for the full operator in the present setting remains **NO_UNIQUE_CRITERION**.

## 5. Zero modes — discrete operator construction

### 5.1 Purpose and setting

The zero-mode task extended the earlier symbolic momentum-space F2 commutator test by constructing the operators in a periodic box.

For each field φ ∈ {f, g}, the construction separates the zero-mode pair (q_{φ,0}, π_{φ,0}) from the right- and left-moving oscillators. The zero-mode momenta enter ∂_tφ. The zero modes do not enter ∂_xφ.

The accepted revision 1.1 run passed ZC1–ZC7. This verdict covers those seven checks only.

### 5.2 Operator algebra and continuum comparisons

Within its construction, the run establishes:

- the finite-cutoff mode algebra and the truncated periodic delta;
- the derivative-field decomposition ∂_tφ = π_{φ,0}/L + R_φ + L_φ and ∂_xφ = L_φ − R_φ;
- the field-level comparison with the registered F1 split and the registered F2 zero-mode list;
- the mode-level comparison with the registered elastic block (L_contact and, separately, Λ) and the registered F3 splitting amplitudes;
- the separate conservation of the zero-mode momenta under H₀ + H_int;
- the specified lattice commutators;
- the registered leading-power L-scaling of all 25 box groups.

The mode normalisation map a(k_n) = √(2|k_n|L)·a_n, with (2π)δ ↔ Lδ_nm, reproduces two sets of registered elements:
- the elastic block;
- the R → R R R elements, together with the parity-related L → L L L elements.

These comparisons hold for distinct labels (§8.3).

Along the secondary floor sequences n(L) = ⌊kL/(2π)⌋, the reported error vanishes as L → ∞ and is O(1/L).

The L-scaling statements concern individual matrix elements. They are not statements about transition rates or probabilities. Summed quantities, including the growth in the number of final collinear states, were not calculated in this task.

### 5.3 Commutators and cutoff comparison

The box construction reproduces all 15 rows of the registered symbolic F2 table:

- eight rows have the registered value True;
- seven rows have the registered value False.

The word **confirmed** means that the box construction agrees with the registered True/False value. It does not mean that all 15 groups commute with H₀.

This confirmation holds within the stated construction and for the checked operator groups. It is not a universal statement for arbitrary states and labels. (Editorial clarification to the zero-mode memorandum `1b4fdce`, §2.2.)

For the cutoff-independence test:

- identically vanishing commutators have zero blocks to compare, so these cases are trivial;
- the 11 nonzero commutators have 192 comparable blocks in total;
- between 4 and 48 blocks are compared per group;
- all compared blocks are independent of N_max.

The 11 nonzero commutators are a different grouping from the seven False rows of the F2 table. The two counts should not be compared.

This is block-wise independence only. The image of a sector can depend on N_max, because a larger cutoff admits further output sectors. It is not a statement that whole sectors are cutoff-independent.

### 5.4 Z4 — equal-E₀ first-order matrix elements

The Z4 analysis classifies first-order matrix elements between states with equal free energy in a fixed (π_f, π_g) sector.

It identifies:

- number-conserving bilinears, including f ↔ g species mixing in the π_fπ_g subgroup;
- four number-changing groups:
  - L_contact (3,0,π_f);
  - L_contact (3,0,π_g);
  - L_contact (0,3,π_f);
  - L_contact (0,3,π_g).

These are four groups by their notation, not four separate transitions. (Editorial clarification to the zero-mode memorandum `1b4fdce`, §3.) They contain R → R R and L → L L transitions with collinear momentum conservation.

The matrix elements are nonzero between free states degenerate in E₀. They are proportional to π_{φ,0}/L. They vanish in the zero-mode vacuum limit at fixed π_{φ,0}, and are suppressed by L^(−1/2) relative to the registered elastic and R → R R R elements.

No probability, summed rate or secular coefficient for these π-dependent groups was calculated.

Therefore Z4 is a classification result. It is not evidence of physical decay, instability or unitarity violation.

### 5.5 Z6 — Wick remainders

The zero-mode run distinguishes the two operator structures.

**L_contact.** The chiral monomials have even degree in each family. Their Wick remainders reduce to:
- same-family bilinears;
- c-numbers;
- single fields, whose spatial integral vanishes.

The reported remainders commute with H₀ for the checked cutoffs N_max = 1, 2, 3.

**The separate Λ term.** The chiral decomposition contains the odd-family monomials Λ L_f³R_f and Λ L_fR_f³. One contraction leaves the mixed bilinear ∫:L_f R_f:. Its commutator with H₀ is nonzero for N_max = 1, 2, 3. The (2,2) Λ group behaves like L_contact: its reported remainders commute with H₀.

This is consistent with ZC7 PASS. ZC7 verifies the prescribed Wick-remainder identities and reports their commutators. It does not require every remainder to commute with H₀.

The nonzero mixed Λ remainder concerns the difference between unordered and normal-ordered chiral parts. With H_int normal-ordered by construction, this remainder is not part of H_int. It matters only if the chiral Λ parts are used without normal ordering, for example in a definition of chiral charges.

No renormalisation prescription and no relation between N_max and a physical cutoff has been selected.

## 6. Synthesis: what the combined results say

### 6.1 Established structural results

1. The Proca construction matches the rank-one subfamily of B4 (s = κ²/p) at order M⁻², not general B4.
2. The Q2 derivation establishes the specified generic-branch perturbative identities. It explicitly excludes finite-k/M reality and light–heavy separation.
3. In the registered Q3 free-Fock setting, the chiral part of L_contact produces first-order collinear splitting channels with a spectator.
4. For the specified packet, those channels give first-order secular probabilities proportional to t².
5. The periodic-box construction confirms the registered symbolic operator pattern in its stated scope. It also identifies further equal-E₀ structures at finite zero-mode momentum.
6. The zero-mode Wick analysis distinguishes the L_contact remainders from the mixed Λ-term remainder.

### 6.2 The central unresolved issue

The principal unresolved issue is not whether the registered symbolic checks pass. It is whether the present free massless Fock asymptotic setting is sufficient for a complete scattering and unitarity analysis of the interacting theory.

Each part of the work adds to this question:
- F3/F4 show a secular limitation of first-order perturbation theory around that free-Fock setting.
- Z4 adds degenerate first-order structures at finite zero-mode momentum.
- Z6 shows that the normal-ordering behaviour of the separate Λ term is not identical to that of L_contact.

These results motivate further analysis of the asymptotic construction. They do not determine which construction is correct.

The relationship between the Q2 background-dependent light branches and the Q3 chiral splitting structure has not been established as a common physical mechanism. Both arise from the registered operator framework, but a joint derivation is required before asserting more.

## 7. Open options: O1, O2 and O3

No option is selected in this report.

- **O1 — incorporate the chiral interactions into the asymptotic Hamiltonian.** The F3/F4 secular result motivates investigating this direction. It does not establish that a non-perturbative construction is necessary or sufficient.
- **O2 — reduced S-matrix prescription.** This would require a separate analysis of the relevant commutators, including terms of the type [Q_R, L_RL], and of their perturbative order. These were not calculated in the registered F3/F4 run or in the zero-mode task.
- **O3 — restrict the conclusions to the present free-Fock setting.** This is the conservative conclusion directly supported by the registered calculations. It does not rule out O1 or O2.

The zero-mode calculation adds information that any future prescription must address at finite π_{f,0} and π_{g,0}. It does not choose the prescription.

## 8. Coverage limits binding on this report

The following limits apply to every conclusion above.

1. **F3/F4 setting:** vacuum, massless f and g in 1+1D, free massless Fock asymptotic states, first order in the contact couplings. The packet probabilities are established for the specified packet, not for arbitrary states.
2. **F3 transition statement:** the three allowed first-order transition types are stated within the registered F3 transition analysis from |a_R b_L⟩. They are not a universal statement about all possible states and operators.
3. **Equal final mode labels:** cases with equal final mode labels have not been checked and are not part of the coverage. They need separate checks of the Fock-state normalisations and the symmetry factors.
4. **Route 2:** the explicit finite-matrix verification was performed only for N_max ∈ {1, 2, 3} and initial sectors 0 ≤ N_R, N_L ≤ 2. Wider statements rest on route 1.
5. **Cutoff independence:** established only for the 192 comparable blocks in 11 nonzero commutators. The zero-block comparisons for identically vanishing commutators are trivial. Sector images can depend on N_max.
6. **Z4:** a classification of first-order equal-E₀ matrix elements only. No probability, rate or physical instability follows.
7. **Z6:** the nonzero mixed Λ-remainder commutators are separate from the ZC7 PASS verdict. No renormalisation or physical-cutoff interpretation is supplied.
8. **Q2:** generic simple branches and generic non-resonant heavy-root continuation only. D4(a) and D4(c) remain not derived. D4(b) is necessary only.
9. **Verification:** `ALL_CHECKS_PASS` refers to the mandatory checks of the respective registered scripts (Q2: 12; Q3: 14; zero modes: 7). Package integrity and internal consistency were reviewed. No independent second-code re-derivation of the full calculations was performed.

## 9. What is not concluded

This report does not conclude:

- that the theory is non-unitary;
- that no unitary completion or scattering description exists;
- that O1, O2 or O3 is necessary, sufficient or preferred;
- that the Q2 light-branch results and the Q3 collinear splitting are the same physical effect;
- that the Λ term has a computed unitarity threshold;
- that the full L_contact has a determined Λ_unit;
- that the zero-mode groups have ceased to be candidates;
- that all mode-label configurations or all cutoff sectors have been checked;
- that B4.12 or Q3 is closed;
- that T7a is unblocked;
- that any Zenodo registration is authorised.

## 10. Recommended next steps

The next steps should be decided only after review of this report.

1. Obtain an independent technical review of this consolidation, especially of the separation between the Q2, F3/F4 and zero-mode conclusions.
2. Resolve any wording corrections without changing registered specifications or run outputs.
3. Prepare a separate decision document for the O1/O2/O3 question only once the missing operator and asymptotic analyses have been specified.
4. Keep T7a blocked until its prerequisites are explicitly satisfied.
5. Do not initiate new calculations, alter registered thresholds, or publish a new scientific claim on the basis of this status report alone.

---

**Current status**

- **Q2:** registered derivation complete within its stated scope; limitations remain.
- **Q3:** OPEN.
- **Zero-mode task:** accepted and registered; coverage limits remain.
- **B4.12:** OPEN.
- **O1/O2/O3:** none selected.
- **T7a:** BLOCKED.
- **Zenodo:** no action authorised.

**This document is a status consolidation, not a closure report or a new proof of unitarity.**
