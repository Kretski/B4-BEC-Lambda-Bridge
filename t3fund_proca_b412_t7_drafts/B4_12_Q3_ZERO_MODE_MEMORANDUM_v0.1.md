# B4.12 Q3 — Analytic memorandum on the zero-mode run (proposal v0.1)

**Status: PROPOSAL FOR REVIEW. Interpretation of an accepted run only. No new computation, no code change, no option selected.** Q3 and B4.12 remain **OPEN**. O1, O2 and O3 are not selected. T7a remains blocked. No commit or Zenodo action is implied.

## 0. Inputs

| Item | Identifier |
|---|---|
| Specification | `B4_12_Q3_ZERO_MODE_SPEC_v0.2.md`, commit `fcd7ace`, sha256 `2a198f7f1e74360ca909f2e1b4d0a0db826acf6e0c19494491d6cb06f962657a` |
| Script | `b412_q3_zero_modes.py` revision 1.1, sha256 `1e2f51ade19aecbca331727bda72ce48f36b9f4044e305c47663f3465212e509` (66 936 B), approved for execution 2026-11-15 |
| Run | `results_b412_q3_zero_modes_rev1_1/`, 2026-11-15 19:56:17–19:59:19 (Europe/Sofia), 181.4 s, exit code 0, Python 3.13.18, SymPy 1.14.0 |
| Results JSON | `b412_q3_zero_modes_results.json`, sha256 `5057c0cdd4fd36eec117d5a121f89d234c6f8503f014b50a7c0750870fa66683` (111 184 B) |
| Output | `b412_q3_zero_modes_output.txt`, sha256 `c876e412e5340fd43d4073a175bf40972bfcfbeba8d701e6146f8e281f2fdcc7` (3 441 B) |
| Archive | `results_b412_q3_zero_modes_rev1_1.tar.gz`, sha256 `efce73f3aafd98272fc97941ac72903e1f7df449f50f3cf10b2de49e67103f63` (31 094 B) |
| Review of the run | reviewer decision of 2026-11-15: run accepted for scientific analysis; no rerun, no revision 1.2 |

The script and the run are not yet committed. Their commit identifiers will be added when they are registered.

**Verdict of the run:** ZC1–ZC7 all PASS (`ALL_CHECKS_PASS`). This covers the seven listed checks only. It confirms the execution of the registered checks. It is not an independent re-derivation of every result and not an independent proof, and it selects none of O1, O2, O3.

## 1. Coverage and limits (stated first, binding on every section below)

1. **Equal final mode labels.** Cases with equal final mode labels have not been checked and are not part of the coverage. They need separate checks of the Fock-state normalisations and the symmetry factors. The mode-level comparisons (ZC3, ZC6) assume distinct labels for the splitting elements.
2. **Route 2** was evaluated only for N_max ∈ {1, 2, 3} and initial sectors 0 ≤ N_R, N_L ≤ 2. Within that range it is exact.
3. **Wider statements** (general N_max, general labels) rest on route 1, not on a universal finite-matrix check.
4. **Cutoff independence** is established only for the comparable blocks, under the rule of §2.2. It is not a statement about whole sectors.
5. **Z4** is a classification of first-order matrix elements. It does not settle the asymptotic construction.
6. **The nonzero Λ-term commutators of Z6** are a separate result (§4). They are not covered by the statement "ZC7 PASS".

## 2. Operator algebra in the box (Z1–Z3, Z5)

### 2.1 Confirmed structure

- **Derivative fields (ZC2).** From the mode expansion: ∂_tφ = π_{φ,0}/L + R_φ + L_φ and ∂_xφ = L_φ − R_φ for f and g. Mode coefficients are ±i√(π|n|)/L. The zero mode enters ∂_t only.
- **Mode algebra (ZC1).** [φ(x), ∂_tφ(y)] at finite N_max is i times the truncated periodic delta. Every Fourier coefficient with |m| ≤ N_max is i/L, and the m = 0 term comes from the zero mode. The f–g cross commutator is zero. Each coefficient is N_max-independent once N_max ≥ |m|.
- **Field level (Z1, ZC3).** The π-free part equals the registered F1 split L_contact = L_RL + L_ch monomial by monomial, for both L_contact and the Λ term. The π-dependent groups, aggregated by total zero-mode degree, equal the registered symbolic F2 list. The Λ term has no π-dependent group (∂_x has no zero mode).
- **Mode level (ZC3, ZC6).** The map a(k_n) = √(2|k_n|L)·a_n is derived from the two commutators. With (2π)δ ↔ Lδ_nm it reproduces two things exactly, for symbolic labels:
  - the registered elastic 4 × 4 block (L_contact and, separately, Λ);
  - the registered R → RRR elements and, by parity, the L → LLL elements.

  Along the primary commensurate sequences the agreement is exact at every L_j. Along the secondary floor sequences n(L) = ⌊kL/2π⌋ the error vanishes as L → ∞ and is O(1/L) (the limits of L × error are finite and listed in the JSON).
- **Zero-mode momenta (Z2, ZC4).** q_{f,0} and q_{g,0} are absent from H₀ + H_int. Hence [π_{f,0}, H₀ + H_int] = 0 and [π_{g,0}, H₀ + H_int] = 0 separately. The problem decomposes into joint (π_{f,0}, π_{g,0}) eigensectors, in which π_{φ,0}/L are c-numbers.

### 2.2 Lattice commutators (Z3, ZC5)

- [Q_R^osc, H₀] = 0, [Q_L^osc, H₀] = 0 and [Q_R^osc, Q_L^osc] = 0 hold exactly as operators at N_max = 1, 2, 3 (route 1). They agree with route 2 on every evaluated sector, and no NOT_COVERED entry occurred.
- For every group G, [G, H₀] = −ΔE(M)·M monomial by monomial. Here ΔE is the oscillator energy created minus the energy annihilated.
- **Comparison with the registered symbolic F2 table.** All 15 rows of `commutes_with_H0_by_group` are reproduced by the box construction. That is an agreement between the registered and the new result, not a statement that all 15 groups commute with H₀. Eight rows are True and seven are False:
  - L_contact (1,1,π), (1,2,π), (2,1,π) and (2,2,0);
  - Λ (1,3), (2,2) and (3,1).

  (Seven box sub-groups with π-split degrees map to the three π-carrying L_contact rows among these. With L_contact (2,2,0) and the three Λ groups they make the 11 nonzero commutators of the cutoff test below.) The registered F2 test was a symbolic momentum-space test. The box construction **confirms** its pattern at operator level on the lattice.
- **Cutoff independence (rule).** A block (in, out) is compared between N_max and N_max′ > N_max only if N_max ≥ N_req. Here N_req = max(level(in), level(out), the largest level of the intermediate sectors reached at N_max′), with level = max(N_R, N_L). The outcome:
  - For commutators that vanish identically, 0 blocks are compared. Independence there is trivial.
  - For the 11 nonzero commutators, 4 to 48 blocks are compared per group, 192 in total. All are independent of N_max.
  - **The image of a sector does depend on N_max** for these groups. Example: the (1,1,π) groups create pairs a†_{φ,n} a†_{φ′,−n} for every n ≤ N_max, so new output sectors appear as the cutoff grows. Independence holds block by block, not sector by sector.

### 2.3 L-scaling in the zero-mode vacuum limit (Z5)

- For box-normalised states, fixed physical momenta and fixed π_{φ,0}, a matrix element of a group with n_osc oscillator factors and n_π zero-mode factors scales as L^{1 − n_osc/2 − n_π}. All 25 box groups agree with the registered leading-power table.
- Relative to the elastic and the chiral-splitting elements (L^{−1}), the π-dependent groups are suppressed by L^{−1/2} (n_osc = 3, n_π = 1) or L^{−1} (n_osc = 2, n_π = 2).
- **Limitation.** This is a statement about individual matrix elements only. It is not a statement about rates or probabilities. For example, the number of final states of a collinear splitting grows with L, and summed quantities are not computed here (spec §7).

## 3. Z4 — first-order elements between states of equal E₀

- **Cancellation shown.** In a fixed (π_f, π_g) sector, E₀(out) − E₀(in) = (2π/L)(ΔN_R + ΔN_L). Momentum conservation gives ΔN_R = ΔN_L for every monomial, checked for every monomial. Equal E₀ is therefore ΔN_R = ΔN_L = 0.
- **Groups without equal-E₀ elements:** (1,1,π), (1,2,π), (2,1,π). Every monomial has ΔN_R = ΔN_L = n ≠ 0, which means pair-type creation or annihilation across R and L.
- **Number-conserving groups with equal-E₀ elements:** (2,0,π²) and (0,2,π²), in their (π_f², π_fπ_g, π_g²) sub-groups. In each π sector they are bilinears a†_{φ,n} a_{φ′,n}.
  - For φ = φ′ the bilinear is a rescaling of the free one-chirality Hamiltonian. Example: −(p/2)(π_f/L)² Σ_{n>0} k_n N_{f,n}, with element −πpπ_f²/L³ on |f,+1⟩.
  - The mixed sub-group (π_fπ_g) contains a†_{f,n} a_{g,n}. That is **species mixing between f and g at equal energy**, e.g. the monomial a†_{f,+1} a_{g,+1}.
- **Number-changing groups with equal-E₀ elements (finding of the run):**
  - L_contact (3,0,π_f), (3,0,π_g): R → RR;
  - L_contact (0,3,π_f), (0,3,π_g): L → LL.

  Examples, with collinear momentum conservation:
  - |f,+2⟩ → |f,+1²⟩ with element −6iπ^{3/2} p π_f/L³;
  - |g,+2⟩ → |f,+1²⟩ with element 2iπ^{3/2} κ π_g/L³.

**What Z4 does not establish.** These are nonzero first-order matrix elements between free states that are degenerate in E₀. They are not evidence of an instability, a physical decay or a violation of unitarity. Three facts bound their weight:
- they are proportional to π_{φ,0}/L, so they vanish in the zero-mode vacuum limit at fixed π (§2.3);
- they are suppressed by L^{−1/2} relative to the registered elastic and R → RRR elements;
- no probability, secular coefficient or summed rate was computed (spec §7).

Nor does the classification by itself call for O1, O2 or O3. It adds a new degenerate structure (1 → 2 collinear, and f ↔ g mixing) that any future asymptotic construction would have to treat at finite π. That consequence is noted here, not analysed.

## 4. Z6 — Wick remainders: L_contact versus the Λ term

- **Closed form.** C(N_max, L) = Σ_{n=1}^{N_max} k_n/(2L) = πN_max(N_max+1)/(2L²), derived from the mode coefficients. It is the same for R and L. Cross-family contractions vanish because R and L use disjoint modes and f and g modes commute. The unordered = normal-ordered + remainder identity was verified explicitly for every chiral monomial at N_max = 1, 2. No N_max → ∞ limit is taken and no renormalisation is chosen.
- **L_contact.** Every chiral monomial has even degree in each family (R_f, R_g, L_f, L_g). Its remainders are therefore only:
  - same-family bilinears ∫:X²: = Σ |k_n| a†a, i.e. rescalings of one-chirality free Hamiltonians;
  - c-numbers;
  - single fields, whose integral vanishes.

  All of these commute with H₀ at N_max = 1, 2, 3. Example: from (p/2)R_f⁴ in L_ch, the coefficient of ∫:R_f²: is 3πpN_max(N_max+1)/(2L²), and the c-number term is proportional to N_max²(N_max+1)²/L⁴.
- **Λ term (D-a), separate result.** −(Λ/4)(L_f − R_f)⁴ contains the odd-family monomials Λ L_f³R_f (group (1,3)) and Λ L_fR_f³ (group (3,1)). One contraction leaves ∫:L_f R_f:, with coefficient 3πΛN_max(N_max+1)/(2L²) from each group.
  - Momentum conservation makes ∫:L_f R_f: a pure pair operator, a†_{f,n} a†_{f,−n} + h.c.
  - **Its commutator with H₀ is nonzero at N_max = 1, 2, 3.** This is consistent with the creation and annihilation of a pair of opposite chirality and equal |k|.
  - The (2,2) Λ group behaves like L_contact: its remainders commute with H₀.

  Therefore the statement "the Wick remainders commute with H₀" holds for L_contact in the checked range. **It does not hold for every Λ-term remainder.**
- **Consistency with ZC7.** ZC7 PASS requires the closed forms, their explicit verification and the reported commutators. It does not require every remainder to commute with H₀. The nonzero Λ commutators are therefore consistent with ZC7 PASS.
- **Scope.** This concerns the difference between unordered and normal-ordered chiral parts. With H_int normal-ordered by construction (spec §2), the pair remainder is not part of H_int. It matters only if the Λ chiral parts are written without normal ordering, for instance in a definition of chiral charges for the light theory. The Λ term remains outside the primary matrix and outside Λ_unit (D-a).

## 5. Relation to the registered Q3 artefacts

- The registered F2 table stays as recorded. The box construction **confirms** its 15 rows (§2.2).
- The zero-mode groups remain **candidates** for an additional chiral sector. This memorandum does not change that status. The run adds two refinements a later reviewed report may use:
  - the operator form in the box, including the decomposition into π eigensectors;
  - the equal-E₀ structure of Z4.
- The memorandum v0.2 (`f7962fe`) conclusions A1–A5 are not revisited. The NO_UNIQUE_CRITERION outcome for the full L_contact and the conditional L_RL result are unchanged. Λ_unit for the full L_contact is not computed.

## 6. Items for the future Q3 status report

1. The six coverage limits of §1, with the equal-label sentence verbatim.
2. Cutoff independence stated block-wise (192 comparable blocks, 11 groups), with the N_max-dependent sector image named explicitly.
3. Z4: the four number-changing equal-E₀ groups and the f ↔ g mixing sub-groups, as classification only.
4. Z6: the L_contact / Λ distinction as a separate result.
5. The reviewer's clarification on E3 ("within the registered F3 transition analysis"), deferred earlier.

## 7. Not concluded here

- any physical instability, decay or violation of unitarity;
- any need for, or choice of, a new asymptotic construction (O1, O2, O3);
- any probability, rate or secular coefficient for the π-dependent groups;
- the limit at fixed π/L = w (the background setting of Q2);
- any renormalisation of the Wick remainders, or any relation between N_max and a physical cutoff;
- Λ_unit for the full L_contact;
- closing Q3 or B4.12; T7a; Zenodo.
