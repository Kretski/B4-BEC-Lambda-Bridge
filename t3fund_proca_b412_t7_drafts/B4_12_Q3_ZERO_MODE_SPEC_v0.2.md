# B4.12 Q3 — Specification: zero modes and discrete operator construction in a periodic box (proposal v0.2)

**Status: PROPOSAL FOR REVIEW. No code has been written and nothing has been computed.**

This specification defines the zero-mode task named as the next step in:
- `B4_12_Q3_RUN_REGISTRATION.md` (`17f53fe`), "Next steps", item 2;
- `B4_12_Q3_F3_F4_MEMORANDUM_v0.2.md` (`f7962fe`), §5.

It does not change any registered artefact:
- the Q3 specification v0.3.1 (`2a411d0`);
- its registration record (`5169551`);
- the chiral-sector note v0.2 (`ee1e0c7`);
- the revision 2.1 script and run (`59dd10d`, `df4d8b9`);
- the run registration (`17f53fe`);
- the memorandum v0.2 (`f7962fe`).

**No asymptotic prescription (O1/O2/O3) is selected or prepared here.** Q3 and B4.12 remain **OPEN**.

**Changes from v0.1** (review of v0.1, sha256 13cb9318…, revision required, items A–H and conditions on D-z1–D-z3):
1. Zero modes are indexed per field throughout: (q_{f,0}, π_{f,0}), (q_{g,0}, π_{g,0}) (A).
2. Momentum conservation is defined by the sign of each operator, and zero modes carry no spatial momentum. The F1 comparison is done at field level and at mode level after an explicit normalisation map (B, G).
3. Z2 requires the two zero-mode commutators separately (C).
4. Route 2 must include every intermediate sector, and truncation may never count as an identity (D).
5. Equal energy in Z4 is defined with the full free Hamiltonian, zero-mode energies included (E).
6. One regulator N_max is used for everything, with a prescribed order of limits. ZC1 is checked at finite N_max and in the distributional limit. Z6 gives the exact N_max dependence (F, D-z2).
7. The external-momentum limit in Z5 is defined by explicit sequences n(L) (H).
8. "Vacuum limit" becomes the "zero-mode vacuum limit" with the reviewer's wording (D-z1).
9. Route 2 is constructively independent of route 1 (D-z3).
10. H₀ is normal-ordered with respect to the oscillator vacuum, so E_c = 0 (minor edit 1.2 of the v0.2 review).
11. In Z3, cutoff independence is tested against the complete set of external and intermediate sectors (minor edit 1.1 of the v0.2 review).

---

## 1. Purpose

The registered F2 result is a **symbolic momentum-space commutator test**. It uses continuous signed momentum symbols and formal momentum conservation. It does not construct the Fock operators in the periodic box and does not impose the lattice k_n = 2πn/L. The zero-mode groups it lists are recorded only as **candidates** for an additional chiral sector.

This task replaces that test by an explicit operator construction in a periodic box. It answers, at operator level and at first order in the couplings, what the zero modes contribute.

## 2. Setting, regulator and definitions

**Box.** Interval of length L > 0 with periodic boundary conditions. Fields f and g are non-compact, so there are no winding sectors. Metric (−,+), plane waves e^{i(kx − ωt)}.

**Regulator (D-z2).**
- A single mode cutoff 1 ≤ |n| ≤ N_max applies to the **field expansions themselves**, and therefore to every composite operator built from them in Z1–Z6. No other regulator is used.
- The zero modes are not affected by N_max.

**Mode expansion, for each φ ∈ {f, g} separately:**
- φ(t, x) = q_{φ,0} + π_{φ,0} t / L + Σ_{1 ≤ |n| ≤ N_max} (2|k_n| L)^{−1/2} [a_{φ,n} e^{i(k_n x − |k_n| t)} + h.c.], with k_n = 2πn/L.
- [a_{φ,n}, a_{φ′,m}†] = δ_{φφ′} δ_nm.
- [q_{φ,0}, π_{φ′,0}] = i δ_{φφ′}.
- Zero-mode operators commute with all oscillators, and the f and g zero modes commute with each other.

**Directions.**
- n > 0 are R oscillators and n < 0 are L oscillators. Energy |k_n|, spatial momentum k_n. **Zero modes carry no spatial momentum**, because q_{φ,0} and π_{φ,0} are x-independent.
- Sector labels: N_R := Σ_φ Σ_{n>0} n·N_{φ,n} ≥ 0 and N_L := Σ_φ Σ_{n<0} |n|·N_{φ,n} ≥ 0, where N_{φ,n} = a_{φ,n}† a_{φ,n}.

**State space.**
- (zero-mode space of f) ⊗ (zero-mode space of g) ⊗ Fock over all oscillators.
- Each zero mode is a free particle with no normalisable ground state. It is described in generalised eigensectors of π_{f,0} and π_{g,0}. These are used only if their conservation is derived (Z2).

**Free Hamiltonian.**
- H₀ = Σ_φ π_{φ,0}²/(2L) + Σ_φ Σ_{1≤|n|≤N_max} |k_n| a_{φ,n}† a_{φ,n}.
- H₀ is normal-ordered with respect to the oscillator vacuum, so its normal-ordering constant is E_c(N_max, L) = 0. This is a convention of the specification.
- On an oscillator Fock state with zero-mode eigenvalues (π_{f,0}, π_{g,0}), the free energy is E₀ = (π_{f,0}² + π_{g,0}²)/(2L) + (2π/L)(N_R + N_L).

**Derivative fields** (derived from the expansion, not assumed):
- ∂_t φ = π_{φ,0}/L + R_φ + L_φ, where R_φ collects the n > 0 oscillators and L_φ the n < 0 ones;
- ∂_x φ = L_φ − R_φ, with no zero mode.

**Interaction.**
- The operators are L_contact (spec v0.3.1 §2) with p > 0 and κ real (spec §1.1), and, separately, the Λ term (decision D-a).
- The interaction Hamiltonian is H_int = −∫₀^L dx :L_int:, valid at first order only (registered B0). Higher-order Legendre terms are out of scope.
- Normal ordering is with respect to the oscillators. The zero-mode factors π_{φ,0}/L commute with the oscillators and are not reordered.

**Momentum conservation (sign convention).**
- For a monomial Π_j a†_{φ_j, m_j} · Π_i a_{φ_i, n_i}, the spatial integral ∫₀^L dx gives L·δ(Σ_j m_j − Σ_i n_i, 0).
- Creation operators enter with +, annihilation operators with −, and zero-mode factors contribute 0.

**Box-to-continuum map (G).**
- The registered continuum convention is relativistic: [a(k), a†(k′)] = 2|k| (2π) δ(k − k′), as in registered C9.
- With a(k_n) := √(2|k_n| L)·a_n, the box commutator becomes [a(k_n), a†(k_m)] = 2|k_n| L δ_nm. This is the map (2π)δ(k − k′) ↔ L δ_nm.
- This map is derived in the script from the two commutators. It is the only map used to compare coefficients with registered continuum results.

**Zero-mode vacuum limit (D-z1).**
- This means L → ∞ at fixed finite π_{f,0} and π_{g,0}, so that the zero-mode density vanishes: π_{φ,0}/L → 0 and π_{φ,0}²/(2L) → 0.
- It does not assert the existence of a normalisable zero-mode vacuum state.

**Order of limits (D-z2).** The steps are taken in this order only:
1. All primary results are exact at finite L and finite N_max.
2. Then N_max → ∞ at fixed L, for each quantity that has such a limit (for example commutators on fixed sectors, and ZC1 in the distributional sense).
3. Then L → ∞ in the zero-mode vacuum limit, with external momenta fixed as in Z5.

The two limits are not assumed to commute, and the reverse order is not used. Quantities without an N_max → ∞ limit, such as the Wick remainders of Z6, are reported at finite N_max and L only. The physical cutoff k_max = 2πN_max/L is given alongside N_max, and no prescription relating them is chosen.

## 3. Questions

**Z1 — Operator form of the interaction.**
- Write H_int as a sum of normal-ordered monomials in {a_{φ,n}, a_{φ,n}†}, π_{f,0} and π_{g,0}, with exact Kronecker momentum conservation as defined in §2 over integer labels.
- Group the monomials by (degree in R oscillators, degree in L oscillators, degree in π_{f,0}, degree in π_{g,0}).
- **Comparison with F1, at two levels:**
  - at field level, the π-free part, written in R_φ and L_φ, must equal the registered F1 split L_contact = L_RL + L_ch monomial by monomial. This level is normalisation-independent.
  - at mode level, its coefficients, mapped by the §2 box-to-continuum map, must equal the registered continuum vertex coefficients.
- List all π-dependent groups with their coefficients. This is the operator version of the registered candidate list.

**Z2 — Zero-mode momenta, separately.**
- Derive whether q_{f,0} and q_{g,0} appear anywhere in H₀ + H_int. They are expected to be absent, because only derivative fields enter, but this must be shown.
- Derive [π_{f,0}, H₀ + H_int] and [π_{g,0}, H₀ + H_int] separately, at first order. If π_{φ,0} enters H_int but q_{φ,0} does not, the commutator can still vanish. The result is reported as derived.
- If both commutators vanish, the problem decomposes into joint (π_{f,0}, π_{g,0}) eigensectors, and π_{φ,0}/L become c-numbers in each sector. That decomposition is used only if derived.

**Z3 — Exact commutators on the lattice.** On sectors (N_R, N_L), at finite N_max and L, compute exactly:
- [Q_R^osc, H₀], [Q_L^osc, H₀] and [Q_R^osc, Q_L^osc], where Q_R^osc and Q_L^osc are the oscillator-only chiral parts, defined as in the registered F2;
- [G, H₀] for each π-dependent group G of Z1, in each (π_{f,0}, π_{g,0}) sector.

The result for each is either exactly zero or a nonzero operator with its structure stated. Nothing is assumed to vanish. It must also be shown whether the results on a given sector are independent of N_max. Cutoff independence must be tested against the largest mode label required by the complete set of external and intermediate sectors involved in the operator product, not merely by the initial sector (N_R, N_L).

**Z4 — Role of the π-dependent groups.** For each group, within a fixed (π_{f,0}, π_{g,0}) sector, state:
- whether it changes the oscillator particle number;
- whether it conserves N_R and N_L separately or only in total;
- whether it has nonzero first-order matrix elements between free states of **equal full free energy** E₀ (§2).

In the same zero-mode sector the term (π_{f,0}² + π_{g,0}²)/(2L) cancels, and equality reduces to N_R + N_L = N_R′ + N_L′. This cancellation is shown, not assumed. The task classifies only. It does not compute transition probabilities or secular coefficients for these groups.

**Z5 — Zero-mode vacuum limit and external momenta (H).** Following the order of limits in §2:
- Derive the L-scaling of the matrix elements of each group between box-normalised states. The registered leading-power table is a symbolic scaling statement; here it is to be derived.
- **External momenta.**
  - Fixed physical momenta k are realised by sequences of integers n(L) with k_{n(L)} → k.
  - **Primary:** commensurate external momenta, with L restricted to a sequence L_j such that every external k lies exactly on the lattice (k = 2π n_j/L_j with integers n_j → ∞).
  - **Secondary,** reported only as a check of the error: n(L) = ⌊kL/2π⌋, with |k_{n(L)} − k| ≤ 2π/L.
- State which zero-mode contributions vanish relative to the elastic and chiral-splitting elements of the registered run.
- Show that the box elastic block and the R → R R R / L → L L L elements, mapped by §2, reduce to the registered continuum results along the primary sequences. This is a consistency check, not a new result.

**Z6 — Normal-ordering remainders.**
- Derive the exact dependence on N_max and L of the Wick terms separating the chiral parts written without normal ordering from the normal-ordered ones. For example, the closed form of C(N_max, L) = Σ_{n=1}^{N_max} k_n/(2L) and of every remainder coefficient.
- Derive the commutators of the remainder terms with H₀.
- No N_max → ∞ limit is taken for divergent quantities, and no renormalisation prescription is chosen.

## 4. Method requirements (D-z3)

Two routes are required. They must be constructively independent.

**Route 1, symbolic mode algebra.**
- Monomials are carried with explicit integer mode labels and field labels.
- Momentum conservation is the Kronecker delta of §2.
- Commutators come from the canonical relations by normal-ordering algebra.

**Route 2, explicit finite matrices, built independently.**
- **Basis.** Occupation-number states over all oscillator modes 1 ≤ |n| ≤ N_max of both fields.
- **Operators.** a_{φ,n} and a_{φ,n}† are matrices defined by their action on this basis (√N factors).
- **Interaction operators** are built directly from the mode expansion:
  - each derivative field is a sum of a and a† matrices with the mode-function factors of §2;
  - each product of field factors is expanded into creation and annihilation parts, with all creation parts placed to the left of all annihilation parts in every term (normal ordering by construction);
  - the x-integral is evaluated by lattice orthogonality.
- **Independence.** Route 2 does **not** take the monomial table, the coefficients or any intermediate expression of route 1 as input.
- **Intermediate sectors (D).** Matrix products must include every intermediate momentum sector reached by the operator factors. The script must determine the set of reachable sectors from the operator content and include all of them. No truncation boundary may be interpreted as an operator identity or as a vanishing commutator. If a sector reached by a product lies outside the constructed basis, the corresponding result is reported as NOT_COVERED, not as zero.
- **Symbols.** p, κ, Λ, L, π_{f,0} and π_{g,0} are kept symbolic. Only the integer labels (N_R, N_L, N_max) are concrete.

**Agreement and coverage.**
- The two routes must agree on every evaluated sector.
- Route 2 is exact on the sectors it covers, but it covers only those sectors and finite N_max. Statements beyond them rest on route 1, and the report states this.

## 5. Mandatory checks and outputs

| # | Check | PASS condition |
|---|---|---|
| ZC1 | Mode algebra | at finite N_max: [φ(x), ∂_tφ(y)] = i δ_L^{(N_max)}(x − y), the truncated periodic delta (1/L) Σ_{\|n\|≤N_max} e^{i k_n (x−y)}, with the n = 0 term supplied by the zero mode. In the distributional limit N_max → ∞, every Fourier coefficient equals 1/L, i.e. i δ_L(x − y). Both levels are checked for f and g separately |
| ZC2 | Derivative fields | ∂_tφ = π_{φ,0}/L + R_φ + L_φ and ∂_xφ = L_φ − R_φ derived from the expansion, for f and g |
| ZC3 | Continuum consistency | π-free part of Z1 = registered F1 split at field level; mode-level coefficients equal the registered continuum vertex coefficients after the §2 map |
| ZC4 | Zero-mode momenta | q_{f,0}, q_{g,0} absent; [π_{f,0}, H₀ + H_int] and [π_{g,0}, H₀ + H_int] derived separately (zero, or the nonzero result reported) |
| ZC5 | Lattice commutators, two routes | route 1 = route 2 for every Z3 commutator on every evaluated sector; intermediate sectors complete; no NOT_COVERED entry counted as PASS |
| ZC6 | Box-to-continuum limit | map derived from the two commutators; registered elastic entries and R → R R R / L → L L L elements recovered along the primary sequences of Z5 |
| ZC7 | Wick remainders | exact closed-form N_max and L dependence of every remainder coefficient, and their commutators with H₀, reported |

Each check has a defined FAIL or NOT_COVERED outcome. `ALL_CHECKS_PASS` requires PASS on all seven and covers only these checks. It is not an independent proof and does not select O1, O2 or O3.

**Outputs:**
- standalone script `b412_q3_zero_modes.py`;
- `b412_q3_zero_modes_output.txt`;
- `b412_q3_zero_modes_results.json`.

They are written to a new folder `results_b412_q3_zero_modes/`. The registered `results_b412_q3/` is not touched.

## 6. Compatibility with the registered Q3 artefacts

- **Specification v0.3.1.**
  - Same operators, same domain (p > 0, κ real), and the Λ term kept separate (D-a).
  - The relativistic normalisation of v0.3.1 §1 is the continuum side of the §2 map.
  - The vacuum setting of v0.3.1 corresponds to the zero-mode vacuum limit (D-z1).
  - The fibre and phase-space results (C1–C9) are not recomputed. They enter only through ZC3 and ZC6.
- **Revision 2.1 and the run registration.**
  - This task tests the registered symbolic F2 table (`commutes_with_H0_by_group`). It does not replace it.
  - If the box construction confirms, refines or contradicts that pattern, the report states so explicitly. The registered result stays as recorded.
  - The status "candidates" changes only in a later report, after review.
- **Memorandum v0.2.** The conclusions A1–A5 are not revisited here. The zero modes are the open item named there.

## 7. Out of scope

- Any choice among O1, O2 and O3. Any construction of asymptotic states or Møller operators. Any degenerate perturbation theory in the collinear subspaces.
- Second-order commutators such as [Q_R, L_RL] (the O2 analysis), higher-order Legendre terms, and one-loop computations.
- Transition probabilities or secular coefficients for the π-dependent groups. Z4 classifies only.
- **The limit L → ∞ at fixed π_{φ,0}/L = w_φ ≠ 0.** This is a homogeneous background ∂_t φ = w_φ, i.e. the setting of the B4 background and of Q2, not the Q3 vacuum. It is named and not analysed.
- Any renormalisation prescription for the Wick remainders, and any relation between N_max and a physical cutoff.
- Λ_unit for the full L_contact, and any change of the conditional L_RL result.
- Numerical evaluation on registered parameter points.
- Changes to registered documents. Closing Q3 or B4.12. T7a. Zenodo.

## 8. Decisions

| Item | Status in v0.2 |
|---|---|
| D-z1 | Accepted by the reviewer; wording "zero-mode vacuum limit" adopted (§2) |
| D-z2 | One cutoff N_max on the field expansions; prescribed order of limits; no N_max → ∞ for divergent quantities (§2, Z6) |
| D-z3 | Two routes with constructive independence and complete intermediate sectors (§4) |

No further decisions are requested. Execution remains forbidden until this specification is approved and registered and a script has passed static review.
