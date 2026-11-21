# Λ: map of independent extraction routes — v0.4

**Status:** a map and planning note. No computation and no run.
- It is registered together with `N0_HAMILTONIAN_QNM_INVENTORY_v0.2.md`, in the commit that adds both, on top of f18065a.
- **The registration fixes documentary results, i.e. what the files say and what was found in them. It does not raise any of them to independent physical confirmation.**

**Supersedes:** v0.3 (6fbf9972…), v0.2 (b6a3e41b…) and v0.1 (050f6ec9…). None of them is registered.

**Changes from v0.3:**
- N0 is completed: inventory v0.2, from files read and checked by SHA-256.
- The factor-2 Hamiltonian row is withdrawn.
- Finding A0 is added.
- BEC is reconciled.
- The Steinhauer table is added.
- The photon bounds are updated.
- The priorities are revised (review of 2026-11-21, 16:30).

**Criterion (unchanged).** Agreement between two formulas from the same theoretical construction is not independent physical support. Support comes only from measured data.

## §0 Decision (project lead, 2026-11-21), unchanged

| direction | status |
|---|---|
| **(B) dispersive Λ**, ω² = c²k²(1 + Λk²) | **primary target** |
| (A) nonlinear Λ (Λ_B4; Q in SOC) | a separate line |
| bridge through c₄ (SOC) | **unverified hypothesis** (§7) |
| one universal Λ across sectors | **not established** and not assumed |

(A) and (B) are not declared to be the same physical quantity.

## §1 Hamiltonian normalisation: canonical; the factor-2 row is withdrawn

- Paper 1 (`vacuum_oscillator_v0.4.tex` l.121–123; preprint v4 l.224, 249, 724): ω² = c²α²k²(1 + Λk²), from the mass shell g^{μν}k_μk_ν + Λk_loc⁴ = 0.
- Paper 2 code (`A1v3_zamo_v04.py` l.188–213): **H = ½[g^{μν}k_μk_ν + Λ(k_loc²)²]**. The factor ½ multiplies both terms, and H = 0 gives the canonical law.
- H_V4 = αk_loc√(1 + Λk_loc²) is canonical. "1 + 2Λ…" appears there only in derived quantities: the numerator of v_g, and the einbein variable change k_loc = w(1 + 2Λw²).
- **Withdrawn:** the v0.3 row "mass-shell H = ½g kk + Λ_H k⁴ ⇒ Λ = 2Λ_H". That form does not occur in the files.
- **The only source of "1 + 2Λk²"** is the obsolete "Kretski" normalisation, ω² = c²k² + 2Λc²k⁴, i.e. **Λ_old = Λ/2**.
  - It is used in `lambda_wave_equation.py` (outside the git repo) and in a few Lambda-model scripts.
  - It is labelled obsolete in `PAPER3_VALIDATION_STATUS.md` l.62, l.187.

Details and file hashes: N0 v0.2, §1 and §3.

## §2 Finding A0: the original nonlinear term does not generate (B)

**Source.** `THEORY_STATUS_GW.md` §2 and §9 of the Lambda-model repository (96987ae7…).
- Paper 1 derives its mass shell from the action term **(Λ/2)(h^{μν}∂_μφ ∂_νφ)²**.
- That term is **quartic in the field**: it drops out at linear order, and its mass dimension is 8, not 6. **It does not produce the linear k⁴ dispersion** (item A0).
- A **linear**, dimension-6 operator built from the Laplacian of the preferred slices does give ω² = k²(1 + Λk²). On Kerr its eikonal limit reproduces the Paper 2 Hamiltonian (A1). That repository lists it as a **candidate**.

**Reading in the terms of this map.** The original construction mixed two mechanisms:
- a type-(A) nonlinear self-interaction (the same class as Λ_B4 (∂φ)⁴);
- a type-(B) linear dispersive correction, which (A) does not generate at linear order.

A type-(B) Λ needs a linear higher-derivative operator. The dimension-6 candidate is **not** an established physical mechanism.

**Scope of A0 (required wording).** A0 is the result of a **theoretical audit of the existing construction**. It is **not** an experimental refutation of any realisation of (B), and it does not exclude a type-(B) Λ produced by another mechanism.

**Lab remark (not a result).** In a BEC, the k⁴ term of the Bogoliubov spectrum, ħ²k⁴/(4m²), comes from the kinetic term (quantum pressure). The self-interaction sets only c.

## §3 Canonical convention and the BEC reconciliation

**Canonical:** ω = ck√(1 + Λk²) ⇔ ω² = c²k²(1 + Λk²), [Λ] = m².

| written form | canonical Λ |
|---|---|
| ω² = c²k²(1 + 2Λ″k²) (obsolete "Kretski") | Λ = 2Λ″ |
| ω = ck(1 + Λ′k²) | Λ = 2Λ′ |
| ω = c₀k + βk³ | Λ = 2β/c₀ |
| LVK α = 4 | Λ = A₄(ħc)² |

**BEC: one physical value, two symbol conventions.**

| source | normalisation | ξ | stated | canonical Λ |
|---|---|---|---|---|
| `lambda_wave_equation.py` | obsolete | ħ/(√2 m c_s) | Λ_old = ξ²/4 | **ħ²/(4m²c_s²)** |
| `example_BEC.py`, `paper3_final.tex` | canonical | ħ/(m c_s) | Λ = ξ²/4 | **ħ²/(4m²c_s²)** |

**Λ_BEC = ħ²/(4m²c_s²).** "ξ²/4" is meaningful only together with its ξ, and this map uses ξ = ħ/(m c_s). This line is kept explicitly so that the apparent conflict does not come back.

## §4 Lab BEC: Steinhauer null test, numbers as recorded

**Source.** `steinhauer_test_v3.py` and `steinhauer_kmin_scan.csv` (Lambda-model repo). Data: Steinhauer et al., PRL 88, 120407 (2002), 13 points extracted from the vector EPS of Fig. 3.
- H0 is the local-density Bogoliubov spectrum. H1 is H0 with the k⁴ coefficient scaled by (1 + g_extra). So **g_extra = ΔΛ_lab/Λ_Bog**.
- A Gaussian prior is used on μ/h = 1.91 ± 0.09 kHz, with a profile likelihood. The k_min scan uses a null test of 1000 pseudo-experiments.
- **The numbers are quoted from the recorded output, not re-run in this repo.**

| k_min [µm⁻¹] | N | g_extra ± err | 95% interval | Δχ² | p_emp | H0 goodness of fit p |
|---|---|---|---|---|---|---|
| 0.0 | 13 | 0.031 ± 0.023 | [−0.01, 0.07] | 1.79 | 0.19 | 5.7e-5 |
| 0.5 | 12 | 0.032 ± 0.023 | [−0.01, 0.07] | 1.88 | 0.18 | 3.5e-5 |
| 0.8 | 11 | 0.026 ± 0.023 | [−0.01, 0.07] | 1.23 | 0.29 | 2.8e-3 |
| 1.2 | 10 | 0.019 ± 0.023 | [−0.02, 0.06] | 0.62 | 0.42 | 0.039 |
| 1.5 | 9 | 0.012 ± 0.024 | [−0.03, 0.05] | 0.28 | 0.59 | 0.28 |
| 2.0 | 8 | 0.019 ± 0.024 | [−0.02, 0.06] | 0.62 | 0.44 | 0.37 |

**Result.** There is **no statistically established excess** over standard Bogoliubov theory. ΔΛ_lab/Λ_Bog is consistent with 0, and the 95% upper end is about 0.05–0.07. The interval edges are quantised by a profile grid of step 0.01.

**Limit.** H0 is a poor fit when the low-k points are included (p = 5.7e-5 at k_min = 0) and becomes acceptable only from k_min ≈ 1.5 µm⁻¹. **The bound is therefore not independent of how the low-k points are treated.** It is model-dependent there (local density, points flagged as upper-bound errors). The source script states that g_extra is not an estimate of the GW / strong-field Λ.

## §5 Hamiltonian / eikonal / QNM route: developed, with no observational channel

- **What exists:** a scalar-field and ray framework, with photon-ring and shadow responses, Paper 2 Table 3 / Eq. 20, the retrograde-ring D_Λ^eik = 0.0127 Λ̂ l³ ("negligible at natural Λ"), and eikonal QNM estimates.
- **Status, per the Lambda-model repository** (`THEORY_STATUS_GW.md` §9):
  - **A2(a):** with the Boyer–Lindquist / ZAMO foliation the QNM problem is **not well posed**;
  - the Paper 2 QNMs remain eikonal estimates;
  - decision: no new QNM computations.
- **GW link:** none. A tensor equation needs dynamics with an independent coefficient Λ_g. Λ_GW = Λ is a universality assumption.
- **The giant-vortex analogue** is an implementation check (model vs model; `CORRECTION_QNMS_CSV.md`, 19 Oct 2026). EXP C compares with markers digitised from a figure, not with the original data. Either way it tests the vortex's own model, not Λ.
- **Decision (review of 2026-11-21):** this route is **closed for now as a path to independent extraction**. Its theory and code stay in the project history. No new QNM computations are started.

## §6 Photon bounds: the old Fermi number and the newer LHAASO numbers (unconfirmed)

| bound | source as recorded | status here |
|---|---|---|
| Λ_γ ≲ 2e-53 m² | Paper 1, Fermi-LAT GRB 090902B (order of magnitude), via `THEORY_STATUS_GW.md` §4 | the origin of the "1e-53 m²" in the registered `LAMBDA_SECTOR_AND_CONVENTION_MAP.md`; **outdated** |
| 1e-55 to 7e-57 m² | LHAASO GRB 221009A time of flight (Piran & Ofengeim 2024; Yang, Bi & Yin 2024; DisCan), via `notes/SENSITIVITY_MAP.md` | **unconfirmed** (N4) |
| 6e-60 m² | photon decay (review arXiv:2109.09072) | **unconfirmed** (N4) |
| 4e-62 / 1e-64 m² | photon splitting, HAWC / LHAASO (as quoted in arXiv:2109.09072) | **unconfirmed** (N4) |

All of these are superluminal bounds, converted with Λ = (ħc/E_QG,2)². SENSITIVITY_MAP calls them orders of magnitude "to be checked against each source before citation". **None is cited as a result until N4 is done.**

The GW row is unchanged: Λ_GW ∈ [−2.4e-11, +7.4e-12] m², a conversion of the LVK A₄ interval (GWTC-4.0), not an independent analysis.

**Type-(B) status today:**
- no non-zero ΔΛ in any sector;
- GW: a bound consistent with 0;
- photons: bounds, to be confirmed;
- lab BEC: a null result with the stated low-k caveat.

## §7 H-bridge through c₄: unverified hypothesis (unchanged)

- **By hand,** in the single-band SOC EFT (δ = 0, g_↑↓ = g, mean field): Λ_disp = (ξ²/4)(1 + Q), with ξ = 1/(Mc), while the nonlinear observable is R = 1 + 12Q.
- **Status:** no symbolic check yet, no comparison with the two-band Bogoliubov spectrum, and no data. It does not identify (A) with (B). Today the count of independent physical support for Q is zero.

## §8 Priorities after registration (review of 2026-11-21)

1. **N4, photon bounds.** Check the primary sources (LHAASO, HAWC, Fermi), including the sign of each bound and the conversion to the canonical Λ. Literature only.
2. **N3, H-bridge.** A symbolic check of Λ_disp = (ξ²/4)(1 + Q), and a comparison with the full two-band Bogoliubov spectrum at δ = 0. Note → review → permission → one run. No physical conclusion before data.
3. **N2, Steinhauer.** It is inventoried (§4, N0 §4). Bringing the files into this repo is not urgent.
4. **New type-(B) line.** Look for a **linear** dimension-6 operator with an independent physical channel, instead of trying to get dispersion from the nonlinear type-(A) term (§2).

**Closed for now:** new QNM computations and a new GW scanner.

**On hold:** the GP test of R = 1 + 12Q. It is pre-registered (f18065a), with no script.
