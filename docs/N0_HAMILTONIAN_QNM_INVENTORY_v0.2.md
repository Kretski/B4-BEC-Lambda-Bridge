# N0 — inventory of the Hamiltonian / eikonal / QNM method (Lambda-model material) — v0.2

**Status:** inventory only. Files were read; nothing was executed and no computation was made. It is registered together with `LAMBDA_INDEPENDENT_EXTRACTION_MAP_v0.4.md`, in the commit that adds both, on top of f18065a. The registration fixes what the files say; it does not raise any finding to independent physical confirmation.

**Supersedes:** v0.1 (1bd8b1e4…). Every **[needs file]** item of v0.1 is now resolved from the collected files.

**Sources:**
- `lambda_model_inventory_report.txt` (85a604b8…): the search report, matching lines only;
- `n0_files.zip`: 19 files, collected read-only on 2026-11-21 at 16:20 by `collect_n0_files.ps1` (984cdd8b…). The SHA-256 of every file was checked against the zip manifest: 19/19 match (Appendix A).

Statements in Lambda-model documents are reported as that repository states them. Numbers from its outputs are quoted, **not re-run**.

## 1. The five N0 questions

**Q1. Where is the mass-shell Hamiltonian, and which normalisation does it use? Answer: canonical, with no factor 2.**
- Paper 1, "A Hamiltonian Oscillator Extension of Wave Propagation in Schwarzschild Spacetime" (`vacuum_oscillator_v0.4.tex`, l.121–123): ω² = c²k²α²(1 + Λk²).
- The later preprint, "An Oscillator-Motivated Hamiltonian Model of Dispersive Wave Propagation in Schwarzschild Spacetime" (`vacuum_oscillator_preprint_v4.tex`), has the same form:
  - l.224: ω² = c²α²k_loc²(1 + Λk_loc²), with Λ > 0;
  - l.249 and l.724: mass shell g^{μν}k_μk_ν + Λk_loc⁴ = 0, with k_loc² = h^{μν}k_μk_ν.
- Paper 2 code (`A1v3_zamo_v04.py`, l.188–213): **H = ½[g^{μν}k_μk_ν + Λ(k_loc²)²]**. The factor ½ multiplies both terms.
- With g^{μν}k_μk_ν = k_loc² − ω²/(α²c²), the condition H = 0 gives ω² = α²c²k_loc²(1 + Λk_loc²). This is the canonical form.

**Correction to map v0.3, §1.** The row "mass-shell H = ½g kk + Λ_H k⁴ ⇒ Λ = 2Λ_H" describes a form that does **not** occur in these files. It is withdrawn. The Hamiltonian route is canonical throughout.

**Q2. How is H_V4 obtained, and where does "1 + 2Λ" appear? Answer: H_V4 is canonical. "1 + 2Λ" appears there only in derived quantities.**
- H_V4 = αk_loc√(1 + Λk_loc²), with v_g = α(1 + 2Λk²)/√(1 + Λk²).
- In the einbein formulation the variable change is k_loc = w(1 + 2Λw²).
- The later preprint (§"Three Hamiltonians from the quartic structure", l.324–440) gives V4, V6 and the einbein form, and their non-equivalence (eq:nonequiv). `proofs_v4v6_nonequiv.py` is its check.
- None of these is a second dispersion normalisation.

**Q3. Which QNM and photon-ring results exist?** Unchanged from v0.1:
- the Kerr photon-ring response;
- the shadow correction;
- Paper 2 Table 3 / Eq. 20;
- the retrograde-ring D_Λ^eik = 0.0127 Λ̂ l³ at a = 0.9, "negligible at natural Λ";
- eikonal QNM estimates.

**Q4. Regression tests or data? Answer: regression and implementation checks only.**
- `THEORY_STATUS_GW.md` §9:
  - **A2(a):** on Schwarzschild, with the Boyer–Lindquist / ZAMO foliation, "the quasinormal-mode problem is not well posed in this setup";
  - Paper 2 QNMs "remain eikonal estimates";
  - decision (October 2026): "do not start new quasinormal-mode computations".

**Q5. Is there an independent path to Λ_GW? Answer: no.**
- `THEORY_STATUS_GW.md` (header, §3, §6) says:
  - Papers 1–2 derive the quartic dispersion **for a scalar field only**;
  - a GW equation requires tensor dynamics with a coefficient Λ_g that is **independent** of the scalar Λ unless some principle fixes it;
  - the identification Λ_GW = Λ "is an additional universality assumption".
- §8.1: the flat-space form equals the Hořava–Lifshitz tensor dispersion. A GW phase test of this form "is therefore not specific to the Lambda model".

**Conclusion of N0 (unchanged, now on full files).** The Hamiltonian / eikonal / QNM route is a developed scalar-field and ray framework. It has **no observational channel** for an independent extraction of the dispersive Λ at present. Type-(B) information stays with the GW, photon and lab-BEC bounds (§4–§6).

## 2. Finding A0: the original Λ was a type-(A) term

`THEORY_STATUS_GW.md` §2 and §9:
- Paper 1 derives the mass shell from the action term **(Λ/2)(h^{μν}∂_μφ ∂_νφ)²**.
- Item A0 states that this term "is quartic in the field, drops out at linear order, and has mass dimension 8 rather than 6". **It does not produce the stated linear k⁴ dispersion.**
- A **linear**, dimension-6 operator built from the Laplacian of the preferred slices does give ω² = k²(1 + Λk²). On Kerr it reproduces the Paper 2 Hamiltonian in the eikonal limit (A1). It is listed as a **candidate**.

**Meaning for the Λ line (map v0.3, §0–§1).**
- In the terms of the map, the original Lambda-model action term is a **type-(A)** operator: nonlinear, quartic in fields, the same class as Λ_B4 (Λ/4)(∂φ)⁴.
- The dispersion it was meant to produce is **type (B)**.
- The Lambda-model's own audit found that (A) does not generate (B) at linear order. This is the same A/B separation as the 2026-11-21 decision, found independently in the other repository.
- It supports keeping (A) and (B) apart. It also shows that a type-(B) Λ needs a **linear** higher-derivative operator, not a self-interaction.

**Lab analogue.** In a BEC the k⁴ term of the Bogoliubov spectrum, ħ²k⁴/(4m²), comes from the **kinetic** term, the quantum pressure. The self-interaction enters only through c. Paper 1 asks (§9) "What is the gravitational analogue of the GP self-interaction that sets Λ?" Within the GP equation the self-interaction does **not** set the k⁴ coefficient; it sets the scale c. This is a remark for the open question, not a result.

## 3. Conventions (N1): resolved for the files collected

- **The canonical form** ω² = c²k²(1 + Λk²) is used by Paper 1 (v0.4, l.123; preprint v4, l.224), Paper 2 / A1v3_zamo_v04, THEORY_STATUS_GW (§2, §5), SENSITIVITY_MAP, `example_BEC.py` and `paper3_final.tex`.
- **The obsolete form** is "Kretski: ω² = c²k² + 2Λc²k⁴", i.e. 1 + 2Λk², with **Λ_old = Λ/2**.
  - It is used in `lambda_wave_equation.py` (outside the git repo, dated 2026-08-23, l.175–186) and in the files listed in v0.1 §3.
  - It is labelled obsolete in `PAPER3_VALIDATION_STATUS.md` (l.62, l.187).

**BEC, "Λ = ξ²/4": the physics is consistent and only the symbol differs.**

| source | normalisation | ξ | stated | canonical Λ |
|---|---|---|---|---|
| `lambda_wave_equation.py` l.162–186 | old (1 + 2Λk²) | ħ/(√2 m c_s) | Λ_old = ξ²/4 | ξ²/2 = **ħ²/(4m²c_s²)** |
| `example_BEC.py` l.13–21; `paper3_final.tex` l.343–360 | canonical | ħ/(m c_s) | Λ = ξ²/4 | **ħ²/(4m²c_s²)** |

Both give the same physical value, **Λ_BEC = ħ²/(4m²c_s²)**.
- "ξ²/4" is correct only together with its ξ definition.
- The decision text of 2026-11-21 ("Λ_BEC = ξ²/4") is correct with ξ = ħ/(m c_s), the Paper 3 convention.

`gwosc_zero_crossing_injection_recovery.py` (l.108, 205: f_obs = f√(1 + 2Λf²)) is a legacy injection script in the obsolete normalisation, applied to frequency. It is not used by the matched-filter results. Converting it is part of the Lambda-model to-do list (`PAPER3_VALIDATION_STATUS.md` l.187).

## 4. Lab BEC null test (N2): numbers as recorded

Source: `steinhauer_test_v3.py` and its output `steinhauer_kmin_scan.csv`; data from Steinhauer et al., PRL 88, 120407 (2002), 13 points extracted from the vector EPS of Fig. 3.
- **H0** is the local-density Bogoliubov spectrum. **H1** is H0 with the k⁴ coefficient scaled by (1 + g_extra). So **g_extra = ΔΛ_lab / Λ_Bog**, a relative excess over the standard value.
- A Gaussian prior is used on μ/h = 1.91 ± 0.09 kHz (published), with a profile likelihood and, in the k_min scan, a null test of 1000 pseudo-experiments.

| k_min [µm⁻¹] | N | g_extra ± err | 95% interval | Δχ² | p_emp | H0 goodness of fit p |
|---|---|---|---|---|---|---|
| 0.0 | 13 | 0.031 ± 0.023 | [−0.01, 0.07] | 1.79 | 0.19 | 5.7e-5 |
| 0.8 | 11 | 0.026 ± 0.023 | [−0.01, 0.07] | 1.23 | 0.29 | 2.8e-3 |
| 1.5 | 9 | 0.012 ± 0.024 | [−0.03, 0.05] | 0.28 | 0.59 | 0.28 |
| 2.0 | 8 | 0.019 ± 0.024 | [−0.02, 0.06] | 0.62 | 0.44 | 0.37 |

**Reading:**
- There is **no significant excess**: ΔΛ_lab/Λ_Bog is consistent with 0 at about the 5–7% level (95%).
- The interval edges are quantised by a profile grid of step 0.01.
- **Caveat:** H0 fits poorly when the low-k points are included (p ≈ 6e-5 at k_min = 0). The model is adequate from k_min ≈ 1.5 µm⁻¹. So the bound is model-dependent at low k (local density, the error flags marked "sigma_is_upper_bound").
- The script states that g_extra "is NOT an estimate of the fundamental Λ … no theoretical bridge between the two exists".

## 5. Analogue giant vortex and EXP C

- `CORRECTION_QNMS_CSV.md` (dated 19 Oct 2026, confirmed by P. Švančara):
  - every "experimental" EXP B comparison is **model vs model**, i.e. two implementations of the same WKB model;
  - the agreement of 0.01–0.06% (|m| ≥ 10) measures implementation agreement;
  - "a comparison [of EXP C] with measured spectra has not been made".
- **EXP C digitisation (uncommitted).** `compare_expC_digitization.py` compares the frozen predictions with markers digitised from published Fig. S5. It states that this is "NOT a comparison with the original numeric experimental dataset".
  - Rank matching without tolerance gives 38/58 matched and MAE = 0.022 Hz.
  - **Whether the Fig. S5 markers are measured peaks or model lines is not stated in these files.**
  - Either way, it tests the vortex's own capillary–gravity model, not Λ.
- The stale "measured / experimental" wording listed in v0.1 §2 predates this correction. That remains a correction task for the Lambda-model repo.

## 6. Photon bounds (N4): source of the numbers

- **The registered map of this repo (Λ_γ ≲ 1e-53 m²)** traces to Paper 1: "Λ ≲ 2×10⁻⁵³ m² from Fermi-LAT (GRB 090902B)" (`THEORY_STATUS_GW.md` §4, order of magnitude).
- **The newer superluminal bounds** in `notes/SENSITIVITY_MAP.md`, using the conversion Λ = (ħc/E_QG,2)² (orders of magnitude only, "to be checked against each source before citation"):
  - photon time of flight, LHAASO GRB 221009A: 1e-55 to 7e-57 m²;
  - photon decay: 6e-60 m²;
  - photon splitting: 4e-62 (HAWC) and 1e-64 m² (LHAASO);
  - neutrinos: about 3e-62 m²;
  - GW: 7e-12 m².
- SENSITIVITY_MAP's own reading: if Λ were universal, photon splitting would already exclude the GUT-scale natural value. A gravity-only Λ is tested only by GW, about 58 orders above the Planck-scale value.
- **For this repo (N4):** the registered map's 1e-53 m² is outdated. Before any citation, each newer bound must be checked against its source, as SENSITIVITY_MAP itself requires.

## 7. Changes to map v0.3

| map item | change |
|---|---|
| §1, Hamiltonian row (Λ = 2Λ_H) | **withdrawn**: the Hamiltonian is canonical |
| §1, source of "1 + 2Λk²" | the obsolete "Kretski" normalisation, Λ_old = Λ/2 |
| §1, BEC row | the physics is consistent: Λ_BEC = ħ²/(4m²c_s²) = ξ²/4 with ξ = ħ/(m c_s) |
| §3, BEC row | ΔΛ_lab/Λ_Bog = 0.03 ± 0.02 (all points), about 95% ≤ 0.07; consistent with 0; H0 fit adequate only for k_min ≥ 1.5 µm⁻¹ |
| §3, photon row | 1e-53 m² is the Paper 1 Fermi value; newer LHAASO-based bounds are 1e-55 to 1e-64 m² (to be source-checked) |
| §3a, Hamiltonian / QNM route | no observational channel (A2(a): QNM problem not well posed; eikonal effects negligible) |
| new: A0 | the original Lambda-model Λ was a type-(A) term that does not give type-(B) dispersion; the type-(B) candidate is a linear dimension-6 operator |

## Appendix A — files read (SHA-256 verified against the zip manifest, 19/19)

| # | file | sha256 (first 16) |
|---|---|---|
| 01, 13 | A1v3_zamo_v04.py (root copy and repo copy, identical) | b9da944ffe143ff0 |
| 02 | lambda_wave_equation.py (outside the repo) | 6bc43a033c2ffa60 |
| 03 | proofs_v4v6_nonequiv.py (outside the repo) | 3468f8b11d02cb1b |
| 04 | repo/examples/example_BEC.py | 534a6108a08c309b |
| 05 | repo/notes/EIKONAL_D_LAMBDA.md | 05863d538f79b01d |
| 06 | repo/notes/SENSITIVITY_MAP.md | 4d7e0fad8a85daa0 |
| 07 | repo/paper3/bec_validation/steinhauer2002/steinhauer_kmin_scan.csv | f9a007638d068ce4 |
| 08 | repo/paper3/bec_validation/steinhauer2002/steinhauer_test_v3.py | 95461853dfc41486 |
| 09 | repo/paper3/bec_validation/steinhauer2002/steinhauer2002_fig3_extracted.csv | 8162274b82a16a0a |
| 10 | repo/paper3/CORRECTION_QNMS_CSV.md | b5a02dfa1486b434 |
| 11 | repo/paper3/expC_prediction/compare_expC_digitization.py (uncommitted) | c1a6d70ca7201c90 |
| 12 | repo/paper3/expC_prediction/comparison_expC_digitization_summary.csv (uncommitted) | 853f282ddb6af188 |
| 14 | repo/paper3/gwosc_zero_crossing_injection_recovery.py | 4472c4525b18f675 |
| 15 | repo/paper3/paper3_final.tex | 2e38e3a9ed3e42ce |
| 16 | repo/paper3/PAPER3_VALIDATION_STATUS.md | 354e80ca7ebddec6 |
| 17 | repo/THEORY_STATUS_GW.md | 96987ae7ad12b5d6 |
| 18 | vacuum_oscillator_preprint_v4.tex (later Paper 1 preprint, outside the repo) | b135051c9ee57db1 |
| 19 | vacuum_oscillator_v0.4.tex (Paper 1, Zenodo version v0.4, outside the repo) | 543bf7463efcc43f |
