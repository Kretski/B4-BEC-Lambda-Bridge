# EXT-1 — Two-component BEC model versus published measurements: specification (v0.2.1)

**Changes from v0.2 (review of 2026-11-16, approved with one clarification).** One sentence is added to the primary metric in §4: the statistical independence of γ_exp and γ_pred is stated as an assumption, unless a covariance is available. No other change.

**Changes from v0.1 (review of 2026-11-16, REVISION REQUIRED).**
- (A) The published ratio γ_exp = 0.193(22) of D1 is the primary observable as published. It is never recomputed from c_s/c_n.
- (B) The correlation of the uncertainties of a and a₁₂ is treated explicitly.
- (C) The density uncertainty is not taken from printed digits.
- (D) The density convention is defined.
- (E) D2 is restricted to Ω_R = 0.
- D3 becomes a separate candidate module.
- Journal references are updated.

**Status: PROPOSAL FOR REVIEW.** Specification only. No data have been extracted, no prediction has been computed and no comparison has been made. No code. Evidence level of every EXT-1 output: **level 3, analogue-system validation** (`UV_DYNAMIC_ETA_SCRIPT_SPEC_v0.3.2.md` §1a). No EXT-1 result may be described as level 5 (B4-specific empirical support). B4.12, the dynamic-η programme and T7a/T7b are not affected.

## 0. What EXT-1 can and cannot test (stated first)

**The model under test.** `theory/bec_dispersion.py`: the two-component mean-field Bogoliubov–de Gennes (BdG) dispersion with quantum pressure, without LHY corrections. Its characteristic polynomial is

    H(c, k²) = [(c − V₁)² − a₁ − b₁k²] [(c − V₂)² − a₂ − b₂k²] − γ,
    a_i = g_ii n_i / m_i,   b_i = ħ²/(4m_i²),   γ = g₁₂² n₁ n₂ / (m₁ m₂).

**What it can test.**
- (i) Whether this module reproduces measured density-mode and spin-mode sound velocities, and their ratio, from independently fixed microscopic parameters.
- (ii) Where a dataset resolves finite k, whether it reproduces the curvature of the branches. The curvature is the quantum-pressure term b_i, the source of β and of Λ_eff = 2β/c₀ in this repository.

**What it cannot test.**
- **For the symmetric mixtures available (m₁ = m₂, n₁ = n₂, g₁₁ = g₂₂, V₁ = V₂ = 0), the module reduces to the standard textbook mean-field result:**

      c_d,s² = (g ± g₁₂) n_i / m = (g ± g₁₂) n_tot / (2m),   so   c_s/c_d = √[(a − a₁₂)/(a + a₁₂)],

  with **n_i the density of each component** and n_tot = 2n_i the total density. This matches `bec_dispersion.py`, where a_i = g_ii n_i/m_i uses the component density n_i. The ratio does not depend on the density. Which density D1 uses in its formula for the absolute speeds (component or total, and the effective factor) is transcribed exactly at extraction, and the conversion to n_i is recorded.

  The first published dataset (§1) already states this prediction and its agreement. A successful comparison therefore confirms that the repository module **implements the standard mean-field description correctly**, and that this description holds on that platform. It does **not** discriminate the module from the standard BEC description. For these data the two coincide.
- **No dataset identified so far has the two components in relative motion (V₁ ≠ V₂).** The counterflow parameter dV is the one that matters for the B4 analogy (Stage 2, `docs/STAGE2_PREREGISTRATION.md`). It is therefore **not tested**.
- EXT-1 does not test B4, any BEC-to-B4 mapping, LHY corrections, or the nonlinear operators.

These limits are binding on every EXT-1 conclusion.

## 1. Datasets (selection fixed before any extraction)

| ID | Publication | Observables | Form of the data | Role |
|---|---|---|---|---|
| D1 | J. H. Kim, D. Hong, Y. Shin, "Observation of two sound modes in a binary superfluid gas", arXiv:1907.10289. Journal version as reported by the reviewer: Phys. Rev. A 101, 061601(R) (2020), DOI 10.1103/PhysRevA.101.061601. It is checked against the APS record at extraction, and the version used (arXiv or journal) is recorded. | density sound speed c_n, spin sound speed c_s, ratio γ_exp = c_s/c_n | values with uncertainties **stated in the text**; raw radius-versus-time data only in a figure | **primary** |
| D2 | R. Cominotti et al., "Observation of massless and massive collective excitations with Faraday patterns in a two-component superfluid", Phys. Rev. Lett. 128, 210401 (2022), arXiv:2112.09880 | density and spin dispersion ω(k), including the deviation from linearity of the spin branch, **only the uncoupled data with Ω_R = 0** | **colour plots only**; no tables and no numerical sound speeds | **secondary**, finite-k shape; digitisation required. The Ω_R ≠ 0 data (massive spin branch) are excluded: they need a model with coherent coupling |
| D3 | C. Piekarski, N. Cherroret, T. Aladjidi, Q. Glorieux, "Spin and density modes in a binary fluid of light", Phys. Rev. Lett. 134, 223403 (2025), DOI 10.1103/s58b-3mmx, arXiv:2412.08718 | two-branch dispersion by Bragg-like spectroscopy; a crossing of the branches at high photon density | a public data set is reported by the reviewer; to be verified in the separate module | **separate candidate module, not part of EXT-1 v0.2.** It is a fluid of light in a nonlinear atomic vapour, not an atomic BEC. At high intensity, saturation of the nonlinearity changes the ordering of the branches, and the standard two-component mean-field module does not describe that. It requires its own specification, including the processing of the intensity images and the treatment of saturation. Agreement in a low-intensity regime would not be accepted as confirmation of the saturated behaviour |

The selection is fixed here. A dataset is not added or dropped after any comparison has been seen. A dataset that turns out to be unusable is reported as such, with the reason.

## 2. Data extraction (before any prediction is compared)

**D1, observables. Two different experimental procedures are kept apart.**
- The absolute speeds c_n and c_s come from separate fits of the radius of each propagating pulse.
- The ratio γ_exp = 0.193(22) comes from a separate scale analysis: matching the times at which the two pulses reach equal radii.

γ_exp is therefore **not** computed as c_s/c_n, and its uncertainty is **not** propagated from those of c_n and c_s. The published γ_exp = 0.193(22) is the **primary observed quantity**. c_n and c_s are separate, secondary observations.

**D1, values stated in the text** (to be transcribed exactly, with the page and line recorded):
- c_n, c_s and γ_exp with their uncertainties;
- the scattering lengths a₁₁ = a₂₂ and a₁₂ with uncertainties;
- the gas parameter n₀a³;
- the effective density factor used by the authors;
- the trap frequencies, Thomas–Fermi radii, atom number, thermal fraction, the quadratic Zeeman energy and the authors' LHY estimate.

The values are entered in `ext1_data.json`, each with: source, location in the paper, the value as printed, its uncertainty as printed, and whether it is **measured** or **model-derived by the authors**. For example, the authors' γ₀ is model-derived and is recorded only as a reference, not as data.

**D2 (only if used).** The dispersion is available only as colour plots. Then:
- the digitisation tool and version;
- the image source (the arXiv PDF version), resolution and axis calibration points;
- the digitised points;
- the digitisation uncertainty in both axes, estimated by repeated independent digitisation (at least 3 passes);

are all recorded. Digitised points are labelled **DIGITISED**, never as primary data. Parameters not given numerically in D2 (n₀, a₁₂) are taken only as the paper states them (μ_d, μ_s, g₁₂ ≈ 0.93g). They are not fitted.

## 3. Reference hypothesis and model predictions (fixed before comparison)

**R0, reference hypothesis.** The standard mean-field two-fluid description, as used by the authors, with the authors' own stated inputs.

**P, prediction of the repository module.** `theory/bec_dispersion.py`, called through its public functions, with:
- m = m(²³Na);
- g = 4πħ²a/m and g₁₂ = 4πħ²a₁₂/m from the D1 scattering lengths;
- equal densities;
- V₁ = V₂ = 0.

**Observables predicted:**
- **PD1-a, the ratio** γ_pred = √[(a − a₁₂)/(a + a₁₂)] from the module. It needs no density.
  - **Uncertainty:** from the uncertainties of a and a₁₂, by linear propagation and by Monte Carlo (10⁴ samples, seed fixed in the script). a and a₁₂ enter both the numerator and the denominator, so their correlation matters.
  - **Correlation:** if D1 or its cited source publishes a covariance, it is used. If not, independence (ρ = 0) is the **stated assumption** for the primary result, not an established fact.
  - **Sensitivity:** the result is also reported for ρ ∈ {−0.9, −0.5, +0.5, +0.9}, with the z value for each.
- **PD1-b, the absolute speeds** c_n and c_s. These need a density n_i (component density, §0), from the authors' n₀ and effective factor.
  - The density uncertainty is the **published** uncertainty, if D1 gives one, or a budget reconstructed from the authors' stated measurements and calibrations, documented step by step.
  - The precision to which a number is printed is **not** used as an experimental uncertainty.
  - If neither a published uncertainty nor a reconstruction is possible, PD1-b is reported with **incomplete uncertainty** and is **excluded from the quantitative conclusion**. It is shown for information only.
- **PD2, only if D2 is used, and only for the Ω_R = 0 data:** the curvature of the spin branch, from the full H(c, k²) with the stated μ_s and g₁₂/g, over the digitised k range.

The predictions are computed and stored in `ext1_predictions.json`, with SHA-256, **before** the comparison script is run. The comparison script reads both files and does not recompute the predictions.

**Honest limitation on blinding.** D1 itself publishes the reference prediction (γ₀) and its agreement. A blind prediction is therefore not possible for D1. EXT-1 records this. Its value for D1 is the independent reproduction by the repository code, not a blind test.

## 4. Comparison and criteria (fixed before comparison)

- **Primary metric (D1, ratio):** z = (γ_exp − γ_pred)/√(σ_exp² + σ_pred²).
  - **PASS** if |z| ≤ 2;
  - **TENSION** if 2 < |z| ≤ 3;
  - **FAIL** if |z| > 3.
  - The primary z-score assumes that the experimental scale-analysis result γ_exp and the prediction derived from the independently measured scattering lengths are statistically independent. If a covariance is available, it is included explicitly; otherwise independence is stated as an assumption.
- **Secondary metric (D1, absolute speeds):** the same z for c_n and c_s separately, with the density uncertainty included. These are reported, not decisive, because of the density calibration.
- **Both branches:** c_n and c_s are always reported together. No branch is dropped.
- **D2 (if used), Ω_R = 0 only:** residuals of the digitised spin-branch points against PD2, normalised by the digitisation uncertainty. χ²/dof is reported with the number of points. No PASS/FAIL threshold is set for D2 in this version. The result is descriptive, because the data are DIGITISED.
- **Fitted versus fixed parameters.** In D1 no parameter is fitted. If any parameter is ever fitted, it is listed as fitted, and the result is reported separately from the fixed-parameter comparison.

**Sensitivity tests** (fixed list):
- (a) a₁₂ varied over ±1σ and ±2σ;
- (b) the LHY estimate quoted by the authors applied as a systematic shift to c_n;
- (c) for PD1-b, f_eff varied between the bulk value 1 and the stated value.

**Exclusions.** None planned. Any exclusion is documented with its reason before the comparison is rerun, and both results are reported.

## 5. Code-level checks (internal)

- **E1.** For the symmetric case, the module's c_d and c_s equal the closed form √[(g ± g₁₂)n/m] to machine precision.
- **E2.** The ratio γ_pred from the module equals √[(a − a₁₂)/(a + a₁₂)] to machine precision.
- **E3.** Unit test of the conversion a → g in SI units, with CODATA constants stated.

E1–E3 check the implementation only. They are not part of the empirical verdict.

## 6. Outputs and order of work

**Order:** this specification, then review, then registration, then `ext1_data.json` (extraction, reviewed), then `ext1_predictions.json` (computed and registered by hash), then the comparison script, then static review, then one authorised run, then review of the results, then registration.

**Outputs:** `results_ext1/` with the comparison JSON, the console output, `RUN_CONTEXT.txt` and `SHA256SUMS.txt`, as for the other runs. Every number is labelled measured, model-derived (authors), DIGITISED or predicted (repository).

## 7. Out of scope

- B4, any BEC-to-B4 mapping, Stage 2 (dV ≠ 0), LHY beyond the authors' quoted estimate, nonlinear operators.
- Any claim of evidence level 4 or 5.
- Fitting the module to the data to obtain agreement.
- The dynamic-η programme, B4.12, T7a, T7b, Zenodo.

## 8. Requested from the reviewer

1. Approval of §0 as the binding statement of what EXT-1 can and cannot test. In particular: for the available symmetric data, the module coincides with the standard mean-field description, so EXT-1 validates the implementation and that description, not a discriminating prediction.
2. Approval of D1 as primary (ratio, no density calibration), and of D2 as secondary and DIGITISED, descriptive only.
3. Confirmation that D3 is a separate candidate module with its own specification, outside EXT-1 v0.2 (reviewer's decision of 2026-11-16).
4. Approval of the thresholds |z| ≤ 2 (PASS), 2 < |z| ≤ 3 (TENSION) and |z| > 3 (FAIL) for the D1 ratio, fixed before any comparison.
5. Agreement that the absence of counterflow data (V₁ ≠ V₂) is recorded as an open gap. A separate search for such data could be a later task.

**B4.12: OPEN. Dynamic-η script A: specification approved, no code yet. No run authorised. Zenodo: paused.**
