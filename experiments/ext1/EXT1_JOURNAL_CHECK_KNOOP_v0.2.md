# EXT-1 — journal-version check: Knoop et al., PRA 83, 042704 (2011) — v0.2

**Status (unambiguous):**

- **Knoop sub-check: CLOSED.**
- **EXT-1 journal-version check: OPEN, pending Kim et al.**
- **EXT-1 comparison: NOT AUTHORISED.**

This is a separate provenance note. It is not a validation of EXT-1. v0.2 applies the precision points of the review of 2026-11-17.

The registered files `experiments/ext1/ext1_data_v3.json` and `ext1_predictions.json` (fa19637) are **not** changed. This note records only where the journal version agrees with or differs from the arXiv version that was used.

**Comparison:** still NOT authorised. The journal version of Kim et al. (PRA 101, 061601(R)) is still missing.

## Sources

| Version | File | SHA-256 | Pages |
|---|---|---|---|
| Journal | `PhysRevA.83.042704.pdf` (supplied by the user 2026-11-17) | `272c1a1ed6784282dae99fa6da31348172b4423b6e9bebd5b76c7962deac4268` | 10 |
| arXiv | `knoop2011.pdf`, arXiv:1102.0572v1 (used for the v3 extraction) | `d7fae796…` | 11 |

How the journal version was identified:

- page header "PHYSICAL REVIEW A 83, 042704 (2011)";
- "©2011 American Physical Society" and "1050-2947/2011/83(4)/042704(10)" on p. 042704-1;
- "Received 2 February 2011; published 8 April 2011".

**Method.** The text of all 10 journal pages was extracted. The pages containing the values were found by a search for 54.54, 50.78, 52.66 and 48.91; these are pp. 9–10 (042704-9 and 042704-10). Those two pages were read in full and compared sentence by sentence with the arXiv text (arXiv pp. 9–11).

## Values used by EXT-1

| Entry in v3 | arXiv | Journal | Result |
|---|---|---|---|
| `knoop_a`: \|1,1⟩+\|1,1⟩, B = 0 | 54.54(20)a₀ | 54.54(20)a₀, p. 042704-9 | **identical** |
| equality with \|1,−1⟩+\|1,−1⟩ and \|1,±1⟩+\|1,0⟩ | stated | stated, same sentence, p. 042704-9 | **identical** |
| `knoop_a_updn`: \|1,1⟩+\|1,−1⟩ | 50.78(40)a₀ | 50.78(40)a₀, p. 042704-9 | **identical** |
| \|1,0⟩+\|1,0⟩ (not used) | 52.66(40)a₀ | 52.66(40)a₀ | identical |
| `knoop_coupled_basis`: (11)2 and (11)0 | 54.54(20)a₀ and 48.91(40)a₀, footnote [46] | the same values, footnote **[47]**, p. 042704-10 | values identical; footnote number changed |
| `knoop_correlation_statement` | "strong correlation between this scattering lengths and the one for \|1,1⟩+\|1,1⟩" (triplet) | same sentence, p. 042704-9 | **identical** |

What the check establishes, and what it does not:

- the **values and uncertainties** used by EXT-1 are identical in both versions;
- the **channel assignments** are identical;
- the **footnote number and the page locations** differ (see below), so the two versions are **not** identical as sources;
- **neither version publishes a covariance matrix or a numerical correlation coefficient** for a and a↑↓.

Flag 2 of v3 (possible positive correlation; ρ = 0 primary; four-value ρ analysis mandatory) therefore **cannot be removed** on the basis of this check, and stays unchanged.

## Differences that affect provenance only

1. **Reference numbers.** Journal [46] = van Abeelen & Verhaar, and [47] = the coupled-basis footnote. The arXiv version has them the other way round. A citation of "Knoop [46]" for the coupled-basis values must read "[47]" for the journal version.
2. **Location of the a₀ value.** The arXiv text gives "(a₀ = 0.529177 × 10⁻¹⁰ m)" next to 54.54(20)a₀. In the journal it appears in the potential-parameter list (R₀ = 11.0a₀, a₀ = 0.529 177 × 10⁻¹⁰ m). EXT-1 uses CODATA 2018 a₀ in any case, so nothing changes.
3. **Page and line numbers.** v3 cites arXiv page/line (p. 9 l. 125; p. 10 l. 100). The journal equivalents are pp. 042704-9 (main text) and 042704-10 (footnote [47]). Journal line numbers depend on the extraction tool and are not given.

## Consequence

- **Knoop sub-check: CLOSED**, with no change to any value. This closes only the Knoop part. It does not close the EXT-1 journal-version check.
- The provenance differences (reference number, pages) are recorded here rather than by editing v3. Editing v3 would change its SHA-256, which `ext1_predict.py` fixes.
- **Open, and the real blocker:** the journal version of Kim et al., PRA 101, 061601(R) (2020). It must be checked directly, from the published PDF:
  1. γ_exp = 0.193(22);
  2. the definitions and values of c_n and c_s;
  3. the density convention (n, n_eff = (2/3)n₀, n₀a³);
  4. for each quantity, whether it is measured, obtained from a fit, or calculated from other quantities. The same symbol in two versions does not guarantee the same physical definition.

  The abstract of the published version (as listed on the group's website) is worded differently from arXiv v1 on how the velocity ratio was determined. This is a further reason not to assume the two versions agree. The second file supplied (`QGL - RESEARCH.pdf`, SHA-256 `9d119161…`, one page) is a printout of the research page of the Shin group (Quantum Gases Lab). It lists the paper with its abstract. It is not the paper and cannot replace the check of γ_exp, c_n, c_s or the density convention.
