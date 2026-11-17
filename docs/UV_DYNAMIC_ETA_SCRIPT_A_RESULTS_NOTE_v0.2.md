# Dynamic η — Script A rev 1.1: results note v0.2

**Status:** v0.2, with the editorial points of the review of 2026-11-16 applied (v0.1 approved in substance). This note interprets the registered run. It does not change any output of the run.

**Evidence level:** the physical claim remains at evidence level 1; DC0–DC6 provide internal numerical consistency at level 2. This is not experimental validation, and not evidence for B4. C4 (Script B) was not evaluated, so "passes T3-fund" cannot be stated for any case.

## 1. The run

| Item | Value |
|---|---|
| Code | `experiments/uv_candidates/uv_dynamic_eta_check.py` rev 1.1, commit `20efdce`, SHA-256 `fa46712b…` |
| Specification | `docs/UV_DYNAMIC_ETA_SCRIPT_SPEC_v0.3.2.md` (aece78a), unchanged |
| Inputs | B3 `99a6b04b…`, manifest `1bae892a…`, verified before the start and by the script |
| Execution | one authorised run, 2026-11-16 18:31:42–18:45:45 UTC, exit code 0, no retry |
| Environment | Anthropic cloud workspace (Linux x86_64), clean GitHub checkout with LF line endings. **Not** the local laptop |
| Outputs | `results_uv_dynamic_eta_rev1_1/` (JSON `de758935…`, log `7c931d1a…`, `RUN_CONTEXT`, `SHA256SUMS`) |
| Completeness | 41 backgrounds; 31 have a real η̄, 10 do not (the B3 anchor is among the 10); 39 simple, separated roots; 117 (root, f) combinations |
| Internal checks | DC0–DC6 PASS, so ALL_DC_PASS = true. This covers the internal checks only; it is not a physics verdict |

Root numbering: this note uses the JSON field `root_index`, which starts at 0. The review of 2026-11-16 counts roots from 1, so its "root 1, 1, 2, 3" are `root_index` 0, 0, 1, 2 here.

## 2. T3-fund (C0–C3)

The counts below are identical for f = 0.1, 1 and 10.

| Outcome | Count | Detail |
|---|---|---|
| no real η̄ (INCONCLUSIVE, C0) | 10 | |
| FAIL-UNSTABLE | 35 | 21 roots: P2a FAIL-UNSTABLE at k = 0 (all TACHYONIC). 14 roots: P2a ZERO MODE, but C3 FAIL-UNSTABLE at k > 0 (6 HEAVY, 8 TACHYONIC) |
| PASS (C0–C3; C4 not evaluated) | 4 | listed in §3 |

So the failure is not explained by the tachyonic class of η alone. Six roots with a heavy, positive η mass still have a confirmed growing mode at k > 0.

C1 and C2 pass analytically for all 39 roots. The principal symbol is block-diagonal and strongly hyperbolic. The failures are lower-order instabilities (C3), not ill-posedness.

**Reading.** For the registered Lagrangian, the 41 backgrounds and the fixed criteria, dynamic η does not show general stability. Four roots out of 39 are consistent by C0–C3. This is a negative result for this specific realisation. It does not show that every dynamical realisation of a preferred field is impossible.

## 3. The four cases that pass C0–C3

Each of these four specific roots is a minimum of V_eff, has P2a ZERO MODE, and the full system confirms its η class at k = 0 (`conf` true). This is a statement about these four solutions, not a general property of the model: 29 of the 39 roots are maxima of V_eff (TACHYONIC). T6 is STATED-DYNAMICAL for all four.

| Background | `root_index` | η̄ | m²_η,eff at f = 0.1 / 1 / 10 | Class at f = 10 | T5-static, u frame | T5-static, lab frame | η-exchange rel. norm |
|---|---|---|---|---|---|---|---|
| manifest_20_index_6052 | 0 | −0.0855 | 7.58 / 0.758 / 0.0758 | HEAVY | PASS | PASS | 1.06 |
| manifest_21_index_522 | 0 | 0.7841 | 7.77 / 0.777 / 0.0777 | HEAVY | PASS | **FAIL** | 0.83 |
| manifest_28_index_523 | 1 | 0.1859 | 2.70 / 0.270 / 0.0270 | HEAVY | PASS | PASS | 0.144 |
| manifest_30_index_825 | 2 | 1.6212 | 0.718 / 0.0718 / 0.00718 | **LIGHT** | PASS | **FAIL** | 1.20 |

In the u frame, for all four roots and all three values of f:

- the velocity block A and the large-k block C₂ are positive;
- the Schur complement at k → 0 is positive;
- the grid minimum is positive;
- H₂ has no negative eigenvalue at k = 0.01, 0.1, 1 and 10.

### The two lab-frame T5 failures

In the lab frame the velocity block A has a negative eigenvalue: −0.512 for manifest_21 and −0.630 for manifest_30. The value does not depend on f, so it comes from the matter sector, not the η kinetic term. These are the two roots with the largest rapidity, η̄ = 0.78 and 1.62.

Because A is not positive definite there, H₂ and the Krein test are NOT EVALUATED in the lab frame for these two roots. This matches the code rule.

This is consistent with the registered interpretation, which uses the preferred foliation as the basis for T3-fund and T5. It shows that positivity of the energy depends on the frame. It must not be summarised as "positive energy in every frame".

## 4. T6

- **Set level: PARTIAL.** η̄ exists on 31 of the 41 backgrounds. The 10 without η̄ include the B3 anchor.
- **Per root:** 4 STATED-DYNAMICAL and 35 FAIL, the same for every f.

## 5. T7

The result is NOT-REPRODUCED for all 39 roots and all three values of f. The reasons differ by class:

- **TACHYONIC (29 roots).** Under rule F2(ii′), eliminating η at a maximum of V_eff is not admissible in this regime. This is a statement about the elimination procedure. It is not a proof that no UV origin exists.
- **HEAVY (10 roots; 9 at f = 10).** Tree-level elimination is admissible, but the η-exchange correction to the quadratic operator does not vanish. Its relative norm ranges from 0.144 to 36.8.
- **LIGHT (1 root, manifest_30 `root_index` 2, only at f = 10).** Not eliminable (F2(ii)).

The frozen-η comparison is exact:

- operator-by-operator identity with rank-one B4, s = κ²/p;
- identical characteristic polynomial;
- the Λ term gives 1 + 3Λq².

The rank-one B4 structure is therefore reproduced with η frozen. It is **not** preserved once η is treated dynamically and eliminated.

## 6. Diagnostics

- **P1.** det K_u is irreducible, of degree 8, over the field of rational functions of the background parameters. Particular numerical backgrounds where it factorises are not excluded.
- **DC1.** 9360 matrix entries agree. There are no disagreements and no R2 non-convergence.
- **X1.** The maximum relative frame-consistency deviation is 1.9 × 10⁻⁴⁵.
- **DC5.** The historical count of 17/41 is reproduced. This checks the earlier classification; it is not a new physical confirmation.
- **X2.** For 6 roots (at each value of f) the stored X2 is between 0.73 and 4.27. The registered X2 metric pairs the two root lists after sorting by real part only. Complex-conjugate pairs share a real part, so they can appear in opposite order in the two lists.

  A post-hoc, order-independent comparison of the stored roots gives a difference of 0 in all 117 cases. The method is a two-sided nearest-neighbour (Hausdorff-type) distance between the two root lists. For lists of equal length with distinct roots, as here, this is adequate for the present diagnostic. It is not a general algorithm for proving equality of multisets with possibly repeated roots. This is in `x2_posthoc_order_independent.txt`, produced by `x2_posthoc_order_independent.py`.

  That check is a separate diagnostic artifact and was not part of the run. The original JSON and log are unchanged. X2 remains outside ALL_DC_PASS, as registered.

  Recorded defect of the original metric: it is sensitive to the order of complex-conjugate roots.

## 7. What this result does and does not change

- It is the first complete Script A result on the full pre-registered set.
- It does **not** support dynamic η as a generally stable UV completion of B4.
- It gives four bounded cases of consistency by C0–C3 and T5 in the u frame. These are not evidence for the whole model.
- The 35 T3-fund failures stand. No threshold or model parameter is to be retuned.
- The status of T7a, Script B/C4, T7b and B4.12 is unchanged. No next stage is authorised by this result.
