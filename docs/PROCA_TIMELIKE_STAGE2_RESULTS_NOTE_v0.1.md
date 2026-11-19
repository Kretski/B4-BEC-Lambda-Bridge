# Proca with a timelike current — Stage 2 results note v0.1

**Specification:** `docs/PROCA_TIMELIKE_CURRENT_ANALYTIC_SPEC_v0.5.md` (87fc31b), unchanged.

**Status:** Stage 2 (arithmetic only) accepted by the reviewer on 2026-11-19. Model class: **OPEN**.

**Evidence level:** 1 (arithmetic of derived formulas on the registered backgrounds; no new criterion).

**Not run:** roots of Q_new(c; 1), V3, NOT-VIABLE (d), Stage 3. Each needs a separate decision.

## 1. Run

| Item | Value |
|---|---|
| Script | `proca_timelike_stage2.py`, SHA-256 `585143fa50d390b27e40b6aa1d8a2114a32e2b881600d6ee21d0cc22e27a5ff6` |
| Inputs | `data/B3_max_speed_7_315096.json` `99a6b04b…`; `data/e3_manifest.json` `1bae892a…` (41 backgrounds) |
| Repository state | fresh clone, master = origin/master = 87fc31b, clean, LF |
| Command | `python3 -I proca_timelike_stage2.py <outdir>` from the repository root; one run, rc = 0, 3.6 s |
| Environment | Python 3.13.18, SymPy 1.14.0, mpmath 1.3.0 |
| Result | `proca_timelike_stage2.json` `789e5cc4…`; `SHA256SUMS` `aa65d115…` |

Conventions: M = 1; g_φ = +√p, g_ψ = −κ/√p; γ_u = σ√(Λ/2), σ = ±1 evaluated separately; kernel with s_eff = κ²/p.

**Script revisions before the run (all statically reviewed):**

- `16fa3506` → `ff359c11`: exact R-UNDEFINED identification (rational isolating intervals, Sturm counts, gcd), V2 and H3 = None when there is no real branch.
- `ff359c11` → `bd6eb9fc`: every D3 sign decision is certified by outward-rounded interval arithmetic on exact rational root enclosures, with adaptive precision; an unseparated sign gives None. Rank-one input check: exact δ_rank1 = |κ²/(ps) − 1|, ABORT above 10⁻⁹ (the registered K4 threshold).
- `bd6eb9fc` → `585143fa`: s < 0, or s = 0 with κ ≠ 0, is rejected.

Synthetic smoke test: 18/18.

## 2. Controls

- c_max against the registered common-M CSV: PASS (41/41).
- Spectral-norm ratio from the exact kernels against the N formula at x = 0.5: maximum relative difference 3.4·10⁻⁴⁹.
- max δ_rank1 = 2.96·10⁻¹⁶ (exact on the rationalised inputs; the registered float K4 gave 2.2·10⁻¹⁶).

## 3. D4 (common-M, sets and BORDERLINE rule unchanged)

| Set | T7a PASS / FAIL | Extended (c′_max) PASS / FAIL | Smallest \|margin\| |
|---|---|---|---|
| A | 0 / 41 | 0 / 41 | 0.257 |
| B | 0 / 41 | 0 / 41 | 0.457 |
| C | 21 / 20 (registered result reproduced) | 19 / 22 | 0.011 |

- No BORDERLINE.
- c_Λ dominates on 9 backgrounds.
- Two transitions PASS → FAIL in set C: `rank1_minus_stable_7761` (+0.130 → −1.172) and `rank1_plus_stable_522` (+0.060 → −0.130).
- There is no transition FAIL → PASS. Adding Λ_eff makes the condition stricter on two backgrounds and improves none.

## 4. D3 (per sign)

The leading kernel has 140 real branches in total: 29 backgrounds with 4 and 12 with 2. All roots are simple, and gcd(F_eff, tr K̂) = 1 on every background, so R-UNDEFINED = 0. All 280 branch decisions are certified at the first precision level. The smallest separations from 0 are: S − 1: 0.015; V2 margin: 0.023; H3: 4·10⁻⁴.

| | σ = +1 | σ = −1 |
|---|---|---|
| Crossing / no crossing | 114 / 26 | 103 / 37 |
| Poles (\|c\| > 1) | 87 | 87 |
| k_der status | 41 defined | 41 defined |
| k_der/M, min–max | 0.096–0.871 | 0.102–1.148 |
| V2 true / false / undetermined | **0 / 41 / 0** | **2 / 39 / 0** |
| H3 on all branches: true / false / undetermined | 0 / 41 / 0 | 6 / 35 / 0 |

V2 holds only for σ = −1, on two backgrounds:

- `rank1_minus_boundary_9745`: max(N + c²) = 1.759, k_der/M = 1.148, H3 on all branches, extended set C PASS;
- `rank1_plus_boundary_7900`: max(N + c²) = 1.936, extended set C FAIL.

On the anchor B3, V2 fails for both signs, with max(N + c²) ≈ 97 (σ = −1) and ≈ 109 (σ = +1). This is far from the bound 2, not a numerical edge effect.

## 5. Decision (reviewer, 2026-11-19)

- **Stage 2: accepted** as the authorised arithmetic.
- D4 gives no overall PASS for any of the three sets.
- D3: the global V2 criterion fails for both signs (σ = +1: 0/41; σ = −1: 2/41).
- Neither sign can be classified VIABLE-FOR-NUMERICS.
- NOT-VIABLE is not declared, because none of the registered NOT-VIABLE criteria is met by these results.
- **Model class: OPEN.**

**Scientific reading.** σ = −1 behaves better than σ = +1 on V2 and H3, but only on a few backgrounds. This does not overcome the limitation on the registered population. The sign of γ_u relative to g_φ is physical (Stage 1), and it matters here beyond leading order.

## 6. Files

- `proca_timelike_stage2.json`, `run_log.txt`, `RUN_CONTEXT.txt`, `SHA256SUMS`;
- post-hoc (read only, no recomputation): `stage2_per_background.csv` `98b05c66…`, `stage2_summary_tables.py` `de7698a9…`, `stage2_summary_tables.txt` `df849ecf…`;
- smoke test: `smoke_stage2_fix2.py` `d2ea9272…`, output `e2fd161a…`.
