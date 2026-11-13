# T3-fund HS2 — C4 result (registered spec, SPEC_HASH 9d52c927a9acd3bd)

Backgrounds (pre-designated): K3 (w,q,v,r)=(-1.0192715089, 0.2166884864, -0.3685933381, -0.5311944604); K4 (-0.1690732891, 0.2054735758, 1.2243113249, -0.7372513395).
Box/amplitude/T per spec (T=1000, eps=1e-3·max(|q|,|r|)); resolutions (N,dt) = (16384,0.025) and (32768,0.0125). Two light-field packets (one in each field), free control run on identical data.

| bg | run | norm_ratio_max | norm_ratio_final | total-energy drift max | fluct. energy drift max | min eig Z |
|---|---|---|---|---|---|---|
| K3 | control 16384 | 1.0 | 0.9999999988 | 1.5e-16 | 6.0e-8 | 1 |
| K3 | control 32768 | 1.0 | 0.99999999996 | 1.5e-16 | 6.0e-8 | 1 |
| K3 | coupled 16384 | 1.2983309173 | 1.0473439325 | 1.2e-16 | 5.8e-8 | 0.84583 |
| K3 | coupled 32768 | 1.2983309098 | 1.0473439342 | 1.2e-16 | 5.8e-8 | 0.84583 |
| K4 | control 16384 | 1.0 | 0.9999999988 | 2.1e-16 | 6.2e-8 | 1 |
| K4 | control 32768 | 1.0 | 0.99999999996 | 2.1e-16 | 6.2e-8 | 1 |
| K4 | coupled 16384 | 1.1738381094 | 0.9959509373 | 3.8e-16 | 1.1e-7 | 1.03974 |
| K4 | coupled 32768 | 1.1738381110 | 0.9959509076 | 1.9e-16 | 5.7e-8 | 1.03974 |

Registered rule: norm ratio ≤3 ✓ (max 1.298); total-energy drift ≤1e-9 ✓ (≤4e-16); control passes ✓; resolution agreement ≤1% ✓ (relative difference ≤1.3e-8 in norm_ratio_max); min eig Z > 1e-3 ✓.

**Outcome: K3 PASS, K4 PASS.**

Scope: two backgrounds only, one amplitude, T=1000. Total energy is dominated by the background energy, so the norm criterion is the substantive one. This does not establish T7 or a UV completion of B4. T7-a/T7-b not yet run.
