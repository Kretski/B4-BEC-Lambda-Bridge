# T3-fund C4 (time domain) — result
Spec: t3fund_c4.py, SPEC_HASH f2f43852bf62ce62 (registered in d7c33d4). Nothing changed after registration.
Backgrounds: SM = A1 seed 11 (J0=0); HS = B2 (w=1,q=0.3,g=0.4,M=1). Amplitude 1e-3, T=1000, Lbox 1024.

| Model | run | N, dt | norm_ratio_max | drift_total_max | drift_fl_max (diag.) |
|---|---|---|---|---|---|
| SM | coupled | 16384, 0.025 | 1.10778288 | 1.3e-16 | 2.5e-8 |
| SM | coupled | 32768, 0.0125 | 1.10778287 | 2.6e-16 | 5.1e-8 |
| SM | free control | 16384 / 32768 | 1.0 / 1.0 | 1.3e-16 / 1.3e-16 | 2.5e-8 |
| HS | coupled | 16384, 0.025 | 1.02132411 | 2.9e-16 | 5.8e-7 |
| HS | coupled | 32768, 0.0125 | 1.02132413 | 2.9e-16 | 5.8e-7 |
| HS | free control | 16384 / 32768 | 1.0 / 1.0 | 2.0e-16 / 2.0e-16 | 3.8e-7 |

Rule check (registered): norm ratio <= 3: yes (1.108, 1.021). Drift <= 1e-9: yes, all runs (~1e-16). Control meets target: yes.
Convergence of norm_ratio_max between resolutions: relative difference ~1e-8 (SM), ~1e-8 (HS) <= 1%: yes.
Outcome: SM on A1/seed 11: PASS. HS on B2: PASS.

Limits of this PASS (stated, not new criteria):
- E_total is dominated by the background energy (E0 = 1761 and 788) while the fluctuation energy is ~1e-5 and ~4e-7. A drift of 1e-9 in E_total therefore constrains the fluctuation energy only to ~1e5 resp. ~1e3 relative; drift_fl is at round-off level (E_total round-off / E_fl), so it is not a meaningful measurement. The energy criterion mainly tests the integrator on the background; the norm criterion (bounded, converged) is the substantive result.
- One background per model, one perturbation shape, amplitude 1e-3, T = 1000. This is evidence for those backgrounds only.
- C4 PASS does not establish T7 or UV completion of B4. The IR-instability result of Stage 1 (sets A2, A3, B1 failures) stands.
