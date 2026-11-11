# UV candidate assessment (T1-T8)

Criteria: see UV_ACCEPTANCE_PROTOCOL.md (T1-T6 fixed before the scalar-mediator runs; T7-T8 added later). Every verdict holds only for the named construction and the sampled backgrounds.
Evidence codes: **A** derived analytically; **N** checked numerically on a stated sample; **C** true by construction; **-** not checked; **X** negative result in that class.

| Candidate | T1 (p,s,kappa) | T2 (real, hyperbolic) | T3 (front <= 1) | T4 (perturbative) | T5 (energy) | T6 (frame) | T7 (operator match) | T8 (observable) | Scope of evidence |
|---|---|---|---|---|---|---|---|---|---|
| Proca mediator (B4.10-B4.12) | A: exact factorisation, kappa^2 = ps | N: scanned backgrounds (6 become complex at Lambda=0) | X, N: front > 1 on all sampled backgrounds (min 1.348 at Lambda=0, anchor 4.84 / 8.26) | see B4_12_REPORT.md | see B4.10 audit | external lab frame | - | - | the rank-one Proca construction, 41 backgrounds + population cells |
| Scalar mediator, lab frame, Lambda = 0 | A: same F (det K/Q = -M^6 X^4) | N: 399/400 on J0 = 0 sample | X on manifest (0/41, min 1.008); C/A on J0 = 0 (front = 1, 400/400) | C on sample (sqrtX <= 0.3) | N: H2(k) >= 0 on 400 J0 = 0 backgrounds | external frame | - | - | J0 = 0 is a background condition, not a property of the theory |
| Scalar mediator + dynamical eta (u^mu) | - (Lambda term not included) | N: 17/17 real | X: A (front = 1 only if sigma_bar = 0, incompatible with stationary eta except degenerate) and N (0/17, min 1.007) | - | - | eta stationary on 17/41 backgrounds only | - | - | this action (kinetic f/2 (d eta)^2, no potential) |
| HeavyScalar_LorentzInvariant (1 light field) | n/a (one field) | N: scanned X in [1e-3,1e3], 4 backgrounds | A: principal symbol Lorentz invariant (front = 1) ; N roots | - | A: Z > 0 on tested backgrounds | not needed (Lorentz invariant action) | - | - | one light field, tree level, coupling chi (d phi)^2 |
| HeavyScalar_TwoField | X, N: best-fit B4 constants leave a residual of 18-60 % of the deviation from free; fitted kappa^2/(ps) = 9.7-29 | A: det K = (c^2-1) x quartic; T5-static implies real frequencies. N: 474 / 500 / 500 of 500 backgrounds (seed 21), 478 / 500 / 500 (seed 22) | A: principal symbol diag(Z,1)(omega^2-k^2), front = 1 iff Z > 0. N: Z > 0 on 500/500 in all sets | official T4 **not checked**; T4-proxy (max abs(eig Z - 1) <= 0.3): 432 / 500 / 500 (seed 21) | T5-static (A + N): 465 / 500 / 500; T5-Hamiltonian: H2 positive definite exactly on the T5-static backgrounds; T5-Krein: no negative-energy mode inside T5-static, 10 backgrounds with one outside (seed 21). Time-domain: energy conserved to 5e-13 of total on one background | not needed | X, A: Lorentz-invariant quartics have boost weight 0, B4 phi_x^4 and (phi_t phi_x)^2 do not | - | three coupling sets, backgrounds uniform in [-1,1]^4, M = 1, k in {0.01, 0.1, 1, 10} for H2; not a statement about all backgrounds |
| BEC mean field (control) | - | - | - | - | - | n/a | A: GP has a spatial quartic with Lambda_BEC < 0, no kappa operator, tied to the cubic term | - | single component symbolic; two-component Stage 2 deferred |

## Findings that are not tied to one candidate
- Three different things must not be conflated: scattering (time delay ~ A^2, no lasting change of dispersion, A/N), a stationary-background dispersion correction (present, measured), and B4 operator correspondence (not achieved).
- Lambda_eff = 2 beta / c0 > 0 appears for all the light modes of the J0 = 0 scalar mediator (800/800) and of the heavy scalar (two modes). The same sign follows for any normal dispersion (omega ~ ck + beta k^3 with |c| increasing toward its front), e.g. the Bogoliubov quantum-pressure term, so the sign alone does not discriminate between mechanisms and is not the bare B4 Lambda (frame dependent, undefined at c0 = 0).
- Hypothesis B (dynamics of the overlap): no lasting dispersion in vacuum for HeavyScalar_TwoField (shift -1.42, dk 1.2e-6).
- The lab frame is not fixed by the backgrounds of the manifest: X_psi > 0 on 24/41, X_phi > 0 on 27/41, no constant lambda gives u from d phi + lambda d psi.

## Open items
T8 (observable), T7 for any Lorentz-violating mechanism, T4 for the heavy scalar in the sqrtX sense, other couplings and backgrounds for the heavy scalar, quantum/loop effects (positivity bounds recalled from literature, not re-derived).
