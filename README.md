# B4-BEC-Lambda-Bridge

Effective Lambda Generation and Operator-Level Matching in Two-Component Quantum Fluids.

## Research objective

This repository investigates whether a two-component quantum fluid can provide a physically consistent analogue or realization of the effective nonlinear correction represented by Lambda in the B4 effective field theory (EFT).

The central question is whether an effective spatial quartic interaction can be derived from the BEC effective action and consistently matched to the B4 coefficient.

The project does not assume that such a correspondence exists.

## Main research questions

1. Can an effective coefficient Lambda_eff be derived from explicitly defined BEC coefficients?
2. Under what field normalization and background conventions can Lambda_eff be compared with the bare B4 parameter Lambda?
3. Are the cubic and quartic fluctuation operators compatible between the two descriptions?
4. Can one fixed set of microscopic parameters reproduce a family of B4 backgrounds?
5. Which results are analytical, numerical, approximate, or inconclusive?

## Theoretical components

### B4 EFT

The working two-field Lagrangian is:

L = 1/2 w^2 - 1/2 q^2 - Lambda/4 q^4
  + 1/2 v^2 - 1/2 r^2 - 1/2 m^2 psi^2
  + p/2 w^2 q^2 - kappa w q v r
  + s/2 v^2 r^2.

Definitions:

- w = d_t(phi)
- q = d_x(phi)
- v = d_t(psi)
- r = d_x(psi)
- p = lambda_phi
- s = lambda_psi

### BEC model

The starting point is the two-component Gross-Pitaevskii model, with possible Lee-Huang-Yang (LHY) corrections.

Mean-field and LHY contributions will be kept separate.

## Research stages

1. Derive Lambda_eff from explicitly defined BEC coefficients.
2. Establish field normalization and dimensional consistency.
3. Compare Lambda_eff with the B4 bare Lambda.
4. Derive and compare the complete cubic and quartic fluctuation sectors.
5. Test a family of backgrounds with fixed couplings and a predefined parameter dictionary.
6. Document limitations and reproducibility.

## Methodological rules

- No assumption of a positive matching result.
- No background-dependent redefinition of fixed microscopic couplings to force agreement.
- Explicit field normalization and reference frame.
- Dimensional checks for every matched coefficient.
- Separate analytical identities from numerical evidence and approximations.
- Preserve negative and inconclusive results.

## Current status

Repository initialized. The exact definitions of the BEC coefficients A_i and the corresponding Lambda_eff mapping remain to be specified and independently checked.

No BEC-to-B4 operator matching is claimed.
## UV candidate assessment (exploratory)
See `docs/UV_ACCEPTANCE_PROTOCOL.md` (criteria T1-T8), `docs/UV_CANDIDATE_ASSESSMENT.md` (evidence table with scope of each verdict) and `docs/SCALAR_MEDIATOR_RESULT.md`.
Proca is kept as a negative control. No candidate reproduces the B4 operator structure; BEC Stage 2 is deferred. Verdicts are limited to the named constructions.
