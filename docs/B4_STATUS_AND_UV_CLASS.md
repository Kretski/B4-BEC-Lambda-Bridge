# Status of B4 and the UV class to search (decision note)

Purpose: fix what kind of theory B4 is before testing further UV candidates. Nothing here is a new result; each line cites what was established and where.
No option is chosen in this note.

## Established inputs
1. Operator level: phi_x^4 and (phi_t phi_x)^2 contain boost weights +-4; every Lorentz-invariant quartic has weight 0. B4 with constant p, s, kappa, Lambda
   in a fixed frame is therefore not a Lorentz-invariant operator set (SCALAR_MEDIATOR_RESULT.md, HeavyScalar_TwoField).
2. B4 can be written covariantly with a frame vector u (n = eps u). With u taken from the backgrounds (d psi, or d phi + lambda d psi) it does not reduce to
   lab-frame B4 on the manifest: X_psi > 0 on 24/41, no constant lambda (B4 Lorentz classification).
3. Bare Lambda > 0 gives a local speed sqrt(1 + 3 Lambda q^2) > 1 in B4 units, independent of any mediator (b4_causality_decomposition.md).
4. Lorentz-invariant heavy scalar: front = 1, stable on most sampled backgrounds, effective operator +O^2/(2M^2) (subluminal sign), T1/T7 negative
   (HEAVY_SCALAR_POPULATION_REPORT.md).
5. Frame-dependent mediators tested so far (Proca, scalar, scalar + dynamical eta) fail T3 against the Minkowski cone on the manifest backgrounds.

## Option 1: B4 is a Lorentz-violating effective theory (the lab frame is physical)
Meaning: a third physical agent (matter rest frame, reservoir, medium, detector) fixes u. The unit speed in B4 is then a property of that medium.
Requirements on a UV mechanism:
- T6 becomes the central test: the origin and dynamics of u must be stated and u must have its own stability and causality check.
- T3 needs a decision: which cone is the causal one? If the unit speed of B4 is the fundamental light cone, input 3 already violates T3 for Lambda > 0, before any UV completion.
  If the unit speed is an emergent speed of a medium (as the sound speed in a BEC, where Bogoliubov dispersion is supersonic and physical), T3 must be
  tested against the fundamental cone of the underlying system, and "front > 1" in B4 units is not by itself a failure.
- Consequence for existing results: the Proca and scalar-mediator T3 verdicts were obtained against the B4 unit cone. They stay valid as stated; under the
  second reading they would have to be re-evaluated against the fundamental cone, with the ratio c_fundamental / c_B4 as a new parameter fixed in advance.
- T7 is in principle attainable (the operators are allowed once u exists); T8 would be a frame-dependent, direction-dependent dispersion.

## Option 2: B4 is a Lorentz-invariant theory seen around a background
Meaning: no external frame; p, s, kappa, Lambda are not constants but functions of the invariants X_phi, X_psi, Y of the background.
Requirements on a UV mechanism:
- T1/T7 must be reformulated: match background-dependent coefficient functions, not constants. The B3 anchor numbers are then valid only on its own background.
- T3 holds automatically for Lorentz-invariant principal symbols with Z > 0 (input 4).
- Sign constraint: for the heavy scalar the generated operator has the subluminal sign; a local speed above 1 (input 3, Lambda > 0) is not produced.
  Positivity bounds from the literature (not re-derived here) suggest this is general for standard Lorentz-invariant UV completions. If so, Option 2 is
  compatible with a Lambda-like dispersion (Lambda_eff from beta) but not with bare Lambda > 0 as a local superluminal speed.
- T8: sign and size of the cubic dispersion term are predictions; the sign is fixed.

## What each choice implies for the next step
| Choice | Next UV class | First task |
|---|---|---|
| Option 1, B4 unit speed = fundamental cone | none can pass T3 with Lambda > 0 (input 3) | decide whether to drop Lambda > 0 or the Minkowski cone |
| Option 1, B4 unit speed = emergent medium speed | medium / analogue-type mechanisms (BEC line becomes central again) | fix c_fundamental / c_B4 and re-state T3 before any run |
| Option 2 | Lorentz-invariant UV (heavy scalar as reference) | reformulate T1/T7 with coefficient functions; test against Lambda-model dispersion sign |

## Not decided here
Which option describes the intended physics of the Lambda-model (phi = Lambda-field, psi = matter). This needs the physical identification of phi, psi and of the frame.
