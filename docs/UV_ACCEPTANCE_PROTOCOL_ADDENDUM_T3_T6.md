# Addendum to the UV acceptance protocol: causal cone (T3) and preferred frame (T6)

Date: 2026-11-12. The original T1-T6 text in UV_ACCEPTANCE_PROTOCOL.md is not modified. Existing numerical results are not changed, only relabelled.
Working interpretation chosen by the project lead: option (b) of B4_CAUSALITY_ANALYSIS.md, i.e. B4 is a Lorentz-violating effective theory with a
preferred observer field (static observer / ZAMO) and Lambda > 0 is kept. This is a working hypothesis, not an established fact.

## Three speeds
- c_B4: the reference (unit) speed of the effective theory.
- c_fund: the fundamental causal limit.
- c_front: the front velocity of a given UV realisation (large-k limit of omega/k).

## T3 is split
- T3-B4: c_front <= c_B4. All T3 verdicts obtained so far (Proca, scalar mediator, scalar + dynamical eta, heavy scalar) are T3-B4 verdicts and keep their values.
- T3-fund: causality with respect to the fundamental structure. Its definition depends on one open decision (below). Until that decision is registered,
  T3-fund is "not evaluated" for every candidate, and no candidate is accepted.

## Registered decision (project lead, 2026-11-12): alternative (iii), preferred-foliation causality
> Causality Convention: The B4-Lambda framework adopts preferred-foliation causality defined by the observer field u^mu. No finite Lorentz-invariant
> upper bound on spatial front velocity is imposed at the effective-theory level. T3-fund evaluates the well-posedness and stability of evolution with
> respect to the preferred time. The physical origin, admissibility and domain of definition of the observer field are assessed independently under T6.

Accepted with it: c_B4 = c (speed of light in the normalisation used); no separate finite c_fund > c is introduced; unbounded front or group velocity is
allowed provided evolution is well posed in the preferred time. This is a working hypothesis for the interpretation, not an established fact.
Alternatives considered and not adopted: (i) finite cone c_fund = c_B4 (incompatible with Lambda > 0 for every k > 0); (ii) finite cone c_fund = R c_B4, R > 1
(requires a fundamental speed above the speed of light).

## T3-fund: what is checked
1. Real spectrum of linear excitations.
2. Well-posed Cauchy problem with respect to the preferred time.
3. Hyperbolicity, or another appropriate notion of well-posedness of the evolution system.
4. Bounded growth of solutions at fixed time (no growth rate increasing without bound in k).
5. No instability that makes the evolution physically unacceptable.
Unbounded propagation speed and unbounded amplitude growth are different properties; only the second fails T3-fund.
Note: items 1 and 3 overlap with T2. Under this convention T3-fund no longer bounds speeds, so the discriminating weight moves to T2, T5 and T6.
Numerical tolerances for T3-fund are not fixed yet; they must be registered before the first T3-fund run. Status for all candidates: not evaluated.

## T6 as a stand-alone task (now the central open physical problem)
It is not sufficient to specify u^mu mathematically. To be established:
1. What the observer field is physically.
2. How it is determined, dynamically or geometrically.
3. What its time foliation is (existence of a global preferred time; hypersurface orthogonality of the congruence).
4. What happens where a given congruence is not defined. Static observers do not exist inside the ergoregion of Kerr; ZAMOs remain defined outside the
   horizon. One frame cannot be assumed to cover the whole spacetime without an additional construction.
5. Its own stability and well-posedness if it is made dynamical. The dynamical-eta test already done is one negative data point for one construction.

## Acceptance rule
No UV candidate is called accepted until T3-fund (in the registered form) and T6 are evaluated for it.
