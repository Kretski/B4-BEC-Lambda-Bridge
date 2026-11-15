# Which energy? — naming correction and three separate results (v0.3.1)

## Naming correction
`E = ¼ ω z†H'z` (also ¼u†hu from the Legendre transform) is the Noether energy of the quadratic Lagrangian L2, which is
an expansion about a background linear in time (φ0 = wt+qx). It is therefore the second variation of **H − μ_φQ_φ − μ_ψQ_ψ**
(μ_φ=w, μ_ψ=v), not the lab-frame T00. Earlier text in `energy.py`/README that called it ⟨T00⟩ was imprecise; fixed.

## Independent consistency check (`b410/charges.py`, `tests/test_audit.py`)
From the NONLINEAR action, expanding to ε² with first-order waves: E_lab = ε²[T00], Q2_φ = ε²[∂L/∂φ_t], Q2_ψ = ε²[∂L/∂ψ_t].
Verified: E_lab − wQ2_φ − vQ2_ψ = Whitham energy = Legendre/h energy (rel 1e-8 in tests; ≤7.6e-7 over all 41 backgrounds at the
sampled nodes, QZ-limited at large x). So the three descriptions are consistent, and E_lab is a *different* quantity.
Linear terms in second-order fields are dropped (total derivatives with constant coefficients).

## Three separate results per background (`results_e3/e3_energy_triplet.csv`, 41 backgrounds, 121 nodes each)
| | lab δ²H (E_lab) | δ²(H−μQ) (E_shift = h) | dynamic reality |
|---|---|---|---|
| 5 backgrounds with negative branch (5713, 4487, 2020 minus; 7900, 3437 plus) | positive for every sample (min E_lab/|E_shift| 0.75–0.95) | negative on one slow scalar branch (54–122 of ~120 nodes) | real at all x, except #4487 (complex for x>0.1887 via Krein collision) |
| all other real-root backgrounds (23 + anchor) | positive | positive | real |
| 12 negative controls | – | h indefinite | complex at every x |

## My reading (partly different from the proposed table — please check)
For waves with k≠0 the first-order charge change vanishes, and conservation of Q requires the wave's second-order charge to be
compensated by the homogeneous (k=0) mode, which costs μ·ΔQ. Hence the energy at fixed charges equals δ²(H−μQ) restricted to δ¹Q=0,
i.e. the same E_shift (this is the usual Bogoliubov statement that excitation energies are ħω in the μ-frame). Under this reading the
"fixed charges" and "fixed μ" rows coincide at quadratic order, and the lab energy E_lab is the relevant one only for a reservoir that
can absorb the wave's charge at zero chemical potential. This argument is not tested against a full (nonlinear/mean-field) calculation.

## Decision recorded
The five backgrounds are NOT rejected and NOT declared stable: they are marked "negative shift-symmetric energy"; physical
admissibility depends on the charge-conservation conditions of the intended system. Negative energy ≠ complex frequency.
