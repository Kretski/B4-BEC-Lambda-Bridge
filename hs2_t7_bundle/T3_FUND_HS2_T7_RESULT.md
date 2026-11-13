# T3-fund HS2 — T7 result (registered spec, SPEC_HASH 9d52c927a9acd3bd; stage commit 0c6f8b8)

## T7-a (analytic, no fit)
As registered in spec §6: the generated operators O²/(2M²) are Lorentz scalar contractions (boost weight 0); the B4 operators Λφ_x⁴ and p(φ_tφ_x)² are not Lorentz invariant (light-cone weights φ_x⁴: +4,+2,0,−2,−4; (φ_tφ_x)²: +4,0,−4). Hence no operator-level match is expected from this tree-level Lorentz-invariant class without extra symmetry breaking. This is an argument about operator structure, not a proof about the full UV theory.

## T7-b (numerical), thresholds as registered
HS2 target: k→0 local polynomial P(c)=D(c)·F(c), D=1−c², F=M²D·detZ̄−4sᵀadj(Z̄)s, monic-normalised (reduces to (c²−1)² for C'=0). B4 polynomial as in heavy_scalar_two_field_matching.py. Per set: first 15 H1 draws with F real; joint fit of (p,s,κ,Λ), 60 random starts (seed 71). Residual per background = ‖(P_B4−P_HS2)[:4]‖/‖(P_HS2−free)[:4]‖; MATCH requires ≤1e-6 on all 15.

| set | max residual | rms residual | min per-bg | control B4→B4 (≤1e-12) | outcome |
|---|---|---|---|---|---|
| K1 (rank one) | 0.943 | 0.833 | 0.462 | 1.8e-15 ✓ | NO-MATCH |
| K3 | 0.816 | 0.519 | 0.242 | 3.6e-16 ✓ | NO-MATCH |
| K4 | 0.875 | 0.702 | 0.284 | 4.8e-16 ✓ | NO-MATCH |

Pre-registered expectation (NO-MATCH for all three sets) is confirmed. Residuals are of order the deviation itself (not near the 1e-6 threshold by any margin), the control recovers to machine precision, so the fit machinery is not the limiting factor. No degenerate or special-background MATCH occurred, so no extra check was triggered.

Notes: (i) the residual normalisation is per background (my implementation of "normalised by the deviation from the free polynomial"); the earlier script used an rms ratio; both are O(1), far above 1e-6. (ii) K1 is rank one and equivalent to one field with g=2a, a test of the implementation, not an independent case. (iii) The 60-start fit is a numerical search over 4 parameters, not a proof that no parameters exist; together with T7-a it is the registered evidence.

Scope: tree level, this Lagrangian, sampled backgrounds only. This does not establish or exclude a UV completion of B4 in general.
