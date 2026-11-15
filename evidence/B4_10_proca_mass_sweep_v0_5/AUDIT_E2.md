# Audit of Experiment 2 (v0.2.x) — independent derivation

Everything below is in `b410/audit_hamiltonian.py` + `tests/test_audit.py` (59 tests total).
It does **not** use H_B49, H', the null-vector Noether code, or `energy.terms`.

**Source of truth** (T1.7 text): L = ½φ_t²−½φ_x²−(Λ/4)φ_x⁴ + ½ψ_t²−½ψ_x² + ½(A1_t−A0_x)² + ½M²(A0²−A1²) + A1(g_φφ_tφ_x + g_ψψ_tψ_x),
i.e. J_φ^1 = φ_tφ_x (so (1/2M²)J·J = (p/2)(φ_tφ_x)², κ=−g_φg_ψ/M²). Background φ0=wt+qx, ψ0=vt+rx, A1=J0/M². The identification q=φ0_x, w=φ0_t follows from the 3Λq² gradient term.
That the action is the T1.7 one is taken from the T1.7 text, not re-proven.

| Audit point | Result |
|---|---|
| Hamiltonian independent or reconstructed? | Now independent: sympy second variation of the nonlinear L reproduces every coefficient of the reconstructed L2 (rel 1e-12); the linear term has no undifferentiated-field part (background stationary, A1_bg = J0/M² confirmed). |
| A0 equation | ∂L/∂A0=0 ⇒ A0 = −P_a'/M² (Hamiltonian form); reduced A-sector ℋ_A = ½(1+k²/M²)P_a² + ½M²a². No extra root. |
| Reduced Hamiltonian & Hamilton equations | 6×6 Hermitian h(k) built by Legendre transform with A0 eliminated. Eigenfrequencies of J·h equal the six roots to ≤6e-11 relative (6e-16 for x≲10). |
| Energy factor and sign | ¼u†hu equals E=¼ωz†H'z to 1e-13 (3e-9 at x=10³, QZ-limited). Factor ¼ = ½ (Hamiltonian density) × ½ (cycle average). |
| Krein / positivity | h(k) is positive definite on a 400-point log grid x∈[1e-4,1e4] (min eigenvalue ∝k² as k→0). Positive-definite h with simple roots ⇒ all six E_n>0 and real roots — this is a stronger statement than the sign of z†H'z. Negative control (a background with complex roots) gives indefinite h and complex Hamilton frequencies. Grid check, not an analytic proof for all x. |
| Which energy? | E is the second variation of H−wQ_φ−vQ_ψ (relative equilibrium of the shift symmetry), the right object for linear stability. The lab-frame second-order T00 differs by w g_φ⟨a1 φ_x⟩+v g_ψ⟨a1 ψ_x⟩; it is also positive (min ratio 0.90 over 97 x-points × 6 branches). Mean-flow second-order terms are neglected. |
| Invariance of energy share | The A-sector share computed from Noether terms and from the Hamiltonian representation agree to 1e-8. The scalar/interaction split between φ,ψ and the cross terms is representation dependent and is not used for conclusions. The Euclidean amplitude share is not equivalent (documented). |
| Tracker with close roots | Roots are simple for all x>0 (exact discriminant); smallest relative adjacent gap 0.33 (x≈0.28). Overlap tracking = ω-ordering at every node. Euclidean overlap margins are small (non-orthogonal eigenvectors of a polynomial pencil) so tracking rests on assigned overlap ≥0.98 + agreement with ω-order. Subspace tracker exercised on a synthetic exact degeneracy. |
| Resolution | Exchange branches {0,1,3,5} and crossing x identical for n0=241 and 481 (≤1e-3 relative) and within 3% for n0=61. |
| eft_domain | Reconstructed (four X<0.95, T1.7.1 `rank_one_grid`); original T1.7.3 scan code unavailable. B3 row passes 0.95, fails 0.5. |

## Not settled by this audit
* Causality: phase/group velocities exceed 1 (max |v_g|≈8.3 in these units); real ω and E>0 do not imply causality. Front velocity not computed (E4).
* UV: g∝M_A; perturbativity not examined (E5).
* One background only (E3).
