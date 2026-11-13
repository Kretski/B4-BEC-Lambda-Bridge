# T3-fund specification: HeavyScalar_TwoLight (HS2) — C0–C4 and T7 test
Status: PROPOSED for registration. No HS2 population evaluation, no C3 run and no C4 run has been performed. Only done: (i) the algebra identity check `t3fund_hs2_check_poly.py` (sympy, one random point + 20 random points, max rel. error 6e-14), (ii) manifest generation `t3fund_hs2_manifest.py` (random draws and analytic pre-classification, no root computation of det K). After registration nothing below changes; deviations are reported as deviations.

## 1. Model (M=1, signature +−, (∂f)·(∂g)=f_t g_t − f_x g_x)
L = ½(∂φ)² + ½(∂ψ)² + ½(∂χ)² − ½M²χ² + χ·O,  O = a(∂φ)² + b(∂ψ)² + e ∂φ·∂ψ = Σ C'_ij ∂f_i·∂f_j,  C' = [[a, e/2],[e/2, b]], f=(φ,ψ).
Equations: Z_ij □f_j + 2C'_ij ∂_μχ ∂^μ f_j = 0, Z = I + 2χC';  □χ + M²χ = O.   Hamiltonian: H = ½Z_ij(f_i,t f_j,t + f_i,x f_j,x) + ½χ_t² + ½χ_x² + ½M²χ².
Background: f0=(w t+q x, v t+r x); χ̄ = O(w,q,v,r)/M² (exact constant solution, C0 residual of all three equations ≤1e-14 in double); Z̄ = I+2χ̄C'.
Coupling sets (fixed): K1=(a,b,e)=(0.2,0.2,0.4) [C' rank one: O=a(∂(φ+ψ))², reduces to one light field with g=2a; kept for continuity and as implementation test]; K3=(0.3,−0.2,0.1) [full rank, indefinite C']; K4=(0.15,0.25,0.2) [full rank, positive definite C'].

## 2. C1/C2 (analytic, registered before results)
Fluctuation (α,β,ξ), plane wave e^{i(kx−ωt)}, c=ω/k, ε=1/k², σ̃_j=c w_j+q_j (w=(w,v), q=(q,r)), s=C'σ̃, D=1−c². Symbol K(ω,k) 3×3 (light block Z̄(k²−ω²), light–χ entries ∓2i k s, χχ: k²−ω²+M²). Closed form (verified):
  det K = k⁶ · D · Q(c,ε),  Q = (D+εM²)·D·det Z̄ − 4ε·sᵀ adj(Z̄) s.
- Principal part: diag(Z̄,1)(ω²−k²). Two light fields and χ share the luminal cone; couplings enter at first derivative order only.
- First-order reduction: U=(α,α_t,α_x,β,β_t,β_x,ξ,ξ_t,ξ_x) (9 components): U_t + A U_x = B U, A real symmetric constant (blocks [[0,−1],[−1,0]] for each field, zero rows for amplitudes), B constant (contains Z̄⁻¹ and the first-derivative couplings; requires Z̄ invertible). Hence C2: real eigenvalues, complete orthogonal eigenvectors (also at repeated roots); C1: γ(k) ≤ λmax((B+Bᵀ)/2) for all k incl. k=0 (computed in mp, 50 digits, per background).
- Ghost condition: Z̄ must be positive definite (H positive definite for the principal part). Z̄ not positive definite ⇒ class FAIL-GHOST (reported separately, not part of the T3 PASS fractions, no C3/C4 attempted).
- k→∞ (analytic): near c=±1, c = ±1 + η, η = ±√ε·√(s_±ᵀ Z̄⁻¹ s_±), s_±=C'(±w+q); i.e. ω = ±k ± √(s_±ᵀ Z̄⁻¹ s_±) + o(1). For Z̄>0 these are real; the exactly luminal pair c=±1 belongs to the factor D (all k). The two heavy roots → ±1 too.
- k→0 (IR, analytic): F(c) = M²(1−c²) detZ̄ − 4 sᵀadj(Z̄)s, quadratic in c (coefficients A2,A1,A0 in the manifest generator); complex roots ⇔ F_disc<0.
- k=0 (P2a): ω⁴·detZ̄·[(M²−ω²) + 4 s_wᵀ Z̄⁻¹ s_w] with s_w=C'w; four zero modes are the ZERO MODE class (expected), heavy mode ω_h² = M²+4 s_wᵀZ̄⁻¹s_w must be ≥0 (else FAIL-UNSTABLE at k=0). Computed numerically from the 3×3 matrix; the formula is a cross-check.

## 3. Populations (manifest `t3fund_hs2_manifest.json`, MANIFEST_HASH in §9)
Draws (w,q,v,r) ~ U[−1.5,1.5]⁴, numpy default_rng(31) for K1, 32 for K3, 33 for K4, 300 raw per set, M=1, in generation order.
- H1: raw draws with Z̄ positive definite (896: K1 296, K3 300, K4 300) — population for C3.
- H2: raw draws with Z̄ not positive definite (4, all K1) — expected FAIL-GHOST.
- H3 (control): first 50 draws in H1 order with F_disc<0 (analytic IR complex roots) — expected C3 FAIL-UNSTABLE. This is an analytic pre-selection; it tests detection, not a physical claim.
- Anchors (w,q,v,r)=(1,0.3,0.8,−0.2): B2h K1, B2k K4, B2s K3, B2z (a,b,e)=(0.2,0.2,0). Predicted C3 PASS (F_disc>0 in all four, Z̄min 0.90–1.12).

## 4. C3 (same as Stage 1)
P2 grid 61 pts k/M∈[1e-3,1e3] + 31 pts in [1e3,1e6]; roots of Q (and the exact factor D) in mp at 50 digits; classes: r≤1e-7 real, 1e-7<r≤1e-5 unresolved, >1e-5 candidate complex, confirmation at 50 digits for r>1e-12; P2a at k=0 as in §2; P6 three-band rule with u(k)=|γ50−γ100| in arbitrary precision (mp.dps=120 for the difference), margin 10, ≥3 consecutive refined points. Verdicts: PASS / FAIL-UNSTABLE / UNRESOLVED / FAIL-GHOST. The exactly luminal factor D=1−c² (roots c=±1 at all k) is classified separately from Q, so that the luminal pair is not misclassified. The same normalisation of U used to build B is used in the energy bound. The polynomial implementation is checked against the explicit K matrix (as in t3fund_stage1_check_poly.py) and the first-order matrices A,B against it before use.
Registered predictions (falsifiable, results will be reported whether or not they hold):
- Pred-1: on H1, C3 FAIL-UNSTABLE ⇔ F has complex roots (F_disc<0). Hypothesis from Stage 1 (740/740 there).
- Pred-2: no H1 background has complex roots only for k/M>1 (every unstable one is already unstable in the IR); in the band k/M∈[1e3,1e6] all roots real for all of H1 (§2 expansion).
- Pred-3: H3 50/50 FAIL-UNSTABLE; H2 4/4 FAIL-GHOST; anchors 4/4 PASS.
- Pred-4 (implementation test): for K1 the roots coincide with one-field HS(g=2a) on background ((w+v)/√2,(q+r)/√2) (rotated combination θ=(φ+ψ)/√2) plus the free luminal pair of (φ−ψ)/√2, to 1e-12.
Count of Pred-1 exceptions is reported; no criterion changes.

## 5. C4 (time domain)
Designated backgrounds (manifest `C4_designated`, rule: first draw in H1 order with F real, per coupling set; run only if C3 PASS, else the next on the same rule, recorded): 
 K3: (w,q,v,r)=(−1.0192715089, 0.2166884864, −0.3685933381, −0.5311944604) (full values in manifest);  K4: (−0.1690732891, 0.2054735758, 1.2243113249, −0.7372513395).
Method as registered C4: periodic box 1024, pseudo-spectral D, RK4, D(D(f)), T=1000, records every 5, base (16384,0.025), convergence (32768,0.0125). State: fluctuations (a,pa,b,pb,xi,ct); Pt_j=w_j+p_j, Px_j=q_j+D(f_j), χ=χ̄+xi, Z=I+2χC'.
 f_tt = D(D f) − Z⁻¹·[2C'_ij (ct·Pt_j − D(xi)·Px_j)],   ct_t = D(D xi) − M²χ + Σ C'_ij(Pt_iPt_j − Px_iPx_j).
Energy E=dx Σ[½Z_ij(Pt_iPt_j+Px_iPx_j)+½ct²+½D(xi)²+½M²χ²]; diagnostic E_fl=E−E_bg.
Perturbation: Gaussian s=15, k0p=0.5, centre 0, amplitude eps=1e-3·max(|q|,|r|); a right-moving (pa=−D a), b left-moving (pb=+D b); xi=ct=0. Norm N=√(dx Σ(a²+b²+xi²)); criterion max N/N(0) ≤ 3. Free control: a=b=e=0, χ̄=0, same gradients/perturbation/resolutions. Additional monitor: min eigenvalue of Z at each record; if ≤1e-3 the run is INCONCLUSIVE (Z-degenerate) and reported.
Outcome rule as in registered C4 (incl. Clarification 1): PASS if at both resolutions N-ratio≤3, drift target (≤1e-9 for total energy) reached, control reaches target, max N-ratio agrees within 1%; INCONCLUSIVE if control misses target, the 1% agreement fails, drift target is missed with control passing, or Z-degenerate; FAIL-UNSTABLE if N-ratio>3 with agreement. No FAIL on energy alone.
Order: K3 coupled+control (both resolutions), K4 coupled+control, then drift/norm/convergence check, report, commit. Estimated ~4 h with two parallel lanes.
Stated limitation (as in the C4 report): the total-energy drift is dominated by the background energy and constrains the fluctuation energy only weakly; the norm criterion is the substantive one.

## 6. T7 (operator matching) — what is tested; existence of a match is NOT assumed
Hypothesis H_T7: the k→0 local polynomial F(c) of HS2 on a coupling set equals the B4 local polynomial F_B4(c;p,s,κ,Λ) (as implemented in heavy_scalar_two_field_matching.py) with constants (p,s,κ,Λ) independent of the background.
- T7-a (analytic, no fit): At O(1/M²), the generated operators are O²/(2M²), constructed from Lorentz scalar contractions (∂f_i·∂f_j), and therefore have boost weight zero. The B4 operators Λφ_x⁴ and p(φ_tφ_x)² are not Lorentz invariant: their light-cone decompositions (φ_x ∝ ∂₊φ−∂₋φ) contain nonzero boost-weight components (φ_x⁴: weights +4,+2,0,−2,−4; (φ_tφ_x)²: +4,0,−4). Consequently, no exact operator-level matching to these B4 terms is expected from this Lorentz-invariant tree-level UV class without additional symmetry-breaking structure. This is recorded as the expected, falsifiable result.
- T7-b (numerical): per coupling set (K1,K3,K4) take the first 15 draws of H1 with F real; fit (p,s,κ,Λ) jointly (60 random starts) to minimise the residual normalised by the deviation from the free polynomial (same residual as the earlier test); control B4→B4 must give ≤1e-12. Registered threshold: MATCH if residual ≤1e-6 on all 15 for a set; otherwise NO-MATCH (value reported). Pre-registered expectation: NO-MATCH for all three sets.
- A MATCH, if it occurred, would be reported as unexpected and checked (fit degeneracy, rank-one case K1) before any claim.

## 7. Scope and non-claims
Tree level, this Lagrangian, sampled backgrounds only. PASS in C3/C4 does not establish a UV completion of B4 or T7. C3 failures on the IR quadratic F are IR effective-theory instabilities, not UV pathologies (Stage 1 finding, to be tested here by Pred-1/2).

## 8. Files and untouched items
t3fund_hs2_check_poly.py, t3fund_hs2_manifest.py, t3fund_hs2_manifest.json, this document. tests/test_dimensions.py and theory/operator_matching.py are not touched. Implementation scripts (t3fund_hs2_stage1.py, t3fund_hs2_c4.py) will be written after registration to implement §2–§6 exactly; any difference is declared in the report.

## 9. Hashes
SPEC_HASH (sha256, first 16 hex, of this file up to this heading): 9d52c927a9acd3bd
MANIFEST_HASH (manifest generator, canonical JSON): 6c312456ac894009; sha256 of the file t3fund_hs2_manifest.json as written: bd638f8f76c21e96
