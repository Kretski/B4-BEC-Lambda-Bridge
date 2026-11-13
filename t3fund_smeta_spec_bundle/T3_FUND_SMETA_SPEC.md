# T3-fund specification: Scalar mediator with dynamical preferred-frame field η (SM+η)
Registered under UV_ACCEPTANCE_PROTOCOL_ADDENDUM_T3_T6.md (+ Clarification 1). All criteria, grids, tolerances and outcome classes of that addendum apply unchanged; this document only fixes the model, populations, closed forms, predictions and the C4 setup *before* any evaluation. Previous results for the static scalar mediator and for HS/HS2 are not carried over: C1–C4 are evaluated afresh for this model.

## 1. Model and conventions (M=1, signature +−)
Fields φ, ψ (light), σ (mediator, mass M), η (dynamical boost parameter of the preferred frame). With ab(f_t,f_x;η) = (cosh η f_t + sinh η f_x)(sinh η f_t + cosh η f_x) = ½ sinh2η (f_t²+f_x²) + cosh2η f_t f_x:
 L = ½(∂φ)² + ½(∂ψ)² + ½(∂σ)² − ½M²σ² + ½ f (∂η)² + σ [g_φ ab(φ_t,φ_x;η) + g_ψ ab(ψ_t,ψ_x;η)],  (∂Φ)² ≡ Φ_t² − Φ_x².
Sign convention: with signature (+−) and the same convention as the φ, ψ, σ kinetic terms, ½f(η_t² − η_x²) has the correct (positive-energy) sign for f>0; f<0 is a ghost η; f=0 is a degenerate kinetic term. (Under (−+++) the same statement reads ½f(∂η)²→−½f(∂η)²; the sign test used everywhere below is the sign of the coefficient of η_t², i.e. of the time-derivative Hessian ∂²L/∂η_t² = f, not the symbol in front of (∂η)².)
Definitions: α=(w,v), β=(q,r) the background gradients (φ̄=wt+qx, ψ̄=vt+rx); A = g_φ w q + g_ψ v r;  B = ½[g_φ(w²+q²) + g_ψ(v²+r²)].

## 2. C0 exact background and closed forms
Constant jets with constant σ̄, η̄. The φ, ψ equations are identically satisfied. The σ equation: M²σ̄ = J(η̄) with J(η)=A cosh2η + B sinh2η; the η equation: σ̄ J′(η̄)=0, J′ = 2(B cosh2η + A sinh2η).
 Branch S (σ̄≠0): tanh2η̄ = −B/A, requires |B|<|A|; σ̄ = sign(A)√(A²−B²)/M², hence σ̄ sinh2η̄ = −B, σ̄ cosh2η̄ = A (exact).
 Branch Z (σ̄=0): J(η̄)=0, tanh2η̄ = −A/B, requires |A|<|B|; here Z_t=X=1, Y=0 below.
 Boundary |A|=|B| has no stationary branch and is excluded. Note |A|>|B| is impossible for g_φ g_ψ>0 (|ab-coupling| ≤ ½Σ g(a²+b²)), so branch S needs g_φ g_ψ<0; this is why the coupling sets below have mixed signs.
Branch-S closed forms (verified, §8): with Z_t,i = 1 − g_i B, X_i = 1 + g_i B, Y_i = g_i A (i=φ,ψ),
 Hessian blocks of the light fields on the (t,x) jets: [[Z_t,i, Y_i],[Y_i, −X_i]];  Y_i² + Z_t,i X_i = 1 + g_i²σ̄² ≡ R_i² > 0;
 η–η value Hessian 4σ̄²; σ–η value Hessian 0; N_i ≡ ∂²L/∂φ_{i,t}∂η = 2g_i(Aα_i − Bβ_i); m_i ≡ ∂²L/∂σ∂φ_{i,t} = g_i(Aβ_i − Bα_i)/σ̄.
 C0 acceptance: residuals of all four Euler–Lagrange equations on the background ≤1e-14 (double) and ≤1e-40 (50 digits).

## 3. C1/C2 analytic
Linearised operator K(ω,k) (4×4, fields φ,ψ,σ,η; plane wave e^{i(kx−ωt)}, c=ω/k): K_φφ-type entries contain the principal part k²·diag(Z_t,i c² − 2Y_i c − X_i, c²−1, f(c²−1)) and lower-order (order k and k⁰) couplings between different fields only. Therefore (principal symbol, analytic):
 (i) characteristic speeds of the light-field blocks: c_{i,±} = (Y_i ± R_i)/Z_t,i, real for every background with Z_t,i>0 because R_i² = 1+g_i²σ̄² > 0 (independent of the sign of X_i); σ and η blocks are luminal (c=±1 for f≠0).
 (ii) Kinetic classification (time-derivative Hessian values κ = (Z_t,φ, Z_t,ψ, 1 (σ), f), tolerance τ=1e-12 absolute; all κ are O(1) in the manifest). Priority order:
   1. any κ < −τ → FAIL-GHOST (negative kinetic eigenvalue; takes priority even if another κ is zero);
   2. else any |κ| ≤ τ → DEGENERATE (kinetic matrix singular; not a ghost, not counted in PASS/FAIL statistics);
   3. else all κ > τ → the kinetic test passes.
   The leading coefficient of det K in c⁸ is Z_t,φ Z_t,ψ f (verified, §8).
 (iii) X_i>0 (gradient-energy positivity of the fluctuation Hamiltonian ½Z_t φ_t² + ½X φ_x² per light field) is recorded as a diagnostic flag; X_i≤0 with Z_t,i>0 is NOT classified as ghost (population S2x, evaluated by C3 only, reported separately).
 C1 (well-posedness): for Z_t,i>0 and f>0 the principal part is block diagonal with real speeds and a constant block symmetrizer built from characteristic variables, lower-order terms are constant; hence γ(k) ≤ γ_max for all k with γ_max computed per background from the first-order reduction (U = (Φ,Φ_t,Φ_x) for all fields) in characteristic variables. C1 consistency checks: eigenvalues of the first-order reduction match the roots of det K to relative 1e-9 at three k; max_k γ(k) over the C3 grid ≤ γ_max(1+1e-9). Violation → INCONCLUSIVE (implementation inconsistency), not FAIL.
 Dispersion polynomial: det K(ck,k) is a degree-8 real-coefficient polynomial in c (Hermitian K ⇒ real); built by exact polynomial arithmetic in arbitrary precision (Leibniz expansion with polynomial entries in c), not by sampling.

## 4. Populations (manifest `t3fund_smeta_manifest.json`, MANIFEST_HASH in §9)
Coupling sets (g_φ,g_ψ): G1=(0.3,−0.2) seed 41; G2=(0.25,−0.35) seed 42; G3=(0.8,−0.6) (strong coupling) seed 43. Per set 200 draws (w,q,v,r)~U[−1.5,1.5]⁴ (numpy default_rng(seed)).
 S1: branch S, Z_t,i>0 and X_i>0 (main population): G1 138, G2 133, G3 117; each evaluated for f ∈ {0.1, 1, 10} (1164 evaluations).
 S2g: branch S with some Z_t,i ≤ 0: G3 2 (f=1). Expected FAIL-GHOST if the offending Z_t,i < −τ (DEGENERATE if |Z_t,i| ≤ τ, by the §3(ii) priority rule).
 S2x: branch S, Z_t,i>0, some X_i ≤ 0: empty in this manifest (reported as empty).
 Z0: branch Z (σ̄=0): G1 62, G2 67, G3 81 (f=1), evaluated by C3 numerics, no closed-form prediction, reported separately and never mixed with S1.
 H3 (negative control, ghost detection): first 40 S1 draws in manifest order with f=−1 (separate group). Expected FAIL-GHOST 40/40 (f=−1<−τ).
 H3b (boundary): first 10 S1 draws with f=0 → class DEGENERATE (kinetic matrix singular, §3(ii)); diagnostic only, not counted.
 Anchors: A_G1, A_G2, A_G3 at (w,q,v,r)=(1,0.3,0.8,−0.2), f=1 (all branch S): implementation check, no prediction beyond the closed forms.
The sign of f in H3 is the sign of ∂²L/∂η_t² (§1); the negative control is not a physical candidate and is not mixed into the S1 evaluation.

## 5. C3 (mode stability) and registered predictions
C3 exactly as registered (k grid 61 points k/M∈[1e-3,1e3] + 31 in [1e3,1e6], P2/P3 classes with 1e-7/1e-5 thresholds, 50-digit confirmation r>1e-12, P2a k=0 check, P6 three-band rule with u(k)=|γ50−γ100|, margin 10, ≥3 consecutive refined points). For this model the k=0 polynomial is, with x=ω²: det K(ω,0)=ω⁴ Z_t,φ Z_t,ψ [(x−M²−S1)(f x + 4σ̄² − S2) − S3²], S1=Σm_i²/Z_t,i, S2=ΣN_i²/Z_t,i, S3=Σm_iN_i/Z_t,i (branch S, Z_t>0); the zero-frequency factor ω⁴ is removed before classifying; k=0 class is real iff both roots in x are real and ≥0. A background with complex k=0 roots is FAIL-UNSTABLE (as in the registered rule for P2a).
Predictions (falsifiable; any counterexample is reported without changing a criterion):
 Pred-1 (ghost): the classification follows the §3(ii) priority rule (FAIL-GHOST ⇔ some κ<−τ; DEGENERATE ⇔ none <−τ and some |κ|≤τ) for all evaluated cases; H3 40/40 FAIL-GHOST, S2g 2/2 FAIL-GHOST, H3b 10/10 DEGENERATE (f=0 exactly, other κ>τ); the sign of the numerical c⁸ coefficient equals the sign of Z_t,φ Z_t,ψ f in every non-degenerate case.
 Pred-2 (k=0): the numerical k=0 class (and roots in x to 1e-9 relative) equals the closed form of §5 for every S1 case. (The manifest already stores the closed-form class at f=1; before registration it was seen that G1 and G2 S1 draws are all closed-form-complex and G3 has 85 complex / 32 real. This pre-registration inspection is stated here; Pred-2 tests the agreement of the numerics with the closed form, not that fact.)
 Pred-3 (UV principal speeds): at k=1e6, for S1 cases whose eight roots are real, the sorted real speeds equal the sorted set {±1 (σ), ±1 (η), c_{i,±}} within 1e-4.
 Pred-4 (no UV/intermediate-k instability beyond the IR): every S1 case with real closed-form k=0 class has C3 PASS (no complex root anywhere in k/M∈[1e-3,1e6]). A counterexample would be a non-IR instability.
 Diagnostic (not part of any outcome): front speed max|c| per S1 background; the manifest stores the closed-form value (all 388 S1 backgrounds have front >1+1e-6). Superluminal fronts are allowed under preferred-foliation causality (T3-fund); this was an earlier exploratory hypothesis and is only reported.

## 6. C4 (time domain)
Designated backgrounds (manifest `C4_designated`, rule: the first two S1 draws in manifest order whose closed-form k=0 class at f=1 is real; both fall in G3, f=1). C4 is run on a designated background only if it is C3 PASS; otherwise it is reported NOT-RUN and no replacement is chosen (so C4 for this candidate is then not available).
Method as registered C4: periodic box 1024, pseudo-spectral D, RK4, T=1000, records every 5, base resolution (16384, dt 0.025), convergence (32768, 0.0125). Full nonlinear Euler–Lagrange system for (φ,ψ,σ,η), solved for second time derivatives; the coefficient matrix of second time derivatives is diag(1+σ g_φ sinh2η, 1+σ g_ψ sinh2η, 1, f) (+ off-diagonal terms none), monitored: if min over x,t of the first two entries ≤1e-3 the run is INCONCLUSIVE (Z-degenerate). The implementation is checked against the independent sympy Euler–Lagrange operator at random states (≤1e-12) before any run.
Energy (exactly conserved): E = ∫dx{ Σ_i [½Z_t,i φ_{i,t}² + ½X_i φ_{i,x}²] + ½σ_t² + ½σ_x² + ½M²σ² + ½f(η_t²+η_x²) }, Z_t,i = 1+σ g_i sinh2η, X_i = 1−σ g_i sinh2η (fields σ, η their full values); E_fl = E − E_bg.
Perturbation: all four fields, Gaussian s=15, k0p=0.5, centre 0, amplitude ε=1e-3·max(|q|,|r|); δφ right-moving (δφ_t=−Dδφ), δψ left-moving (δψ_t=+Dδψ), δσ and δη with zero time derivative. Norm 𝒩=√(dx Σ(δφ²+δψ²+δσ²+δη²)); criterion max 𝒩/𝒩(0) ≤ 3. Free control: g_φ=g_ψ=0 (σ̄=0, η̄=0), same gradients, f, perturbation, resolutions.
Outcome rule: as registered C4 incl. Clarification 1 (PASS: ratio ≤3 at both resolutions, total energy drift ≤1e-9 reached, control reaches target, 1% agreement; INCONCLUSIVE: control misses target / 1% agreement fails / drift missed with control passing / Z-degenerate; FAIL-UNSTABLE: ratio >3 with agreement; never FAIL on energy alone). Stated limitation: total-energy drift is dominated by the background energy; the norm criterion is substantive.

## 7. T7 (operator matching) — not part of this document
T7 for SM+η (integrating out σ at η̄ and matching to the B4 operators) is a separate stage to be registered in its own document before it is run, and only if at least one population reaches C3/C4 PASS. Existence of a match is not assumed.

## 8. Scope, verification and files
Tree level, this Lagrangian, sampled backgrounds only; PASS in C3/C4 does not establish a UV completion of B4 or T7. Backgrounds are new random draws (not the earlier B4-derived backgrounds); the earlier exploratory scalar_mediator_dynamicU_check.py is historical and is not reused.
Verified (t3fund_smeta_check.py, output smeta_check_output.txt): EL residual on background ≤7e-17; det K from the jet Hessian equals det of the independently linearised Euler–Lagrange operator (explicit plane waves) to 1e-14; closed forms of §2 to 3e-16; k=0 identity of §5 to 3e-15; principal part and leading coefficient as in §3 (the principal-part difference at k=1e5 is 2.8e-5, i.e. O(1/k) as expected; leading coefficient 2e-14).
Files: T3_FUND_SMETA_SPEC.md, t3fund_smeta_check.py, smeta_check_output.txt, t3fund_smeta_manifest.py, t3fund_smeta_manifest.json. Not touched: tests/test_dimensions.py, theory/operator_matching.py. No simulation or evaluation is run before this document and its hashes are registered by commit.

## 9. Hashes
SPEC_HASH (sha256, first 16 hex, of all bytes of this file before the line '## 9. Hashes'): b2fbbf16e6741ac7
MANIFEST_HASH (sha256, first 16 hex, of canonical JSON json.dumps(sort_keys=True, separators=(',',':'))): 9e9eb854437015e0
