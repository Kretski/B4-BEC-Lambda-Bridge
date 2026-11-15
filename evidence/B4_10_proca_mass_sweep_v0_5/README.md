# B4.10 v0.1 — Proca mediator mass at fixed EFT coefficients (original B3 maximum-speed point)

Separate package. B4.9 is vendored **unchanged** in `vendor/` and checked byte-for-byte
against the hand-off manifest (`tests/test_b49_untouched.py`). This release contains only
the first mandatory experiment (M_A -> infinity) plus exploratory outputs.

## Prescription

Held fixed: `p, s, kappa, Lambda` and the scalar background `(w,q,v,r)` of the recovered
B3 row (`data/B3_max_speed_7_315096.json`). For each `M_A`:

    g_phi = M_A sqrt(p),   g_psi = -M_A kappa / sqrt(p)     (p=g_phi^2/M^2, s=g_psi^2/M^2, kappa=-g_phi g_psi/M^2)

This is not the B4.9 fixed-coupling family; a negative-control test shows the two families
behave differently.

## Structural result (derived, then verified numerically and exactly)

1. The contact coefficients `Abar g_phi = p w q - kappa v r` and `Abar g_psi = -kappa w q + s v r`
   do not depend on `M_A`. Hence the local (B4.9 EFT) polynomial is `M_A`-independent and equals
   the T1.7 B3 characteristic polynomial for every `M_A`.
2. With `x = k/M_A`, `c = omega/k`, `eps = x^2/(1+x^2)`:

       Q(c; eps) = F(c) - eps c^2 d1(c) d2(c) = 0        (degree six in c)

   so the entire six-branch spectrum depends on `x` only. At fixed EFT coefficients a sweep in
   `M_A` at fixed `k` is the same one-parameter family as a sweep in `k` at fixed `M_A`.
3. First-order shift of the four scalar-like roots: `c_i(x) = c_i^B3 + A_i x^2 + O(x^4)`,
   `A_i = c_i^2 d1 d2 / F'(c_i)`.
4. Dressed Proca frequency at `k -> 0`: `omega_P = ± M_A sqrt(1 + p q^2 + s r^2) = ± M_A sqrt(det A)`
   (1.6592 M_A here, not M_A).
5. Exact statement (`b410/exact.py`, rational arithmetic on the exact binary values of the B3
   numbers): the discriminant of `Q` in `eps` has no real root in (0,1], and `Q` has six real
   roots at small `eps`; therefore all six roots are real and simple for every `x > 0`.

## Numerical results (`python run_experiment1.py`, see `results/exp1_summary.json`)

| Check | Result |
|---|---|
| p, s, kappa reconstructed from g(M), M = 1e-3 ... 1e6 | max rel. error 3.4e-16 |
| local roots vs B3 roots and vs independent T1.7 polynomial, M = 1e-3 ... 1e6 | max abs diff 5.3e-15 |
| spectrum collapse in x = k/M (M = 0.3, 1, 7) | max abs dc 2.0e-13 |
| B4.9 polynomial vs QZ (linearized pencil) vs mpmath (60 digits), 6 points x 6 roots | max abs diff 6.2e-15 |
| null-vector residual / second-smallest-to-smallest singular value | <= 2.0e-15 / >= 1.8e14 |
| leading shifts A_i (roots -4.9164, -0.1746, 2.7969, 7.3151) | 17.82, 0.00105, -7.76, -11.49 |

The largest relative root condition number (1.0e4) belongs to the slow root omega = -0.00175
at k = 0.01; it is large because the relative measure divides by |omega|.

## B4.9 limitation found (B4.9 itself left unchanged)

`full_roots` trims the leading coefficient when `|a6| <= 1e-13 max|a_i|` and then silently
returns fewer than six roots. In this family `a6 = -M^2`, `a0 ~ M^4 k^4`, so it happens
for large `M k^2`, e.g. (M=40, k=80) 5 roots, (M=1000, k=140) 4 roots, (M=1e5, k=1e3) 3 roots
(`results/exp1_b49_validity_map.csv`). `b410/spectrum.py` raises `RootCountError` in that
case and provides a scale-robust solver (the exact scaling puts everything at M=1; agreement
with mpmath over the whole map: 6.5e-13 relative).

## Domain note

Recomputed invariants: `p q^2 = 0.9106`, `s r^2 = 0.8423`, `p w^2 = 0.6692`, `s v^2 = 0.0741`.
The original B3 protocol definition of `eft_domain` was not recovered, so no domain claim is
made. `p w^2 > 0.5`, which is where T1.7.2 saw the first non-real roots in random scans; this
specific row is real-rooted (exactly, for all x, by item 5 above).

## Tests

`python -m pytest -q` (19 tests): vendored B4.9 byte-identical; B3 record vs independent T1.7
polynomial; coupling reconstruction; M-independence of the local limit; exact x-collapse;
negative control (fixed-coupling family does not collapse); x^4 residual of the leading shift;
dressed Proca frequency; README matrix vs B4.9 `coupled_matrix`; QZ and mpmath agreement;
simple-root and residual checks; exact absence of complex roots; B4.9 root-drop detection;
extreme-scale stability; overlap tracker stability.
B4.9's own 15 tests were run separately in the recipient copy and pass.

## Not established

Global strong hyperbolicity, causality, UV consistency, the physical meaning of the Proca
background at this large `J0`, and anything for other backgrounds. `exp1b_*` files are
exploratory; `scalar_fraction` uses the plain Euclidean norm of the null vector of an
indefinite matrix and is qualitative only.

---
# Experiment 2 (v0.2.0): energy-weighted mode mixing

Run: `python run_experiment2.py [outdir]` (artifacts in `results/exp2_*`), tests: `pytest` (47 tests).
B4.9 is vendored and hash-checked, unchanged.

## 1. Quadratic Hamiltonian and energy metric (`b410/energy.py`)
Fields (δφ, δψ, A0, A1), 1+1D, metric (−,+). Reconstructed quadratic Lagrangian
`L2 = ½φ_t²+½ψ_t²+½A1_t² − ½(1+3Λq²)φ_x² − ½ψ_x² + ½A0_x² − A1_t A0_x + Ābar-g_φ φ_tφ_x + Ābar-g_ψ ψ_tψ_x
 + g_φ A1(wφ_x+qφ_t) + g_ψ A1(vψ_x+rψ_t) + ½M²A0² − ½M²A1²` (A0 carries no time derivative).
Cycle average on a plane wave: `<L2> = ¼ z†𝓗(ω,k) z`. **B4.9's H equals U†𝓗U with constant U=diag(1,1,i,i)** (checked to 3e-14; u=−i fails).
A0 is non-dynamical: its row of `H z=0` slaves it to (φ,ψ,A1); the Schur complement has K00=M²+k² (constant) so there is no extra root and no extra Krein sector.
Three separate diagnostics:
* eigenvector overlap — tracking only;
* Krein value `K=z†H'(ω)z` (unitary invariant);
* canonical energy `E = <T^00> = Σ coef (n_t−1)<∂Φ∂Φ> = ¼ ω z†H' z` on shell (as <L2>=0). Verified three ways (Whitham, Noether terms, brute-force time-domain average) to ≤4e-9 (1e-13 for x≤10).
`sign E = sign ω · sign K`. Group velocity `v_g = −z†∂_kH z / z†H'z` is checked against finite differences of the roots.

## 2. Results (B3 row, fixed p,s,κ, x=k/M ∈ [1e-3,1e3])
* All six roots real and simple (QZ imag = 0; mpmath dps 60 agrees ≤4.6e-10 relative, growing with x because QZ of the unscaled pencil degrades at k/M=10³; exact discriminant from Exp. 1).
* Krein sign = sign ω on every branch for all x; **no Krein sign change; every branch has E>0**.
* Overlap tracking (adaptive log-x, assigned overlap ≥0.98) equals ω-ordering; roots cannot cross. Hence composition exchange *must* be an avoided-crossing-type exchange between positive-energy branches.
* Composition exchange (energy-share metric: Proca share of the A-group energy crossing 0.5): branches 0,1,3,5. Crossings at x≈0.14 (br5), 0.175 (br0), 0.19 (br3), 0.22 (br1) — |ω|/M at crossing 1.9, 1.8, 0.48, 0.89. Branches 2,4 do not exchange.
* Exchange does **not** coincide with a Krein change (none exists). It lies outside the operational EFT windows (|ω|/M ≤ η and x ≤ η: all four scalar modes inside only for x ≤ 0.013 (η=0.1), 0.040 (η=0.3)).

## 3. Caveats (keep in conclusions)
* `eft_domain`: reconstructed from T1.7.1 `rank_one_grid` (all four X<0.95). The B3 row satisfies it (Xq=0.911, Xr=0.842, Xw=pw²=0.669, Xv=0.074) but fails the strict 0.5 cut. The T1.7.3 scan code is not available — reconstruction, not proof. The operational windows above are an additional, separate criterion.
* g_φ,g_ψ ∝ M_A: M_A→∞ at fixed EFT coefficients is a formal limit; perturbative UV completeness is **not** established. All points are marked `unestablished`.
* Composition metrics are decompositions, not unique mixing angles. The energy-share metric is the justified one (A0 slaved, signed, sums to 1). The Euclidean amplitude share is stored only as a cross-check and is **not** equivalent (it crosses 0.5 for branches 2,4 where the energy metric does not; test documents this).
* Single background point. L2 was first reconstructed from B4.9's H; the audit (`AUDIT_E2.md`) now re-derives it from the T1.7 action independently (see there).

---
# Experiment 3 (v0.3.0): multi-background analysis — see `E3_REPORT.md`; run `python run_experiment3.py`; 68 tests.

**v0.3.1:** energy naming corrected (shift-symmetric energy = δ²(H−μQ), not lab T00); see `E3_ENERGY_NOTE.md`, `b410/charges.py`, `run_energy_triplet.py`; 70 tests.

**v0.4.0 (B4.11):** front velocity / causality — `B4_11_REPORT.md`, `b410/causality.py`, `run_b411.py`, `run_b411_population.py`; 78 tests.
