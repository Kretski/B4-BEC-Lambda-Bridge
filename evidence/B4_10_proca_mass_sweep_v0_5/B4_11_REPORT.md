# B4.11 — front velocity (causality) of the full Proca system (package v0.4.0)

Run: `python run_b411.py` (41 backgrounds incl. time-domain checks), `python run_b411_population.py` (140 000 background/mode samples).
Tests: `pytest` (78). The frozen v0.3.1 results (E3) are untouched.

## Principal symbol (exact)
At fixed M, ω=ck, k→∞, the homogeneous degree-6 part of det H(ω,k) is **−M² k⁶ Q_inf(c)** (sign fixed in v0.5: det H(ck,k)/(k⁶Q_inf)=−M²) with
Q_inf = F − c² d1 d2 = (1−c²)d1d2 + p(w−qc)²d2 + s(v−rc)²d1  (F = local B3 polynomial). Checked: det H(ck,k)/(k⁶Q_inf) = −1 = −M²
for several c at k=1e5 (mpmath, 1e-10). It is M-independent for the fixed-EFT family and equals the x→∞ limit of the E2/E3 sweep.
Front velocity of the (hyperbolic) system = max|c| over the roots of Q_inf. Complex roots of Q_inf ⇒ Im ω ∝ k ⇒ ill-posed Cauchy problem.

## Independent time-domain check
Pseudo-spectral evolution of the reduced Hamiltonian system (N=2¹⁵, L=200, Gaussian φ pulse σ=0.05, threshold 1e-7) on all 28 real-principal backgrounds:
measured front speeds right/left agree with max/min Q_inf roots to ≤0.08 %. For rank1_minus boundary #4487 (complex Q_inf, Im c = 0.166) the evolution
blows up with growth rate = Im c · k_max (to 10 %), amplitude >1e50 — the Cauchy problem is ill posed.

## Results (41 backgrounds; anchor separate)
* 28 of 29 non-control backgrounds have a real principal symbol (anchor + 27); **front velocity from 1.369 to 8.261 (median 2.40); none is subluminal (≤1).** All 12 negative controls and #4487 have complex Q_inf.
* Front velocity ≥ local EFT maximum speed in every case (ratio 1.0001–2.11): the Proca completion does **not** remove the superluminality of the local EFT; it keeps or raises the maximum speed. Anchor: local 7.315 → front 8.261.
* Max group velocity found in E3 equals the front velocity (the supremum is approached as x→∞).
* Propagation is non-reciprocal (|c_max| ≠ |c_min|), as expected with a preferred frame (anchor +8.26 / −6.53).
* Population scan (T1.7.2 generator, 10 000 backgrounds × 2 modes per cut): no subluminal case at any cut. Median front: 1.21 (cut 0.05), 1.30 (0.1), 1.47 (0.2), 2.04 (0.5), 2.92 (0.95); minimum 1.067 at cut 0.05. Real principal fraction: 100 % up to cut 0.5, 99.4/98.6 % at 0.7, 96.2/92.5 % at 0.95 (minus/plus).
  The excess (front−1) grows roughly like √(amplitude cut). Heuristic reading (not proven): at c=±1 the Proca and ψ light cones coincide; the coupling j splits the degenerate pair and the outer root moves beyond 1.

## What this does and does not establish
* Establishes (numerically, plus the exact principal symbol): in this rank-one Proca realisation the UV front velocity exceeds the Minkowski light speed for every sampled non-trivial background; it is a property of the principal symbol, not an artefact of the EFT truncation.
* Does not establish acausality in a global sense: a well-posed system with finite propagation speed 8.26 has its own characteristic cone relative to the preferred frame. Whether that is acceptable (no closed causal curves, coupling to gravity/light cone of other fields) is a separate question. Violation is relative to the Minkowski cone of the A_μ, φ, ψ kinetic terms.
* The complex-principal-symbol cases (#4487, 12 controls) are UV ill-posed — stronger than the local loss of reality found in E3.
* Pointwise backgrounds only; perturbativity (g∝M_A) is B4.12.
