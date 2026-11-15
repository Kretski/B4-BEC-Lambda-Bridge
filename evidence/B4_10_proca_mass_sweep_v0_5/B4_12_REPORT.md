# B4.12 — UV perturbativity / control parameters (v0.5)

Proca action NOT modified. B4.9 untouched (hash-checked). B4.11 frozen (sign fixed: det H|deg6 = −M²k⁶Q_inf).

## Exact identity (M-independent)
With the fixed-EFT prescription g_φ=M√p, |g_ψ|=M√s: **g · (background gradient/M) = √X** for each of q,w (φ) and r,v (ψ).
Raising M_A at fixed EFT therefore raises g ∝ M_A and lowers ρ=gradient/M in proportion; only the product is fixed by the EFT data.
Verified to 2e-16 on all 41 backgrounds (`results_b412/`).
Strong-coupling scale estimate Λ_* ~ 1/√p (dimensional analysis ONLY, not a computed unitarity bound); background gradient = √X · Λ_*.

## Pre-registered control sets (set before the tables were produced)
A: g≤1, ρ≤0.3 (√X≤0.30); B: g≤1, ρ≤0.1 (√X≤0.10); C: g≤√(4π) (α≤1), ρ≤0.3 (√X≤1.06). C was added as the most permissive edge, and α=1 is at the border of perturbativity, not comfortably inside it.

## Results (41 backgrounds incl. B3 anchor)
- √X_max range 0.554 … 0.974 (B3 anchor 0.954). Set A: 0/41, set B: 0/41, set C: 41/41 (29/29 non-controls).
- Minimum g required for ρ≤0.3: 1.85 … 3.25 (anchor 3.18, α≈0.81). I.e. every background in the E3 manifest is outside g≤1.
- Fluctuation-EFT window: scalar speeds deviate 10 % from local ones at x=k/M ≥ 0.19 (anchor 0.20).
- The 41 backgrounds are strong-X by construction (selected from the cut-0.95 population), so this is a statement about THEM.

## Population (10 000 backgrounds per cell, seed 172026)
Principal-symbol front velocity, 0 subluminal in every cell. Median/min of (c_front−1)/√X:
| cut on max X | front min | (c−1)/√X median | min |
|---|---|---|---|
| 0.01 (√X≤0.1) | 1.030 | 1.01 | 0.56 |
| 0.04 | 1.061 | 1.03–1.06 | 0.54 |
| 0.09 (√X≤0.3, set A) | 1.089 | 1.05–1.14 | 0.53 |
| 0.25 | 1.14 | 1.27–1.37 | 0.50 |
| 0.95 | 1.25 | 2.2 | 0.48 |
(rank1_minus/plus; complex principal symbols appear only at cut 0.95: ~4–8 %.)

## Conclusion (limited to the scanned populations)
In the controlled domain (set A, √X≤0.3) the principal-symbol front velocity is still >1 in all 20 000 samples, with excess ≳0.5√X. Superluminality does not disappear as the theory becomes more perturbative; it shrinks only together with the EFT effect itself (c−1 ∝ √X ∝ g ρ). Hence in this Proca realisation "perturbative + hyperbolic + Minkowski-causal" is not simultaneously satisfied on the scanned backgrounds; the effect is ≳ 3 % at √X=0.1 and ≳ 9 % at √X=0.3.

## Caveats
Λ_* is a dimensional estimate; backgrounds are pointwise (no check of dynamical sourcing); front velocity ≠ proof of acausality (preferred-frame cone); the ≳0.5√X bound is an empirical minimum over samples, not a theorem; 1+1D.
