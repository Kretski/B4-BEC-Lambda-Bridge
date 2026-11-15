# B4.9 — Proca eigenvector tracking and mediator-mass robustness v0.1

## Objective

B4.9 strengthens B4.8 in three ways:
1. It constructs the explicit coupled 4x4 Fourier matrix for
   `(delta_phi, delta_psi, delta_A0, delta_A1)`.
2. It tracks the six frequency roots by eigenvector overlap, rather than only
   nearest-frequency assignment.
3. It varies `M_A` while holding the scalar-current background and couplings
   fixed, and checks the avoided-crossing gap and complex-frequency spectrum.

Conventions: 1+1 dimensions, metric `(-,+)`, Fourier `exp(-i omega t+i k x)`.
The sourced-background contact/seagull term is retained.

## Independent matrix check

The coupled matrix is real symmetric for real frequencies:

\[
H(\omega,k)=
\begin{pmatrix}
D_\phi&0&0&j_\phi\\
0&D_\psi&0&j_\psi\\
0&0&M_A^2+k^2&\omega k\\
j_\phi&j_\psi&\omega k&\omega^2-M_A^2
\end{pmatrix},
\]

with

\[
D_\phi=\omega^2-k^2(1+3\Lambda q_0^2)-2\omega k\bar A g_\phi,
\quad
D_\psi=\omega^2-k^2-2\omega k\bar A g_\psi,
\]

\[
\bar A=J_0/M_A^2,\quad
J_0=g_\phi w_0q_0+g_\psi v_0r_0,
\]

\[
j_\phi=g_\phi(w_0k-q_0\omega),\quad
j_\psi=g_\psi(v_0k-r_0\omega).
\]

The implementation verifies numerically that `det(H) = -P_full(omega,k)`,
where `P_full` is the B4.7/B4.8 pole-cleared degree-six polynomial.

At each root, the right null vector is recovered from the smallest right
singular vector. The scalar fraction is
`|delta_phi|^2+|delta_psi|^2` divided by the total four-component norm;
the Proca fraction is its complement.

For complex frequencies, the full complex `omega` is passed to the matrix.
It is not replaced by its real part.

## B3 illustrative avoided-crossing result

For `M_A=1`, `g_phi^2=0.05`, `g_psi=0`, `w0=q0=1`, `v0=r0=0`,
`Lambda=0`:
- local superluminal scalar speed: `c=1.0511838048`;
- uncoupled crossing estimate: `k/M_A=3.0862523`;
- eigenvector-tracked minimum gap: `0.0366612 M_A`;
- at the endpoints `k/M_A=2.7` and `3.6`, the two upper positive-frequency
  tracked modes exchange scalar/Proca composition;
- minimum mean adjacent-step eigenvector overlap: about `0.999996`;
- maximum null-vector residual in this crossing window: below `5e-13`;
- no complex frequencies above the stated `1e-6` relative threshold.

The composition exchange is the key avoided-crossing signature: the branch
that begins scalar-dominated becomes more Proca-dominated, while its partner
gains scalar character. This is an adiabatic eigenmode continuation result,
not a claim that a field species retains the same identity through the crossing.

## Fixed-background mediator-mass sweep

All quantities except `M_A` are held fixed:
`g_phi^2=0.05`, `g_psi=0`, `w0=q0=1`, `v0=r0=0`, `Lambda=0`.
The scan window is centered on the uncoupled crossing estimate for each mass.

| M_A | local c | predicted k_cross | minimum gap |
|---:|---:|---:|---:|
| 0.50 | 1.215250 | 0.724080 | 0.080342 |
| 0.75 | 1.092453 | 1.705188 | 0.050188 |
| 1.00 | 1.051184 | 3.086252 | 0.036661 |
| 1.50 | 1.022463 | 7.037399 | 0.023961 |
| 2.00 | 1.012577 | 12.570799 | 0.017843 |
| 3.00 | 1.005571 | 28.381773 | 0.011834 |

The gap decreases as `M_A` increases for this fixed-coupling, fixed-background
family. The local Proca-induced EFT coefficients scale as `g^2/M_A^2`, so
the interaction becomes weaker as the mediator is made heavier. This is a
result for this benchmark family, not a universal monotonicity theorem.

## Complex-frequency checks

The B3 illustrative and B1 boundary scans remain real within tolerance.
The deliberately selected stress case outside the B1 sufficient domain
contains complex-conjugate roots and is retained as a regression test for
the instability detector.

The B1 boundary remains defective in the local first-order reduction.
Real frequencies and small numerical null residuals do not prove strong
hyperbolicity at that boundary.

## Tests and outputs

`pytest -q`: 15 tests pass.

The runner writes:
- `results/scan_summary.json`
- `results/B3_illustrative_tracked_roots.csv`
- `results/B1_boundary_tracked_roots.csv`
- `results/complex_frequency_stress_outside_B1_domain_tracked_roots.csv`
- `results/B3_avoided_crossing_eigenvectors.csv`
- `results/mediator_mass_sweep.csv`

Run:

```bash
python -m pip install -e ".[test]"
python -m pytest -q
python run_scan.py
```

## Limitations

- Constant backgrounds and linear perturbations; classical tree-level mediator.
- B3 is the controlled illustrative point, not the original B3 scan point
  with reported maximum `|c|≈7.315`. Exact reproduction requires that point's
  full parameter tuple.
- Eigenvector overlap tracking is reliable here because the grid is dense
  and adjacent overlaps remain close to one. In a future larger parameter
  space, adaptive step control and subspace tracking at degeneracies should
  be added.
- This is not a proof of global strong hyperbolicity or UV consistency.
