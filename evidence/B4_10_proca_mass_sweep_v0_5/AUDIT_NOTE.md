# Audit note — B4.10 v0.1

* B4.9 files are byte-identical to `SHA256SUMS.txt` of the hand-off (16 files). The hand-off
  states that two corrections were made to `eigenmodes.py` after the original archive
  (`B4_9_ORIGINAL_ARCHIVE_SHA256.txt`); the original archive was not available, so that
  statement is taken from the hand-off and not independently verified.
* The README matrix `H(omega,k)` was written again independently (`reference.H_blocks`) and
  agrees with B4.9 `coupled_matrix` to 1e-10 on random complex omega.
* New B4.9 defect found, not fixed here: silent loss of roots in `full_roots` (see README).
* No result in this package is a proof of strong hyperbolicity or causality.
* The exploratory tracked-composition CSV should not be quoted as a physical mixing result.

## v0.2.0 (Experiment 2) additions
* B4.9 untouched (hash manifest verified on every run).
* `eft_domain` reconstructed from T1.7.1 text (cut 0.95 for all four X); original T1.7.3 scan generator not available. The cut-0.95 sampling of T1.7.2 `generate_samples(seed 172026)` does NOT reproduce row 8899, so provenance of that row remains unverified.
* QZ vs mpmath: 6e-16 at x=1e-3..1, 4e-14 at x=10, 6e-12 at 100, 6e-11 at 1e3 (relative to largest root). mpmath is the more reliable value at large x.
* Max relative root condition (5e7) belongs to the near-zero scalar root (relative measure divides by |ω|).
