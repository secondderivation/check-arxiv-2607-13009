# Release validation protocol — 2026-10-09

This is a post-hoc package protocol, written before the new release runs. It does not replace the original review gates, count as a blind review, or retroactively fix an unsuccessful reproduction.

Scope: A1 checks the exact GR / spatially flat FLRW identities. N1 is a new, implementation-conditional test of the NEC statement, using the reconstruction defined in METHODS.md. N1 is not a replacement for the original C7 author-pipeline claim. The original headline significance and the printed mock table remain unresolved and are excluded from reproduced claims.

Run quartic closures I4 and I2, with and without the negative-pressure prior (REIRA and RA), using seeds 731009 and 20261009, 50,000 Student-t proposals per configuration. Both variants impose positive distance derivative. Save ESS, maximum weight, exact-likelihood spot checks, H/(1+z) at z=0.47, and the pointwise NEC quantiles. The claimed agreement is absence of a negative 97.5th percentile at any tested grid point; this is not evidence that NEC holds universally, a simultaneous confidence band, or an author-chain reproduction.

Engineering acceptance, fixed before these runs: ESS >= 1000, maximum normalized weight <= 0.01, finite reported numbers, maximum absolute exact-minus-compressed SN chi-square difference <= 0.01 in 40 posterior-weighted spot checks. Mean and SD at z=0.47 must agree between seeds within 0.03 km/s/Mpc. Each configuration must satisfy the stated NEC criterion; failure remains visible and blocks its N1 inclusion. Replay acceptance against the saved same-seed scalar references: mean/SD tolerance 0.002 km/s/Mpc, minimum upper NEC percentile 0.02, ESS 5% relative and max-weight 10% relative; diagnostics also pass their absolute gates. These are numerical portability tolerances, not paper-reproduction tolerances.

Controls: exact matter-only/de Sitter/LambdaCDM limits; positive-pressure and negative-NEC toy cases; covariance-vs-precision whitening and constant-magnitude marginalization; damaged data and reference files; intentionally shifted replay values; failed diagnostics; private or unexpected payload files. All four input data hashes must pass. Archive contents, manifest and internal links are checked. Reference outputs are immutable during replay.

One importance proposal and two seeds can miss remote modes. Acceptance is a bounded numerical diagnostic, not proof of posterior completeness. The h-r_d prior is implemented approximately after integrating q=h*r_d; finite redshift-grid constraints and unconfirmed endpoint conventions remain limitations. No headline-significance calculation is certified here.

## rc2 implementation correction — 2026-10-28

The scientific and replay tolerances above remain unchanged. The actual input matrix must be finite and symmetric; a triangular Cholesky factor must reconstruct it (rtol=atol=1e-10; symmetry rtol=atol=1e-12). These are numerical factorization invariants. Every reported numeric field and full NEC grid/quantile array must be finite, shape-consistent and ordered, with summaries exactly matching the curves. Invalid outputs fail before persistence or acceptance.
