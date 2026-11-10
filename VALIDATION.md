# Validation

The acceptance rules are in [PROTOCOL.md](PROTOCOL.md). They are post-hoc engineering checks and conditional scientific diagnostics, not the original paper-reproduction tolerances.

`python3 reproduce.py --all` verifies every manifest-listed file, rejects unexpected payload files, fetches and hashes the four pinned data inputs, runs the falsification controls, derives the exact identities, and executes all eight seeded numerical configurations. Each result is compared with immutable reference JSON files. It checks diagnostics and seed stability, then verifies the package again. A failed condition causes a nonzero exit. Outputs go into a timestamped `_runs/` directory. `_runs/`, `_verify/`, `.venv/`, `.git/` and Python bytecode caches are excluded from payload checks; the distributor must export only the manifest allowlist.

`python3 -m unittest test_validation.py` independently exercises damaged files, invalid numerical results, shifted outputs, covariance/precision orientation and analytical nuisance-magnitude projection. `checks/algebra.py` includes physically distinct exact limits and opposite-sign controls.

The reference table in README is generated from the same JSON files used for replay. It is a compact representation of the saved results, not separate evidence. All NEC curves are provided as numerical JSON arrays. No paper figures or upstream data are redistributed.

This candidate has a local clean-environment replay. That is not a fresh-context scientific review. Remote CI has not run because the candidate has not been pushed. The included GitHub workflow runs tests only and has read-only permissions. A full independent P2 review and an explicit owner approval for the exact archive, version and destination are separate release gates. No script in this package sends mail or publishes anything.

The rc2 correction adds real-matrix Cholesky reconstruction checks, complete numeric finiteness checks and full NEC profile/summary validation. Regression controls exercise unused-triangle contamination, wrong factorization, nonfinite/malformed profiles and auxiliary diagnostics. The rc1 references are historical and must not be used to judge rc2 replay.
