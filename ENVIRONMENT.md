# Environment

Reference calculations: Python 3.14.3, NumPy 2.5.3, SciPy 1.18.1, SymPy 1.14.0, mpmath 1.3.0, macOS arm64. BLAS thread counts are set to one by the runner. A full separate-environment replay and its timing are recorded in the private delivery report. The suite is designed for a workstation; allow several minutes plus installation/download time and approximately 1 GB RAM. Cross-platform replay remains subject to the declared numerical tolerances; the supplied Linux CI has not yet executed. Installation is tested with the exact pins, not with unbounded latest versions.

The rc2 informed correction reused the previously clean pinned Python environment for reference regeneration and a manifest-only scratch replay. It did not reinstall dependencies or claim a new clean-context review. Factorization and complete-output controls were added; numerical/scientific settings and tolerances remain unchanged.
