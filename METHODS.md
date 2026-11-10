# Methods and scope

Paper: [arXiv:2607.13009v1](https://arxiv.org/abs/2607.13009v1). The adopted spline method is described in [arXiv:2505.19083v2](https://arxiv.org/abs/2505.19083v2). The reference cosmology is from [Planck 2018, arXiv:1807.06209v4](https://arxiv.org/abs/1807.06209v4).

## A1: exact identities

With GR, spatial flatness, the Friedmann constraint, total density and pressure, and dimensionless distance D=H0*d_M/c, H/H0=1/D'. Pressure in units of rho_c0*c^2 is

`p_bar = -(2/3)*(1+z)*D''/(D')^3 - 1/(D')^2`.

The total-fluid NEC combination is

`NEC_bar = -(2/3)*(1+z)*D''/(D')^3`.

Thus negative total pressure requires `D' + (2/3)*(1+z)*D'' > 0` when D'>0. Negative pressure is an additional stated prior; it does not follow from GR and flatness alone. The matter-only, de Sitter and flat LambdaCDM pressure limits are 0, -1 and Omega_m-1. These exact identities, not every prose assertion associated with the paper's original derivation claim, constitute A1.

## N1: conditional numerical agreement

N1 asks whether a specified independent reconstruction has a negative pointwise upper 97.5% NEC percentile anywhere on its grid. The absence of such a point is consistent with the paper's no-violation-evidence statement. It does not reproduce the authors' unavailable posterior, prove the NEC, or provide a simultaneous band.

Distance residuals use six free amplitudes plus the fixed origin, at z=[0, .35, .512, .782, 1.236, 1.626, 2.33], about a flat LambdaCDM baseline Omega_m=.3153. A degree-four spline is used (order five). Closure I4 fixes residual derivatives 1 and 3 at the origin and derivative 2 at the last knot to zero. Closure I2 instead fixes derivatives 1 and 2 at the origin, and derivative 2 at the last knot. These choices are ours; neither is identified as the authors' exact closure.

DESI DR2 contributes 13 BAO measurements with the full covariance. DES-Dovekie contributes 1820 supernovae with STAT+SYS. The NPZ `cov` array is packed precision, not a covariance: it is expanded symmetrically and whitened with the transpose of its lower Cholesky factor. A constant supernova magnitude mode is projected out analytically. The luminosity-distance redshift prefactor uses zHEL, with zHD for the comoving distance. Cubic log-distance expansion and QR compression accelerate evaluation; exact likelihood spot checks bound the observed error at sampled points, not the entire prior.

The sampled calibration coordinate is q=h*r_d. A Gaussian r_d=147.09 +/- .27 Mpc and nominal flat h range [.5,1.2] motivate q bounds (73.5,176.5). Conditional r_d sampling uses the first-order 1/r_d Jacobian mean shift, .27^2/147.09. This is an approximation to the joint prior and uses hard q bounds rather than exact h truncation for each draw. N1 depends on distance derivatives rather than this sampled calibration, but the approximate marginal prior and joint likelihood can affect the distance posterior. The saved H0 and H/(1+z) diagnostics are not certified reproductions of the paper's values.

Both variants bound residual amplitudes to +/- .5 and require positive D' on 468 redshift points spanning [0,2.33], including z=.47. REIRA additionally requires negative pressure on that grid; RA does not. Off-grid violations are not ruled out. The baseline distance is integrated on 48001 points through z=2.4. A Student-t importance proposal with 8 degrees of freedom is fitted locally, with 50,000 draws per run. Direct normalized weights determine means and NEC percentiles. There is no redshift smoothing. Two seeds and two closures test limited numerical stability. A well-behaved importance diagnostic is not proof against missed modes; these checks do not establish author-pipeline fidelity.

## Claims outside the reproduced scope

The headline 3.5-sigma maximum, its location, the reported H0 and h*r_d calibration, the numerical RA-to-REIRA contrast, all robustness fits, the printed mock table, and global calibration of the selected peak are outside the reproduced scope of this package. Qualitative figure agreement or acceptable calibration of a different mock ensemble is insufficient to label those claims reproduced. Adding redshift averaging changes the target quantity and cannot substitute for the authors' sample-to-band recipe.

No major error finding is asserted. Implementation questions should first be resolved with the authors under the project's response process. The candidate is prepared while that process is pending; it has not been approved for release.

## rc2 numerical integrity correction

The lower Cholesky factor is explicitly triangular and checked against the full source matrix before whitening. This corrects a pinned local library result containing unused-triangle entries; it does not change the intended likelihood, data, priors or spline model. Full numeric outputs and NEC arrays must be finite and internally consistent.
