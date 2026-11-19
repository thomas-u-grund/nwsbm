# TODO nwsbm
See dev/DESIGN.md.

## Milestone 2
- Poisson and Gaussian SBM (valued networks): M step and E step per `blockmodels` zpoisson / zgaussian; the engine in `native/sbm_vem.c` is written around the two sufficient statistics (A.Z and the pair counts), so the new models mainly change `m_step` and the E-step log terms and the ICL parameter counts (Gaussian: Q^2 (or Q(Q+1)/2) means plus one variance).
- Degree-corrected SBM (Karrer and Newman 2011; greed DcSbm) and greed's exact ICL with conjugate priors.
- Latent block model for two-mode networks (`blocks(#1 #2)`, blockmodels LBM).
- Compare recovery with `greed` on the same fixtures.
- Covariates (blockmodels `*_covariates`).
## Smaller
- Posterior memberships and block order: option to order blocks by size.
- Backward/forward exploration is a fixed pass (all Q fitted); R's early stopping is not needed but a `fast` option (fewer starts) may help for n above 3000.
- Threads (pthreads) for the restarts; not done for portability.
- Linux and Windows plugin build and test via CI after the first push; Stata 18/19 tests.
- Integration into the core (`nwinstall`, `nwversion`, website, book chapter) is a core release.
