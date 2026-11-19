# TODO nwsbm
See dev/DESIGN.md.

## Milestone 3
- Covariates on the ties (blockmodels `bernoulli_covariates`, `poisson_covariates`, `gaussian_covariates`; R `sbm` `covariates`): not done in 0.2.0. Reason: with a covariate the log-likelihood of a pair depends on the pair, so the block statistics (sum of y and number of cells per block pair) are no longer sufficient. The E step then needs the dense n x n x Q^2 evaluation of the pairwise terms each sweep (no sparse adjacency), and the M step a joint Newton or gradient step for block parameters and regression coefficients (blockmodels uses a numerical optimiser). Plan: dense n x n covariate matrices in the plugin (memory 8 n^2 p bytes), Newton with analytic gradient for the logit/log link, same starts; check against `sbm::estimateSimpleSBM(..., covariates = )` including ICL (the covariate coefficients add `p` parameters).
- Degree-corrected two-mode model (greed `DcLbm`) and the exact ICL for the Poisson (not degree-corrected) and Gaussian models (no greed model exists for them).
- Gaussian model with block-specific variances; ZIP / zero-inflated models of `blockmodels`.
- greed-style greedy merge/swap search on the exact criterion (`icl(exact)` evaluates the exact ICL at the VEM partition; greed's greedy search can be one or two nodes better).
- Check the plugin timing against 0.1.0 on an idle machine (the single-thread times of 0.2.0 were measured with a load average of about 6, see DESIGN).
## Smaller
- Posterior memberships and block order: option to order blocks by size.
- Linux and Windows plugin build (pthreads and `_beginthreadex` paths are in the source) and test via CI after the first push; Stata 18/19 tests.
- Menu entry for nwsbmplot in the core (nwinstall.ado, Graphics or Statistical Models; see the hand-over note of 0.2.0).
- Integration into the core (`nwinstall`, `nwversion`, website, book chapter) is a core release.
