# TODO nwsbm
See dev/DESIGN.md.

## Milestone 4
- Standard errors of the covariate coefficients (conditional on the memberships from the final Newton Hessian, or a bootstrap); R does not give them either.
- Faster Bernoulli covariate E and M steps (a tabulated softplus instead of one log1p per pair and block pair; caching the pair terms between sweeps); sparse treatment is not possible with covariates.
- Dynamic / longitudinal block model (design note, not to implement yet): networks observed at T time points with a hidden Markov chain on the block of every actor (transition matrix between blocks, or a stochastic block model with block-pair parameters per period); variational EM with a mean-field posterior over the actor paths (forward-backward on the chain for the E step, block parameters per period or shared); the ICL gets the transition parameters Q(Q-1); the input would be a list of network names (one per period) with the same actors; R references: `dynsbm` (Matias and Miele), `dynnet`. Needs a design decision on constant versus period-specific block parameters and on how nwcommands represents the periods.
- Degree-corrected two-mode model (greed `DcLbm`) and the exact ICL for the Poisson (not degree-corrected) and Gaussian models (no greed model exists for them).
- Gaussian model with block-specific variances; ZIP / zero-inflated models of `blockmodels`.
- greed-style greedy merge/swap search on the exact criterion (`icl(exact)` evaluates the exact ICL at the VEM partition; greed's greedy search can be one or two nodes better).
## Smaller
- Posterior memberships and block order: option to order blocks by size.
- Linux and Windows plugin build (pthreads and `_beginthreadex` paths are in the source) and test via CI after the first push; Stata 18/19 tests.
- Menu entry for nwsbmplot in the core (nwinstall.ado, Graphics or Statistical Models; see the hand-over note of 0.2.0).
- Integration into the core (`nwinstall`, `nwversion`, website, book chapter) is a core release.
