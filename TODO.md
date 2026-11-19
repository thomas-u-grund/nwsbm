# TODO nwsbm
See dev/DESIGN.md.

## Open
- Zero-inflated models (ZIP) of `blockmodels`.
- Dynamic / longitudinal block model (design note, not to implement yet): networks observed at T time points with a hidden Markov chain on the block of every actor (transition matrix between blocks, or a stochastic block model with block-pair parameters per period); variational EM with a mean-field posterior over the actor paths (forward-backward on the chain for the E step, block parameters per period or shared); the ICL gets the transition parameters Q(Q-1); the input would be a list of network names (one per period) with the same actors; R references: `dynsbm` (Matias and Miele), `dynnet`. Needs a design decision on constant versus period-specific block parameters and on how nwcommands represents the periods.
- Exact ICL for the Poisson (not degree-corrected), Gaussian and covariate models (no greed model exists for them); `greedy` for networks with missing ties.
- Standard errors that account for the uncertainty of the memberships (bootstrap); the current ones are conditional on the memberships.
- Integration into the core (`nwinstall`, `nwversion`, website, book chapter) is a core release.
