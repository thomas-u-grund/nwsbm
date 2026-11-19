# Changelog

## [0.2.0] - 2026-11-19
- Edge models: `model(poisson)` for counts and `model(gaussian)` for real-valued ties (shared variance estimated, or fixed with `sigma2(#)`), for undirected and directed one-mode networks; objective and ICL as in R `blockmodels`/`sbm` (Poisson includes the log-factorial term, Gaussian one extra parameter for the variance).
- `degreecorrected` (with `model(poisson)`): degree-corrected block model (Karrer and Newman 2011) with the degree parameters profiled out (proportional to the observed degrees).
- Two-mode networks: the latent block model, `blocks(#1 #2)` or `maxblocks(#1 #2)` / `minblocks(#1 #2)` (blocks of mode 1 and mode 2), Bernoulli, Poisson and Gaussian; results `e(Q1)`, `e(Q2)`, `e(alpha2)`, `e(blocksize2)`, `e(pi)` (block of mode 1 by block of mode 2).
- `icl(exact)`: the exact ICL of R `greed` (conjugate priors) for the Bernoulli model and the degree-corrected Poisson model of one-mode networks without missing ties; the table keeps the criterion of `blockmodels` in `e(ICLbm)`.
- `threads(#)`: the starts of one number of blocks are evaluated by independent threads (default one per core); the result does not depend on the number of threads. `fast`: at most two random starts and no merge pass.
- `estat icl` and `estat blocks` after a fit; `nwsbmplot` (blocked adjacency matrix, block graph) with dialog `nwsbmplot.dlg`; `nwsbm.dlg` has the new options.
- Not in this version: covariates (see TODO), multiplex and dynamic networks.
- Run so far (Stata 17 BE on macOS only): agreement with R sbm 1.0.0 for Poisson and Gaussian (6 fixtures, undirected, directed, missing ties) and for two-mode networks (4 fixtures), with greed 0.6.2 for the degree-corrected model and the exact ICL, recovery of planted partitions, error paths, timing. Not run: Stata 18/19, Windows, Linux.

## [0.1.0] - 2026-11-19
- First version: `nwsbm netname, blocks() maxblocks() minblocks() model(bernoulli) seed() restarts() generate() replace` fits the Bernoulli stochastic block model (undirected = blockmodels `SBM_sym`, directed = `SBM`) by variational EM in the native plugin `sbm_vem`; missing ties are masked. ICL as in R `blockmodels`/`sbm` (expected complete-data log likelihood minus the BIC-type penalties, the entropy is not added). Starts: spectral clustering of the residual matrix (subspace iteration, k-means), splits of the (Q-1) solution, merges of the (Q+1) solution, `restarts()` random starts. Results in `e()` (ICL table, connectivity matrix `e(pi)`, block proportions and sizes), memberships and probabilities by `generate()`.
- Dialog box `nwsbm.dlg` (menu entry Statistical Models in the core menu, opens with `db nwsbm`): network, blocks, maxblocks, minblocks, model, restarts, seed, generate and replace.
- Not in this version: Poisson/Gaussian edges, degree correction, latent block model (two-mode), covariates, plots, `estat`.
- Run so far (Stata 17 BE on macOS only): agreement with R sbm 1.0.0 on 8 fixtures (same Q, ARI 1.0, ICL equal to 1e-5 at the selected Q), recovery of planted partitions, error paths, missing plugin, timing against R. Not run: Stata 18/19, Windows, Linux (plugins for Linux and Windows are built by the workflow in `.github/workflows`).
