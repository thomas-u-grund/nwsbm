# nwsbm

<!--version-->**Version 1.0.2** (released 20 November 2026)<!--/version-->

Stochastic block models for networks in Stata, fitted by variational EM with a native C plugin. Part of the nwcommands family; requires the core package
`nwcommands` (3.26.0 or later). Free software, GPL-3. Follows the R packages `blockmodels` and `sbm` (objective, ICL, model names).

```
nwsbm netname, [blocks(#) | maxblocks(#) [minblocks(#)]] [model(bernoulli|poisson|gaussian) degreecorrected sigma2(#) blockvar icl(exact) greedy orderblocks
                covariates(netlist) absdiff(varlist) same(varlist) signeddiff(varlist)
                seed(#) restarts(#) threads(#) fast generate(stub) replace]
nwsbm twomodenet, blocks(#1 #2) | maxblocks(#1 #2) [minblocks(#1 #2)] [model(...)]     (latent block model)
estat icl | estat blocks
nwsbmplot [, type(matrix|blocks) shade nolines missing nolabels]
```

Status: development, not released. Implemented: Bernoulli, Poisson and Gaussian stochastic block models for undirected and directed one-mode networks (missing ties allowed),
the degree-corrected Poisson model, the latent block model for two-mode networks, ICL selection of the number of blocks (also the exact ICL of R `greed`), spectral / split /
merge / random starts evaluated by threads, `estat` output, plots (`nwsbmplot`: blocked adjacency matrix, block graph) and dialogs (`db nwsbm`, `db nwsbmplot`). The ICL, the number
of blocks and the memberships agree with R `sbm` / `blockmodels` (Bernoulli, Poisson, Gaussian, two-mode) and R `greed` (degree-corrected model, exact ICL) on the test fixtures
(`dev/ref_sbm.R`, `dev/ref_valued.R`, `dev/ref_dc.R`, `cscripts/`). Dyadic covariates (networks and node-variable differences) for all three edge models, one-mode and two-mode, agree with R `sbm` on the fixtures of `dev/ref_cov.R`. The degree-corrected two-mode model and its exact ICL agree with R `greed` (`dev/ref_dclbm.R`); the covariate coefficients come with standard errors (conditional on the memberships, which the R packages do not give); `greedy` improves the partitions on the exact ICL as `greed` does; `blockvar` fits a variance for every pair of blocks. Not yet implemented: multiplex and dynamic networks, zero-inflated models. See `nwsbm.sthlp`, `TODO.md` and `dev/DESIGN.md`.

## Manual

A user manual (PDF in `manual/`, with the LaTeX source and the do-files of its examples) is archived at Zenodo: <https://doi.org/10.5281/zenodo.24398296> (all versions). Cite the manual by that address together with the package version printed by `nwsbm version`.
