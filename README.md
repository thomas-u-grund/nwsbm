# nwsbm

<!--version-->**Version 0.1.0** (released 19 November 2026)<!--/version-->

Stochastic block models for networks in Stata, fitted by variational EM with a native C plugin. Part of the nwcommands family; requires the core package
`nwcommands` (3.26.0 or later). Free software, GPL-3. Follows the R packages `blockmodels` and `sbm` (objective, ICL, model names).

```
nwsbm netname, [blocks(#) | maxblocks(#) [minblocks(#)]] [model(bernoulli) seed(#) restarts(#) generate(stub) replace]
```

Status: development, not released. Implemented: Bernoulli SBM for binary undirected and directed one-mode networks with missing ties, ICL selection of the number of
blocks, spectral / split / merge / random starts, memberships and block connectivity returned (see `nwsbm.sthlp` and `cscripts/`). The ICL, the number of blocks
and the memberships agree with R `sbm` on the test fixtures (`dev/ref_sbm.R`, `cscripts/test_nwsbm_rsbm.do`). Not yet implemented: Poisson and Gaussian
models, degree correction, two-mode latent block model, covariates. See `TODO.md` and `dev/DESIGN.md`.
