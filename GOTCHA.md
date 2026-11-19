# GOTCHA - known pitfalls of nwsbm

- **ICL is the blockmodels/sbm one, not greed's.** `ICL = PL - 1/2 (npar log(Ndata) + (Q-1) log(n))` with PL the expected complete-data log likelihood (no entropy), npar = Q(Q+1)/2 (undirected) or Q^2 (directed), Ndata = observed unordered (undirected) or ordered pairs. It is the value `sbm::estimateSimpleSBM(...)$ICL` returns. Do not compare it with the exact ICL of `greed`.
- **blockmodels clamps the memberships** into [0.1/n, 1 - 0.1/n] (TOL_P/n) and renormalises after every E step. The plugin does the same; without it PL (and hence ICL) differs from R in the 3rd to 4th decimal.
- **R's local optima.** Away from the planted Q (over-fitted Q) R and nwsbm can end in different local optima (either one better, 0.1 percent in ICL). The tests therefore require equality of ICL only at the selected Q and "not worse than R by 0.3 percent" elsewhere.
- **`drop` of the only variable removes the observations** (Stata resets _N to 0). `_nwsbm_generate` drops old result variables first and calls `set obs` afterwards.
- **Test fixtures come from `dev/ref_sbm.R`** (R sbm 1.0.0, blockmodels 1.2.1). Re-running it changes nothing (fixed seeds) except the timing values in `time_u_n*_fit.csv` and `timing.csv`.
- **Mata `{ a; b }` on one line** is a syntax error in a do-file (`invalid expression`): write one statement per line inside braces.
- **Mata library is built with Stata 17** (`lib/build.do` refuses newer Stata).
- **Plugin protocol check:** the Mata layer asks the plugin for its protocol number (code 99); after changing the layout of the input block raise `SBM_PROTOCOL` in `native/sbm_vem.c` and `SBM_PROTOCOL_NEEDED()` in `unw_sbm.do` together. Never shown to users.
- **R sbm has no `dcsbm` and no covariates-free fixed variance.** `blockmodels` 1.2.1 / `sbm` 1.0.0 offer bernoulli, poisson, gaussian (one free common variance), ZIP and covariate versions; degree correction comes from `greed::DcSbm` (exact ICL, not comparable with the blockmodels ICL). `sigma2(#)` is an extension of nwsbm.
- **greed's printed ICL uses the Dirichlet(1) marginal of the labels** (`lgamma(K) + sum lgamma(n_k + 1) - lgamma(N + K)`), not the `-log K + sum lgamma(n_k) - lgamma(N)` of the C++ code that is used internally for the merge path. `icl(exact)` reproduces the printed value.
- **nwset puts the columns of a wide two-mode matrix into mode 1** and the rows into mode 2; `blocks(#1 #2)` follows this order (mode 1 first). R `estimateBipartiteSBM` has rows first. The test of the two-mode fixtures swaps the modes.
- **Missing ties are Stata missing values in the plugin's input block** (not -1, which is a legal Gaussian value). A tie value is "absent" (0) only for Bernoulli/Poisson sparsity; Gaussian zeros are observed zeros and stay in the counts.
- **Two-mode square matrices**: in test helpers a two-mode fixture with n1 == n2 must not be symmetrised or have its diagonal removed (`sbmt_loadv` with dir = -1).
- **Mata classes of the core are not available while building lib/lnwsbm.mlib**: pass plain vectors (`get_modes()`) from the ado instead of a pointer to the network class.
- **Machine load distorts the timing test.** Compare timings only on an idle machine (uptime load below 2).
