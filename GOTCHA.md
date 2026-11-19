# GOTCHA - known pitfalls of nwsbm

- **ICL is the blockmodels/sbm one, not greed's.** `ICL = PL - 1/2 (npar log(Ndata) + (Q-1) log(n))` with PL the expected complete-data log likelihood (no entropy), npar = Q(Q+1)/2 (undirected) or Q^2 (directed), Ndata = observed unordered (undirected) or ordered pairs. It is the value `sbm::estimateSimpleSBM(...)$ICL` returns. Do not compare it with the exact ICL of `greed`.
- **blockmodels clamps the memberships** into [0.1/n, 1 - 0.1/n] (TOL_P/n) and renormalises after every E step. The plugin does the same; without it PL (and hence ICL) differs from R in the 3rd to 4th decimal.
- **R's local optima.** Away from the planted Q (over-fitted Q) R and nwsbm can end in different local optima (either one better, 0.1 percent in ICL). The tests therefore require equality of ICL only at the selected Q and "not worse than R by 0.3 percent" elsewhere.
- **`drop` of the only variable removes the observations** (Stata resets _N to 0). `_nwsbm_generate` drops old result variables first and calls `set obs` afterwards.
- **Test fixtures come from `dev/ref_sbm.R`** (R sbm 1.0.0, blockmodels 1.2.1). Re-running it changes nothing (fixed seeds) except the timing values in `time_u_n*_fit.csv` and `timing.csv`.
- **Mata `{ a; b }` on one line** is a syntax error in a do-file (`invalid expression`): write one statement per line inside braces.
- **Mata library is built with Stata 17** (`lib/build.do` refuses newer Stata).
- **Plugin protocol check:** the Mata layer asks the plugin for its protocol number (code 99); after changing the layout of the input block raise `SBM_PROTOCOL` in `native/sbm_vem.c` and `SBM_PROTOCOL_NEEDED()` in `unw_sbm.do` together. Never shown to users.
