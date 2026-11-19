# nwsbm design (internal; not shipped)

Owner: Thomas U. Grund (sole author). GPL-3. Part of the nwcommands family; depends only on the core. Native C plugin for all heavy computation (no Mata fallback). Mata library built with Stata 17.

## References (match as closely as possible; document deviations)
- R `blockmodels` 1.2.1 (Leger) and `sbm` 1.0.0 (Chiquet, Donnet, Barbillon): variational EM, ICL, model names bernoulli / poisson / gaussian, membership types SBM (directed), SBM_sym (undirected), LBM (two-mode). `sbm::estimateSimpleSBM` is a wrapper of blockmodels, so the two share every number (checked in sources: `ICL <- PL - .5*(n_parameters*log(data_number) + (Q-1)*log(n))`).
- R `greed` 0.6.2 (Come, Randriamanamihaga, Ouattara, Rey): exact ICL with conjugate priors and greedy merge/swap search; degree-corrected SBM. A later milestone.

## Models (stages)
| Milestone | Model | Status |
|---|---|---|
| 0.1.0 | Bernoulli SBM, undirected (blockmodels `SBM_sym`) and directed (`SBM`) one-mode networks, missing ties allowed | this version |
| 0.2.0 | Poisson and Gaussian SBM (valued networks), same engine | TODO |
| later | degree-corrected SBM (greed's DcSbm, Karrer-Newman), latent block model (two-mode, `blocks(#1 #2)`), covariates | TODO |

Out of scope for 0.1.0: valued edges (error with a pointer to dichotomising), two-mode networks, degree correction, covariates, multiplex/dynamic, self-loops (ignored, as blockmodels), exact ICL of greed, plots, `estat`, parallel threads.

## Model and objective (Bernoulli, one block structure over ordered pairs)
Latent labels Z_i ~ Multinomial(alpha); y_ij | Z_i=q, Z_j=l ~ Bernoulli(pi_ql); undirected: y_ij = y_ji and pi symmetric, pairs i<j only. Observed pairs only (missing ties are masked, as blockmodels: they contribute neither to the likelihood nor to the counts).
Variational EM with the mean-field family q(Z) = prod_i Mult(tau_i):
- M step: alpha_q = mean_i tau_iq; pi_ql = sum_{i!=j} tau_iq tau_jl y_ij / sum_{i!=j} tau_iq tau_jl m_ij (m = observed mask). Undirected networks use the symmetric full matrix and the objective is halved (as blockmodels).
- E step (fixed point, at most 10 sweeps, stop at max change < 0.1 like blockmodels TOL_F): log tau_iq = log alpha_q + sum_j sum_l tau_jl [y_ij log pi_ql + (m_ij - y_ij) log(1-pi_ql)] (+ the same with transposed roles for directed networks); tau is clamped into [0.1/n, 1-0.1/n] and renormalised (blockmodels TOL_P/n).
- PL = E_q log p(y, Z; theta) = sum_iq tau_iq log alpha_q + sum_ql [A_ql log pi_ql + (N_ql - A_ql) log(1-pi_ql)] (halved if undirected); H = -sum tau log tau; lower bound J = PL + H. Iterate until J improves by less than 1e-5 (TOL_EM) or 1000 iterations.
- ICL (blockmodels/sbm definition, NOT greed's exact ICL): ICL = PL - 1/2 [ npar log(Ndata) + (Q-1) log n ], npar = Q^2 (directed) or Q(Q+1)/2 (undirected), Ndata = number of observed ordered pairs (directed) or unordered pairs (undirected), n = number of actors. Note that PL (not J) enters, i.e. the entropy is not added.

## Q selection and initialisation
`blocks(Q)` fits one Q; `maxblocks(Qmax)` (with optional `minblocks()`, default 1) fits every Q in minblocks..Qmax and reports the ICL table; the reported solution is the Q with maximal ICL (ties: smaller Q). Default maxblocks is min(10, n).
Starts per Q (all run through the VEM, the one with the largest ICL is kept, as blockmodels chooses by ICL among its candidate starts):
1. Spectral clustering: eigenvectors with the largest |eigenvalue| of the symmetrised residual matrix of the one-block model, (A + A')/2 - p, via orthogonal subspace iteration in C (cost O(m k) per iteration on the sparse adjacency lists; no dense n x n eigen decomposition); rows of the first Q eigenvectors clustered by k-means with farthest-point start (deterministic, like blockmodels `blockmodelskmeans`).
2. Splits of the best (Q-1)-solution: each block bisected by k-means(2) in the embedding (the forward exploration of blockmodels).
3. `restarts(#)` extra starts (default 10): k-means++ on the embedding (half) and random partitions (half), from the RNG seeded by `seed()`.
4. After the forward pass, a backward pass: merges of the best (Q+1)-solution (all pairs, or the 20 most similar rows/columns of pi), then one more sweep of splits. The search always runs one level above the largest Q asked for so that merges also improve that Q.
Deviations from blockmodels: no exploration_factor stopping rule (we fit all Q in the requested range); the embedding is of the residual matrix by subspace iteration, not of blockmodels' sigmoid-of-residual-Gram Laplacian (an n^3 dense decomposition); fixed passes instead of the forward-backward walk with a stopping rule.

## Syntax
`nwsbm netname, [blocks(#) | maxblocks(#) [minblocks(#)]] [model(bernoulli)] [seed(#) restarts(#) generate(stub) replace]`; `nwsbm` alone replays; `nwsbm version`. `blocks(#1 #2)` (LBM) and model(poisson|gaussian) are accepted by the parser but stop with "not available in this version" until milestone 2+. Directedness comes from the network (no `directed` option): the network's flag decides SBM vs SBM_sym.

## Output
- generate(stub): `stub` = hard block (1..Q, MAP of tau), `stub_p1..stub_pQ` = posterior membership probabilities tau (replace allowed); needs the dataset to have >= n rows (rows are added, as other nwcommands, `_nwsyncrows`).
- e(): cmd, model, membership (SBM / SBM_sym), network, directed (0/1), N (actors), Q, ICL (selected), J, PL, H, npar, n_missing, seed, restarts; matrices e(pi) (Q x Q block connectivity), e(alpha) (1 x Q block proportions), e(blocksize) (1 x Q), e(ICLtable) (Qs x [Q ICL J PL H npar]). No e(b)/e(V) (no ereturn post).
- Mata: tau of the last fit is kept in the global `__nwsbm_tau` until the next fit, so `generate()` also works on a replay. Stata matrices only for small objects (BE matrix limit 800 x 800), never for tau.

## Plugin interface (internal)
Frame `__nwsbm_native`; input block in column v1 (header, then the n x n adjacency with -1 for missing), output in v2 (per Q: ICL, J, PL, H, npar, tau n x Q (row major), pi Q x Q, alpha Q). A version handshake (code 99) protects against stale plugins; the number never appears in user text.

## Validation plan (milestone 1)
1. Agreement with R sbm/blockmodels on fixtures: planted-partition networks (n = 60, 150; Q = 2, 3, 4; undirected and directed) and a network with missing ties; compare: Q selected, adjusted Rand index of the MAP memberships vs R (>= 0.95), ICL (relative tolerance stated in the test), pi after matching blocks.
2. Recovery of planted partitions across contrast strengths (ARI vs truth, by strength).
3. Error paths (two-mode, valued, bad options, no plugin).
4. Timing vs R at n = 200, 500, 1000.

## Results of milestone 1 (macOS, Stata 17 BE)
- 8 fixtures vs R sbm 1.0.0: same Q, ARI 1.0000, ICL at the selected Q equal to the 3rd decimal (relative 1e-5), connectivity matrices identical to 4 decimals; other Q: never worse than R by more than 0.3 percent, sometimes better.
- Timing (Q = 1..6, 10 extra starts per Q): n = 200: 0.4 s (R 0.56 s), n = 500: 3.0 s (R 4.2 s), n = 1000: 8.3 s (R 18.1 s); n = 2000 sparse: 12 s with restarts(2).
- Recovery: see cscripts/test_nwsbm_recovery.do (planted n = 100, 3 blocks, p_between 0.10: perfect from p_within 0.5, near-perfect directed from 0.3; below about 0.2 no method recovers the blocks).
