# nwsbm design (internal; not shipped)

Owner: Thomas U. Grund (sole author). GPL-3. Part of the nwcommands family; depends only on the core. Native C plugin for all heavy computation (no Mata fallback). Mata library built with Stata 17.

## References (match as closely as possible; document deviations)
- R `blockmodels` 1.2.1 (Leger) and `sbm` 1.0.0 (Chiquet, Donnet, Barbillon): variational EM, ICL, model names bernoulli / poisson / gaussian, membership types SBM (directed), SBM_sym (undirected), LBM (two-mode). `sbm::estimateSimpleSBM` is a wrapper of blockmodels, so the two share every number (checked in sources: `ICL <- PL - .5*(n_parameters*log(data_number) + (Q-1)*log(n))`).
- R `greed` 0.6.2 (Come, Randriamanamihaga, Ouattara, Rey): exact ICL with conjugate priors and greedy merge/swap search; degree-corrected SBM. A later milestone.

## Models (stages)
| Milestone | Model | Status |
|---|---|---|
| 0.1.0 | Bernoulli SBM, undirected (blockmodels `SBM_sym`) and directed (`SBM`) one-mode networks, missing ties allowed | this version |
| 0.2.0 | Poisson and Gaussian SBM (valued networks), degree-corrected Poisson SBM, latent block model (two-mode), exact ICL, threads, estat, nwsbmplot | this version |
| later | covariates, degree-corrected two-mode, ZIP | TODO |

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

## Milestone 2 (0.2.0)
Same engine for all models: the plugin works with the sufficient statistics ZA (sum of y) and ZO (number of observed cells; degree-weighted for the degree-corrected model) per block pair, so a model is a triple (M step, E-step coefficients, parameter count). The E step is log tau_iq += sum_l [a_il U_ql + n_il V_ql] with U, V from the parameters (Bernoulli: U = logit p, V = log(1-p); Poisson: U = log lambda, V = -lambda; Gaussian: U = mu/s2, V = -mu^2/(2 s2)); a = sum_j tau_jl y_ij from the sparse adjacency, n = observed cells from the missing list.
- Poisson (blockmodels zpoisson): lambda = ZA/ZO, PL = sum(ZA log lambda - ZO lambda) - sum log y! (halved for undirected; the log-factorial constant is in PL, as in R, so ICL equals R's). npar as Bernoulli.
- Gaussian (blockmodels zgaussian): one free variance common to all pairs, sigma2 = (sum y^2 + sum(mu^2 ZO - 2 mu ZA)) / (observed ordered cells), PL = -nobs/2 log(2 pi sigma2) - resid/(2 sigma2) (halved if undirected); npar + 1. R has no fixed-variance option; `sigma2(#)` is an extension (no extra parameter, variance not updated).
- Degree correction (not in blockmodels/sbm; R greed `DcSbm`): y_ij ~ Poisson(theta_i theta_j omega_ql) with theta profiled out, theta_i = observed degree (out-degree times in-degree when directed); omega = ZA / ZD, ZD = sum tau_iq tau_jl d_i d_j m_ij; PL adds sum y log(d_i d_j) - sum log y!. npar as the plain model (the degree parameters are data, same for every Q). Starts: spectral embedding of the degree-corrected residual with rows normalised to the unit sphere (as in SCORE).
- Latent block model (blockmodels LBM): tau1 (n1 x Q1), tau2 (n2 x Q2), rectangular adjacency and its transpose as two sparse structures; E step updates both sides from the old tau, entropy and alpha terms of both sides; ICL = PL - 1/2 [npar log(Ndata) + (Q1-1) log n1 + (Q2-1) log n2], npar = Q1 Q2 (+1 Gaussian), Ndata = observed cells. Starts: SVD-like embedding (subspace iteration on B B', columns C = B'U), k-means of rows and columns separately; the grid (Q1 x Q2) is fitted completely (splits of a row block from (Q1-1,Q2), of a column block from (Q1,Q2-1), merges from the larger neighbours, random starts); R explores a path with a stopping rule.
- Exact ICL (greed 0.6.2; reproduced to 1e-3 absolute on the fixtures): ICL = emission term + lgamma(K) + sum lgamma(n_k + 1) - lgamma(N + K) (Dirichlet(1)); Bernoulli SBM emission sum over block pairs of log B(1 + x, 1 + cc - x) (undirected: pairs k<l with cc = n_k n_l, diagonal x/2 with cc = n_k (n_k-1)/2; directed: all ordered pairs with cc = n_k n_l); DcSbm emission: sum_k [lgamma(n_k) - lgamma(n_k + d_k) + d_k log n_k] (twice for directed, in and out degrees) + sum over block pairs [lgamma(x + 1) - (x + 1) log(p cc + 1)] + sum_i lgamma(d_i + 1) - sum lgamma(y + 1) (/2 undirected) + links log p, p = sum y / n^2. Evaluated at the hard partition (empty blocks dropped).
- Threads: starts of one Q (spectral, splits, random) are a batch; each is deterministic (the starts are built sequentially from one RNG stream, then fitted in parallel with a per-thread workspace); results are collected in start order, ties go to the earlier start, so any thread count gives the same output. Merge and split passes are batches too. `fast` = at most 2 random starts and no merge/second split pass.

### Results of milestone 2 (macOS, Stata 17 BE; R sbm 1.0.0, blockmodels 1.2.1, greed 0.6.2)
- Poisson and Gaussian vs sbm: 6 fixtures (undirected, directed, missing ties; n = 80 and 100; Q = 2, 3): same Q, ARI 1.0000, ICL at every stored Q equal to R to a relative 5e-8 (never worse), connectivity matrices and the Gaussian variance equal to 4 decimals.
- Two-mode (LBM) vs sbm: 4 fixtures (Bernoulli, Poisson, Gaussian 40 x 60; Bernoulli with missing 50 x 50): same (Q1, Q2), ARI 1.0000, ICL at every stored cell equal to a relative 5e-8, connectivity equal to 1e-4.
- Degree-corrected vs greed DcSbm (3 fixtures): same number of blocks, ARI 1.0000 / 1.0000 / 0.987 against greed; the plain Poisson SBM selects 6 blocks on the same networks (degree heterogeneity). Exact ICL of nwsbm vs greed at identical partitions: equal to a relative 1e-6 (5 fixtures incl. directed binary SBM); in dc_u_n200_Q4 the partitions differ by two nodes and greed's ICL is larger by 0.05 (7e-6 relative).
- Recovery: Poisson, Gaussian, degree-corrected (strong degree heterogeneity: ARI 0.90), two-mode Bernoulli/Poisson/Gaussian: right Q and ARI >= 0.90 (cscripts/test_nwsbm_recovery2.do).
- Timing (Stata 17, one thread, same Q range 1..6 with 10 extra starts, load average about 6 on the machine): n = 200: 0.64 s (R 0.56 s), 500: 5.7 s (R 4.2 s), 1000: 12.6 s (R 18.1 s); threads(4): 0.27 / 2.9 / 6.3 s; `fast` one thread: 0.19 / 2.3 / 4.2 s. In 0.1.0 (idle machine) the single-thread times were 0.4 / 3.0 / 8.3 s; the difference is probably the load, to be re-measured on an idle machine.
