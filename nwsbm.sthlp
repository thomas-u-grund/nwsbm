{smcl}
{* *! version 0.3.0  19nov2026 author: Thomas Grund}{...}
{title:Title}

{p2colset 9 17 24 2}{...}
{p2col :nwsbm  {hline 2}}Stochastic block models for networks, fitted by variational EM{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:nwsbm} {it:netname} [{cmd:,} {it:options}]

{p 8 16 2}
{cmd:nwsbm}  [{cmd:,} {opt gen:erate(stub)} {opt replace}]  (replay after a fit)

{p 8 16 2}
{cmd:nwsbm version}

{pstd}
{it:netname} is a network of the nwcommands family: a one-mode network (undirected or directed, binary or valued) or a two-mode network. Self-loops of a one-mode network are ignored; ties with a missing value are left out of the likelihood. The network's own flag decides between the undirected model (one symmetric connectivity matrix) and the directed model, and between the one-mode and the two-mode (latent block) model.

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt block:s(#)}}fit exactly this number of blocks; two numbers, {cmd:blocks(}{it:#1 #2}{cmd:)}, for a two-mode network (blocks of mode 1 and of mode 2){p_end}
{synopt :{opt maxb:locks(#)}}fit every number of blocks from {opt minblocks()} to this value and select by ICL; default {cmd:min(10, actors)} when {opt blocks()} is not given; two numbers for a two-mode network (all pairs are fitted){p_end}
{synopt :{opt minb:locks(#)}}smallest number of blocks fitted with {opt maxblocks()}; default {cmd:1}; two numbers for a two-mode network{p_end}
{synopt :{opt mod:el(name)}}tie model: {cmd:bernoulli} (default; ties 0 or 1), {cmd:poisson} (non-negative integer counts) or {cmd:gaussian} (real values){p_end}
{synopt :{opt degree:corrected}}degree-corrected block model; one-mode networks, with {cmd:model(poisson)}{p_end}
{synopt :{opt sig:ma2(#)}}fix the variance of the Gaussian model; default: one common variance, estimated{p_end}
{synopt :{opt icl(type)}}{cmd:blockmodels} (default) or {cmd:exact}, the exact ICL with conjugate priors of R {cmd:greed}{p_end}
{synopt :{opt cov:ariates(netlist)}}dyadic covariates: networks of the nwcommands session with the same actors (for a two-mode network, the same two modes); the coefficients are estimated with the block parameters{p_end}
{synopt :{opt absd:iff(varlist)}}covariates {it:|x_i - x_j|} from numeric node variables (one covariate per variable){p_end}
{synopt :{opt sam:e(varlist)}}covariates {it:1(x_i = x_j)} from node variables (numeric or string){p_end}
{synopt :{opt sign:eddiff(varlist)}}covariates {it:x_i - x_j} from numeric node variables; directed networks only{p_end}
{synopt :{opt seed(#)}}random-number seed of the starts; default random{p_end}
{synopt :{opt rest:arts(#)}}extra random starts per number of blocks, in addition to the spectral, split and merge starts; default {cmd:10} for networks with up to 500 actors, {cmd:4} up to 1500, {cmd:2} above (and at most {cmd:4} for the Bernoulli model with covariates){p_end}
{synopt :{opt thr:eads(#)}}number of threads; default {cmd:0}, one per core. The result does not depend on it{p_end}
{synopt :{opt fast}}at most two random starts per number of blocks and no merge pass; for large networks{p_end}
{synopt :{opt gen:erate(stub)}}save the memberships: {it:stub} holds the most probable block, {it:stub}{cmd:_p1}, {it:stub}{cmd:_p2}, ... the posterior membership probabilities (for a two-mode network the block and probabilities of the node's own mode){p_end}
{synopt :{opt replace}}allow {opt generate()} to overwrite existing variables{p_end}
{synopt :{opt nodis:play}}suppress the output{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:nwsbm} fits stochastic block models (Holland, Laskey and Leinhardt 1983; Nowicki and Snijders 2001): every actor belongs to one of Q blocks, and the distribution of a tie depends only on the blocks of the two actors. For two-mode networks it fits the latent block model (Govaert and Nadif 2003), with its own blocks for the actors of each mode. The models are fitted by variational EM with a mean-field posterior for the memberships, as in the R packages {cmd:blockmodels} and {cmd:sbm} (Leger; Chiquet, Donnet and Barbillon), whose definitions of the objective and of the integrated classification likelihood (ICL) it follows. All computation runs in a compiled plugin.

{pstd}
{bf:Edge models.} With {cmd:model(bernoulli)} the tie probability of a pair depends on the two blocks; with {cmd:model(poisson)} the tie value is a count with a rate that depends on the two blocks; with {cmd:model(gaussian)} it is a normal variable whose mean depends on the two blocks and whose variance is common to all pairs (estimated, or fixed with {opt sigma2()}). Zero ties of a valued network count as observed zeros.

{pstd}
{bf:Covariates.} With {opt covariates()}, {opt absdiff()}, {opt same()} or {opt signeddiff()} the expected tie of a pair depends on the blocks and on pair-level covariates {it:x_ij1, ..., x_ijp}: the link function of the mean (logit for {cmd:bernoulli}, log for {cmd:poisson}, identity for {cmd:gaussian}) equals a block parameter plus the sum of the coefficients times the covariates, as in the covariate models of {cmd:blockmodels} and {cmd:sbm}. The coefficients are common to all blocks; each adds one parameter to the ICL. Node variables are read from the dataset: row {it:i} is actor {it:i} of the network (as with {opt generate()}); for a two-mode network the variable is compared between an actor of mode 1 and an actor of mode 2. Covariate networks must contain no missing values, must be symmetric when the network is undirected and must not be constant over the pairs. The block parameters {cmd:e(pi)} are then those of a pair with all covariates equal to 0 (response scale); the coefficients are in {cmd:e(beta)} and are displayed after the block table (no standard errors). Missing ties are left out of the covariate models as of the models without covariates (the R covariate models for the Bernoulli and Poisson edges count them as zeros). Degree correction and {opt icl(exact)} are not available with covariates.

{pstd}
{it:Computing time and memory.} Covariates make the pair terms depend on the pair, so the plugin holds dense {it:n x n} matrices (about {it:8 (3 + p + 3T) n^2} bytes for {it:p} covariates and {it:T} threads) instead of the sparse adjacency lists. The Poisson and Gaussian models need {it:O(n^2 Q)} operations per sweep. The Bernoulli model needs {it:O(n^2 Q^2)} (one logarithm per pair and pair of blocks) and is the slowest: on one thread, 100 actors with {cmd:maxblocks(5)} take about 10 seconds and 200 actors about 40 seconds on a recent laptop; use more threads, {opt fast} or a smaller {opt maxblocks()}. A note is shown above 1000 actors; networks above about 3000 actors are refused.

{pstd}
{bf:Degree correction.} With {opt degreecorrected} the Poisson rate of a pair is the product of a degree parameter of each actor and a rate that depends on the two blocks (Karrer and Newman 2011), so that actors with very different degrees can share a block. The degree parameters are profiled out (they are proportional to the observed degrees); a network with strongly varying degrees is otherwise split into blocks by degree. A directed network has an out-degree and an in-degree parameter.

{pstd}
{bf:Choosing the number of blocks.} With {opt maxblocks()} every Q in the range is fitted and the Q with the largest ICL is reported (ties go to the smaller Q); the table shows ICL, the variational lower bound J, the expected complete-data log likelihood PL and the number of parameters for every Q (every pair of block numbers for a two-mode network). Each Q is started from a spectral clustering of the residual of the one-block model, from splits of the best solution with one block fewer, from merges of the best solution with one block more, and from {opt restarts()} random starts; the fit with the largest ICL is kept. The starts are evaluated by {opt threads()} threads; each start is deterministic given {opt seed()}, so the result is the same for any number of threads.

{pstd}
The ICL of {cmd:blockmodels} is
{it:ICL = PL - 1/2 [ P log(N) + (Q - 1) log(n) ]}
for a one-mode network and
{it:ICL = PL - 1/2 [ P log(N) + (Q1 - 1) log(n1) + (Q2 - 1) log(n2) ]}
for a two-mode network, where {it:PL} is the expected complete-data log likelihood under the fitted membership probabilities, {it:P} the number of parameters of the edge model (Q*Q for a directed network, Q(Q+1)/2 for an undirected one, Q1*Q2 for a two-mode network, plus one for the estimated variance of the Gaussian model), {it:N} the number of observed pairs (ordered pairs when directed, unordered pairs when undirected, cells when two-mode) and {it:n} the number of actors of the mode. The entropy of the memberships is not added. Larger is better.

{pstd}
{opt icl(exact)} replaces it by the exact ICL with conjugate priors (Dirichlet(1) for the block proportions, Beta(1,1) for the tie probabilities of the Bernoulli model, the priors of R {cmd:greed} for the degree-corrected Poisson model), evaluated at the most probable partition of every fitted Q, and selects Q by it. It is available for {cmd:model(bernoulli)} and for {cmd:model(poisson)} with {opt degreecorrected}, one-mode networks without missing ties. The table then holds the exact ICL in the ICL column and the criterion of {cmd:blockmodels} in {cmd:e(ICLbm)}.

{pstd}
The block labels are arbitrary; blocks are not ordered by size. Without covariates the cost of an iteration grows with the number of ties and blocks, not with the square of the number of actors, but each Q is started many times: reduce {opt restarts()} or use {opt fast} for large networks.

{title:Results}

{pstd}
Scalars: {cmd:e(N)} actors, {cmd:e(Q)} selected blocks (for a two-mode network {cmd:e(Q1)} and {cmd:e(Q2)}, {cmd:e(N1)} and {cmd:e(N2)} actors of the two modes; {cmd:e(Q)} equals {cmd:e(Q1)}), {cmd:e(ICL)}, {cmd:e(J)} (lower bound), {cmd:e(PL)}, {cmd:e(H)} (entropy of the memberships), {cmd:e(npar)}, {cmd:e(directed)}, {cmd:e(n_missing)} pairs with a missing tie, {cmd:e(sigma2)} (Gaussian), {cmd:e(ICLbm)} (with {cmd:icl(exact)}), {cmd:e(seed)}, {cmd:e(restarts)} (the number used), {cmd:e(ncov)} number of covariates.
Macros: {cmd:e(cmd)}, {cmd:e(model)}, {cmd:e(network)}, {cmd:e(degreecorrected)}, {cmd:e(icltype)}, {cmd:e(covariates)} (names of the covariates), {cmd:e(membership)} ({cmd:SBM} directed, {cmd:SBM_sym} undirected, {cmd:LBM} two-mode, the names used by {cmd:blockmodels}).
Matrices: {cmd:e(pi)} (block connectivity from a row block to a column block: tie probability, expected tie value or mean; Q x Q, or Q1 x Q2 for a two-mode network; for the degree-corrected model the rate that multiplies the degree parameters), {cmd:e(alpha)} (block proportions; {cmd:e(alpha2)} for mode 2), {cmd:e(blocksize)} (actors by most probable block; {cmd:e(blocksize2)}), {cmd:e(ICLtable)} (one row per fitted Q or pair of block numbers), {cmd:e(beta)} (1 x p coefficients of the covariates; with the names {it:netname}, {cmd:abs_}{it:var}, {cmd:same_}{it:var}, {cmd:diff_}{it:var}).

{pstd}
The membership probabilities are not stored in a Stata matrix (a matrix holds at most 800 rows in Stata BE); use {opt generate()}, also after a fit by typing {cmd:nwsbm, generate(}{it:stub}{cmd:)}. In a two-mode network the nwcommands convention puts the columns of the wide matrix into mode 1.

{title:Postestimation}

{p 8 16 2}{cmd:estat icl}{p_end}
{p 8 16 2}{cmd:estat blocks}{p_end}

{pstd}
{cmd:estat icl} lists the ICL, J, PL, H and number of parameters of every fitted Q (returned in {cmd:r(table)}); {cmd:estat blocks} lists size, proportion and within-block connectivity (two-mode: mean tie value) of the blocks (returned in {cmd:r(blocks)}, {cmd:r(blocks2)}). {help nwsbmplot} draws the blocked adjacency matrix and the block graph.

{title:Examples}

{phang2}{cmd:. nwrandom 120, prob(0.1) name(g) undirected}{p_end}
{phang2}{cmd:. nwsbm g, maxblocks(6) seed(1) generate(blk)}{p_end}
{phang2}{cmd:. nwsbm g, blocks(3) seed(1) generate(b3)}{p_end}
{phang2}{cmd:. matrix list e(pi)}{p_end}
{phang2}{cmd:. nwsbm g, model(poisson) degreecorrected maxblocks(6) icl(exact)}{p_end}
{phang2}{cmd:. nwsbm g, model(poisson) maxblocks(5) covariates(distance) absdiff(age) same(gender)}{p_end}
{phang2}{cmd:. matrix list e(beta)}{p_end}
{phang2}{cmd:. estat blocks}{p_end}
{phang2}{cmd:. nwsbmplot, shade}{p_end}

{title:Not in this version}

{pstd}
Multiplex and dynamic networks, the degree-corrected two-mode model, covariates with the degree-corrected model and a Gaussian model with block-specific variances are not available.

{title:References}

{phang}
Holland, P. W., K. B. Laskey and S. Leinhardt. 1983. Stochastic blockmodels: first steps. {it:Social Networks} 5: 109-137.

{phang}
Nowicki, K. and T. A. B. Snijders. 2001. Estimation and prediction for stochastic blockstructures. {it:Journal of the American Statistical Association} 96: 1077-1087.

{phang}
Govaert, G. and M. Nadif. 2003. Clustering with block mixture models. {it:Pattern Recognition} 36: 463-473.

{phang}
Karrer, B. and M. E. J. Newman. 2011. Stochastic blockmodels and community structure in networks. {it:Physical Review E} 83: 016107.

{phang}
Daudin, J.-J., F. Picard and S. Robin. 2008. A mixture model for random graphs. {it:Statistics and Computing} 18: 173-183.

{phang}
Biernacki, C., G. Celeux and G. Govaert. 2000. Assessing a mixture model for clustering with the integrated completed likelihood. {it:IEEE Transactions on Pattern Analysis and Machine Intelligence} 22: 719-725.

{title:Author}

{pstd}
Thomas U. Grund, RWTH Aachen University. nwsbm is part of the nwcommands family and requires the core package nwcommands.

{title:Also see}

{pstd}
{help nwsbmplot}, {help nwset}, {help nwrandom}

{hline}
{pstd}{it:nwsbm version 0.3.0, 19 November 2026}{p_end}{* nwversion}
