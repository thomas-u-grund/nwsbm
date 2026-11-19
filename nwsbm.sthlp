{smcl}
{* *! version 0.1.0  19nov2026 author: Thomas Grund}{...}
{title:Title}

{p2colset 9 17 24 2}{...}
{p2col :nwsbm  {hline 2}}Stochastic block model for the ties of a binary network, fitted by variational EM{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:nwsbm} {it:netname} [{cmd:,} {it:options}]

{p 8 16 2}
{cmd:nwsbm}  [{cmd:,} {opt gen:erate(stub)} {opt replace}]  (replay after a fit)

{p 8 16 2}
{cmd:nwsbm version}

{pstd}
{it:netname} is a one-mode network of the nwcommands family, undirected or directed and binary. Self-loops are ignored; ties with a missing value are left out of the likelihood. The network's own flag decides between the undirected model (one connectivity matrix, symmetric) and the directed model.

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt block:s(#)}}fit exactly this number of blocks{p_end}
{synopt :{opt maxb:locks(#)}}fit every number of blocks from {opt minblocks()} to this value and select by ICL; default {cmd:min(10, actors)} when {opt blocks()} is not given{p_end}
{synopt :{opt minb:locks(#)}}smallest number of blocks fitted with {opt maxblocks()}; default {cmd:1}{p_end}
{synopt :{opt mod:el(bernoulli)}}tie model; {cmd:bernoulli} is the only one in this version{p_end}
{synopt :{opt seed(#)}}random-number seed of the starts; default random{p_end}
{synopt :{opt rest:arts(#)}}extra random starts per number of blocks, in addition to the spectral and the split and merge starts; default {cmd:10}{p_end}
{synopt :{opt gen:erate(stub)}}save the memberships: {it:stub} holds the most probable block, {it:stub}{cmd:_p1}, {it:stub}{cmd:_p2}, ... the posterior membership probabilities{p_end}
{synopt :{opt replace}}allow {opt generate()} to overwrite existing variables{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:nwsbm} fits the stochastic block model (Holland, Laskey and Leinhardt 1983; Nowicki and Snijders 2001): every actor belongs to one of Q blocks, and the probability of a tie depends only on the blocks of the two actors. The model is fitted by variational EM with a mean-field posterior for the memberships, as in the R packages {cmd:blockmodels} and {cmd:sbm} (Leger; Chiquet, Donnet and Barbillon), whose definitions of the objective and of the integrated classification likelihood (ICL) it follows. All computation runs in a compiled plugin.

{pstd}
With {opt maxblocks()} the model is fitted for every Q in the range and the Q with the largest ICL is reported (ties go to the smaller Q); the table shows ICL, the variational lower bound J, the expected complete-data log likelihood PL and the number of tie parameters for every Q. Each Q is started from a spectral clustering of the residual of the one-block model, from splits of the best solution with one block fewer, from merges of the best solution with one block more, and from {opt restarts()} random starts; the fit with the largest ICL is kept. {opt seed()} makes the result reproducible.

{pstd}
The ICL is
{it:ICL = PL - 1/2 [ P log(N) + (Q - 1) log(n) ]},
where {it:PL} is the expected complete-data log likelihood under the fitted membership probabilities, {it:P} the number of tie parameters (Q*Q for a directed network, Q(Q+1)/2 for an undirected one), {it:N} the number of observed pairs (ordered pairs when directed, unordered pairs when undirected) and {it:n} the number of actors. This is the criterion of {cmd:blockmodels} and {cmd:sbm}; it is not the exact ICL with conjugate priors of other packages, and the entropy of the memberships is not added. Larger is better.

{pstd}
The block labels are arbitrary; blocks are not ordered by size. Networks with several thousand actors are fitted in seconds to minutes: the cost of an iteration grows with the number of ties and the number of blocks, not with the square of the number of actors, but each Q is started many times, so reduce {opt restarts()} for large networks.

{title:Results}

{pstd}
Scalars: {cmd:e(N)} actors, {cmd:e(Q)} selected blocks, {cmd:e(ICL)}, {cmd:e(J)} (lower bound), {cmd:e(PL)}, {cmd:e(H)} (entropy of the memberships), {cmd:e(npar)} tie parameters, {cmd:e(directed)}, {cmd:e(n_missing)} pairs with a missing tie, {cmd:e(seed)}, {cmd:e(restarts)}.
Macros: {cmd:e(cmd)}, {cmd:e(model)}, {cmd:e(network)}, {cmd:e(membership)} ({cmd:SBM} directed, {cmd:SBM_sym} undirected, the names used by {cmd:blockmodels}).
Matrices: {cmd:e(pi)} (Q x Q block connectivity: probability of a tie from a row block to a column block), {cmd:e(alpha)} (block proportions), {cmd:e(blocksize)} (actors by most probable block), {cmd:e(ICLtable)} (one row per fitted Q: Q, ICL, J, PL, H, npar).

{pstd}
The membership probabilities are not stored in a Stata matrix (a matrix holds at most 800 rows in Stata BE); use {opt generate()}, also after a fit by typing {cmd:nwsbm, generate(}{it:stub}{cmd:)}.

{title:Examples}

{phang2}{cmd:. nwrandom 120, prob(0.1) name(g) undirected}{p_end}
{phang2}{cmd:. nwsbm g, maxblocks(6) seed(1) generate(blk)}{p_end}
{phang2}{cmd:. nwsbm g, blocks(3) seed(1) generate(b3)}{p_end}
{phang2}{cmd:. matrix list e(pi)}{p_end}

{title:Not in this version}

{pstd}
Poisson and Gaussian edge models, degree-corrected blocks, the latent block model for two-mode networks and covariates are not available; {cmd:model(poisson)}, {cmd:model(gaussian)} and {cmd:blocks(}{it:#1 #2}{cmd:)} stop with a message. Valued networks must be dichotomised first.

{title:References}

{phang}
Holland, P. W., K. B. Laskey and S. Leinhardt. 1983. Stochastic blockmodels: first steps. {it:Social Networks} 5: 109-137.

{phang}
Nowicki, K. and T. A. B. Snijders. 2001. Estimation and prediction for stochastic blockstructures. {it:Journal of the American Statistical Association} 96: 1077-1087.

{phang}
Daudin, J.-J., F. Picard and S. Robin. 2008. A mixture model for random graphs. {it:Statistics and Computing} 18: 173-183.

{phang}
Biernacki, C., G. Celeux and G. Govaert. 2000. Assessing a mixture model for clustering with the integrated completed likelihood. {it:IEEE Transactions on Pattern Analysis and Machine Intelligence} 22: 719-725.

{title:Author}

{pstd}
Thomas U. Grund, RWTH Aachen University. nwsbm is part of the nwcommands family and requires the core package nwcommands.

{title:Also see}

{pstd}
{help nwset}, {help nwrandom}

{hline}
{pstd}{it:nwsbm version 0.1.0, 19 November 2026}{p_end}{* nwversion}
