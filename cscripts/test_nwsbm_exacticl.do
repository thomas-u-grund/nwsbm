* icl(exact): the exact ICL of R greed 0.6.2 (conjugate priors, Dirichlet(1) on the proportions) evaluated at the partition nwsbm finds.
* Fixtures of dev/ref_dc.R: binary SBM (undirected, directed; greed Sbm) and degree-corrected Poisson (undirected, directed; greed DcSbm).
* PASS per fixture: nwsbm selects greed's number of blocks by the exact ICL, its partition equals greed's (ARI >= 0.95), and its exact ICL equals
* greed's within 1e-6 (relative) when the partitions are identical, otherwise within 1e-4 (greed's greedy search on the exact criterion can end one or two
* nodes away from the VEM partition, and its exact ICL is then slightly the larger).
do cscripts/nwsbm_testlib.do
local fix "gr_sbm_u_n100_Q3 100 0 bernoulli gr_sbm_d_n100_Q2 100 1 bernoulli dc_u_n150_Q3 150 0 dc dc_d_n120_Q2 120 1 dc dc_u_n200_Q4 200 0 dc"
local fail 0
local nfix : word count `fix'
forvalues f = 1(4)`nfix' {
	local nm : word `f' of `fix'
	local n : word `=`f'+1' of `fix'
	local dir : word `=`f'+2' of `fix'
	local mod : word `=`f'+3' of `fix'
	if "`mod'" == "dc" {
		sbmt_mknetv g, ref(`nm') n(`n') dir(`dir')
		local opt "model(poisson) degreecorrected"
	}
	else {
		sbmt_mknet g, ref(`nm') n(`n') dir(`dir')
		local opt "model(bernoulli)"
	}
	preserve
	import delimited "$SBM_REPO/dev/ref/`nm'_fit.csv", clear varnames(1)
	local rq = q[1]
	local ricl = icl[1]
	import delimited "$SBM_REPO/dev/ref/`nm'_membership.csv", clear varnames(1)
	rename block rb
	tempfile m
	save `m'
	restore
	nwsbm g, `opt' maxblocks(7) icl(exact) seed(5) generate(mb) replace nodisplay
	local ok 1
	capture drop id rb
	gen id = _n
	merge 1:1 id using `m', nogen
	sbmt_ari mb rb
	local ari = r(ari)
	local rel = abs(e(ICL) - `ricl') / abs(`ricl')
	local tol = cond(`ari' > 0.9999, 1e-6, 1e-4)
	if e(Q) != `rq' | `ari' < 0.95 | `rel' > `tol' {
		di as error "`nm': Q `=e(Q)' (greed `rq'), ARI `ari', exact ICL `=e(ICL)' vs greed `ricl' (relative difference `rel')"
		local ok 0
	}
	di as text "`nm': Q greed=`rq' nwsbm=`=e(Q)'  ARI=" %6.4f `ari' "  exact ICL nwsbm=" %10.3f e(ICL) " greed=" %10.3f `ricl' "  rel.diff=" %9.2e `rel'
	if `ok' di as result "PASS `nm'"
	else local fail = `fail' + 1
}
* errors
capture nwsbm g, model(gaussian) icl(exact)
if _rc != 198 {
	di as error "icl(exact) with the Gaussian model must fail with r(198)"
	local fail = `fail' + 1
}
if `fail' == 0 di as result "PASS all exact-ICL fixtures agree with greed"
else {
	di as error "FAIL: `fail' checks differ"
	exit 9
}
