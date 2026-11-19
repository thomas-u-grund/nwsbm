* Degree-corrected stochastic block model (model(poisson) degreecorrected) against R greed 0.6.2 DcSbm on the fixtures of dev/ref_dc.R.
* greed's ICL is the exact ICL (priors) and differs from the blockmodels-type ICL used here, so the numbers of ICL are not compared; compared are
* the selected number of blocks and the memberships (ARI >= 0.95 against greed, >= 0.9 against the planted blocks). The model without degree
* correction must select more blocks than the planted ones on these degree-heterogeneous networks (the reason for the correction).
do cscripts/nwsbm_testlib.do
local fix "dc_u_n150_Q3 150 0 dc_d_n120_Q2 120 1 dc_u_n200_Q4 200 0"
local fail 0
local nfix : word count `fix'
forvalues f = 1(3)`nfix' {
	local nm : word `f' of `fix'
	local n : word `=`f'+1' of `fix'
	local dir : word `=`f'+2' of `fix'
	sbmt_mknetv g, ref(`nm') n(`n') dir(`dir')
	preserve
	import delimited "$SBM_REPO/dev/ref/`nm'_fit.csv", clear varnames(1)
	local rq = q[1]
	import delimited "$SBM_REPO/dev/ref/`nm'_membership.csv", clear varnames(1)
	rename block rb
	tempfile m
	save `m'
	import delimited "$SBM_REPO/dev/ref/`nm'_truth.csv", clear varnames(1)
	merge 1:1 id using `m', nogen
	tempfile m2
	save `m2'
	restore
	nwsbm g, maxblocks(7) model(poisson) degreecorrected seed(5) generate(mb) replace
	local ok 1
	local qdc = e(Q)
	if `qdc' != `rq' {
		di as error "`nm': Q selected `qdc' but greed selected `rq'"
		local ok 0
	}
	capture drop id
	capture drop rb
	capture drop truth
	gen id = _n
	merge 1:1 id using `m2', nogen
	sbmt_ari mb rb
	local a1 = r(ari)
	sbmt_ari mb truth
	local a2 = r(ari)
	if missing(`a1') | `a1' < 0.95 | missing(`a2') | `a2' < 0.9 {
		di as error "`nm': ARI vs greed `a1', vs truth `a2'"
		local ok 0
	}
	nwsbm g, maxblocks(7) model(poisson) seed(5)
	local qp = e(Q)
	if `qp' <= `rq' {
		di as error "`nm': plain Poisson model selects `qp' blocks, expected more than `rq'"
		local ok 0
	}
	di as text "`nm': Q greed=`rq' DC=`qdc' plain Poisson=`qp'  ARI vs greed=" %6.4f `a1' "  vs truth=" %6.4f `a2'
	if `ok' di as result "PASS `nm'"
	else local fail = `fail' + 1
}
if `fail' == 0 di as result "PASS all degree-corrected fixtures agree with greed"
else {
	di as error "FAIL: `fail' fixtures differ"
	exit 9
}
