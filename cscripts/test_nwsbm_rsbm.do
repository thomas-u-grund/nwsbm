* Agreement of nwsbm with R sbm 1.0.0 / blockmodels 1.2.1 on the fixtures of dev/ref_sbm.R.
* Per fixture: Q selected equal; adjusted Rand index of the memberships (nwsbm vs R) >= 0.95; ICL at the selected Q within 1e-5 and at every
* other Q that R stored not worse than R by more than 0.3% (relative; away from the truth the two programs can end in different local optima; nwsbm
* also searches merges of the next larger solution and is sometimes better); block connectivity matrix (blocks matched by
* overlap) within 0.02 absolute.
do cscripts/nwsbm_testlib.do
local fix "pp_u_n60_Q2 60 0 pp_u_n60_Q3 60 0 pp_u_n150_Q3 150 0 pp_u_n150_Q4 150 0 pp_d_n100_Q3 100 1 pp_d_n150_Q2 150 1 miss_u_n80_Q3 80 0 miss_d_n80_Q2 80 1"
local fail 0
local nfix : word count `fix'
forvalues f = 1(3)`nfix' {
	local nm : word `f' of `fix'
	local n : word `=`f'+1' of `fix'
	local dir : word `=`f'+2' of `fix'
	sbmt_mknet g, ref(`nm') n(`n') dir(`dir')
	preserve
	import delimited "$SBM_REPO/dev/ref/`nm'_fit.csv", clear varnames(1)
	local rq = q[1]
	local ricl = icl[1]
	import delimited "$SBM_REPO/dev/ref/`nm'_icl.csv", clear varnames(1)
	local nrow = _N
	forvalues r = 1/`nrow' {
		local rQ`=nbblocks[`r']' = icl[`r']
	}
	local rmaxq = nbblocks[`nrow']
	import delimited "$SBM_REPO/dev/ref/`nm'_membership.csv", clear varnames(1)
	rename block rblock
	tempfile rm
	save `rm'
	import delimited "$SBM_REPO/dev/ref/`nm'_pi.csv", clear varnames(1)
	mkmat _all, matrix(RPI)
	restore
	nwsbm g, maxblocks(`rmaxq') seed(7) generate(mb) replace
	local ok 1
	* Q selected
	if e(Q) != `rq' {
		di as error "`nm': Q selected `=e(Q)' but R selected `rq'"
		local ok 0
	}
	* ARI
	capture drop id rblock
	gen id = _n
	merge 1:1 id using `rm', nogen
	sbmt_ari mb rblock
	local ari = r(ari)
	if `ari' < 0.95 {
		di as error "`nm': ARI `ari' < 0.95"
		local ok 0
	}
	* ICL at every stored Q
	matrix T = e(ICLtable)
	local maxrel 0
	local relsel 0
	forvalues q = 1/`rmaxq' {
		local rel = abs(T[`q', 2] - `rQ`q'') / abs(`rQ`q'')
		local worse = (`rQ`q'' - T[`q', 2]) / abs(`rQ`q'')
		if `worse' > `maxrel' local maxrel = `worse'
		if `q' == `rq' local relsel = `rel'
		di as text "   `nm' Q=`q': ICL nwsbm " %10.3f T[`q', 2] "  R " %10.3f `rQ`q''
	}
	if `maxrel' > 0.003 | `relsel' > 0.00001 {
		di as error "`nm': ICL worse than R by `maxrel' (relative; at the selected Q: difference `relsel')"
		local ok 0
	}
	* pi after matching blocks by overlap (only when Q equal)
	local maxpi .
	if e(Q) == `rq' {
		mata: sbmt_pidiff("mb", "rblock", `n')
	}
	di as text "`nm': n=`n' dir=`dir'  Q R=`rq' nwsbm=`=e(Q)'  ARI=" %6.4f `ari' "  ICL R=" %10.3f `ricl' " nwsbm=" %10.3f e(ICL) "  max rel ICL shortfall vs R (Q=1..`rmaxq')=" %9.2e `maxrel' "  max |pi diff|=" %7.4f scalar(__pidiff)
	if scalar(__pidiff) > 0.02 {
		di as error "`nm': connectivity differs by " scalar(__pidiff)
		local ok 0
	}
	if `ok' di as result "PASS `nm'"
	else local fail = `fail' + 1
}
if `fail' == 0 di as result "PASS all fixtures agree with R sbm"
else {
	di as error "FAIL: `fail' fixtures differ"
	exit 9
}
