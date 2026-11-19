* Agreement of nwsbm on two-mode networks (latent block model) with R sbm 1.0.0 estimateBipartiteSBM on the fixtures of dev/ref_valued.R.
* nwset puts the columns of the wide matrix into mode 1 and the rows into mode 2, so R's row blocks are nwsbm's blocks of mode 2.
* Per fixture: selected (Q1, Q2) equal; ARI of both memberships >= 0.95; ICL at the selected cell within 1e-5 (relative) and at every other cell R stored
* not worse than R by 0.3 percent; connectivity matrix within 0.05 after matching the blocks by overlap.
do cscripts/nwsbm_testlib.do
local fix "lbm_ber_40x60 40 60 bernoulli lbm_po_40x60 40 60 poisson lbm_ga_40x60 40 60 gaussian lbm_ber_miss_50x50 50 50 bernoulli"
local fail 0
local nfix : word count `fix'
forvalues f = 1(4)`nfix' {
	local nm : word `f' of `fix'
	local n1 : word `=`f'+1' of `fix'
	local n2 : word `=`f'+2' of `fix'
	local mod : word `=`f'+3' of `fix'
	sbmt_mknet2 g, ref(`nm') n1(`n1') n2(`n2') valued
	preserve
	import delimited "$SBM_REPO/dev/ref/`nm'_fit.csv", clear varnames(1)
	local rq1 = q1[1]
	local rq2 = q2[1]
	local ricl = icl[1]
	if "`mod'" == "gaussian" local rs2 = sigma2[1]
	import delimited "$SBM_REPO/dev/ref/`nm'_icl.csv", clear varnames(1)
	local nrow = _N
	local mx1 0
	local mx2 0
	forvalues r = 1/`nrow' {
		local a = rowblocks[`r']
		local b = colblocks[`r']
		local rI_`a'_`b' = icl[`r']
		if `a' > `mx1' local mx1 = `a'
		if `b' > `mx2' local mx2 = `b'
	}
	import delimited "$SBM_REPO/dev/ref/`nm'_pi.csv", clear varnames(1)
	mkmat _all, matrix(RPI)
	restore
	* nwsbm: mode 1 = columns (R), mode 2 = rows (R)
	nwsbm g, model(`mod') maxblocks(`mx2' `mx1') seed(7) generate(mb) replace
	local ok 1
	if e(Q1) != `rq2' | e(Q2) != `rq1' {
		di as error "`nm': (Q1,Q2) selected `=e(Q1)' x `=e(Q2)' but R selected `rq2' x `rq1' (columns x rows)"
		local ok 0
	}
	_nwsyntax g, max(1)
	mata: sbm_setidx(`netobj'->get_modes())
	mata: st_local("n1c", strofreal(rows(__nwsbm_i1)))
	if `n1c' != `n2' {
		di as error "`nm': mode 1 has `n1c' nodes, expected `n2'"
		local ok 0
	}
	* memberships: mode 1 nodes against R's columns, mode 2 nodes against R's rows
	mata: sbmt_ari_csv("mb", "$SBM_REPO/dev/ref/`nm'_membership2.csv", __nwsbm_i1)
	local ari2 = scalar(__ari)
	mata: sbmt_ari_csv("mb", "$SBM_REPO/dev/ref/`nm'_membership.csv", __nwsbm_i2)
	local ari1 = scalar(__ari)
	if missing(`ari1') | missing(`ari2') | `ari1' < 0.95 | `ari2' < 0.95 {
		di as error "`nm': ARI rows `ari1' columns `ari2' < 0.95"
		local ok 0
	}
	matrix T = e(ICLtable)
	local maxrel 0
	local relsel 0
	forvalues a = 1/`mx1' {
		forvalues b = 1/`mx2' {
			if "`rI_`a'_`b''" == "" continue
			local rv = `rI_`a'_`b''
			* nwsbm cell (Q1 = b, Q2 = a)
			local rr = (`b' - 1) * `mx1' + `a'
			local worse = (`rv' - T[`rr', 3]) / abs(`rv')
			if `worse' > `maxrel' local maxrel = `worse'
			if `a' == `rq1' & `b' == `rq2' local relsel = abs(T[`rr', 3] - `rv') / abs(`rv')
		}
	}
	if `maxrel' > 0.003 | `relsel' > 0.00001 {
		di as error "`nm': ICL worse than R by `maxrel' (relative; at the selected cell: difference `relsel')"
		local ok 0
	}
	if "`mod'" == "gaussian" {
		if abs(e(sigma2) - `rs2') / `rs2' > 1e-4 {
			di as error "`nm': variance `=e(sigma2)' vs R `rs2'"
			local ok 0
		}
	}
	mata: sbmt_pidiff2("mb", "$SBM_REPO/dev/ref/`nm'_membership2.csv", "$SBM_REPO/dev/ref/`nm'_membership.csv", __nwsbm_i1, __nwsbm_i2)
	di as text "`nm': `n1'x`n2'  cell R=`rq2'x`rq1' (cols x rows) nwsbm=`=e(Q1)'x`=e(Q2)'  ARI rows=" %6.4f `ari1' " cols=" %6.4f `ari2' "  ICL R=" %10.3f `ricl' " nwsbm=" %10.3f e(ICL) "  max rel ICL shortfall=" %9.2e `maxrel' "  max |pi diff|=" %7.4f scalar(__pidiff)
	if scalar(__pidiff) > 0.05 {
		di as error "`nm': connectivity differs by " scalar(__pidiff)
		local ok 0
	}
	if `ok' di as result "PASS `nm'"
	else local fail = `fail' + 1
}
if `fail' == 0 di as result "PASS all two-mode fixtures agree with R sbm"
else {
	di as error "FAIL: `fail' fixtures differ"
	exit 9
}
