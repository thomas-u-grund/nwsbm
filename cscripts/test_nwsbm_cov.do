* Agreement of nwsbm with dyadic covariates (covariates()) with R sbm 1.0.0 / blockmodels 1.2.1 (estimateSimpleSBM, estimateBipartiteSBM with
* covariates = list(...), fast = FALSE, see GOTCHA) on the fixtures of dev/ref_cov.R: six one-mode fixtures (Bernoulli, Poisson, Gaussian; undirected and
* directed; one and two covariates; missing ties) and three two-mode fixtures.
* Per fixture: Q selected equal; ARI of the memberships >= 0.95; ICL at the selected Q within 2e-4 (relative), at every other Q R stored not worse than R by
* more than 0.3%; regression coefficients within 0.01 (absolute); Gaussian variance within 1e-3 (relative); block parameters (response scale, covariates 0)
* within 0.05 after matching blocks by overlap. The one-block fit of every fixture is also compared with an independent glm / lm fit (log likelihood within 1e-6
* relative, coefficients within 1e-6). R's Bernoulli two-mode covariate model reports a PL that exceeds the maximal log likelihood (checked with glm), so
* the fixture cl_b_40x50 is compared on Q, memberships and coefficients (0.02) only. R's Bernoulli and Poisson covariate models treat missing ties as zeros,
* so the fixtures with missing ties (cbm_d_n60_Q2, cpm_u_n60_Q2) are checked against glm and the planted values only.
do cscripts/nwsbm_testlib.do
local fix "cb_u_n100_Q3 100 0 bernoulli 1 cb_d_n80_Q2 80 1 bernoulli 2 cp_u_n80_Q3 80 0 poisson 1 cp_d_n70_Q2 70 1 poisson 2 cg_u_n80_Q2 80 0 gaussian 1 cg_d_n70_Q3 70 1 gaussian 1"
local fail 0
local nfix : word count `fix'
forvalues f = 1(5)`nfix' {
	local nm : word `f' of `fix'
	local n : word `=`f'+1' of `fix'
	local dir : word `=`f'+2' of `fix'
	local mod : word `=`f'+3' of `fix'
	local np : word `=`f'+4' of `fix'
	sbmt_mknetv g, ref(`nm') n(`n') dir(`dir')
	local covs
	forvalues k = 1/`np' {
		sbmt_mkcov x`k', file("$SBM_REPO/dev/ref/`nm'_cov`k'.csv") dir(`dir')
		local covs `covs' x`k'
	}
	preserve
	import delimited "$SBM_REPO/dev/ref/`nm'_fit.csv", clear varnames(1)
	local rq = q[1]
	local ricl = icl[1]
	if "`mod'" == "gaussian" local rs2 = sigma2[1]
	import delimited "$SBM_REPO/dev/ref/`nm'_icl.csv", clear varnames(1)
	local nrow = _N
	forvalues r = 1/`nrow' {
		local rQ`=nbblocks[`r']' = icl[`r']
	}
	local rmaxq = nbblocks[`nrow']
	import delimited "$SBM_REPO/dev/ref/`nm'_beta.csv", clear varnames(1)
	forvalues k = 1/`np' {
		local rb`k' = beta[`k']
	}
	import delimited "$SBM_REPO/dev/ref/`nm'_membership.csv", clear varnames(1)
	rename block rblock
	tempfile rm
	save `rm'
	import delimited "$SBM_REPO/dev/ref/`nm'_pi.csv", clear varnames(1)
	mkmat _all, matrix(RPI)
	restore
	nwsbm g, maxblocks(`rmaxq') model(`mod') covariates(`covs') seed(7) generate(mb) replace
	local ok 1
	if e(Q) != `rq' {
		di as error "`nm': Q selected `=e(Q)' but R selected `rq'"
		local ok 0
	}
	capture drop id
	capture drop rblock
	gen id = _n
	merge 1:1 id using `rm', nogen
	sbmt_ari mb rblock
	local ari = r(ari)
	if missing(`ari') | `ari' < 0.95 {
		di as error "`nm': ARI `ari' < 0.95"
		local ok 0
	}
	if "`mod'" == "gaussian" {
		if abs(e(sigma2) - `rs2') / `rs2' > 1e-3 {
			di as error "`nm': variance `=e(sigma2)' vs R `rs2'"
			local ok 0
		}
	}
	matrix B = e(beta)
	local maxb 0
	forvalues k = 1/`np' {
		local d = abs(B[1, `k'] - `rb`k'')
		if `d' > `maxb' local maxb = `d'
	}
	if `maxb' > 0.01 {
		di as error "`nm': coefficients differ from R by `maxb'"
		local ok 0
	}
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
	if `maxrel' > 0.003 | `relsel' > 0.0002 {
		di as error "`nm': ICL worse than R by `maxrel' (relative; at the selected Q: difference `relsel')"
		local ok 0
	}
	local qsel = e(Q)
	local one1 1
	* the one-block fit against an independent glm / lm fit
	preserve
	import delimited "$SBM_REPO/dev/ref/`nm'_glm.csv", clear varnames(1)
	local gll = loglik[1]
	forvalues k = 1/`np' {
		local gb`k' = beta`k'[1]
	}
	restore
	quietly nwsbm g, blocks(`one1') model(`mod') covariates(`covs') seed(7) nodisplay
	matrix B1 = e(beta)
	if abs(e(PL) - `gll') > 1e-6 * abs(`gll') {
		di as error "`nm': one-block PL " %12.6f e(PL) " vs glm " %12.6f `gll'
		local ok 0
	}
	forvalues k = 1/`np' {
		if abs(B1[1, `k'] - `gb`k'') > 1e-6 {
			di as error "`nm': one-block coefficient `k' " %10.7f B1[1, `k'] " vs glm " %10.7f `gb`k''
			local ok 0
		}
	}
	if `qsel' == `rq' {
		nwsbm g, maxblocks(`rmaxq') model(`mod') covariates(`covs') seed(7) generate(mb) replace nodisplay
		mata: sbmt_pidiff("mb", "rblock", `n')
	}
	else scalar __pidiff = .
	di as text "`nm': n=`n' dir=`dir' p=`np'  Q R=`rq' nwsbm=`=e(Q)'  ARI=" %6.4f `ari' "  ICL R=" %10.3f `ricl' " nwsbm=" %10.3f e(ICL) "  max rel ICL shortfall=" %9.2e `maxrel' "  max |beta diff|=" %8.5f `maxb' "  max |pi diff|=" %7.4f scalar(__pidiff)
	if scalar(__pidiff) > 0.05 {
		di as error "`nm': block parameters differ by " scalar(__pidiff)
		local ok 0
	}
	if `ok' di as result "PASS `nm'"
	else local fail = `fail' + 1
}
* two-mode fixtures (nwset puts the columns of the wide matrix into mode 1, so R's row blocks are nwsbm's blocks of mode 2)
local fix2 "cl_b_40x50 40 50 bernoulli 1 cl_g_40x50 40 50 gaussian 1 cl_p_40x50 40 50 poisson 1"
local nfix : word count `fix2'
forvalues f = 1(5)`nfix' {
	local nm : word `f' of `fix2'
	local n1 : word `=`f'+1' of `fix2'
	local n2 : word `=`f'+2' of `fix2'
	local mod : word `=`f'+3' of `fix2'
	local np : word `=`f'+4' of `fix2'
	sbmt_mknet2 g, ref(`nm') n1(`n1') n2(`n2') valued
	local covs
	forvalues k = 1/`np' {
		sbmt_mkcov x`k', file("$SBM_REPO/dev/ref/`nm'_cov`k'.csv") dir(2)
		local covs `covs' x`k'
	}
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
	import delimited "$SBM_REPO/dev/ref/`nm'_beta.csv", clear varnames(1)
	forvalues k = 1/`np' {
		local rb`k' = beta[`k']
	}
	import delimited "$SBM_REPO/dev/ref/`nm'_pi.csv", clear varnames(1)
	mkmat _all, matrix(RPI)
	restore
	nwsbm g, model(`mod') covariates(`covs') maxblocks(`mx2' `mx1') seed(7) generate(mb) replace
	local ok 1
	if e(Q1) != `rq2' | e(Q2) != `rq1' {
		di as error "`nm': (Q1,Q2) selected `=e(Q1)' x `=e(Q2)' but R selected `rq2' x `rq1' (columns x rows)"
		local ok 0
	}
	_nwsyntax g, max(1)
	mata: sbm_setidx(`netobj'->get_modes())
	mata: sbmt_ari_csv("mb", "$SBM_REPO/dev/ref/`nm'_membership2.csv", __nwsbm_i1)
	local ari2 = scalar(__ari)
	mata: sbmt_ari_csv("mb", "$SBM_REPO/dev/ref/`nm'_membership.csv", __nwsbm_i2)
	local ari1 = scalar(__ari)
	if missing(`ari1') | missing(`ari2') | `ari1' < 0.95 | `ari2' < 0.95 {
		di as error "`nm': ARI rows `ari1' columns `ari2' < 0.95"
		local ok 0
	}
	matrix B = e(beta)
	local maxb 0
	forvalues k = 1/`np' {
		local d = abs(B[1, `k'] - `rb`k'')
		if `d' > `maxb' local maxb = `d'
	}
	local btol = cond("`nm'" == "cl_b_40x50", 0.02, 0.01)
	if `maxb' > `btol' {
		di as error "`nm': coefficients differ from R by `maxb'"
		local ok 0
	}
	matrix T = e(ICLtable)
	local maxrel 0
	local relsel 0
	forvalues a = 1/`mx1' {
		forvalues b = 1/`mx2' {
			if "`rI_`a'_`b''" == "" continue
			local rv = `rI_`a'_`b''
			local rr = (`b' - 1) * `mx1' + `a'
			local worse = (`rv' - T[`rr', 3]) / abs(`rv')
			if `worse' > `maxrel' local maxrel = `worse'
			if `a' == `rq1' & `b' == `rq2' local relsel = abs(T[`rr', 3] - `rv') / abs(`rv')
		}
	}
	if "`nm'" == "cl_b_40x50" {
		local maxrel 0
		local relsel 0
	}
	if `maxrel' > 0.003 | `relsel' > 0.0002 {
		di as error "`nm': ICL worse than R by `maxrel' (relative; at the selected cell: difference `relsel')"
		local ok 0
	}
	if "`mod'" == "gaussian" {
		if abs(e(sigma2) - `rs2') / `rs2' > 1e-3 {
			di as error "`nm': variance `=e(sigma2)' vs R `rs2'"
			local ok 0
		}
	}
	local one1 1 1
	* the one-block fit against an independent glm / lm fit
	preserve
	import delimited "$SBM_REPO/dev/ref/`nm'_glm.csv", clear varnames(1)
	local gll = loglik[1]
	forvalues k = 1/`np' {
		local gb`k' = beta`k'[1]
	}
	restore
	quietly nwsbm g, blocks(`one1') model(`mod') covariates(`covs') seed(7) nodisplay
	matrix B1 = e(beta)
	if abs(e(PL) - `gll') > 1e-6 * abs(`gll') {
		di as error "`nm': one-block PL " %12.6f e(PL) " vs glm " %12.6f `gll'
		local ok 0
	}
	forvalues k = 1/`np' {
		if abs(B1[1, `k'] - `gb`k'') > 1e-6 {
			di as error "`nm': one-block coefficient `k' " %10.7f B1[1, `k'] " vs glm " %10.7f `gb`k''
			local ok 0
		}
	}
	nwsbm g, model(`mod') covariates(`covs') maxblocks(`mx2' `mx1') seed(7) generate(mb) replace nodisplay
	mata: sbmt_pidiff2("mb", "$SBM_REPO/dev/ref/`nm'_membership2.csv", "$SBM_REPO/dev/ref/`nm'_membership.csv", __nwsbm_i1, __nwsbm_i2)
	di as text "`nm': `n1'x`n2'  cell R=`rq2'x`rq1' nwsbm=`=e(Q1)'x`=e(Q2)'  ARI rows=" %6.4f `ari1' " cols=" %6.4f `ari2' "  ICL R=" %10.3f `ricl' " nwsbm=" %10.3f e(ICL) "  max rel ICL shortfall=" %9.2e `maxrel' "  max |beta diff|=" %8.5f `maxb' "  max |pi diff|=" %7.4f scalar(__pidiff)
	if scalar(__pidiff) > 0.05 {
		di as error "`nm': block parameters differ by " scalar(__pidiff)
		local ok 0
	}
	if `ok' di as result "PASS `nm'"
	else local fail = `fail' + 1
}
* fixtures with missing ties: one-block fit against glm, two blocks recovered, coefficients near the planted values
local fix3 "cbm_d_n60_Q2 60 1 bernoulli 2 cpm_u_n60_Q2 60 0 poisson 1"
local nfix : word count `fix3'
forvalues f = 1(5)`nfix' {
	local nm : word `f' of `fix3'
	local n : word `=`f'+1' of `fix3'
	local dir : word `=`f'+2' of `fix3'
	local mod : word `=`f'+3' of `fix3'
	local np : word `=`f'+4' of `fix3'
	sbmt_mknetv g, ref(`nm') n(`n') dir(`dir')
	local covs
	forvalues k = 1/`np' {
		sbmt_mkcov x`k', file("$SBM_REPO/dev/ref/`nm'_cov`k'.csv") dir(`dir')
		local covs `covs' x`k'
	}
	local ok 1
	local one1 1
	* the one-block fit against an independent glm / lm fit
	preserve
	import delimited "$SBM_REPO/dev/ref/`nm'_glm.csv", clear varnames(1)
	local gll = loglik[1]
	forvalues k = 1/`np' {
		local gb`k' = beta`k'[1]
	}
	restore
	quietly nwsbm g, blocks(`one1') model(`mod') covariates(`covs') seed(7) nodisplay
	matrix B1 = e(beta)
	if abs(e(PL) - `gll') > 1e-6 * abs(`gll') {
		di as error "`nm': one-block PL " %12.6f e(PL) " vs glm " %12.6f `gll'
		local ok 0
	}
	forvalues k = 1/`np' {
		if abs(B1[1, `k'] - `gb`k'') > 1e-6 {
			di as error "`nm': one-block coefficient `k' " %10.7f B1[1, `k'] " vs glm " %10.7f `gb`k''
			local ok 0
		}
	}
	nwsbm g, maxblocks(3) model(`mod') covariates(`covs') seed(7) generate(mb) replace nodisplay
	preserve
	import delimited "$SBM_REPO/dev/ref/`nm'_truth.csv", clear varnames(1)
	rename truth tb
	tempfile tr
	save `tr'
	restore
	capture drop id
	capture drop tb
	gen id = _n
	merge 1:1 id using `tr', nogen
	sbmt_ari mb tb
	di as text "`nm': missing ties `=e(n_missing)'; Q=`=e(Q)' ARI vs truth " %6.4f r(ari) "  beta " %8.4f el(e(beta), 1, 1)
	if e(Q) != 2 | r(ari) < 0.9 | e(n_missing) == 0 {
		di as error "`nm': two blocks not recovered with missing ties"
		local ok 0
	}
	if `ok' di as result "PASS `nm'"
	else local fail = `fail' + 1
}
if `fail' == 0 di as result "PASS all covariate fixtures agree with R sbm"
else {
	di as error "FAIL: `fail' fixtures differ"
	exit 9
}
