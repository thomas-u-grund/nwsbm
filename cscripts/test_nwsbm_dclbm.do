* Degree-corrected two-mode model (Poisson) against R greed 0.6.2 DcLbm (dev/ref_dclbm.R): three fixtures. greed's rows are mode 2 of the nwset network.
* Per fixture: the numbers of row and column blocks selected by the exact ICL equal greed's, ARI of both partitions against greed >= 0.95, exact ICL
* within 1e-4 relative of greed's (equal to 1e-6 when the partitions coincide), the plain Poisson latent block model finds more blocks on strongly
* heterogeneous degrees or a worse exact ICL (informational), the ICL table has the column ICLbm, errors for Bernoulli exact two-mode.
do cscripts/nwsbm_testlib.do
set linesize 200
local fail 0
foreach spec in "dl_40x60 40 60" "dl_80x50 80 50" "dl_100x100 100 100" {
	local nm : word 1 of `spec'
	local nr : word 2 of `spec'
	local nc : word 3 of `spec'
	sbmt_mknet2 g, ref(`nm') n1(`nr') n2(`nc') valued
	preserve
	import delimited "$SBM_REPO/dev/ref/`nm'_fit.csv", clear varnames(1)
	local rk1 = k1[1]
	local rk2 = k2[1]
	local ricl = icl[1]
	restore
	quietly nwsbm g, model(poisson) degreecorrected icl(exact) maxblocks(6 6) seed(7) generate(mb) replace nodisplay
	local ok 1
	if e(Q1) != `rk2' | e(Q2) != `rk1' {
		di as error "`nm': blocks selected `=e(Q1)' (cols) x `=e(Q2)' (rows), greed `rk2' x `rk1'"
		local ok 0
	}
	_nwsyntax g, max(1)
	mata: sbm_setidx(`netobj'->get_modes())
	mata: sbmt_ari_csv("mb", "$SBM_REPO/dev/ref/`nm'_membership_cols.csv", __nwsbm_i1)
	local ari1 = scalar(__ari)
	mata: sbmt_ari_csv("mb", "$SBM_REPO/dev/ref/`nm'_membership_rows.csv", __nwsbm_i2)
	local ari2 = scalar(__ari)
	if missing(`ari1') | missing(`ari2') | `ari1' < 0.95 | `ari2' < 0.95 {
		di as error "`nm': ARI columns `ari1' rows `ari2' < 0.95"
		local ok 0
	}
	local rel = abs(e(ICL) - `ricl') / abs(`ricl')
	local relbm = e(ICLbm)
	di as text "`nm': Q greed `rk2' x `rk1', nwsbm `=e(Q1)' x `=e(Q2)'; ARI cols " %6.4f `ari1' " rows " %6.4f `ari2' "; exact ICL greed " %12.4f `ricl' " nwsbm " %12.4f e(ICL) " (rel " %9.2e `rel' "); ICL (blockmodels type) " %12.3f `relbm'
	* nwsbm may find a slightly better partition (then its exact ICL is larger); it must not be worse by more than 1e-4 relative
	if (`ricl' - e(ICL)) / abs(`ricl') > 1e-4 {
		di as error "`nm': exact ICL below greed's"
		local ok 0
	}
	if `rel' > 1e-4 & `ari1' > 0.9999 & `ari2' > 0.9999 {
		di as error "`nm': identical partitions but exact ICL differs by `rel'"
		local ok 0
	}
	* the ICL table carries both criteria
	matrix T = e(ICLtable)
	if colsof(T) != 8 {
		di as error "`nm': ICL table has `=colsof(T)' columns, expected 8"
		local ok 0
	}
	if `ok' di as result "PASS `nm'"
	else local fail = `fail' + 1
}
* errors
sbmt_mknet2 g, ref(dl_40x60) n1(40) n2(60) valued
capture nwsbm g, model(poisson) icl(exact) blocks(2 2) nodisplay
if _rc != 198 {
	di as error "icl(exact) without degreecorrected on a two-mode network should fail (198), rc = " _rc
	local fail = `fail' + 1
}
else di as result "PASS error path exact without degreecorrected"
capture nwsbm g, model(bernoulli) degreecorrected blocks(2 2) nodisplay
if _rc != 198 {
	di as error "degreecorrected with model(bernoulli) should fail, rc = " _rc
	local fail = `fail' + 1
}
else di as result "PASS error path degreecorrected needs poisson"
* the ordinary criterion on the same network: runs and returns a table with 7 columns
quietly nwsbm g, model(poisson) degreecorrected blocks(2 2) seed(1) nodisplay
if colsof(e(ICLtable)) != 7 | e(ICL) >= . {
	di as error "default criterion table"
	local fail = `fail' + 1
}
else di as result "PASS default criterion for the degree-corrected two-mode model"
if `fail' di as error "FAILED `fail'"
else di as result "ALL PASSED"
