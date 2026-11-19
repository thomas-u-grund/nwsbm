* Greedy search on the exact ICL (option greedy with icl(exact)): swaps of single nodes and merges, as the hybrid algorithm of R greed.
* Fixtures of dev/ref_dc.R and dev/ref_greedy.R (binary SBM, degree-corrected Poisson SBM; the gw_ ones have weak block structure) and dev/ref_dclbm.R (degree-corrected two-mode) with the exact ICL printed by greed.
* Checks: (1) the exact ICL with greedy is never smaller than without it (the search starts from the VEM partition) for every number of blocks of the table;
* (2) at the selected number of blocks it is within 1e-7 (relative) of greed's or larger; for dc_u_n200_Q4, where the VEM partition lies two nodes from greed's
* optimum, it reaches greed's value (the plain exact ICL does not); (3) the partition is hard (every membership probability 0 or 1); (4) greedy without
* icl(exact) is refused; (5) blocks() with greedy works for a fixed number of blocks.
do cscripts/nwsbm_testlib.do
set linesize 200
local fail 0
local fix "gr_sbm_u_n100_Q3 100 0 bernoulli gr_sbm_d_n100_Q2 100 1 bernoulli dc_u_n150_Q3 150 0 dc dc_d_n120_Q2 120 1 dc dc_u_n200_Q4 200 0 dc gw_sbm_u_n200_Q3 200 0 bernoulli gw_dc_u_n150_Q3 150 0 dc gw_sbm_d_n120_Q3 120 1 bernoulli"
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
	restore
	quietly nwsbm g, `opt' maxblocks(7) icl(exact) seed(5) nodisplay
	matrix T0 = e(ICLtable)
	local icl0 = e(ICL)
	timer clear 1
	timer on 1
	quietly nwsbm g, `opt' maxblocks(7) icl(exact) greedy seed(5) generate(mb) replace nodisplay
	timer off 1
	quietly timer list 1
	local tsec = r(t1)
	matrix T1 = e(ICLtable)
	local ok 1
	local worst 0
	forvalues r = 1/7 {
		local d = T1[`r', 2] - T0[`r', 2]
		if `d' < `worst' local worst = `d'
	}
	if `worst' < -1e-6 {
		di as error "`nm': greedy exact ICL below the plain one by `worst'"
		local ok 0
	}
	* hard memberships
	quietly summarize mb_p1 in 1/`n'
	local hard 1
	forvalues k = 1/`=e(Q)' {
		quietly count if mb_p`k' > 1e-9 & mb_p`k' < 1 - 1e-9 in 1/`n'
		if r(N) > 0 local hard 0
	}
	if !`hard' {
		di as error "`nm': greedy memberships are not hard"
		local ok 0
	}
	* selected number of blocks: the plain search finds greed's number; the greedy search must too and its exact ICL reach greed's
	if e(Q) != `rq' {
		di as error "`nm': greedy selects `=e(Q)' blocks, greed `rq'"
		local ok 0
	}
	local gap = (`ricl' - e(ICL)) / abs(`ricl')
	local gap0 = (`ricl' - `icl0') / abs(`ricl')
	di as text "`nm': exact ICL greed " %12.4f `ricl' "  nwsbm plain " %12.4f `icl0' " (rel. gap " %9.2e `gap0' ")  nwsbm greedy " %12.4f e(ICL) " (rel. gap " %9.2e `gap' ")  time " %5.2f `tsec' " s"
	if `gap' > 1e-7 {
		di as error "`nm': greedy exact ICL below greed's"
		local ok 0
	}
	if `ok' di as result "PASS `nm'"
	else local fail = `fail' + 1
}
* two-mode degree-corrected fixtures
foreach spec in "dl_40x60 40 60" "dl_80x50 80 50" "dl_100x100 100 100" {
	local nm : word 1 of `spec'
	local nr : word 2 of `spec'
	local nc : word 3 of `spec'
	sbmt_mknet2 g, ref(`nm') n1(`nr') n2(`nc') valued
	preserve
	import delimited "$SBM_REPO/dev/ref/`nm'_fit.csv", clear varnames(1)
	local ricl = icl[1]
	restore
	quietly nwsbm g, model(poisson) degreecorrected icl(exact) maxblocks(5 5) seed(7) nodisplay
	local icl0 = e(ICL)
	quietly nwsbm g, model(poisson) degreecorrected icl(exact) greedy maxblocks(5 5) seed(7) nodisplay
	local gap = (`ricl' - e(ICL)) / abs(`ricl')
	di as text "`nm': exact ICL greed " %12.4f `ricl' "  plain " %12.4f `icl0' "  greedy " %12.4f e(ICL) " (rel. gap " %9.2e `gap' ")"
	if `gap' > 1e-7 | e(ICL) < `icl0' - 1e-6 {
		di as error "`nm': greedy two-mode exact ICL"
		local fail = `fail' + 1
	}
	else di as result "PASS `nm' (two-mode)"
}
* errors
sbmt_mknet g, ref(gr_sbm_u_n100_Q3) n(100) dir(0)
capture nwsbm g, greedy blocks(3)
if _rc != 198 {
	di as error "greedy without icl(exact) should fail (198), rc = " _rc
	local fail = `fail' + 1
}
else di as result "PASS error path greedy without icl(exact)"
capture nwsbm g, greedy icl(exact) fast blocks(3)
if _rc != 198 {
	di as error "greedy with fast should fail (198), rc = " _rc
	local fail = `fail' + 1
}
else di as result "PASS error path greedy with fast"
quietly nwsbm g, greedy icl(exact) blocks(3) nodisplay seed(2)
if e(Q) != 3 | "`e(greedy)'" != "yes" {
	di as error "blocks(3) with greedy"
	local fail = `fail' + 1
}
else di as result "PASS blocks(3) with greedy"
if `fail' di as error "FAILED `fail'"
else di as result "ALL PASSED"
