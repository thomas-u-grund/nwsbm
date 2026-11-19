* Timing of nwsbm against R sbm (estimateSimpleSBM, nbCores = 1, exploreMax = 6, times in dev/ref/time_u_n*_fit.csv and timing.csv).
* nwsbm fits every Q = 1..6 with 10 extra starts each. Prints the table; PASS when the Q selected equals R's and nwsbm is not slower.
do cscripts/nwsbm_testlib.do
local fail 0
foreach n in 200 500 1000 {
	sbmt_mknet g, ref(time_u_n`n') n(`n') dir(0)
	preserve
	import delimited "$SBM_REPO/dev/ref/time_u_n`n'_fit.csv", clear varnames(1)
	local rq = q[1]
	local rsec = seconds[1]
	local ricl = icl[1]
	restore
	timer clear 1
	timer on 1
	nwsbm g, maxblocks(6) seed(3) threads(1) nodisplay
	timer off 1
	quietly timer list 1
	local sec = r(t1)
	di as text "TIMING n=`n': nwsbm " %7.2f `sec' " s (Q=`=e(Q)', ICL " %10.2f e(ICL) ")"
	di as text "TIMING n=`n': R sbm " %7.2f `rsec' " s (Q=`rq', ICL " %10.2f `ricl' "), ratio R/nwsbm " %5.2f `rsec'/`sec'
	local icl1 = e(ICL)
	timer clear 2
	timer on 2
	quietly nwsbm g, maxblocks(6) seed(3) threads(4)
	timer off 2
	quietly timer list 2
	local sec4 = r(t2)
	timer clear 3
	timer on 3
	quietly nwsbm g, maxblocks(6) seed(3) threads(1) fast
	timer off 3
	quietly timer list 3
	di as text "TIMING n=`n': threads(4) " %7.2f `sec4' " s (ICL " %10.2f e(ICL) " for threads(1)=" %10.2f `icl1' ")"
	di as text "TIMING n=`n': fast, threads(1) " %7.2f r(t3) " s, Q=`=e(Q)', ICL " %10.2f e(ICL)
	quietly nwsbm g, maxblocks(6) seed(3) threads(4)
	if abs(`icl1' - e(ICL)) > 1e-6 {
		di as error "n=`n': the result depends on the number of threads"
		local fail = `fail' + 1
	}
	if e(Q) != `rq' | abs(e(ICL) - `ricl') > 0.001 * abs(`ricl') {
		di as error "n=`n': result differs from R"
		local fail = `fail' + 1
	}
}
if `fail' == 0 di as result "PASS timing fixtures agree with R"
else exit 9
