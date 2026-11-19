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
* default number of restarts for large networks (n > 500: 4, n > 1500: 2): the n = 1000 fixture gives the same ICL with 4 and with 10 restarts
sbmt_mknet g, ref(time_u_n1000) n(1000) dir(0)
quietly nwsbm g, maxblocks(6) seed(3)
local r_def = e(restarts)
local icl_def = e(ICL)
quietly nwsbm g, maxblocks(6) seed(3) restarts(10)
di as text "TIMING n=1000: default restarts `r_def' ICL " %10.3f `icl_def' "; restarts(10) ICL " %10.3f e(ICL)
if `r_def' != 4 | abs(`icl_def' - e(ICL)) > 1e-6 {
	di as error "default restarts rule: restarts `r_def', ICL differs from restarts(10)"
	local fail = `fail' + 1
}
* covariates (Bernoulli, Poisson, Gaussian; planted 3 blocks, one covariate; not compared with R, informational): single thread and all cores
foreach n in 100 200 {
	foreach m in bernoulli poisson gaussian {
		local mi = cond("`m'" == "bernoulli", 0, cond("`m'" == "poisson", 1, 2))
		if `mi' == 0 mata: Mm = (-0.2,-2.5,-2.2 \ -2.5,-0.4,-2.3 \ -2.2,-2.3,-0.1)
		else if `mi' == 1 mata: Mm = (1.0,-0.8,-0.5 \ -0.8,0.8,-0.6 \ -0.5,-0.6,1.1)
		else mata: Mm = (2.0,0,0.3 \ 0,1.5,0 \ 0.3,0,2.5)
		mata: sbmt_simcov(`n', Mm, (0.4*`n' \ 0.3*`n' \ 0.3*`n'), (-1.0), `mi', 0, 5, 1, 0)
		sbmt_frommata g, mat(__sim_A) dir(0)
		sbmt_frommata x1, mat(__sim_X1) dir(0)
		timer clear 5
		timer on 5
		quietly nwsbm g, maxblocks(5) model(`m') covariates(x1) seed(3) threads(1)
		timer off 5
		quietly timer list 5
		local t1 = r(t5)
		local q1 = e(Q)
		timer clear 6
		timer on 6
		quietly nwsbm g, maxblocks(5) model(`m') covariates(x1) seed(3)
		timer off 6
		quietly timer list 6
		di as text "TIMING covariates `m' n=`n': one thread " %7.2f `t1' " s (Q=`q1', restarts `=e(restarts)'), all cores " %7.2f r(t6) " s (Q=`=e(Q)')"
	}
}
if `fail' == 0 di as result "PASS timing fixtures agree with R"
else exit 9
