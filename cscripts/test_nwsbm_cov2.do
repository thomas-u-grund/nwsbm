* Covariates: recovery of planted partitions and coefficients (simulated), node-variable options against the equivalent covariate networks, error paths.
do cscripts/nwsbm_testlib.do
local fail 0
* ---- recovery: Q = 3 planted blocks, two covariates (|x_i - x_j| and a same-group indicator) ----
* model n dir pmiss beta1 beta2
local cases "bernoulli 150 0 0 -1.5 0.8 poisson 120 0 0 -1.0 0.5 gaussian 120 0 0 -1.5 0.8 poisson 120 1 0.05 -1.0 0.5 gaussian 120 1 0.05 -1.5 0.8 bernoulli 120 1 0 -1.5 0.8"
local nc : word count `cases'
forvalues f = 1(6)`nc' {
	local mod : word `f' of `cases'
	local n : word `=`f'+1' of `cases'
	local dir : word `=`f'+2' of `cases'
	local pm : word `=`f'+3' of `cases'
	local b1 : word `=`f'+4' of `cases'
	local b2 : word `=`f'+5' of `cases'
	if "`mod'" == "bernoulli" {
		mata: M = (-0.2,-2.5,-2.2 \ -2.5,-0.4,-2.3 \ -2.2,-2.3,-0.1)
		local m 0
	}
	else if "`mod'" == "poisson" {
		mata: M = (1.0,-0.8,-0.5 \ -0.8,0.8,-0.6 \ -0.5,-0.6,1.1)
		local m 1
	}
	else {
		mata: M = (2.0,0,0.3 \ 0,1.5,0 \ 0.3,0,2.5)
		local m 2
	}
	mata: sbmt_simcov(`n', M, (0.4*`n' \ 0.3*`n' \ 0.3*`n'), (`b1', `b2'), `m', `dir', 11, 1, `pm')
	sbmt_frommata g, mat(__sim_A) dir(`dir')
	sbmt_frommata x1, mat(__sim_X1) dir(0)
	sbmt_frommata x2, mat(__sim_X2) dir(0)
	* the covariate networks are undirected; in a directed network they are used as they are
	nwsbm g, maxblocks(5) model(`mod') covariates(x1 x2) seed(5) generate(mb) replace nodisplay
	capture drop tz
	mata: st_addvar("double", "tz")
	mata: st_store((1::`n'), "tz", __sim_z)
	sbmt_ari mb tz
	local ari = r(ari)
	matrix B = e(beta)
	local e1 = abs(B[1, 1] - `b1')
	local e2 = abs(B[1, 2] - `b2')
	di as text "recovery `mod' n=`n' dir=`dir' miss=`pm': Q=`=e(Q)' ARI=" %6.4f `ari' " beta " %7.3f B[1,1] " (" `b1' ") " %7.3f B[1,2] " (" `b2' ")"
	if e(Q) != 3 | `ari' < 0.9 | `e1' > 0.45 | `e2' > 0.45 {
		di as error "recovery failed: `mod' n=`n' dir=`dir'"
		local fail = `fail' + 1
	}
	else di as result "PASS recovery `mod' n=`n' dir=`dir' miss=`pm'"
}
* ---- node variables give the same fit as the covariate networks built from them ----
mata: sbmt_simcov(100, (-0.2,-2.5,-2.2 \ -2.5,-0.4,-2.3 \ -2.2,-2.3,-0.1), (40 \ 30 \ 30), (-1.5, 0.8), 0, 0, 12, 1, 0)
sbmt_frommata g, mat(__sim_A) dir(0)
sbmt_frommata x1, mat(__sim_X1) dir(0)
sbmt_frommata x2, mat(__sim_X2) dir(0)
capture drop age grp
mata: st_addvar("double", ("age", "grp"))
mata: st_store((1::100), ("age", "grp"), (__sim_x, __sim_g))
nwsbm g, blocks(3) covariates(x1 x2) seed(4) nodisplay
local icl1 = e(ICL)
matrix B1 = e(beta)
nwsbm g, blocks(3) absdiff(age) same(grp) seed(4) nodisplay
matrix B2 = e(beta)
if abs(e(ICL) - `icl1') > 1e-8 | abs(B1[1, 1] - B2[1, 1]) > 1e-8 | abs(B1[1, 2] - B2[1, 2]) > 1e-8 | "`: colnames B2'" != "abs_age same_grp" {
	di as error "absdiff()/same() differ from the equivalent covariate networks"
	local fail = `fail' + 1
}
else di as result "PASS node variables equal covariate networks"
* string variable for same()
gen str3 sgrp = "a" + string(grp)
nwsbm g, blocks(3) same(sgrp) seed(4) nodisplay
nwsbm g, blocks(3) same(grp) seed(4) nodisplay
local iclg = e(ICL)
nwsbm g, blocks(3) same(sgrp) seed(4) nodisplay
if abs(e(ICL) - `iclg') > 1e-8 {
	di as error "same() with a string variable differs from the numeric one"
	local fail = `fail' + 1
}
else di as result "PASS same() with a string variable"
* directed: signeddiff
mata: sbmt_simcov(80, (-0.2,-2.5 \ -2.5,-0.4), (40 \ 40), (-1.5), 0, 1, 13, 1, 0)
sbmt_frommata gd, mat(__sim_A) dir(1)
capture drop age
gen double age = .
mata: st_store((1::80), "age", __sim_x[1::80])
nwsbm gd, blocks(2) signeddiff(age) seed(4) nodisplay
if missing(el(e(beta), 1, 1)) | "`: colnames e(beta)'" != "diff_age" {
	di as error "signeddiff() failed"
	local fail = `fail' + 1
}
else di as result "PASS signeddiff on a directed network"
* ---- error paths (all must stop with the stated return code) ----
mata: sbmt_simcov(60, (-0.2,-2.5 \ -2.5,-0.4), (30 \ 30), (-1.5), 0, 0, 14, 1, 0)
sbmt_frommata g6, mat(__sim_A) dir(0)
sbmt_frommata c6, mat(__sim_X1) dir(0)
mata: sbmt_simcov(50, (-0.2,-2.5 \ -2.5,-0.4), (25 \ 25), (-1.5), 0, 0, 15, 1, 0)
sbmt_frommata c50, mat(__sim_X1) dir(0)
mata: sbmt_simcov(60, (-0.2,-2.5 \ -2.5,-0.4), (30 \ 30), (-1.5), 0, 1, 16, 1, 0)
sbmt_frommata cdir, mat(__sim_A) dir(1)
capture drop v60 vmiss vconst vstr
mata: st_addvar("double", ("v60", "vmiss", "vconst"))
mata: st_store((1::60), ("v60", "vmiss", "vconst"), (runiform(60, 1), (.\ runiform(59, 1)), J(60, 1, 2)))
gen str2 vstr = "x"
local e1 "nwsbm g6, blocks(2) covariates(c50)"
local c1 198
local e2 "nwsbm g6, blocks(2) covariates(cdir)"
local c2 198
local e3 "nwsbm g6, blocks(2) absdiff(vmiss)"
local c3 416
local e4 "nwsbm g6, blocks(2) same(vconst)"
local c4 198
local e5 "nwsbm g6, blocks(2) model(poisson) degreecorrected covariates(c6)"
local c5 198
local e6 "nwsbm g6, blocks(2) icl(exact) covariates(c6)"
local c6 198
local e7 "nwsbm g6, blocks(2) signeddiff(v60)"
local c7 198
local e8 "nwsbm g6, blocks(2) covariates(nosuchnet)"
local c8 0
local e9 "nwsbm g6, blocks(2) absdiff(vstr)"
local c9 198
local e10 "nwsbm g6, blocks(2) absdiff(nosuchvar)"
local c10 111
forvalues k = 1/10 {
	capture `e`k''
	local rc = _rc
	if (`c`k'' > 0 & `rc' != `c`k'') | (`c`k'' == 0 & `rc' == 0) {
		di as error "error path `k' (`e`k''): return code `rc', expected `c`k''"
		local fail = `fail' + 1
	}
	else di as result "PASS error path `k' (r(`rc'))"
}
if `fail' == 0 di as result "PASS covariate recovery, equivalence and error tests"
else {
	di as error "FAIL: `fail' checks"
	exit 9
}
