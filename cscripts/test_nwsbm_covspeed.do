* Speed and agreement of the Bernoulli covariate model (E and M steps with tabulated softplus). Prints times and the ICL of every Q; run once with the
* plugin of a reference build (NWSBM_PLUGDIR) and once with the current one to compare. PASS when the ICL table agrees with the values of the build
* that computed every softplus with log1p (the version before the tabulation; times there: 5.6 s for n = 100, 18.2 s for n = 200 on the development
* machine, one thread) to a relative 1e-7.
do cscripts/nwsbm_testlib.do
set linesize 200
local fail 0
local tag : env NWSBM_TAG
foreach n in 100 200 {
	mata: Mm = (-0.2,-2.5,-2.2 \ -2.5,-0.4,-2.3 \ -2.2,-2.3,-0.1)
	mata: sbmt_simcov(`n', Mm, (0.4*`n' \ 0.3*`n' \ 0.3*`n'), (-1.0, 0.5), 0, 0, 5, 1, 0)
	sbmt_frommata g, mat(__sim_A) dir(0)
	sbmt_frommata x1, mat(__sim_X1) dir(0)
	sbmt_frommata x2, mat(__sim_X2) dir(0)
	timer clear 5
	timer on 5
	quietly nwsbm g, maxblocks(5) covariates(x1 x2) seed(3) threads(1) nodisplay
	timer off 5
	quietly timer list 5
	matrix T = e(ICLtable)
	matrix B = e(beta)
	di as text "TIMING Bernoulli covariates n=`n' one thread: " %7.2f r(t5) " s, Q=`=e(Q)', beta " %9.6f B[1, 1] " " %9.6f B[1, 2]
	if `n' == 100 local ref "-2377.95774863 -2230.99942229 -2177.59364622 -2193.84946320 -2217.12880113"
	else local ref "-9548.25662943 -8734.82754904 -8280.53650824 -8304.98690362 -8334.64052416"
	forvalues r = 1/5 {
		local rv : word `r' of `ref'
		di as text "ICL n=`n' Q=`r' " %16.8f T[`r', 2] "  reference " %16.8f `rv'
		if abs(T[`r', 2] - `rv') > 1e-7 * abs(`rv') {
			di as error "n=`n' Q=`r': ICL differs from the reference"
			local fail = `fail' + 1
		}
	}
}
if `fail' == 0 di as result "PASS covspeed ICL equals the reference"
else exit 9
