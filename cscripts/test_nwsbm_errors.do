* Error paths of nwsbm (each must stop with the stated return code and must not leave a result behind) and a few result checks.
do cscripts/nwsbm_testlib.do
local fail 0
capture program drop sbmt_expect
program sbmt_expect
	args rc want label
	if `rc' != `want' {
		di as error "FAIL `label': r(`rc') instead of r(`want')"
		c_local fail = `c(fail)'
		global SBM_FAILS = ${SBM_FAILS} + 1
	}
	else di as text "ok   `label' (r(`rc'))"
end
global SBM_FAILS 0
sbmt_mknet g, ref(pp_u_n60_Q2) n(60) dir(0)

* replay without a fit
capture nwsbm
sbmt_expect `=_rc' 301 "replay without a fit"
capture nwsbm g, blocks(2 3)
sbmt_expect `=_rc' 198 "two block numbers (latent block model)"
capture nwsbm g, blocks(2) maxblocks(4)
sbmt_expect `=_rc' 198 "blocks() with maxblocks()"
capture nwsbm g, model(poisson) blocks(2)
sbmt_expect `=_rc' 0 "model(poisson) on a binary network (allowed, as R)"
capture nwsbm g, model(gaussian) blocks(2)
sbmt_expect `=_rc' 0 "model(gaussian) on a binary network (allowed, as R)"
capture nwsbm g, blocks(2) degreecorrected
sbmt_expect `=_rc' 198 "degreecorrected needs the Poisson model"
capture nwsbm g, model(poisson) blocks(2) sigma2(1)
sbmt_expect `=_rc' 198 "sigma2() needs the Gaussian model"
capture nwsbm g, model(gaussian) blocks(2) sigma2(-3)
sbmt_expect `=_rc' 198 "negative sigma2()"
capture nwsbm g, model(gaussian) blocks(2) icl(exact)
sbmt_expect `=_rc' 198 "icl(exact) with the Gaussian model"
capture nwsbm g, model(poisson) blocks(2) icl(exact)
sbmt_expect `=_rc' 198 "icl(exact) with the plain Poisson model"
capture nwsbm g, blocks(2) icl(foo)
sbmt_expect `=_rc' 198 "unknown icl()"
capture nwsbm g, blocks(2) threads(300)
sbmt_expect `=_rc' 198 "threads() out of range"
capture nwsbm g, blocks(2) maxblocks(3 3)
sbmt_expect `=_rc' 198 "two maxblocks numbers on a one-mode network"
capture nwsbm g, model(foo) blocks(2)
sbmt_expect `=_rc' 198 "unknown model"
capture nwsbm g, blocks(0)
sbmt_expect `=_rc' 125 "blocks(0)"
capture nwsbm g, blocks(61)
sbmt_expect `=_rc' 198 "more blocks than actors"
capture nwsbm g, minblocks(5) maxblocks(3)
sbmt_expect `=_rc' 198 "minblocks above maxblocks"
capture nwsbm g, blocks(2) seed(-1)
sbmt_expect `=_rc' 198 "negative seed"
capture nwsbm g, blocks(2) seed(abc)
sbmt_expect `=_rc' 198 "non-numeric seed"
capture nwsbm g, blocks(2) restarts(-1)
sbmt_expect `=_rc' 198 "negative restarts"
capture nwsbm nonexistent, blocks(2)
if _rc == 0 {
	di as error "FAIL: a missing network was accepted"
	global SBM_FAILS = ${SBM_FAILS} + 1
}
else di as text "ok   missing network (r(`=_rc'))"

* two-mode and valued networks
clear
set obs 6
gen a = _n
gen b = mod(_n, 3) + 1
nwset a b, twomode name(tm) replace
capture nwsbm tm, blocks(2)
sbmt_expect `=_rc' 198 "two-mode network needs two block numbers"
capture nwsbm tm, blocks(2 2) degreecorrected model(bernoulli)
sbmt_expect `=_rc' 198 "degree correction needs the Poisson model"
capture nwsbm tm, blocks(2 2) model(bernoulli) icl(exact)
sbmt_expect `=_rc' 198 "icl(exact) on a two-mode network without degree correction"
capture nwsbm tm, blocks(2 2) model(poisson) degreecorrected icl(exact) greedy fast
sbmt_expect `=_rc' 198 "greedy with fast"
capture nwsbm tm, blocks(2 2) model(poisson) blockvar
sbmt_expect `=_rc' 198 "blockvar needs the Gaussian model"
capture nwsbm tm, blocks(9 2)
sbmt_expect `=_rc' 198 "more blocks than actors of a mode"
capture nwsbm tm, blocks(2 2) seed(3)
sbmt_expect `=_rc' 0 "two-mode network, blocks(2 2)"
if e(Q1) != 2 | e(Q2) != 2 | "`e(membership)'" != "LBM" {
	di as error "FAIL: two-mode e() results"
	global SBM_FAILS = ${SBM_FAILS} + 1
}
clear
set obs 5
gen double v1 = 0
gen double v2 = 0
gen double v3 = 0
gen double v4 = 0
gen double v5 = 0
replace v2 = 3 in 1
replace v3 = 2 in 1
replace v1 = 3 in 2
replace v4 = 5 in 2
replace v1 = 2 in 3
replace v5 = 1 in 4
replace v4 = 1 in 5
nwset v1-v5, name(vn) undirected valued replace
capture nwsbm vn, blocks(2)
sbmt_expect `=_rc' 198 "valued network with model(bernoulli)"
capture nwsbm vn, blocks(2) model(poisson) seed(1)
sbmt_expect `=_rc' 0 "valued network with model(poisson)"
capture nwsbm vn, blocks(2) model(gaussian) sigma2(2) seed(1)
sbmt_expect `=_rc' 0 "fixed variance"
if abs(e(sigma2) - 2) > 1e-12 {
	di as error "FAIL: fixed variance not kept"
	global SBM_FAILS = ${SBM_FAILS} + 1
}
replace v2 = 0.5 in 1
replace v1 = 0.5 in 2
nwset v1-v5, name(vn2) undirected valued replace
capture nwsbm vn2, blocks(2) model(poisson)
sbmt_expect `=_rc' 198 "Poisson with non-integer ties"
replace v2 = -1 in 1
replace v1 = -1 in 2
nwset v1-v5, name(vn3) undirected valued replace
capture nwsbm vn3, blocks(2) model(poisson)
sbmt_expect `=_rc' 198 "Poisson with negative ties"
capture nwsbm vn3, blocks(2) model(gaussian) seed(1)
sbmt_expect `=_rc' 0 "Gaussian with negative ties"
clear
set obs 2
gen double v1 = 0
gen double v2 = 0
nwset v1 v2, name(tiny) undirected replace
capture nwsbm tiny, blocks(1)
sbmt_expect `=_rc' 198 "two actors"

* generate(): existing variable needs replace
sbmt_mknet g, ref(pp_u_n60_Q2) n(60) dir(0)
clear
set obs 60
gen double mb = 1
capture nwsbm g, blocks(2) generate(mb)
sbmt_expect `=_rc' 110 "generate() over an existing variable"
nwsbm g, blocks(2) generate(mb) replace seed(1)
capture confirm variable mb_p1 mb_p2
sbmt_expect `=_rc' 0 "generate() creates the probability variables"
quietly summarize mb
if r(min) != 1 | r(max) != 2 {
	di as error "FAIL: block variable not in 1..2"
	global SBM_FAILS = ${SBM_FAILS} + 1
}
gen double ps = mb_p1 + mb_p2
quietly summarize ps
if abs(r(min) - 1) > 1e-9 | abs(r(max) - 1) > 1e-9 {
	di as error "FAIL: probabilities do not sum to one"
	global SBM_FAILS = ${SBM_FAILS} + 1
}
else di as text "ok   probabilities sum to one"

* results and replay
nwsbm g, maxblocks(4) seed(5)
if e(cmd) != "nwsbm" | e(Q) != 2 | e(N) != 60 | rowsof(e(ICLtable)) != 4 | colsof(e(pi)) != 2 | abs(e(blocksize)[1,1] + e(blocksize)[1,2] - 60) > 0 {
	di as error "FAIL: e() results"
	global SBM_FAILS = ${SBM_FAILS} + 1
}
else di as text "ok   e() results"
nwsbm
nwsbm, generate(r2) replace
capture confirm variable r2 r2_p1 r2_p2
sbmt_expect `=_rc' 0 "replay with generate()"
* same seed, same result; one block gives density
nwsbm g, maxblocks(4) seed(5)
local i1 = e(ICL)
nwsbm g, maxblocks(4) seed(5)
if e(ICL) != `i1' {
	di as error "FAIL: same seed gives different ICL"
	global SBM_FAILS = ${SBM_FAILS} + 1
}
else di as text "ok   reproducible with seed()"
nwsbm g, blocks(1)
if e(Q) != 1 | e(alpha)[1,1] != 1 | e(pi)[1,1] <= 0 | e(pi)[1,1] >= 1 | e(npar) != 1 {
	di as error "FAIL: one block"
	global SBM_FAILS = ${SBM_FAILS} + 1
}
nwsbm version
if "`r(version)'" == "" {
	di as error "FAIL: nwsbm version"
	global SBM_FAILS = ${SBM_FAILS} + 1
}
if ${SBM_FAILS} == 0 di as result "PASS error paths and result checks"
else {
	di as error "FAIL: ${SBM_FAILS} checks"
	exit 9
}
