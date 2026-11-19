*! version 0.2.0  19nov2026 nwsbm
* Copyright (C) 2026 Thomas U. Grund
* This file is part of nwsbm, a Stata package for stochastic block models. It follows the variational EM and the ICL of the R packages
* blockmodels and sbm and is written independently.
* nwsbm is free software: you can redistribute it and/or modify it under the terms of the GNU General
* Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
* that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
* or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.

capture program drop nwsbm
program nwsbm, eclass
	version 14
	gettoken sub rest : 0, parse(" ,")
	if "`sub'" == "version" {
		nwsbm_version
		exit
	}
	if "`sub'" == "" | substr(trim(`"`0'"'), 1, 1) == "," {
		if "`e(cmd)'" != "nwsbm" error 301
		syntax [, GENerate(string) replace]
		if "`generate'" != "" _nwsbm_generate `generate', `replace'
		_nwsbm_display
		exit
	}
	syntax anything(name=netname) [, BLOCKs(numlist integer >=1 max=2) MAXBlocks(numlist integer >=0 max=2) MINBlocks(numlist integer >=1 max=2) ///
		MODel(string) DEGREEcorrected SIGma2(real -1) SEED(string) RESTarts(integer 10) THReads(integer 0) FAST GENerate(string) replace]
	local __nb : word count `blocks'
	local __nmx : word count `maxblocks'
	local __nmn : word count `minblocks'
	if "`model'" == "" local model "bernoulli"
	if !inlist("`model'", "bernoulli", "poisson", "gaussian") {
		di as error "model() must be bernoulli, poisson or gaussian"
		error 198
	}
	local __model = cond("`model'" == "bernoulli", 0, cond("`model'" == "poisson", 1, 2))
	if "`degreecorrected'" != "" & "`model'" != "poisson" {
		di as error "degreecorrected needs model(poisson)"
		error 198
	}
	if `sigma2' != -1 & ("`model'" != "gaussian" | `sigma2' <= 0) {
		di as error "sigma2() sets a fixed variance (a positive number) and needs model(gaussian)"
		error 198
	}
	if `sigma2' == -1 local sigma2 0
	if `__nb' >= 1 & `__nmx' >= 1 {
		di as error "blocks() and maxblocks() cannot be combined"
		error 198
	}
	if `restarts' < 0 | `restarts' > 1000 {
		di as error "restarts() must be between 0 and 1000"
		error 198
	}
	if `threads' < 0 | `threads' > 256 {
		di as error "threads() must be between 1 and 256 (default: one thread per core)"
		error 198
	}
	if `"`seed'"' != "" {
		capture confirm integer number `seed'
		if _rc | real(`"`seed'"') < 0 | real(`"`seed'"') > 2147483647 {
			di as error "seed() must be an integer between 0 and 2147483647"
			error 198
		}
	}
	_nwsbm_libcheck
	set more off
	_nwsyntax `netname', max(1)
	local __kind = ("`is2mode'" == "true")
	if `__kind' & "`degreecorrected'" != "" {
		di as error "degreecorrected is available for one-mode networks only"
		error 198
	}
	if `__kind' == 0 & (`__nb' == 2 | `__nmx' == 2 | `__nmn' == 2) {
		di as error "{p 0 0 2}blocks(#1 #2), maxblocks(#1 #2) and minblocks(#1 #2) (row and column blocks) are for two-mode networks; `netname' is one-mode.{p_end}"
		error 198
	}
	if `__kind' & (`__nb' == 1 | `__nmx' == 1 | `__nmn' == 1) {
		di as error "{p 0 0 2}`netname' is a two-mode network: give two numbers, one for the blocks of mode 1 and one for mode 2 (for example blocks(3 2)).{p_end}"
		error 198
	}
	if "`valued'" == "true" & "`model'" == "bernoulli" {
		di as error "{p 0 0 2}model(bernoulli) fits binary networks; `netname' is a valued network. Use model(poisson) (counts) or model(gaussian), or dichotomise it with {bf:nwreplace}.{p_end}"
		error 198
	}
	if `nodes' < 3 {
		di as error "the network needs at least 3 actors"
		error 198
	}
	if "`generate'" != "" {
		capture confirm name `generate'
		if _rc {
			di as error "generate(): `generate' is not a valid variable name stub"
			error 198
		}
	}
	local __dir = (`__kind' == 0 & "`directed'" == "true")
	* the data matrix
	capture mata: mata drop __nwsbm_y
	if `__kind' {
		mata: sbm_setidx(`netobj'->get_modes())
		mata: __nwsbm_y = (*(`netobj'->get_matrix()))[__nwsbm_i1, __nwsbm_i2]
		mata: st_local("__n1", strofreal(rows(__nwsbm_y)))
		mata: st_local("__n2", strofreal(cols(__nwsbm_y)))
	}
	else {
		mata: __nwsbm_y = sbm_ymat(*(`netobj'->get_matrix()), `__dir')
		local __n1 = `nodes'
		local __n2 = `nodes'
	}
	mata: st_local("__bad", strofreal(sbm_checkvalues(__nwsbm_y, `__model')))
	if `__bad' == 1 {
		capture mata: mata drop __nwsbm_y
		di as error "model(poisson) needs non-negative integer ties (counts)"
		error 198
	}
	if `__bad' == 2 {
		capture mata: mata drop __nwsbm_y
		di as error "model(bernoulli) needs ties that are 0 or 1"
		error 198
	}
	* ranges of block numbers
	local qmin = 1
	local q2min = 1
	if `__nb' >= 1 {
		local qmin = real(word("`blocks'", 1))
		local qmax = `qmin'
		if `__kind' {
			local q2min = real(word("`blocks'", 2))
			local q2max = `q2min'
		}
	}
	else {
		if `__nmn' >= 1 local qmin = real(word("`minblocks'", 1))
		if `__nmn' == 2 local q2min = real(word("`minblocks'", 2))
		local qmax = cond(`__nmx' >= 1, real(word("`maxblocks'", 1)), 0)
		if `qmax' == 0 local qmax = min(10, `__n1')
		if `__kind' {
			local q2max = cond(`__nmx' == 2, real(word("`maxblocks'", 2)), 0)
			if `q2max' == 0 local q2max = min(10, `__n2')
		}
	}
	if !`__kind' {
		local q2min = 0
		local q2max = 0
	}
	if `qmax' < `qmin' | (`__kind' & `q2max' < `q2min') {
		capture mata: mata drop __nwsbm_y
		di as error "maxblocks() must not be smaller than minblocks()"
		error 198
	}
	if `qmax' > `__n1' | `qmax' > 100 | (`__kind' & (`q2max' > `__n2' | `q2max' > 100)) {
		capture mata: mata drop __nwsbm_y
		di as error "{p 0 0 2}the number of blocks cannot exceed the number of actors (of the mode) or 100.{p_end}"
		error 198
	}
	if `"`seed'"' == "" local __seed = floor(runiform() * 2147483646) + 1
	else local __seed = `seed'
	local __gq = max(`qmax', cond(`__kind', `q2max', 0))
	if "`generate'" != "" & "`replace'" == "" {
		local __vars "`generate'"
		forvalues k = 1/`__gq' {
			local __vars `__vars' `generate'_p`k'
		}
		foreach __v of local __vars {
			capture confirm variable `__v', exact
			if !_rc {
				capture mata: mata drop __nwsbm_y
				di as error "variable {bf:`__v'} already exists; specify {bf:replace}"
				error 110
			}
		}
	}
	tempname TAB PI AL BS AL2 BS2
	local __fast = ("`fast'" != "")
	local __dc = ("`degreecorrected'" != "")
	capture noisily mata: SbmFit(__nwsbm_y, `__dir', `__model', `__kind', `__dc', `qmin', `qmax', `q2min', `q2max', `__seed', `restarts', `sigma2', `threads', `__fast', "`TAB'", "`PI'", "`AL'", "`BS'", "`AL2'", "`BS2'")
	local __rc = _rc
	mata: st_local("__nmiss", strofreal(sum(__nwsbm_y :>= .) / (1 + (!`__dir' & !`__kind'))))
	capture mata: mata drop __nwsbm_y
	if `__rc' exit `__rc'
	local __q = scalar(__nwsbm_Q)
	local __q2 = scalar(__nwsbm_Q2)
	local __sig = scalar(__nwsbm_sigma2)
	capture scalar drop __nwsbm_Q __nwsbm_Q2 __nwsbm_sigma2
	* row of the selected fit in the table
	local __ic = cond(`__kind', 3, 2)
	local __row = 1
	forvalues r = 1/`=rowsof(`TAB')' {
		if `TAB'[`r', 1] == `__q' & (`__kind' == 0 | `TAB'[`r', 2] == `__q2') local __row = `r'
	}
	ereturn clear
	ereturn scalar N = `nodes'
	ereturn scalar Q = `__q'
	if `__kind' {
		ereturn scalar Q1 = `__q'
		ereturn scalar Q2 = `__q2'
		ereturn scalar N1 = `__n1'
		ereturn scalar N2 = `__n2'
	}
	ereturn scalar ICL = `TAB'[`__row', `__ic']
	ereturn scalar J = `TAB'[`__row', `__ic' + 1]
	ereturn scalar PL = `TAB'[`__row', `__ic' + 2]
	ereturn scalar H = `TAB'[`__row', `__ic' + 3]
	ereturn scalar npar = `TAB'[`__row', `__ic' + 4]
	ereturn scalar directed = `__dir'
	ereturn scalar n_missing = `__nmiss'
	ereturn scalar seed = `__seed'
	ereturn scalar restarts = `restarts'
	if "`model'" == "gaussian" ereturn scalar sigma2 = `__sig'
	ereturn local cmd "nwsbm"
	ereturn local model "`model'"
	ereturn local network "`netname'"
	ereturn local degreecorrected = cond(`__dc', "yes", "no")
	ereturn local membership = cond(`__kind', "LBM", cond(`__dir', "SBM", "SBM_sym"))
	if `__kind' matrix colnames `TAB' = Q1 Q2 ICL J PL H npar
	else matrix colnames `TAB' = Q ICL J PL H npar
	forvalues k = 1/`__q' {
		local __bn `__bn' b`k'
	}
	if `__kind' {
		forvalues k = 1/`__q2' {
			local __bn2 `__bn2' c`k'
		}
		matrix rownames `PI' = `__bn'
		matrix colnames `PI' = `__bn2'
		matrix colnames `AL2' = `__bn2'
		matrix colnames `BS2' = `__bn2'
		ereturn matrix alpha2 = `AL2'
		ereturn matrix blocksize2 = `BS2'
	}
	else {
		matrix rownames `PI' = `__bn'
		matrix colnames `PI' = `__bn'
	}
	matrix colnames `AL' = `__bn'
	matrix colnames `BS' = `__bn'
	ereturn matrix pi = `PI'
	ereturn matrix alpha = `AL'
	ereturn matrix blocksize = `BS'
	ereturn matrix ICLtable = `TAB'
	if "`generate'" != "" _nwsbm_generate `generate', `replace'
	_nwsbm_display
end

capture program drop _nwsbm_generate
program _nwsbm_generate
	syntax name(name=stub), [replace]
	local q = e(Q)
	if "`e(membership)'" == "LBM" local q = max(e(Q1), e(Q2))
	local n = e(N)
	capture confirm variable `stub', exact
	if !_rc & "`replace'" == "" {
		di as error "variable {bf:`stub'} already exists; specify {bf:replace}"
		exit 110
	}
	* dropping the only variable of a dataset also removes its observations: drop first, then make sure there are n rows
	capture drop `stub'
	forvalues k = 1/`q' {
		capture drop `stub'_p`k'
	}
	qui if _N < `n' {
		set obs `n'
	}
	qui gen double `stub' = .
	label var `stub' "nwsbm block (MAP)"
	forvalues k = 1/`q' {
		qui gen double `stub'_p`k' = .
		label var `stub'_p`k' "nwsbm P(block `k')"
	}
	if "`e(membership)'" == "LBM" {
		capture mata: mata describe __nwsbm_i1
		if _rc {
			di as error "the network information is no longer available; refit before generate()"
			exit 498
		}
		mata: SbmStore("`stub'", 1, __nwsbm_i1, __nwsbm_i2)
	}
	else mata: SbmStore("`stub'", 1)
end

capture program drop _nwsbm_display
program _nwsbm_display
	local q = e(Q)
	local two = ("`e(membership)'" == "LBM")
	local mname "`e(model)'"
	if "`e(degreecorrected)'" == "yes" local mname "`mname', degree corrected"
	di
	di as text "Stochastic block model (" as result "`mname'" as text ", " as result cond(`two', "two-mode", cond(e(directed), "directed", "undirected")) as text ")" ///
		_col(60) "Actors" _col(70) "= " as result %9.0g e(N)
	di as text "Network: " as result "`e(network)'" as text _col(60) "Blocks" _col(70) "= " as result %9.0g `q' cond(`two', " x " + string(e(Q2)), "")
	di as text _col(60) "ICL" _col(70) "= " as result %9.2f e(ICL)
	di as text _col(60) "Lower bound J" _col(70) "= " as result %9.2f e(J)
	if "`e(model)'" == "gaussian" di as text _col(60) "Variance" _col(70) "= " as result %9.4g e(sigma2)
	tempname T
	matrix `T' = e(ICLtable)
	local c0 = cond(`two', 3, 2)
	if rowsof(`T') > 1 {
		di
		if `two' {
			di as text "{hline 13}{c TT}{hline 61}"
			di as text %12s "Blocks 1 x 2" " {c |}" %12s "ICL" %12s "J" %12s "PL" %12s "parameters"
			di as text "{hline 13}{c +}{hline 61}"
		}
		else {
			di as text "{hline 13}{c TT}{hline 48}"
			di as text %12s "Blocks" " {c |}" %12s "ICL" %12s "J" %12s "PL" %12s "parameters"
			di as text "{hline 13}{c +}{hline 48}"
		}
		forvalues r = 1/`=rowsof(`T')' {
			local sel = (`T'[`r', 1] == `q') & (`two' == 0 | `T'[`r', 2] == e(Q2))
			local mark = cond(`sel', "*", " ")
			local lab = cond(`two', string(`T'[`r', 1]) + " x " + string(`T'[`r', 2]), string(`T'[`r', 1]))
			di as text %11s "`lab'" "`mark' {c |}" as result %12.2f `T'[`r', `c0'] %12.2f `T'[`r', `c0' + 1] %12.2f `T'[`r', `c0' + 2] %12.0f `T'[`r', `c0' + 4]
		}
		di as text "{hline 13}{c BT}{hline 48}" cond(`two', "{hline 13}", "")
		di as text "* selected by the largest ICL"
	}
	di
	di as text "Block sizes" cond(`two', " (mode 1):", ":")
	matrix list e(blocksize), noheader format(%9.0f)
	if `two' {
		di as text "Block sizes (mode 2):"
		matrix list e(blocksize2), noheader format(%9.0f)
	}
	di as text "Block proportions" cond(`two', " (mode 1):", ":")
	matrix list e(alpha), noheader format(%9.4f)
	if `two' {
		di as text "Block proportions (mode 2):"
		matrix list e(alpha2), noheader format(%9.4f)
	}
	if "`e(model)'" == "bernoulli" local what "Block connectivity, tie probability"
	else if "`e(model)'" == "poisson" local what "Block connectivity, expected tie value"
	else local what "Block means, expected tie value"
	di as text "`what', from the row block to the column block:"
	matrix list e(pi), noheader format(%9.4f)
end
