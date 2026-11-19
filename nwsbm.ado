*! version 0.1.0  19nov2026 nwsbm
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
	syntax anything(name=netname) [, BLOCKs(numlist integer >=1 max=2) MAXBlocks(integer 0) MINBlocks(integer 1) MODel(string) ///
		SEED(string) RESTarts(integer 10) GENerate(string) replace]
	local __nb : word count `blocks'
	if `__nb' == 2 {
		di as error "{p 0 0 2}blocks(#1 #2), the latent block model for two-mode networks, is not available in this version.{p_end}"
		error 198
	}
	if "`model'" == "" local model "bernoulli"
	if !inlist("`model'", "bernoulli", "poisson", "gaussian") {
		di as error "model() must be bernoulli, poisson or gaussian"
		error 198
	}
	if "`model'" != "bernoulli" {
		di as error "{p 0 0 2}model(`model') is not available in this version; only model(bernoulli) is.{p_end}"
		error 198
	}
	if `__nb' == 1 & `maxblocks' > 0 {
		di as error "blocks() and maxblocks() cannot be combined"
		error 198
	}
	if `restarts' < 0 | `restarts' > 1000 {
		di as error "restarts() must be between 0 and 1000"
		error 198
	}
	if `minblocks' < 1 {
		di as error "minblocks() must be at least 1"
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
	if "`is2mode'" == "true" {
		di as error "{p 0 0 2}nwsbm fits one-mode networks; `netname' is a two-mode network (the latent block model is not available in this version).{p_end}"
		error 198
	}
	if "`valued'" == "true" {
		di as error "{p 0 0 2}model(bernoulli) fits binary networks; `netname' is a valued network. Use {bf:nwset} or {bf:nwreplace} to dichotomise it first.{p_end}"
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
	local __dir = ("`directed'" == "true")
	if `__nb' == 1 {
		local qmin = real(word("`blocks'", 1))
		local qmax = `qmin'
	}
	else {
		local qmin = `minblocks'
		local qmax = `maxblocks'
		if `qmax' == 0 local qmax = min(10, `nodes')
	}
	if `qmax' < `qmin' {
		di as error "maxblocks() must not be smaller than minblocks()"
		error 198
	}
	if `qmax' > `nodes' | `qmax' > 100 {
		di as error "{p 0 0 2}the number of blocks cannot exceed the number of actors or 100.{p_end}"
		error 198
	}
	if `"`seed'"' == "" local __seed = floor(runiform() * 2147483646) + 1
	else local __seed = `seed'
	if "`generate'" != "" & "`replace'" == "" {
		local __vars "`generate'"
		forvalues k = 1/`qmax' {
			local __vars `__vars' `generate'_p`k'
		}
		foreach __v of local __vars {
			capture confirm variable `__v', exact
			if !_rc {
				di as error "variable {bf:`__v'} already exists; specify {bf:replace}"
				error 110
			}
		}
	}
	tempname TAB PI AL BS
	capture mata: mata drop __nwsbm_y
	mata: __nwsbm_y = sbm_ymat(*(`netobj'->get_matrix()), `__dir')
	capture noisily mata: st_local("__q", strofreal(SbmFit(__nwsbm_y, `__dir', `qmin', `qmax', `__seed', `restarts', "`TAB'", "`PI'", "`AL'", "`BS'")))
	local __rc = _rc
	mata: st_local("__nmiss", strofreal(sum(__nwsbm_y :>= .) / (1 + !`__dir')))
	capture mata: mata drop __nwsbm_y
	if `__rc' exit `__rc'
	local __row = `__q' - `qmin' + 1
	ereturn clear
	ereturn scalar N = `nodes'
	ereturn scalar Q = `__q'
	ereturn scalar ICL = `TAB'[`__row', 2]
	ereturn scalar J = `TAB'[`__row', 3]
	ereturn scalar PL = `TAB'[`__row', 4]
	ereturn scalar H = `TAB'[`__row', 5]
	ereturn scalar npar = `TAB'[`__row', 6]
	ereturn scalar directed = `__dir'
	ereturn scalar n_missing = `__nmiss'
	ereturn scalar seed = `__seed'
	ereturn scalar restarts = `restarts'
	ereturn local cmd "nwsbm"
	ereturn local model "`model'"
	ereturn local network "`netname'"
	ereturn local membership = cond(`__dir', "SBM", "SBM_sym")
	matrix colnames `TAB' = Q ICL J PL H npar
	forvalues k = 1/`__q' {
		local __bn `__bn' b`k'
	}
	matrix rownames `PI' = `__bn'
	matrix colnames `PI' = `__bn'
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
	mata: SbmStore("`stub'", 1)
end

capture program drop _nwsbm_display
program _nwsbm_display
	local q = e(Q)
	di
	di as text "Stochastic block model (" as result "`e(model)'" as text ", " as result cond(e(directed), "directed", "undirected") as text ")" ///
		_col(50) "Actors" _col(60) "= " as result %9.0g e(N)
	di as text "Network: " as result "`e(network)'" as text _col(50) "Blocks" _col(60) "= " as result %9.0g `q'
	di as text _col(50) "ICL" _col(60) "= " as result %9.2f e(ICL)
	di as text _col(50) "Lower bound J" _col(60) "= " as result %9.2f e(J)
	tempname T
	matrix `T' = e(ICLtable)
	if rowsof(`T') > 1 {
		di
		di as text "{hline 13}{c TT}{hline 48}"
		di as text %12s "Blocks" " {c |}" %12s "ICL" %12s "J" %12s "PL" %12s "parameters"
		di as text "{hline 13}{c +}{hline 48}"
		forvalues r = 1/`=rowsof(`T')' {
			local mark = cond(`T'[`r', 1] == `q', "*", " ")
			di as text %11.0f `T'[`r', 1] "`mark' {c |}" as result %12.2f `T'[`r', 2] %12.2f `T'[`r', 3] %12.2f `T'[`r', 4] %12.0f `T'[`r', 6]
		}
		di as text "{hline 13}{c BT}{hline 48}"
		di as text "* selected by the largest ICL"
	}
	di
	di as text "Block sizes:"
	matrix list e(blocksize), noheader format(%9.0f)
	di as text "Block proportions:"
	matrix list e(alpha), noheader format(%9.4f)
	di as text "Block connectivity (tie probability from the row block to the column block):"
	matrix list e(pi), noheader format(%9.4f)
end
