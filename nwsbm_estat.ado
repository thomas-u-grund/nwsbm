*! version 1.0.0  20nov2026 nwsbm_estat
* Copyright (C) 2026 Thomas U. Grund. Part of nwsbm; GNU GPL-3 (see LICENSE).
capture program drop nwsbm_estat
program nwsbm_estat, rclass
	version 14
	if "`e(cmd)'" != "nwsbm" error 301
	gettoken sub rest : 0, parse(" ,")
	local len = length("`sub'")
	if "`sub'" == "icl" | ("`sub'" != "" & `len' >= 3 & "`sub'" == substr("icl", 1, `len')) {
		_nwsbm_estat_icl `rest'
		return add
	}
	else if "`sub'" == "blocks" | ("`sub'" != "" & `len' >= 2 & "`sub'" == substr("blocks", 1, `len')) {
		_nwsbm_estat_blocks `rest'
		return add
	}
	else {
		di as error "estat `sub' is not available after nwsbm; use {bf:estat icl} or {bf:estat blocks}"
		error 198
	}
end

capture program drop _nwsbm_estat_icl
program _nwsbm_estat_icl, rclass
	syntax [, noHEader]
	tempname T
	matrix `T' = e(ICLtable)
	local two = ("`e(membership)'" == "LBM")
	local c0 = cond(`two', 3, 2)
	local ncol = colsof(`T')
	di
	di as text "Model selection by the integrated classification likelihood (ICL), `e(model)' model"
	di as text "{hline 13}{c TT}{hline 14}{hline 12}{hline 12}{hline 12}{hline 12}"
	di as text %12s cond(`two', "Blocks 1 x 2", "Blocks") " {c |}" %13s "ICL" %12s "J" %12s "PL" %12s "H" %12s "parameters"
	di as text "{hline 13}{c +}{hline 14}{hline 12}{hline 12}{hline 12}{hline 12}"
	* the starred row is the fit's own selection (e(Q), e(Q2)); the fit takes the first row of the table whose ICL is within 1e-9 of the largest,
	* i.e. on ties the smallest number of blocks (one-mode) or the smallest Q1 and then the smallest Q2 (two-mode)
	local bi = .
	forvalues r = 1/`=rowsof(`T')' {
		if `T'[`r', 1] == e(Q) & (`two' == 0 | `T'[`r', 2] == e(Q2)) local bi = `r'
	}
	forvalues r = 1/`=rowsof(`T')' {
		local mark = cond(`r' == `bi', "*", " ")
		local lab = cond(`two', string(`T'[`r', 1]) + " x " + string(`T'[`r', 2]), string(`T'[`r', 1]))
		di as text %11s "`lab'" "`mark' {c |}" as result %13.2f `T'[`r', `c0'] %12.2f `T'[`r', `c0' + 1] %12.2f `T'[`r', `c0' + 2] %12.2f `T'[`r', `c0' + 3] %12.0f `T'[`r', `c0' + 4]
	}
	di as text "{hline 13}{c BT}{hline 14}{hline 12}{hline 12}{hline 12}{hline 12}"
	di as text "* selected by the largest ICL (ties within 1e-9: fewest blocks)"
	return matrix table = `T'
end

capture program drop _nwsbm_estat_blocks
program _nwsbm_estat_blocks, rclass
	syntax [, noHEader]
	tempname PI SZ AL S
	local two = ("`e(membership)'" == "LBM")
	matrix `PI' = e(pi)
	matrix `SZ' = e(blocksize)
	matrix `AL' = e(alpha)
	local q = colsof(`SZ')
	local what = cond("`e(model)'" == "bernoulli", "within", cond("`e(model)'" == "poisson", "within", "within"))
	di
	if `two' di as text "Blocks of mode 1 (rows of the connectivity matrix)"
	else di as text "Blocks"
	di as text "{hline 8}{c TT}{hline 36}"
	di as text %7s "Block" " {c |}" %10s "Size" %12s "Proportion" %14s cond(`two', "Mean tie", "Within")
	di as text "{hline 8}{c +}{hline 36}"
	matrix `S' = J(`q', 3, .)
	forvalues k = 1/`q' {
		local w = .
		if `two' {
			local s = 0
			forvalues l = 1/`=colsof(`PI')' {
				local s = `s' + `PI'[`k', `l'] * e(alpha2)[1, `l']
			}
			local w = `s'
		}
		else local w = `PI'[`k', `k']
		matrix `S'[`k', 1] = `SZ'[1, `k']
		matrix `S'[`k', 2] = `AL'[1, `k']
		matrix `S'[`k', 3] = `w'
		di as text %7.0f `k' " {c |}" as result %10.0f `SZ'[1, `k'] %12.4f `AL'[1, `k'] %14.4f `w'
	}
	di as text "{hline 8}{c BT}{hline 36}"
	matrix colnames `S' = size proportion connectivity
	return matrix blocks = `S'
	if `two' {
		tempname SZ2 AL2 S2
		matrix `SZ2' = e(blocksize2)
		matrix `AL2' = e(alpha2)
		local q2 = colsof(`SZ2')
		matrix `S2' = J(`q2', 3, .)
		di
		di as text "Blocks of mode 2 (columns of the connectivity matrix)"
		di as text "{hline 8}{c TT}{hline 36}"
		di as text %7s "Block" " {c |}" %10s "Size" %12s "Proportion" %14s "Mean tie"
		di as text "{hline 8}{c +}{hline 36}"
		forvalues l = 1/`q2' {
			local s = 0
			forvalues k = 1/`q' {
				local s = `s' + `PI'[`k', `l'] * `AL'[1, `k']
			}
			matrix `S2'[`l', 1] = `SZ2'[1, `l']
			matrix `S2'[`l', 2] = `AL2'[1, `l']
			matrix `S2'[`l', 3] = `s'
			di as text %7.0f `l' " {c |}" as result %10.0f `SZ2'[1, `l'] %12.4f `AL2'[1, `l'] %14.4f `s'
		}
		di as text "{hline 8}{c BT}{hline 36}"
		matrix colnames `S2' = size proportion connectivity
		return matrix blocks2 = `S2'
	}
	else {
		di as text "Within: connectivity within the block (the diagonal of the connectivity matrix)."
	}
end
