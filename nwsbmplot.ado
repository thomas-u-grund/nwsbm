*! version 1.0.0  20nov2026 nwsbmplot
* Copyright (C) 2026 Thomas U. Grund. Part of nwsbm; GNU GPL-3 (see LICENSE).
capture program drop nwsbmplot
program nwsbmplot
	version 14
	syntax [, TYpe(string) SHADE NOLInes MISSing NOLABels SORT COLors(string asis) EXPort(string asis) Title(string asis) NAME(string) SAVing(string) SIZEs *]
	if "`e(cmd)'" != "nwsbm" error 301
	if "`type'" == "" local type "matrix"
	local type = lower("`type'")
	if "`type'" == "data" local type "matrix"
	if "`type'" == "meso" local type "blocks"
	if !inlist("`type'", "matrix", "blocks", "icl", "expected", "membership") {
		di as error "type() must be matrix (or data), blocks (or meso), icl, expected or membership"
		error 198
	}
	if "`type'" != "icl" {
		capture mata: mata describe __nwsbm_tau
		if _rc {
			di as error "the memberships of the last nwsbm fit are no longer in memory; refit with nwsbm"
			exit 498
		}
	}
	local two = ("`e(membership)'" == "LBM")
	local nmopt = cond("`name'" != "", "name(`name')", "")
	local svopt = cond(`"`saving'"' != "", `"saving("`saving'", replace)"', "")
	* palette of the blocks
	local pal `"`colors'"'
	if `"`pal'"' == "" | `"`pal'"' == "default" local pal "navy maroon forest_green dkorange teal purple khaki cranberry sand ebblue"
	else if `"`pal'"' == "gray" | `"`pal'"' == "grey" local pal "gs2 gs6 gs9 gs12 gs14"
	tokenize `"`pal'"'
	local i = 1
	while `"``i''"' != "" {
		local pal`i' `"``i''"'
		local ++i
	}
	local npal = `i' - 1
	tempname PI
	capture matrix `PI' = e(pi)
	if "`type'" == "icl" {
		tempname T
		matrix `T' = e(ICLtable)
		local c0 = cond(`two', 3, 2)
		local hasbm = (colsof(`T') > `c0' + 4)
		if rowsof(`T') < 2 {
			di as error "only one number of blocks was fitted; there is no ICL path to plot (fit with maxblocks())"
			exit 198
		}
		preserve
		quietly {
			clear
			svmat double `T', names(__t)
			local qs = e(Q)
			gen double __icl = __t`c0'
			if `hasbm' gen double __iclbm = __t`=colsof(`T')'
			if `two' {
				rename __t1 __q1
				rename __t2 __q2
				gen byte __sel = (__q1 == e(Q1) & __q2 == e(Q2))
			}
			else {
				rename __t1 __q1
				gen byte __sel = (__q1 == e(Q))
			}
			summarize __q1, meanonly
			local qlo = r(min)
			local qhi = r(max)
		}
		local xlab "xlabel(`qlo'(1)`qhi')"
		if `qhi' - `qlo' > 20 local xlab "xlabel(#10)"
		local crit = cond("`e(icltype)'" == "exact", "exact ICL", "ICL")
		local plots ""
		local ord ""
		if `two' {
			quietly levelsof __q2, local(q2s)
			local j = 0
			foreach v of local q2s {
				local ++j
				local c = mod(`j' - 1, `npal') + 1
				local plots `plots' (connected __icl __q1 if __q2 == `v', sort msymbol(O) lcolor(`"`pal`c''"') mcolor(`"`pal`c''"'))
				local ord `ord' `j' "mode 2: `v' blocks"
			}
			local plots `plots' (scatter __icl __q1 if __sel, msymbol(D) msize(large) mfcolor(none) mlcolor(black) mlwidth(medthick))
			local ord `ord' `=`j' + 1' "selected"
			local xt "blocks of mode 1"
			local ttl "`crit' by number of blocks (`e(Q1)' x `e(Q2)' selected)"
			local legopt "legend(order(`ord') size(small) cols(1) position(3) region(margin(l+1 r+6 t+1 b+1)))"
		}
		else {
			local plots (connected __icl __q1, sort msymbol(O) lcolor(`"`pal1'"') mcolor(`"`pal1'"'))
			local plots `plots' (scatter __icl __q1 if __sel, msymbol(D) msize(large) mfcolor(none) mlcolor(black) mlwidth(medthick))
			local ord `"1 "`crit'" 2 "selected (`e(Q)' blocks)""'
			if `hasbm' {
				local plots `plots' (connected __iclbm __q1, sort yaxis(2) msymbol(Oh) lpattern(dash) lcolor(gs8) mcolor(gs8))
				local ord `"`ord' 3 "ICL of blockmodels (right axis)""'
			}
			local xt "number of blocks"
			local ttl "`crit' by number of blocks (`e(Q)' selected)"
			local legopt `"legend(order(`ord') size(small) rows(2) position(6) region(margin(l+1 r+6 t+1 b+1)))"'
		}
		if `"`title'"' == "" local title `""`ttl'""'
		local y2 ""
		if !`two' & `hasbm' local y2 `"ytitle("ICL of blockmodels", axis(2)) ylabel(, angle(0) axis(2))"'
		twoway `plots', `xlab' xtitle("`xt'") ytitle("`crit'") ylabel(, angle(0)) `y2' `legopt' title(`title') `nmopt' `svopt' `options'
		restore
	}
	else if "`type'" == "membership" {
		local nmode = cond(`two', 2, 1)
		forvalues md = 1/`nmode' {
			preserve
			quietly {
				clear
				mata: SbmPlotMember(`md')
			}
			local q = cond(`md' == 1, cond(`two', e(Q1), e(Q)), e(Q2))
			local plots ""
			local ord ""
			forvalues k = 1/`q' {
				local c = mod(`k' - 1, `npal') + 1
				local plots `plots' (rbar __cum`=`k' - 1' __cum`k' __x, barwidth(1) fcolor(`"`pal`c''"') lcolor(`"`pal`c''"') lwidth(none))
				local ord `ord' `k' "`k'"
			}
			local xl ""
			if "`nolines'" == "" & "`__xl'" != "" local xl "xline(`__xl', lcolor(black) lwidth(vthin))"
			local legopt "legend(off)"
			if `q' <= 12 local legopt `"legend(order(`ord') rows(1) size(small) title("block", size(small)) position(6))"'
			local sttl ""
			if `two' local sttl `"Mode `md': `q' blocks"'
			else local sttl `"Posterior membership probabilities, `q' blocks"'
			if `"`title'"' != "" & !`two' local sttl `"`title'"'
			local gopt "`nmopt' `svopt'"
			if `two' local gopt `"nodraw name(__nwsbm_mem`md', replace)"'
			twoway `plots', `xl' xscale(range(0.5 `=`__n' + 0.5')) yscale(range(0 1)) ylabel(0(.25)1, angle(0)) xlabel(none, noticks) ///
				plotregion(margin(zero)) xtitle("actors, ordered by block and certainty") ytitle("probability") `legopt' title(`"`sttl'"') `gopt' `options'
			restore
		}
		if `two' {
			if `"`title'"' == "" local title `""Posterior membership probabilities""'
			graph combine __nwsbm_mem1 __nwsbm_mem2, rows(1) title(`title') `nmopt' `svopt'
			capture graph drop __nwsbm_mem1 __nwsbm_mem2
		}
	}
	else if inlist("`type'", "matrix", "expected") {
		if "`type'" == "matrix" {
			_nwsyntax `e(network)', max(1)
			if `two' {
				mata: sbm_setidx(`netobj'->get_modes())
				mata: __nwsbm_A = (*(`netobj'->get_matrix()))[__nwsbm_i1, __nwsbm_i2]
			}
			else {
				mata: __nwsbm_A = *(`netobj'->get_matrix())
				mata: __nwsbm_i1 = J(0, 1, .)
				mata: __nwsbm_i2 = J(0, 1, .)
			}
		}
		preserve
		quietly {
			clear
			if "`type'" == "matrix" mata: SbmPlotMatrix(__nwsbm_A, `two', __nwsbm_i1, __nwsbm_i2, st_matrix("`PI'"), ("`shade'" != ""), ("`missing'" != ""), ("`sort'" != ""))
			else mata: SbmPlotExpected(`two', st_matrix("`PI'"), ("`sort'" != ""))
		}
		capture mata: mata drop __nwsbm_A
		local nr = `__nr'
		local nc = `__nc'
		local ms = cond(max(`nr', `nc') <= 60, "small", cond(max(`nr', `nc') <= 150, "vsmall", cond(max(`nr', `nc') <= 400, "tiny", "vtiny")))
		local plots ""
		local leg ""
		if "`shade'" != "" {
			local shc "gs15 gs13 gs11 gs9"
			quietly count if __rt == 0
			local nt0 = r(N)
			* one bar per block rectangle (the widths differ)
			forvalues i = `=`nt0' + 1'/`=_N' {
				local col : word `=__sh[`i']' of `shc'
				local plots `plots' (rbar __yl __yh __xc if _n == `i', barwidth(`=__w[`i']') fcolor(`col') lcolor(`col') lwidth(none))
			}
		}
		if "`type'" == "expected" {
			local plots ""
			forvalues lv = 1/8 {
				local plots `plots' (scatter __y __x if __lev == `lv', msymbol(S) msize(`ms') mcolor(gs`=16 - 2 * `lv''))
			}
			local isval = 0
			local noplot = 1
		}
		else {
			local noplot = 0
		}
		if !`noplot' {
		quietly count if __rt == 0 & __v < .
		local isval = 0
		if r(N) > 0 {
			quietly summarize __v if __rt == 0 & __v < ., meanonly
			if r(min) != r(max) local isval = 1
		}
		if `isval' {
			quietly _pctile __v if __rt == 0 & __v < ., percentiles(33 66)
			local p1 = r(r1)
			local p2 = r(r2)
			local plots `plots' (scatter __y __x if __rt == 0 & __v < . & __v <= `p1', msymbol(S) msize(`ms') mcolor(gs10))
			local plots `plots' (scatter __y __x if __rt == 0 & __v < . & __v > `p1' & __v <= `p2', msymbol(S) msize(`ms') mcolor(gs6))
			local plots `plots' (scatter __y __x if __rt == 0 & __v < . & __v > `p2', msymbol(S) msize(`ms') mcolor(black))
		}
		else {
			local plots `plots' (scatter __y __x if __rt == 0 & __v < ., msymbol(S) msize(`ms') mcolor(black))
		}
		if "`missing'" != "" local plots `plots' (scatter __y __x if __rt == 0 & __v >= ., msymbol(S) msize(`ms') mcolor(cranberry))
		}
		local xl ""
		local yl ""
		if "`nolines'" == "" {
			if "`__xl'" != "" local xl "xline(`__xl', lcolor(cranberry) lwidth(thin))"
			if "`__yl'" != "" local yl "yline(`__yl', lcolor(cranberry) lwidth(thin))"
		}
		local xlab "xlabel(none, noticks)"
		local ylab "ylabel(none, noticks)"
		if "`nolabels'" == "" {
			local xlab `"xlabel(`__xlab', noticks labsize(small))"'
			local ylab `"ylabel(`__ylab', noticks angle(0) labsize(small))"'
		}
		local what = cond("`type'" == "expected", "Expected ties", "Blocked adjacency matrix")
		if `"`title'"' == "" {
			local title `""`what', `e(Q)' blocks""'
			if `two' local title `""`what', `e(Q1)' x `e(Q2)' blocks""'
		}
		local subopt ""
		if "`type'" == "expected" local subopt `"subtitle("posterior expected value of the tie, `__vmin' (light) to `__vmax' (dark)", size(small))"'
		twoway `plots', `xl' `yl' `xlab' `ylab' xscale(range(0.5 `=`nc'+0.5') lstyle(none)) yscale(range(0.5 `=`nr'+0.5') lstyle(none)) ///
			legend(off) aspectratio(`=`nr'/`nc'') plotregion(margin(zero) lcolor(gs8)) xtitle("") ytitle("") title(`title') `subopt' `nmopt' `svopt' `options'
		restore
	}
	else {
		preserve
		quietly {
			clear
			tempname SZ1 SZ2
			matrix `SZ1' = e(blocksize)
			if `two' matrix `SZ2' = e(blocksize2)
			local nrw = rowsof(`PI')
			local ncl = colsof(`PI')
			if !`two' local ncl = `nrw'
			* nodes: one-mode blocks on a circle; two-mode blocks of mode 1 on the left and of mode 2 on the right
			local nn = `nrw' + cond(`two', `ncl', 0)
			set obs `nn'
			gen double __x = .
			gen double __y = .
			gen double __sz = .
			gen double __ex1 = .
			gen double __ey1 = .
			gen double __ex2 = .
			gen double __ey2 = .
			gen byte __el = .
			gen byte __rt = 0
			local pi = c(pi)
			forvalues k = 1/`nrw' {
				if `two' {
					replace __x = -1 if _n == `k'
					replace __y = `=(`nrw' - 1) / 2 - (`k' - 1)' if _n == `k'
				}
				else {
					replace __x = cos(2 * `pi' * (`k' - 1) / `nrw' + `pi' / 2) if _n == `k'
					replace __y = sin(2 * `pi' * (`k' - 1) / `nrw' + `pi' / 2) if _n == `k'
				}
				replace __sz = `SZ1'[1, `k'] if _n == `k'
			}
			if `two' {
				forvalues l = 1/`ncl' {
					replace __x = 1 if _n == `nrw' + `l'
					replace __y = `=(`ncl' - 1) / 2 - (`l' - 1)' if _n == `nrw' + `l'
					replace __sz = `SZ2'[1, `l'] if _n == `nrw' + `l'
				}
			}
			gen str12 __lab = cond(_n <= `nrw', string(_n), string(_n - `nrw'))
			local mxp = 0
			forvalues a = 1/`nrw' {
				forvalues b = 1/`ncl' {
					local v = abs(`PI'[`a', `b'])
					if `v' > `mxp' local mxp = `v'
				}
			}
			if `mxp' <= 0 local mxp = 1
			local directed = (!`two' & e(directed) == 1)
			local ne = 0
			forvalues a = 1/`nrw' {
				forvalues b = 1/`ncl' {
					if !`two' & (`a' == `b' | (`b' < `a' & !`directed')) continue
					local v = abs(`PI'[`a', `b'])
					if !`two' & !`directed' local v = max(abs(`PI'[`a', `b']), abs(`PI'[`b', `a']))
					if `v' < 0.02 * `mxp' continue
					local lev = max(1, ceil(5 * `v' / `mxp'))
					local ne = `ne' + 1
					local ib = cond(`two', `nrw' + `b', `b')
					local x1 = __x[`a']
					local y1 = __y[`a']
					local x2 = __x[`ib']
					local y2 = __y[`ib']
					if `directed' {
						* arrows leave and enter at the rim of the nodes; the two directions of a pair run side by side
						local dd = sqrt((`x2' - `x1')^2 + (`y2' - `y1')^2)
						local ux = (`x2' - `x1') / `dd'
						local uy = (`y2' - `y1') / `dd'
						local x1 = `x1' + 0.09 * `ux' - 0.035 * `uy'
						local y1 = `y1' + 0.09 * `uy' + 0.035 * `ux'
						local x2 = `x2' - 0.13 * `ux' - 0.035 * `uy'
						local y2 = `y2' - 0.13 * `uy' + 0.035 * `ux'
					}
					set obs `=_N + 1'
					replace __rt = 1 if _n == _N
					replace __el = `lev' in L
					replace __ex1 = `x1' in L
					replace __ey1 = `y1' in L
					replace __ex2 = `x2' in L
					replace __ey2 = `y2' in L
				}
			}
			local eplots ""
			local lws "vthin thin medthin medium thick"
			forvalues w = 1/5 {
				count if __rt == 1 & __el == `w'
				if r(N) > 0 {
					local lw : word `w' of `lws'
					if `directed' local eplots `eplots' (pcarrow __ey1 __ex1 __ey2 __ex2 if __rt == 1 & __el == `w', lwidth(`lw') lcolor(gs8) mcolor(gs8) msize(small))
					else local eplots `eplots' (pcspike __ey1 __ex1 __ey2 __ex2 if __rt == 1 & __el == `w', lwidth(`lw') lcolor(gs8))
				}
			}
			summarize __sz if __rt == 0, meanonly
			local mxs = r(max)
			local splots ""
			local msz "2 3 4 5 6 7"
			local usec = (`"`colors'"' != "")
			forvalues k = 1/`nn' {
				local cl = ceil(6 * __sz[`k'] / `mxs')
				local c = mod(cond(`two' & `k' > `nrw', `k' - `nrw', `k') - 1, `npal') + 1
				local lbc = cond(`usec', "white", "black")
				local fc = cond(`usec', `"`pal`c''"', "gs12")
				local splots `splots' (scatter __y __x if _n == `k', msymbol(O) msize(*`: word `cl' of `msz'') mcolor(`"`fc'"') mlcolor(gs4) mlabel(__lab) mlabposition(0) mlabcolor(`lbc') mlabsize(small))
			}
		}
		if `"`title'"' == "" local title `""Block graph""'
		local sub `""node size: block size; line width: connectivity (relative to the largest entry of the connectivity matrix)""'
		if `two' local sub `""left: blocks of mode 1, right: blocks of mode 2; node size: block size; line width: connectivity""'
		twoway `eplots' `splots', legend(off) xscale(off range(-1.3 1.3)) yscale(off range(-1.3 1.3)) xlabel(none) ylabel(none) ///
			aspectratio(1) title(`title') subtitle(`sub', size(small)) `nmopt' `svopt' `options'
		restore
	}
	if `"`export'"' != "" {
		if "`name'" != "" {
			gettoken gname : name, parse(",")
			graph export `export', name(`gname') replace
		}
		else graph export `export', replace
	}
end
