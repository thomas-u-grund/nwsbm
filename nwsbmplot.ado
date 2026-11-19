*! version 0.3.0  19nov2026 nwsbmplot
* Copyright (C) 2026 Thomas U. Grund. Part of nwsbm; GNU GPL-3 (see LICENSE).
capture program drop nwsbmplot
program nwsbmplot
	version 14
	syntax [, TYpe(string) SHADE NOLInes MISSing NOLABels Title(string asis) NAME(string) SAVing(string) SIZEs *]
	if "`e(cmd)'" != "nwsbm" error 301
	if "`type'" == "" local type "matrix"
	if !inlist("`type'", "matrix", "blocks") {
		di as error "type() must be matrix or blocks"
		error 198
	}
	capture mata: mata describe __nwsbm_tau
	if _rc {
		di as error "the memberships of the last nwsbm fit are no longer in memory; refit with nwsbm"
		exit 498
	}
	local two = ("`e(membership)'" == "LBM")
	tempname PI
	matrix `PI' = e(pi)
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
		preserve
		quietly {
			clear
			mata: SbmPlotMatrix(__nwsbm_A, `two', __nwsbm_i1, __nwsbm_i2, st_matrix("`PI'"), ("`shade'" != ""), ("`missing'" != ""))
		}
		capture mata: mata drop __nwsbm_A
		local nr = `__nr'
		local nc = `__nc'
		local ms = cond(max(`nr', `nc') <= 60, "small", cond(max(`nr', `nc') <= 150, "vsmall", cond(max(`nr', `nc') <= 400, "tiny", "vtiny")))
		local plots ""
		local leg ""
		if "`shade'" != "" {
			local shc "gs15 gs13 gs11 gs9"
			forvalues s = 1/4 {
				local col : word `s' of `shc'
				local plots `plots' (rbar __yl __yh __xc if __rt == 1 & __sh == `s', barwidth(1) fcolor(`col') lcolor(`col') lwidth(none))
			}
		}
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
		if `"`title'"' == "" {
			local title `""Blocked adjacency matrix, `e(Q)' blocks""'
			if `two' local title `""Blocked adjacency matrix, `e(Q1)' x `e(Q2)' blocks""'
		}
		local nmopt = cond("`name'" != "", "name(`name')", "")
		local svopt = cond(`"`saving'"' != "", `"saving("`saving'", replace)"', "")
		twoway `plots', `xl' `yl' `xlab' `ylab' xscale(range(0.5 `=`nc'+0.5') lstyle(none)) yscale(range(0.5 `=`nr'+0.5') lstyle(none)) ///
			legend(off) aspectratio(`=`nr'/`nc'') plotregion(margin(zero) lcolor(gs8)) xtitle("") ytitle("") title(`title') `nmopt' `svopt' `options'
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
					set obs `=_N + 1'
					replace __rt = 1 if _n == _N
					replace __el = `lev' in L
					replace __ex1 = __x[`a'] in L
					replace __ey1 = __y[`a'] in L
					replace __ex2 = __x[`ib'] in L
					replace __ey2 = __y[`ib'] in L
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
			forvalues c = 1/6 {
				count if __rt == 0 & ceil(6 * __sz / `mxs') == `c'
				if r(N) > 0 local splots `splots' (scatter __y __x if __rt == 0 & ceil(6 * __sz / `mxs') == `c', msymbol(O) msize(*`: word `c' of `msz'') mcolor(gs12) mlcolor(gs4) mlabel(__lab) mlabposition(0) mlabcolor(black) mlabsize(small))
			}
		}
		if `"`title'"' == "" local title `""Block graph""'
		local sub `""node size: block size; line width: connectivity (relative to the largest entry of the connectivity matrix)""'
		if `two' local sub `""left: blocks of mode 1, right: blocks of mode 2; node size: block size; line width: connectivity""'
		local nmopt = cond("`name'" != "", "name(`name')", "")
		local svopt = cond(`"`saving'"' != "", `"saving("`saving'", replace)"', "")
		twoway `eplots' `splots', legend(off) xscale(off range(-1.3 1.3)) yscale(off range(-1.3 1.3)) xlabel(none) ylabel(none) ///
			aspectratio(1) title(`title') subtitle(`sub', size(small)) `nmopt' `svopt' `options'
		restore
	}
end
