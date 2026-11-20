* nwsbmplot: every plot type and option on a one-mode undirected, a directed and a two-mode fit; the graphs exist; error paths
do cscripts/nwsbm_testlib.do
set linesize 200
set graphics off
global nwfail 0
local tmp "$SBM_SCRATCH"

capture program drop chk
program chk
	* chk name, cmd(...) : run nwsbmplot with the options and check that the named graph exists afterwards
	syntax name(name=gn), LABel(string) [OPTS(string asis)]
	capture graph drop `gn'
	capture noisily nwsbmplot, `opts' name(`gn', replace)
	if _rc {
		di as error "FAIL `label': r(" _rc ")"
		global nwfail = $nwfail + 1
		exit
	}
	capture graph describe `gn'
	if _rc {
		di as error "FAIL `label': graph not found"
		global nwfail = $nwfail + 1
		exit
	}
	di as result "PASS `label'"
end

* three fits
sbmt_mknet gu, ref(pp_u_n150_Q3) n(150) dir(0)
sbmt_mknet gd, ref(pp_d_n100_Q3) n(100) dir(1)
sbmt_mknet2 g2, ref(lbm_po_40x60) n1(40) n2(60) valued

foreach fit in u d t {
	if "`fit'" == "u" {
		quietly nwsbm gu, maxblocks(5) seed(3) nodisplay
		local lab "undirected"
	}
	if "`fit'" == "d" {
		quietly nwsbm gd, maxblocks(5) seed(3) nodisplay
		local lab "directed"
	}
	if "`fit'" == "t" {
		quietly nwsbm g2, maxblocks(3 3) model(poisson) seed(3) nodisplay
		local lab "two-mode"
	}
	foreach o in "" "type(matrix)" "type(data)" "shade" "nolines" "nolabels" "sort" "missing" "shade sort title(Test)" ///
		"type(blocks)" "type(meso)" "type(blocks) colors(red blue green orange purple)" "type(blocks) colors(gray)" "type(blocks) colors(default) title(X)" ///
		"type(expected)" "type(expected) sort" "type(expected) nolines nolabels" "type(expected) scheme(s1mono)" ///
		"type(membership)" "type(membership) colors(gray)" "type(membership) colors(red blue green) nolines" "type(membership) title(M)" ///
		"type(icl)" "type(icl) colors(red)" "type(icl) scheme(s1mono)" {
		chk __g, opts(`o') label("`lab': `o'")
	}
	capture graph drop __g
	* saving and export
	capture erase "`tmp'/plot_`fit'.gph"
	capture erase "`tmp'/plot_`fit'.png"
	capture erase "`tmp'/plot_`fit'_m.png"
	nwsbmplot, type(matrix) saving("`tmp'/plot_`fit'.gph") export("`tmp'/plot_`fit'.png")
	capture confirm file "`tmp'/plot_`fit'.gph"
	if _rc {
		di as error "FAIL `lab': saving()"
		global nwfail = $nwfail + 1
	}
	capture confirm file "`tmp'/plot_`fit'.png"
	if _rc {
		di as error "FAIL `lab': export() png"
		global nwfail = $nwfail + 1
	}
	else di as result "PASS `lab': saving and export"
	nwsbmplot, type(membership) name(__gm, replace) export("`tmp'/plot_`fit'_m.png")
	capture confirm file "`tmp'/plot_`fit'_m.png"
	if _rc {
		di as error "FAIL `lab': export() png with name()"
		global nwfail = $nwfail + 1
	}
	else di as result "PASS `lab': export png with name()"
	capture graph drop _all
	* ordering within blocks by certainty (the first arguments of the Mata helper)
	mata: st_numscalar("__bad", sbmt_certcheck())
	if scalar(__bad) {
		di as error "FAIL `lab': sort order within blocks is not by decreasing certainty"
		global nwfail = $nwfail + 1
	}
	else di as result "PASS `lab': certainty order"
}

* ICL path with the exact criterion (adds the blockmodels ICL on the right axis); the selected row is marked
quietly nwsbm gu, maxblocks(5) seed(3) icl(exact) nodisplay
chk __g, opts(type(icl)) label("exact ICL path")
quietly nwsbm g2, maxblocks(3 3) model(poisson) seed(3) icl(exact) degreecorrected nodisplay
chk __g, opts(type(icl)) label("exact ICL path, two-mode degree-corrected")

* error paths
quietly nwsbm gu, blocks(3) seed(3) nodisplay
capture nwsbmplot, type(icl)
if _rc != 198 {
	di as error "FAIL type(icl) with a single fit should give r(198), got " _rc
	global nwfail = $nwfail + 1
}
else di as result "PASS icl: single fit refused"
capture nwsbmplot, type(foo)
if _rc != 198 {
	di as error "FAIL type(foo) should give r(198), got " _rc
	global nwfail = $nwfail + 1
}
else di as result "PASS bad type"
capture nwsbmplot, type(blocks) export("`tmp'/plot_bad.xyz")
if !_rc {
	di as error "FAIL export with an unknown extension should fail"
	global nwfail = $nwfail + 1
}
else di as result "PASS bad export extension"
mata: mata drop __nwsbm_tau __nwsbm_tau2
capture nwsbmplot, type(membership)
if _rc != 498 {
	di as error "FAIL memberships gone should give r(498), got " _rc
	global nwfail = $nwfail + 1
}
else di as result "PASS memberships gone"
capture nwsbmplot, type(icl)
* the ICL path needs only the table: not an error to have lost the memberships (single fit gives 198 here)
ereturn clear
capture nwsbmplot
if _rc != 301 {
	di as error "FAIL no estimates should give r(301), got " _rc
	global nwfail = $nwfail + 1
}
else di as result "PASS no estimates"

* regressions of 1.0.0: shade draws rectangles; two-mode estat icl header; two-mode rows are mode 1
quietly nwsbm gu, maxblocks(4) seed(3) nodisplay
nwsbmplot, type(matrix) export("`tmp'/sh0.png")
nwsbmplot, type(matrix) shade export("`tmp'/sh1.png")
quietly checksum "`tmp'/sh0.png"
local h0 = r(checksum)
quietly checksum "`tmp'/sh1.png"
if "`h0'" == "`r(checksum)'" {
	di as error "FAIL shade has no effect"
	global nwfail = $nwfail + 1
}
else di as result "PASS shade"
* the shaded rectangles: the dark grey fraction of the picture must be large (strips would be tiny); compare file sizes of shade vs unshaded is weak, so test the plot data instead
preserve
quietly {
	clear
	mata: SbmPlotMatrix(J(150, 150, 0), 0, J(0, 1, .), J(0, 1, .), st_matrix("e(pi)"), 1, 0, 0)
	summarize __w if __rt == 1, meanonly
}
if r(max) <= 1 {
	di as error "FAIL shade rectangles have unit width"
	global nwfail = $nwfail + 1
}
else di as result "PASS shade rectangle widths"
restore
quietly nwsbm g2, maxblocks(3 3) model(poisson) seed(3) nodisplay
tempname lg
quietly log using "`tmp'/icl2.log", text replace name(icl2)
estat icl
quietly log close icl2
quietly filefilter "`tmp'/icl2.log" "`tmp'/icl2b.log", from("Blocks 1 x 2") to("STALE") replace
if r(occurrences) > 0 {
	di as error "FAIL stale header Blocks 1 x 2"
	global nwfail = $nwfail + 1
}
else di as result "PASS estat icl header"
nwsbmplot, type(matrix) name(__gt, replace)
* mode 1 on the rows: the plot data has as many rows as actors of mode 1
preserve
quietly {
	clear
	mata: SbmPlotMatrix(J(`=e(N1)', `=e(N2)', 1), 1, (1::`=e(N1)'), (1::`=e(N2)'), st_matrix("e(pi)"), 0, 0, 0)
}
if `__nr' != e(N1) | `__nc' != e(N2) {
	di as error "FAIL two-mode rows are not mode 1"
	global nwfail = $nwfail + 1
}
else di as result "PASS two-mode rows are mode 1"
restore

if $nwfail di as error "FAILED $nwfail"
else di as result "ALL PASSED"
