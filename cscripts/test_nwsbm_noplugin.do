* Without the compiled plugin nwsbm stops with r(601) and a message; there is no Mata fallback. The package files are copied
* to a scratch directory that has no plugins/ folder.
capture program drop nwinstall
local core : env NWCORE
if "`core'" == "" local core "/Users/tgrund/FILES_NEW/SOFTWARE/nwcommands"
local scratch : env NWSBM_SCRATCH
if "`scratch'" == "" local scratch "`c(tmpdir)'/nwsbm_noplugin_test"
capture mkdir "`scratch'"
local d "`scratch'/noplugin"
capture mkdir "`d'"
capture mkdir "`d'/plus"
sysdir set PLUS "`d'/plus"
local repo "`c(pwd)'"
foreach f in nwsbm.ado _nwsbm_libcheck.ado nwsbm_version.ado lib/lnwsbm.mlib {
	local b = subinstr("`f'", "lib/", "", 1)
	copy "`repo'/`f'" "`d'/`b'", replace
}
adopath ++ "`core'"
adopath ++ "`core'/lib"
adopath ++ "`d'"
set obs 5
forvalues j = 1/5 {
	gen double __a`j' = (`j' == mod(_n, 5) + 1)
}
nwset __a*, name(g) undirected
capture nwsbm g, blocks(2)
if _rc == 601 di as result "PASS missing plugin stops with r(601)"
else {
	di as error "FAIL: r(" _rc ") instead of r(601)"
	exit 9
}
