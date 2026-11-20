* The printed output of nwsbm, estat icl and estat blocks does not wrap at linesize 80 (a wrapped line starts with "> " in the log).
do cscripts/nwsbm_testlib.do
global nwfail 0
set linesize 80
local lg "$SBM_SCRATCH/linesize_`c(pid)'.log"
capture log close _all
log using "`lg'", replace text name(lsz)
nwwebuse lazega, nwclear
nwsbm lazega_adv, maxblocks(8)
estat icl
estat blocks
nwsbm lazega_adv, maxblocks(5) icl(exact)
estat icl
nwsbm lazega_adv, blocks(3) model(poisson) degreecorrected
nwsbm lazega_adv, blocks(2) model(gaussian) blockvar
log close lsz
set linesize 200
tempname fh
local nwrap 0
local nok 0
file open `fh' using "`lg'", read text
file read `fh' line
while r(eof) == 0 {
	if substr(`"`macval(line)'"', 1, 2) == "> " & strpos(`"`macval(line)'"', "linesize_") == 0 & strpos(`"`macval(line)'"', "scratchpad") == 0 & strpos(`"`macval(line)'"', "nwsbm_tests") == 0 {
		local ++nwrap
		di as error "wrapped line: `macval(line)'"
	}
	if strpos(`"`macval(line)'"', "Lower bound J") > 0 local ++nok
	file read `fh' line
}
file close `fh'
if `nwrap' > 0 {
	di as error "FAIL: `nwrap' wrapped lines"
	global nwfail 1
}
if `nok' < 3 {
	di as error "FAIL: the header lines were not found"
	global nwfail 1
}
if $nwfail == 0 di as result "PASS no wrapped lines at linesize 80"
else exit 9
