* Ties in the ICL table: the fit and estat icl must star the same cell, by a documented rule.
* Rule: the selected cell is the first row of the table (smallest Q, or smallest Q1 and then smallest Q2) whose ICL is within 1e-9 of the largest ICL.
* Data: Davis, Gardner and Gardner (1941), Southern women, 18 x 14 (dev/ref/davis_matrix.csv). The degree-corrected Poisson model with icl(exact) gives
* the exact same ICL for 1 x 1, 1 x 2, ..., 4 x 1 (fits with empty or merged blocks have the partition of the one-block fit). Before the fix estat icl
* starred 4 x 1 (its running maximum went through a rounded macro) while the fit's table starred 1 x 1.
do cscripts/nwsbm_testlib.do
set linesize 200
local fail 0
sbmt_mkcov dv, file("$SBM_REPO/dev/ref/davis_matrix.csv") dir(2)
quietly nwsbm dv, maxblocks(4 4) seed(20261120) model(poisson) degreecorrected icl(exact) nodisplay
matrix T = e(ICLtable)
* expected row by the rule
local mx = -1e300
forvalues r = 1/`=rowsof(T)' {
	if T[`r', 3] > `mx' local mx = T[`r', 3]
}
local want = .
forvalues r = `=rowsof(T)'(-1)1 {
	if T[`r', 3] >= `mx' - 1e-9 local want = `r'
}
local nt 0
forvalues r = 1/`=rowsof(T)' {
	if abs(T[`r', 3] - `mx') <= 1e-9 local ++nt
}
di as text "two-mode: " `nt' " cells tie, expected selection row `want' = " T[`want', 1] " x " T[`want', 2] "; fit selected " e(Q1) " x " e(Q2)
if `nt' < 2 {
	di as error "the network no longer produces tied cells"
	local ++fail
}
if e(Q1) != T[`want', 1] | e(Q2) != T[`want', 2] {
	di as error "fit does not select the first tied cell"
	local ++fail
}
else di as result "PASS fit selects the first of the tied cells (two-mode)"
tempfile lg
quietly log using "`lg'", text name(icltie) replace
estat icl
quietly log close icltie
tempname fh
file open `fh' using "`lg'", read text
local star ""
file read `fh' line
while r(eof) == 0 {
	if strpos(`"`line'"', "*") & strpos(`"`line'"', "x") & strpos(`"`line'"', "{c |}") + strpos(`"`line'"', "|") > 0 & strpos(`"`line'"', "largest") == 0 {
		local star `"`line'"'
	}
	file read `fh' line
}
file close `fh'
di as text `"starred line of estat icl: `star'"'
local lab = string(T[`want', 1]) + " x " + string(T[`want', 2]) + "*"
if strpos(`"`star'"', "`lab'") == 0 {
	di as error "estat icl does not star `lab'"
	local ++fail
}
else di as result "PASS estat icl stars the same cell as the fit (two-mode)"

* one-mode network with tied cells is not needed; the rule is the same code. Check estat icl on a one-mode fit stars e(Q).
sbmt_mknet g, ref(gr_sbm_u_n100_Q3) n(100) dir(0)
quietly nwsbm g, model(bernoulli) maxblocks(6) icl(exact) seed(5) nodisplay
quietly log using "`lg'", text name(icltie) replace
estat icl
quietly log close icltie
file open `fh' using "`lg'", read text
local star ""
file read `fh' line
while r(eof) == 0 {
	if strpos(`"`line'"', "*") & strpos(`"`line'"', "|") & strpos(`"`line'"', "largest") == 0 local star `"`line'"'
	file read `fh' line
}
file close `fh'
if strpos(`"`star'"', "`=e(Q)'*") == 0 {
	di as error "one-mode: estat icl does not star Q = `=e(Q)': `star'"
	local ++fail
}
else di as result "PASS estat icl stars e(Q) (one-mode)"
if `fail' di as error "FAILED `fail'"
else di as result "ALL PASSED"
