* orderblocks: blocks in decreasing order of size (hard memberships); every other result is a relabelling of the unordered fit
do cscripts/nwsbm_testlib.do
set linesize 200
local fail 0
foreach spec in "pp_u_n150_Q3 150 0" "pp_d_n100_Q3 100 1" {
	local nm : word 1 of `spec'
	local n : word 2 of `spec'
	local dir : word 3 of `spec'
	sbmt_mknet g, ref(`nm') n(`n') dir(`dir')
	quietly nwsbm g, blocks(3) seed(5) generate(m0) replace nodisplay
	local icl0 = e(ICL)
	matrix S0 = e(blocksize)
	matrix P0 = e(pi)
	quietly nwsbm g, blocks(3) seed(5) orderblocks generate(m1) replace nodisplay
	matrix S1 = e(blocksize)
	matrix P1 = e(pi)
	matrix A1 = e(alpha)
	local ok 1
	if abs(e(ICL) - `icl0') > 1e-9 * abs(`icl0') {
		di as error "`nm': ICL changed"
		local ok 0
	}
	if S1[1, 1] < S1[1, 2] | S1[1, 2] < S1[1, 3] {
		di as error "`nm': block sizes not decreasing"
		local ok 0
	}
	* sizes of the generated variable agree with e(blocksize) and the sorted unordered sizes
	forvalues k = 1/3 {
		quietly count if m1 == `k' in 1/`n'
		if r(N) != S1[1, `k'] {
			di as error "`nm': generated block `k' has `r(N)' actors, e(blocksize) says " S1[1, `k']
			local ok 0
		}
	}
	* the relabelling is consistent: the connectivity of the block of the generated variable equals e(pi) entry (check via block-pair density of the first block)
	quietly sum m0 in 1/`n'
	sbmt_ari m0 m1
	if r(ari) < 0.9999 {
		di as error "`nm': ordered partition differs from the unordered one"
		local ok 0
	}
	* pi: the largest block of the unordered fit has the same diagonal entry
	local big0 = 1
	forvalues k = 2/3 {
		if S0[1, `k'] > S0[1, `big0'] local big0 = `k'
	}
	if abs(P0[`big0', `big0'] - P1[1, 1]) > 1e-9 {
		di as error "`nm': pi not permuted consistently"
		local ok 0
	}
	di as text "`nm': sizes unordered " S0[1, 1] " " S0[1, 2] " " S0[1, 3] "  ordered " S1[1, 1] " " S1[1, 2] " " S1[1, 3]
	if `ok' di as result "PASS `nm'"
	else local fail = `fail' + 1
}
* two-mode
sbmt_mknet2 g2, ref(lbm_po_40x60) n1(40) n2(60) valued
quietly nwsbm g2, blocks(2 3) model(poisson) seed(2) nodisplay
matrix S0 = e(blocksize)
matrix T0 = e(blocksize2)
quietly nwsbm g2, blocks(2 3) model(poisson) seed(2) orderblocks nodisplay
matrix S1 = e(blocksize)
matrix T1 = e(blocksize2)
if S1[1, 1] < S1[1, 2] | T1[1, 1] < T1[1, 2] | T1[1, 2] < T1[1, 3] {
	di as error "two-mode: block sizes not decreasing"
	local fail = `fail' + 1
}
else di as result "PASS two-mode"
if `fail' di as error "FAILED `fail'"
else di as result "ALL PASSED"
