* Standard errors of the covariate coefficients (e(se_beta), e(V_beta)): inverse of the information of the coefficients with the block parameters profiled
* out, conditional on the estimated memberships. R (sbm / blockmodels) reports no standard errors, so the references are direct computations:
*  (a) one block (Q = 1): the model is an ordinary logit / Poisson / linear regression of the pairs; the standard errors must equal those of
*      logit / poisson / regress (the Gaussian ones rescaled from the unbiased to the maximum likelihood variance), to 1e-5 relative;
*  (b) Q blocks with well separated memberships: the same regression with a dummy for every block pair at the hard memberships; relative
*      difference at most 1 percent (the memberships are fuzzy to a small degree and the likelihood is the expected one).
* Missing ties are dropped from the pair data (the plugin masks them).
do cscripts/nwsbm_testlib.do
set linesize 200
mata:
// pair data: y, x1..xp, block of i and of j (z1, z2: memberships of rows and columns; one-mode: same vector); kind 1 = two-mode
void sbmt_pairs(string scalar Aname, real scalar np, real scalar dir, real scalar kind, real colvector z1, real colvector z2)
{
	real matrix A, D, X
	real scalar i, j, k, r
	A = st_matrix(Aname)
	D = J(0, 3 + np, .)
	for (i = 1; i <= rows(A); i++) for (j = 1; j <= cols(A); j++) {
		if (!kind & i == j) continue
		if (!kind & !dir & j < i) continue
		if (A[i, j] >= .) continue
		X = J(1, np, .)
		for (k = 1; k <= np; k++) X[k] = st_matrix("XM" + strofreal(k))[i, j]
		D = D \ (A[i, j], X, z1[i], z2[j])
	}
	st_addobs(rows(D))
	(void) st_addvar("double", ("y", "x" :+ strofreal(1..np), "z1", "z2"))
	st_store(., ., D)
}
end
local fix "cb_u_n100_Q3 100 0 bernoulli 1 3 cb_d_n80_Q2 80 1 bernoulli 2 2 cp_u_n80_Q3 80 0 poisson 1 3 cp_d_n70_Q2 70 1 poisson 2 2 cg_u_n80_Q2 80 0 gaussian 1 2 cg_d_n70_Q3 70 1 gaussian 1 3 cbm_d_n60_Q2 60 1 bernoulli 2 2 cpm_u_n60_Q2 60 0 poisson 1 2"
local fail 0
local nfix : word count `fix'
forvalues f = 1(6)`nfix' {
	local nm : word `f' of `fix'
	local n : word `=`f'+1' of `fix'
	local dir : word `=`f'+2' of `fix'
	local mod : word `=`f'+3' of `fix'
	local np : word `=`f'+4' of `fix'
	local Q : word `=`f'+5' of `fix'
	sbmt_mknetv g, ref(`nm') n(`n') dir(`dir')
	mata: sbmt_loadv(`n', `n', `dir', "$SBM_REPO/dev/ref/`nm'_edges.csv", "$SBM_REPO/dev/ref/`nm'_missing.csv")
	matrix AA = __sbmt_A
	local covs
	forvalues k = 1/`np' {
		sbmt_mkcov x`k', file("$SBM_REPO/dev/ref/`nm'_cov`k'.csv") dir(`dir')
		local covs `covs' x`k'
		preserve
		import delimited "$SBM_REPO/dev/ref/`nm'_cov`k'.csv", clear varnames(nonames)
		mkmat _all, matrix(XM`k')
		restore
	}
	local ok 1
	* (a) one block
	quietly nwsbm g, blocks(1) model(`mod') covariates(`covs') seed(3) nodisplay
	matrix SE1 = e(se_beta)
	matrix V1 = e(V_beta)
	matrix B1 = e(beta)
	preserve
	clear
	mata: sbmt_pairs("AA", `np', `dir', 0, J(`n', 1, 1), J(`n', 1, 1))
	local xs
	forvalues k = 1/`np' {
		local xs `xs' x`k'
	}
	if "`mod'" == "bernoulli" quietly logit y `xs'
	else if "`mod'" == "poisson" quietly poisson y `xs'
	else {
		quietly regress y `xs'
	}
	local fac = 1
	if "`mod'" == "gaussian" local fac = sqrt((e(N) - e(df_m) - 1) / e(N))
	forvalues k = 1/`np' {
		local rel = abs(SE1[1, `k'] - _se[x`k'] * `fac') / (_se[x`k'] * `fac')
		di as text "   `nm' (a) coef `k': se nwsbm " %10.6f SE1[1, `k'] "  reference " %10.6f _se[x`k'] * `fac' "  rel diff " %9.2e `rel'
		if `rel' > 1e-5 | missing(`rel') {
			di as error "`nm': one-block standard error of covariate `k' differs from the regression by `rel'"
			local ok 0
		}
	}
	restore
	* e(V_beta) symmetric, diagonal = se^2
	if abs(V1[1, 1] - SE1[1, 1]^2) > 1e-12 | (`np' > 1 & abs(V1[1, 2] - V1[2, 1]) > 1e-12) {
		di as error "`nm': V_beta inconsistent with se_beta"
		local ok 0
	}
	* (b) Q blocks
	quietly nwsbm g, blocks(`Q') model(`mod') covariates(`covs') seed(7) generate(mb) replace nodisplay
	matrix SEQ = e(se_beta)
	mata: st_matrix("__ZZ", st_data((1::`n'), "mb"))
	preserve
	clear
	mata: sbmt_pairs("AA", `np', `dir', 0, st_matrix("__ZZ"), st_matrix("__ZZ"))
	if `dir' == 0 {
		gen zlo = min(z1, z2)
		gen zhi = max(z1, z2)
		egen bp = group(zlo zhi)
	}
	else egen bp = group(z1 z2)
	if "`mod'" == "bernoulli" quietly logit y `xs' i.bp
	else if "`mod'" == "poisson" quietly poisson y `xs' i.bp
	else quietly regress y `xs' i.bp
	local fac = 1
	if "`mod'" == "gaussian" local fac = sqrt((e(N) - e(df_m) - 1) / e(N))
	local maxrel 0
	forvalues k = 1/`np' {
		local rel = abs(SEQ[1, `k'] - _se[x`k'] * `fac') / (_se[x`k'] * `fac')
		if `rel' > `maxrel' local maxrel = `rel'
		di as text "   `nm' (b) coef `k': se nwsbm " %10.6f SEQ[1, `k'] "  hard-membership regression " %10.6f _se[x`k'] * `fac' "  rel diff " %9.2e `rel'
	}
	restore
	if `maxrel' > 0.01 | missing(`maxrel') {
		di as error "`nm': standard errors at Q = `Q' differ from the hard-membership regression by `maxrel'"
		local ok 0
	}
	if `ok' di as result "PASS `nm'"
	else local fail = `fail' + 1
}
* two-mode, one block against the regression
local fix2 "cl_b_40x50 40 50 bernoulli cl_p_40x50 40 50 poisson cl_g_40x50 40 50 gaussian"
local nfix : word count `fix2'
forvalues f = 1(4)`nfix' {
	local nm : word `f' of `fix2'
	local n1 : word `=`f'+1' of `fix2'
	local n2 : word `=`f'+2' of `fix2'
	local mod : word `=`f'+3' of `fix2'
	sbmt_mknet2 g, ref(`nm') n1(`n1') n2(`n2') valued
	mata: sbmt_loadv(`n1', `n2', -1, "$SBM_REPO/dev/ref/`nm'_edges.csv", "$SBM_REPO/dev/ref/`nm'_missing.csv")
	matrix AA = __sbmt_A
	sbmt_mkcov x1, file("$SBM_REPO/dev/ref/`nm'_cov1.csv") dir(2)
	preserve
	import delimited "$SBM_REPO/dev/ref/`nm'_cov1.csv", clear varnames(nonames)
	mkmat _all, matrix(XM1)
	restore
	quietly nwsbm g, blocks(1 1) model(`mod') covariates(x1) seed(3) nodisplay
	matrix SE1 = e(se_beta)
	preserve
	clear
	mata: sbmt_pairs("AA", 1, -1, 1, J(`n1', 1, 1), J(`n2', 1, 1))
	if "`mod'" == "bernoulli" quietly logit y x1
	else if "`mod'" == "poisson" quietly poisson y x1
	else quietly regress y x1
	local fac = 1
	if "`mod'" == "gaussian" local fac = sqrt((e(N) - 2) / e(N))
	local rel = abs(SE1[1, 1] - _se[x1] * `fac') / (_se[x1] * `fac')
	di as text "   `nm' two-mode one block: se nwsbm " %10.6f SE1[1, 1] "  reference " %10.6f _se[x1] * `fac' "  rel diff " %9.2e `rel'
	restore
	if `rel' > 1e-5 | missing(`rel') {
		di as error "`nm': two-mode one-block standard error differs by `rel'"
		local fail = `fail' + 1
	}
	else di as result "PASS `nm' (two-mode)"
}
if `fail' di as error "FAILED `fail'"
else di as result "ALL PASSED"
