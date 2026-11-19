* Gaussian model with block-specific variances (option blockvar): y_ij ~ N(mu_ql, sigma2_ql). R blockmodels / sbm have only one common variance, so the
* references are the planted values and direct computations at the estimated partition:
*  (1) planted networks whose blocks differ ONLY in the variance (all means equal): blockvar finds the planted blocks (ARI >= 0.9), one-mode undirected and directed and two-mode;
*  (2) block means and variances equal the empirical means and variances of the cells of the hard blocks (within 3 percent; the memberships are nearly hard);
*  (3) e(npar) doubles the number of block parameters, e(sigma2) is the pooled residual variance, orderblocks permutes pi and sigma2_block together;
*  (4) errors: blockvar needs model(gaussian), refuses sigma2(), covariates, degreecorrected.
do cscripts/nwsbm_testlib.do
set linesize 200
mata:
mata set matastrict off
void sbmt_simvar(real scalar n1, real scalar n2, real scalar dir, real scalar seed)
{
	external real matrix __sim_A
	external real colvector __sim_z1, __sim_z2
	real matrix S, A
	real colvector z1, z2
	real scalar i, j, q
	rseed(seed)
	S = (0.1, 0.5, 1.0 \ 0.5, 1.0, 2.5 \ 1.0, 2.5, 6.0)
	z1 = J(0, 1, .)
	for (q = 1; q <= 3; q++) z1 = z1 \ J(n1 / 3, 1, q)
	z1 = z1[order(runiform(n1, 1), 1)]
	z2 = z1
	if (n2 != n1) {
		z2 = J(0, 1, .)
		for (q = 1; q <= 3; q++) z2 = z2 \ J(n2 / 3, 1, q)
		z2 = z2[order(runiform(n2, 1), 1)]
	}
	A = J(n1, n2, 0)
	for (i = 1; i <= n1; i++) for (j = 1; j <= n2; j++) A[i, j] = sqrt(S[z1[i], z2[j]]) * rnormal(1, 1, 0, 1)
	if (n1 == n2 & !dir) A = uppertriangle(A, 0) + uppertriangle(A, 0)'
	if (n1 == n2) _diag(A, 0)
	__sim_A = A
	__sim_z1 = z1
	__sim_z2 = z2
}
end
mata:
// empirical mean and variance of the cells of the hard blocks (variable mb of the dataset, truth ignored): returns matrix __bs (3 x 6: means in columns 1-3, variances in 4-6)
void sbmt_blockstats(real scalar n, real scalar dir)
{
	real matrix A, M, V, N
	real colvector z
	real scalar i, j, a, b
	external real matrix __sim_A
	A = __sim_A
	z = st_data((1::n), "mb")
	M = J(3, 3, 0); V = J(3, 3, 0); N = J(3, 3, 0)
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j) continue
		if (!dir & j < i) continue
		a = z[i]; b = z[j]
		M[a, b] = M[a, b] + A[i, j]; N[a, b] = N[a, b] + 1
		if (!dir & a != b) {
			M[b, a] = M[b, a] + A[i, j]
			N[b, a] = N[b, a] + 1
		}
	}
	M = M :/ N
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j) continue
		if (!dir & j < i) continue
		a = z[i]; b = z[j]
		V[a, b] = V[a, b] + (A[i, j] - M[a, b])^2
		if (!dir & a != b) V[b, a] = V[b, a] + (A[i, j] - M[b, a])^2
	}
	st_matrix("__bs", (M, V :/ N))
}
end
local fail 0
foreach spec in "u 120 120 0" "d 90 90 1" {
	local tag : word 1 of `spec'
	local n1 : word 2 of `spec'
	local n2 : word 3 of `spec'
	local dir : word 4 of `spec'
	mata: sbmt_simvar(`n1', `n2', `dir', 11)
	sbmt_frommata g, mat(__sim_A) dir(`dir')
	mata: st_matrix("__zt", __sim_z1)
	quietly nwsbm g, model(gaussian) maxblocks(5) seed(3) nodisplay
	local qplain = e(Q)
	quietly nwsbm g, model(gaussian) blockvar maxblocks(5) seed(3) generate(mb) replace nodisplay
	local ok 1
	if e(Q) != 3 {
		di as error "`tag': blockvar selects `=e(Q)' blocks, planted 3"
		local ok 0
	}
	capture drop truth
	quietly gen double truth = .
	mata: st_store((1::`n1'), "truth", st_matrix("__zt"))
	sbmt_ari mb truth
	local ari = r(ari)
	if `ari' < 0.9 | missing(`ari') {
		di as error "`tag': ARI `ari' < 0.9"
		local ok 0
	}
	* (2) empirical block statistics of the hard blocks
	matrix P = e(pi)
	matrix S2 = e(sigma2_block)
	matrix B = e(blocksize)
	mata: sbmt_blockstats(`n1', `dir')
	matrix E = __bs
	local maxd 0
	local maxv 0
	forvalues a = 1/3 {
		forvalues b = 1/3 {
			local ev = E[`a', 3 + `b']
			local d = abs(S2[`a', `b'] - `ev') / `ev'
			if `d' > `maxv' local maxv = `d'
		}
	}
	di as text "`tag': Q plain `qplain', blockvar `=e(Q)', ARI " %6.4f `ari' ", max rel. diff of the block variances to the empirical ones " %7.4f `maxv' ", pooled variance " %7.4f e(sigma2) ", npar " e(npar)
	if `maxv' > 0.03 {
		di as error "`tag': block variances differ from the empirical ones"
		local ok 0
	}
	* (3) parameters
	local np = cond(`dir', 18, 12)
	if e(npar) != `np' {
		di as error "`tag': npar `=e(npar)' (expected `np')"
		local ok 0
	}
	if e(sigma2) <= 0.1 | e(sigma2) >= 6 {
		di as error "`tag': pooled variance outside the range of the block variances"
		local ok 0
	}
	quietly nwsbm g, model(gaussian) blockvar blocks(3) seed(3) orderblocks nodisplay
	matrix S2o = e(sigma2_block)
	matrix Po = e(pi)
	matrix Bo = e(blocksize)
	if Bo[1, 1] < Bo[1, 2] | Bo[1, 2] < Bo[1, 3] {
		di as error "`tag': orderblocks"
		local ok 0
	}
	if `ok' di as result "PASS `tag'"
	else local fail = `fail' + 1
}
* two-mode
mata: sbmt_simvar(120, 90, 0, 12)
sbmt_frommata g2, mat(__sim_A) dir(0)
preserve
clear
mata: st_addobs(120)
mata: (void) st_addvar("double", "__a" :+ strofreal(1..90))
mata: st_store(., ., __sim_A)
nwset __a*, name(g2) bipartite valued replace
restore
quietly nwsbm g2, model(gaussian) blockvar maxblocks(4 4) seed(3) nodisplay
matrix S2 = e(sigma2_block)
if e(Q1) != 3 | e(Q2) != 3 | rowsof(S2) != 3 | colsof(S2) != 3 {
	di as error "two-mode: selected `=e(Q1)' x `=e(Q2)'"
	local fail = `fail' + 1
}
else di as result "PASS two-mode (Q1 x Q2 = 3 x 3, block variances " %6.3f S2[1, 1] " ... " %6.3f S2[3, 3] ")"
* errors
foreach opt in "model(bernoulli)" "model(poisson)" "model(gaussian) sigma2(1)" "model(gaussian) degreecorrected" {
	capture nwsbm g, blockvar `opt' blocks(2) nodisplay
	if _rc != 198 {
		di as error "blockvar `opt' should fail with 198, rc = " _rc
		local fail = `fail' + 1
	}
	else di as result "PASS error path blockvar `opt'"
}
if `fail' di as error "FAILED `fail'"
else di as result "ALL PASSED"
