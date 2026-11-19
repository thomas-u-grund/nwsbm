* Recovery of planted partitions by the new models: Poisson, Gaussian, degree-corrected Poisson (heterogeneous degrees), two-mode latent block model.
* For each setting: the selected number of blocks equals the planted one and the ARI of the memberships against the planted blocks is at least the stated minimum.
do cscripts/nwsbm_testlib.do
mata:
mata set matastrict off
// simulate: kind 0 one-mode undirected / 1 two-mode; model 0 bernoulli 1 poisson 2 gaussian; Pm Qr x Qc mean matrix; dcsd: lognormal sd of degree propensities
void sbmr_sim(real scalar kind, real scalar n1, real scalar n2, real matrix Pm, real scalar model, real scalar sd, real scalar dcsd, real scalar seed)
{
	real scalar i, j, qr, qc
	real colvector z1, z2, th
	real matrix A, M
	rseed(seed)
	qr = rows(Pm)
	qc = cols(Pm)
	z1 = mod((1::n1) :+ floor(runiform(1, 1) * 0), qr) :+ 1
	z1 = jumble(z1)
	if (kind) {
		z2 = jumble(mod(1::n2, qc) :+ 1)
	}
	else {
		z2 = z1
		n2 = n1
	}
	th = exp(dcsd :* rnormal(n1, 1, 0, 1))
	th = th :/ mean(th)
	M = Pm[z1, z2]
	if (dcsd > 0 & !kind) M = M :* (th * th')
	A = J(n1, n2, 0)
	for (i = 1; i <= n1; i++) {
		for (j = 1; j <= n2; j++) {
			if (!kind & j <= i) continue
			if (model == 0) A[i, j] = (runiform(1, 1) < M[i, j])
			else if (model == 1) A[i, j] = rpoisson(1, 1, M[i, j])
			else A[i, j] = M[i, j] + sd * rnormal(1, 1, 0, 1)
		}
	}
	if (!kind) A = A + A'
	st_matrix("__sbmr_A", A)
	st_matrix("__sbmr_z1", z1)
	st_matrix("__sbmr_z2", z2)
}
end

capture program drop sbmr_net
program sbmr_net
	syntax name(name=nm), KIND(integer) N1(integer) N2(integer) MODEL(integer) SD(real) DCSD(real) SEED(integer) PM(string)
	tempname P
	matrix `P' = `pm'
	mata: sbmr_sim(`kind', `n1', `n2', st_matrix("`P'"), `model', `sd', `dcsd', `seed')
	preserve
	clear
	local nc = cond(`kind', `n2', `n1')
	set obs `n1'
	forvalues j = 1/`nc' {
		quietly gen double __a`j' = .
	}
	mata: st_store(., ., st_matrix("__sbmr_A"))
	if `kind' nwset __a*, name(`nm') bipartite valued replace
	else nwset __a*, name(`nm') undirected valued replace
	restore
end

local fail 0
capture program drop sbmr_check
program sbmr_check
	args label q ari minari wantq
	if `q' != `wantq' | `ari' < `minari' {
		di as error "FAIL `label': Q `q' (planted `wantq'), ARI " %6.4f `ari' " (minimum `minari')"
		global SBM_RFAILS = ${SBM_RFAILS} + 1
	}
	else di as result "PASS `label': Q `q', ARI " %6.4f `ari'
end
global SBM_RFAILS 0

* Poisson, 3 blocks, undirected
sbmr_net p1, kind(0) n1(120) n2(120) model(1) sd(0) dcsd(0) seed(11) pm(3,0.5,0.4\0.5,2.5,0.5\0.4,0.5,2)
nwsbm p1, model(poisson) maxblocks(5) seed(1) generate(mb) replace nodisplay
local q = e(Q)
gen double zt = .
mata: st_store((1::120), "zt", st_matrix("__sbmr_z1"))
sbmt_ari mb zt
sbmr_check "Poisson 3 blocks" `q' `r(ari)' 0.95 3
drop mb* zt
* Gaussian, 2 blocks, directed-like strength, unit variance
sbmr_net g1, kind(0) n1(100) n2(100) model(2) sd(1) dcsd(0) seed(12) pm(1,0\0,1)
matrix PM = (1.2, 0 \ 0, 1.2)
sbmr_net g1, kind(0) n1(100) n2(100) model(2) sd(1) dcsd(0) seed(12) pm(1.2,0\0,1.2)
nwsbm g1, model(gaussian) maxblocks(4) seed(1) generate(mb) replace nodisplay
local q = e(Q)
gen double zt = .
mata: st_store((1::100), "zt", st_matrix("__sbmr_z1"))
sbmt_ari mb zt
sbmr_check "Gaussian 2 blocks" `q' `r(ari)' 0.95 2
drop mb* zt
* degree corrected: heterogeneous degrees (lognormal sd 0.8), 3 blocks; the plain Poisson SBM over-splits
sbmr_net d1, kind(0) n1(150) n2(150) model(1) sd(0) dcsd(0.8) seed(13) pm(0.30,0.03,0.03\0.03,0.26,0.03\0.03,0.03,0.22)
nwsbm d1, model(poisson) degreecorrected maxblocks(6) seed(1) generate(mb) replace nodisplay
local q = e(Q)
gen double zt = .
mata: st_store((1::150), "zt", st_matrix("__sbmr_z1"))
sbmt_ari mb zt
sbmr_check "degree-corrected Poisson 3 blocks" `q' `r(ari)' 0.9 3
nwsbm d1, model(poisson) maxblocks(6) seed(1) nodisplay
if e(Q) <= 3 {
	di as error "FAIL plain Poisson SBM expected to select more than 3 blocks on heterogeneous degrees, selected `=e(Q)'"
	global SBM_RFAILS = ${SBM_RFAILS} + 1
}
else di as result "PASS plain Poisson SBM over-splits (selects `=e(Q)' blocks) where the degree-corrected one selects 3"
drop mb* zt
* two-mode: 60 x 90, 2 x 3 blocks, Bernoulli and Poisson and Gaussian
foreach mod in 0 1 2 {
	local nm : word `=`mod'+1' of bernoulli poisson gaussian
	if `mod' == 0 local pm "0.55,0.08,0.15\0.1,0.5,0.08"
	if `mod' == 1 local pm "4,0.5,1.2\0.7,3.5,0.6"
	if `mod' == 2 local pm "2,0,0.8\0.3,1.8,0"
	sbmr_net l`mod', kind(1) n1(60) n2(90) model(`mod') sd(1) dcsd(0) seed(`=20+`mod'') pm(`pm')
	nwsbm l`mod', model(`nm') maxblocks(4 4) seed(1) generate(mb) replace nodisplay
	* mode 1 of the nwset network = columns (90), mode 2 = rows (60)
	local q1 = e(Q1)
	local q2 = e(Q2)
	_nwsyntax l`mod', max(1)
	mata: sbm_setidx(`netobj'->get_modes())
	gen double zt = .
	mata: st_store(__nwsbm_i1, "zt", st_matrix("__sbmr_z2"))
	mata: st_store(__nwsbm_i2, "zt", st_matrix("__sbmr_z1"))
	* ARI separately for each mode
	mata: st_numscalar("__a1", sbmt_ari(st_data(__nwsbm_i1, "mb"), st_data(__nwsbm_i1, "zt")))
	mata: st_numscalar("__a2", sbmt_ari(st_data(__nwsbm_i2, "mb"), st_data(__nwsbm_i2, "zt")))
	if `q1' != 3 | `q2' != 2 | scalar(__a1) < 0.9 | scalar(__a2) < 0.9 {
		di as error "FAIL two-mode `nm': blocks `q1' x `q2' (planted 3 x 2 for columns x rows), ARI " %6.4f scalar(__a1) " / " %6.4f scalar(__a2)
		global SBM_RFAILS = ${SBM_RFAILS} + 1
	}
	else di as result "PASS two-mode `nm': blocks `q1' x `q2', ARI " %6.4f scalar(__a1) " / " %6.4f scalar(__a2)
	drop mb* zt
}
if ${SBM_RFAILS} == 0 di as result "PASS recovery of planted partitions (Poisson, Gaussian, degree-corrected, two-mode)"
else {
	di as error "FAIL: ${SBM_RFAILS} settings"
	exit 9
}
