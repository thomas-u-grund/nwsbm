* Default number of random starts. The plain one-mode Bernoulli model with up to 500 actors uses 2 extra random starts, every other model the old default (10 up to 500 actors).
* Agreement with restarts(10) on simulated networks (n = 100, directed and undirected, planted 3 and 5 blocks, three signal strengths, 5 seeds):
* the same number of blocks in at least 95% of the networks and an ICL that is never lower than that of restarts(10) by more than 1 (the validation grid in dev/restarts found 0.56 as the largest loss in 960 networks).
do cscripts/nwsbm_testlib.do
mata:
mata set matastrict off
void rst_sim(real scalar n, real scalar K, real scalar pin, real scalar pout, real scalar dir)
{
	real colvector z
	real matrix P, A
	z = mod((1::n), K) :+ 1
	z = z[order(uniform(n, 1), 1)]
	P = (z * J(1, n, 1) :== J(n, 1, 1) * z') * (pin - pout) :+ pout
	A = (uniform(n, n) :< P)
	if (!dir) A = (lowertriangle(A, 0) + lowertriangle(A, 0)')
	_diag(A, 0)
	st_matrix("__rst_A", A)
}
end
local fail 0
local n 100
local nn 0
local nq 0
local worst 0
foreach dir in 0 1 {
	foreach K in 3 5 {
		foreach pp in "0.40 0.05" "0.25 0.10" "0.25 0.08" {
			tokenize "`pp'"
			forvalues rep = 1/5 {
				set seed `=1000 * `dir' + 100 * `K' + `rep' + 10 * real("`1'") * 100'
				mata: rst_sim(`n', `K', `1', `2', `dir')
				clear
				quietly set obs `n'
				forvalues j = 1/`n' {
					quietly gen double __a`j' = .
				}
				mata: st_store(., ., st_matrix("__rst_A"))
				if `dir' quietly nwset __a*, name(g) directed replace
				else quietly nwset __a*, name(g) undirected replace
				quietly nwsbm g, maxblocks(`=`K' + 3') seed(`rep') threads(1)
				local r_d = e(restarts)
				local q_d = e(Q)
				local i_d = e(ICL)
				quietly nwsbm g, maxblocks(`=`K' + 3') seed(`rep') threads(1) restarts(10)
				local nn = `nn' + 1
				if `q_d' == e(Q) local nq = `nq' + 1
				local worst = min(`worst', `i_d' - e(ICL))
				if `r_d' != 2 {
					di as error "default restarts is `r_d', expected 2"
					local fail = `fail' + 1
				}
			}
		}
	}
}
di as text "RESTARTS: `nn' networks, same Q in `nq', worst ICL difference (default - restarts(10)) " %8.4f `worst'
if `nq' < 0.95 * `nn' | `worst' < -1 {
	di as error "default restarts: agreement with restarts(10) too weak"
	local fail = `fail' + 1
}
* the old default stays for the other models and above 500 actors; restarts() overrides
foreach m in poisson gaussian {
	quietly nwsbm g, blocks(2) model(`m') seed(1) threads(1)
	if e(restarts) != 10 {
		di as error "model(`m'): default restarts `=e(restarts)', expected 10"
		local fail = `fail' + 1
	}
}
quietly nwsbm g, blocks(2) degreecorrected model(poisson) seed(1) threads(1)
if e(restarts) != 10 {
	di as error "degreecorrected: default restarts `=e(restarts)', expected 10"
	local fail = `fail' + 1
}
quietly nwsbm g, blocks(2) seed(1) threads(1) restarts(10)
if e(restarts) != 10 {
	di as error "restarts(10) not honoured"
	local fail = `fail' + 1
}
if `fail' == 0 di as result "PASS default number of restarts"
else exit 9
