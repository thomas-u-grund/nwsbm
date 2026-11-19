* Recovery of planted partitions: n = 100 actors, 3 blocks of 40/35/25, within-block tie probability pin, between-block 0.10,
* undirected and directed, 10 simulated networks per cell. Reports the mean adjusted Rand index (ARI) of the memberships against the planted
* blocks (with blocks(3)) and the share of replications in which ICL selects Q = 3 (maxblocks(6)).
* Assertions: pin = 0.50: mean ARI >= 0.95 and Q = 3 selected in >= 9 of 10; pin = 0.30: mean ARI >= 0.85; the mean ARI never decreases with pin by more than 0.05.
do cscripts/nwsbm_testlib.do
mata:
mata set matastrict off
void sbmt_sim(real scalar n, real scalar pin, real scalar pout, real scalar dir)
{
	real colvector z
	real matrix P, A
	real scalar i, j
	z = J(40, 1, 1) \ J(35, 1, 2) \ J(25, 1, 3)
	z = z[order(uniform(n, 1), 1)]
	A = J(n, n, 0)
	for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) {
		if (i == j) continue
		if (!dir & j < i) {
			A[i, j] = A[j, i]
			continue
		}
		A[i, j] = (runiform(1, 1) < (z[i] == z[j] ? pin : pout))
	}
	st_matrix("__sbmt_A", A)
	st_matrix("__sbmt_z", z)
}
end
local fail 0
local prev .
foreach dir in 0 1 {
	local prev = -1
	foreach pin in 0.15 0.20 0.30 0.50 {
		local sari 0
		local nq3 0
		local sarq 0
		forvalues rep = 1/10 {
			set seed `=1000 * `dir' + 100 * `pin' * 100 + `rep''
			mata: sbmt_sim(100, `pin', 0.10, `dir')
			clear
			quietly set obs 100
			forvalues j = 1/100 {
				quietly gen double __a`j' = .
			}
			mata: st_store(., ., st_matrix("__sbmt_A"))
			quietly gen truth = .
			mata: st_store(., "truth", st_matrix("__sbmt_z"))
			if `dir' nwset __a*, name(g) directed replace
			else nwset __a*, name(g) undirected replace
			quietly nwsbm g, blocks(3) seed(`rep') generate(m3)
			sbmt_ari m3 truth
			local sari = `sari' + r(ari)
			quietly nwsbm g, maxblocks(6) seed(`rep') generate(ms) replace
			if e(Q) == 3 local nq3 = `nq3' + 1
			sbmt_ari ms truth
			local sarq = `sarq' + r(ari)
		}
		local mari = `sari' / 10
		local marq = `sarq' / 10
		di as text "RECOVERY dir=`dir' pin=" %4.2f `pin' "  mean ARI (blocks(3)) = " %6.3f `mari' "  mean ARI (selected Q) = " %6.3f `marq' "  Q=3 selected in `nq3'/10"
		if `pin' == 0.5 & (`mari' < 0.95 | `nq3' < 9) {
			di as error "dir=`dir' pin=0.5: recovery too weak"
			local fail = `fail' + 1
		}
		if `pin' == 0.3 & `mari' < 0.85 {
			di as error "dir=`dir' pin=0.3: recovery too weak"
			local fail = `fail' + 1
		}
		if `mari' < `prev' - 0.05 {
			di as error "dir=`dir': ARI decreases with the strength"
			local fail = `fail' + 1
		}
		local prev = `mari'
	}
}
if `fail' == 0 di as result "PASS recovery of planted partitions"
else exit 9
