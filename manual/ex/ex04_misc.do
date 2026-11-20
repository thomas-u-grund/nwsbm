* options: seeds and threads, errors, size and timing
clear all
set linesize 90
nwwebuse zachary, nwclear
log using out/misc_seed.txt, text replace
quietly nwsbm zachary_bin, maxblocks(5) seed(1) threads(1)
display "threads(1): ICL " %12.6f e(ICL) ", Q = " e(Q)
quietly nwsbm zachary_bin, maxblocks(5) seed(1) threads(4)
display "threads(4): ICL " %12.6f e(ICL) ", Q = " e(Q)
quietly nwsbm zachary_bin, maxblocks(5) seed(2) threads(4)
display "seed(2):    ICL " %12.6f e(ICL) ", Q = " e(Q) ", restarts used = " e(restarts)
log close
nwwebuse lazega, nwclear
nwgenerate w = lazega_adv + lazega_cow
log using out/misc_errors.txt, text replace
capture noisily nwsbm zachary_bin, model(gaussian) degreecorrected
capture noisily nwsbm zachary_bin, greedy
capture noisily nwsbm zachary_bin, blocks(2) maxblocks(4)
capture noisily nwsbm zachary_bin, blocks(2 2)
capture noisily nwsbm zachary_bin, icl(exact) model(poisson)
capture noisily nwsbm w, blocks(2)
capture noisily nwsbm w, model(poisson) blocks(2) covariates(lazega_fr) degreecorrected
capture noisily nwsbm w, model(poisson) blocks(2) icl(exact) greedy fast
log close
log using out/misc_timing.txt, text replace
foreach n in 200 500 {
	nwrandom `n', prob(0.05) name(g`n') undirected
	timer clear 1
	timer on 1
	quietly nwsbm g`n', maxblocks(6) seed(1) threads(1)
	timer off 1
	quietly timer list 1
	display "n = `n', maxblocks(6), threads(1):   " %6.2f r(t1) " seconds"
	timer clear 2
	timer on 2
	quietly nwsbm g`n', maxblocks(6) seed(1) threads(1) fast
	timer off 2
	quietly timer list 2
	display "n = `n', maxblocks(6), threads(1), fast: " %6.2f r(t2) " seconds"
	timer clear 3
	timer on 3
	quietly nwsbm g`n', maxblocks(6) seed(1)
	timer off 3
	quietly timer list 3
	display "n = `n', maxblocks(6), all cores:    " %6.2f r(t3) " seconds"
}
log close
