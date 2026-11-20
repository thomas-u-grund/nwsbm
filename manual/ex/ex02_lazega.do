* Lazega's law firm (71 lawyers, directed friendship network), loaded from the nwcommands website
clear all
set linesize 90
log using out/lazega_data.txt, text replace
nwwebuse lazega, nwclear
nwsummarize lazega_fr
describe, short
log close
log using out/lazega_plain.txt, text replace
nwsbm lazega_fr, maxblocks(6) seed(20261120) generate(b0) replace
log close
log using out/lazega_cov.txt, text replace
nwsbm lazega_fr, maxblocks(4) seed(20261120) same(office practice gender) absdiff(seniority) generate(b1) replace
matrix list e(beta)
matrix list e(se_beta)
display "z of same_office: " %6.2f el(e(beta), 1, colnumb(e(beta), "same_office")) / el(e(se_beta), 1, colnumb(e(beta), "same_office"))
display "odds ratio of same_office: " %6.2f exp(el(e(beta), 1, colnumb(e(beta), "same_office")))
log close
log using out/lazega_exact.txt, text replace
nwsbm lazega_fr, maxblocks(6) seed(20261120) icl(exact) generate(e1) replace nodisplay
estat icl
display "exact ICL: " %9.2f e(ICL) "   criterion of blockmodels at the same fit: " %9.2f e(ICLbm)
nwsbm lazega_fr, maxblocks(6) seed(20261120) icl(exact) greedy generate(e2) replace nodisplay
estat icl
tabulate e1 e2
log close
log using out/lazega_valued.txt, text replace
nwgenerate ties3 = lazega_adv + lazega_cow + lazega_fr
nwsummarize ties3
nwsbm ties3, model(poisson) maxblocks(5) seed(20261120) nodisplay
estat icl
nwsbm ties3, model(poisson) degreecorrected maxblocks(5) seed(20261120) nodisplay
estat icl
nwsbm ties3, model(gaussian) maxblocks(5) seed(20261120) generate(g1) replace nodisplay
display "one variance: ICL " %9.2f e(ICL) ", " e(Q) " blocks, variance " %6.3f e(sigma2)
nwsbm ties3, model(gaussian) blockvar maxblocks(5) seed(20261120) orderblocks generate(g2) replace nodisplay
matrix list e(sigma2_block)
display "pooled residual variance: " %6.3f e(sigma2) ", Q = " e(Q) ", ICL " %9.2f e(ICL)
tabulate g1 g2
log close
