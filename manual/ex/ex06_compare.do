* comparing block models: ICL of the fits of one network; nwcompare lists the active fit
clear all
set linesize 90
nwwebuse lazega, nwclear
log using out/compare_icl.txt, text replace
quietly nwsbm lazega_fr, maxblocks(6) seed(20261120)
display "no covariates:       Q = " e(Q) ", ICL = " %9.2f e(ICL) ", parameters = " e(npar)
quietly nwsbm lazega_fr, maxblocks(6) seed(20261120) same(office practice gender) absdiff(seniority)
display "with covariates:     Q = " e(Q) ", ICL = " %9.2f e(ICL) ", parameters = " e(npar)
nwcompare .
log close
