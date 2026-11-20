* Davis, Gardner and Gardner (1941): 18 Southern women and 14 social events, a two-mode network
clear all
set linesize 90
log using out/davis_data.txt, text replace
nwwebuse davis, nwclear
nwsummarize davis
log close
log using out/davis_plain.txt, text replace
nwsbm davis, maxblocks(4 4) seed(20261120) generate(d0) replace
log close
log using out/davis_dc.txt, text replace
nwsbm davis, maxblocks(4 4) seed(20261120) model(poisson) degreecorrected icl(exact) generate(d1) replace nodisplay
estat icl
nwsbm davis, blocks(2 2) seed(20261120) model(poisson) degreecorrected icl(exact) greedy orderblocks generate(d3) replace
estat blocks
log close
