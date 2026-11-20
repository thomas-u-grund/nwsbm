* Zachary's karate club (34 members, undirected): the data are loaded from the nwcommands website
clear all
set linesize 90
log using out/karate_data.txt, text replace
nwwebuse zachary, nwclear
nwsummarize zachary_bin
log close
log using out/karate_fit.txt, text replace
nwsbm zachary_bin, maxblocks(5) seed(20261120) generate(block) replace
log close
log using out/karate_estat.txt, text replace
estat icl
estat blocks
display "selected: " e(Q) " blocks, ICL " %9.3f e(ICL)
matrix list e(pi)
matrix list e(alpha)
log close
log using out/karate_cross.txt, text replace
tabulate block group
nwsbm zachary_bin, blocks(3) seed(20261120) generate(block3) replace nodisplay
estat blocks
tabulate block3 group
log close
log using out/karate_replay.txt, text replace
nwsbm zachary_bin, blocks(2) seed(20261120) nodisplay
nwsbm, generate(b2) replace
summarize b2 b2_p1 b2_p2
nwsbm version
log close
