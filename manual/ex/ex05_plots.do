* the figures of chapter 7 (PDF files in ../fig)
clear all
set linesize 90
set scheme s1network
nwwebuse zachary, nwclear
quietly nwsbm zachary_bin, maxblocks(5) seed(20261120)
nwsbmplot, type(matrix) shade blocklabels(periphery core) name(f_matrix, replace) export(../fig/fig_matrix.pdf) xsize(5.5) ysize(5.5)
nwsbmplot, type(icl) name(f_icl, replace) export(../fig/fig_icl.pdf) xsize(5.5) ysize(4)
nwsbmplot, type(membership) blocklabels(periphery core) title("") colors(gray) name(f_member, replace) export(../fig/fig_member.pdf) xsize(5.5) ysize(3.5)
quietly nwsbm zachary_bin, blocks(3) seed(20261120)
nwsbmplot, type(blocks) colors(gray) subtitle("node size: block size; line width: connectivity") name(f_blocks, replace) export(../fig/fig_blocks.pdf) xsize(5.5) ysize(4.5)
nwwebuse lazega, nwclear
quietly nwsbm lazega_fr, maxblocks(6) seed(20261120) icl(exact)
nwsbmplot, type(expected) sort name(f_expected, replace) export(../fig/fig_expected.pdf) xsize(5.5) ysize(5.5)
nwsbmplot, type(icl) name(f_iclexact, replace) export(../fig/fig_iclexact.pdf) xsize(5.5) ysize(4)
nwwebuse davis, nwclear
quietly nwsbm davis, maxblocks(4 4) seed(20261120)
nwsbmplot, type(matrix) name(f_davis, replace) export(../fig/fig_davis.pdf) xsize(5.5) ysize(4.5)
nwsbmplot, type(icl) colors(gray) name(f_davisicl, replace) export(../fig/fig_davisicl.pdf) xsize(5.5) ysize(4)
