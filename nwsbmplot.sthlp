{smcl}
{* *! version 1.1.0  23nov2026 author: Thomas Grund}{...}
{title:Title}

{p2colset 9 20 27 2}{...}
{p2col :nwsbmplot  {hline 2}}Plots after nwsbm: blocked matrix, expected ties, memberships, block graph, ICL path{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:nwsbmplot} [{cmd:,} {it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt ty:pe(type)}}{cmd:matrix} (default, or {cmd:data}), {cmd:expected}, {cmd:membership}, {cmd:blocks} (or {cmd:meso}) or {cmd:icl}{p_end}
{synopt :{opt shade}}matrix: shade the block rectangles by their connectivity (four levels, relative to the largest entry of {cmd:e(pi)}){p_end}
{synopt :{opt nolin:es}}matrix, expected, membership: no lines between the blocks{p_end}
{synopt :{opt miss:ing}}matrix: show pairs with a missing tie{p_end}
{synopt :{opt nolab:els}}matrix, expected: no block numbers on the axes{p_end}
{synopt :{opt blocklabels(list)}}names of the blocks (mode 1 first, then mode 2; as many as there are blocks) for the axes of matrix and expected, the legend of membership and the nodes of blocks; not used by icl{p_end}
{synopt :{opt actorlabels(varname)}}matrix, expected, membership: names of the actors (a string variable, or a numeric one with value labels, one observation per actor of the network) on the axes or under the bars; only for up to 40 actors, otherwise ignored with a note; on the axes they replace the block names{p_end}
{synopt :{opt linec:olor(color)}}matrix, expected: colour of the lines between the blocks; default {cmd:gs8}{p_end}
{synopt :{opt nolegend}}icl, membership: no legend box{p_end}
{synopt :{opt sort}}matrix, expected: within a block order the actors by decreasing certainty (largest membership probability){p_end}
{synopt :{opt col:ors(list)}}colours of the blocks (blocks, membership, icl): {cmd:default}, {cmd:gray}, or a list of Stata colour names, used in turn{p_end}
{synopt :{opt title(text)}}title of the graph; {cmd:title("")} draws no title{p_end}
{synopt :{opt name(name)}}name of the graph in memory{p_end}
{synopt :{opt sav:ing(filename)}}save the graph as a {cmd:.gph} file{p_end}
{synopt :{opt exp:ort(filename)}}export the graph with {helpb graph export} ({cmd:.png}, {cmd:.pdf}, {cmd:.eps}, ... by the extension; an existing file is replaced){p_end}
{synopt :{it:twoway_options}}other options of {helpb twoway}, such as {cmd:xtitle()} or {cmd:ytitle()}; they come last and override the defaults. Titles and axis titles are black by default; the axes of matrix and expected are titled "block of sender" and "block of receiver" ("block of mode 1" and "block of mode 2" for a two-mode network){p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:nwsbmplot} uses the last {helpb nwsbm} fit (the network and the membership probabilities must still be in memory). The {cmd:matrix} plot draws every tie of the network as a point in the adjacency matrix, with the actors ordered by their most probable block, so that the block structure shows as dense and sparse rectangles; for a valued network the points are darker for larger values. For a two-mode network the rows are the actors of mode 1 and the columns those of mode 2, each ordered by its blocks.

{pstd}
The {cmd:expected} plot draws the posterior expected value of every tie, E = T1 {it:pi} T2', where {it:pi} is the connectivity matrix {cmd:e(pi)} and T1, T2 are the membership probabilities of the actors (of the two modes for a two-mode network); the actors are ordered as in the {cmd:matrix} plot and eight grey levels span the range of the values (the range is in the subtitle). It shows the block part of the model: degree parameters of {cmd:degreecorrected} and the covariate terms are not included. One-mode networks: the diagonal is left out. At most 2.5 million cells are drawn.

{pstd}
The {cmd:membership} plot draws the membership probabilities of the actors as stacked bars (one colour per block), the actors ordered by their most probable block and within it by decreasing certainty; blocks that are not clear-cut show as mixed bars. Two-mode networks: one panel per mode, combined.

{pstd}
The {cmd:blocks} plot draws the blocks as nodes (size proportional to the number of actors) and the connectivity between blocks as lines whose width follows the entries of {cmd:e(pi)}, relative to the largest entry; entries below 2 percent of it are not drawn. One-mode blocks sit on a circle; for a directed network the two directions of a pair of blocks are drawn as two arrows side by side; the blocks of a two-mode network sit in two columns. Connectivity within a block is not drawn. With {opt colors()} the nodes are coloured by block.

{pstd}
The {cmd:icl} plot draws the ICL of every fitted number of blocks and marks the selected one (a diamond). With {cmd:icl(exact)} it draws the exact ICL and the ICL of {cmd:blockmodels} ({cmd:e(ICLbm)}, right axis). For a two-mode network there is one line per number of blocks of mode 2 over the number of blocks of mode 1. It needs at least two fitted numbers of blocks ({cmd:maxblocks()}) and works without the network or the memberships in memory.

{title:Examples}

{phang2}{cmd:. nwsbm g, maxblocks(6) seed(1)}{p_end}
{phang2}{cmd:. nwsbmplot, shade}{p_end}
{phang2}{cmd:. nwsbmplot, type(blocks) colors(navy maroon forest_green)}{p_end}
{phang2}{cmd:. nwsbmplot, type(expected) sort}{p_end}
{phang2}{cmd:. nwsbmplot, type(membership)}{p_end}
{phang2}{cmd:. nwsbmplot, type(icl) export(icl.png)}{p_end}

{title:Also see}

{pstd}
{help nwsbm}, {help nwplotmatrix}

{hline}
{pstd}{it:nwsbm version 1.1.0, 23 November 2026}{p_end}{* nwversion}
