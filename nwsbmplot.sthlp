{smcl}
{* *! version 0.4.0  19nov2026 author: Thomas Grund}{...}
{title:Title}

{p2colset 9 20 27 2}{...}
{p2col :nwsbmplot  {hline 2}}Plots after nwsbm: blocked adjacency matrix and block graph{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}
{cmd:nwsbmplot} [{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt ty:pe(matrix|blocks)}}{cmd:matrix} (default): the adjacency matrix with rows and columns ordered by block; {cmd:blocks}: the block graph{p_end}
{synopt :{opt shade}}matrix: shade the block rectangles by their connectivity (four levels, relative to the largest entry of {cmd:e(pi)}){p_end}
{synopt :{opt nolines}}matrix: no lines between the blocks{p_end}
{synopt :{opt missing}}matrix: show pairs with a missing tie{p_end}
{synopt :{opt nolabels}}matrix: no block numbers on the axes{p_end}
{synopt :{opt title(text)}}title of the graph{p_end}
{synopt :{opt name(name)}}name of the graph in memory{p_end}
{synopt :{opt sav:ing(filename)}}save the graph as a {cmd:.gph} file{p_end}
{synopt :{it:twoway_options}}other options of {helpb twoway}{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:nwsbmplot} uses the last {helpb nwsbm} fit (the network and the membership probabilities must still be in memory). The {cmd:matrix} plot draws every tie of the network as a point in the adjacency matrix, with the actors ordered by their most probable block, so that the block structure shows as dense and sparse rectangles; for a valued network the points are darker for larger values. For a two-mode network the rows are the actors of mode 2 and the columns those of mode 1, each ordered by its blocks.

{pstd}
The {cmd:blocks} plot draws the blocks as nodes (size proportional to the number of actors) and the connectivity between blocks as lines whose width follows the entries of {cmd:e(pi)}, relative to the largest entry; entries below 2 percent of it are not drawn. One-mode blocks sit on a circle (directed networks have arrows); the blocks of a two-mode network sit in two columns.

{title:Examples}

{phang2}{cmd:. nwsbm g, maxblocks(6) seed(1)}{p_end}
{phang2}{cmd:. nwsbmplot, shade}{p_end}
{phang2}{cmd:. nwsbmplot, type(blocks)}{p_end}

{title:Also see}

{pstd}
{help nwsbm}, {help nwplotmatrix}

{hline}
{pstd}{it:nwsbm version 0.4.0, 19 November 2026}{p_end}{* nwversion}
