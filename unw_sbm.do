/*
 * Copyright (C) 2026 Thomas U. Grund
 * This file is part of nwsbm, a Stata package for stochastic block models. It follows the variational EM of the R packages
 * blockmodels and sbm and is written independently.
 * nwsbm is free software: you can redistribute it and/or modify it under the terms of the GNU General
 * Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.
 */
/*
	unw_sbm.do -- Mata glue of nwsbm (built into lib/lnwsbm.mlib). The estimation runs in the compiled plugin native/sbm_vem.c;
	Mata writes the adjacency matrix into the plugin's frame and unpacks the results. There is no Mata fallback.
*/
capture mata: mata drop SbmPlatform()
capture mata: mata drop SbmPluginFilename()
capture mata: mata drop SbmPluginPath()
capture mata: mata drop SbmNativeAvailable()
capture mata: mata drop SBM_PROTOCOL_NEEDED()
capture mata: mata drop SbmNat()
capture mata: mata drop SbmFit()
capture mata: mata drop sbm_ymat()
capture mata: mata drop SbmStore()
capture mata: mata drop sbm_reclen()
capture mata: mata drop SbmPlotMatrix()
capture mata: mata drop sbm_hard()
capture mata: mata drop sbm_setidx()
capture mata: mata drop sbm_checkvalues()

mata:
mata set matastrict off

string scalar SbmPlatform()
{
	string scalar os

	os = st_global("c(os)")
	if (os == "Windows") return("windows")
	if (os == "MacOSX") return("macos")
	if (substr(st_global("c(machine_type)"), 1, 3) == "Mac") return("macos")
	return("unix")
}

string scalar SbmPluginFilename() return("sbm_vem_" + SbmPlatform() + ".plugin")

string scalar SbmPluginPath()
{
	string scalar fname, found, full, dir, fn

	fname = SbmPluginFilename()
	found = findfile(fname)
	if (found != "") return(found)
	full = findfile("nwsbm.ado")
	if (full == "") return("")
	pathsplit(full, dir, fn)
	return(pathjoin(pathjoin(dir, "lib"), pathjoin("plugins", pathjoin(SbmPlatform(), fname))))
}

real scalar SbmNativeAvailable()
{
	string scalar p

	p = SbmPluginPath()
	if (p == "") return(0)
	return(fileexists(p))
}

real scalar SBM_PROTOCOL_NEEDED() return(2)

class SbmNat {
	real scalar isopen, L
	string scalar origframe
	void open()
	void close()
	real colvector run()
}

void SbmNat::open()
{
	string scalar p, fname, q
	real scalar rc, ver

	isopen = 0
	fname = SbmPluginFilename()
	p = SbmPluginPath()
	if (p == "" | !fileexists(p)) {
		errprintf("{p 0 0 2}nwsbm: the compiled plugin {bf:%s} was not found. nwsbm fits every model with this plugin and has no Mata fallback; reinstall the package with {bf:net install nwsbm, replace}.{p_end}\n", fname)
		exit(601)
	}
	origframe = st_framecurrent()
	stata("capture frame drop __nwsbm_native")
	stata("frame create __nwsbm_native")
	st_framecurrent("__nwsbm_native")
	isopen = 1
	st_addobs(4)
	(void) st_addvar("double", "v1")
	(void) st_addvar("double", "v2")
	q = char(34)
	stata("capture program drop sbmnative")
	stata("capture program sbmnative, plugin using(" + q + p + q + ")")
	rc = _stata("plugin call sbmnative v1 v2, " + q + "99 0" + q, 1)
	ver = .
	if (rc == 0) ver = st_data(1, 2)
	if (rc != 0 | ver != SBM_PROTOCOL_NEEDED()) {
		close()
		errprintf("{p 0 0 2}nwsbm: the installed plugin {bf:%s} does not match this version of nwsbm (or it cannot be loaded). Reinstall the package with {bf:net install nwsbm, replace}.{p_end}\n", fname)
		exit(3495)
	}
}

void SbmNat::close()
{
	if (isopen) {
		st_framecurrent(origframe)
		stata("capture frame drop __nwsbm_native")
		isopen = 0
	}
}

// y: n x n (one-mode) or n1 x n2 (two-mode) matrix (missing = .); returns the stacked output of the plugin
real colvector SbmNat::run(real matrix y, real scalar directed, real scalar model, real scalar kind, real scalar dc, real scalar qmin, real scalar qmax,
	real scalar q2min, real scalar q2max, real scalar seed, real scalar restarts, real scalar s2fix, real scalar threads, real scalar fast, real scalar tot)
{
	real scalar n, need
	real colvector blob

	n = rows(y)
	blob = (SBM_PROTOCOL_NEEDED() \ n \ directed \ model \ qmin \ qmax \ seed \ restarts \ s2fix \ threads \ kind \ cols(y) \ q2min \ q2max \ dc \ fast) \ vec(y')
	L = rows(blob)
	need = max((L, tot, 4))
	if (need > st_nobs()) st_addobs(need - st_nobs())
	st_store((1::L), 1, blob)
	stata("plugin call sbmnative v1 v2, " + char(34) + "10 " + strofreal(L) + char(34))
	return(st_data((1::tot), 2))
}

// adjacency matrix of a one-mode network with the diagonal removed; an undirected network is symmetrised by the upper triangle
real matrix sbm_ymat(real matrix A, real scalar directed)
{
	real matrix y

	y = A
	_diag(y, 0)
	if (!directed) y = uppertriangle(y, 0) + uppertriangle(y, 0)'
	return(y)
}

// positions of the mode-1 and mode-2 nodes of a two-mode network (globals __nwsbm_i1, __nwsbm_i2)
void sbm_setidx(string rowvector md)
{
	external real colvector __nwsbm_i1, __nwsbm_i2
	__nwsbm_i1 = select(1..cols(md), md :== "1")'
	__nwsbm_i2 = select(1..cols(md), md :== "2")'
}

// 0 = fine, 1 = Poisson needs non-negative integers, 2 = Bernoulli needs 0/1
real scalar sbm_checkvalues(real matrix y, real scalar model)
{
	real matrix o
	o = select(vec(y), vec(y) :< .)
	if (rows(o) == 0) return(0)
	if (model == 1) return((min(o) < 0 | any(o :!= floor(o))) ? 1 : 0)
	if (model == 0) return(any((o :!= 0) :& (o :!= 1)) ? 2 : 0)
	return(0)
}

// length of one plugin record
real scalar sbm_reclen(real scalar n1, real scalar n2, real scalar q1, real scalar q2)
{
	real scalar l
	l = 8 + n1 * q1 + q1 * (q2 > 0 ? q2 : q1) + q1
	if (q2 > 0) l = l + n2 * q2 + q2
	return(l)
}

// Fits Qmin..Qmax (one-mode) or the grid (Q1min..Q1max) x (Q2min..Q2max) (two-mode); sets Stata matrices
//   tabname (one row per fit: Q [Q2] ICL J PL H npar) and, for the selected fit, pi, alpha, blocksize (and alpha2 / blocksize2 for the columns of a
//   two-mode network); the membership probabilities stay in __nwsbm_tau (and __nwsbm_tau2); scalars __nwsbm_Q, __nwsbm_Q2, __nwsbm_sigma2
// returns the selected Q (rows block count for a two-mode network)
real scalar SbmFit(real matrix y, real scalar directed, real scalar model, real scalar kind, real scalar dc, real scalar qmin, real scalar qmax,
	real scalar q2min, real scalar q2max, real scalar seed, real scalar restarts, real scalar s2fix, real scalar threads, real scalar fast,
	string scalar tabname, string scalar piname, string scalar alphaname, string scalar sizename, string scalar alpha2name, string scalar size2name)
{
	class SbmNat scalar nat
	real colvector out
	real scalar n1, n2, k, l, off, nfit, best, i, qq, qq2, tot, bestoff
	real matrix tab, pm
	real rowvector hs, mx
	external real matrix __nwsbm_tau, __nwsbm_tau2

	n1 = rows(y)
	n2 = cols(y)
	if (!kind) {
		q2min = 0
		q2max = 0
	}
	tot = 0
	if (kind) {
		for (k = qmin; k <= qmax; k++) for (l = q2min; l <= q2max; l++) tot = tot + sbm_reclen(n1, n2, k, l)
	}
	else {
		for (k = qmin; k <= qmax; k++) tot = tot + sbm_reclen(n1, n1, k, 0)
	}
	nat.open()
	out = nat.run(y, directed, model, kind, dc, qmin, qmax, q2min, q2max, seed, restarts, s2fix, threads, fast, tot)
	nat.close()
	nfit = kind ? (qmax - qmin + 1) * (q2max - q2min + 1) : (qmax - qmin + 1)
	tab = J(nfit, kind ? 7 : 6, .)
	off = 0
	best = 1
	bestoff = 0
	for (i = 1; i <= nfit; i++) {
		k = out[off + 7]
		l = out[off + 8]
		if (kind) tab[i, .] = (k, l, out[off + 1], out[off + 2], out[off + 3], out[off + 4], out[off + 5])
		else tab[i, .] = (k, out[off + 1], out[off + 2], out[off + 3], out[off + 4], out[off + 5])
		if (i == 1 | out[off + 1] > tab[best, kind ? 3 : 2] + 1e-9) {
			best = i
			bestoff = off
		}
		off = off + sbm_reclen(n1, n2, k, l)
	}
	qq = out[bestoff + 7]
	qq2 = out[bestoff + 8]
	st_numscalar("__nwsbm_Q", qq)
	st_numscalar("__nwsbm_Q2", qq2)
	st_numscalar("__nwsbm_sigma2", out[bestoff + 6])
	off = bestoff + 8
	__nwsbm_tau = colshape(out[(off + 1)::(off + n1 * qq)], qq)
	off = off + n1 * qq
	if (kind) {
		__nwsbm_tau2 = colshape(out[(off + 1)::(off + n2 * qq2)], qq2)
		off = off + n2 * qq2
		pm = colshape(out[(off + 1)::(off + qq * qq2)], qq2)
		off = off + qq * qq2
	}
	else {
		__nwsbm_tau2 = J(0, 0, .)
		pm = colshape(out[(off + 1)::(off + qq * qq)], qq)
		off = off + qq * qq
	}
	st_matrix(tabname, tab)
	st_matrix(piname, pm)
	st_matrix(alphaname, out[(off + 1)::(off + qq)]')
	off = off + qq
	if (kind) st_matrix(alpha2name, out[(off + 1)::(off + qq2)]')
	hs = J(1, qq, 0)
	for (k = 1; k <= n1; k++) {
		mx = __nwsbm_tau[k, .]
		i = (select(1..qq, mx :== max(mx)))[1]
		hs[i] = hs[i] + 1
	}
	st_matrix(sizename, hs)
	if (kind) {
		hs = J(1, qq2, 0)
		for (k = 1; k <= n2; k++) {
			mx = __nwsbm_tau2[k, .]
			i = (select(1..qq2, mx :== max(mx)))[1]
			hs[i] = hs[i] + 1
		}
		st_matrix(size2name, hs)
	}
	return(qq)
}

// stores the hard memberships (and optionally the probabilities) of the last fit in variables. One-mode: rows 1..n. Two-mode: idx1 and idx2 give
// the rows (positions of the network's nodes) of the mode-1 and mode-2 nodes; the probabilities are those of the node's own mode.
void SbmStore(string scalar stub, real scalar probs, | real colvector idx1, real colvector idx2)
{
	external real matrix __nwsbm_tau, __nwsbm_tau2
	real scalar n, q, k
	real colvector hard
	real rowvector mx

	n = rows(__nwsbm_tau)
	q = cols(__nwsbm_tau)
	hard = J(n, 1, 1)
	for (k = 1; k <= n; k++) {
		mx = __nwsbm_tau[k, .]
		hard[k] = (select(1..q, mx :== max(mx)))[1]
	}
	if (args() < 3) {
		st_store((1::n), stub, hard)
		if (probs) for (k = 1; k <= q; k++) st_store((1::n), stub + "_p" + strofreal(k), __nwsbm_tau[., k])
		return
	}
	st_store(idx1, stub, hard)
	if (probs) for (k = 1; k <= q; k++) st_store(idx1, stub + "_p" + strofreal(k), __nwsbm_tau[., k])
	n = rows(__nwsbm_tau2)
	q = cols(__nwsbm_tau2)
	hard = J(n, 1, 1)
	for (k = 1; k <= n; k++) {
		mx = __nwsbm_tau2[k, .]
		hard[k] = (select(1..q, mx :== max(mx)))[1]
	}
	st_store(idx2, stub, hard)
	if (probs) for (k = 1; k <= q; k++) st_store(idx2, stub + "_p" + strofreal(k), __nwsbm_tau2[., k])
}

// hard (MAP) block of every row of a tau matrix
real colvector sbm_hard(real matrix tau)
{
	real colvector h
	real scalar k
	real rowvector mx

	h = J(rows(tau), 1, 1)
	for (k = 1; k <= rows(tau); k++) {
		mx = tau[k, .]
		h[k] = (select(1..cols(tau), mx :== max(mx)))[1]
	}
	return(h)
}

// data of the blocked adjacency matrix plot: writes the variables of the current dataset (ties: __x __y __v; shaded block rectangles: __yl __yh __xc __w
// __sh) and the axis specifications in locals. A: the matrix (one-mode: n x n, two-mode: rows i1, columns i2 of the full matrix); pi: connectivity
void SbmPlotMatrix(real matrix A, real scalar kind, real colvector i1, real colvector i2, real matrix pim, real scalar shade, real scalar showmiss)
{
	external real matrix __nwsbm_tau, __nwsbm_tau2
	real colvector h1, h2, p1, p2, c1, c2, sz1, sz2, ctr1, ctr2
	real matrix M, T
	real scalar nr, nc, q1, q2, k, l, nt, mxp, lev, r, cc
	real colvector i, j, idx, vv
	real matrix Msk
	string scalar sx, sy, lx, ly

	h1 = sbm_hard(__nwsbm_tau)
	q1 = cols(__nwsbm_tau)
	if (kind) {
		h2 = sbm_hard(__nwsbm_tau2)
		q2 = cols(__nwsbm_tau2)
		M = A
	}
	else {
		h2 = h1
		q2 = q1
		M = A
	}
	p1 = order((h1, (1::rows(h1))), (1, 2))
	p2 = order((h2, (1::rows(h2))), (1, 2))
	M = M[p1, p2]
	h1 = h1[p1]
	h2 = h2[p2]
	nr = rows(M)
	nc = cols(M)
	sz1 = J(q1, 1, 0)
	sz2 = J(q2, 1, 0)
	for (k = 1; k <= q1; k++) sz1[k] = sum(h1 :== k)
	for (k = 1; k <= q2; k++) sz2[k] = sz2[k] + sum(h2 :== k)
	c1 = runningsum(sz1)
	c2 = runningsum(sz2)
	ctr1 = c1 - sz1 :/ 2 :+ 0.5
	ctr2 = c2 - sz2 :/ 2 :+ 0.5
	// ties: observed nonzero cells (missing cells with showmiss as value .)
	Msk = (M :!= 0)
	if (!showmiss) Msk = Msk :& (M :< .)
	if (!kind) _diag(Msk, 0)
	idx = selectindex(vec(Msk))
	vv = vec(M)[idx]
	j = ceil(idx :/ nr)
	i = idx - (j :- 1) :* nr
	T = (j, nr :- i :+ 1, vv)
	nt = rows(T)
	// rectangles of the block pairs, shaded by the relative connectivity
	r = 0
	mxp = max(abs(pim))
	if (mxp <= 0) mxp = 1
	if (shade) r = q1 * q2
	st_addobs(nt + r)
	(void) st_addvar("double", ("__x", "__y", "__v", "__yl", "__yh", "__xc", "__w", "__sh", "__rt"))
	if (nt > 0) st_store((1::nt), ("__x", "__y", "__v"), T)
	st_store((1::(nt + r)), "__rt", (J(nt, 1, 0) \ J(r, 1, 1)))
	if (shade) {
		cc = nt
		for (k = 1; k <= q1; k++) {
			for (l = 1; l <= q2; l++) {
				cc++
				lev = ceil(4 * abs(pim[k, l]) / mxp)
				if (lev < 1) lev = 1
				if (lev > 4) lev = 4
				st_store(cc, ("__yl", "__yh", "__xc", "__w", "__sh"), (nr - c1[k] + 0.5, nr - (c1[k] - sz1[k]) + 0.5, c2[l] - sz2[l] / 2 + 0.5, sz2[l], lev))
			}
		}
	}
	sx = ""
	sy = ""
	lx = ""
	ly = ""
	for (k = 1; k < q2; k++) sx = sx + " " + strofreal(c2[k] + 0.5)
	for (k = 1; k < q1; k++) sy = sy + " " + strofreal(nr - c1[k] + 0.5)
	for (k = 1; k <= q2; k++) lx = lx + " " + strofreal(ctr2[k]) + " " + char(34) + strofreal(k) + char(34)
	for (k = 1; k <= q1; k++) ly = ly + " " + strofreal(nr - ctr1[k] + 1) + " " + char(34) + strofreal(k) + char(34)
	st_local("__xl", sx)
	st_local("__yl", sy)
	st_local("__xlab", lx)
	st_local("__ylab", ly)
	st_local("__nr", strofreal(nr))
	st_local("__nc", strofreal(nc))
}
end
