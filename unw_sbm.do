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

real scalar SBM_PROTOCOL_NEEDED() return(1)

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

// y: n x n (0/1, missing = .); returns the stacked output of the plugin
real colvector SbmNat::run(real matrix y, real scalar directed, real scalar qmin, real scalar qmax, real scalar seed, real scalar restarts)
{
	real scalar n, tot, k, need
	real colvector blob

	n = rows(y)
	blob = (SBM_PROTOCOL_NEEDED() \ n \ directed \ 0 \ qmin \ qmax \ seed \ restarts \ J(8, 1, 0)) \ vec(editmissing(y, -1)')
	L = rows(blob)
	tot = 0
	for (k = qmin; k <= qmax; k++) tot = tot + 5 + n * k + k * k + k
	need = max((L, tot, 4))
	if (need > st_nobs()) st_addobs(need - st_nobs())
	st_store((1::L), 1, blob)
	stata("plugin call sbmnative v1 v2, " + char(34) + "10 " + strofreal(L) + char(34))
	return(st_data((1::tot), 2))
}

// adjacency matrix with the diagonal removed; an undirected network is symmetrised by the upper triangle
real matrix sbm_ymat(real matrix A, real scalar directed)
{
	real matrix y

	y = A
	_diag(y, 0)
	if (!directed) y = uppertriangle(y, 0) + uppertriangle(y, 0)'
	return(y)
}

// Fits Qmin..Qmax; leaves the results in the Mata globals __nwsbm_out (the plugin's output) and sets Stata matrices:
//   tabname (nq x 6: Q, ICL, J, PL, H, npar), and for the selected Q: pi, alpha, blocksize; the tau matrix stays in __nwsbm_tau
// returns the selected Q
real scalar SbmFit(real matrix y, real scalar directed, real scalar qmin, real scalar qmax, real scalar seed, real scalar restarts,
	string scalar tabname, string scalar piname, string scalar alphaname, string scalar sizename)
{
	class SbmNat scalar nat
	real colvector out
	real scalar n, k, off, nq, best, i, qq
	real matrix tab
	real matrix pm
	real rowvector hs, mx

	n = rows(y)
	nat.open()
	out = nat.run(y, directed, qmin, qmax, seed, restarts)
	nat.close()
	nq = qmax - qmin + 1
	tab = J(nq, 6, .)
	off = 0
	external real matrix __nwsbm_tau
	best = 1
	for (k = qmin; k <= qmax; k++) {
		i = k - qmin + 1
		tab[i, .] = (k, out[off + 1], out[off + 2], out[off + 3], out[off + 4], out[off + 5])
		off = off + 5 + n * k + k * k + k
	}
	best = 1
	for (i = 2; i <= nq; i++) if (tab[i, 2] > tab[best, 2] + 1e-9) best = i
	qq = qmin + best - 1
	off = 0
	for (k = qmin; k < qq; k++) off = off + 5 + n * k + k * k + k
	off = off + 5
	__nwsbm_tau = colshape(out[(off + 1)::(off + n * qq)], qq)
	off = off + n * qq
	pm = colshape(out[(off + 1)::(off + qq * qq)], qq)
	off = off + qq * qq
	st_matrix(tabname, tab)
	st_matrix(piname, pm)
	st_matrix(alphaname, out[(off + 1)::(off + qq)]')
	hs = J(1, qq, 0)
	for (k = 1; k <= n; k++) {
		mx = __nwsbm_tau[k, .]
		i = (select(1..qq, mx :== max(mx)))[1]
		hs[i] = hs[i] + 1
	}
	st_matrix(sizename, hs)
	return(qq)
}

// stores the hard memberships (and optionally the probabilities) of the last fit in variables; rows 1..n
void SbmStore(string scalar stub, real scalar probs)
{
	external real matrix __nwsbm_tau
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
	st_store((1::n), stub, hard)
	if (probs) for (k = 1; k <= q; k++) st_store((1::n), stub + "_p" + strofreal(k), __nwsbm_tau[., k])
}
end
