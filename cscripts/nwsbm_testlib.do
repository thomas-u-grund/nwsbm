* Shared set-up of the nwsbm tests: scratch PLUS directory, ad-paths, networks from the R fixtures.
* Run the tests from the root of the repository. NWCORE: directory of the core package (default: the owner's checkout).
capture program drop nwinstall
local core : env NWCORE
if "`core'" == "" local core "/Users/tgrund/FILES_NEW/SOFTWARE/nwcommands"
local scratch : env NWSBM_SCRATCH
if "`scratch'" == "" local scratch "/private/tmp/claude-501/-Users-tgrund-FILES-NEW-SOFTWARE-nwcommands/9d5243e9-ee49-4558-835d-02925b410730/scratchpad/sbm/work"
local tag : env NWSBM_TAG
local plus "`scratch'/plus_`tag'_`=floor(runiform()*1e9)'"
capture mkdir "`scratch'"
capture mkdir "`plus'"
sysdir set PLUS "`plus'"
local repo "`c(pwd)'"
adopath ++ "`core'"
adopath ++ "`core'/lib"
adopath ++ "`repo'"
adopath ++ "`repo'/lib"
adopath ++ "`repo'/lib/plugins/macos"
global SBM_REPO "`repo'"
global SBM_SCRATCH "`scratch'"

capture program drop sbmt_mknet
* sbmt_mknet name, ref(fixture) n(#) dir(0|1): network of n actors from the fixture's edge list (and its missing list)
program sbmt_mknet
	syntax name(name=nm), REF(string) N(integer) DIR(integer)
	mata: sbmt_load(`n', `dir', "$SBM_REPO/dev/ref/`ref'_edges.csv", "$SBM_REPO/dev/ref/`ref'_missing.csv")
	preserve
	clear
	set obs `n'
	forvalues j = 1/`n' {
		quietly gen double __a`j' = .
	}
	mata: st_store(., ., st_matrix("__sbmt_A"))
	if `dir' nwset __a*, name(`nm') directed replace
	else nwset __a*, name(`nm') undirected replace
	restore
	capture matrix drop __sbmt_A
end

capture program drop sbmt_ari
* sbmt_ari var1 var2: adjusted Rand index of two partitions (r(ari))
program sbmt_ari, rclass
	args v1 v2
	mata: st_numscalar("__sbmt_ari", sbmt_ari(st_data(., "`v1'"), st_data(., "`v2'")))
	return scalar ari = scalar(__sbmt_ari)
end

mata:
mata set matastrict off
void sbmt_load(real scalar n, real scalar dir, string scalar fe, string scalar fm)
{
	real matrix A, E
	real scalar k
	A = J(n, n, 0)
	E = sbmt_readcsv(fe)
	for (k = 1; k <= rows(E); k++) {
		A[E[k, 1], E[k, 2]] = 1
		if (!dir) A[E[k, 2], E[k, 1]] = 1
	}
	E = sbmt_readcsv(fm)
	for (k = 1; k <= rows(E); k++) {
		A[E[k, 1], E[k, 2]] = .
		if (!dir) A[E[k, 2], E[k, 1]] = .
	}
	_diag(A, 0)
	st_matrix("__sbmt_A", A)
}
real matrix sbmt_readcsv(string scalar f)
{
	real scalar fh, k
	string scalar ln
	string rowvector r
	real matrix E
	string matrix L
	fh = fopen(f, "r")
	L = J(0, 1, "")
	while ((ln = fget(fh)) != J(0, 0, "")) L = L \ ln
	fclose(fh)
	if (rows(L) <= 1) return(J(0, 2, .))
	E = J(rows(L) - 1, 2, .)
	for (k = 2; k <= rows(L); k++) {
		r = tokens(subinstr(L[k], ",", " "))
		E[k - 1, .] = strtoreal(r[1..2])
	}
	return(E)
}
real scalar sbmt_ari(real colvector a, real colvector b)
{
	real matrix T
	real scalar n, sij, sa, sb, e, mx
	real colvector ua, ub
	ua = uniqrows(a); ub = uniqrows(b)
	T = J(rows(ua), rows(ub), 0)
	for (n = 1; n <= rows(a); n++) T[selectindex(ua :== a[n]), selectindex(ub :== b[n])] = T[selectindex(ua :== a[n]), selectindex(ub :== b[n])] + 1
	n = rows(a)
	sij = sum(T :* (T :- 1) :/ 2)
	sa = sum(rowsum(T) :* (rowsum(T) :- 1) :/ 2)
	sb = sum(colsum(T)' :* (colsum(T)' :- 1) :/ 2)
	e = sa * sb / (n * (n - 1) / 2)
	mx = (sa + sb) / 2
	if (mx == e) return(1)
	return((sij - e) / (mx - e))
}
void sbmt_pidiff(string scalar vm, string scalar vr, real scalar n)
{
	real matrix PM, PR, M
	real colvector a, b
	real scalar q, k, j
	real colvector map
	a = st_data(., vm)[1::n]
	b = st_data(., vr)[1::n]
	PM = st_matrix("e(pi)")
	PR = st_matrix("RPI")
	q = rows(PM)
	map = J(q, 1, 0)
	for (k = 1; k <= q; k++) {
		M = J(q, 1, 0)
		for (j = 1; j <= q; j++) M[j] = sum((b :== k) :* (a :== j))
		map[k] = (select(1::q, M :== max(M)))[1]
	}
	// PM is stored lower-triangular-symmetric in Stata only by display; e(pi) is the full matrix
	st_numscalar("__pidiff", max(abs(PM[map, map] - PR)))
}
end

capture program drop sbmt_mknetv
* sbmt_mknetv name, ref(fixture) n(#) dir(0|1): valued network (edge list with values; absent = 0, missing list = .)
program sbmt_mknetv
	syntax name(name=nm), REF(string) N(integer) DIR(integer)
	mata: sbmt_loadv(`n', `n', `dir', "$SBM_REPO/dev/ref/`ref'_edges.csv", "$SBM_REPO/dev/ref/`ref'_missing.csv")
	preserve
	clear
	set obs `n'
	forvalues j = 1/`n' {
		quietly gen double __a`j' = .
	}
	mata: st_store(., ., st_matrix("__sbmt_A"))
	if `dir' nwset __a*, name(`nm') directed valued replace
	else nwset __a*, name(`nm') undirected valued replace
	restore
	capture matrix drop __sbmt_A
end

capture program drop sbmt_mknet2
* sbmt_mknet2 name, ref(fixture) n1(#) n2(#): two-mode network from an n1 x n2 matrix (edge list with values; missing list)
program sbmt_mknet2
	syntax name(name=nm), REF(string) N1(integer) N2(integer) [valued]
	mata: sbmt_loadv(`n1', `n2', -1, "$SBM_REPO/dev/ref/`ref'_edges.csv", "$SBM_REPO/dev/ref/`ref'_missing.csv")
	preserve
	clear
	set obs `n1'
	forvalues j = 1/`n2' {
		quietly gen double __a`j' = .
	}
	mata: st_store(., ., st_matrix("__sbmt_A"))
	nwset __a*, name(`nm') bipartite `valued' replace
	restore
	capture matrix drop __sbmt_A
end

mata:
void sbmt_loadv(real scalar n1, real scalar n2, real scalar dir, string scalar fe, string scalar fm)
{
	real matrix A, E
	real scalar k
	A = J(n1, n2, 0)
	E = sbmt_readcsv3(fe)
	for (k = 1; k <= rows(E); k++) {
		A[E[k, 1], E[k, 2]] = E[k, 3]
		if (dir == 0) A[E[k, 2], E[k, 1]] = E[k, 3]
	}
	E = sbmt_readcsv(fm)
	for (k = 1; k <= rows(E); k++) {
		A[E[k, 1], E[k, 2]] = .
		if (dir == 0) A[E[k, 2], E[k, 1]] = .
	}
	if (dir >= 0) _diag(A, 0)
	st_matrix("__sbmt_A", A)
}
real matrix sbmt_readcsv3(string scalar f)
{
	real scalar fh, k
	string scalar ln
	real matrix E
	string matrix L
	fh = fopen(f, "r")
	L = J(0, 1, "")
	while ((ln = fget(fh)) != J(0, 0, "")) L = L \ ln
	fclose(fh)
	if (rows(L) <= 1) return(J(0, 3, .))
	E = J(rows(L) - 1, 3, .)
	for (k = 2; k <= rows(L); k++) E[k - 1, .] = strtoreal(tokens(subinstr(L[k], ",", " "))[1..3])
	return(E)
}
end

mata:
real colvector sbmt_readcol(string scalar f)
{
	real scalar fh, k
	string scalar ln
	string matrix L
	real colvector v
	fh = fopen(f, "r")
	L = J(0, 1, "")
	while ((ln = fget(fh)) != J(0, 0, "")) L = L \ ln
	fclose(fh)
	v = J(rows(L) - 1, 1, .)
	for (k = 2; k <= rows(L); k++) v[k - 1] = strtoreal(tokens(subinstr(L[k], ",", " "))[2])
	return(v)
}
void sbmt_ari_csv(string scalar var, string scalar f, real colvector idx)
{
	st_numscalar("__ari", sbmt_ari(st_data(idx, var), sbmt_readcol(f)))
}
// largest |difference| between e(pi) and R's pi (R rows x R columns), blocks matched by overlap on both modes
void sbmt_pidiff2(string scalar var, string scalar fcol, string scalar frow, real colvector i1, real colvector i2)
{
	real matrix PM, PR, M
	real colvector a1, b1, a2, b2, map1, map2
	real scalar q1, q2, k, j
	a1 = st_data(i1, var); b1 = sbmt_readcol(fcol)
	a2 = st_data(i2, var); b2 = sbmt_readcol(frow)
	PM = st_matrix("e(pi)")
	PR = st_matrix("RPI")'
	q1 = rows(PM); q2 = cols(PM)
	map1 = J(q1, 1, 0); map2 = J(q2, 1, 0)
	for (k = 1; k <= q1; k++) {
		M = J(q1, 1, 0)
		for (j = 1; j <= q1; j++) M[j] = sum((b1 :== k) :* (a1 :== j))
		map1[k] = (select(1::q1, M :== max(M)))[1]
	}
	for (k = 1; k <= q2; k++) {
		M = J(q2, 1, 0)
		for (j = 1; j <= q2; j++) M[j] = sum((b2 :== k) :* (a2 :== j))
		map2[k] = (select(1::q2, M :== max(M)))[1]
	}
	st_numscalar("__pidiff", max(abs(PM[map1, map2] - PR)))
}
end
