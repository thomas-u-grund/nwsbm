/*
 * Copyright (C) 2026 Thomas U. Grund
 * This file is part of nwsbm, a Stata package for stochastic block models. It follows the variational EM of the R
 * packages blockmodels (Leger) and sbm (Chiquet, Donnet, Barbillon) and is written independently from the published method.
 * nwsbm is free software: you can redistribute it and/or modify it under the terms of the GNU General
 * Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.
 */
/*
 * sbm_vem -- Bernoulli stochastic block model by variational EM, with spectral / split / random starts and ICL.
 *
 * Input block (column v1 of the frame, one number per row):
 *   [0] protocol  [1] n  [2] directed (0/1)  [3] model (0 = bernoulli)  [4] Qmin  [5] Qmax  [6] seed  [7] restarts
 *   [8..15] reserved  then the n x n adjacency, row major, values 0, 1 or -1 (missing); the diagonal is ignored.
 * Output (column v2, from row 1): for Q = Qmin..Qmax: ICL, J, PL, H, npar, tau (n x Q row major), pi (Q x Q row major), alpha (Q).
 * Call: plugin call sbmvem v1 v2, "10 <L>" fits; "99 0" returns the protocol number in v2 row 1.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stplugin.h"

#define SBM_PROTOCOL 1
#define HDR 16
#define TOL_P 1e-1
#define TOL_F 1e-1
#define TOL_EM 1e-5
#define MAXIT 1000
#define LOGFLOOR 1e-12

typedef struct { int *p, *j; } Csr;

typedef struct {
	int n, directed;
	Csr ao, ai, mo, mi;      /* adjacency out/in, missing out/in (for undirected, in == out) */
	int nmiss;               /* number of missing off-diagonal ordered pairs */
	double ndata;            /* observed pairs (ordered if directed, unordered if not) */
} Net;

/* ---------- random numbers (xorshift64*) ---------- */
static unsigned long long rng_s;
static void rng_seed(unsigned long long s) { rng_s = s * 2685821657736338717ULL + 88172645463325252ULL; if (!rng_s) rng_s = 1; }
static double rng_u(void) {
	rng_s ^= rng_s >> 12; rng_s ^= rng_s << 25; rng_s ^= rng_s >> 27;
	return ((rng_s * 2685821657736338717ULL) >> 11) * (1.0 / 9007199254740992.0);
}
static int rng_int(int k) { int r = (int)(rng_u() * k); return r >= k ? k - 1 : r; }

/* ---------- sparse structure ---------- */
static int csr_build(Csr *c, int n, const int *ri, const int *ci, long m) {
	long k; int i;
	c->p = calloc(n + 1, sizeof(int)); c->j = malloc(sizeof(int) * (m > 0 ? m : 1));
	if (!c->p || !c->j) return 0;
	for (k = 0; k < m; k++) c->p[ri[k] + 1]++;
	for (i = 0; i < n; i++) c->p[i + 1] += c->p[i];
	{
		int *fill = malloc(sizeof(int) * (n + 1));
		if (!fill) return 0;
		memcpy(fill, c->p, sizeof(int) * (n + 1));
		for (k = 0; k < m; k++) c->j[fill[ri[k]]++] = ci[k];
		free(fill);
	}
	return 1;
}
static void csr_free(Csr *c) { free(c->p); free(c->j); c->p = c->j = NULL; }

/* out = S X, X n x k row major */
static void spmm(const Csr *c, int n, const double *X, int k, double *out) {
	int i, l; long e;
	for (i = 0; i < n; i++) {
		double *o = out + (long)i * k;
		for (l = 0; l < k; l++) o[l] = 0;
		for (e = c->p[i]; e < c->p[i + 1]; e++) {
			const double *x = X + (long)c->j[e] * k;
			for (l = 0; l < k; l++) o[l] += x[l];
		}
	}
}

/* N_i = sum over observed j != i of X_j */
static void obs_sum(const Net *g, const Csr *m, const double *X, int k, double *out) {
	int i, l, n = g->n;
	double *cs = calloc(k, sizeof(double)), *tmp = malloc(sizeof(double) * (long)n * k);
	for (i = 0; i < n; i++) for (l = 0; l < k; l++) cs[l] += X[(long)i * k + l];
	if (g->nmiss) spmm(m, n, X, k, tmp); else memset(tmp, 0, sizeof(double) * (long)n * k);
	for (i = 0; i < n; i++) for (l = 0; l < k; l++) out[(long)i * k + l] = cs[l] - X[(long)i * k + l] - tmp[(long)i * k + l];
	free(cs); free(tmp);
}

/* ---------- one fitted model ---------- */
typedef struct {
	int Q;
	double *tau, *pi, *alpha;
	double PL, H, J, icl, npar;
} Fit;

static void fit_alloc(Fit *f, int n, int Q) {
	f->Q = Q;
	f->tau = malloc(sizeof(double) * (long)n * Q); f->pi = malloc(sizeof(double) * Q * Q); f->alpha = malloc(sizeof(double) * Q);
}
static void fit_free(Fit *f) { free(f->tau); free(f->pi); free(f->alpha); f->tau = f->pi = f->alpha = NULL; }

static void clamp_tau(double *tau, int n, int Q) {
	double tol = TOL_P / n; int i, q;
	for (i = 0; i < n; i++) {
		double s = 0;
		double *t = tau + (long)i * Q;
		for (q = 0; q < Q; q++) { if (t[q] < tol) t[q] = tol; else if (t[q] > 1 - tol) t[q] = 1 - tol; s += t[q]; }
		for (q = 0; q < Q; q++) t[q] /= s;
	}
}

/* workspace for the VEM */
typedef struct { double *AZ, *NZ, *AZi, *NZi, *lz, *tmp; } Work;
static int work_alloc(Work *w, int n, int Q) {
	long s = (long)n * Q;
	w->AZ = malloc(sizeof(double) * s); w->NZ = malloc(sizeof(double) * s); w->AZi = malloc(sizeof(double) * s);
	w->NZi = malloc(sizeof(double) * s); w->lz = malloc(sizeof(double) * s); w->tmp = malloc(sizeof(double) * s);
	return w->AZ && w->NZ && w->AZi && w->NZi && w->lz && w->tmp;
}
static void work_free(Work *w) { free(w->AZ); free(w->NZ); free(w->AZi); free(w->NZi); free(w->lz); free(w->tmp); }

/* M step: alpha, pi from tau; returns PL; sets H */
static double m_step(const Net *g, Fit *f, Work *w) {
	int n = g->n, Q = f->Q, i, q, l; double pl = 0;
	double *ZA = calloc(Q * Q, sizeof(double)), *ZO = calloc(Q * Q, sizeof(double));
	spmm(&g->ao, n, f->tau, Q, w->AZ);
	obs_sum(g, &g->mo, f->tau, Q, w->NZ);
	for (i = 0; i < n; i++) for (q = 0; q < Q; q++) {
		double t = f->tau[(long)i * Q + q];
		for (l = 0; l < Q; l++) { ZA[q * Q + l] += t * w->AZ[(long)i * Q + l]; ZO[q * Q + l] += t * w->NZ[(long)i * Q + l]; }
	}
	for (q = 0; q < Q; q++) { double s = 0; for (i = 0; i < n; i++) s += f->tau[(long)i * Q + q]; f->alpha[q] = s / n; }
	for (q = 0; q < Q * Q; q++) {
		double p = ZO[q] > 0 ? ZA[q] / ZO[q] : 0.0;
		if (p < 0) p = 0; if (p > 1) p = 1;
		f->pi[q] = p;
		if (ZA[q] > 0) pl += ZA[q] * log(p);
		if (ZO[q] - ZA[q] > 1e-12) pl += (ZO[q] - ZA[q]) * log(1 - p);
	}
	if (!g->directed) pl /= 2;
	for (i = 0; i < n; i++) for (q = 0; q < Q; q++) pl += f->tau[(long)i * Q + q] * log(f->alpha[q] > LOGFLOOR ? f->alpha[q] : LOGFLOOR);
	{ double h = 0; for (i = 0; i < n * Q; i++) if (f->tau[i] > 0) h -= f->tau[i] * log(f->tau[i]); f->H = h; }
	free(ZA); free(ZO);
	return pl;
}

/* E step (fixed point iterations) */
static void e_step(const Net *g, Fit *f, Work *w) {
	int n = g->n, Q = f->Q, i, q, l, it;
	double *lp = malloc(sizeof(double) * Q * Q), *lq = malloc(sizeof(double) * Q * Q);
	for (q = 0; q < Q * Q; q++) {
		double p = f->pi[q]; if (p < LOGFLOOR) p = LOGFLOOR; if (p > 1 - LOGFLOOR) p = 1 - LOGFLOOR;
		lp[q] = log(p); lq[q] = log(1 - p);
	}
	for (it = 0; it < 10; it++) {
		double step = 0;
		spmm(&g->ao, n, f->tau, Q, w->AZ);
		obs_sum(g, &g->mo, f->tau, Q, w->NZ);
		if (g->directed) { spmm(&g->ai, n, f->tau, Q, w->AZi); obs_sum(g, &g->mi, f->tau, Q, w->NZi); }
		for (i = 0; i < n; i++) {
			double mx = -1e300, s = 0, *z = w->lz + (long)i * Q;
			for (q = 0; q < Q; q++) {
				double v = log(f->alpha[q] > LOGFLOOR ? f->alpha[q] : LOGFLOOR);
				for (l = 0; l < Q; l++) {
					double a = w->AZ[(long)i * Q + l], nn = w->NZ[(long)i * Q + l];
					v += a * lp[q * Q + l] + (nn - a) * lq[q * Q + l];
					if (g->directed) {
						double a2 = w->AZi[(long)i * Q + l], n2 = w->NZi[(long)i * Q + l];
						v += a2 * lp[l * Q + q] + (n2 - a2) * lq[l * Q + q];
					}
				}
				z[q] = v; if (v > mx) mx = v;
			}
			for (q = 0; q < Q; q++) { z[q] = exp(z[q] - mx); s += z[q]; }
			for (q = 0; q < Q; q++) z[q] /= s;
		}
		clamp_tau(w->lz, n, Q);
		for (i = 0; i < n * Q; i++) { double d = fabs(w->lz[i] - f->tau[i]); if (d > step) step = d; }
		memcpy(f->tau, w->lz, sizeof(double) * (long)n * Q);
		if (step <= TOL_F) break;
	}
	free(lp); free(lq);
}

static double icl_of(const Net *g, const Fit *f) {
	int Q = f->Q;
	double npar = g->directed ? (double)Q * Q : (double)Q * (Q + 1) / 2;
	return f->PL - 0.5 * (npar * log(g->ndata) + (Q - 1) * log((double)g->n));
}

/* full VEM from the start tau in f */
static void vem(const Net *g, Fit *f, Work *w) {
	int Q = f->Q, it;
	double J, old;
	clamp_tau(f->tau, g->n, Q);
	f->PL = m_step(g, f, w);
	J = f->PL + f->H;
	if (Q > 1) {
		for (it = 0; it < MAXIT; it++) {
			e_step(g, f, w);
			f->PL = m_step(g, f, w);
			old = J; J = f->PL + f->H;
			if (!(J - old > TOL_EM)) break;
		}
	}
	f->J = J;
	f->npar = g->directed ? (double)Q * Q : (double)Q * (Q + 1) / 2;
	f->icl = icl_of(g, f);
}

/* ---------- spectral embedding ---------- */
typedef struct { int k; double *V; /* n x k row major, columns ordered by |eigenvalue| */ } Embed;

static void jacobi_eig(int k, double *A, double *E, double *lam) { /* A symmetric k x k (destroyed), eigenvectors in columns of E */
	int i, j, p, q, sweep;
	for (i = 0; i < k * k; i++) E[i] = 0;
	for (i = 0; i < k; i++) E[i * k + i] = 1;
	for (sweep = 0; sweep < 60; sweep++) {
		double off = 0;
		for (p = 0; p < k; p++) for (q = p + 1; q < k; q++) off += A[p * k + q] * A[p * k + q];
		if (off < 1e-26) break;
		for (p = 0; p < k; p++) for (q = p + 1; q < k; q++) {
			double apq = A[p * k + q], th, t, c, s;
			if (fabs(apq) < 1e-300) continue;
			th = (A[q * k + q] - A[p * k + p]) / (2 * apq);
			t = (th >= 0 ? 1.0 : -1.0) / (fabs(th) + sqrt(th * th + 1));
			c = 1 / sqrt(t * t + 1); s = t * c;
			for (j = 0; j < k; j++) { double x = A[j * k + p], y = A[j * k + q]; A[j * k + p] = c * x - s * y; A[j * k + q] = s * x + c * y; }
			for (j = 0; j < k; j++) { double x = A[p * k + j], y = A[q * k + j]; A[p * k + j] = c * x - s * y; A[q * k + j] = s * x + c * y; }
			for (j = 0; j < k; j++) { double x = E[j * k + p], y = E[j * k + q]; E[j * k + p] = c * x - s * y; E[j * k + q] = s * x + c * y; }
		}
	}
	for (i = 0; i < k; i++) lam[i] = A[i * k + i];
}

static void orthonormalise(double *V, int n, int k) { /* modified Gram-Schmidt on columns of V (n x k row major) */
	int a, b, i;
	for (a = 0; a < k; a++) {
		for (b = 0; b < a; b++) {
			double d = 0; for (i = 0; i < n; i++) d += V[(long)i * k + a] * V[(long)i * k + b];
			for (i = 0; i < n; i++) V[(long)i * k + a] -= d * V[(long)i * k + b];
		}
		{ double s = 0; for (i = 0; i < n; i++) s += V[(long)i * k + a] * V[(long)i * k + a];
		  s = sqrt(s); if (s < 1e-300) { for (i = 0; i < n; i++) V[(long)i * k + a] = rng_u() - 0.5; s = 1; for (i = 0, s = 0; i < n; i++) s += V[(long)i*k+a]*V[(long)i*k+a]; s = sqrt(s); }
		  for (i = 0; i < n; i++) V[(long)i * k + a] /= s; }
	}
}

/* B V for B = (A + A')/2 - p (J - I) restricted to observed pairs */
static void resid_apply(const Net *g, double p, const double *V, int k, double *out, double *t1, double *t2) {
	int n = g->n, i, l;
	spmm(&g->ao, n, V, k, out);
	if (g->directed) {
		spmm(&g->ai, n, V, k, t1);
		for (i = 0; i < n * k; i++) out[i] = 0.5 * (out[i] + t1[i]);
	}
	obs_sum(g, &g->mo, V, k, t2);
	if (g->directed) { obs_sum(g, &g->mi, V, k, t1); for (i = 0; i < n * k; i++) t2[i] = 0.5 * (t2[i] + t1[i]); }
	for (i = 0; i < n; i++) for (l = 0; l < k; l++) out[(long)i * k + l] -= p * t2[(long)i * k + l];
}

static int embed_compute(const Net *g, double p, int k, Embed *em) {
	int n = g->n, it, a, b, i;
	double *V = malloc(sizeof(double) * (long)n * k), *W = malloc(sizeof(double) * (long)n * k);
	double *t1 = malloc(sizeof(double) * (long)n * k), *t2 = malloc(sizeof(double) * (long)n * k);
	double *H = malloc(sizeof(double) * k * k), *E = malloc(sizeof(double) * k * k), *lam = malloc(sizeof(double) * k), *prev = calloc(k, sizeof(double));
	if (!V || !W || !t1 || !t2 || !H || !E || !lam || !prev) return 0;
	for (i = 0; i < n * k; i++) V[i] = rng_u() - 0.5;
	orthonormalise(V, n, k);
	for (it = 0; it < 500; it++) {
		double ch = 0;
		resid_apply(g, p, V, k, W, t1, t2);
		for (a = 0; a < k; a++) for (b = 0; b < k; b++) { double d = 0; for (i = 0; i < n; i++) d += V[(long)i * k + a] * W[(long)i * k + b]; H[a * k + b] = d; }
		for (a = 0; a < k; a++) for (b = a + 1; b < k; b++) { double m = 0.5 * (H[a * k + b] + H[b * k + a]); H[a * k + b] = H[b * k + a] = m; }
		jacobi_eig(k, H, E, lam);
		for (a = 0; a < k; a++) { double d = fabs(fabs(lam[a]) - prev[a]); if (d > ch * fabs(lam[a]) ) ch = d / (fabs(lam[a]) + 1e-12); prev[a] = fabs(lam[a]); }
		/* V <- W (not rotated) then orthonormalise; the Ritz rotation is applied at the end */
		memcpy(V, W, sizeof(double) * (long)n * k);
		orthonormalise(V, n, k);
		if (it > 3 && ch < 1e-7) break;
	}
	/* Rayleigh-Ritz on the final V */
	resid_apply(g, p, V, k, W, t1, t2);
	for (a = 0; a < k; a++) for (b = 0; b < k; b++) { double d = 0; for (i = 0; i < n; i++) d += V[(long)i * k + a] * W[(long)i * k + b]; H[a * k + b] = d; }
	for (a = 0; a < k; a++) for (b = a + 1; b < k; b++) { double m = 0.5 * (H[a * k + b] + H[b * k + a]); H[a * k + b] = H[b * k + a] = m; }
	jacobi_eig(k, H, E, lam);
	{
		int *ord = malloc(sizeof(int) * k); double *R = malloc(sizeof(double) * (long)n * k);
		for (a = 0; a < k; a++) ord[a] = a;
		for (a = 0; a < k; a++) for (b = a + 1; b < k; b++) if (fabs(lam[ord[b]]) > fabs(lam[ord[a]])) { int t = ord[a]; ord[a] = ord[b]; ord[b] = t; }
		for (i = 0; i < n; i++) for (a = 0; a < k; a++) {
			double s = 0; for (b = 0; b < k; b++) s += V[(long)i * k + b] * E[b * k + ord[a]];
			R[(long)i * k + a] = s;
		}
		em->k = k; em->V = R; free(ord);
	}
	free(V); free(W); free(t1); free(t2); free(H); free(E); free(lam); free(prev);
	return 1;
}

/* ---------- k-means on the first d columns of the embedding (row stride k) ---------- */
static double sqd(const double *x, const double *y, int d) { double s = 0; int a; for (a = 0; a < d; a++) { double t = x[a] - y[a]; s += t * t; } return s; }

/* mode 0: farthest point start (deterministic); mode 1: k-means++ (random) ; idx: optional subset of rows (m of them) */
static void kmeans(const double *V, int k, int d, const int *idx, int m, int Q, int mode, int *lab) {
	int i, c, it, changed;
	double *cen = malloc(sizeof(double) * Q * d), *mind = malloc(sizeof(double) * m);
	int *cnt = malloc(sizeof(int) * Q);
	#define ROW(r) (V + (long)(idx ? idx[r] : (r)) * k)
	if (Q == 1) { for (i = 0; i < m; i++) lab[i] = 0; free(cen); free(mind); free(cnt); return; }
	if (mode == 0) {
		int a = 0, b = 0; double best = -1, mu[64];
		int dd = d > 64 ? 64 : d;
		for (c = 0; c < dd; c++) { mu[c] = 0; for (i = 0; i < m; i++) mu[c] += ROW(i)[c]; mu[c] /= m; }
		for (i = 0; i < m; i++) { double s = sqd(ROW(i), mu, dd); if (s > best) { best = s; a = i; } }
		best = -1; for (i = 0; i < m; i++) { double s = sqd(ROW(i), ROW(a), d); if (s > best) { best = s; b = i; } }
		memcpy(cen, ROW(a), sizeof(double) * d); memcpy(cen + d, ROW(b), sizeof(double) * d);
		for (i = 0; i < m; i++) { double s1 = sqd(ROW(i), cen, d), s2 = sqd(ROW(i), cen + d, d); mind[i] = s1 < s2 ? s1 : s2; }
		for (c = 2; c < Q; c++) {
			int far = 0; best = -1;
			for (i = 0; i < m; i++) if (mind[i] > best) { best = mind[i]; far = i; }
			memcpy(cen + c * d, ROW(far), sizeof(double) * d);
			for (i = 0; i < m; i++) { double s = sqd(ROW(i), cen + c * d, d); if (s < mind[i]) mind[i] = s; }
		}
	} else {
		int first = rng_int(m);
		memcpy(cen, ROW(first), sizeof(double) * d);
		for (i = 0; i < m; i++) mind[i] = sqd(ROW(i), cen, d);
		for (c = 1; c < Q; c++) {
			double tot = 0, u; int pick = m - 1;
			for (i = 0; i < m; i++) tot += mind[i];
			u = rng_u() * tot;
			if (tot <= 0) pick = rng_int(m);
			else { double acc = 0; for (i = 0; i < m; i++) { acc += mind[i]; if (acc >= u) { pick = i; break; } } }
			memcpy(cen + c * d, ROW(pick), sizeof(double) * d);
			for (i = 0; i < m; i++) { double s = sqd(ROW(i), cen + c * d, d); if (s < mind[i]) mind[i] = s; }
		}
	}
	for (i = 0; i < m; i++) lab[i] = -1;
	for (it = 0; it < 100; it++) {
		changed = 0;
		for (i = 0; i < m; i++) {
			int bc = 0; double bs = 1e300;
			for (c = 0; c < Q; c++) { double s = sqd(ROW(i), cen + c * d, d); if (s < bs) { bs = s; bc = c; } }
			if (lab[i] != bc) { lab[i] = bc; changed = 1; }
		}
		if (!changed) break;
		memset(cen, 0, sizeof(double) * Q * d); memset(cnt, 0, sizeof(int) * Q);
		for (i = 0; i < m; i++) { cnt[lab[i]]++; for (c = 0; c < d; c++) cen[lab[i] * d + c] += ROW(i)[c]; }
		for (c = 0; c < Q; c++) {
			if (cnt[c] == 0) { int r = rng_int(m); memcpy(cen + c * d, ROW(r), sizeof(double) * d); }
			else for (i = 0; i < d; i++) cen[c * d + i] /= cnt[c];
		}
	}
	#undef ROW
	free(cen); free(mind); free(cnt);
}

static void tau_from_labels(double *tau, int n, int Q, const int *lab) {
	int i, q;
	for (i = 0; i < n; i++) for (q = 0; q < Q; q++) tau[(long)i * Q + q] = (lab[i] == q) ? 1.0 : 0.0;
}

/* ---------- driver ---------- */
static void keep_if_better(Fit *best, int *have, Fit *cand, int n) {
	if (!*have || cand->icl > best->icl) {
		if (*have) fit_free(best);
		*best = *cand; *have = 1;
		cand->tau = cand->pi = cand->alpha = NULL;
	} else fit_free(cand);
}

static void hard_labels(const Fit *f, int n, int *lab) {
	int i, k, Q = f->Q;
	for (i = 0; i < n; i++) {
		int bq = 0; double bv = -1;
		for (k = 0; k < Q; k++) if (f->tau[(long)i * Q + k] > bv) { bv = f->tau[(long)i * Q + k]; bq = k; }
		lab[i] = bq;
	}
}

/* forward step: starts for Q made by splitting each block of the Q-1 solution `src` in two (k-means on the embedding) */
static void try_splits(const Net *g, const Embed *em, int kk, const Fit *src, int Q, Fit *best, int *have, Work *w) {
	int n = g->n, i, q, d = Q < kk ? Q : kk;
	int *idx = malloc(sizeof(int) * n), *hard = malloc(sizeof(int) * n), *hs = malloc(sizeof(int) * n), *l2 = malloc(sizeof(int) * n);
	hard_labels(src, n, hard);
	for (q = 0; q < Q - 1; q++) {
		int m = 0;
		for (i = 0; i < n; i++) if (hard[i] == q) idx[m++] = i;
		if (m < 2) continue;
		memcpy(hs, hard, sizeof(int) * n);
		kmeans(em->V, kk, d, idx, m, 2, 0, l2);
		for (i = 0; i < m; i++) if (l2[i] == 1) hs[idx[i]] = Q - 1;
		{ Fit c; fit_alloc(&c, n, Q); tau_from_labels(c.tau, n, Q, hs); vem(g, &c, w); keep_if_better(best, have, &c, n); }
	}
	free(idx); free(hard); free(hs); free(l2);
}

/* backward step: starts for Q made by merging two blocks of the Q+1 solution `src` (at most 20 pairs: the most similar ones) */
static void try_merges(const Net *g, const Fit *src, int Q, Fit *best, int *have, Work *w) {
	int n = g->n, Q1 = Q + 1, a, b, i, np = 0, k, maxp = 20;
	int *lab = malloc(sizeof(int) * n), *hard = malloc(sizeof(int) * n);
	int *pa = malloc(sizeof(int) * Q1 * Q1), *pb = malloc(sizeof(int) * Q1 * Q1);
	double *pd = malloc(sizeof(double) * Q1 * Q1);
	hard_labels(src, n, hard);
	for (a = 0; a < Q1; a++) for (b = a + 1; b < Q1; b++) {
		double d = 0;
		for (k = 0; k < Q1; k++) { double x = src->pi[a * Q1 + k] - src->pi[b * Q1 + k], y = src->pi[k * Q1 + a] - src->pi[k * Q1 + b]; d += x * x + y * y; }
		pa[np] = a; pb[np] = b; pd[np] = d; np++;
	}
	if (np > maxp) { /* keep the maxp smallest distances (partial selection sort) */
		int t, u;
		for (t = 0; t < maxp; t++) {
			int mi = t;
			for (u = t + 1; u < np; u++) if (pd[u] < pd[mi]) mi = u;
			{ double dd = pd[t]; int ta = pa[t], tb = pb[t]; pd[t] = pd[mi]; pa[t] = pa[mi]; pb[t] = pb[mi]; pd[mi] = dd; pa[mi] = ta; pb[mi] = tb; }
		}
		np = maxp;
	}
	for (k = 0; k < np; k++) {
		a = pa[k]; b = pb[k];
		for (i = 0; i < n; i++) {
			int l = hard[i];
			if (l == b) l = a;
			if (l > b) l--;
			lab[i] = l;
		}
		{ Fit c; fit_alloc(&c, n, Q); tau_from_labels(c.tau, n, Q, lab); vem(g, &c, w); keep_if_better(best, have, &c, n); }
	}
	free(lab); free(hard); free(pa); free(pb); free(pd);
}

static int run(Net *g, int Qmin, int Qout, int restarts, double *out) {
	/* the search runs one level above the largest Q asked for (Qout) so that merges of the (Qout+1) solution also improve Qout */
	int Qmax = Qout < g->n ? Qout + 1 : Qout;
	int n = g->n, Q, r, i;
	Fit *best; int *have; Work w; Embed em; double p1;
	long off = 0;
	int *lab = malloc(sizeof(int) * n);
	int kk = (Qmax + 3 < n) ? Qmax + 3 : n;
	best = calloc(Qmax + 2, sizeof(Fit)); have = calloc(Qmax + 2, sizeof(int));
	if (!lab || !best || !have || !work_alloc(&w, n, Qmax > 1 ? Qmax : 1)) return 909;
	{ /* one-block density for the residual matrix */
		Fit f1; Work w1; fit_alloc(&f1, n, 1); work_alloc(&w1, n, 1);
		for (i = 0; i < n; i++) f1.tau[i] = 1.0;
		(void)m_step(g, &f1, &w1);
		p1 = f1.pi[0];
		fit_free(&f1); work_free(&w1);
	}
	if (!embed_compute(g, p1, kk, &em)) return 909;
	/* forward pass: spectral start, splits of the best (Q-1) solution, extra random starts */
	for (Q = 1; Q <= Qmax; Q++) {
		if (Q == 1) {
			Fit c; fit_alloc(&c, n, 1);
			for (i = 0; i < n; i++) c.tau[i] = 1.0;
			vem(g, &c, &w); keep_if_better(&best[1], &have[1], &c, n);
		} else {
			int d = Q < kk ? Q : kk;
			{ Fit c; fit_alloc(&c, n, Q); kmeans(em.V, kk, d, NULL, n, Q, 0, lab); tau_from_labels(c.tau, n, Q, lab); vem(g, &c, &w); keep_if_better(&best[Q], &have[Q], &c, n); }
			try_splits(g, &em, kk, &best[Q - 1], Q, &best[Q], &have[Q], &w);
			for (r = 0; r < restarts; r++) {
				Fit c; fit_alloc(&c, n, Q);
				if (r % 2 == 0) kmeans(em.V, kk, d, NULL, n, Q, 1, lab);
				else for (i = 0; i < n; i++) lab[i] = rng_int(Q);
				tau_from_labels(c.tau, n, Q, lab);
				vem(g, &c, &w); keep_if_better(&best[Q], &have[Q], &c, n);
			}
		}
	}
	/* backward pass: merges of the best (Q+1) solution; then one more forward sweep of splits from the improved solutions */
	for (Q = Qmax - 1; Q >= 2; Q--) try_merges(g, &best[Q + 1], Q, &best[Q], &have[Q], &w);
	for (Q = 2; Q <= Qmax; Q++) try_splits(g, &em, kk, &best[Q - 1], Q, &best[Q], &have[Q], &w);
	for (Q = Qmin; Q <= Qout; Q++) {
		Fit *b = &best[Q];
		out[off++] = b->icl; out[off++] = b->J; out[off++] = b->PL; out[off++] = b->H; out[off++] = b->npar;
		memcpy(out + off, b->tau, sizeof(double) * (long)n * Q); off += (long)n * Q;
		memcpy(out + off, b->pi, sizeof(double) * Q * Q); off += Q * Q;
		memcpy(out + off, b->alpha, sizeof(double) * Q); off += Q;
	}
	for (Q = 1; Q <= Qmax; Q++) if (have[Q]) fit_free(&best[Q]);
	free(best); free(have); free(em.V); free(lab); work_free(&w);
	return 0;
}

STDLL stata_call(int argc, char *argv[]) {
	char *buf, *tok;
	int code, rc = 0;
	long L = 0;
	if (argc < 1) { SF_error("sbm_vem: missing argument string\n"); return 198; }
	buf = malloc(strlen(argv[0]) + 1);
	if (!buf) return 909;
	strcpy(buf, argv[0]);
	tok = strtok(buf, " \t");
	code = tok ? (int)atof(tok) : 0;
	if (code == 99) { free(buf); SF_vstore(2, 1, (ST_double)SBM_PROTOCOL); return 0; }
	tok = strtok(NULL, " \t");
	L = tok ? (long)atof(tok) : 0;
	free(buf);
	if (code != 10) { SF_error("sbm_vem: unknown code (plugin and package versions do not match)\n"); return 198; }
	{
		ST_double v; int n, directed, model, Qmin, Qmax, restarts; unsigned long long seed;
		long i, j, ma = 0, mm = 0, cap, k, outn = 0, tot;
		int *ar, *ac, *mr, *mc; double *out; Net g;
		SF_vdata(1, 1, &v); if (v != SBM_PROTOCOL) { SF_error("sbm_vem: the model block is invalid (plugin and package versions do not match?)\n"); return 198; }
		SF_vdata(1, 2, &v); n = (int)v;
		SF_vdata(1, 3, &v); directed = (int)v;
		SF_vdata(1, 4, &v); model = (int)v;
		SF_vdata(1, 5, &v); Qmin = (int)v;
		SF_vdata(1, 6, &v); Qmax = (int)v;
		SF_vdata(1, 7, &v); seed = (unsigned long long)v;
		SF_vdata(1, 8, &v); restarts = (int)v;
		if (model != 0 || n < 2 || Qmin < 1 || Qmax < Qmin || Qmax > n || L != HDR + (long)n * n) { SF_error("sbm_vem: the model block is invalid\n"); return 198; }
		rng_seed(seed);
		cap = 1024;
		ar = malloc(sizeof(int) * cap); ac = malloc(sizeof(int) * cap); mr = malloc(sizeof(int) * cap); mc = malloc(sizeof(int) * cap);
		memset(&g, 0, sizeof(g)); g.n = n; g.directed = directed;
		{
			long obs = 0;
			for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
				if (i == j) continue;
				if (!directed && j < i) continue;
				SF_vdata(1, (ST_int)(HDR + i * n + j + 1), &v);
				if (!directed) { /* use the upper triangle value; the Mata layer symmetrises */ }
				if (v < 0) {
					if (ma + 2 >= cap || mm + 2 >= cap) { cap *= 2; ar = realloc(ar, sizeof(int) * cap); ac = realloc(ac, sizeof(int) * cap); mr = realloc(mr, sizeof(int) * cap); mc = realloc(mc, sizeof(int) * cap); }
					mr[mm] = (int)i; mc[mm] = (int)j; mm++;
					if (!directed) { mr[mm] = (int)j; mc[mm] = (int)i; mm++; }
				} else {
					obs++;
					if (v > 0) {
						if (ma + 2 >= cap || mm + 2 >= cap) { cap *= 2; ar = realloc(ar, sizeof(int) * cap); ac = realloc(ac, sizeof(int) * cap); mr = realloc(mr, sizeof(int) * cap); mc = realloc(mc, sizeof(int) * cap); }
						ar[ma] = (int)i; ac[ma] = (int)j; ma++;
						if (!directed) { ar[ma] = (int)j; ac[ma] = (int)i; ma++; }
					}
				}
			}
			g.ndata = (double)obs; g.nmiss = (int)mm;
		}
		if (!ar || !ac || !mr || !mc) return 909;
		if (g.ndata < 1) { SF_error("sbm_vem: no observed ties\n"); return 198; }
		csr_build(&g.ao, n, ar, ac, ma); csr_build(&g.mo, n, mr, mc, mm);
		if (directed) { csr_build(&g.ai, n, ac, ar, ma); csr_build(&g.mi, n, mc, mr, mm); }
		else { g.ai = g.ao; g.mi = g.mo; }
		tot = 0; for (k = Qmin; k <= Qmax; k++) tot += 5 + (long)n * k + k * k + k;
		out = malloc(sizeof(double) * tot);
		if (!out) return 909;
		rc = run(&g, Qmin, Qmax, restarts, out);
		if (!rc) for (outn = 0; outn < tot; outn++) SF_vstore(2, (ST_int)(1 + outn), out[outn]);
		free(out); free(ar); free(ac); free(mr); free(mc);
		csr_free(&g.ao); csr_free(&g.mo);
		if (directed) { csr_free(&g.ai); csr_free(&g.mi); }
	}
	return rc;
}
