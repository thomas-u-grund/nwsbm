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
 * sbm_vem -- stochastic block models by variational EM, with spectral / split / merge / random starts and ICL.
 * Models: Bernoulli, Poisson, Gaussian (free shared variance or fixed variance) edges; one-mode networks (undirected / directed
 * stochastic block model, optionally degree corrected for Poisson edges) and two-mode networks (latent block model).
 *
 * Input block (column v1 of the frame, one number per row, rows 1-based):
 *   row 1 protocol, 2 n (rows of the matrix), 3 directed (0/1), 4 model (0 bernoulli, 1 poisson, 2 gaussian), 5 Qmin, 6 Qmax, 7 seed,
 *   8 restarts, 9 fixed variance (<= 0: estimated), 10 threads (0 = one per core), 11 kind (0 one-mode, 1 two-mode), 12 n2 (columns,
 *   two-mode), 13 Q2min, 14 Q2max, 15 degree correction (0/1), 16 fast (0/1), 17 number of dyadic covariates p; then the n x n2 matrix, row major (missing ties are
 *   Stata missing values; the diagonal of a one-mode matrix is ignored), then the p covariate matrices (n x n2, no missing values).
 * Output (column v2, from row 1): one record per fit (one-mode: Q = Qmin..Qmax; two-mode: Q1 = Qmin..Qmax outer, Q2 = Q2min..Q2max inner):
 *   ICL, J, PL, H, npar, sigma2, Q1, Q2 (0 for one-mode), exact ICL of the hard partition (-1e300 if not available), tau (n x Q1, row major), [tau2 (n2 x Q2)], pi (Q1 x Q2 or Q x Q; with covariates on the response scale for covariate value 0), alpha (Q1), [alpha2 (Q2)], beta (p).
 * Call: plugin call sbmvem v1 v2, "10 <L>" fits; "99 0" returns the protocol number in v2 row 1.
 * The starts of one Q are evaluated by independent threads; every start is deterministic given the seed, so results do not depend
 * on the number of threads.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "stplugin.h"
#if SYSTEM == STWIN32
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#include <unistd.h>
#endif

#define SBM_PROTOCOL 5
#define HDR 17
#define MAXCOV 32
#define TOL_P 1e-1
#define TOL_F 1e-1
#define TOL_EM 1e-5
#define MAXIT 1000
#define SCREEN_TOL 5e-3   /* covariate Bernoulli model: starts are fitted to this tolerance (one Newton step per M step); the selected fits are polished to TOL_EM */
#define LOGFLOOR 1e-12

/* ---------- threads ---------- */
#if SYSTEM == STWIN32
typedef HANDLE sbm_thread;
#define SBM_THREAD_FN unsigned __stdcall
#define SBM_THREAD_RET return 0
static int thread_start(sbm_thread *t, unsigned (__stdcall *fn)(void *), void *arg) { *t = (HANDLE)_beginthreadex(NULL, 0, fn, arg, 0, NULL); return *t == 0 ? -1 : 0; }
static void thread_join(sbm_thread t) { WaitForSingleObject(t, INFINITE); CloseHandle(t); }
static int ncores(void) { SYSTEM_INFO si; GetSystemInfo(&si); return si.dwNumberOfProcessors > 0 ? (int)si.dwNumberOfProcessors : 1; }
#else
typedef pthread_t sbm_thread;
#define SBM_THREAD_FN void *
#define SBM_THREAD_RET return NULL
static int thread_start(sbm_thread *t, void *(*fn)(void *), void *arg) { return pthread_create(t, NULL, fn, arg) == 0 ? 0 : -1; }
static void thread_join(sbm_thread t) { pthread_join(t, NULL); }
static int ncores(void) { long c = sysconf(_SC_NPROCESSORS_ONLN); return c > 0 ? (int)c : 1; }
#endif

typedef struct { int *p, *j; double *v; } Csr;

typedef struct {
	int kind;                /* 0 one-mode (SBM), 1 two-mode (LBM) */
	int n, n2, directed, model, dc;
	Csr ao, ai, mo, mi;      /* one-mode: adjacency out/in, missing out/in (undirected: in == out); two-mode: rows -> columns and columns -> rows */
	int nmiss;
	double ndata;            /* observed pairs (one-mode: ordered if directed, unordered if not; two-mode: observed cells) */
	double nobs;             /* observed ordered cells */
	double sumy2;            /* sum of y^2 over the observed ordered cells (Gaussian) */
	double cst;              /* constant of the complete-data log likelihood (-sum log y! ; degree terms), in the units of PL */
	double sigma2fix;        /* fixed variance (Gaussian), <= 0: estimated */
	double ybar;             /* overall mean of the observed cells */
	double sumw;             /* sum of the out-weights (degree correction) */
	double *wo, *wi;         /* degree-correction weights: out and in degree (undirected: the degree) */
	int p;                   /* number of dyadic covariates (0: block statistics suffice, sparse path) */
	double *Y, *Mf;          /* covariate models: dense tie values (0 where unobserved) and observation indicator (diagonal of a one-mode matrix: 0) */
	double **X;              /* p dense covariate matrices (0 at unobserved cells) */
	double *sxy, *sxx;       /* sum x_k y, sum x_k x_k' over observed ordered cells */
	double *Rres;            /* residual matrix of the one-block covariate model (starts) */
} Net;

/* ---------- random numbers (xorshift64*), one stream per owner so that threads never share one ---------- */
typedef struct { unsigned long long s; } Rng;
static void rng_seed(Rng *r, unsigned long long s) { r->s = s * 2685821657736338717ULL + 88172645463325252ULL; if (!r->s) r->s = 1; }
static double rng_u(Rng *r) {
	r->s ^= r->s >> 12; r->s ^= r->s << 25; r->s ^= r->s >> 27;
	return ((r->s * 2685821657736338717ULL) >> 11) * (1.0 / 9007199254740992.0);
}
static int rng_int(Rng *r, int k) { int x = (int)(rng_u(r) * k); return x >= k ? k - 1 : x; }

/* ---------- sparse structure ---------- */
static int csr_build(Csr *c, int nr, const int *ri, const int *ci, const double *val, long m) {
	long k; int i;
	c->p = calloc(nr + 1, sizeof(int)); c->j = malloc(sizeof(int) * (m > 0 ? m : 1)); c->v = malloc(sizeof(double) * (m > 0 ? m : 1));
	if (!c->p || !c->j || !c->v) return 0;
	for (k = 0; k < m; k++) c->p[ri[k] + 1]++;
	for (i = 0; i < nr; i++) c->p[i + 1] += c->p[i];
	{
		int *fill = malloc(sizeof(int) * (nr + 1));
		if (!fill) return 0;
		memcpy(fill, c->p, sizeof(int) * (nr + 1));
		for (k = 0; k < m; k++) { int pos = fill[ri[k]]++; c->j[pos] = ci[k]; c->v[pos] = val ? val[k] : 1.0; }
		free(fill);
	}
	return 1;
}
static void csr_free(Csr *c) { free(c->p); free(c->j); free(c->v); c->p = c->j = NULL; c->v = NULL; }

/* out = S X (nr rows), optional values and optional column weights w; X row major with k columns */
static void spmm(const Csr *c, int nr, const double *X, int k, double *out, int useval, const double *w) {
	int i, l; long e;
	for (i = 0; i < nr; i++) {
		double *o = out + (long)i * k;
		for (l = 0; l < k; l++) o[l] = 0;
		for (e = c->p[i]; e < c->p[i + 1]; e++) {
			const double *x = X + (long)c->j[e] * k;
			double f = (useval ? c->v[e] : 1.0) * (w ? w[c->j[e]] : 1.0);
			for (l = 0; l < k; l++) o[l] += f * x[l];
		}
	}
}

/* one-mode: out_i = sum over observed j != i of w_j X_j */
static void obs_sum(const Net *g, const Csr *m, const double *X, int k, double *out, const double *w) {
	int i, l, n = g->n;
	double *cs = calloc(k, sizeof(double)), *tmp = malloc(sizeof(double) * (long)n * k);
	for (i = 0; i < n; i++) { double f = w ? w[i] : 1.0; for (l = 0; l < k; l++) cs[l] += f * X[(long)i * k + l]; }
	if (g->nmiss) spmm(m, n, X, k, tmp, 0, w); else memset(tmp, 0, sizeof(double) * (long)n * k);
	for (i = 0; i < n; i++) { double f = w ? w[i] : 1.0; for (l = 0; l < k; l++) out[(long)i * k + l] = cs[l] - f * X[(long)i * k + l] - tmp[(long)i * k + l]; }
	free(cs); free(tmp);
}

/* two-mode: out_i = sum over observed cells (i, j) of w_j X_j (w optional: degree correction); X has nx rows */
static void obs_sum_rect(const Net *g, const Csr *m, int nr, const double *X, int nx, int k, double *out, const double *w) {
	int i, l;
	double *cs = calloc(k, sizeof(double));
	for (i = 0; i < nx; i++) { double f = w ? w[i] : 1.0; for (l = 0; l < k; l++) cs[l] += f * X[(long)i * k + l]; }
	if (g->nmiss) spmm(m, nr, X, k, out, 0, w); else memset(out, 0, sizeof(double) * (long)nr * k);
	for (i = 0; i < nr; i++) for (l = 0; l < k; l++) out[(long)i * k + l] = cs[l] - out[(long)i * k + l];
	free(cs);
}

/* ---------- one fitted model ---------- */
typedef struct {
	int Q, Q2;                     /* Q2 = 0 for one-mode networks */
	double *tau, *tau2, *pi, *alpha, *alpha2;
	double *mq, *beta, *info;      /* covariate models: info = information matrix of beta (p x p, block parameters profiled out, conditional on the memberships); : block parameters on the link scale, regression coefficients */
	int warm, nnewton;             /* nnewton: Newton steps per M step once warm (covariate Bernoulli model) */
	double sigma2, PL, H, J, icl, npar;
} Fit;

static void fit_alloc(Fit *f, const Net *g, int Q, int Q2) {
	int c2 = Q2 ? Q2 : Q;
	f->Q = Q; f->Q2 = Q2;
	f->tau = malloc(sizeof(double) * (long)g->n * Q);
	f->tau2 = Q2 ? malloc(sizeof(double) * (long)g->n2 * Q2) : NULL;
	f->pi = malloc(sizeof(double) * Q * c2); f->alpha = malloc(sizeof(double) * Q);
	f->alpha2 = Q2 ? malloc(sizeof(double) * Q2) : NULL;
	f->sigma2 = 0; f->warm = 0; f->nnewton = 40;
	f->mq = g->p ? calloc((size_t)Q * c2, sizeof(double)) : NULL;
	f->beta = g->p ? calloc(g->p, sizeof(double)) : NULL;
	f->info = g->p ? calloc((size_t)g->p * g->p, sizeof(double)) : NULL;
}
static void fit_free(Fit *f) { free(f->tau); free(f->tau2); free(f->pi); free(f->alpha); free(f->alpha2); free(f->mq); free(f->beta); free(f->info); f->tau = f->tau2 = f->pi = f->alpha = f->alpha2 = f->mq = f->beta = f->info = NULL; }

static void clamp_tau(double *tau, int n, int Q) {
	double tol = TOL_P / n; int i, q;
	for (i = 0; i < n; i++) {
		double s = 0;
		double *t = tau + (long)i * Q;
		for (q = 0; q < Q; q++) { if (t[q] < tol) t[q] = tol; else if (t[q] > 1 - tol) t[q] = 1 - tol; s += t[q]; }
		for (q = 0; q < Q; q++) t[q] /= s;
	}
}

typedef struct { double *a1, *n1, *a2, *n2, *l1, *l2; double *Cc, *Ec, *Rm, *tmp; } Work;
static int work_alloc(Work *w, long cap, const Net *g) {
	long N = (long)g->n * g->n2;
	w->a1 = malloc(sizeof(double) * cap); w->n1 = malloc(sizeof(double) * cap); w->a2 = malloc(sizeof(double) * cap);
	w->n2 = malloc(sizeof(double) * cap); w->l1 = malloc(sizeof(double) * cap); w->l2 = malloc(sizeof(double) * cap);
	w->tmp = malloc(sizeof(double) * cap);
	w->Cc = w->Ec = w->Rm = NULL;
	if (g->p) { w->Cc = malloc(sizeof(double) * N); w->Ec = malloc(sizeof(double) * N); w->Rm = malloc(sizeof(double) * N); }
	return w->a1 && w->n1 && w->a2 && w->n2 && w->l1 && w->l2 && w->tmp && (!g->p || (w->Cc && w->Ec && w->Rm));
}
static void work_free(Work *w) { free(w->a1); free(w->n1); free(w->a2); free(w->n2); free(w->l1); free(w->l2); free(w->tmp); free(w->Cc); free(w->Ec); free(w->Rm); }

/* model part of PL at the maximiser for the block statistics ZA (sum of y) and ZO (number of observed cells, degree weighted for the degree-corrected model);
   fills the parameter matrix pi (probability / rate / mean) and the variance */
static double mpart(const Net *g, int cells, const double *ZA, const double *ZO, double *pi, double *sigma2) {
	int k; double pl = 0;
	if (g->model == 0) {
		for (k = 0; k < cells; k++) {
			double p = ZO[k] > 0 ? ZA[k] / ZO[k] : 0.0;
			if (p < 0) p = 0; if (p > 1) p = 1;
			pi[k] = p;
			if (ZA[k] > 0) pl += ZA[k] * log(p);
			if (ZO[k] - ZA[k] > 1e-12) pl += (ZO[k] - ZA[k]) * log(1 - p);
		}
	} else if (g->model == 1) {
		for (k = 0; k < cells; k++) {
			double lam = ZO[k] > 0 ? ZA[k] / ZO[k] : 0.0;
			pi[k] = lam;
			if (ZA[k] > 0) pl += ZA[k] * log(lam);
			pl -= lam * ZO[k];
		}
	} else {
		double resid = g->sumy2, s2;
		for (k = 0; k < cells; k++) {
			double mu = ZO[k] > 0 ? ZA[k] / ZO[k] : 0.0;
			pi[k] = mu;
			resid += mu * mu * ZO[k] - 2 * mu * ZA[k];
		}
		if (resid < 0) resid = 0;
		s2 = g->sigma2fix > 0 ? g->sigma2fix : resid / g->nobs;
		if (s2 < 1e-12) s2 = 1e-12;
		*sigma2 = s2;
		pl = -0.5 * g->nobs * log(2 * M_PI * s2) - resid / (2 * s2);
	}
	return pl;
}

/* E-step coefficients: log f(a, nn; q, l) = a U + nn V for a = sum of y, nn = number of observed cells (weighted) */
static void uv_fill(const Net *g, const Fit *f, int cells, double *U, double *V) {
	int k;
	for (k = 0; k < cells; k++) {
		double p = f->pi[k];
		if (g->model == 0) {
			if (p < LOGFLOOR) p = LOGFLOOR; if (p > 1 - LOGFLOOR) p = 1 - LOGFLOOR;
			U[k] = log(p) - log(1 - p); V[k] = log(1 - p);
		} else if (g->model == 1) {
			if (p < LOGFLOOR) p = LOGFLOOR;
			U[k] = log(p); V[k] = -p;
		} else {
			U[k] = p / f->sigma2; V[k] = -p * p / (2 * f->sigma2);
		}
	}
}

static double lg(double a) { return log(a > LOGFLOOR ? a : LOGFLOOR); }

/* ---- one-mode ---- */
static double m_step_sbm(const Net *g, Fit *f, Work *w) {
	int n = g->n, Q = f->Q, i, q, l; double pl = 0;
	double *ZA = calloc(Q * Q, sizeof(double)), *ZO = calloc(Q * Q, sizeof(double));
	spmm(&g->ao, n, f->tau, Q, w->a1, 1, NULL);
	obs_sum(g, &g->mo, f->tau, Q, w->n1, g->dc ? g->wi : NULL);
	for (i = 0; i < n; i++) {
		double own = g->dc ? g->wo[i] : 1.0;
		for (q = 0; q < Q; q++) {
			double t = f->tau[(long)i * Q + q];
			for (l = 0; l < Q; l++) { ZA[q * Q + l] += t * w->a1[(long)i * Q + l]; ZO[q * Q + l] += own * t * w->n1[(long)i * Q + l]; }
		}
	}
	for (q = 0; q < Q; q++) { double s = 0; for (i = 0; i < n; i++) s += f->tau[(long)i * Q + q]; f->alpha[q] = s / n; }
	pl = mpart(g, Q * Q, ZA, ZO, f->pi, &f->sigma2);
	if (!g->directed) pl /= 2;
	pl += g->cst;
	for (i = 0; i < n; i++) for (q = 0; q < Q; q++) pl += f->tau[(long)i * Q + q] * lg(f->alpha[q]);
	{ double h = 0; for (i = 0; i < n * Q; i++) if (f->tau[i] > 0) h -= f->tau[i] * log(f->tau[i]); f->H = h; }
	free(ZA); free(ZO);
	return pl;
}

static void e_step_sbm(const Net *g, Fit *f, Work *w) {
	int n = g->n, Q = f->Q, i, q, l, it;
	double *U = malloc(sizeof(double) * Q * Q), *V = malloc(sizeof(double) * Q * Q);
	uv_fill(g, f, Q * Q, U, V);
	for (it = 0; it < 10; it++) {
		double step = 0;
		spmm(&g->ao, n, f->tau, Q, w->a1, 1, NULL);
		obs_sum(g, &g->mo, f->tau, Q, w->n1, g->dc ? g->wi : NULL);
		if (g->directed) { spmm(&g->ai, n, f->tau, Q, w->a2, 1, NULL); obs_sum(g, &g->mi, f->tau, Q, w->n2, g->dc ? g->wo : NULL); }
		for (i = 0; i < n; i++) {
			double mx = -1e300, s = 0, *z = w->l1 + (long)i * Q;
			double own_o = g->dc ? g->wo[i] : 1.0, own_i = g->dc ? g->wi[i] : 1.0;
			for (q = 0; q < Q; q++) {
				double v = lg(f->alpha[q]);
				for (l = 0; l < Q; l++) {
					v += w->a1[(long)i * Q + l] * U[q * Q + l] + own_o * w->n1[(long)i * Q + l] * V[q * Q + l];
					if (g->directed) v += w->a2[(long)i * Q + l] * U[l * Q + q] + own_i * w->n2[(long)i * Q + l] * V[l * Q + q];
				}
				z[q] = v; if (v > mx) mx = v;
			}
			for (q = 0; q < Q; q++) { z[q] = exp(z[q] - mx); s += z[q]; }
			for (q = 0; q < Q; q++) z[q] /= s;
		}
		clamp_tau(w->l1, n, Q);
		for (i = 0; i < n * Q; i++) { double d = fabs(w->l1[i] - f->tau[i]); if (d > step) step = d; }
		memcpy(f->tau, w->l1, sizeof(double) * (long)n * Q);
		if (step <= TOL_F) break;
	}
	free(U); free(V);
}

/* ---- two-mode ---- */
static double m_step_lbm(const Net *g, Fit *f, Work *w) {
	int n1 = g->n, n2 = g->n2, Q1 = f->Q, Q2 = f->Q2, i, q, l; double pl;
	double *ZA = calloc(Q1 * Q2, sizeof(double)), *ZO = calloc(Q1 * Q2, sizeof(double));
	spmm(&g->ao, n1, f->tau2, Q2, w->a1, 1, NULL);
	obs_sum_rect(g, &g->mo, n1, f->tau2, n2, Q2, w->n1, g->dc ? g->wi : NULL);
	for (i = 0; i < n1; i++) for (q = 0; q < Q1; q++) {
		double t = f->tau[(long)i * Q1 + q], own = g->dc ? g->wo[i] : 1.0;
		for (l = 0; l < Q2; l++) { ZA[q * Q2 + l] += t * w->a1[(long)i * Q2 + l]; ZO[q * Q2 + l] += own * t * w->n1[(long)i * Q2 + l]; }
	}
	for (q = 0; q < Q1; q++) { double s = 0; for (i = 0; i < n1; i++) s += f->tau[(long)i * Q1 + q]; f->alpha[q] = s / n1; }
	for (l = 0; l < Q2; l++) { double s = 0; for (i = 0; i < n2; i++) s += f->tau2[(long)i * Q2 + l]; f->alpha2[l] = s / n2; }
	pl = mpart(g, Q1 * Q2, ZA, ZO, f->pi, &f->sigma2) + g->cst;
	for (i = 0; i < n1; i++) for (q = 0; q < Q1; q++) pl += f->tau[(long)i * Q1 + q] * lg(f->alpha[q]);
	for (i = 0; i < n2; i++) for (l = 0; l < Q2; l++) pl += f->tau2[(long)i * Q2 + l] * lg(f->alpha2[l]);
	{ double h = 0; for (i = 0; i < n1 * Q1; i++) if (f->tau[i] > 0) h -= f->tau[i] * log(f->tau[i]);
	  for (i = 0; i < n2 * Q2; i++) if (f->tau2[i] > 0) h -= f->tau2[i] * log(f->tau2[i]); f->H = h; }
	free(ZA); free(ZO);
	return pl;
}

static void e_step_lbm(const Net *g, Fit *f, Work *w) {
	int n1 = g->n, n2 = g->n2, Q1 = f->Q, Q2 = f->Q2, i, q, l, it;
	double *U = malloc(sizeof(double) * Q1 * Q2), *V = malloc(sizeof(double) * Q1 * Q2);
	uv_fill(g, f, Q1 * Q2, U, V);
	for (it = 0; it < 10; it++) {
		double step = 0, s, mx;
		spmm(&g->ao, n1, f->tau2, Q2, w->a1, 1, NULL);
		obs_sum_rect(g, &g->mo, n1, f->tau2, n2, Q2, w->n1, g->dc ? g->wi : NULL);
		spmm(&g->ai, n2, f->tau, Q1, w->a2, 1, NULL);
		obs_sum_rect(g, &g->mi, n2, f->tau, n1, Q1, w->n2, g->dc ? g->wo : NULL);
		for (i = 0; i < n1; i++) {
			double *z = w->l1 + (long)i * Q1, own_r = g->dc ? g->wo[i] : 1.0; mx = -1e300; s = 0;
			for (q = 0; q < Q1; q++) {
				double v = lg(f->alpha[q]);
				for (l = 0; l < Q2; l++) v += w->a1[(long)i * Q2 + l] * U[q * Q2 + l] + own_r * w->n1[(long)i * Q2 + l] * V[q * Q2 + l];
				z[q] = v; if (v > mx) mx = v;
			}
			for (q = 0; q < Q1; q++) { z[q] = exp(z[q] - mx); s += z[q]; }
			for (q = 0; q < Q1; q++) z[q] /= s;
		}
		for (i = 0; i < n2; i++) {
			double *z = w->l2 + (long)i * Q2, own_c = g->dc ? g->wi[i] : 1.0; mx = -1e300; s = 0;
			for (l = 0; l < Q2; l++) {
				double v = lg(f->alpha2[l]);
				for (q = 0; q < Q1; q++) v += w->a2[(long)i * Q1 + q] * U[q * Q2 + l] + own_c * w->n2[(long)i * Q1 + q] * V[q * Q2 + l];
				z[l] = v; if (v > mx) mx = v;
			}
			for (l = 0; l < Q2; l++) { z[l] = exp(z[l] - mx); s += z[l]; }
			for (l = 0; l < Q2; l++) z[l] /= s;
		}
		clamp_tau(w->l1, n1, Q1); clamp_tau(w->l2, n2, Q2);
		for (i = 0; i < n1 * Q1; i++) { double d = fabs(w->l1[i] - f->tau[i]); if (d > step) step = d; }
		for (i = 0; i < n2 * Q2; i++) { double d = fabs(w->l2[i] - f->tau2[i]); if (d > step) step = d; }
		memcpy(f->tau, w->l1, sizeof(double) * (long)n1 * Q1);
		memcpy(f->tau2, w->l2, sizeof(double) * (long)n2 * Q2);
		if (step <= TOL_F) break;
	}
	free(U); free(V);
}

/* ====================================================================================================================
   Dyadic covariates (R sbm / blockmodels *_covariates): link(E[y_ij | Z_i = q, Z_j = l]) = m_ql + sum_k beta_k x_ijk with the logit link
   (Bernoulli), log link (Poisson) or identity link (Gaussian, one common variance). The pair terms no longer reduce to block statistics, so the
   plugin keeps dense n x n2 matrices (tie values, observation indicator, covariates). Poisson and Gaussian still need only n x n2 x Q work per
   sweep (the pair term factorises); the Bernoulli model needs n x n2 x Q1 x Q2 (one log1p per pair and block pair).
   M step: Gaussian: closed form (profile of the block means, normal equations for beta); Poisson: Newton on beta with the block rates profiled out;
   Bernoulli: joint Newton on all block parameters and beta (diagonal block Hessian of the m part, Schur complement for beta), with a line search.
   ==================================================================================================================== */
static void dmm(const double *D, int nr, int nc, const double *X, int k, double *out) {   /* out (nr x k) = D (nr x nc) X (nc x k) */
	int i, j, l;
	for (i = 0; i < nr; i++) {
		double *o = out + (long)i * k;
		const double *d = D + (long)i * nc;
		for (l = 0; l < k; l++) o[l] = 0;
		for (j = 0; j < nc; j++) {
			double v = d[j];
			if (v != 0) { const double *x = X + (long)j * k; for (l = 0; l < k; l++) o[l] += v * x[l]; }
		}
	}
}
static void dmmt(const double *D, int nr, int nc, const double *X, int k, double *out) {  /* out (nc x k) = D' X, X (nr x k) */
	int i, j, l;
	memset(out, 0, sizeof(double) * (long)nc * k);
	for (i = 0; i < nr; i++) {
		const double *d = D + (long)i * nc, *x = X + (long)i * k;
		for (j = 0; j < nc; j++) {
			double v = d[j];
			if (v != 0) { double *o = out + (long)j * k; for (l = 0; l < k; l++) o[l] += v * x[l]; }
		}
	}
}
/* S (Q1 x Q2) = T1' D T2 with T1 (nr x Q1), T2 (nc x Q2); tmp has nr * Q2 doubles */
static void blockstat(const double *D, int nr, int nc, const double *T1, int Q1, const double *T2, int Q2, double *tmp, double *S) {
	int i, q, l;
	dmm(D, nr, nc, T2, Q2, tmp);
	memset(S, 0, sizeof(double) * Q1 * Q2);
	for (i = 0; i < nr; i++) for (q = 0; q < Q1; q++) { double t = T1[(long)i * Q1 + q]; for (l = 0; l < Q2; l++) S[q * Q2 + l] += t * tmp[(long)i * Q2 + l]; }
}
static int solve_lin(int p, double *a, double *b, double *x) {   /* a (p x p, destroyed), b (destroyed): x = a^-1 b */
	int i, j, k;
	for (i = 0; i < p; i++) {
		int piv = i; double mx = fabs(a[i * p + i]), t;
		for (j = i + 1; j < p; j++) if (fabs(a[j * p + i]) > mx) { mx = fabs(a[j * p + i]); piv = j; }
		if (mx < 1e-300) return 0;
		if (piv != i) { for (k = 0; k < p; k++) { t = a[i * p + k]; a[i * p + k] = a[piv * p + k]; a[piv * p + k] = t; } t = b[i]; b[i] = b[piv]; b[piv] = t; }
		for (j = i + 1; j < p; j++) { double f = a[j * p + i] / a[i * p + i]; for (k = i; k < p; k++) a[j * p + k] -= f * a[i * p + k]; b[j] -= f * b[i]; }
	}
	for (i = p - 1; i >= 0; i--) { double s = b[i]; for (k = i + 1; k < p; k++) s -= a[i * p + k] * x[k]; x[i] = s / a[i * p + i]; }
	return 1;
}
/* linear predictor part c_ij = sum_k beta_k x_ijk and the model's working matrices (Bernoulli, Poisson: exp(c) on observed cells; Gaussian: y - c) */
static void cov_refresh(const Net *g, const double *beta, Work *w) {
	long N = (long)g->n * g->n2, t; int k, p = g->p;
	for (t = 0; t < N; t++) {
		double c = 0;
		for (k = 0; k < p; k++) c += beta[k] * g->X[k][t];
		w->Cc[t] = c;
		if (g->model == 2) w->Rm[t] = g->Mf[t] != 0 ? g->Y[t] - c : 0;
		else w->Ec[t] = (g->model == 1 && g->Mf[t] == 0) ? 0 : exp(c);
	}
}
/* Tabulated softplus f(x) = log(1 + e^x) and logistic s(x) = f'(x) on a grid of 1/64 with cubic Hermite interpolation (f from f and f', s from f' and
   f'' = s (1 - s)); the interpolation error is below 1e-10. Outside [-36, 36]: f = x, s = 1, respectively f = s = e^x. Built once before any thread starts. */
#define SP_PER 64
#define SP_LIM 36
typedef struct { double f, s, d; } SpNode;   /* d = s (1 - s) */
static SpNode *sp_tab = NULL;
static int sp_init(void) {
	int k, N = 2 * SP_LIM * SP_PER + 2;
	if (sp_tab) return 1;
	sp_tab = malloc(sizeof(SpNode) * N);
	if (!sp_tab) return 0;
	for (k = 0; k < N; k++) {
		double x = -SP_LIM + (double)k / SP_PER, e = exp(-fabs(x)), sg;
		sp_tab[k].f = (x > 0 ? x : 0) + log1p(e);
		sg = x >= 0 ? 1 / (1 + e) : e / (1 + e);
		sp_tab[k].s = sg; sp_tab[k].d = sg * (1 - sg);
	}
	return 1;
}
#ifdef _MSC_VER
#define SBM_INLINE __inline
#else
#define SBM_INLINE inline
#endif
static SBM_INLINE void sp_eval(double x, double *f, double *s) {
	if (x >= SP_LIM) { *f = x; *s = 1; return; }
	if (x <= -SP_LIM) { double e = exp(x); *f = e; *s = e; return; }
	{
		double u = (x + SP_LIM) * SP_PER, t, t2, t3, h = 1.0 / SP_PER;
		int k = (int)u; const SpNode *a, *b;
		t = u - k; t2 = t * t; t3 = t2 * t; a = sp_tab + k; b = a + 1;
		*f = (2 * t3 - 3 * t2 + 1) * a->f + (t3 - 2 * t2 + t) * h * a->s + (3 * t2 - 2 * t3) * b->f + (t3 - t2) * h * b->s;
		*s = (2 * t3 - 3 * t2 + 1) * a->s + (t3 - 2 * t2 + t) * h * a->d + (3 * t2 - 2 * t3) * b->s + (t3 - t2) * h * b->d;
	}
}

/* Bernoulli covariate model, E step scores of both sides in one pass over the observed cells: z1 (n1 x Q1, nodes as rows) and, for a two-mode
   network, z2 (n2 x Q2); a one-mode network has z2 = z1 (cell (i, j) feeds node i as the row and node j as the column; an undirected network visits
   every pair once because its block matrix and its covariates are symmetric). score(i, q) = sum_j sum_l tau_jl [y_ij m_ql - softplus(m_ql + c_ij)]. */
static void e_scores_bern(const Net *g, const Fit *f, Work *w, double *z1, double *z2) {
	int n1 = g->n, n2 = g->n2, Q1 = f->Q, Q2 = g->kind ? f->Q2 : f->Q, i, j, q, l;
	const double *ta = f->tau, *tb = g->kind ? f->tau2 : f->tau, *m = f->mq;
	double *zb = g->kind ? z2 : z1;
	int und = (!g->kind && !g->directed);
	memset(z1, 0, sizeof(double) * (size_t)n1 * Q1);
	if (g->kind) memset(z2, 0, sizeof(double) * (size_t)n2 * Q2);
	for (i = 0; i < n1; i++) {
		const double *ti = ta + (size_t)i * Q1; double *zi = z1 + (size_t)i * Q1;
		for (j = und ? i + 1 : 0; j < n2; j++) {
			size_t t = (size_t)i * n2 + j; double y, c0; const double *tj; double *zj;
			if (g->Mf[t] == 0) continue;
			y = g->Y[t]; c0 = w->Cc[t]; tj = tb + (size_t)j * Q2; zj = zb + (size_t)j * Q2;
			for (q = 0; q < Q1; q++) {
				double acc = 0, tq = ti[q]; const double *mr = m + (size_t)q * Q2;
				for (l = 0; l < Q2; l++) {
					double sp, sg, F;
					sp_eval(mr[l] + c0, &sp, &sg);
					F = y * mr[l] - sp;
					acc += tj[l] * F; zj[l] += tq * F;
				}
				zi[q] += acc;
			}
		}
	}
	(void)n2;
}

/* E-step scores of one side: role 0 = nodes are rows, role 1 = nodes are columns; tauO: memberships of the other side; z (nS x Qs) = sum over the other side of the expected log density */
static void cov_scores(const Net *g, const Fit *f, Work *w, int role, const double *tauO, double *z) {
	int n1 = g->n, n2 = g->n2, Q1 = f->Q, Q2 = g->kind ? f->Q2 : f->Q;
	int nS = role ? n2 : n1, Qs = role ? Q2 : Q1, QO = role ? Q1 : Q2;
	int i, q, l;
	double *mm = malloc(sizeof(double) * Qs * QO), *EE = malloc(sizeof(double) * Qs * QO);
	for (q = 0; q < Qs; q++) for (l = 0; l < QO; l++) { int idx = role ? l * Q2 + q : q * Q2 + l; mm[q * QO + l] = f->mq[idx]; EE[q * QO + l] = g->model == 2 ? 0 : exp(f->mq[idx]); }
	if (g->model == 1) {
		spmm(role ? &g->ai : &g->ao, nS, tauO, QO, w->a1, 1, NULL);
		if (role) dmmt(w->Ec, n1, n2, tauO, QO, w->n1); else dmm(w->Ec, n1, n2, tauO, QO, w->n1);
		for (i = 0; i < nS; i++) for (q = 0; q < Qs; q++) {
			double acc = 0;
			for (l = 0; l < QO; l++) acc += w->a1[(long)i * QO + l] * mm[q * QO + l] - EE[q * QO + l] * w->n1[(long)i * QO + l];
			z[(long)i * Qs + q] = acc;
		}
	} else {
		double s2 = f->sigma2;
		if (role) { dmmt(w->Rm, n1, n2, tauO, QO, w->a1); dmmt(g->Mf, n1, n2, tauO, QO, w->n1); }
		else { dmm(w->Rm, n1, n2, tauO, QO, w->a1); dmm(g->Mf, n1, n2, tauO, QO, w->n1); }
		for (i = 0; i < nS; i++) for (q = 0; q < Qs; q++) {
			double acc = 0;
			for (l = 0; l < QO; l++) { double m = mm[q * QO + l]; acc += w->a1[(long)i * QO + l] * m - 0.5 * w->n1[(long)i * QO + l] * m * m; }
			z[(long)i * Qs + q] = acc / s2;
		}
	}
	free(mm); free(EE);
}

static void softmax_rows(double *z, int n, int Q, const double *alpha) {
	int i, q;
	for (i = 0; i < n; i++) {
		double mx = -1e300, s = 0, *zi = z + (long)i * Q;
		for (q = 0; q < Q; q++) { zi[q] += lg(alpha[q]); if (zi[q] > mx) mx = zi[q]; }
		for (q = 0; q < Q; q++) { zi[q] = exp(zi[q] - mx); s += zi[q]; }
		for (q = 0; q < Q; q++) zi[q] /= s;
	}
}

static void e_step_cov(const Net *g, Fit *f, Work *w) {
	int n1 = g->n, n2 = g->n2, Q1 = f->Q, Q2 = g->kind ? f->Q2 : f->Q, it, i;
	for (it = 0; it < 10; it++) {
		double step = 0;
		if (g->model == 0) {
			e_scores_bern(g, f, w, w->l1, w->l2);
		} else {
			cov_scores(g, f, w, 0, g->kind ? f->tau2 : f->tau, w->l1);
			if (!g->kind && g->directed) {
				cov_scores(g, f, w, 1, f->tau, w->l2);
				for (i = 0; i < n1 * Q1; i++) w->l1[i] += w->l2[i];
			}
		}
		softmax_rows(w->l1, n1, Q1, f->alpha);
		clamp_tau(w->l1, n1, Q1);
		if (g->kind) {
			if (g->model != 0) cov_scores(g, f, w, 1, f->tau, w->l2);
			softmax_rows(w->l2, n2, Q2, f->alpha2);
			clamp_tau(w->l2, n2, Q2);
		}
		for (i = 0; i < n1 * Q1; i++) { double d = fabs(w->l1[i] - f->tau[i]); if (d > step) step = d; }
		if (g->kind) for (i = 0; i < n2 * Q2; i++) { double d = fabs(w->l2[i] - f->tau2[i]); if (d > step) step = d; }
		memcpy(f->tau, w->l1, sizeof(double) * (long)n1 * Q1);
		if (g->kind) memcpy(f->tau2, w->l2, sizeof(double) * (long)n2 * Q2);
		if (step <= TOL_F) break;
	}
}

/* Bernoulli: objective sum_ij sum_ql tau_iq tau_jl [y eta - log(1 + e^eta)] and (optionally) its gradient and Hessian blocks */
static double bern_eval(const Net *g, const Fit *f, Work *w, const double *mq, const double *beta, int deriv,
	double *gm, double *hm, double *gb, double *hbm, double *hbb) {
	int n1 = g->n, n2 = g->n2, Q1 = f->Q, Q2 = g->kind ? f->Q2 : f->Q, p = g->p, cells = Q1 * Q2, i, j, q, l, k, kk;
	const double *ta = f->tau, *tb = g->kind ? f->tau2 : f->tau;
	double o = 0, xs[MAXCOV];
	int und = (!g->kind && !g->directed);
	cov_refresh(g, beta, w);
	if (deriv) { memset(gm, 0, sizeof(double) * cells); memset(hm, 0, sizeof(double) * cells); memset(gb, 0, sizeof(double) * p); memset(hbm, 0, sizeof(double) * cells * p); memset(hbb, 0, sizeof(double) * p * p); }
	for (i = 0; i < n1; i++) {
		const double *t1 = ta + (long)i * Q1;
		for (j = und ? i + 1 : 0; j < n2; j++) {
			long t = (long)i * n2 + j; double y, c0, sg = 0, sh = 0; const double *t2;
			if (g->Mf[t] == 0) continue;
			y = g->Y[t]; c0 = w->Cc[t]; t2 = tb + (long)j * Q2;
			if (deriv) for (k = 0; k < p; k++) xs[k] = g->X[k][t];
			for (q = 0; q < Q1; q++) {
				double tq = t1[q];
				for (l = 0; l < Q2; l++) {
					int idx = q * Q2 + l; double ww = und ? tq * t2[l] + t2[q] * t1[l] : tq * t2[l], eta = mq[idx] + c0, sp, s;
					sp_eval(eta, &sp, &s);
					o += ww * (y * eta - sp);
					if (deriv) {
						double rr = ww * (y - s), hh = ww * s * (1 - s);
						gm[idx] += rr; hm[idx] += hh; sg += rr; sh += hh;
						for (k = 0; k < p; k++) hbm[idx * p + k] += hh * xs[k];
					}
				}
			}
			if (deriv) for (k = 0; k < p; k++) { gb[k] += sg * xs[k]; for (kk = 0; kk <= k; kk++) hbb[k * p + kk] += sh * xs[k] * xs[kk]; }
		}
	}
	if (deriv) for (k = 0; k < p; k++) for (kk = k + 1; kk < p; kk++) hbb[k * p + kk] = hbb[kk * p + k];
	return o;
}

static double m_step_cov(const Net *g, Fit *f, Work *w) {
	int n1 = g->n, n2 = g->n2, Q1 = f->Q, kind = g->kind, Q2 = kind ? f->Q2 : f->Q, p = g->p, cells = Q1 * Q2, i, q, l, k, kk, c, it;
	const double *t1 = f->tau, *t2 = kind ? f->tau2 : f->tau;
	double *ZA = calloc(cells, sizeof(double)), *ZO = calloc(cells, sizeof(double)), pl = 0, half = (!kind && !g->directed) ? 0.5 : 1.0;
	blockstat(g->Y, n1, n2, t1, Q1, t2, Q2, w->tmp, ZA);
	blockstat(g->Mf, n1, n2, t1, Q1, t2, Q2, w->tmp, ZO);
	if (g->model == 2) {
		double *ZX = malloc(sizeof(double) * p * cells), *S = malloc(sizeof(double) * p * p), *tv = malloc(sizeof(double) * p), *bx = malloc(sizeof(double) * p);
		double sumres, resid, s2;
		for (k = 0; k < p; k++) blockstat(g->X[k], n1, n2, t1, Q1, t2, Q2, w->tmp, ZX + k * cells);
		for (k = 0; k < p; k++) {
			for (kk = 0; kk < p; kk++) {
				double v = g->sxx[k * p + kk];
				for (c = 0; c < cells; c++) if (ZO[c] > 0) v -= ZX[k * cells + c] * ZX[kk * cells + c] / ZO[c];
				S[k * p + kk] = v;
			}
			tv[k] = g->sxy[k];
			for (c = 0; c < cells; c++) if (ZO[c] > 0) tv[k] -= ZX[k * cells + c] * ZA[c] / ZO[c];
			S[k * p + k] += 1e-12 * (fabs(S[k * p + k]) + 1);
		}
		for (k = 0; k < p * p; k++) f->info[k] = S[k];
		if (!solve_lin(p, S, tv, bx)) for (k = 0; k < p; k++) bx[k] = f->beta[k];
		for (k = 0; k < p; k++) f->beta[k] = bx[k];
		sumres = g->sumy2;
		for (k = 0; k < p; k++) { sumres -= 2 * f->beta[k] * g->sxy[k]; for (kk = 0; kk < p; kk++) sumres += f->beta[k] * f->beta[kk] * g->sxx[k * p + kk]; }
		resid = sumres;
		for (c = 0; c < cells; c++) {
			double zr = ZA[c];
			for (k = 0; k < p; k++) zr -= f->beta[k] * ZX[k * cells + c];
			f->mq[c] = ZO[c] > 0 ? zr / ZO[c] : 0.0;
			resid += -2 * f->mq[c] * zr + f->mq[c] * f->mq[c] * ZO[c];
		}
		if (resid < 0) resid = 0;
		s2 = g->sigma2fix > 0 ? g->sigma2fix : resid / g->nobs;
		if (s2 < 1e-12) s2 = 1e-12;
		f->sigma2 = s2;
		for (k = 0; k < p * p; k++) f->info[k] *= half / s2;
		pl = -0.5 * g->nobs * log(2 * M_PI * s2) - resid / (2 * s2);
		for (c = 0; c < cells; c++) f->pi[c] = f->mq[c];
		cov_refresh(g, f->beta, w);
		free(ZX); free(S); free(tv); free(bx);
	} else if (g->model == 1) {
		double *ZE = malloc(sizeof(double) * cells), *ZEX = malloc(sizeof(double) * p * cells), *ZEXX = malloc(sizeof(double) * p * p * cells);
		double *gb = malloc(sizeof(double) * p), *H = malloc(sizeof(double) * p * p), *dl = malloc(sizeof(double) * p), *bt = malloc(sizeof(double) * p), *Ht = malloc(sizeof(double) * p * p);
		double obj = 0, objt;
		long N = (long)n1 * n2, t;
		for (it = 0; it < 40; it++) {
			double mxs = 0, tt = 1; int ok = 0, ls;
			cov_refresh(g, f->beta, w);
			blockstat(w->Ec, n1, n2, t1, Q1, t2, Q2, w->tmp, ZE);
			obj = 0;
			for (c = 0; c < cells; c++) if (ZA[c] > 1e-12 && ZE[c] > 0) obj -= ZA[c] * log(ZE[c]);
			for (k = 0; k < p; k++) obj += f->beta[k] * g->sxy[k];
			for (k = 0; k < p; k++) {
				for (t = 0; t < N; t++) w->Rm[t] = w->Ec[t] * g->X[k][t];
				blockstat(w->Rm, n1, n2, t1, Q1, t2, Q2, w->tmp, ZEX + k * cells);
				for (kk = 0; kk <= k; kk++) {
					for (t = 0; t < N; t++) w->Rm[t] = w->Ec[t] * g->X[k][t] * g->X[kk][t];
					blockstat(w->Rm, n1, n2, t1, Q1, t2, Q2, w->tmp, ZEXX + (k * p + kk) * cells);
				}
			}
			for (k = 0; k < p; k++) {
				gb[k] = g->sxy[k];
				for (c = 0; c < cells; c++) if (ZA[c] > 1e-12 && ZE[c] > 0) gb[k] -= ZA[c] * ZEX[k * cells + c] / ZE[c];
				for (kk = 0; kk <= k; kk++) {
					double v = 0;
					for (c = 0; c < cells; c++) if (ZA[c] > 1e-12 && ZE[c] > 0) v += ZA[c] * (ZEXX[(k * p + kk) * cells + c] / ZE[c] - ZEX[k * cells + c] * ZEX[kk * cells + c] / (ZE[c] * ZE[c]));
					H[k * p + kk] = H[kk * p + k] = v;
				}
				H[k * p + k] += 1e-12 * (fabs(H[k * p + k]) + 1);
			}
			memcpy(Ht, H, sizeof(double) * p * p);
			if (!solve_lin(p, Ht, gb, dl)) break;
			for (k = 0; k < p; k++) if (fabs(dl[k]) > mxs) mxs = fabs(dl[k]);
			if (mxs < 1e-9) break;
			for (ls = 0; ls < 30 && !ok; ls++) {
				for (k = 0; k < p; k++) bt[k] = f->beta[k] + tt * dl[k];
				cov_refresh(g, bt, w);
				blockstat(w->Ec, n1, n2, t1, Q1, t2, Q2, w->tmp, ZE);
				objt = 0;
				for (c = 0; c < cells; c++) if (ZA[c] > 1e-12 && ZE[c] > 0) objt -= ZA[c] * log(ZE[c]);
				for (k = 0; k < p; k++) objt += bt[k] * g->sxy[k];
				if (objt >= obj - 1e-10 * (fabs(obj) + 1)) ok = 1; else tt *= 0.5;
			}
			if (!ok) break;
			for (k = 0; k < p; k++) f->beta[k] = bt[k];
			if (mxs * tt < 1e-9) break;
		}
		cov_refresh(g, f->beta, w);
		blockstat(w->Ec, n1, n2, t1, Q1, t2, Q2, w->tmp, ZE);
		for (k = 0; k < p; k++) {
			for (t = 0; t < N; t++) w->Rm[t] = w->Ec[t] * g->X[k][t];
			blockstat(w->Rm, n1, n2, t1, Q1, t2, Q2, w->tmp, ZEX + k * cells);
			for (kk = 0; kk <= k; kk++) {
				double v = 0;
				for (t = 0; t < N; t++) w->Rm[t] = w->Ec[t] * g->X[k][t] * g->X[kk][t];
				blockstat(w->Rm, n1, n2, t1, Q1, t2, Q2, w->tmp, ZEXX + (k * p + kk) * cells);
				for (c = 0; c < cells; c++) if (ZA[c] > 1e-12 && ZE[c] > 0) v += ZA[c] * (ZEXX[(k * p + kk) * cells + c] / ZE[c] - ZEX[k * cells + c] * ZEX[kk * cells + c] / (ZE[c] * ZE[c]));
				f->info[k * p + kk] = f->info[kk * p + k] = v * half;
			}
		}
		for (c = 0; c < cells; c++) {
			double za = ZA[c] > 1e-12 ? ZA[c] : 1e-12;
			f->mq[c] = ZE[c] > 0 ? log(za / ZE[c]) : 0.0;
			f->pi[c] = exp(f->mq[c]);
			if (ZA[c] > 0) pl += ZA[c] * f->mq[c];
			pl -= f->pi[c] * ZE[c];
		}
		for (k = 0; k < p; k++) pl += f->beta[k] * g->sxy[k];
		free(ZE); free(ZEX); free(ZEXX); free(gb); free(H); free(dl); free(bt); free(Ht);
	} else {
		double *gm = malloc(sizeof(double) * cells), *hm = malloc(sizeof(double) * cells), *gb = malloc(sizeof(double) * p), *hbm = malloc(sizeof(double) * cells * p), *hbb = malloc(sizeof(double) * p * p);
		double *gm2 = malloc(sizeof(double) * cells), *hm2 = malloc(sizeof(double) * cells), *gb2 = malloc(sizeof(double) * p), *hbm2 = malloc(sizeof(double) * cells * p), *hbb2 = malloc(sizeof(double) * p * p);
		double *S = malloc(sizeof(double) * p * p), *rhs = malloc(sizeof(double) * p), *db = malloc(sizeof(double) * p), *dm = malloc(sizeof(double) * cells);
		double *mt = malloc(sizeof(double) * cells), *bt = malloc(sizeof(double) * p), obj, objt;
		if (!f->warm) {
			for (c = 0; c < cells; c++) { double pp = (ZA[c] + 0.5) / (ZO[c] + 1); f->mq[c] = log(pp / (1 - pp)); }
			for (k = 0; k < p; k++) f->beta[k] = 0;
		}
		obj = bern_eval(g, f, w, f->mq, f->beta, 1, gm, hm, gb, hbm, hbb);
		for (it = 0; it < (f->warm ? f->nnewton : 40); it++) {
			double mxs = 0, tt = 1; int ok = 0, ls;
			memcpy(S, hbb, sizeof(double) * p * p);
			for (k = 0; k < p; k++) rhs[k] = gb[k];
			for (c = 0; c < cells; c++) {
				double d = hm[c] + 1e-9;
				for (k = 0; k < p; k++) { rhs[k] -= hbm[c * p + k] * gm[c] / d; for (kk = 0; kk < p; kk++) S[k * p + kk] -= hbm[c * p + k] * hbm[c * p + kk] / d; }
			}
			for (k = 0; k < p; k++) S[k * p + k] += 1e-10 * (fabs(S[k * p + k]) + 1);
			{ double *St = malloc(sizeof(double) * p * p), *rt = malloc(sizeof(double) * p); memcpy(St, S, sizeof(double) * p * p); memcpy(rt, rhs, sizeof(double) * p);
			  if (!solve_lin(p, St, rt, db)) { for (k = 0; k < p; k++) db[k] = 0; }
			  free(St); free(rt); }
			for (c = 0; c < cells; c++) {
				double d = hm[c] + 1e-9, v = gm[c];
				for (k = 0; k < p; k++) v -= hbm[c * p + k] * db[k];
				dm[c] = v / d;
				if (fabs(dm[c]) > mxs) mxs = fabs(dm[c]);
			}
			for (k = 0; k < p; k++) if (fabs(db[k]) > mxs) mxs = fabs(db[k]);
			if (mxs < 1e-8) break;
			for (ls = 0; ls < 30 && !ok; ls++) {
				for (c = 0; c < cells; c++) { mt[c] = f->mq[c] + tt * dm[c]; if (mt[c] > 30) mt[c] = 30; if (mt[c] < -30) mt[c] = -30; }
				for (k = 0; k < p; k++) bt[k] = f->beta[k] + tt * db[k];
				objt = bern_eval(g, f, w, mt, bt, 1, gm2, hm2, gb2, hbm2, hbb2);
				if (objt >= obj - 1e-10 * (fabs(obj) + 1)) ok = 1; else tt *= 0.5;
			}
			if (!ok) break;
			{ double *tp; double ob = objt - obj;
			  memcpy(f->mq, mt, sizeof(double) * cells); memcpy(f->beta, bt, sizeof(double) * p); obj = objt;
			  tp = gm; gm = gm2; gm2 = tp; tp = hm; hm = hm2; hm2 = tp; tp = gb; gb = gb2; gb2 = tp; tp = hbm; hbm = hbm2; hbm2 = tp; tp = hbb; hbb = hbb2; hbb2 = tp;
			  if (ob < 1e-11 * (fabs(obj) + 1) && mxs * tt < 1e-6) break; }
		}
		{ /* information of beta at the final point (gm, hm, hbm, hbb belong to it): Schur complement of the block part */
		  memcpy(S, hbb, sizeof(double) * p * p);
		  for (c = 0; c < cells; c++) { double d = hm[c] + 1e-9; for (k = 0; k < p; k++) for (kk = 0; kk < p; kk++) S[k * p + kk] -= hbm[c * p + k] * hbm[c * p + kk] / d; }
		  for (k = 0; k < p * p; k++) f->info[k] = S[k] * half; }
		cov_refresh(g, f->beta, w);
		pl = obj;
		for (c = 0; c < cells; c++) f->pi[c] = 1 / (1 + exp(-f->mq[c]));
		f->warm = 1;
		free(gm); free(hm); free(gb); free(hbm); free(hbb); free(gm2); free(hm2); free(gb2); free(hbm2); free(hbb2); free(S); free(rhs); free(db); free(dm); free(mt); free(bt);
	}
	if (g->model != 0) f->warm = 1;
	/* proportions, entropy, label terms */
	for (q = 0; q < Q1; q++) { double s = 0; for (i = 0; i < n1; i++) s += f->tau[(long)i * Q1 + q]; f->alpha[q] = s / n1; }
	if (kind) for (l = 0; l < Q2; l++) { double s = 0; for (i = 0; i < n2; i++) s += f->tau2[(long)i * Q2 + l]; f->alpha2[l] = s / n2; }
	pl = pl * half + g->cst;
	for (i = 0; i < n1; i++) for (q = 0; q < Q1; q++) pl += f->tau[(long)i * Q1 + q] * lg(f->alpha[q]);
	if (kind) for (i = 0; i < n2; i++) for (l = 0; l < Q2; l++) pl += f->tau2[(long)i * Q2 + l] * lg(f->alpha2[l]);
	{ double h = 0; for (i = 0; i < n1 * Q1; i++) if (f->tau[i] > 0) h -= f->tau[i] * log(f->tau[i]);
	  if (kind) for (i = 0; i < n2 * Q2; i++) if (f->tau2[i] > 0) h -= f->tau2[i] * log(f->tau2[i]);
	  f->H = h; }
	free(ZA); free(ZO);
	return pl;
}

static double npar_of(const Net *g, int Q, int Q2) {
	double np = g->kind ? (double)Q * Q2 : (g->directed ? (double)Q * Q : (double)Q * (Q + 1) / 2);
	if (g->model == 2 && g->sigma2fix <= 0) np += 1;
	return np + g->p;
}
static double icl_of(const Net *g, const Fit *f) {
	double pen = (f->Q - 1) * log((double)g->n);
	if (g->kind) pen += (f->Q2 - 1) * log((double)g->n2);
	return f->PL - 0.5 * (f->npar * log(g->ndata) + pen);
}

/* full VEM from the start tau (and tau2) in f */
static void vem_tol(const Net *g, Fit *f, Work *w, double tol) {
	int it, trivial = g->kind ? (f->Q == 1 && f->Q2 == 1) : (f->Q == 1);
	double J, old;
	clamp_tau(f->tau, g->n, f->Q);
	if (g->kind) clamp_tau(f->tau2, g->n2, f->Q2);
	f->PL = g->p ? m_step_cov(g, f, w) : g->kind ? m_step_lbm(g, f, w) : m_step_sbm(g, f, w);
	J = f->PL + f->H;
	if (!trivial) {
		for (it = 0; it < MAXIT; it++) {
			if (g->p) e_step_cov(g, f, w); else if (g->kind) e_step_lbm(g, f, w); else e_step_sbm(g, f, w);
			f->PL = g->p ? m_step_cov(g, f, w) : g->kind ? m_step_lbm(g, f, w) : m_step_sbm(g, f, w);
			old = J; J = f->PL + f->H;
			if (!(J - old > tol)) break;
		}
	}
	f->J = J;
	f->npar = npar_of(g, f->Q, f->Q2);
	f->icl = icl_of(g, f);
}
static void vem(const Net *g, Fit *f, Work *w) { vem_tol(g, f, w, TOL_EM); }

/* residual matrix of the one-block covariate model (the starts use its spectral embedding, as the residual of the one-block model without covariates) */
static int cov_residual(Net *g) {
	Fit f; Work w; long N = (long)g->n * g->n2, t; double m0;
	long mx = g->n > g->n2 ? g->n : g->n2;
	if (!work_alloc(&w, mx + 1, g)) return 0;
	fit_alloc(&f, g, 1, g->kind ? 1 : 0);
	{ long i; for (i = 0; i < g->n; i++) f.tau[i] = 1; if (g->kind) for (i = 0; i < g->n2; i++) f.tau2[i] = 1; }
	vem(g, &f, &w);
	m0 = f.mq[0];
	g->Rres = malloc(sizeof(double) * N);
	if (!g->Rres) return 0;
	for (t = 0; t < N; t++) {
		double eta = m0 + w.Cc[t], mu = g->model == 0 ? 1 / (1 + exp(-eta)) : g->model == 1 ? exp(eta) : eta;
		g->Rres[t] = g->Mf[t] != 0 ? g->Y[t] - mu : 0;
	}
	fit_free(&f); work_free(&w);
	return 1;
}

/* ---------- spectral embedding ---------- */
typedef struct { int k; double *V; double *C; /* V: rows (n x k); C: columns of a two-mode matrix (n2 x k) */ } Embed;

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

static void orthonormalise(double *V, int n, int k, Rng *rng) { /* modified Gram-Schmidt on columns of V (n x k row major) */
	int a, b, i;
	for (a = 0; a < k; a++) {
		double s;
		for (b = 0; b < a; b++) {
			double d = 0; for (i = 0; i < n; i++) d += V[(long)i * k + a] * V[(long)i * k + b];
			for (i = 0; i < n; i++) V[(long)i * k + a] -= d * V[(long)i * k + b];
		}
		s = 0; for (i = 0; i < n; i++) s += V[(long)i * k + a] * V[(long)i * k + a];
		s = sqrt(s);
		if (s < 1e-300) { for (i = 0; i < n; i++) V[(long)i * k + a] = rng_u(rng) - 0.5; s = 0; for (i = 0; i < n; i++) s += V[(long)i * k + a] * V[(long)i * k + a]; s = sqrt(s); }
		for (i = 0; i < n; i++) V[(long)i * k + a] /= s;
	}
}

/* residual operator on the row space: one-mode B = (A + A')/2 - p M (degree corrected: minus the expected degree product); two-mode B B' with B = A - p M */
static void emb_apply(const Net *g, double p, const double *V, int k, double *out) {
	int n = g->n, i, l;
	if (g->Rres) {
		if (g->kind) {
			int n2 = g->n2;
			double *mid = malloc(sizeof(double) * (long)n2 * k);
			dmmt(g->Rres, n, n2, V, k, mid);
			dmm(g->Rres, n, n2, mid, k, out);
			free(mid);
		} else {
			double *t1 = malloc(sizeof(double) * (long)n * k);
			dmm(g->Rres, n, n, V, k, out);
			if (g->directed) { dmmt(g->Rres, n, n, V, k, t1); for (i = 0; i < n * k; i++) out[i] = 0.5 * (out[i] + t1[i]); }
			free(t1);
		}
		return;
	}
	if (g->kind) {
		int n2 = g->n2;
		double *mid = malloc(sizeof(double) * (long)n2 * k), *t = malloc(sizeof(double) * (long)n2 * k);
		double *mm = malloc(sizeof(double) * (long)n * k);
		spmm(&g->ai, n2, V, k, mid, 1, NULL);                    /* A' V */
		if (g->dc) {                                             /* degree corrected: B = A - wo wi' / sum(w) */
			double *cv = calloc(k, sizeof(double)), *cm = calloc(k, sizeof(double));
			for (i = 0; i < n; i++) for (l = 0; l < k; l++) cv[l] += g->wo[i] * V[(long)i * k + l];
			for (i = 0; i < n2; i++) for (l = 0; l < k; l++) mid[(long)i * k + l] -= g->wi[i] * cv[l] / g->sumw;
			spmm(&g->ao, n, mid, k, out, 1, NULL);
			for (i = 0; i < n2; i++) for (l = 0; l < k; l++) cm[l] += g->wi[i] * mid[(long)i * k + l];
			for (i = 0; i < n; i++) for (l = 0; l < k; l++) out[(long)i * k + l] -= g->wo[i] * cm[l] / g->sumw;
			free(cv); free(cm); free(mid); free(t); free(mm);
			return;
		}
		obs_sum_rect(g, &g->mi, n2, V, n, k, t, NULL);             /* M' V */
		for (i = 0; i < n2 * k; i++) mid[i] -= p * t[i];
		spmm(&g->ao, n, mid, k, out, 1, NULL);                    /* A (B' V) */
		obs_sum_rect(g, &g->mo, n, mid, n2, k, mm, NULL);         /* M (B' V) */
		for (i = 0; i < n * k; i++) out[i] -= p * mm[i];
		free(mid); free(t); free(mm);
		return;
	}
	{
		double *t1 = malloc(sizeof(double) * (long)n * k), *t2 = malloc(sizeof(double) * (long)n * k);
		spmm(&g->ao, n, V, k, out, 1, NULL);
		if (g->directed) {
			spmm(&g->ai, n, V, k, t1, 1, NULL);
			for (i = 0; i < n * k; i++) out[i] = 0.5 * (out[i] + t1[i]);
		}
		if (g->dc) {
			double *cw = calloc(k, sizeof(double)), *cv = calloc(k, sizeof(double));
			for (i = 0; i < n; i++) for (l = 0; l < k; l++) { cw[l] += g->wi[i] * V[(long)i * k + l]; cv[l] += g->wo[i] * V[(long)i * k + l]; }
			for (i = 0; i < n; i++) for (l = 0; l < k; l++) out[(long)i * k + l] -= 0.5 * (g->wo[i] * cw[l] + g->wi[i] * cv[l]) / g->sumw;
			free(cw); free(cv);
		} else {
			obs_sum(g, &g->mo, V, k, t2, NULL);
			if (g->directed) { obs_sum(g, &g->mi, V, k, t1, NULL); for (i = 0; i < n * k; i++) t2[i] = 0.5 * (t2[i] + t1[i]); }
			for (i = 0; i < n * k; i++) out[i] -= p * t2[i];
		}
		free(t1); free(t2);
	}
}

static int embed_compute(const Net *g, double p, int k, Embed *em, Rng *rng) {
	int n = g->n, it, a, b, i;
	double *V = malloc(sizeof(double) * (long)n * k), *W = malloc(sizeof(double) * (long)n * k);
	double *H = malloc(sizeof(double) * k * k), *E = malloc(sizeof(double) * k * k), *lam = malloc(sizeof(double) * k), *prev = calloc(k, sizeof(double));
	if (!V || !W || !H || !E || !lam || !prev) return 0;
	for (i = 0; i < n * k; i++) V[i] = rng_u(rng) - 0.5;
	orthonormalise(V, n, k, rng);
	for (it = 0; it < 500; it++) {
		double ch = 0;
		emb_apply(g, p, V, k, W);
		for (a = 0; a < k; a++) for (b = 0; b < k; b++) { double d = 0; for (i = 0; i < n; i++) d += V[(long)i * k + a] * W[(long)i * k + b]; H[a * k + b] = d; }
		for (a = 0; a < k; a++) for (b = a + 1; b < k; b++) { double m = 0.5 * (H[a * k + b] + H[b * k + a]); H[a * k + b] = H[b * k + a] = m; }
		jacobi_eig(k, H, E, lam);
		for (a = 0; a < k; a++) { double d = fabs(fabs(lam[a]) - prev[a]); if (d > ch * fabs(lam[a])) ch = d / (fabs(lam[a]) + 1e-12); prev[a] = fabs(lam[a]); }
		memcpy(V, W, sizeof(double) * (long)n * k);
		orthonormalise(V, n, k, rng);
		if (it > 3 && ch < 1e-7) break;
	}
	emb_apply(g, p, V, k, W);
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
		em->k = k; em->V = R; em->C = NULL; free(ord);
	}
	if (g->kind) { /* column embedding: C = B' U, columns scaled to unit length */
		int n2 = g->n2; double *t = malloc(sizeof(double) * (long)n2 * k);
		em->C = malloc(sizeof(double) * (long)n2 * k);
		if (g->Rres) { dmmt(g->Rres, n, n2, em->V, k, em->C); memset(t, 0, sizeof(double) * (long)n2 * k); }
		else {
			spmm(&g->ai, n2, em->V, k, em->C, 1, NULL);
			if (g->dc) {
				double *cv = calloc(k, sizeof(double)); int l;
				for (i = 0; i < n; i++) for (l = 0; l < k; l++) cv[l] += g->wo[i] * em->V[(long)i * k + l];
				for (i = 0; i < n2; i++) for (l = 0; l < k; l++) em->C[(long)i * k + l] -= g->wi[i] * cv[l] / g->sumw;
				free(cv);
			} else {
				obs_sum_rect(g, &g->mi, n2, em->V, n, k, t, NULL);
				for (i = 0; i < n2 * k; i++) em->C[i] -= p * t[i];
			}
		}
		for (a = 0; a < k; a++) {
			double s = 0; for (i = 0; i < n2; i++) s += em->C[(long)i * k + a] * em->C[(long)i * k + a];
			s = sqrt(s); if (s < 1e-300) s = 1;
			for (i = 0; i < n2; i++) em->C[(long)i * k + a] /= s;
		}
		free(t);
	}
	if (g->dc) for (i = 0; i < n; i++) { /* degree correction: rows of the embedding on the unit sphere */
		double s = 0; for (a = 0; a < k; a++) s += em->V[(long)i * k + a] * em->V[(long)i * k + a];
		s = sqrt(s); if (s > 1e-300) for (a = 0; a < k; a++) em->V[(long)i * k + a] /= s;
	}
	if (g->dc && g->kind) for (i = 0; i < g->n2; i++) {
		double s = 0; for (a = 0; a < k; a++) s += em->C[(long)i * k + a] * em->C[(long)i * k + a];
		s = sqrt(s); if (s > 1e-300) for (a = 0; a < k; a++) em->C[(long)i * k + a] /= s;
	}
	free(V); free(W); free(H); free(E); free(lam); free(prev);
	return 1;
}

/* ---------- k-means on the first d columns of the embedding (row stride k) ---------- */
static double sqd(const double *x, const double *y, int d) { double s = 0; int a; for (a = 0; a < d; a++) { double t = x[a] - y[a]; s += t * t; } return s; }

/* mode 0: farthest point start (deterministic); mode 1: k-means++ (random); idx: optional subset of rows (m of them) */
static void kmeans(const double *V, int k, int d, const int *idx, int m, int Q, int mode, int *lab, Rng *rng) {
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
			int farj = 0; best = -1;
			for (i = 0; i < m; i++) if (mind[i] > best) { best = mind[i]; farj = i; }
			memcpy(cen + c * d, ROW(farj), sizeof(double) * d);
			for (i = 0; i < m; i++) { double s = sqd(ROW(i), cen + c * d, d); if (s < mind[i]) mind[i] = s; }
		}
	} else {
		int first = rng_int(rng, m);
		memcpy(cen, ROW(first), sizeof(double) * d);
		for (i = 0; i < m; i++) mind[i] = sqd(ROW(i), cen, d);
		for (c = 1; c < Q; c++) {
			double tot = 0, u; int pick = m - 1;
			for (i = 0; i < m; i++) tot += mind[i];
			u = rng_u(rng) * tot;
			if (tot <= 0) pick = rng_int(rng, m);
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
			if (cnt[c] == 0) { int r = rng_int(rng, m); memcpy(cen + c * d, ROW(r), sizeof(double) * d); }
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

static void hard_labels(const double *tau, int n, int Q, int *lab) {
	int i, k;
	for (i = 0; i < n; i++) {
		int bq = 0; double bv = -1;
		for (k = 0; k < Q; k++) if (tau[(long)i * Q + k] > bv) { bv = tau[(long)i * Q + k]; bq = k; }
		lab[i] = bq;
	}
}

/* ---------- starts and their parallel evaluation ---------- */
typedef struct { int **l1, **l2; int cnt, cap; } Starts;
static void starts_add(Starts *s, int *l1, int *l2) {
	if (s->cnt == s->cap) { s->cap = s->cap ? 2 * s->cap : 16; s->l1 = realloc(s->l1, sizeof(int *) * s->cap); s->l2 = realloc(s->l2, sizeof(int *) * s->cap); }
	s->l1[s->cnt] = l1; s->l2[s->cnt] = l2; s->cnt++;
}
static void starts_free(Starts *s) {
	int i;
	for (i = 0; i < s->cnt; i++) { free(s->l1[i]); free(s->l2[i]); }
	free(s->l1); free(s->l2); s->l1 = s->l2 = NULL; s->cnt = s->cap = 0;
}
static int *dupl(const int *a, int n) { int *b = malloc(sizeof(int) * n); memcpy(b, a, sizeof(int) * n); return b; }

typedef struct { const Net *g; const Starts *st; Fit *res; int Q1, Q2, T, tid; } Job;
static SBM_THREAD_FN job_main(void *arg) {
	Job *j = (Job *)arg; Work w; int t;
	long mx = j->g->n > j->g->n2 ? j->g->n : j->g->n2, mq = j->Q1 > j->Q2 ? j->Q1 : j->Q2;
	work_alloc(&w, mx * mq + 1, j->g);
	for (t = j->tid; t < j->st->cnt; t += j->T) {
		Fit *c = &j->res[t];
		fit_alloc(c, j->g, j->Q1, j->Q2);
		tau_from_labels(c->tau, j->g->n, j->Q1, j->st->l1[t]);
		if (j->Q2) tau_from_labels(c->tau2, j->g->n2, j->Q2, j->st->l2[t]);
		if (j->g->p && j->g->model == 0) { c->nnewton = 1; vem_tol(j->g, c, &w, SCREEN_TOL); } else vem(j->g, c, &w);
	}
	work_free(&w);
	SBM_THREAD_RET;
}
/* fits all starts (threads share nothing but the read-only network); the results are in the order of the starts */
static Fit *eval_starts(const Net *g, const Starts *st, int Q1, int Q2, int nthr) {
	Fit *res = calloc(st->cnt > 0 ? st->cnt : 1, sizeof(Fit));
	int T = nthr > 0 ? nthr : ncores(), t, started = 0;
	Job *jobs; sbm_thread *th;
	if (T > st->cnt) T = st->cnt;
	if (T < 1) T = 1;
	jobs = malloc(sizeof(Job) * T); th = malloc(sizeof(sbm_thread) * T);
	for (t = 0; t < T; t++) { jobs[t].g = g; jobs[t].st = st; jobs[t].res = res; jobs[t].Q1 = Q1; jobs[t].Q2 = Q2; jobs[t].T = T; jobs[t].tid = t; }
	for (t = 1; t < T; t++) { if (thread_start(&th[t], job_main, &jobs[t]) == 0) started |= (1 << (t < 30 ? t : 30)); else { job_main(&jobs[t]); th[t] = 0; } }
	job_main(&jobs[0]);
	for (t = 1; t < T; t++) if (started & (1 << (t < 30 ? t : 30))) thread_join(th[t]);
	free(jobs); free(th);
	return res;
}
/* keeps the start with the largest ICL (earliest wins ties); consumes res */
static void reduce_best(Fit *res, int cnt, Fit *best, int *have) {
	int t;
	for (t = 0; t < cnt; t++) {
		if (!*have || res[t].icl > best->icl) { if (*have) fit_free(best); *best = res[t]; *have = 1; }
		else fit_free(&res[t]);
	}
	free(res);
}

/* ---------- one-mode search ---------- */
static void sbm_splits(const Net *g, const Embed *em, int kk, const Fit *src, int Q, Starts *st, Rng *rng) {
	int n = g->n, i, q, d = Q < kk ? Q : kk;
	int *idx = malloc(sizeof(int) * n), *hard = malloc(sizeof(int) * n), *l2 = malloc(sizeof(int) * n);
	hard_labels(src->tau, n, src->Q, hard);
	for (q = 0; q < Q - 1; q++) {
		int m = 0, *hs;
		for (i = 0; i < n; i++) if (hard[i] == q) idx[m++] = i;
		if (m < 2) continue;
		hs = dupl(hard, n);
		kmeans(em->V, kk, d, idx, m, 2, 0, l2, rng);
		for (i = 0; i < m; i++) if (l2[i] == 1) hs[idx[i]] = Q - 1;
		starts_add(st, hs, NULL);
	}
	free(idx); free(hard); free(l2);
}

/* the maxp most similar pairs of rows (side 0) or columns (side 1) of the parameter matrix pi (R x C) */
static int similar_pairs(const double *pi, int R, int C, int side, int maxp, int *pa, int *pb, int symmetric) {
	int a, b, k, np = 0, cnt = side ? C : R;
	double *pd = malloc(sizeof(double) * cnt * cnt);
	for (a = 0; a < cnt; a++) for (b = a + 1; b < cnt; b++) {
		double d = 0;
		if (!side) { for (k = 0; k < C; k++) { double x = pi[a * C + k] - pi[b * C + k]; d += x * x; } }
		else { for (k = 0; k < R; k++) { double x = pi[k * C + a] - pi[k * C + b]; d += x * x; } }
		if (symmetric) { for (k = 0; k < R; k++) { double y = pi[k * C + a] - pi[k * C + b]; d += y * y; } }
		pa[np] = a; pb[np] = b; pd[np] = d; np++;
	}
	if (np > maxp) {
		int t, u;
		for (t = 0; t < maxp; t++) {
			int mi = t;
			for (u = t + 1; u < np; u++) if (pd[u] < pd[mi]) mi = u;
			{ double dd = pd[t]; int ta = pa[t], tb = pb[t]; pd[t] = pd[mi]; pa[t] = pa[mi]; pb[t] = pb[mi]; pd[mi] = dd; pa[mi] = ta; pb[mi] = tb; }
		}
		np = maxp;
	}
	free(pd);
	return np;
}
static void merge_labels(const int *hard, int n, int a, int b, int *lab) {
	int i;
	for (i = 0; i < n; i++) { int l = hard[i]; if (l == b) l = a; if (l > b) l--; lab[i] = l; }
}

static void sbm_merges(const Net *g, const Fit *src, Starts *st) {
	int n = g->n, Q1 = src->Q, k, np;
	int *hard = malloc(sizeof(int) * n), *pa = malloc(sizeof(int) * Q1 * Q1), *pb = malloc(sizeof(int) * Q1 * Q1);
	hard_labels(src->tau, n, Q1, hard);
	np = similar_pairs(src->pi, Q1, Q1, 0, 20, pa, pb, 1);
	for (k = 0; k < np; k++) { int *lab = malloc(sizeof(int) * n); merge_labels(hard, n, pa[k], pb[k], lab); starts_add(st, lab, NULL); }
	free(hard); free(pa); free(pb);
}

/* Exact ICL of the hard partition (conjugate priors, as R greed 0.6.2): Bernoulli stochastic block model (Beta(1,1) on the connectivities) and
   Poisson degree-corrected block model (uniform degree parameters in each block, rate prior as greed's default p); Dirichlet(1) on the proportions.
   Returns NaN for the other models. */
/* exact ICL of the degree-corrected latent block model (R greed DcLbm; Dirichlet(1) on both sets of proportions, p = mean tie value; greed leaves the
   terms that do not depend on the partition out: sum log (d_i)!, sum log y!, links log p, and so does this function), at the hard partition */
static double exact_icl_dclbm(const Net *g, const Fit *f) {
	int n1 = g->n, n2 = g->n2, Q1 = f->Q, Q2 = f->Q2, i, k, l, K1 = 0, K2 = 0, *l1, *l2, *m1, *m2;
	long e;
	double *X, *c1, *c2, *dr, *dcl, s = 0, nl = 0, p;
	l1 = malloc(sizeof(int) * n1); l2 = malloc(sizeof(int) * n2); m1 = malloc(sizeof(int) * Q1); m2 = malloc(sizeof(int) * Q2);
	hard_labels(f->tau, n1, Q1, l1); hard_labels(f->tau2, n2, Q2, l2);
	c1 = calloc(Q1, sizeof(double)); c2 = calloc(Q2, sizeof(double));
	for (i = 0; i < n1; i++) c1[l1[i]] += 1;
	for (i = 0; i < n2; i++) c2[l2[i]] += 1;
	for (k = 0; k < Q1; k++) m1[k] = c1[k] > 0 ? K1++ : -1;
	for (l = 0; l < Q2; l++) m2[l] = c2[l] > 0 ? K2++ : -1;
	free(c1); free(c2);
	X = calloc((size_t)K1 * K2, sizeof(double)); c1 = calloc(K1, sizeof(double)); c2 = calloc(K2, sizeof(double)); dr = calloc(K1, sizeof(double)); dcl = calloc(K2, sizeof(double));
	for (i = 0; i < n1; i++) c1[m1[l1[i]]] += 1;
	for (i = 0; i < n2; i++) c2[m2[l2[i]]] += 1;
	for (i = 0; i < n1; i++) for (e = g->ao.p[i]; e < g->ao.p[i + 1]; e++) { X[(size_t)m1[l1[i]] * K2 + m2[l2[g->ao.j[e]]]] += g->ao.v[e]; nl += g->ao.v[e]; }
	for (k = 0; k < K1; k++) for (l = 0; l < K2; l++) { dr[k] += X[(size_t)k * K2 + l]; dcl[l] += X[(size_t)k * K2 + l]; }
	p = nl / ((double)n1 * n2);
	for (k = 0; k < K1; k++) s += lgamma(c1[k]) - lgamma(c1[k] + dr[k]) + dr[k] * log(c1[k]);
	for (l = 0; l < K2; l++) s += lgamma(c2[l]) - lgamma(c2[l] + dcl[l]) + dcl[l] * log(c2[l]);
	for (k = 0; k < K1; k++) for (l = 0; l < K2; l++) { double x = X[(size_t)k * K2 + l]; s += lgamma(x + 1) - (x + 1) * log(p * c1[k] * c2[l] + 1); }
	s += lgamma((double)K1) - lgamma((double)n1 + K1) + lgamma((double)K2) - lgamma((double)n2 + K2);
	for (k = 0; k < K1; k++) s += lgamma(c1[k] + 1);
	for (l = 0; l < K2; l++) s += lgamma(c2[l] + 1);
	free(l1); free(l2); free(m1); free(m2); free(X); free(c1); free(c2); free(dr); free(dcl);
	return s;
}

static double exact_icl(const Net *g, const Fit *f) {
	int n = g->n, Q = f->Q, i, k, l, K = 0, *lab, *map;
	long e;
	double *X, *cnt, *din, *dout, s = 0, icl = 0;
	if (g->kind) return (g->dc && g->model == 1 && !g->nmiss && !g->p) ? exact_icl_dclbm(g, f) : NAN;
	if (g->nmiss || g->p || !((g->model == 0 && !g->dc) || (g->model == 1 && g->dc))) return NAN;
	lab = malloc(sizeof(int) * n); map = malloc(sizeof(int) * Q);
	hard_labels(f->tau, n, Q, lab);
	cnt = calloc(Q, sizeof(double));
	for (i = 0; i < n; i++) cnt[lab[i]] += 1;
	for (k = 0; k < Q; k++) map[k] = cnt[k] > 0 ? K++ : -1;
	free(cnt);
	X = calloc((size_t)K * K, sizeof(double)); cnt = calloc(K, sizeof(double)); din = calloc(K, sizeof(double)); dout = calloc(K, sizeof(double));
	for (i = 0; i < n; i++) cnt[map[lab[i]]] += 1;
	for (i = 0; i < n; i++) for (e = g->ao.p[i]; e < g->ao.p[i + 1]; e++) X[(size_t)map[lab[i]] * K + map[lab[g->ao.j[e]]]] += g->ao.v[e];
	for (k = 0; k < K; k++) for (l = 0; l < K; l++) { dout[k] += X[(size_t)k * K + l]; din[l] += X[(size_t)k * K + l]; }
	if (g->model == 0) {
		#define LB(x, cc) (lgamma(1 + (x)) + lgamma(1 + (cc) - (x)) + lgamma(2.0) - lgamma(1.0) - lgamma(1.0) - lgamma(2 + (cc)))
		for (k = 0; k < K; k++) for (l = 0; l < K; l++) {
			double x = X[(size_t)k * K + l];
			if (g->directed) s += LB(x, cnt[k] * cnt[l]);
			else if (k < l) s += LB(x, cnt[k] * cnt[l]);
			else if (k == l && cnt[k] > 1) s += LB(x / 2, cnt[k] * (cnt[k] - 1) / 2);
		}
		#undef LB
	} else {
		double nl = 0, cden = 0, p, cst = 0;
		for (i = 0; i < n; i++) for (e = g->ao.p[i]; e < g->ao.p[i + 1]; e++) { nl += g->ao.v[e]; cden += lgamma(g->ao.v[e] + 1); }
		p = nl / ((double)n * n);
		for (k = 0; k < K; k++) s += lgamma(cnt[k]) - lgamma(cnt[k] + din[k]) + din[k] * log(cnt[k]);
		if (g->directed) {
			for (k = 0; k < K; k++) s += lgamma(cnt[k]) - lgamma(cnt[k] + dout[k]) + dout[k] * log(cnt[k]);
			for (k = 0; k < K; k++) for (l = 0; l < K; l++) { double x = X[(size_t)k * K + l]; s += lgamma(x + 1) - (x + 1) * log(p * cnt[k] * cnt[l] + 1); }
			for (i = 0; i < n; i++) cst += lgamma(g->wo[i] + 1) + lgamma(g->wi[i] + 1);
			cst += -cden + nl * log(p);
		} else {
			for (k = 0; k < K; k++) for (l = 0; l <= k; l++) {
				double x = X[(size_t)k * K + l];
				if (k == l) { if (cnt[k] > 1) s += lgamma(x / 2 + 1) - (x / 2 + 1) * log(p * cnt[k] * (cnt[k] - 1) / 2 + 1); }
				else s += lgamma(x + 1) - (x + 1) * log(p * cnt[k] * cnt[l] + 1);
			}
			for (i = 0; i < n; i++) cst += lgamma(g->wo[i] + 1);
			cst += -cden / 2 + (nl / 2) * log(p);
		}
		s += cst;
	}
	icl = s + lgamma((double)K);
	for (k = 0; k < K; k++) icl += lgamma(cnt[k] + 1);
	icl -= lgamma((double)n + K);
	free(lab); free(map); free(X); free(cnt); free(din); free(dout);
	return icl;
}

/* polishes a screened fit: full Newton M steps and the tight EM tolerance, starting from its own memberships */
static void polish_fit(const Net *g, Fit *f, Work *w) {
	if (!(g->p && g->model == 0)) return;
	f->nnewton = 40;
	vem(g, f, w);
}

static void write_fit(const Net *g, const Fit *b, double *out, long *off) {
	long o = *off;
	int n = g->n, Q1 = b->Q, Q2 = b->Q2, c2 = Q2 ? Q2 : Q1;
	out[o++] = b->icl; out[o++] = b->J; out[o++] = b->PL; out[o++] = b->H; out[o++] = b->npar; out[o++] = b->sigma2; out[o++] = Q1; out[o++] = Q2;
	{ double ex = exact_icl(g, b); out[o++] = (ex == ex) ? ex : -1e300; }
	memcpy(out + o, b->tau, sizeof(double) * (long)n * Q1); o += (long)n * Q1;
	if (Q2) { memcpy(out + o, b->tau2, sizeof(double) * (long)g->n2 * Q2); o += (long)g->n2 * Q2; }
	memcpy(out + o, b->pi, sizeof(double) * Q1 * c2); o += (long)Q1 * c2;
	memcpy(out + o, b->alpha, sizeof(double) * Q1); o += Q1;
	if (Q2) { memcpy(out + o, b->alpha2, sizeof(double) * Q2); o += Q2; }
	if (g->p) { memcpy(out + o, b->beta, sizeof(double) * g->p); o += g->p; memcpy(out + o, b->info, sizeof(double) * g->p * g->p); o += (long)g->p * g->p; }
	*off = o;
}
static long fit_len(const Net *g, int Q1, int Q2) {
	long l = 9 + (long)g->n * Q1 + (long)Q1 * (Q2 ? Q2 : Q1) + Q1;
	if (Q2) l += (long)g->n2 * Q2 + Q2;
	return l + g->p + (long)g->p * g->p;
}

static int run_sbm(Net *g, int Qmin, int Qout, int restarts, int nthr, int fast, Rng *rng, double *out) {
	int Qmax = Qout < g->n ? Qout + 1 : Qout;
	int n = g->n, Q, r, i;
	Fit *best; int *have; Embed em; double p1;
	long off = 0;
	int *lab = malloc(sizeof(int) * n);
	int kk = (Qmax + 3 < n) ? Qmax + 3 : n;
	best = calloc(Qmax + 2, sizeof(Fit)); have = calloc(Qmax + 2, sizeof(int));
	if (!lab || !best || !have) return 909;
	p1 = g->ybar;
	if (g->p && !cov_residual(g)) return 909;
	if (!embed_compute(g, p1, kk, &em, rng)) return 909;
	for (Q = 1; Q <= Qmax; Q++) {
		Starts st; Fit *res; memset(&st, 0, sizeof(st));
		if (Q == 1) { int *l = calloc(n, sizeof(int)); starts_add(&st, l, NULL); }
		else {
			int d = Q < kk ? Q : kk;
			kmeans(em.V, kk, d, NULL, n, Q, 0, lab, rng); starts_add(&st, dupl(lab, n), NULL);
			sbm_splits(g, &em, kk, &best[Q - 1], Q, &st, rng);
			for (r = 0; r < restarts; r++) {
				if (r % 2 == 0) kmeans(em.V, kk, d, NULL, n, Q, 1, lab, rng);
				else for (i = 0; i < n; i++) lab[i] = rng_int(rng, Q);
				starts_add(&st, dupl(lab, n), NULL);
			}
		}
		res = eval_starts(g, &st, Q, 0, nthr);
		reduce_best(res, st.cnt, &best[Q], &have[Q]);
		starts_free(&st);
	}
	if (!fast) {
		for (Q = Qmax - 1; Q >= 2; Q--) {
			Starts st; Fit *res; memset(&st, 0, sizeof(st));
			sbm_merges(g, &best[Q + 1], &st);
			res = eval_starts(g, &st, Q, 0, nthr); reduce_best(res, st.cnt, &best[Q], &have[Q]); starts_free(&st);
		}
		for (Q = 2; Q <= Qmax; Q++) {
			Starts st; Fit *res; memset(&st, 0, sizeof(st));
			sbm_splits(g, &em, kk, &best[Q - 1], Q, &st, rng);
			res = eval_starts(g, &st, Q, 0, nthr); reduce_best(res, st.cnt, &best[Q], &have[Q]); starts_free(&st);
		}
	}
	if (g->p && g->model == 0) {
		Work pw; long mx = g->n > g->n2 ? g->n : g->n2;
		if (work_alloc(&pw, mx * Qmax + 1, g)) { for (Q = Qmin; Q <= Qout; Q++) polish_fit(g, &best[Q], &pw); work_free(&pw); }
	}
	for (Q = Qmin; Q <= Qout; Q++) write_fit(g, &best[Q], out, &off);
	for (Q = 1; Q <= Qmax; Q++) if (have[Q]) fit_free(&best[Q]);
	free(best); free(have); free(em.V); free(lab);
	return 0;
}

/* ---------- two-mode search over the grid of (Q1, Q2) ---------- */
static int run_lbm(Net *g, int Q1min, int Q1max, int Q2min, int Q2max, int restarts, int nthr, int fast, Rng *rng, double *out) {
	int n1 = g->n, n2 = g->n2, a, b, r, i, q;
	int kk = (((Q1max > Q2max ? Q1max : Q2max) + 3) < (n1 < n2 ? n1 : n2)) ? ((Q1max > Q2max ? Q1max : Q2max) + 3) : (n1 < n2 ? n1 : n2);
	Fit *best = calloc((long)Q1max * Q2max, sizeof(Fit)); int *have = calloc((long)Q1max * Q2max, sizeof(int));
	Embed em; long off = 0;
	int *lab1 = malloc(sizeof(int) * n1), *lab2 = malloc(sizeof(int) * n2);
	int *idx = malloc(sizeof(int) * (n1 > n2 ? n1 : n2)), *l2 = malloc(sizeof(int) * (n1 > n2 ? n1 : n2));
	int pass;
	#define CELL(x, y) ((long)((x) - 1) * Q2max + ((y) - 1))
	if (g->p && !cov_residual(g)) return 909;
	if (!embed_compute(g, g->ybar, kk, &em, rng)) return 909;
	for (pass = 0; pass < (fast ? 1 : 3); pass++) {
		int fwd = (pass != 1);
		for (a = fwd ? 1 : Q1max; fwd ? a <= Q1max : a >= 1; a += fwd ? 1 : -1)
		for (b = fwd ? 1 : Q2max; fwd ? b <= Q2max : b >= 1; b += fwd ? 1 : -1) {
			Starts st; Fit *res; long c = CELL(a, b);
			memset(&st, 0, sizeof(st));
			if (fwd) {
				if (pass == 0) {
					if (a == 1 && b == 1) starts_add(&st, calloc(n1, sizeof(int)), calloc(n2, sizeof(int)));
					else {
						kmeans(em.V, kk, a < kk ? a : kk, NULL, n1, a, 0, lab1, rng); kmeans(em.C, kk, b < kk ? b : kk, NULL, n2, b, 0, lab2, rng);
						starts_add(&st, dupl(lab1, n1), dupl(lab2, n2));
					}
				}
				if (a > 1 && have[CELL(a - 1, b)]) { /* split a row block of the (a-1, b) solution */
					const Fit *s = &best[CELL(a - 1, b)]; int *hr = malloc(sizeof(int) * n1), *hc = malloc(sizeof(int) * n2);
					hard_labels(s->tau, n1, s->Q, hr); hard_labels(s->tau2, n2, s->Q2, hc);
					for (q = 0; q < a - 1; q++) {
						int m = 0, *hs;
						for (i = 0; i < n1; i++) if (hr[i] == q) idx[m++] = i;
						if (m < 2) continue;
						hs = dupl(hr, n1);
						kmeans(em.V, kk, a < kk ? a : kk, idx, m, 2, 0, l2, rng);
						for (i = 0; i < m; i++) if (l2[i] == 1) hs[idx[i]] = a - 1;
						starts_add(&st, hs, dupl(hc, n2));
					}
					free(hr); free(hc);
				}
				if (b > 1 && have[CELL(a, b - 1)]) {
					const Fit *s = &best[CELL(a, b - 1)]; int *hr = malloc(sizeof(int) * n1), *hc = malloc(sizeof(int) * n2);
					hard_labels(s->tau, n1, s->Q, hr); hard_labels(s->tau2, n2, s->Q2, hc);
					for (q = 0; q < b - 1; q++) {
						int m = 0, *hs;
						for (i = 0; i < n2; i++) if (hc[i] == q) idx[m++] = i;
						if (m < 2) continue;
						hs = dupl(hc, n2);
						kmeans(em.C, kk, b < kk ? b : kk, idx, m, 2, 0, l2, rng);
						for (i = 0; i < m; i++) if (l2[i] == 1) hs[idx[i]] = b - 1;
						starts_add(&st, dupl(hr, n1), hs);
					}
					free(hr); free(hc);
				}
				if (pass == 0 && !(a == 1 && b == 1)) for (r = 0; r < restarts; r++) {
					if (r % 2 == 0) { kmeans(em.V, kk, a < kk ? a : kk, NULL, n1, a, 1, lab1, rng); kmeans(em.C, kk, b < kk ? b : kk, NULL, n2, b, 1, lab2, rng); }
					else { for (i = 0; i < n1; i++) lab1[i] = rng_int(rng, a); for (i = 0; i < n2; i++) lab2[i] = rng_int(rng, b); }
					starts_add(&st, dupl(lab1, n1), dupl(lab2, n2));
				}
			} else { /* backward: merges of the (a+1, b) and (a, b+1) solutions */
				int np, k; int *pa = malloc(sizeof(int) * 5000), *pb = malloc(sizeof(int) * 5000);
				if (a < Q1max && have[CELL(a + 1, b)]) {
					const Fit *s = &best[CELL(a + 1, b)]; int *hr = malloc(sizeof(int) * n1), *hc = malloc(sizeof(int) * n2);
					hard_labels(s->tau, n1, s->Q, hr); hard_labels(s->tau2, n2, s->Q2, hc);
					np = similar_pairs(s->pi, s->Q, s->Q2, 0, 20, pa, pb, 0);
					for (k = 0; k < np; k++) { int *lr = malloc(sizeof(int) * n1); merge_labels(hr, n1, pa[k], pb[k], lr); starts_add(&st, lr, dupl(hc, n2)); }
					free(hr); free(hc);
				}
				if (b < Q2max && have[CELL(a, b + 1)]) {
					const Fit *s = &best[CELL(a, b + 1)]; int *hr = malloc(sizeof(int) * n1), *hc = malloc(sizeof(int) * n2);
					hard_labels(s->tau, n1, s->Q, hr); hard_labels(s->tau2, n2, s->Q2, hc);
					np = similar_pairs(s->pi, s->Q, s->Q2, 1, 20, pa, pb, 0);
					for (k = 0; k < np; k++) { int *lc = malloc(sizeof(int) * n2); merge_labels(hc, n2, pa[k], pb[k], lc); starts_add(&st, dupl(hr, n1), lc); }
					free(hr); free(hc);
				}
				free(pa); free(pb);
			}
			if (st.cnt > 0) { res = eval_starts(g, &st, a, b, nthr); reduce_best(res, st.cnt, &best[c], &have[c]); }
			starts_free(&st);
		}
	}
	if (g->p && g->model == 0) {
		Work pw; long mx = g->n > g->n2 ? g->n : g->n2, mq = Q1max > Q2max ? Q1max : Q2max;
		if (work_alloc(&pw, mx * mq + 1, g)) { for (a = Q1min; a <= Q1max; a++) for (b = Q2min; b <= Q2max; b++) polish_fit(g, &best[CELL(a, b)], &pw); work_free(&pw); }
	}
	for (a = Q1min; a <= Q1max; a++) for (b = Q2min; b <= Q2max; b++) write_fit(g, &best[CELL(a, b)], out, &off);
	for (i = 0; i < Q1max * Q2max; i++) if (have[i]) fit_free(&best[i]);
	free(best); free(have); free(em.V); free(em.C); free(lab1); free(lab2); free(idx); free(l2);
	#undef CELL
	return 0;
}

/* ---------- plugin entry ---------- */
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
		ST_double v; int n, n2, directed, model, Qmin, Qmax, Q2min, Q2max, restarts, nthr, kind, dc, fast, ncov; unsigned long long seed; double s2fix;
		long i, j, ma = 0, mm = 0, cap, k, outn = 0, tot;
		int *ar, *ac, *mr, *mc; double *av, *out; Net g; Rng rng; int ncol;
		SF_vdata(1, 1, &v); if (v != SBM_PROTOCOL) { SF_error("sbm_vem: the model block is invalid (plugin and package versions do not match?)\n"); return 198; }
		SF_vdata(1, 2, &v); n = (int)v;
		SF_vdata(1, 3, &v); directed = (int)v;
		SF_vdata(1, 4, &v); model = (int)v;
		SF_vdata(1, 5, &v); Qmin = (int)v;
		SF_vdata(1, 6, &v); Qmax = (int)v;
		SF_vdata(1, 7, &v); seed = (unsigned long long)v;
		SF_vdata(1, 8, &v); restarts = (int)v;
		SF_vdata(1, 9, &v); s2fix = (double)v;
		SF_vdata(1, 10, &v); nthr = (int)v;
		SF_vdata(1, 11, &v); kind = (int)v;
		SF_vdata(1, 12, &v); n2 = (int)v;
		SF_vdata(1, 13, &v); Q2min = (int)v;
		SF_vdata(1, 14, &v); Q2max = (int)v;
		SF_vdata(1, 15, &v); dc = (int)v;
		SF_vdata(1, 16, &v); fast = (int)v;
		SF_vdata(1, 17, &v); ncov = (int)v;
		if (!kind) { n2 = n; Q2min = Q2max = 0; }
		if (model < 0 || model > 2 || n < 2 || Qmin < 1 || Qmax < Qmin || Qmax > n || ncov < 0 || ncov > MAXCOV || L != HDR + (long)(1 + ncov) * n * n2) { SF_error("sbm_vem: the model block is invalid\n"); return 198; }
		if (kind && (n2 < 2 || Q2min < 1 || Q2max < Q2min || Q2max > n2)) { SF_error("sbm_vem: the model block is invalid\n"); return 198; }
		if (dc && ncov > 0) { SF_error("sbm_vem: degree correction cannot be combined with covariates\n"); return 198; }
		if (dc && model != 1) { SF_error("sbm_vem: degree correction needs the Poisson model\n"); return 198; }
		if (!sp_init()) return 909;
		rng_seed(&rng, seed);
		cap = 1024;
		ar = malloc(sizeof(int) * cap); ac = malloc(sizeof(int) * cap); mr = malloc(sizeof(int) * cap); mc = malloc(sizeof(int) * cap); av = malloc(sizeof(double) * cap);
		memset(&g, 0, sizeof(g)); g.kind = kind; g.n = n; g.n2 = n2; g.directed = kind ? 0 : directed; g.model = model; g.dc = dc; g.sigma2fix = s2fix; g.p = ncov;
		ncol = n2;
		if (ncov > 0) {
			long NN = (long)n * n2;
			g.Y = calloc(NN, sizeof(double)); g.Mf = calloc(NN, sizeof(double)); g.X = calloc(ncov, sizeof(double *));
			g.sxy = calloc(ncov, sizeof(double)); g.sxx = calloc((size_t)ncov * ncov, sizeof(double));
			if (!g.Y || !g.Mf || !g.X || !g.sxy || !g.sxx) { SF_error("sbm_vem: not enough memory for the covariate matrices\n"); return 909; }
			for (k = 0; k < ncov; k++) { g.X[k] = calloc(NN, sizeof(double)); if (!g.X[k]) { SF_error("sbm_vem: not enough memory for the covariate matrices\n"); return 909; } }
		}
		{
			long obs = 0; double sy = 0;
			int und = (!kind && !directed);
			for (i = 0; i < n; i++) for (j = 0; j < ncol; j++) {
				if (!kind && i == j) continue;
				if (und && j < i) continue;
				SF_vdata(1, (ST_int)(HDR + i * ncol + j + 1), &v);
				if (ma + 2 >= cap || mm + 2 >= cap) { cap *= 2; ar = realloc(ar, sizeof(int) * cap); ac = realloc(ac, sizeof(int) * cap); mr = realloc(mr, sizeof(int) * cap); mc = realloc(mc, sizeof(int) * cap); av = realloc(av, sizeof(double) * cap); }
				if (!ar || !ac || !mr || !mc || !av) return 909;
				if (SF_is_missing(v)) {
					mr[mm] = (int)i; mc[mm] = (int)j; mm++;
					if (und) { mr[mm] = (int)j; mc[mm] = (int)i; mm++; }
				} else {
					double x = (model == 0) ? (v > 0 ? 1.0 : 0.0) : (double)v;
					obs++;
					if (ncov > 0) { g.Y[i * ncol + j] = x; g.Mf[i * ncol + j] = 1; if (und) { g.Y[j * ncol + i] = x; g.Mf[j * ncol + i] = 1; } }
					sy += x * (und ? 2 : 1);
					g.sumy2 += x * x * (und ? 2 : 1);
					if (model == 1) g.cst -= lgamma(x + 1);
					if (x != 0) {
						ar[ma] = (int)i; ac[ma] = (int)j; av[ma] = x; ma++;
						if (und) { ar[ma] = (int)j; ac[ma] = (int)i; av[ma] = x; ma++; }
					}
				}
			}
			if (ncov > 0) {
				long t; int a, b;
				for (a = 0; a < ncov; a++) {
					for (i = 0; i < n; i++) for (j = 0; j < ncol; j++) {
						t = i * ncol + j;
						if (!g.Mf[t]) continue;
						SF_vdata(1, (ST_int)(HDR + (long)(1 + a) * n * ncol + t + 1), &v);
						if (SF_is_missing(v)) { SF_error("sbm_vem: a covariate has missing values\n"); return 198; }
						g.X[a][t] = (double)v;
					}
				}
				for (a = 0; a < ncov; a++) {
					for (t = 0; t < (long)n * ncol; t++) g.sxy[a] += g.X[a][t] * g.Y[t];
					for (b = 0; b < ncov; b++) for (t = 0; t < (long)n * ncol; t++) g.sxx[a * ncov + b] += g.X[a][t] * g.X[b][t];
				}
			}
			g.ndata = (double)obs; g.nmiss = (int)mm; g.nobs = obs * (und ? 2.0 : 1.0);
			g.ybar = g.nobs > 0 ? sy / g.nobs : 0;
		}
		if (g.ndata < 1) { SF_error("sbm_vem: no observed ties\n"); return 198; }
		csr_build(&g.ao, n, ar, ac, av, ma); csr_build(&g.mo, n, mr, mc, NULL, mm);
		if (kind) { csr_build(&g.ai, n2, ac, ar, av, ma); csr_build(&g.mi, n2, mc, mr, NULL, mm); }
		else if (directed) { csr_build(&g.ai, n, ac, ar, av, ma); csr_build(&g.mi, n, mc, mr, NULL, mm); }
		else { g.ai = g.ao; g.mi = g.mo; }
		if (dc) {
			long e;
			int sep = (directed || kind);
			g.wo = calloc(n, sizeof(double)); g.wi = sep ? calloc(kind ? n2 : n, sizeof(double)) : g.wo;
			for (i = 0; i < n; i++) for (e = g.ao.p[i]; e < g.ao.p[i + 1]; e++) { g.wo[i] += g.ao.v[e]; if (sep) g.wi[g.ao.j[e]] += g.ao.v[e]; }
			for (i = 0; i < n; i++) g.sumw += g.wo[i];
			if (g.sumw <= 0) { SF_error("sbm_vem: the network has no ties\n"); return 198; }
			for (i = 0; i < n; i++) for (e = g.ao.p[i]; e < g.ao.p[i + 1]; e++) g.cst += ((directed || kind) ? 1.0 : 0.5) * g.ao.v[e] * log(g.wo[i] * g.wi[g.ao.j[e]]);
		}
		tot = 0;
		if (kind) { for (k = Qmin; k <= Qmax; k++) for (j = Q2min; j <= Q2max; j++) tot += fit_len(&g, (int)k, (int)j); }
		else for (k = Qmin; k <= Qmax; k++) tot += fit_len(&g, (int)k, 0);
		out = malloc(sizeof(double) * tot);
		if (!out) return 909;
		if (fast && restarts > 2) restarts = 2;
		rc = kind ? run_lbm(&g, Qmin, Qmax, Q2min, Q2max, restarts, nthr, fast, &rng, out) : run_sbm(&g, Qmin, Qmax, restarts, nthr, fast, &rng, out);
		if (!rc) for (outn = 0; outn < tot; outn++) SF_vstore(2, (ST_int)(1 + outn), out[outn]);
		free(out); free(ar); free(ac); free(mr); free(mc); free(av);
		csr_free(&g.ao); csr_free(&g.mo);
		if (kind || directed) { csr_free(&g.ai); csr_free(&g.mi); }
		if (dc) { if (directed || kind) free(g.wi); free(g.wo); }
		if (ncov > 0) { int a; free(g.Y); free(g.Mf); for (a = 0; a < ncov; a++) free(g.X[a]); free(g.X); free(g.sxy); free(g.sxx); free(g.Rres); }
	}
	return rc;
}
