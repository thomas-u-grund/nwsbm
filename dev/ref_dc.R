# Reference fits with R greed 0.6.2 (DcSbm and Sbm, exact ICL) for the degree-corrected tests and the exact-ICL option of nwsbm.
# Run: Rscript dev/ref_dc.R dev/ref   (writes <name>_edges.csv, _missing.csv, _truth.csv, _fit.csv (K, ICL), _membership.csv)
suppressMessages(library(greed))
out <- commandArgs(TRUE)[1]; if (is.na(out)) out <- "dev/ref"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
plantdc <- function(n, Q, directed, within, between, sdlog, scale, seed) {
  set.seed(seed); z <- sample(rep(seq_len(Q), length.out = n)); th <- rlnorm(n, 0, sdlog); th2 <- if (directed) rlnorm(n, 0, sdlog) else th
  P <- matrix(between, Q, Q); diag(P) <- within
  M <- outer(th, th2) * P[z, z] * scale
  A <- matrix(rpois(n * n, M), n, n)
  if (!directed) A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  list(A = A, z = z)
}
write_all <- function(name, A, z, directed, model) {
  idx <- which(A != 0, arr.ind = TRUE)
  if (!directed) idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  write.csv(data.frame(from = idx[, 1], to = idx[, 2], value = A[idx]), file.path(out, paste0(name, "_edges.csv")), row.names = FALSE)
  write.csv(data.frame(from = integer(0), to = integer(0)), file.path(out, paste0(name, "_missing.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_along(z), truth = z), file.path(out, paste0(name, "_truth.csv")), row.names = FALSE)
  f <- greed(A, model = model, K = 8, verbose = FALSE)
  write.csv(data.frame(Q = f@K, ICL = f@icl), file.path(out, paste0(name, "_fit.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_along(z), block = clustering(f)), file.path(out, paste0(name, "_membership.csv")), row.names = FALSE)
  cat(name, f@K, f@icl, "\n")
}
d <- plantdc(150, 3, FALSE, 0.12, 0.02, 0.6, 3, 41); write_all("dc_u_n150_Q3", d$A, d$z, FALSE, DcSbm(type = "undirected"))
d <- plantdc(120, 2, TRUE, 0.10, 0.02, 0.5, 3, 42); write_all("dc_d_n120_Q2", d$A, d$z, TRUE, DcSbm(type = "directed"))
d <- plantdc(200, 4, FALSE, 0.10, 0.015, 0.7, 3, 43); write_all("dc_u_n200_Q4", d$A, d$z, FALSE, DcSbm(type = "undirected"))
# binary SBM and Poisson SBM for the exact ICL of greed
set.seed(44); n <- 100; z <- sample(rep(1:3, length.out = n)); P <- matrix(0.05, 3, 3); diag(P) <- c(0.5, 0.4, 0.3)
A <- matrix(rbinom(n * n, 1, P[z, z]), n, n); A[lower.tri(A)] <- t(A)[lower.tri(A)]; diag(A) <- 0
write_all("gr_sbm_u_n100_Q3", A, z, FALSE, Sbm(type = "undirected"))
