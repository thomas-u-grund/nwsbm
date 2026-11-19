# Reference fits with R greed 0.6.2 (DcLbm, exact ICL) for the degree-corrected two-mode tests of nwsbm.
# Run: Rscript dev/ref_dclbm.R dev/ref   (writes <name>_edges.csv (row, col, value), _missing.csv, _fit.csv (K1, K2, ICL), _membership_rows.csv,
# _membership_cols.csv, _truth_rows.csv, _truth_cols.csv; rows of the matrix are mode 2 of the nwset network, see GOTCHA)
suppressMessages(library(greed))
out <- commandArgs(TRUE)[1]; if (is.na(out)) out <- "dev/ref"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
plantdclbm <- function(n1, n2, Q1, Q2, sdlog, scale, seed) {
  set.seed(seed)
  z1 <- sample(rep(seq_len(Q1), length.out = n1)); z2 <- sample(rep(seq_len(Q2), length.out = n2))
  th1 <- rlnorm(n1, 0, sdlog); th2 <- rlnorm(n2, 0, sdlog)
  P <- matrix(0.03, Q1, Q2); for (k in seq_len(min(Q1, Q2))) P[k, k] <- 0.2; P[1, Q2] <- 0.12
  M <- outer(th1, th2) * P[z1, z2] * scale
  list(A = matrix(rpois(n1 * n2, M), n1, n2), z1 = z1, z2 = z2)
}
# exact ICL of the degree-corrected latent block model (Dirichlet(1) on both sets of proportions, p = mean of X), as implemented in nwsbm; like greed it leaves out the terms that do not depend on the partition
# (sum log d_i!, sum log x_ij!, links * log p)
dclbm_icl <- function(A, z1, z2) {
  n1 <- nrow(A); n2 <- ncol(A); K1 <- max(z1); K2 <- max(z2)
  nk <- tabulate(z1, K1); nl <- tabulate(z2, K2)
  X <- matrix(0, K1, K2); for (k in 1:K1) for (l in 1:K2) X[k, l] <- sum(A[z1 == k, z2 == l])
  dr <- rowSums(X); dc <- colSums(X); p <- sum(A) / (n1 * n2)
  s <- sum(lgamma(nk) - lgamma(nk + dr) + dr * log(nk)) + sum(lgamma(nl) - lgamma(nl + dc) + dc * log(nl))
  s <- s + sum(lgamma(X + 1) - (X + 1) * log(p * outer(nk, nl) + 1))
  s + lgamma(K1) + sum(lgamma(nk + 1)) - lgamma(n1 + K1) + lgamma(K2) + sum(lgamma(nl + 1)) - lgamma(n2 + K2)
}
write_all <- function(name, d) {
  A <- d$A; idx <- which(A != 0, arr.ind = TRUE)
  write.csv(data.frame(from = idx[, 1], to = idx[, 2], value = A[idx]), file.path(out, paste0(name, "_edges.csv")), row.names = FALSE)
  write.csv(data.frame(from = integer(0), to = integer(0)), file.path(out, paste0(name, "_missing.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_along(d$z1), truth = d$z1), file.path(out, paste0(name, "_truth_rows.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_along(d$z2), truth = d$z2), file.path(out, paste0(name, "_truth_cols.csv")), row.names = FALSE)
  f <- greed(A, model = DcLbm(), K = 8, verbose = FALSE)
  write.csv(data.frame(K1 = f@Krow, K2 = f@Kcol, ICL = f@icl), file.path(out, paste0(name, "_fit.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_len(nrow(A)), block = f@clrow), file.path(out, paste0(name, "_membership_rows.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_len(ncol(A)), block = f@clcol), file.path(out, paste0(name, "_membership_cols.csv")), row.names = FALSE)
  cat(name, f@Krow, f@Kcol, f@icl, "\n")
  # the exact ICL written out (the formula nwsbm uses), evaluated at greed's partition
  cat("   formula at greed's partition:", dclbm_icl(A, f@clrow, f@clcol), "\n")
}
write_all("dl_40x60", plantdclbm(40, 60, 2, 3, 0.6, 4, 51))
write_all("dl_80x50", plantdclbm(80, 50, 3, 2, 0.5, 3, 52))
write_all("dl_100x100", plantdclbm(100, 100, 3, 3, 0.7, 3, 53))
