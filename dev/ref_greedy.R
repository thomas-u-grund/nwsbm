# Fixtures with weak block structure (the VEM partition and the partition of greed's greedy search can differ in a few nodes) for the greedy option of nwsbm.
# Run: Rscript dev/ref_greedy.R dev/ref   (writes the same files as dev/ref_dc.R: <name>_edges.csv, _missing.csv, _truth.csv, _fit.csv, _membership.csv)
suppressMessages(library(greed))
out <- commandArgs(TRUE)[1]; if (is.na(out)) out <- "dev/ref"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
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
set.seed(61); n <- 200; z <- sample(rep(1:3, length.out = n)); P <- matrix(0.07, 3, 3); diag(P) <- c(0.24, 0.20, 0.18)
A <- matrix(rbinom(n * n, 1, P[z, z]), n, n); A[lower.tri(A)] <- t(A)[lower.tri(A)]; diag(A) <- 0
write_all("gw_sbm_u_n200_Q3", A, z, FALSE, Sbm(type = "undirected"))
set.seed(62); n <- 150; z <- sample(rep(1:3, length.out = n)); th <- rlnorm(n, 0, 0.7); P <- matrix(0.04, 3, 3); diag(P) <- 0.16
M <- outer(th, th) * P[z, z] * 3
A <- matrix(rpois(n * n, M), n, n); A[lower.tri(A)] <- t(A)[lower.tri(A)]; diag(A) <- 0
write_all("gw_dc_u_n150_Q3", A, z, FALSE, DcSbm(type = "undirected"))
set.seed(63); n <- 120; z <- sample(rep(1:3, length.out = n)); P <- matrix(0.08, 3, 3); diag(P) <- c(0.30, 0.25, 0.22)
A <- matrix(rbinom(n * n, 1, P[z, z]), n, n); diag(A) <- 0
write_all("gw_sbm_d_n120_Q3", A, z, TRUE, Sbm(type = "directed"))
