# Reference fits with R sbm 1.0.0 (wrapper of blockmodels 1.2.1) for the nwsbm tests.
# Run: Rscript dev/ref_sbm.R dev/ref   (writes <name>_edges.csv, _missing.csv, _fit.csv, _membership.csv, _pi.csv, _icl.csv, and timing.csv)
suppressMessages(library(sbm))
out <- commandArgs(TRUE)[1]; if (is.na(out)) out <- "dev/ref"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

plant <- function(n, Pm, sizes, directed, seed, pmiss = 0) {
  set.seed(seed)
  Q <- nrow(Pm)
  z <- rep(seq_len(Q), times = sizes); z <- sample(z)
  p <- Pm[z, z]
  A <- matrix(rbinom(n * n, 1, p), n, n)
  if (!directed) { A[lower.tri(A)] <- t(A)[lower.tri(A)] }
  diag(A) <- 0
  M <- matrix(FALSE, n, n)
  if (pmiss > 0) {
    M[upper.tri(M)] <- runif(n * (n - 1) / 2) < pmiss
    if (directed) M[lower.tri(M)] <- runif(n * (n - 1) / 2) < pmiss else M <- M | t(M)
    A[M] <- NA
  }
  list(A = A, z = z, M = M)
}
write_net <- function(name, A, z, directed) {
  n <- nrow(A)
  idx <- which(A == 1, arr.ind = TRUE)
  if (!directed) idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  write.csv(data.frame(from = idx[, 1], to = idx[, 2]), file.path(out, paste0(name, "_edges.csv")), row.names = FALSE)
  mi <- which(is.na(A), arr.ind = TRUE)
  if (!directed) mi <- mi[mi[, 1] < mi[, 2], , drop = FALSE]
  write.csv(data.frame(from = mi[, 1], to = mi[, 2]), file.path(out, paste0(name, "_missing.csv")), row.names = FALSE)
  write.csv(data.frame(id = 1:n, truth = z), file.path(out, paste0(name, "_truth.csv")), row.names = FALSE)
}
fit_ref <- function(name, A, directed, qmax) {
  t0 <- proc.time()[3]
  f <- estimateSimpleSBM(A, "bernoulli", directed = directed,
         estimOptions = list(verbosity = 0, plot = FALSE, nbCores = 1, exploreMax = qmax))
  el <- proc.time()[3] - t0
  bm <- f$storedModels
  Q <- f$nbBlocks
  write.csv(data.frame(Q = Q, ICL = f$ICL, n = nrow(A), seconds = el), file.path(out, paste0(name, "_fit.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_len(nrow(A)), block = apply(f$probMemberships, 1, which.max)), file.path(out, paste0(name, "_membership.csv")), row.names = FALSE)
  write.csv(f$connectParam$mean, file.path(out, paste0(name, "_pi.csv")), row.names = FALSE)
  write.csv(f$storedModels, file.path(out, paste0(name, "_icl.csv")), row.names = FALSE)
  el
}
P3 <- matrix(0.05, 3, 3); diag(P3) <- c(0.5, 0.4, 0.3)
P3a <- matrix(c(0.6, 0.1, 0.05, 0.1, 0.5, 0.15, 0.05, 0.15, 0.4), 3)
P4 <- matrix(0.04, 4, 4); diag(P4) <- 0.35; P4[1, 2] <- P4[2, 1] <- 0.15
P2 <- matrix(c(0.4, 0.05, 0.05, 0.3), 2)
Pd3 <- matrix(c(0.5, 0.05, 0.1, 0.2, 0.4, 0.05, 0.05, 0.3, 0.45), 3, byrow = TRUE)
specs <- list(
  list("pp_u_n60_Q2",  60, P2,  c(30, 30), FALSE, 11, 0, 5),
  list("pp_u_n60_Q3",  60, P3a, c(20, 20, 20), FALSE, 12, 0, 6),
  list("pp_u_n150_Q3", 150, P3, c(60, 50, 40), FALSE, 13, 0, 6),
  list("pp_u_n150_Q4", 150, P4, c(50, 40, 35, 25), FALSE, 14, 0, 7),
  list("pp_d_n100_Q3", 100, Pd3, c(40, 30, 30), TRUE, 15, 0, 6),
  list("pp_d_n150_Q2", 150, matrix(c(0.35, 0.05, 0.12, 0.3), 2), c(75, 75), TRUE, 16, 0, 5),
  list("miss_u_n80_Q3", 80, P3a, c(30, 25, 25), FALSE, 17, 0.08, 6),
  list("miss_d_n80_Q2", 80, matrix(c(0.4, 0.06, 0.1, 0.3), 2), c(40, 40), TRUE, 18, 0.06, 5)
)
for (s in specs) {
  d <- plant(s[[2]], s[[3]], s[[4]], s[[5]], s[[6]], s[[7]])
  write_net(s[[1]], d$A, d$z, s[[5]])
  fit_ref(s[[1]], d$A, s[[5]], s[[8]])
  cat(s[[1]], "done\n")
}
# timing networks (planted Q = 3, undirected)
tim <- data.frame()
for (n in c(200, 500, 1000)) {
  sz <- c(round(n * 0.4), round(n * 0.35)); sz <- c(sz, n - sum(sz))
  Pt <- matrix(0.02, 3, 3); diag(Pt) <- c(0.25, 0.2, 0.15)
  if (n >= 500) { Pt <- Pt * 200 / n * 2.5 }
  d <- plant(n, pmin(Pt, 0.9), sz, FALSE, 100 + n, 0)
  nm <- paste0("time_u_n", n)
  write_net(nm, d$A, d$z, FALSE)
  el <- fit_ref(nm, d$A, FALSE, 6)
  tim <- rbind(tim, data.frame(n = n, seconds = el))
  cat(nm, el, "\n")
}
write.csv(tim, file.path(out, "timing.csv"), row.names = FALSE)
