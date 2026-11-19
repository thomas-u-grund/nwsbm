# Reference fits for the valued (Poisson, Gaussian) and two-mode (latent block model) tests of nwsbm: R sbm 1.0.0 (blockmodels 1.2.1).
# Run: Rscript dev/ref_valued.R dev/ref   (writes <name>_edges.csv (from,to,value), _missing.csv, _truth.csv, _fit.csv, _membership.csv, _pi.csv, _icl.csv)
suppressMessages(library(sbm))
out <- commandArgs(TRUE)[1]; if (is.na(out)) out <- "dev/ref"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
draw <- function(model) function(n, m, sd) if (model == "poisson") rpois(n, m) else rnorm(n, m, sd)
plant1 <- function(n, Pm, sizes, directed, model, seed, pmiss = 0, sd = 1) {
  set.seed(seed); Q <- nrow(Pm); z <- sample(rep(seq_len(Q), times = sizes)); p <- Pm[z, z]
  A <- matrix(draw(model)(n * n, as.vector(p), sd), n, n)
  if (!directed) A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  M <- matrix(FALSE, n, n)
  if (pmiss > 0) { M[upper.tri(M)] <- runif(n * (n - 1) / 2) < pmiss
    if (directed) M[lower.tri(M)] <- runif(n * (n - 1) / 2) < pmiss else M <- M | t(M)
    A[M] <- NA }
  list(A = A, z = z, z2 = NULL)
}
plant2 <- function(n1, n2, Pm, s1, s2, model, seed, pmiss = 0, sd = 1) {
  set.seed(seed); z1 <- sample(rep(seq_len(nrow(Pm)), times = s1)); z2 <- sample(rep(seq_len(ncol(Pm)), times = s2))
  A <- matrix(draw(model)(n1 * n2, as.vector(Pm[z1, z2]), sd), n1, n2)
  if (pmiss > 0) A[matrix(runif(n1 * n2) < pmiss, n1, n2)] <- NA
  list(A = A, z = z1, z2 = z2)
}
write_net <- function(name, A, z, z2 = NULL, directed = FALSE) {
  idx <- which(!is.na(A) & A != 0, arr.ind = TRUE)
  if (!directed && is.null(z2)) idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  write.csv(data.frame(from = idx[, 1], to = idx[, 2], value = A[idx]), file.path(out, paste0(name, "_edges.csv")), row.names = FALSE)
  mi <- which(is.na(A), arr.ind = TRUE)
  if (!directed && is.null(z2)) mi <- mi[mi[, 1] < mi[, 2], , drop = FALSE]
  write.csv(data.frame(from = mi[, 1], to = mi[, 2]), file.path(out, paste0(name, "_missing.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_along(z), truth = z), file.path(out, paste0(name, "_truth.csv")), row.names = FALSE)
  if (!is.null(z2)) write.csv(data.frame(id = seq_along(z2), truth = z2), file.path(out, paste0(name, "_truth2.csv")), row.names = FALSE)
}
fit1 <- function(name, A, model, directed, qmax) {
  t0 <- proc.time()[3]
  f <- estimateSimpleSBM(A, model, directed = directed, estimOptions = list(verbosity = 0, plot = FALSE, nbCores = 1, exploreMax = qmax))
  el <- proc.time()[3] - t0
  s2 <- if (model == "gaussian") f$connectParam$var else NA
  write.csv(data.frame(Q = f$nbBlocks, ICL = f$ICL, n = nrow(A), sigma2 = s2, seconds = el), file.path(out, paste0(name, "_fit.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_len(nrow(A)), block = apply(f$probMemberships, 1, which.max)), file.path(out, paste0(name, "_membership.csv")), row.names = FALSE)
  write.csv(f$connectParam$mean, file.path(out, paste0(name, "_pi.csv")), row.names = FALSE)
  write.csv(f$storedModels, file.path(out, paste0(name, "_icl.csv")), row.names = FALSE)
}
fit2 <- function(name, A, model, qmax) {
  t0 <- proc.time()[3]
  f <- estimateBipartiteSBM(A, model, estimOptions = list(verbosity = 0, plot = FALSE, nbCores = 1, exploreMax = qmax))
  el <- proc.time()[3] - t0
  s2 <- if (model == "gaussian") f$connectParam$var else NA
  write.csv(data.frame(Q1 = f$nbBlocks[1], Q2 = f$nbBlocks[2], ICL = f$ICL, n1 = nrow(A), n2 = ncol(A), sigma2 = s2, seconds = el), file.path(out, paste0(name, "_fit.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_len(nrow(A)), block = apply(f$probMemberships$row, 1, which.max)), file.path(out, paste0(name, "_membership.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_len(ncol(A)), block = apply(f$probMemberships$col, 1, which.max)), file.path(out, paste0(name, "_membership2.csv")), row.names = FALSE)
  write.csv(f$connectParam$mean, file.path(out, paste0(name, "_pi.csv")), row.names = FALSE)
  write.csv(f$storedModels, file.path(out, paste0(name, "_icl.csv")), row.names = FALSE)
}
L3 <- matrix(1, 3, 3); diag(L3) <- c(8, 6, 5)
L3a <- matrix(c(7, 1, 0.5, 1, 6, 1.5, 0.5, 1.5, 5), 3)
L2 <- matrix(c(6, 0.8, 1.2, 4), 2)
sp <- list(
  list("po_u_n100_Q3", "poisson", 100, L3a, c(40, 30, 30), FALSE, 21, 0, 1),
  list("po_d_n100_Q2", "poisson", 100, L2, c(50, 50), TRUE, 22, 0, 1),
  list("po_miss_u_n80_Q2", "poisson", 80, L2, c(40, 40), FALSE, 23, 0.07, 1),
  list("ga_u_n100_Q3", "gaussian", 100, L3a, c(40, 30, 30), FALSE, 24, 0, 1.5),
  list("ga_d_n100_Q2", "gaussian", 100, L2, c(50, 50), TRUE, 25, 0, 1.5),
  list("ga_miss_u_n80_Q2", "gaussian", 80, L2, c(40, 40), FALSE, 26, 0.07, 1.5))
for (s in sp) {
  d <- plant1(s[[3]], s[[4]], s[[5]], s[[6]], s[[2]], s[[7]], s[[8]], s[[9]])
  write_net(s[[1]], d$A, d$z, NULL, s[[6]]); fit1(s[[1]], d$A, s[[2]], s[[6]], 5); cat(s[[1]], "\n")
}
B2 <- matrix(c(0.6, 0.1, 0.2, 0.15, 0.55, 0.1), 2, byrow = TRUE)
P2 <- matrix(c(6, 1, 2, 1.5, 5, 1), 2, byrow = TRUE)
sp2 <- list(
  list("lbm_ber_40x60", "bernoulli", 40, 60, B2, c(20, 20), c(20, 20, 20), 31, 0, 1),
  list("lbm_po_40x60", "poisson", 40, 60, P2, c(20, 20), c(20, 20, 20), 32, 0, 1),
  list("lbm_ga_40x60", "gaussian", 40, 60, P2, c(20, 20), c(20, 20, 20), 33, 0, 1.5),
  list("lbm_ber_miss_50x50", "bernoulli", 50, 50, B2, c(25, 25), c(15, 20, 15), 34, 0.06, 1))
for (s in sp2) {
  d <- plant2(s[[3]], s[[4]], s[[5]], s[[6]], s[[7]], s[[2]], s[[8]], s[[9]], s[[10]])
  if (s[[2]] == "bernoulli") d$A <- (d$A > 0) * 1 + 0 * d$A
  if (s[[2]] == "bernoulli") { set.seed(s[[8]]); z1 <- d$z; z2 <- d$z2; A <- matrix(rbinom(s[[3]] * s[[4]], 1, s[[5]][z1, z2]), s[[3]], s[[4]]); if (s[[9]] > 0) A[matrix(runif(s[[3]] * s[[4]]) < s[[9]], s[[3]], s[[4]])] <- NA; d$A <- A }
  write_net(s[[1]], d$A, d$z, d$z2); fit2(s[[1]], d$A, s[[2]], 5); cat(s[[1]], "\n")
}
