# Reference fits with dyadic covariates for nwsbm: R sbm 1.0.0 (blockmodels 1.2.1), estimateSimpleSBM / estimateBipartiteSBM with covariates = list(...).
# Run: Rscript dev/ref_cov.R dev/ref [fixture-index]   (no index: all fixtures; R needs about a minute per fixture)
# Writes <name>_edges.csv (from,to,value), _missing.csv, _truth.csv (_truth2.csv), _cov<k>.csv (the k-th covariate matrix, no header),
# _fit.csv (Q, ICL, n, sigma2, seconds), _beta.csv, _membership.csv, _pi.csv, _icl.csv.
suppressMessages(library(sbm))
args <- commandArgs(TRUE); out <- args[1]; if (is.na(out)) out <- "dev/ref"
out_dir <- out
only <- if (length(args) > 1) as.integer(args[2]) else NA
dir.create(out, showWarnings = FALSE, recursive = TRUE)
sym <- function(M) { M[lower.tri(M)] <- t(M)[lower.tri(M)]; M }
mkcov <- function(n, directed, kind, seed) {   # kind: "absdiff" (symmetric, from a node attribute), "same" (symmetric, group indicator), "rand" (asymmetric)
  set.seed(seed)
  if (kind == "absdiff") { x <- runif(n); abs(outer(x, x, "-")) }
  else if (kind == "same") { g <- sample(1:3, n, TRUE); 1 * outer(g, g, "==") }
  else { M <- matrix(rnorm(n * n), n, n); if (!directed) M <- sym(M); M }
}
gen <- function(model, eta, sd) {  # eta: linear predictor matrix
  nn <- length(eta)
  y <- switch(model, bernoulli = rbinom(nn, 1, plogis(eta)), poisson = rpois(nn, exp(eta)), gaussian = rnorm(nn, eta, sd))
  matrix(y, nrow(eta), ncol(eta))
}
write_all <- function(name, A, z, covs, z2 = NULL, directed = FALSE) {
  twomode <- !is.null(z2)
  idx <- which(!is.na(A) & A != 0, arr.ind = TRUE)
  if (!directed && !twomode) idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  write.csv(data.frame(from = idx[, 1], to = idx[, 2], value = A[idx]), file.path(out, paste0(name, "_edges.csv")), row.names = FALSE)
  mi <- which(is.na(A), arr.ind = TRUE)
  if (!directed && !twomode) mi <- mi[mi[, 1] < mi[, 2], , drop = FALSE]
  write.csv(data.frame(from = mi[, 1], to = mi[, 2]), file.path(out, paste0(name, "_missing.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_along(z), truth = z), file.path(out, paste0(name, "_truth.csv")), row.names = FALSE)
  if (twomode) write.csv(data.frame(id = seq_along(z2), truth = z2), file.path(out, paste0(name, "_truth2.csv")), row.names = FALSE)
  for (k in seq_along(covs)) write.table(covs[[k]], file.path(out, paste0(name, "_cov", k, ".csv")), sep = ",", row.names = FALSE, col.names = FALSE)
}
write_glm <- function(name, A, covs, model, directed, twomode) {   # independent reference for the one-block model (logistic / Poisson / linear regression on the observed pairs)
  ok <- !is.na(A)
  if (!twomode) { diag(ok) <- FALSE; if (!directed) ok[lower.tri(ok)] <- FALSE }
  d <- data.frame(y = A[ok]); for (k in seq_along(covs)) d[[paste0("x", k)]] <- covs[[k]][ok]
  f <- as.formula(paste("y ~", paste(paste0("x", seq_along(covs)), collapse = " + ")))
  g <- if (model == "gaussian") lm(f, d) else glm(f, data = d, family = if (model == "bernoulli") binomial else poisson)
  ll <- as.numeric(logLik(g))
  out <- data.frame(loglik = ll, nobs = sum(ok), intercept = unname(coef(g)[1]))
  for (k in seq_along(covs)) out[[paste0("beta", k)]] <- unname(coef(g)[k + 1])
  if (model == "gaussian") out$sigma2 <- sum(resid(g)^2) / sum(ok)
  write.csv(out, file.path(out_dir, paste0(name, "_glm.csv")), row.names = FALSE)
}
opt <- function(qmax) list(verbosity = 0, plot = FALSE, nbCores = 1, exploreMax = qmax, fast = FALSE)  # fast = TRUE (the sbm default for Bernoulli with covariates) is an approximation, see GOTCHA
fit_one <- function(name, A, model, directed, covs, qmax) {
  t0 <- proc.time()[3]
  f <- estimateSimpleSBM(A, model, directed = directed, covariates = covs, estimOptions = opt(qmax))
  el <- proc.time()[3] - t0
  s2 <- if (model == "gaussian") f$connectParam$var else NA
  write.csv(data.frame(Q = f$nbBlocks, ICL = f$ICL, n = nrow(A), sigma2 = s2, seconds = el), file.path(out, paste0(name, "_fit.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_len(nrow(A)), block = apply(f$probMemberships, 1, which.max)), file.path(out, paste0(name, "_membership.csv")), row.names = FALSE)
  write.csv(f$connectParam$mean, file.path(out, paste0(name, "_pi.csv")), row.names = FALSE)
  write.csv(data.frame(beta = as.vector(f$covarParam)), file.path(out, paste0(name, "_beta.csv")), row.names = FALSE)
  write.csv(f$storedModels, file.path(out, paste0(name, "_icl.csv")), row.names = FALSE)
}
fit_two <- function(name, A, model, covs, qmax) {
  t0 <- proc.time()[3]
  f <- estimateBipartiteSBM(A, model, covariates = covs, estimOptions = opt(qmax))
  el <- proc.time()[3] - t0
  s2 <- if (model == "gaussian") f$connectParam$var else NA
  write.csv(data.frame(Q1 = f$nbBlocks[1], Q2 = f$nbBlocks[2], ICL = f$ICL, n1 = nrow(A), n2 = ncol(A), sigma2 = s2, seconds = el), file.path(out, paste0(name, "_fit.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_len(nrow(A)), block = apply(f$probMemberships$row, 1, which.max)), file.path(out, paste0(name, "_membership.csv")), row.names = FALSE)
  write.csv(data.frame(id = seq_len(ncol(A)), block = apply(f$probMemberships$col, 1, which.max)), file.path(out, paste0(name, "_membership2.csv")), row.names = FALSE)
  write.csv(f$connectParam$mean, file.path(out, paste0(name, "_pi.csv")), row.names = FALSE)
  write.csv(data.frame(beta = as.vector(f$covarParam)), file.path(out, paste0(name, "_beta.csv")), row.names = FALSE)
  write.csv(f$storedModels, file.path(out, paste0(name, "_icl.csv")), row.names = FALSE)
}
# NOTE: R blockmodels 1.2.1 treats missing ties (NA) as zeros in the Bernoulli and Poisson covariate models (checked against glm), so the fixtures that are compared
# with the R fits have no missing ties, except the Gaussian one; the fixtures in `mis` (missing ties) have only the independent glm reference of the one-block model.
# one-mode: name, model, n, block logit/log/mean matrix M (on the link scale), sizes, directed, seed, pmiss, sd, covariate kinds, beta
one <- list(
  list("cb_u_n100_Q3", "bernoulli", 100, matrix(qlogis(c(.35,.06,.10,.06,.30,.08,.10,.08,.28)), 3), c(40,30,30), FALSE, 41, 0, 1, c("absdiff"), c(-2.0)),
  list("cb_d_n80_Q2", "bernoulli", 80, matrix(qlogis(c(.35,.08,.12,.30)), 2), c(40,40), TRUE, 42, 0, 1, c("rand", "same"), c(0.8, 1.0)),
  list("cp_u_n80_Q3", "poisson", 80, matrix(log(c(2.5,.4,.8,.4,2,.6,.8,.6,1.6)), 3), c(30,25,25), FALSE, 43, 0, 1, c("absdiff"), c(-1.2)),
  list("cp_d_n70_Q2", "poisson", 70, matrix(log(c(2.5,.5,.8,2)), 2), c(35,35), TRUE, 44, 0, 1, c("rand", "same"), c(0.4, 0.6)),
  list("cg_u_n80_Q2", "gaussian", 80, matrix(c(1.5,0,0.3,1.0), 2), c(40,40), FALSE, 45, 0, 1, c("absdiff"), c(-1.5)),
  list("cg_d_n70_Q3", "gaussian", 70, matrix(c(2,0,0.5,0,1.5,0,0.5,0,1), 3), c(25,25,20), TRUE, 46, 0.04, 1, c("same"), c(0.8)))
two <- list(
  list("cl_b_40x50", "bernoulli", 40, 50, matrix(qlogis(c(.45,.08,.12,.4,.1,.1)), 2, byrow = TRUE), c(20,20), c(15,20,15), 51, 0, 1, c("rand"), c(0.7)),
  list("cl_g_40x50", "gaussian", 40, 50, matrix(c(2,0,0.7,0,1.5,0.2), 2, byrow = TRUE), c(20,20), c(15,20,15), 52, 0.05, 1, c("rand"), c(0.8)),
  list("cl_p_40x50", "poisson", 40, 50, matrix(log(c(3,.5,1,.6,2.5,.7)), 2, byrow = TRUE), c(20,20), c(15,20,15), 53, 0, 1, c("rand"), c(0.4)))
k <- 0
for (s in one) {
  k <- k + 1; if (!is.na(only) && only != k) next
  name <- s[[1]]; model <- s[[2]]; n <- s[[3]]; M <- s[[4]]; directed <- s[[6]]; seed <- s[[7]]; pmiss <- s[[8]]
  covs <- lapply(seq_along(s[[10]]), function(i) mkcov(n, directed, s[[10]][i], seed * 10 + i))
  set.seed(seed); z <- sample(rep(seq_len(nrow(M)), times = s[[5]]))
  eta <- M[z, z]; for (i in seq_along(covs)) eta <- eta + s[[11]][i] * covs[[i]]
  A <- gen(model, eta, s[[9]]); if (!directed) A <- sym(A); diag(A) <- 0
  if (pmiss > 0) { Mm <- matrix(FALSE, n, n); Mm[upper.tri(Mm)] <- runif(n * (n - 1) / 2) < pmiss
    if (directed) Mm[lower.tri(Mm)] <- runif(n * (n - 1) / 2) < pmiss else Mm <- Mm | t(Mm); A[Mm] <- NA }
  write_all(name, A, z, covs, NULL, directed); write_glm(name, A, covs, model, directed, FALSE); fit_one(name, A, model, directed, covs, 4); cat(name, "\n")
}
for (s in two) {
  k <- k + 1; if (!is.na(only) && only != k) next
  name <- s[[1]]; model <- s[[2]]; n1 <- s[[3]]; n2 <- s[[4]]; M <- s[[5]]; seed <- s[[8]]; pmiss <- s[[9]]
  covs <- lapply(seq_along(s[[11]]), function(i) mkcov(n1, TRUE, s[[11]][i], seed * 10 + i)[, seq_len(n2) %% n1 + 1, drop = FALSE])
  set.seed(seed); z1 <- sample(rep(seq_len(nrow(M)), times = s[[6]])); z2 <- sample(rep(seq_len(ncol(M)), times = s[[7]]))
  eta <- M[z1, z2]; for (i in seq_along(covs)) eta <- eta + s[[12]][i] * covs[[i]]
  A <- gen(model, eta, s[[10]])
  if (pmiss > 0) A[matrix(runif(n1 * n2) < pmiss, n1, n2)] <- NA
  write_all(name, A, z1, covs, z2, FALSE); write_glm(name, A, covs, model, FALSE, TRUE); fit_two(name, A, model, covs, 4); cat(name, "\n")
}
# fixtures with missing ties: only the glm reference of the one-block model (R's Bernoulli and Poisson covariate models do not mask missing ties)
mis <- list(
  list("cbm_d_n60_Q2", "bernoulli", 60, matrix(qlogis(c(.35,.08,.12,.30)), 2), c(30,30), TRUE, 61, 0.06, 1, c("rand", "absdiff"), c(0.8, -1.5)),
  list("cpm_u_n60_Q2", "poisson", 60, matrix(log(c(2.5,.5,.8,2)), 2), c(30,30), FALSE, 62, 0.06, 1, c("absdiff"), c(-1.0)))
for (s in mis) {
  k <- k + 1; if (!is.na(only) && only != k) next
  name <- s[[1]]; model <- s[[2]]; n <- s[[3]]; M <- s[[4]]; directed <- s[[6]]; seed <- s[[7]]; pmiss <- s[[8]]
  covs <- lapply(seq_along(s[[10]]), function(i) mkcov(n, directed, s[[10]][i], seed * 10 + i))
  set.seed(seed); z <- sample(rep(seq_len(nrow(M)), times = s[[5]]))
  eta <- M[z, z]; for (i in seq_along(covs)) eta <- eta + s[[11]][i] * covs[[i]]
  A <- gen(model, eta, s[[9]]); if (!directed) A <- sym(A); diag(A) <- 0
  Mm <- matrix(FALSE, n, n); Mm[upper.tri(Mm)] <- runif(n * (n - 1) / 2) < pmiss
  if (directed) Mm[lower.tri(Mm)] <- runif(n * (n - 1) / 2) < pmiss else Mm <- Mm | t(Mm); A[Mm] <- NA
  write_all(name, A, z, covs, NULL, directed); write_glm(name, A, covs, model, directed, FALSE); cat(name, "\n")
}
