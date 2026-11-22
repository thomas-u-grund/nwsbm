# Summarise dev/restarts/grid_part*.csv: fast and default against the 20-restart reference (base R).
# Run from the repository root: Rscript dev/restarts/analyse.R [cells]
args <- commandArgs(TRUE)
d <- do.call(rbind, lapply(sort(Sys.glob("dev/restarts/grid_part*.csv")), read.csv))
d$dfd <- d$icl_fast - d$icl_def      # < 0: fast worse than default
d$drd <- d$icl_ref - d$icl_def       # < 0: default better than the 20-restart run
d$dfr <- d$icl_fast - d$icl_ref
d$Qagree <- d$Q_fast == d$Q_def
d$fw <- d$dfd < -0.5; d$fb <- d$dfd > 0.5; d$dw <- d$drd < -0.5; d$dwb <- d$drd > 0.5
summ <- function(x) { m <- function(v) mean(v, na.rm = TRUE); s <- function(v) sum(v, na.rm = TRUE)
  data.frame(runs = nrow(x), Qagree = m(x$Qagree), fast_worse = s(x$fw), fast_better = s(x$fb),
  ref_better_than_def = s(x$dwb), ari_df = m(x$ari_def_fast), ari_df_min = suppressWarnings(min(x$ari_def_fast, na.rm = TRUE)),
  ari_def_true = m(x$ari_def_true), ari_fast_true = m(x$ari_fast_true), ari_ref_true = m(x$ari_ref_true),
  min_dfd = suppressWarnings(min(x$dfd, na.rm = TRUE)), t_def = s(x$t_def), t_fast = s(x$t_fast), speedup = s(x$t_def) / s(x$t_fast)) }
bycols <- function(cols) { k <- unique(d[cols]); k <- k[do.call(order, k), , drop = FALSE]
  do.call(rbind, lapply(seq_len(nrow(k)), function(i) { s <- d; for (c in cols) s <- s[s[[c]] == k[i, c], ]; cbind(k[i, , drop = FALSE], summ(s)) })) }
options(width = 250); rownames(d) <- NULL
cat("runs", nrow(d), " missing", sum(is.na(d$icl_def) | is.na(d$icl_fast)), "\n")
cat("OVERALL\n"); print(round(summ(d), 3), row.names = FALSE)
for (k in list("level", "n", "K", "dir", c("n", "level"))) { cat("\nBY", k, "\n"); print(format(bycols(k), digits = 3), row.names = FALSE) }
if (length(args) && args[1] == "cells") { cat("\nBY CELL\n"); print(format(bycols(c("dir", "n", "K", "level")), digits = 3), row.names = FALSE) }
a <- d[d$Qagree, ]
cat("\nWhere Q agrees: ARI def/fast mean", mean(a$ari_def_fast), " share ARI>=0.99:", mean(a$ari_def_fast >= 0.99), " share |dICL|<0.5:", mean(abs(a$dfd) < 0.5), "\n")
cat("Q ref vs def agree", mean(d$Q_ref == d$Q_def), " fast vs ref", mean(d$Q_fast == d$Q_ref), " selected Q == true K: def", mean(d$Q_def == d$K), "fast", mean(d$Q_fast == d$K), "ref", mean(d$Q_ref == d$K), "\n")
cat("\nWorst 20 runs for fast (dICL = fast - default)\n")
w <- d[order(d$dfd), ][1:20, c("dir", "n", "K", "level", "rep", "Q_def", "Q_fast", "Q_ref", "icl_def", "icl_fast", "icl_ref", "ari_def_fast", "ari_def_true", "ari_fast_true", "t_def", "t_fast")]
print(format(w, digits = 5), row.names = FALSE)
write.csv(d, "dev/restarts/grid_all.csv", row.names = FALSE)
