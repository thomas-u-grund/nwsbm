# restarts(0), (2), (4) against the default (10) and the 20-restart reference. Run from the repository root.
g1 <- read.csv("dev/restarts/grid_all.csv")
g2 <- do.call(rbind, lapply(sort(Sys.glob("dev/restarts/grid2*_part*.csv")), read.csv))
d <- merge(g1, g2, by = c("dir", "n", "K", "level", "pin", "pout", "rep"))
cat("runs", nrow(d), "\n")
arm <- function(nm, Q, icl, t, ari = NULL) {
  dd <- icl - d$icl_def
  data.frame(arm = nm, Qagree_def = mean(Q == d$Q_def), icl_worse_05 = sum(dd < -0.5), icl_worse_05_pct = round(100 * mean(dd < -0.5), 2), icl_better_05 = sum(dd > 0.5),
    worst = min(dd), time_sum = sum(t), speedup = sum(d$t_def) / sum(t))
}
r <- rbind(arm("restarts(0)", d$Q_r0, d$icl_r0, d$t_r0), arm("restarts(2)", d$Q_r2, d$icl_r2, d$t_r2), arm("restarts(4)", d$Q_r4, d$icl_r4, d$t_r4),
           arm("fast", d$Q_fast, d$icl_fast, d$t_fast), arm("restarts(20)", d$Q_ref, d$icl_ref, d$t_ref))
options(width = 200); print(format(r, digits = 4), row.names = FALSE)
for (lv in c("easy", "mid", "hard", "linux")) { cat("\nlevel", lv, "\n"); k <- d$level == lv
  print(format(data.frame(arm = c("r0", "r2", "r4", "fast", "r20"), Qagree = c(mean(d$Q_r0[k] == d$Q_def[k]), mean(d$Q_r2[k] == d$Q_def[k]), mean(d$Q_r4[k] == d$Q_def[k]), mean(d$Q_fast[k] == d$Q_def[k]), mean(d$Q_ref[k] == d$Q_def[k])),
    worse05 = c(sum(d$icl_r0[k] - d$icl_def[k] < -.5), sum(d$icl_r2[k] - d$icl_def[k] < -.5), sum(d$icl_r4[k] - d$icl_def[k] < -.5), sum(d$icl_fast[k] - d$icl_def[k] < -.5), sum(d$icl_ref[k] - d$icl_def[k] < -.5)),
    speedup = c(sum(d$t_def[k]) / sum(d$t_r0[k]), sum(d$t_def[k]) / sum(d$t_r2[k]), sum(d$t_def[k]) / sum(d$t_r4[k]), sum(d$t_def[k]) / sum(d$t_fast[k]), sum(d$t_def[k]) / sum(d$t_ref[k]))), digits = 3), row.names = FALSE) }
cat("\nr4 worse cases\n"); w <- d[d$icl_r4 - d$icl_def < -.5, c("dir", "n", "K", "level", "rep", "Q_def", "Q_r4", "icl_def", "icl_r4", "icl_ref")]; print(w)
cat("\nr2 worse by >0.5, count by n:\n"); print(table(d$n[d$icl_r2 - d$icl_def < -.5]))
cat("\nBy n: share Q agrees with default / ICL worse by >0.5 (r0, r2, r4) and speed-up of r2\n")
for (nn in sort(unique(d$n))) { k <- d$n == nn
  cat(sprintf("n=%d runs=%d  Qagree r0 %.3f r2 %.3f r4 %.3f | worse>0.5 r0 %d r2 %d r4 %d | speedup r2 %.2f  r0 %.2f\n", nn, sum(k), mean(d$Q_r0[k] == d$Q_def[k]), mean(d$Q_r2[k] == d$Q_def[k]), mean(d$Q_r4[k] == d$Q_def[k]),
    sum(d$icl_r0[k] - d$icl_def[k] < -.5), sum(d$icl_r2[k] - d$icl_def[k] < -.5), sum(d$icl_r4[k] - d$icl_def[k] < -.5), sum(d$t_def[k]) / sum(d$t_r2[k]), sum(d$t_def[k]) / sum(d$t_r0[k]))) }
cat("\nMean ARI vs the planted partition: default", mean(d$ari_def_true), " r0", mean(d$ari_r0_true), " r2", mean(d$ari_r2_true), " r4", mean(d$ari_r4_true), " r20", mean(d$ari_ref_true), "\n")
cat("Mean ARI between r2 and r4 partitions:", mean(d$ari_r2_r4, na.rm = TRUE), "; between default and r20:", mean(d$ari_def_ref, na.rm = TRUE), "\n")
q <- d$Q_r2 == d$Q_def
cat("Q agrees (r2 vs default): ARI(r2,true) - ARI(def,true) mean", mean(d$ari_r2_true[q] - d$ari_def_true[q]), " min", min(d$ari_r2_true[q] - d$ari_def_true[q]), "; share ICL within 0.5:", mean(abs(d$icl_r2[q] - d$icl_def[q]) < .5), "\n")
cat("Time per network (sum over runs, seconds): default", sum(d$t_def), " r2", sum(d$t_r2), " r0", sum(d$t_r0), "\n")
