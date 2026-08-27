## Regenerate chapter 07 artifacts from the v8.1 full-probe pipeline run.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
t0 <- Sys.time()
D  <- "ewas_pipeline/run_grady"
out <- "repo/data"
dir.create(file.path(out, "07_bacon"), showWarnings = FALSE, recursive = TRUE)

## The pipeline writes its own BACON diagnostics; copy them in under the names
## the chapter references rather than redrawing them.
file.copy(file.path(D, "all/bacon_plots/PTSD_fit.jpg"),
          file.path(out, "07_bacon/all_fit.jpg"), overwrite = TRUE)
file.copy(file.path(D, "all/bacon_plots/PTSD_qqs.jpg"),
          file.path(out, "07_bacon/all_qqs.jpg"), overwrite = TRUE)
file.copy(file.path(D, "strat/F/bacon_plots/F_PTSD_qqs.jpg"),
          file.path(out, "07_bacon/F_qqs.jpg"), overwrite = TRUE)
file.copy(file.path(D, "strat/M/bacon_plots/M_PTSD_qqs.jpg"),
          file.path(out, "07_bacon/M_qqs.jpg"), overwrite = TRUE)

lam <- function(p) { p <- p[is.finite(p) & p > 0]
                     median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }

A <- fread(file.path(D, "all/PTSD_ewas_bacon_results.csv.gz"))
F <- fread(file.path(D, "strat/F/F_PTSD_ewas_bacon_results.csv.gz"))
M <- fread(file.path(D, "strat/M/M_PTSD_ewas_bacon_results.csv.gz"))
m <- fread(file.path(D, "meta_analysis/PTSD_ewas_meta_analysis_results_1.txt"))
setnames(m, "P-value", "P")

## ---- summary table -------------------------------------------------------
row1 <- function(lab, d) {
  setorder(d, bacon.pval)
  data.table(Analysis = lab, n = d$n[1],
             lambda_raw = lam(d$p.value), lambda_bacon = lam(d$bacon.pval),
             n_P_lt_1e5 = sum(d$bacon.pval < 1e-5, na.rm = TRUE),
             top_cpg = d$cpgid[1], top_P = d$bacon.pval[1])
}
summ <- rbindlist(list(
  row1("Overall", A), row1("Female", F), row1("Male", M),
  data.table(Analysis = "Meta (F+M)", n = NA_integer_,
             lambda_raw = NA_real_, lambda_bacon = lam(m$P),
             n_P_lt_1e5 = sum(m$P < 1e-5, na.rm = TRUE),
             top_cpg = m$MarkerName[which.min(m$P)], top_P = min(m$P))))
fwrite(summ, file.path(out, "07_pipeline_summary.csv"))
print(summ)

## ---- meta QQ -------------------------------------------------------------
qqdt <- function(p, lab) {
  p <- sort(p[is.finite(p) & p > 0])
  data.table(analysis = lab, obs = -log10(p),
             exp = -log10(ppoints(length(p))))
}
qm <- qqdt(m$P, "Meta (F+M)")
pq <- ggplot(qm, aes(exp, obs)) +
  geom_abline(slope = 1, intercept = 0, color = "gray60", linetype = 2) +
  geom_point(size = 0.5, color = "#4C72B0", alpha = 0.6) +
  labs(x = expression(Expected~-log[10](p)), y = expression(Observed~-log[10](p)),
       title = "Meta-analysis QQ (BACON-adjusted, inverse-variance)",
       subtitle = sprintf("%s CpGs, lambda = %.3f", format(nrow(m), big.mark = ","), lam(m$P))) +
  theme_minimal(base_size = 11)
ggsave(file.path(out, "07_meta_qq.png"), plot = pq, width = 6, height = 5, dpi = 200)

## ---- four-way QQ overlay -------------------------------------------------
## Thin to a log-spaced index set: 756k points x 4 panels is unreadable and slow.
thin <- function(d) { n <- nrow(d)
  idx <- unique(round(c(1:2000, exp(seq(log(2001), log(n), length.out = 4000)))))
  d[idx[idx <= n]] }
ov <- rbindlist(lapply(list(
  thin(qqdt(A$bacon.pval, "Overall (n=87)")),
  thin(qqdt(F$bacon.pval, "Female (n=45)")),
  thin(qqdt(M$bacon.pval, "Male (n=42)")),
  thin(qm)), identity))
ov[, analysis := factor(analysis, levels = c("Overall (n=87)", "Female (n=45)",
                                             "Male (n=42)", "Meta (F+M)"))]
po <- ggplot(ov, aes(exp, obs, color = analysis)) +
  geom_abline(slope = 1, intercept = 0, color = "gray60", linetype = 2) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05 / nrow(m)), color = "#C44E52", linetype = 3) +
  scale_color_manual(values = c("Overall (n=87)" = "gray30", "Female (n=45)" = "#DD8452",
                                 "Male (n=42)" = "#4C72B0", "Meta (F+M)" = "#C44E52"),
                      name = NULL) +
  labs(x = expression(Expected~-log[10](p)), y = expression(Observed~-log[10](p)),
       title = "All four analyses, BACON-adjusted",
       subtitle = "dotted line = Bonferroni threshold") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(out, "07_qq_overlay.png"), plot = po, width = 7, height = 5.5, dpi = 200)

## ---- effect concordance --------------------------------------------------
J <- merge(F[, .(cpgid, F_es = bacon.es, F_p = bacon.pval)],
           M[, .(cpgid, M_es = bacon.es, M_p = bacon.pval)], by = "cpgid")
bonf <- 0.05 / nrow(m)
hits <- m[P < bonf, .(cpgid = MarkerName, Direction, P)]
J <- merge(J, hits, by = "cpgid", all.x = TRUE)
J[, sig := !is.na(P)]
r_all <- cor(J$F_es, J$M_es)
## Plot a random subsample of the null cloud; keep every hit.
set.seed(42)
bg <- J[sig == FALSE][sample(.N, min(60000, .N))]
pc <- ggplot() +
  geom_hline(yintercept = 0, color = "gray80") + geom_vline(xintercept = 0, color = "gray80") +
  geom_point(data = bg, aes(F_es, M_es), color = "gray75", size = 0.3, alpha = 0.4) +
  geom_point(data = J[sig == TRUE], aes(F_es, M_es, color = Direction), size = 2.6) +
  scale_color_manual(values = c(`++` = "#C44E52", `--` = "#4C72B0",
                                 `+-` = "#DD8452", `-+` = "#8172B3"), name = "Direction") +
  labs(x = "Female stratum effect (BACON-adjusted)",
       y = "Male stratum effect (BACON-adjusted)",
       title = "Per-CpG PTSD effect sizes, female vs male",
       subtitle = sprintf("Pearson r = %.3f genome-wide; colored = %d Bonferroni meta hits",
                          r_all, nrow(hits))) +
  theme_minimal(base_size = 11)
ggsave(file.path(out, "07_effect_concordance.png"), plot = pc, width = 6.5, height = 5.5, dpi = 200)

## ---- numbers the prose needs --------------------------------------------
tt <- merge(hits, J[, .(cpgid, F_es, F_p, M_es, M_p)], by = "cpgid")
setorder(tt, P)
cat("\n--- Bonferroni meta hits ---\n"); print(tt)
cat("\nbonf thresh:", bonf, "| n bonf:", nrow(hits),
    "| n FDR<0.05:", sum(p.adjust(m$P, "BH") < 0.05), "\n")
cat("one-arm-only (other arm p>0.05):", sum(tt$F_p > 0.05 | tt$M_p > 0.05), "of", nrow(tt), "\n")
cat("cor F/M effects:", round(r_all, 4), "\n")
cat("sign agreement:", round(mean(sign(J$F_es) == sign(J$M_es)), 4), "\n")
cat("het P<0.05:", sum(m$HetPVal < 0.05, na.rm = TRUE),
    sprintf("(%.2f%%)", 100 * mean(m$HetPVal < 0.05, na.rm = TRUE)), "\n")
cat("median HetISq:", median(m$HetISq, na.rm = TRUE), "| mean:", round(mean(m$HetISq, na.rm = TRUE), 2), "\n")
cat("median SE F:", round(median(F$std.error), 4), "M:", round(median(M$std.error), 4),
    "overall:", round(median(A$std.error), 4), "\n")
saveRDS(list(summ = summ, hits = tt, r_all = r_all, bonf = bonf,
             het_frac = mean(m$HetPVal < 0.05, na.rm = TRUE),
             n_bonf = nrow(hits), n_fdr = sum(p.adjust(m$P, "BH") < 0.05)),
        file.path(out, "07_pipeline_extra.rds"))
cat("\nelapsed:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
