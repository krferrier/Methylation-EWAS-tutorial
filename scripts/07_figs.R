## 07_figs.R -- draws the five figures chapter 07 displays from the EWAS
## pipeline's result files. No analysis is re-run.
##   tutorial/data/07_meta_qq.png            meta-analysis QQ
##   tutorial/data/07_qq_overlay.png         QQ of all four analyses, overlaid
##   tutorial/data/07_volcano_facets.png     one volcano per analysis
##   tutorial/data/07_effect_concordance.png female vs male effect, meta hits
##   tutorial/data/07_foothills.png          top CpGs across all four analyses
##
## Run from the repository root after the pipeline run, with its output in
## ewas_pipeline/ (the folder chapter 07 reads as ../ewas_pipeline):
##   Rscript scripts/07_figs.R
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel)
})
D   <- "ewas_pipeline/run_grady"
DA  <- "ewas_pipeline/run_grady_all"
out <- "tutorial/data"

A <- fread(file.path(DA, "PTSD_ewas_bacon_results.csv.gz"))      # overall
F <- fread(file.path(D, "F/F_PTSD_ewas_bacon_results.csv.gz"))   # female
M <- fread(file.path(D, "M/M_PTSD_ewas_bacon_results.csv.gz"))   # male
m <- fread(file.path(D, "PTSD_ewas_meta_analysis_results_1.txt"))  # METAL
setnames(m, "P-value", "P")

# Labels carry each arm's sample size, read from the pipeline's n column;
# every figure uses the same label and colour for each analysis.
LV  <- c(sprintf("Overall (n=%d)", A$n[1]), sprintf("Female (n=%d)", F$n[1]),
         sprintf("Male (n=%d)", M$n[1]), "Meta (F+M)")
PAL <- setNames(c("gray30", "#DD8452", "#4C72B0", "#C44E52"), LV)

# Bonferroni threshold over all tested CpGs, the only threshold drawn
BONF <- 0.05 / nrow(m)

# genomic inflation factor
lam <- function(p) { p <- p[is.finite(p) & p > 0]
                     median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }

## ================= FIGURE 1: meta-analysis QQ ==============================
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

## ================= FIGURE 2: QQ overlay of all four analyses ==============
## Thin to a log-spaced index set: 756k points x 4 panels is unreadable and slow.
thin <- function(d) { n <- nrow(d)
  idx <- unique(round(c(1:2000, exp(seq(log(2001), log(n), length.out = 4000)))))
  d[idx[idx <= n]] }
ov <- rbindlist(lapply(list(
  thin(qqdt(A$bacon.pval, LV[1])),
  thin(qqdt(F$bacon.pval, LV[2])),
  thin(qqdt(M$bacon.pval, LV[3])),
  thin(qm)), identity))
ov[, analysis := factor(analysis, levels = LV)]
po <- ggplot(ov, aes(exp, obs, color = analysis)) +
  geom_abline(slope = 1, intercept = 0, color = "gray60", linetype = 2) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(BONF), color = "#C44E52", linetype = 3) +
  scale_color_manual(values = PAL, name = NULL) +
  labs(x = expression(Expected~-log[10](p)), y = expression(Observed~-log[10](p)),
       title = "All four analyses, BACON-adjusted",
       subtitle = "dotted line = Bonferroni threshold") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(out, "07_qq_overlay.png"), plot = po, width = 7, height = 5.5, dpi = 200)

## ---- long-format table of all four analyses ------------------------------
L <- rbindlist(list(
  A[, .(analysis = LV[1], cpg = cpgid, es = bacon.es, se = bacon.se, p = bacon.pval)],
  F[, .(analysis = LV[2], cpg = cpgid, es = bacon.es, se = bacon.se, p = bacon.pval)],
  M[, .(analysis = LV[3], cpg = cpgid, es = bacon.es, se = bacon.se, p = bacon.pval)],
  m[, .(analysis = LV[4], cpg = MarkerName, es = Effect, se = StdErr, p = P)]))
L <- L[is.finite(p) & p > 0 & is.finite(es)]
L[, analysis := factor(analysis, levels = LV)]
L[, nlp := -log10(p)]

## ================= FIGURE 3: faceted volcano ==============================
## Keep every CpG with p < 1e-4 so no point near the top of any panel is lost,
## then thin the null cloud below that so four dense panels stay legible and
## the PNG stays a reasonable size. 1e-4 is a plotting cutoff for point
## retention only -- it is not drawn and it is not a significance threshold.
set.seed(42)
keep <- L[p < 1e-4]
bgv  <- L[p >= 1e-4][, .SD[sample(.N, min(60000, .N))], by = analysis]
V    <- rbind(keep, bgv)
V[, tier := ifelse(p < BONF, "Bonferroni", "Not significant")]
V[, tier := factor(tier, levels = c("Bonferroni", "Not significant"))]

nsig <- L[, .(n_bonf = sum(p < BONF)), by = analysis]
nsig[, lab := sprintf("Bonferroni: %d", n_bonf)]
print(nsig)
xr <- range(V$es); yr <- range(V$nlp)

p_volc <- ggplot(V, aes(es, nlp)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "gray45") +
  geom_hline(yintercept = -log10(BONF), colour = "#C44E52") +
  geom_point(data = V[tier == "Not significant"], colour = "gray72",
             size = 0.35, alpha = 0.45) +
  geom_point(data = V[tier == "Bonferroni"], shape = 21, fill = "#C44E52",
             colour = "gray10", stroke = 0.35, size = 2.6) +
  ggrepel::geom_text_repel(data = V[tier == "Bonferroni"], aes(label = cpg),
                           size = 2.4, min.segment.length = 0, segment.size = 0.3,
                           segment.colour = "gray40", box.padding = 0.4,
                           max.overlaps = 30, seed = 42) +
  geom_label(data = nsig, aes(x = xr[2], y = yr[2], label = lab),
             hjust = 1, vjust = 1, size = 2.9, label.size = 0.25,
             fill = "white") +
  facet_wrap(~ analysis, nrow = 1) +
  scale_x_continuous(expand = expansion(mult = 0.06)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = "PTSD effect on M-value (BACON-adjusted)",
       y = expression(-log[10](p)),
       title = "Where each analysis puts its signal",
       subtitle = "red line and red points = Bonferroni; the null cloud is thinned to 60,000 points per panel") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "gray90", colour = "gray60"),
        strip.text = element_text(size = 10),
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9, colour = "gray30"))
ggsave(file.path(out, "07_volcano_facets.png"), plot = p_volc,
       width = 11, height = 4.2, dpi = 200)

## ================= FIGURE 4: F-vs-M effect concordance ====================
## One point per CpG the meta-analysis calls Bonferroni-significant: its female
## effect against its male effect, with each stratum's standard error drawn.
J <- merge(F[, .(cpg = cpgid, F_es = bacon.es, F_se = bacon.se, F_p = bacon.pval)],
           M[, .(cpg = cpgid, M_es = bacon.es, M_se = bacon.se, M_p = bacon.pval)],
           by = "cpg")
hits <- m[P < BONF, .(cpg = MarkerName, meta_es = Effect, meta_p = P,
                      dirs = Direction, hetisq = HetISq)]
S <- merge(hits, J, by = "cpg")[order(meta_p)]
S[, concordant := sign(F_es) == sign(M_es)]
cat("\nmeta-Bonferroni CpGs:", nrow(S),
    "| concordant in sign:", sum(S$concordant), "\n")
print(S[, .(cpg, meta_p = signif(meta_p, 3), dirs,
            F_es = round(F_es, 3), M_es = round(M_es, 3),
            hetisq = round(hetisq, 1))])
rho <- if (nrow(S) > 2) cor(S$F_es, S$M_es) else NA_real_
cat("Pearson r (F vs M effect, meta hits):", round(rho, 3), "\n")
lim <- range(c(S$F_es - S$F_se, S$F_es + S$F_se,
               S$M_es - S$M_se, S$M_es + S$M_se)) * 1.10

p_scat <- ggplot(S, aes(F_es, M_es)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "gray55") +
  geom_vline(xintercept = 0, linetype = 2, colour = "gray55") +
  geom_abline(slope = 1, intercept = 0, colour = "gray20") +
  geom_errorbar(aes(ymin = M_es - M_se, ymax = M_es + M_se),
                width = 0, colour = "gray65", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = F_es - F_se, xmax = F_es + F_se),
                 height = 0, colour = "gray65", linewidth = 0.4) +
  geom_point(aes(fill = concordant), shape = 21, colour = "gray15",
             stroke = 0.4, size = 3) +
  geom_text_repel(aes(label = cpg), size = 2.6, min.segment.length = 0,
                  segment.size = 0.3, segment.colour = "gray40",
                  box.padding = 0.55, max.overlaps = 30, seed = 42) +
  scale_fill_manual(values = c(`TRUE` = "#4C72B0", `FALSE` = "#C44E52"),
                    labels = c(`TRUE` = "same sign", `FALSE` = "opposite sign"),
                    name = "Effect direction across strata") +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "Female stratum effect (BACON-adjusted)",
       y = "Male stratum effect",
       title = sprintf("The %d meta-analysis hits are driven by one stratum, not both",
                       nrow(S)),
       subtitle = paste0("bars are \u00b11 SE; the solid line is y = x, where both strata ",
                         "would have equal effects (Pearson r = ",
                         sprintf("%.2f", rho), ")")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9, colour = "gray30"))
ggsave(file.path(out, "07_effect_concordance.png"), plot = p_scat,
       width = 6.8, height = 6.4, dpi = 200)

## ================= FIGURE 5: foothills ====================================
## Union of the top 8 CpGs from each analysis, ordered by strongest evidence
## anywhere, with every analysis' estimate for that CpG shown.
top8 <- L[order(p), head(.SD, 8), by = analysis]$cpg
top8 <- unique(top8)
Ft <- L[cpg %in% top8]
ord <- Ft[, .(best = min(p)), by = cpg][order(best)]
Ft[, cpg := factor(cpg, levels = ord$cpg)]
Ft[, direction := ifelse(es > 0, "Hypermethylated in cases",
                                 "Hypomethylated in cases")]
cat("\nfoothills CpGs:", length(top8), "| rows:", nrow(Ft), "\n")

p_foot <- ggplot(Ft, aes(cpg, nlp, colour = analysis, shape = direction)) +
  geom_hline(yintercept = -log10(BONF), colour = "#C44E52") +
  geom_point(size = 2.4, alpha = 0.9) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = c("Hypermethylated in cases" = 24,
                                "Hypomethylated in cases"  = 25), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  labs(x = NULL, y = expression(-log[10](p)),
       title = "Every top CpG is carried by one analysis, not all four",
       subtitle = "union of each analysis' top 8 CpGs; the red line is Bonferroni") +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape  = guide_legend(order = 2, nrow = 2)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "gray92"),
        axis.text.x = element_text(angle = 55, hjust = 1, size = 7.5),
        legend.position = "right", legend.key.height = unit(0.9, "lines"),
        plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9, colour = "gray30"))
ggsave(file.path(out, "07_foothills.png"), plot = p_foot,
       width = 10, height = 5.2, dpi = 200)
