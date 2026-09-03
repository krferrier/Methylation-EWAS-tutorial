## ---------------------------------------------------------------------------
## restyle_figs.R -- redraw every referenced tutorial figure in the site theme.
##
## FIGURES ONLY. Nothing here re-runs an analysis: every panel is drawn from a
## checkpoint that already exists on disk, so no frozen number can move. The
## only figures that are *authored* rather than restyled are the three that
## never had a generating script (00b_catalog_overview, 04_validation_scatter,
## 08_dmr_locus); those are rebuilt from their source data, and the DMR panel is
## additionally corrected from the stale 10-probe region to the v8.1 6-probe one.
##
## Not touched: data/CGI.png (user-supplied schematic), 06_bacon/*.jpg and
## 07_bacon/*.jpg (base-graphics diagnostics from the bacon object / pipeline,
## which is not checkpointed), and the nine figures no chapter references.
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(patchwork); library(scales)
})
setwd("repo"); source("_setup.R")
t0 <- Sys.time()
set.seed(42)

DPI   <- 300
DENSE <- 220            # dense point clouds: 300 dpi buys nothing but bytes
out   <- "data"
outd  <- file.path(out, "08_annotation")

## Every ggsave in this script goes through here so the device, background and
## font are identical across all fifteen figures.
save_fig <- function(p, path, width, height, dpi = DPI) {
  ggsave(path, plot = p, width = width, height = height, dpi = dpi,
         device = grDevices::png, type = "cairo", bg = "white")
  cat(sprintf("  %-46s %5.1f x %-4.1f  %7.0f KB\n", path, width, height,
              file.size(path) / 1024))
}

## Palette shorthands (values come from _setup.R, which the site theme shares).
teal_d <- unname(ewas_col["teal_dark"]);  teal  <- unname(ewas_col["teal"])
teal_l <- unname(ewas_col["teal_light"]); sand  <- unname(ewas_col["sand"])
plum   <- unname(ewas_col["plum"]);       gry_d <- unname(ewas_col["gray_dark"])
gry    <- unname(ewas_col["gray"]);       gry_l <- unname(ewas_col["gray_light"])

cat("font family:", if (nzchar(.ewas_family)) .ewas_family else "(device default)", "\n")
cat("figures ->", normalizePath(out), "\n\n")

## ===========================================================================
## 00b -- catalog overview (AUTHORED: no generating script existed)
## Source: methylation_geo_catalog.csv, the 70-study screen described in ch00b.
## ===========================================================================
cat("[00b] catalog overview\n")
cat0 <- as.data.table(read.csv("methylation_geo_catalog.csv", check.names = FALSE))
stopifnot(nrow(cat0) == 70)

tier_lv <- c("Strong", "Moderate", "Limited")
tiers <- cat0[, .N, by = .(tier = factor(Suitability_tier, levels = tier_lv))][order(tier)]
pa <- ggplot(tiers, aes(tier, N, fill = tier)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = N), vjust = -0.4, size = 3.2, colour = gry_d,
            family = .ewas_family) +
  scale_fill_manual(values = c(Strong = teal, Moderate = teal_l, Limited = gry_l),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Teaching-suitability tier", x = NULL, y = "Studies")

fields <- data.table(
  field = c("Sex", "Age", "Tissue / cell", "Phenotype or exposure",
            "Ancestry / race", "Cell-type composition"),
  n = c(sum(cat0$Has_sex == "True"), sum(cat0$Has_age == "True"),
        sum(cat0$Has_tissue == "True"), sum(cat0$N_phenotype_fields > 0),
        sum(cat0$Has_ancestry == "True"), sum(cat0$Has_cell_composition == "True")))
fields[, field := factor(field, levels = rev(field))]
## Accent the two fields a first EWAS most wants and least often gets.
fields[, key := n < 20]
pb <- ggplot(fields, aes(n, field, fill = key)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = n), hjust = -0.28, size = 3.2, colour = gry_d,
            family = .ewas_family) +
  scale_fill_manual(values = c(`TRUE` = sand, `FALSE` = teal), guide = "none") +
  scale_x_continuous(limits = c(0, 70), expand = expansion(mult = c(0, 0.10))) +
  labs(title = "Metadata fields present", x = "Studies (of 70)", y = NULL)

pc <- ggplot(cat0, aes(N_samples)) +
  geom_histogram(bins = 18, fill = teal_l, colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = median(cat0$N_samples), colour = plum,
             linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = median(cat0$N_samples) * 1.08, y = Inf,
           label = sprintf("median %s", format(round(median(cat0$N_samples)), big.mark = ",")),
           hjust = 0, vjust = 1.6, size = 3.1, colour = plum, family = .ewas_family) +
  scale_x_log10(breaks = c(150, 300, 1000, 3000),
                labels = c("150", "300", "1k", "3k")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = sprintf("Study size (%s-%s samples)",
                       min(cat0$N_samples), format(max(cat0$N_samples), big.mark = ",")),
       x = "Samples per study (log scale)", y = "Studies")

arr_lv <- c("EPIC (850K)", "450K+EPIC (850K)", "EPIC (850K)+EPIC v2")
arrs <- cat0[, .N, by = .(arr = factor(Meth_array, levels = arr_lv))][order(arr)]
arrs[, lab := c("EPIC v1\nonly", "450K +\nEPIC", "EPIC +\nEPIC v2")[match(arr, arr_lv)]]
arrs[, lab := factor(lab, levels = lab)]
pd <- ggplot(arrs, aes(lab, N)) +
  geom_col(width = 0.6, fill = teal) +
  geom_text(aes(label = N), vjust = -0.4, size = 3.2, colour = gry_d,
            family = .ewas_family) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Array generation",
       subtitle = sprintf("%d studies are multi-omic super-series",
                          sum(nzchar(cat0$Other_platforms))),
       x = NULL, y = "Studies")

p00b <- (pa | pb) / (pc | pd) +
  plot_annotation(
    title = "A catalogue of 70 EPIC-array GEO series screened for EWAS teaching",
    tag_levels = "a",
    theme = theme_ewas() + theme(
      plot.title = element_text(size = rel(1.15), face = "bold", colour = teal_d))) &
  theme(plot.tag = element_text(face = "bold", size = rel(1.0), colour = gry_d))
save_fig(p00b, data_path("00b_catalog_overview.png"), 11, 7.6)

## ---------------------------------------------------------------------------
## [04] validation scatter -- our estimates vs the study's deposited values.
## Authored: this figure never had a generating script.
## ---------------------------------------------------------------------------
val <- readRDS(data_path("04_validation.rds"))
mg  <- as.data.table(val$merged); vv <- as.data.table(val$validation)
cells <- vv$cell
long <- rbindlist(lapply(cells, function(cl)
  data.table(cell = cl, est = mg[[cl]], dep = mg[[tolower(cl)]])))
## Facet strips carry r and MAE so the reader never has to cross-reference the
## table above it in the chapter.
vv[, lab := sprintf("%s\nr = %.3f, MAE = %.1f pp", cell, r, 100 * mae)]
long[, lab := vv$lab[match(cell, vv$cell)]]
long[, lab := factor(lab, levels = vv$lab[order(-vv$r)])]

p04 <- ggplot(long, aes(100 * dep, 100 * est)) +
  geom_abline(slope = 1, intercept = 0, colour = gry_l, linetype = 2,
              linewidth = 0.45) +
  geom_point(colour = teal, size = 1.5, alpha = 0.75) +
  facet_wrap(~ lab, scales = "free", nrow = 2) +
  scale_x_continuous(labels = label_number(suffix = "%")) +
  scale_y_continuous(labels = label_number(suffix = "%")) +
  labs(x = "Deposited proportion (study's own estimate)",
       y = "Our estimate (IDOL / preprocessNoob)",
       title = "Our cell-type estimates reproduce the study's deposited proportions",
       subtitle = sprintf("n = %d samples; dashed line is perfect agreement", nrow(mg))) +
  theme_ewas(base_size = 11) +
  theme(panel.spacing = unit(1.05, "lines"))
save_fig(p04, data_path("04_validation_scatter.png"), 9.5, 6.0)

## ---------------------------------------------------------------------------
## [05a] PC association heatmap
## ---------------------------------------------------------------------------
pca <- readRDS(data_path("05_batch_pca.rds"))
r2  <- pca$r2; pve <- pca$pve
nice <- c(chip_slide = "Chip (slide)", array_position = "Array position",
          PTSD = "PTSD", sex = "Sex", age = "Age",
          Neu = "Neutrophil", CD4T = "CD4 T", CD8T = "CD8 T")
hm <- as.data.table(as.table(r2)); setnames(hm, c("variable", "PC", "r2"))
hm[, variable := factor(nice[as.character(variable)], levels = rev(unname(nice)))]
hm[, PC := factor(PC, levels = paste0("PC", 1:10))]

p05a <- ggplot(hm, aes(PC, variable, fill = r2)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sub("^0", "", sprintf("%.2f", r2)),
                colour = r2 > 0.45), size = 2.9, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = gry_d)) +
  scale_fill_gradient(low = "#F4F2EE", high = teal_d, limits = c(0, 1),
                      name = expression(R^2),
                      guide = guide_colourbar(barheight = unit(3.2, "cm"))) +
  scale_x_discrete(position = "top", labels = function(x)
    sprintf("%s\n%.1f%%", x, 100 * pve[match(x, paste0("PC", 1:10))])) +
  labs(x = NULL, y = NULL,
       title = "What the leading principal components are made of",
       subtitle = "Cell composition owns PC1; chip and sex are entangled on PC5; PTSD is nowhere",
       caption = "Percentages under each label are the variance that component explains.") +
  theme_ewas(base_size = 11) +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.text.x.top = element_text(colour = gry_d, lineheight = 1.15))
save_fig(p05a, data_path("05_pc_heatmap.png"), 9.6, 4.6)

## ---------------------------------------------------------------------------
## [05b] PC1-PC2 coloured by neutrophil fraction
## ---------------------------------------------------------------------------
sc <- data.table(PC1 = pca$pcs[, 1], PC2 = pca$pcs[, 2],
                 Neu = 100 * pca$props[, "Neu"])
p05b <- ggplot(sc, aes(PC1, PC2, colour = Neu)) +
  geom_point(size = 2.7, alpha = 0.92) +
  scale_colour_gradient(low = "#DCE9EA", high = teal_d,
                        name = "Neutrophil\nfraction",
                        labels = label_number(suffix = "%")) +
  labs(x = sprintf("PC1 (%.1f%% of variance)", 100 * pve[1]),
       y = sprintf("PC2 (%.1f%%)", 100 * pve[2]),
       title = "PC1 of a blood methylome is immune-cell mixture",
       subtitle = "Each point is one array, coloured by its estimated neutrophil fraction") +
  theme_ewas(base_size = 11)
save_fig(p05b, data_path("05_pc_scatter_neu.png"), 7.4, 5.0)

## ---------------------------------------------------------------------------
## [05c] what each surrogate variable captures
## ---------------------------------------------------------------------------
sva_o <- readRDS(data_path("05_sva.rds"))
SV <- sva_o$SV; mdk <- sva_o$mdk
colnames(SV) <- paste0("SV", seq_len(ncol(SV)))
sv_levels <- colnames(SV)          # hold the labels outside the data.table: inside
                                   # bars[...], the name `SV` resolves to the COLUMN,
                                   # so colnames(SV) would be NULL and every level NA.
r2_on <- function(x) sapply(seq_len(ncol(SV)),
                            function(j) summary(lm(SV[, j] ~ x))$r.squared)
bars <- rbindlist(list(
  data.table(SV = sv_levels, r2 = r2_on(factor(mdk$slide)), what = "Chip (technical)"),
  data.table(SV = sv_levels, r2 = r2_on(factor(mdk$sex)),   what = "Sex (biological)"),
  data.table(SV = sv_levels, r2 = r2_on(factor(mdk$ptsd)),  what = "PTSD (exposure)")))
bars[, SV := factor(SV, levels = sv_levels)]
bars[, what := factor(what, levels = c("Chip (technical)", "Sex (biological)",
                                       "PTSD (exposure)"))]
p05c <- ggplot(bars, aes(SV, r2, fill = what)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  scale_fill_manual(values = c("Chip (technical)" = teal_d,
                               "Sex (biological)" = sand,
                               "PTSD (exposure)"  = plum), name = NULL) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL, y = expression(R^2),
       title = "The surrogate variables absorb chip, not the exposure",
       subtitle = sprintf(paste("Buja-Eyuboglu suggested %d SVs on the full filtered matrix;",
                                "the %d carried into the model are shown"),
                          sva_o$prev_n_sv, ncol(SV))) +
  theme_ewas(base_size = 11) + theme(legend.position = "top")
save_fig(p05c, data_path("05_sva_bars.png"), 8.4, 4.4)

## ---------------------------------------------------------------------------
## [06] Manhattan + volcano.
## Coordinates come from annotated.rds (Zhou hg38), not the hg19 minfi
## annotation the first draft used -- the rest of the tutorial is hg38.
## ---------------------------------------------------------------------------
ew  <- readRDS(data_path("06_ewas.rds"))
tt  <- as.data.table(ew$tt)
ann <- readRDS(file.path(outd, "annotated.rds"))
pl  <- merge(tt[, .(probe, P.Value, delta_beta)],
             ann[, .(probe, chrm = CpG_chrm, pos = CpG_beg)],
             by = "probe")
pl[, chrn := suppressWarnings(as.integer(sub("^chr", "", chrm)))]
pl <- pl[!is.na(chrn) & !is.na(pos)]
setorder(pl, chrn, pos)
chrmax <- pl[, .(m = max(pos)), by = chrn][order(chrn)]
chrmax[, off := cumsum(as.numeric(m)) - m]
pl[chrmax, x := pos + i.off, on = "chrn"]
ctr  <- pl[, .(x = (min(x) + max(x)) / 2), by = chrn][order(chrn)]
bonf <- 0.05 / ew$n_tested

p06m <- ggplot(pl, aes(x, -log10(P.Value), colour = factor(chrn %% 2))) +
  geom_point(size = 0.3, alpha = 0.7) +
  geom_hline(yintercept = -log10(bonf), colour = plum, linewidth = 0.5) +
  annotate("text", x = max(pl$x), y = -log10(bonf), vjust = -0.6, hjust = 1,
           label = "Bonferroni", colour = plum, size = 3,
           family = if (nzchar(.ewas_family)) .ewas_family else NULL) +
  scale_colour_manual(values = c(`0` = teal_d, `1` = teal_l), guide = "none") +
  scale_x_continuous(breaks = ctr$x, labels = ctr$chrn, expand = c(0.008, 0)) +
  scale_y_continuous(limits = c(0, max(7, -log10(min(pl$P.Value)) * 1.12)),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Chromosome", y = expression(-log[10](italic(p))),
       title = "PTSD EWAS, uncorrected limma p-values",
       subtitle = sprintf("%s CpGs tested in %d samples; nothing clears Bonferroni before BACON rescaling",
                          format(ew$n_tested, big.mark = ","), ew$n)) +
  theme_ewas(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = rel(0.8)))
save_fig(p06m, data_path("06_manhattan.png"), 10.5, 4.6, dpi = DENSE)

## Volcano: label only what the prose actually discusses.
tt[, dbp := 100 * delta_beta]
lab <- head(tt[order(P.Value)], 3)
p06v <- ggplot(tt, aes(dbp, -log10(P.Value))) +
  geom_vline(xintercept = 0, colour = gry_l, linewidth = 0.4) +
  geom_point(size = 0.3, alpha = 0.45, colour = teal) +
  geom_point(data = lab, colour = plum, size = 2) +
  ggrepel::geom_text_repel(data = lab, aes(label = probe), size = 3,
                           colour = plum, min.segment.length = 0,
                           box.padding = 0.6, seed = 42,
                           family = if (nzchar(.ewas_family)) .ewas_family else "") +
  scale_x_continuous(labels = label_number(style_positive = "plus",
                                           suffix = " pp")) +
  labs(x = "Methylation difference, case \u2212 control",
       y = expression(-log[10](italic(p))),
       title = "Effect size against significance",
       subtitle = sprintf("Largest shift among the top hits is %.1f percentage points \u2014 EWAS effects are small",
                          max(abs(lab$dbp)))) +
  theme_ewas(base_size = 11)
save_fig(p06v, data_path("06_volcano.png"), 7.6, 5.2, dpi = DENSE)

## ---------------------------------------------------------------------------
## [07] pipeline figures. Read-only over the pipeline's result CSVs -- the
## Snakemake run itself is not re-executed, so every p-value is as-run.
## ---------------------------------------------------------------------------
D <- "../ewas_pipeline/run_grady"
DA <- "../ewas_pipeline/run_grady_all"
lam <- function(p) { p <- p[is.finite(p) & p > 0]
                     median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
qqdt <- function(p, lab) { p <- sort(p[is.finite(p) & p > 0])
  data.table(analysis = lab, obs = -log10(p), exp = -log10(ppoints(length(p)))) }

A <- fread(file.path(DA, "PTSD_ewas_bacon_results.csv.gz"))
Fm <- fread(file.path(D, "F/F_PTSD_ewas_bacon_results.csv.gz"))
Mm <- fread(file.path(D, "M/M_PTSD_ewas_bacon_results.csv.gz"))
me <- fread(file.path(D, "PTSD_ewas_meta_analysis_results_1.txt"))
setnames(me, "P-value", "P")

## [07a] meta QQ, with a 95% concentration band so "on the line" is quantified.
qm  <- qqdt(me$P, "Meta (F+M)")
nq  <- nrow(qm)
ci  <- data.table(exp = qm$exp, i = seq_len(nq))
ci[, `:=`(lo = -log10(qbeta(0.975, i, nq - i + 1)),
          hi = -log10(qbeta(0.025, i, nq - i + 1)))]
p07a <- ggplot() +
  geom_ribbon(data = ci[seq(1, nq, length.out = 4000)],
              aes(exp, ymin = lo, ymax = hi), fill = gry_l, alpha = 0.35) +
  geom_abline(slope = 1, intercept = 0, colour = gry, linetype = 2, linewidth = 0.4) +
  geom_point(data = qm, aes(exp, obs), size = 0.5, colour = teal_d, alpha = 0.7) +
  labs(x = expression(Expected~-log[10](italic(p))),
       y = expression(Observed~-log[10](italic(p))),
       title = "Meta-analysis QQ, BACON-adjusted inverse-variance",
       subtitle = sprintf("%s CpGs, \u03bb = %.3f; shaded band is the 95%% null concentration",
                          format(nrow(me), big.mark = ","), lam(me$P))) +
  theme_ewas(base_size = 11)
save_fig(p07a, data_path("07_meta_qq.png"), 6.4, 5.2, dpi = DENSE)

## [07b] four-way overlay. Thin to a log-spaced index set: 756k x 4 is unreadable.
thin <- function(d) { n <- nrow(d)
  i <- unique(round(c(1:2000, exp(seq(log(2001), log(n), length.out = 4000)))))
  d[i[i <= n]] }
lv <- c(sprintf("Overall (n=%d)", A$n[1]), sprintf("Female (n=%d)", Fm$n[1]),
        sprintf("Male (n=%d)", Mm$n[1]), "Meta (F+M)")
ov <- rbindlist(list(thin(qqdt(A$bacon.pval,  lv[1])),
                     thin(qqdt(Fm$bacon.pval, lv[2])),
                     thin(qqdt(Mm$bacon.pval, lv[3])), thin(qm)))
ov[, analysis := factor(analysis, levels = lv)]
p07b <- ggplot(ov, aes(exp, obs, colour = analysis)) +
  geom_abline(slope = 1, intercept = 0, colour = gry, linetype = 2, linewidth = 0.4) +
  geom_point(size = 0.55, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05 / nrow(me)), colour = plum,
             linetype = 3, linewidth = 0.5) +
  annotate("text", x = 0, y = -log10(0.05 / nrow(me)), hjust = 0, vjust = -0.6,
           label = "Bonferroni", colour = plum, size = 3,
           family = if (nzchar(.ewas_family)) .ewas_family else NULL) +
  scale_colour_manual(values = setNames(c(gry_d, sand, teal_l, plum), lv),
                      name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.4))) +
  labs(x = expression(Expected~-log[10](italic(p))),
       y = expression(Observed~-log[10](italic(p))),
       title = "All four analyses, BACON-adjusted",
       subtitle = "Splitting by sex halves the sample and costs power; the meta-analysis recovers it") +
  theme_ewas(base_size = 11) + theme(legend.position = "top")
save_fig(p07b, data_path("07_qq_overlay.png"), 7.2, 5.6, dpi = DENSE)

## [07c] female vs male effect concordance.
J <- merge(Fm[, .(cpgid, F_es = bacon.es, F_p = bacon.pval)],
           Mm[, .(cpgid, M_es = bacon.es, M_p = bacon.pval)], by = "cpgid")
bonf7 <- 0.05 / nrow(me)
hits  <- me[P < bonf7, .(cpgid = MarkerName, Direction, P)]
J <- merge(J, hits, by = "cpgid", all.x = TRUE)
J[, sig := !is.na(P)]
r_all <- cor(J$F_es, J$M_es)
set.seed(42)
bg <- J[sig == FALSE][sample(.N, min(60000, .N))]
dir_lab <- c(`++` = "Both up", `--` = "Both down",
             `+-` = "Opposite (F up)", `-+` = "Opposite (M up)")
hi <- J[sig == TRUE]; hi[, dl := dir_lab[Direction]]
p07c <- ggplot() +
  geom_hline(yintercept = 0, colour = gry_l, linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = gry_l, linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, colour = gry_l, linetype = 2,
              linewidth = 0.4) +
  geom_point(data = bg, aes(F_es, M_es), colour = gry_l, size = 0.3, alpha = 0.35) +
  geom_point(data = hi, aes(F_es, M_es, colour = dl), size = 2.8) +
  scale_colour_manual(values = c("Both up" = plum, "Both down" = teal_d,
                                 "Opposite (F up)" = sand,
                                 "Opposite (M up)" = teal_l), name = NULL) +
  labs(x = "Female stratum effect (BACON-adjusted)",
       y = "Male stratum effect (BACON-adjusted)",
       title = "Per-CpG PTSD effects, female versus male",
       subtitle = sprintf("Pearson r = %.3f genome-wide; %d Bonferroni meta hits highlighted",
                          r_all, nrow(hits)),
       caption = "Gray cloud is a 60,000-CpG random subsample of the non-significant probes; dashed line is y = x.") +
  theme_ewas(base_size = 11) + theme(legend.position = "top")
save_fig(p07c, data_path("07_effect_concordance.png"), 7.0, 6.0, dpi = DENSE)

## ---------------------------------------------------------------------------
## [08a/b] genomic context of the top 1,000 CpGs vs everything tested.
## ---------------------------------------------------------------------------
setorder(ann, bacon.p)
TOPN <- 1000
topd <- head(ann, TOPN)
mk <- function(col, lev) {
  d <- rbind(data.table(set = "Top 1,000",  grp = factor(topd[[col]], levels = lev)),
             data.table(set = "All tested", grp = factor(ann[[col]],  levels = lev))
  )[!is.na(grp), .N, by = .(set, grp)]
  d[, pct := 100 * N / sum(N), by = set]
  d[, set := factor(set, levels = c("Top 1,000", "All tested"))][]
}
ctx_plot <- function(d, title, sub) {
  ggplot(d, aes(grp, pct, fill = set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    geom_text(aes(label = sprintf("%.1f", pct)),
              position = position_dodge(width = 0.8), vjust = -0.45,
              size = 2.9, colour = gry_d,
              family = if (nzchar(.ewas_family)) .ewas_family else "") +
    scale_fill_manual(values = c("Top 1,000" = plum, "All tested" = teal), name = NULL) +
    scale_y_continuous(labels = label_number(suffix = "%"),
                       expand = expansion(mult = c(0, 0.13))) +
    labs(x = NULL, y = "CpGs in category", title = title, subtitle = sub) +
    theme_ewas(base_size = 11) +
    theme(legend.position = "top", panel.grid.major.x = element_blank())
}
fd <- mk("feature", c("Promoter", "Gene body", "Intergenic"))
idl <- c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea")
id <- mk("island", idl)
id[, grp := factor(c(Island = "Island", N_Shore = "N shore", S_Shore = "S shore",
                     N_Shelf = "N shelf", S_Shelf = "S shelf",
                     OpenSea = "Open sea")[as.character(grp)],
                   levels = c("Island", "N shore", "S shore", "N shelf",
                              "S shelf", "Open sea"))]
isl_top <- id[set == "Top 1,000" & grp == "Island", pct]
isl_all <- id[set == "All tested" & grp == "Island", pct]

save_fig(ctx_plot(fd, "Gene-feature context of the strongest signals",
                  sprintf("Zhou promoter (\u00b11.5 kb of TSS) / gene-body / intergenic scheme; %s CpGs tested",
                          format(nrow(ann), big.mark = ","))),
         file.path(outd, "08_feature_distribution.png"), 8.2, 4.4)
save_fig(ctx_plot(id, "CpG-island context of the strongest signals",
                  sprintf("The top set leans toward islands (%.1f%% vs %.1f%%) \u2014 a lean this size is not evidence on its own",
                          isl_top, isl_all)),
         file.path(outd, "08_island_distribution.png"), 8.8, 4.4)

## ---------------------------------------------------------------------------
## [08c] gometh GO + KEGG
## ---------------------------------------------------------------------------
gm <- readRDS(file.path(outd, "08_gometh.rds"))
go <- as.data.table(gm$go); kg <- as.data.table(gm$kegg)
pd <- rbind(
  data.table(collection = "GO",   term = go$TERM[1:10],        P = go$P.DE[1:10]),
  data.table(collection = "KEGG", term = kg$Description[1:10], P = kg$P.DE[1:10]))
pd[, term := ifelse(nchar(term) > 48, paste0(substr(term, 1, 45), "\u2026"), term)]
pd[, lab := factor(term, levels = rev(term))]
p08c <- ggplot(pd, aes(-log10(P), lab, fill = collection)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, colour = plum,
             linewidth = 0.45) +
  annotate("text", x = -log10(0.05), y = 0.6, label = "nominal p = 0.05",
           hjust = -0.08, size = 2.9, colour = plum,
           family = if (nzchar(.ewas_family)) .ewas_family else NULL) +
  scale_fill_manual(values = c(GO = teal_d, KEGG = sand), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
  facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression(-log[10]~italic(p)~"(nominal)"), y = NULL,
       title = "gometh enrichment of the top 1,000 CpGs",
       subtitle = sprintf("Nothing survives multiple testing: smallest FDR is %.2f (GO) and %.2f (KEGG)",
                          min(go$FDR), min(kg$FDR))) +
  theme_ewas(base_size = 10) +
  theme(legend.position = "none", panel.grid.major.y = element_blank())
save_fig(p08c, file.path(outd, "08_gometh_enrichment.png"), 9.2, 6.0)

## ---------------------------------------------------------------------------
## [08d] KYCG
## ---------------------------------------------------------------------------
ky <- as.data.table(readRDS(file.path(outd, "08_kycg.rds")))
kp <- head(ky[is.finite(estimate)][order(FDR)], 12)
kp[, grp := sub("^KYCG\\.EPIC\\.", "", sub("\\.[0-9]+$", "", group))]
kp[, lab := factor(dbname, levels = rev(dbname))]
p08d <- ggplot(kp, aes(estimate, lab, fill = grp)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(chromHMM = teal_d, TFBSconsensus = plum,
                               HMconsensus = sand, CGI = teal_l),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
  labs(x = expression(log[2]~"odds ratio"), y = NULL,
       title = "KYCG enrichment of the top 1,000 CpGs",
       subtitle = sprintf(paste("Twelve most significant features, drawn from %d of the four EPIC",
                                "knowledgebases; all clear FDR < 0.05 (max %.1e)"),
                          uniqueN(kp$grp), max(kp$FDR))) +
  theme_ewas(base_size = 10) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())
save_fig(p08d, file.path(outd, "08_kycg_enrichment.png"), 8.8, 5.2)

## ---------------------------------------------------------------------------
## [08e] the DMR locus. Authored, and corrected: the file on disk was the stale
## ten-probe region from the pre-v8.1 run. The v8.1 comb-p call is
## chr6:28,633,534-28,633,601 with six probes.
## ---------------------------------------------------------------------------
reg <- fread(cmd = sprintf("zcat %s",
             file.path(outd, "dmr/PTSD_dmr_demo.regions-p.bed.gz")))
setnames(reg, c("chrom", "start", "end", "min_p", "n_probes", "z_p", "z_sidak_p"))
PAD <- 260
loc <- ann[CpG_chrm == reg$chrom[1] &
           CpG_beg >= reg$start[1] - PAD & CpG_end <= reg$end[1] + PAD,
           .(probe, pos = CpG_beg, p = bacon.p, db = delta_beta)]
loc[, inreg := pos >= reg$start[1] & pos <= reg$end[1]]
stopifnot(sum(loc$inreg) == reg$n_probes[1])
sidak <- reg$z_sidak_p[1]

p08e <- ggplot(loc, aes(pos, -log10(p))) +
  annotate("rect", xmin = reg$start[1], xmax = reg$end[1], ymin = -Inf, ymax = Inf,
           fill = teal_l, alpha = 0.16) +
  geom_hline(yintercept = -log10(sidak), colour = plum, linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05 / nrow(ann)), colour = gry,
             linetype = 2, linewidth = 0.4) +
  geom_segment(aes(xend = pos, yend = 0, colour = inreg), linewidth = 0.6) +
  geom_point(aes(colour = inreg, shape = db < 0), size = 2.4) +
  annotate("text", x = min(loc$pos), y = -log10(sidak), hjust = 0, vjust = -0.55,
           size = 3, colour = plum,
           label = sprintf("region \u0160id\u00e1k p = %.2e", sidak),
           family = if (nzchar(.ewas_family)) .ewas_family else NULL) +
  annotate("text", x = min(loc$pos), y = -log10(0.05 / nrow(ann)), hjust = 0,
           vjust = 1.5, size = 3, colour = gry,
           label = "single-CpG genome-wide threshold",
           family = if (nzchar(.ewas_family)) .ewas_family else NULL) +
  scale_colour_manual(values = c(`TRUE` = teal_d, `FALSE` = gry_l), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 25, `FALSE` = 24),
                     labels = c(`TRUE` = "Hypomethylated in cases",
                                `FALSE` = "Hypermethylated in cases"), name = NULL) +
  scale_y_continuous(limits = c(0, -log10(sidak) * 1.12),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(breaks = scales::breaks_width(50),
                     expand = expansion(mult = c(0.05, 0.05)),
                     labels = function(x) format(x, big.mark = ",", trim = TRUE,
                                                 scientific = FALSE)) +
  labs(x = sprintf("Position on %s (hg38)", reg$chrom[1]),
       y = expression(-log[10](italic(p))),
       title = sprintf("Six adjacent CpGs make a region no single CpG could make"),
       subtitle = sprintf("%s:%s\u2013%s; the shaded %d bp window is what comb-p reports.\nThe strongest single CpG in it reaches only \u2212log\u2081\u2080 p = %.2f",
                          reg$chrom[1], format(reg$start[1], big.mark = ","),
                          format(reg$end[1], big.mark = ","),
                          reg$end[1] - reg$start[1],
                          max(-log10(loc[inreg == TRUE]$p)))) +
  theme_ewas(base_size = 11) +
  theme(legend.position = "top",
        plot.margin = margin(5.5, 22, 5.5, 5.5))
save_fig(p08e, file.path(outd, "08_dmr_locus.png"), 9.0, 5.2)

cat("\nall figures done in",
    round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
