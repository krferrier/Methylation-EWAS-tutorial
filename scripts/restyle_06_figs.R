## ---------------------------------------------------------------------------
## restyle_06_figs.R -- redraw ONLY data/06_manhattan.png and data/06_volcano.png.
##
## FIGURES ONLY. Both panels are drawn from checkpoints that already exist
## (06_ewas.rds, 08_annotation/annotated.rds), so no published number moves.
## Changes vs restyle_figs.R: larger base font and point size on both panels,
## and the volcano now encodes the sign of the methylation difference by colour
## (neg = teal, pos = plum, the pheno_pal keys already in _setup.R).
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(scales)
})
setwd("repo"); source("_setup.R")
set.seed(42)

DENSE <- 220
outd  <- file.path("data", "08_annotation")

save_fig <- function(p, path, width, height, dpi = DENSE) {
  ggsave(path, plot = p, width = width, height = height, dpi = dpi,
         device = grDevices::png, type = "cairo", bg = "white")
  cat(sprintf("  %-34s %5.1f x %-4.1f  %7.0f KB\n", path, width, height,
              file.size(path) / 1024))
}

teal_d <- unname(ewas_col["teal_dark"]);  teal  <- unname(ewas_col["teal"])
teal_l <- unname(ewas_col["teal_light"]); plum  <- unname(ewas_col["plum"])
gry_d  <- unname(ewas_col["gray_dark"]);  gry   <- unname(ewas_col["gray"])
gry_l  <- unname(ewas_col["gray_light"])
fam    <- if (nzchar(.ewas_family)) .ewas_family else NULL
cat("font family:", if (is.null(fam)) "(device default)" else fam, "\n")

BASE <- 16          # was 11
PT   <- 1.3         # Manhattan; was 0.3
PTV  <- 1.3         # volcano: same point size as the Manhattan panel

## ---- [06] Manhattan -------------------------------------------------------
ew  <- readRDS(data_path("06_ewas.rds"))
bs  <- readRDS(data_path("06_bacon_summary.rds"))
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
  geom_point(size = PT, alpha = 0.5) +
  geom_hline(yintercept = -log10(bonf), colour = plum, linewidth = 0.7) +
  annotate("text", x = max(pl$x), y = -log10(bonf), vjust = -0.55, hjust = 1,
           label = "Bonferroni", colour = plum, size = 4.6, fontface = "bold",
           family = fam) +
  scale_colour_manual(values = c(`0` = teal_d, `1` = teal_l), guide = "none") +
  scale_x_continuous(breaks = ctr$x, labels = ctr$chrn, expand = c(0.008, 0)) +
  scale_y_continuous(limits = c(0, max(-log10(bonf) * 1.08,
                                       -log10(min(pl$P.Value)) * 1.12)),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Chromosome", y = expression(-log[10](italic(p))),
       title = "PTSD EWAS, uncorrected limma p-values",
       subtitle = sprintf("%s CpGs tested in %d samples; nothing clears Bonferroni before BACON rescaling",
                          format(ew$n_tested, big.mark = ","), ew$n)) +
  theme_ewas(base_size = BASE) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = rel(0.8)))
save_fig(p06m, data_path("06_manhattan.png"), 11.0, 5.2)

## ---- [06] Volcano ---------------------------------------------------------
## x-axis is the limma model coefficient (logFC, M-value scale) -- the EWAS
## effect estimate, adjusted for every covariate in the design. The earlier
## draft used delta_beta, a raw case-minus-control group mean difference with
## no covariate adjustment; the two disagree in sign for 171,484 of the
## 756,251 probes, so the colour split moves with the axis.
##
## Colour carries the direction of effect; the three labelled top hits are
## ringed in neutral dark gray so the hue channel belongs to direction alone.
tt[, dir := factor(ifelse(logFC < 0, "neg", "pos"), levels = c("neg", "pos"))]
lab <- head(tt[order(P.Value)], 3)

p06v <- ggplot(tt, aes(logFC, -log10(P.Value))) +
  geom_vline(xintercept = 0, colour = gry_l, linewidth = 0.5) +
  geom_point(aes(colour = dir), size = PTV, alpha = 0.40) +
  geom_point(data = lab, shape = 21, fill = NA, colour = gry_d,
             size = 4.6, stroke = 1.2) +
  ggrepel::geom_text_repel(data = lab, aes(label = probe), size = 4.6,
                           colour = gry_d, min.segment.length = 0,
                           box.padding = 0.7, seed = 42, family = fam) +
  scale_colour_manual(
    values = c(neg = unname(pheno_pal["neg"]), pos = unname(pheno_pal["pos"])),
    labels = c(neg = "Lower in cases (hypomethylated)",
               pos = "Higher in cases (hypermethylated)"),
    name = NULL,
    guide = guide_legend(override.aes = list(size = 3.4, alpha = 1))) +
  scale_x_continuous(labels = label_number(style_positive = "plus",
                                           accuracy = 0.5),
                     expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = "EWAS effect size (logFC, M-value scale)",
       y = expression(-log[10](italic(p))),
       title = "Effect size against significance",
       subtitle = sprintf("Largest effect among the top hits is %.2f on the M-value scale \u2014 EWAS effects are small",
                          max(abs(lab$logFC)))) +
  theme_ewas(base_size = BASE) +
  theme(legend.position = "top",
        legend.text = element_text(size = rel(0.9)),
        legend.margin = margin(b = -2))
save_fig(p06v, data_path("06_volcano.png"), 9.8, 6.2)

## ---- [06] Manhattan, BACON-adjusted ---------------------------------------
## Same coordinates and the same Bonferroni threshold as the uncorrected panel,
## drawn on the BACON-rescaled p-values instead. This is the version that
## matches what the chapter tells you to report: the rescaling is what moves
## the top CpG across the line, so the two panels differ only in the y value.
plb <- merge(tt[, .(probe, bacon.p, bacon.adj.P)],
             pl[, .(probe, chrn, x)], by = "probe")
hit  <- plb[bacon.p    <  bonf]                    # Bonferroni after BACON
fdr  <- plb[bacon.adj.P < 0.05 & bacon.p >= bonf]  # FDR only
cat("\nBACON Manhattan: n_bonf =", nrow(hit), " n_fdr_only =", nrow(fdr),
    " max -log10(bacon.p) =", sprintf("%.2f", max(-log10(plb$bacon.p))), "\n")

p06mb <- ggplot(plb, aes(x, -log10(bacon.p), colour = factor(chrn %% 2))) +
  geom_point(size = PT, alpha = 0.5) +
  geom_hline(yintercept = -log10(bonf), colour = plum, linewidth = 0.7) +
  annotate("text", x = max(plb$x), y = -log10(bonf), vjust = -0.55, hjust = 1,
           label = "Bonferroni", colour = plum, size = 4.6, fontface = "bold",
           family = fam) +
  geom_point(data = fdr, shape = 21, fill = NA, colour = gry_d,
             size = 4.4, stroke = 1.1) +
  geom_point(data = hit, shape = 21, fill = plum, colour = "gray10",
             size = 4.4, stroke = 1.1) +
  ggrepel::geom_text_repel(data = rbind(hit[, .(probe, x, bacon.p)],
                                        fdr[, .(probe, x, bacon.p)]),
                           aes(x, -log10(bacon.p), label = probe),
                           inherit.aes = FALSE, size = 4.4, colour = gry_d,
                           min.segment.length = 0, box.padding = 0.8,
                           seed = 42, family = fam) +
  scale_colour_manual(values = c(`0` = teal_d, `1` = teal_l), guide = "none") +
  scale_x_continuous(breaks = ctr$x, labels = ctr$chrn, expand = c(0.008, 0)) +
  scale_y_continuous(limits = c(0, max(-log10(bonf) * 1.08,
                                       -log10(min(plb$bacon.p)) * 1.12)),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Chromosome", y = expression(-log[10](italic(p)))) +
  ggtitle("PTSD EWAS, BACON-adjusted p-values",
          subtitle = sprintf(
            "\u03bb %.3f \u2192 %.3f after rescaling; %d CpG clears Bonferroni, %d more passes FDR < 0.05 (open ring)",
            bs$lambda_raw, bs$lambda_bacon, nrow(hit), nrow(fdr))) +
  theme_ewas(base_size = BASE) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = rel(0.8)))
save_fig(p06mb, data_path("06_manhattan_bacon.png"), 11.0, 5.2)

cat("\nn_neg =", tt[dir == "neg", .N], " n_pos =", tt[dir == "pos", .N],
    " n_tested =", ew$n_tested, "\n")
cat("labelled:", paste(lab$probe, collapse = ", "), "\n")
cat("max |logFC| among labelled:", sprintf("%.4f", max(abs(lab$logFC))), "\n")
cat("logFC range:", sprintf("%.4f %.4f", min(tt$logFC), max(tt$logFC)), "\n")
cat("bonf =", format(bonf, digits = 4), " -log10 =", sprintf("%.2f", -log10(bonf)), "\n")
