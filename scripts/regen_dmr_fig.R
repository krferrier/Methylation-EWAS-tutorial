## Redraw 08_dmr_locus.png against the v8.1 refit and the corrected comb-p run
## (seed 1e-4). Source of p-values is 06_ewas.rds, joined to YAME v8.1
## coordinates -- the same BED comb-p was actually given.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(scales)})
setwd("repo"); source("_setup.R")
set.seed(42)
DPI <- 300
outd <- "data/08_annotation"

save_fig <- function(p, path, width, height, dpi = DPI) {
  ggsave(path, plot = p, width = width, height = height, dpi = dpi,
         device = grDevices::png, type = "cairo", bg = "white")
  cat(sprintf("  %-46s %5.1f x %-4.1f  %7.0f KB\n", path, width, height,
              file.size(path) / 1024))
}
teal_d <- unname(ewas_col["teal_dark"]);  teal  <- unname(ewas_col["teal"])
teal_l <- unname(ewas_col["teal_light"]); sand  <- unname(ewas_col["sand"])
plum   <- unname(ewas_col["plum"]);       gry_d <- unname(ewas_col["gray_dark"])
gry    <- unname(ewas_col["gray"]);       gry_l <- unname(ewas_col["gray_light"])

ord <- fread("../yame_data/EPIC/EPIC.ordering.tsv.gz", select = "Probe_ID")
crd <- fread("../yame_data/EPIC/EPIC.hg38.coord.tsv.gz")
map <- data.table(probe = ord$Probe_ID, CpG_chrm = crd$CpG_chrm, CpG_beg = crd$CpG_beg)

e   <- readRDS("data/06_ewas.rds")
ann <- merge(as.data.table(e$tt), map, by = "probe")
N_TESTED <- nrow(as.data.table(e$tt))

reg <- fread(cmd = sprintf("zcat %s", file.path(outd, "dmr/PTSD_dmr.regions-p.bed.gz")))
setnames(reg, c("chrom", "start", "end", "min_p", "n_probes", "z_p", "z_sidak_p"))
reg <- reg[z_sidak_p < 0.05]
stopifnot(nrow(reg) == 1)

PAD <- 260
loc <- ann[CpG_chrm == reg$chrom[1] &
           CpG_beg >= reg$start[1] - PAD & CpG_beg <= reg$end[1] + PAD,
           .(probe, pos = CpG_beg, p = bacon.p, db = delta_beta)]
loc[, inreg := pos >= reg$start[1] & pos <= reg$end[1]]
stopifnot(sum(loc$inreg) == reg$n_probes[1])
sidak <- reg$z_sidak_p[1]

p08e <- ggplot(loc, aes(pos, -log10(p))) +
  annotate("rect", xmin = reg$start[1], xmax = reg$end[1], ymin = -Inf, ymax = Inf,
           fill = teal_l, alpha = 0.16) +
  geom_hline(yintercept = -log10(sidak), colour = plum, linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05 / N_TESTED), colour = gry,
             linetype = 2, linewidth = 0.4) +
  geom_segment(aes(xend = pos, yend = 0, colour = inreg), linewidth = 0.6) +
  geom_point(aes(colour = inreg, shape = db < 0), size = 2.4) +
  annotate("text", x = min(loc$pos), y = -log10(sidak), hjust = 0, vjust = -0.55,
           size = 3, colour = plum,
           label = sprintf("region \u0160id\u00e1k p = %.2e", sidak),
           family = if (nzchar(.ewas_family)) .ewas_family else NULL) +
  annotate("text", x = min(loc$pos), y = -log10(0.05 / N_TESTED), hjust = 0,
           vjust = 1.5, size = 3, colour = gry,
           label = "single-CpG genome-wide threshold",
           family = if (nzchar(.ewas_family)) .ewas_family else NULL) +
  scale_colour_manual(values = c(`TRUE` = teal_d, `FALSE` = gry_l), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 25, `FALSE` = 24),
                     labels = c(`TRUE` = "Hypomethylated in cases",
                                `FALSE` = "Hypermethylated in cases"), name = NULL) +
  scale_y_continuous(limits = c(0, -log10(sidak) * 1.12),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(breaks = scales::breaks_width(100),
                     expand = expansion(mult = c(0.05, 0.05)),
                     labels = function(x) format(x, big.mark = ",", trim = TRUE,
                                                 scientific = FALSE)) +
  labs(x = sprintf("Position on %s (hg38)", reg$chrom[1]),
       y = expression(-log[10](italic(p))),
       title = sprintf("%d adjacent CpGs make a region no single CpG could make",
                       reg$n_probes[1]),
       subtitle = sprintf("%s:%s\u2013%s; the shaded %d bp window is what comb-p reports.\nThe strongest single CpG in it reaches only \u2212log\u2081\u2080 p = %.2f",
                          reg$chrom[1], format(reg$start[1], big.mark = ","),
                          format(reg$end[1], big.mark = ","),
                          reg$end[1] - reg$start[1],
                          max(-log10(loc[inreg == TRUE]$p)))) +
  theme_ewas(base_size = 11) +
  theme(legend.position = "top",
        plot.margin = margin(5.5, 22, 5.5, 5.5))
save_fig(p08e, file.path(outd, "08_dmr_locus.png"), 9.0, 5.2)

cat("n_tested:", N_TESTED, "\n")
cat("n_in_window(plotted):", nrow(loc), " in_region:", sum(loc$inreg), "\n")
cat("sidak:", signif(sidak,3), " span:", reg$end[1]-reg$start[1], "\n")
cat("max_neglog10_in_region:", round(max(-log10(loc[inreg==TRUE]$p)),2), "\n")
cat("bonf_threshold_neglog10:", round(-log10(0.05/N_TESTED),2), "\n")
