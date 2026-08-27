## Two-panel diagnostic for choosing how many SVs to carry into the model.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
setwd("repo"); source("_setup.R")

d <- readRDS(data_path("05_sv_ranking.rds"))
setorder(d, -tech_max)
d[, rank := .I]
d[, lab := factor(sv, levels = sv)]
K <- 8L   # decided from the elbow in panel A; see prose

## Panel A: technical association, SVs in descending order
pa <- ggplot(d, aes(lab, tech_max)) +
  geom_col(aes(fill = rank <= K), width = 0.72) +
  geom_hline(yintercept = d[rank == K]$tech_max, linetype = 3, colour = "grey45") +
  geom_vline(xintercept = K + 0.5, linetype = 2, colour = "#C44E52") +
  annotate("text", x = K + 0.75, y = 0.80, hjust = 0, size = 3.1, colour = "#C44E52",
           label = sprintf("keep k = %d", K)) +
  scale_fill_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "grey78"), guide = "none") +
  labs(x = NULL, y = expression(R^2~"with chip / array position"),
       title = "A. Each SV against known technical factors",
       subtitle = "SVs ordered by technical association, not by index: SVA does not order its output by variance") +
  theme_minimal(base_size = 10) + theme(panel.grid.major.x = element_blank())

## Panel B: variance vs technical association — the two criteria disagree
pb <- ggplot(d, aes(tech_max, pve)) +
  geom_point(aes(fill = rank <= K), shape = 21, size = 3.4, stroke = 0.4) +
  ggrepel::geom_text_repel(aes(label = sv), size = 2.9, seed = 42, max.overlaps = 20) +
  scale_fill_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "grey78"), guide = "none") +
  labs(x = expression(R^2~"with chip / array position"), y = "% of residual variance",
       title = "B. Variance rank and technical rank disagree",
       subtitle = "SV2 carries the most variance but the least chip signal; SV4 the reverse") +
  theme_minimal(base_size = 10)

p <- pa / pb + plot_layout(heights = c(1, 1.05))
ggsave(data_path("05_sv_selection.png"), plot = p, width = 8.6, height = 7.6, dpi = 200)

cat("kept:", paste(d[rank <= K]$sv, collapse = ", "), "\n")
cat("dropped:", paste(d[rank > K]$sv, collapse = ", "), "\n")
cat("technical R2 captured by kept set:", round(d[rank == K]$cum_tech, 1), "%\n")
cat("max PTSD R2 among kept:", round(max(d[rank <= K]$PTSD), 4), "\n")
cat("resid df -> Overall:", 87 - (10 + K), "| Female:", 45 - (9 + K), "| Male:", 42 - (9 + K), "\n")
saveRDS(list(K = K, kept = d[rank <= K]$sv, ranking = d), data_path("05_sv_selection.rds"))
