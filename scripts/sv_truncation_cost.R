## What does truncating the SV set actually cost?
## Evaluate every k under three orderings, tracking BOTH nuisance axes.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
setwd("repo"); source("_setup.R")

inp <- readRDS(data_path("05_sv_selection_inputs.rds"))
d <- copy(inp$ranking)

tot_tech <- sum(d$tech); tot_smk <- sum(d$smoking)

## three candidate orderings of the SVs
ords <- list(
  `SV index`        = d[order(k)]$sv,
  `technical R2`    = d[order(-tech)]$sv,
  `combined rank`   = d[order(-(tech / max(tech) + smoking / max(smoking)))]$sv
)

## NB: do the lookup with plain vectors, not d[...]. Inside a data.table `[`,
## a loop variable named `k` is shadowed by d's own `k` column.
res <- rbindlist(lapply(names(ords), function(nm) {
  o <- ords[[nm]]
  tech_v <- d$tech[match(o, d$sv)]
  smk_v  <- d$smoking[match(o, d$sv)]
  data.table(ordering = nm, k = seq_along(o),
             tech_pct = 100 * cumsum(tech_v) / tot_tech,
             smk_pct  = 100 * cumsum(smk_v)  / tot_smk)
}))

cat("smoking signal retained (% of total across all 15 SVs):\n")
print(dcast(res[k %in% c(4,6,8,10,12,15)], k ~ ordering, value.var = "smk_pct")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)])
cat("\ntechnical signal retained:\n")
print(dcast(res[k %in% c(4,6,8,10,12,15)], k ~ ordering, value.var = "tech_pct")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)])

## the cost of k=8 under a technical ordering, spelled out
k8 <- d[order(-tech)][1:8]
drop8 <- d[order(-tech)][9:15]
cat("\nk=8 by technical rank drops:", paste(drop8$sv, collapse = ", "), "\n")
cat("  those dropped SVs carry", round(100*sum(drop8$smoking)/tot_smk, 1),
    "% of the total smoking signal\n")
cat("  including SV3 (smoking R2", round(d[sv=="SV3"]$smoking,3),
    ") and SV2 (", round(d[sv=="SV2"]$smoking,3), ")\n")

## residual df at each k
dfs <- data.table(k = 0:15)[, `:=`(Overall = 87 - (10 + k),
                                   Female = 45 - (9 + k),
                                   Male   = 42 - (9 + k))]
cat("\nresidual df by k:\n"); print(dfs[k %in% c(0,4,6,8,10,12,15)])

## ---- figure -----------------------------------------------------------------
long <- melt(res, id.vars = c("ordering","k"),
             measure.vars = c("tech_pct","smk_pct"),
             variable.name = "axis", value.name = "pct")
long[, axis := factor(axis, levels = c("tech_pct","smk_pct"),
                      labels = c("technical (chip / array position)","smoking proxy"))]

pa <- ggplot(long, aes(k, pct, colour = ordering)) +
  geom_hline(yintercept = 80, linetype = 3, colour = "grey50") +
  geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
  facet_wrap(~ axis) +
  scale_colour_manual(values = c(`SV index` = "grey55",
                                 `technical R2` = "#C44E52",
                                 `combined rank` = "#4C72B0")) +
  scale_x_continuous(breaks = seq(0, 15, 3)) +
  labs(x = "number of SVs retained (k)", y = "% of nuisance signal retained",
       colour = "SVs added in order of",
       title = "A. Truncating the SV set trades away smoking adjustment faster than technical adjustment",
       subtitle = "Ranking on technical association alone (red, right panel) discards the SVs carrying smoking signal") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")

pb <- ggplot(d, aes(tech, smoking)) +
  geom_point(aes(size = pve), shape = 21, fill = "#4C72B0", alpha = 0.75, stroke = 0.3) +
  ggrepel::geom_text_repel(aes(label = sv), size = 2.9, seed = 42, max.overlaps = 20) +
  scale_size_continuous(range = c(2, 7), name = "% resid. var") +
  labs(x = expression(R^2~"with chip / array position"),
       y = expression(R^2~"with smoking proxy"),
       title = "B. The two nuisance axes are carried by different SVs",
       subtitle = "SV4 is almost purely technical; SV5 and SV3 lead on smoking. No SV captures both.") +
  theme_minimal(base_size = 10)

p <- pa / pb + plot_layout(heights = c(1, 1.1))
ggsave(data_path("05_sv_truncation.png"), plot = p, width = 9.2, height = 8.2, dpi = 200)
fwrite(res, data_path("05_sv_truncation_curves.csv"))
