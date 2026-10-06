## catalog_overview_fig.R -- draws tutorial/data/00b_catalog_overview.png, the
## four-panel overview of the 70 GEO Series screened in the dataset catalogue
## chapter (00b_dataset_catalog.qmd).
##
## The figure describes all 70 Series as screened, before the chapter groups
## Series from the same publication into 60 studies, so it reads the raw screen
## (00b_geo_catalog_source.csv) rather than the published 60-row catalogue.
##
## Run from the repository root:   Rscript scripts/catalog_overview_fig.R
## Needs ggplot2, patchwork and scales. Colours, font and theme come from
## tutorial/_setup.R, the same file every chapter sources.

suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(patchwork); library(scales)
})
source("tutorial/_setup.R")

teal_d <- unname(ewas_col["teal_dark"]);  teal  <- unname(ewas_col["teal"])
teal_l <- unname(ewas_col["teal_light"]); sand  <- unname(ewas_col["sand"])
plum   <- unname(ewas_col["plum"]);       gry_d <- unname(ewas_col["gray_dark"])
gry_l  <- unname(ewas_col["gray_light"])

cat0 <- as.data.table(read.csv("tutorial/00b_geo_catalog_source.csv", check.names = FALSE))
stopifnot(nrow(cat0) == 70)

## a: teaching-suitability tier
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

## b: metadata fields present; the two a first EWAS most wants and least often
## gets are accented
fields <- data.table(
  field = c("Sex", "Age", "Tissue / cell", "Phenotype or exposure",
            "Ancestry / race", "Cell-type composition"),
  n = c(sum(cat0$Has_sex == "True"), sum(cat0$Has_age == "True"),
        sum(cat0$Has_tissue == "True"), sum(cat0$N_phenotype_fields > 0),
        sum(cat0$Has_ancestry == "True"), sum(cat0$Has_cell_composition == "True")))
fields[, field := factor(field, levels = rev(field))]
fields[, key := n < 20]
pb <- ggplot(fields, aes(n, field, fill = key)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = n), hjust = -0.28, size = 3.2, colour = gry_d,
            family = .ewas_family) +
  scale_fill_manual(values = c(`TRUE` = sand, `FALSE` = teal), guide = "none") +
  scale_x_continuous(limits = c(0, 70), expand = expansion(mult = c(0, 0.10))) +
  labs(title = "Metadata fields present", x = "Studies (of 70)", y = NULL)

## c: study size
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

## d: array generation
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

ggsave("tutorial/data/00b_catalog_overview.png", plot = p00b, width = 11, height = 7.6,
       dpi = 300, device = grDevices::png, type = "cairo", bg = "white")
