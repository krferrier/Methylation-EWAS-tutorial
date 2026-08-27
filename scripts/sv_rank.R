## Rank SVs by what they are for: capturing known technical structure.
## SVA does not order SVs by variance, so a scree/80%-variance rule is not
## well-defined. Rank instead by association with measured technical factors,
## and check that no SV being kept or dropped is entangled with PTSD.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setwd("repo"); source("_setup.R")

sva <- readRDS(data_path("05_sva.rds"))
pca <- readRDS(data_path("05_batch_pca.rds"))
svv <- readRDS(data_path("05_sv_variance.rds"))
SV  <- sva$SV; mdk <- sva$mdk
cat("SV:", paste(dim(SV), collapse=" x "), "| mdk cols:", paste(names(mdk), collapse=", "), "\n")

r2 <- function(y, x) {
  d <- data.frame(y = y, x = x); d <- d[complete.cases(d), ]
  if (length(unique(d$x)) < 2) return(NA_real_)
  summary(lm(y ~ x, data = d))$r.squared
}

## slide is integer64 -> factor via character (a direct as.numeric corrupts it);
## array_pos "R01C01" -> row, matching chapter 05's construction
tech <- list(chip_slide     = factor(as.character(mdk$slide)),
             array_position = factor(substr(mdk$array_pos, 1, 3)))
bio  <- list(PTSD = mdk$ptsd, sex = mdk$sex, age = mdk$age)

d <- data.table(k = seq_len(ncol(SV)), sv = paste0("SV", seq_len(ncol(SV))))
d[, pve := 100 * svv$pve]
for (nm in names(tech)) d[[nm]] <- vapply(seq_len(ncol(SV)), function(i) r2(SV[, i], tech[[nm]]), 0)
for (nm in names(bio))  d[[nm]] <- vapply(seq_len(ncol(SV)), function(i) r2(SV[, i], bio[[nm]]), 0)
d[, tech_max := pmax(chip_slide, array_position, na.rm = TRUE)]
setorder(d, -tech_max)

print(d[, .(sv, pve = round(pve,2), chip = round(chip_slide,3),
            arraypos = round(array_position,3), tech_max = round(tech_max,3),
            PTSD = round(PTSD,3), sex = round(sex,3), age = round(age,3))])

cat("\nSVs with tech R2 >= 0.5:", paste(d[tech_max >= 0.5]$sv, collapse=", "), "\n")
cat("SVs with tech R2 >= 0.3:", paste(d[tech_max >= 0.3]$sv, collapse=", "), "\n")
cat("max PTSD R2 across all SVs:", round(max(d$PTSD, na.rm=TRUE), 4),
    "on", d[which.max(PTSD)]$sv, "\n")
cat("SVs with PTSD R2 > 0.05:", paste(d[PTSD > 0.05]$sv, collapse=", "), "\n")

## cumulative technical capture as SVs are added in tech-rank order
d[, rank := .I]
d[, cum_tech := cumsum(tech_max) / sum(tech_max) * 100]
cat("\ncumulative share of total technical R2 by rank:\n")
print(d[, .(rank, sv, tech_max = round(tech_max,3), cum_pct = round(cum_tech,1))])

fwrite(d, data_path("05_sv_ranking.csv"))
saveRDS(d, data_path("05_sv_ranking.rds"))
