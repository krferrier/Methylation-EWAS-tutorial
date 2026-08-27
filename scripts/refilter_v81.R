## Refilter the normalized GenomicRatioSet under the Zhou v8.1 eight-code mask.
## Mask read from the bit-packed EPIC.hg38.mask.cm via YAME (see DECISIONS.md).
suppressPackageStartupMessages({
  library(minfi); library(data.table)
})
setwd("repo")
source("_setup.R")

t0 <- Sys.time()
grs  <- readRDS(data_path("02_funnorm_grs.rds"))
detP <- readRDS("../bundle_extract/data/01_detP.rds")
cat("grs:", nrow(grs), "probes x", ncol(grs), "samples\n")
cat("detP:", nrow(detP), "x", ncol(detP), "\n")

## ---- v8.1 mask: the eight codes, union ----
terms  <- readLines("../term_order.txt")
probes <- readLines("../probe_order.txt")
want <- c("M_mapping","M_nonuniq",
          "M_SNP_AFR_1pt","M_1baseSwitchSNP_AFR_1pt","M_2extBase_SNP_AFR_1pt",
          "M_SNP_EUR_1pt","M_1baseSwitchSNP_EUR_1pt","M_2extBase_SNP_EUR_1pt")
bits <- fread("../all24.txt", header = FALSE)
setnames(bits, terms)
stopifnot(nrow(bits) == length(probes))

per_code <- sapply(want, function(w) sum(bits[[w]] == 1L))
mask_any <- Reduce(`|`, lapply(want, function(w) bits[[w]] == 1L))
drop8    <- probes[mask_any]
cat("union of 8 codes:", length(drop8), "\n")

## also keep the general mask and a EUR-only comparator for the qmd figures
m_general <- probes[bits$M_general == 1L]
eur_only  <- c("M_mapping","M_nonuniq","M_SNP_EUR_1pt",
               "M_1baseSwitchSNP_EUR_1pt","M_2extBase_SNP_EUR_1pt")
drop_eur  <- probes[Reduce(`|`, lapply(eur_only, function(w) bits[[w]] == 1L))]
snp_by_pop <- sapply(c("AFR","EUR","EAS","SAS","AMR"),
                     function(p) sum(bits[[paste0("M_SNP_", p, "_1pt")]] == 1L))

## ---- the funnel ----
rn <- rownames(grs)
n0 <- length(rn)

drop_mask <- intersect(rn, drop8)
frac_fail <- rowMeans(detP[intersect(rn, rownames(detP)), , drop = FALSE] > 0.01)
drop_detp <- names(frac_fail)[frac_fail > 0.10]
ann       <- getAnnotation(grs)
drop_sex  <- rn[ann$chr %in% c("chrX", "chrY")]

keep <- setdiff(rn, Reduce(union, list(drop_mask, drop_detp, drop_sex)))
grs_filt <- grs[keep, ]

## EUR-only comparator (mask step only, same detp/sex rules)
keep_eur <- setdiff(rn, Reduce(union, list(intersect(rn, drop_eur), drop_detp, drop_sex)))

funl <- list(
  n0 = n0,
  drop_mask       = length(drop_mask),
  drop_detp       = length(drop_detp),
  drop_sex        = length(drop_sex),
  retained        = length(keep),
  retained_eur    = length(keep_eur),
  afr_not_eur_kept = length(setdiff(keep_eur, keep)),
  per_code        = per_code,
  union8          = length(drop8),
  m_general_n     = length(m_general),
  snp_by_pop      = snp_by_pop,
  mask_version    = "Zhou EPIC.hg38.mask.cm v8.1 (via YAME v1.40)",
  retained_probes = keep
)
mp <- list(
  probe_order   = probes,
  drop8         = drop8,
  per_code_sets = setNames(lapply(want, function(w) probes[bits[[w]] == 1L]), want),
  general_mask  = m_general,
  drop_eur      = drop_eur,
  snp_by_pop    = snp_by_pop
)

saveRDS(grs_filt, data_path("03_grs_filtered.rds"))
saveRDS(funl,     data_path("03_filter_funnel.rds"))
saveRDS(mp,       data_path("03_mask_pieces.rds"))

cat("\n=== FUNNEL (v8.1) ===\n")
cat("n0                :", n0, "\n")
cat("drop mask (union8):", length(drop_mask), "\n")
cat("drop detp         :", length(drop_detp), "\n")
cat("drop sex          :", length(drop_sex), "\n")
cat("RETAINED          :", length(keep), "\n")
cat("retained (EUR-only mask):", length(keep_eur), "\n")
cat("kept by EUR but dropped by union8:", funl$afr_not_eur_kept, "\n")
cat("\nper-code:\n"); print(per_code)
cat("\nSNP by super-population:\n"); print(snp_by_pop)
cat("\nelapsed:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
