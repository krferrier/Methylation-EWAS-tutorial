## Rebuild the two v8.1 render-tree checkpoints that were lost:
##   1. repo/data/03_grs_filtered.rds  -> 756,273 probes x 96 samples
##   2. repo/data/05_mvals_combat.rds  -> stratified-ComBat M-value matrix
##
## Every parameter is taken from the frozen checkpoints already on disk
## (03_filter_funnel.rds, 05_sva.rds, 05_smoking_proxy.rds) rather than
## recomputed, so no published number can move. 05_sva.rds is READ ONLY here:
## the SVs, k=6, mdk, propk and the smoking proxy are reused exactly as saved.
suppressPackageStartupMessages({
  library(minfi); library(sva)
})
setwd("repo"); source("_setup.R")
set.seed(42)
t0 <- Sys.time()

funl <- readRDS(data_path("03_filter_funnel.rds"))
sva_o <- readRDS(data_path("05_sva.rds"))
prox <- readRDS(data_path("05_smoking_proxy.rds"))

## ---- 1. refilter the normalized set to the frozen v8.1 probe list -----------
g <- readRDS(data_path("02_funnorm_grs.rds"))
cat("02_funnorm_grs:", nrow(g), "x", ncol(g), "\n")
stopifnot(all(funl$retained_probes %in% rownames(g)))
grs_filt <- g[funl$retained_probes, ]
rm(g); invisible(gc())
cat("03_grs_filtered:", nrow(grs_filt), "x", ncol(grs_filt),
    "| matches funnel retained:", nrow(grs_filt) == funl$retained, "\n")
saveRDS(grs_filt, data_path("03_grs_filtered.rds"))

## ---- 2. stratified ComBat, exactly as recompute_05_v3.R ---------------------
keep <- sva_o$keep; mdk <- sva_o$mdk; propk <- sva_o$propk
smoke <- sva_o$smoke                      # frozen proxy, not recomputed
pv2 <- prox$present                       # 20 panel CpGs excluded from testing

beta <- minfi::getBeta(grs_filt)[, keep, drop = FALSE]
rm(grs_filt); invisible(gc())
cat("beta:", paste(dim(beta), collapse = " x "), "\n")

M_all <- log2(beta / (1 - beta))
excl <- rownames(M_all) %in% pv2
M <- M_all[!excl, , drop = FALSE]
M <- M[rowSums(!is.finite(M)) == 0, , drop = FALSE]
rm(M_all, beta); invisible(gc())
cat("excluded panel probes:", sum(excl), "| tested probes:", nrow(M),
    "| matches 05_sva tested_probes:", nrow(M) == sva_o$tested_probes, "\n")

ptsd <- factor(mdk$ptsd); sex <- factor(mdk$sex); age <- as.numeric(mdk$age)
slide <- factor(as.character(mdk$slide))
levels(slide) <- paste0("s", seq_len(nlevels(slide)))
P <- as.data.frame(propk)
cell <- c("Neu", "NK", "CD4T", "CD8T", "Bcell", "Mono")

cat("\n--- stratified ComBat on slide (mean.only = TRUE) ---\n")
Mcb <- M
for (s in levels(sex)) {
  i <- which(sex == s)
  bs <- droplevels(slide[i])
  d <- cbind(P[i, cell, drop = FALSE],
             ptsd = ptsd[i], age = age[i], smoke = smoke[i])
  cb_mod <- model.matrix(~ ptsd + age + smoke + Neu + NK + CD4T + CD8T + Bcell + Mono,
                         data = d)
  tb <- table(bs)
  cat(sprintf("  %s: n = %d | slides = %d | min batch = %d | cb_mod cols = %d\n",
              s, length(i), nlevels(bs), min(tb), ncol(cb_mod)))
  Mcb[, i] <- ComBat(dat = M[, i, drop = FALSE], batch = bs, mod = cb_mod,
                     par.prior = TRUE, mean.only = TRUE, prior.plots = FALSE)
}
rm(M); invisible(gc())
saveRDS(Mcb, data_path("05_mvals_combat.rds"))
cat("\n05_mvals_combat:", paste(dim(Mcb), collapse = " x "), "\n")
cat("elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
