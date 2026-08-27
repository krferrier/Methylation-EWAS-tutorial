## Chapter 05 rebuild, v2.
##
## Two changes driven by user decisions:
##   1. The methylation-derived smoking proxy enters the SVA null model, so
##      smoking is an EXPLICITLY MODELLED covariate and the SVs estimate only
##      the variation that remains unmeasured after it.
##   2. The CpGs used to build the proxy are EXCLUDED from the tested matrix,
##      removing the circularity of testing probes that define a covariate.
##
## Also quantifies whether ComBat on slide reduces the number of SVs needed,
## and whether ComBat is even admissible given slide/sex confounding.
suppressPackageStartupMessages({
  library(minfi); library(data.table); library(ggplot2); library(sva)
})
setwd("repo"); source("_setup.R")
set.seed(42)

g     <- readRDS(data_path("03_grs_filtered.rds"))
sva_o <- readRDS(data_path("05_sva.rds"))
keep  <- sva_o$keep; mdk <- sva_o$mdk; propk <- sva_o$propk
beta  <- minfi::getBeta(g)[, keep, drop = FALSE]
cat("beta:", paste(dim(beta), collapse = " x "), "\n")

## ---- 1. smoking proxy (unchanged construction) ------------------------------
panel <- list(
  AHRR      = c("cg05575921","cg21161138","cg23576855","cg25648203","cg26703534","cg11902777"),
  F2RL3     = c("cg03636183","cg21911711"),
  GPR15     = c("cg19859270"),
  `2q37.1`  = c("cg05951221","cg21566642","cg01940273","cg03329539","cg06126421"),
  GFI1      = c("cg09935388","cg12876356","cg18146737","cg06338710"),
  MYO1G     = c("cg12803068","cg22132788","cg04180046","cg19089201"),
  PRSS23    = c("cg14391737"),
  RARA      = c("cg17739917")
)
pv      <- unlist(panel, use.names = FALSE)
present <- pv %in% rownames(beta)
pv2     <- pv[present]
cat("panel probes:", length(pv), "| present:", sum(present), "\n")

Bp <- beta[pv2, , drop = FALSE]
Z  <- t(scale(t(Bp))); Z <- Z[rowSums(is.na(Z)) == 0, , drop = FALSE]
pc <- prcomp(t(Z), center = TRUE, scale. = FALSE)
smoke <- pc$x[, 1]
pve_panel <- pc$sdev^2 / sum(pc$sdev^2)
if (cor(smoke, Bp["cg05575921", ], use = "complete.obs") > 0) smoke <- -smoke
smoke <- as.numeric(scale(smoke))
cat("panel PC1 pve:", round(100 * pve_panel[1], 1), "% | anchor cor:",
    round(cor(smoke, Bp["cg05575921", ], use = "complete.obs"), 3), "\n")

## ---- 2. exclude panel probes from the tested matrix -------------------------
M_all <- log2(beta / (1 - beta))
excl  <- rownames(M_all) %in% pv2
M     <- M_all[!excl, , drop = FALSE]
ok    <- rowSums(!is.finite(M)) == 0
M     <- M[ok, , drop = FALSE]
cat("excluded panel probes:", sum(excl), "| non-finite dropped:", sum(!ok),
    "| tested probes:", nrow(M), "\n")

## ---- 3. SVA with the proxy protected in both models ------------------------
ptsd <- factor(mdk$ptsd); sex <- factor(mdk$sex); age <- as.numeric(mdk$age)
P    <- as.data.frame(propk)
mod  <- model.matrix(~ ptsd + sex + age + smoke + Neu + NK + CD4T + CD8T + Bcell + Mono,
                     data = cbind(P, ptsd = ptsd, sex = sex, age = age, smoke = smoke))
mod0 <- model.matrix(~ sex + age + smoke + Neu + NK + CD4T + CD8T + Bcell + Mono,
                     data = cbind(P, sex = sex, age = age, smoke = smoke))
cat("mod cols:", ncol(mod), "| mod0 cols:", ncol(mod0), "\n")

nsv_new <- num.sv(M, mod, method = "be")
cat("num.sv (be) WITH smoking proxy in the model:", nsv_new, "\n")
cat("num.sv previously (no proxy, 15):", ncol(sva_o$SV), "\n")

svo <- sva(M, mod, mod0, n.sv = nsv_new)
SV  <- svo$sv; colnames(SV) <- paste0("SV", seq_len(ncol(SV)))
cat("SV dims:", paste(dim(SV), collapse = " x "), "\n")

## how much of the proxy do the new SVs still explain?
r2 <- function(y, X) summary(lm(y ~ X))$r.squared
cat("R2(smoke ~ new SVs):", round(r2(smoke, SV), 3), "\n")
slide <- factor(as.character(mdk$slide)); pos <- factor(substr(mdk$array_pos, 1, 3))
chip_r2 <- apply(SV, 2, function(s) r2(s, slide))
cat("chip R2 per SV:", paste(sprintf("SV%d=%.2f", seq_along(chip_r2), chip_r2), collapse = " "), "\n")
cat("max PTSD R2 across SVs:",
    round(max(apply(SV, 2, function(s) r2(s, model.matrix(~ptsd)[, -1]))), 4), "\n")

## ---- 4. is ComBat admissible here? -----------------------------------------
cat("\n--- ComBat feasibility ---\n")
cat("R2(sex ~ slide):", round(summary(lm(as.numeric(sex) ~ slide))$r.squared, 3),
    " R2(PTSD ~ slide):", round(summary(lm(as.numeric(ptsd) ~ slide))$r.squared, 3), "\n")
cb_mod <- model.matrix(~ ptsd + sex + age + smoke + Neu + NK + CD4T + CD8T + Bcell + Mono,
                       data = cbind(P, ptsd = ptsd, sex = sex, age = age, smoke = smoke))
combat_ok <- TRUE
Mcb <- tryCatch(
  ComBat(dat = M, batch = slide, mod = cb_mod, par.prior = TRUE, prior.plots = FALSE),
  error = function(e) { combat_ok <<- FALSE; cat("ComBat FAILED (sex-protected):", conditionMessage(e), "\n"); NULL })

## fallback: protect everything except sex
if (is.null(Mcb)) {
  cb_mod2 <- model.matrix(~ ptsd + age + smoke + Neu + NK + CD4T + CD8T + Bcell + Mono,
                          data = cbind(P, ptsd = ptsd, age = age, smoke = smoke))
  Mcb <- tryCatch(
    ComBat(dat = M, batch = slide, mod = cb_mod2, par.prior = TRUE, prior.plots = FALSE),
    error = function(e) { cat("ComBat FAILED (sex unprotected too):", conditionMessage(e), "\n"); NULL })
  if (!is.null(Mcb)) cat("ComBat ran only WITHOUT sex protected.\n")
}

if (!is.null(Mcb)) {
  nsv_cb <- num.sv(Mcb, mod, method = "be")
  cat("num.sv after ComBat on slide:", nsv_cb, " (vs", nsv_new, "uncorrected)\n")
  sex_r2_before <- mean(apply(M[1:20000, ], 1, function(y) summary(lm(y ~ sex))$r.squared))
  sex_r2_after  <- mean(apply(Mcb[1:20000, ], 1, function(y) summary(lm(y ~ sex))$r.squared))
  cat("mean sex R2 over 20k probes  before:", round(sex_r2_before, 4),
      " after ComBat:", round(sex_r2_after, 4), "\n")
  saveRDS(list(n_sv = nsv_cb, sex_r2_before = sex_r2_before, sex_r2_after = sex_r2_after,
               ran_with_sex = combat_ok),
          data_path("05_combat_diag.rds"))
}

## ---- 5. save ---------------------------------------------------------------
saveRDS(list(SV = SV, n.sv = nsv_new, keep = keep, mdk = mdk, propk = propk,
             smoke = smoke, panel_probes = pv2, panel_pve = pve_panel,
             excluded_n = sum(excl), tested_probes = nrow(M),
             prev_n_sv = ncol(sva_o$SV), r2_smoke_sv = r2(smoke, SV),
             chip_r2 = chip_r2),
        data_path("05_sva.rds"))
saveRDS(list(smoke = smoke, panel = panel, present = pv2, pve = pve_panel,
             loadings = pc$rotation[, 1]),
        data_path("05_smoking_proxy.rds"))
cat("\nsaved 05_sva.rds and 05_smoking_proxy.rds\n")
