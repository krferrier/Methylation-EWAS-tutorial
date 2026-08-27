## Chapter 05 rebuild, v3: STRATIFIED ComBat -> recombine -> SVA.
##
## Design (user ruling): slide is confounded with sex across the whole cohort
## (R2 = 0.874), so a cohort-wide ComBat on slide must protect sex against a
## variable it is nearly collinear with. Running ComBat WITHIN each sex stratum
## dissolves the confounding -- inside a stratum, slide is just a batch -- and
## because ComBat preserves the mean of the data it is handed, between-sex
## differences pass through untouched. The corrected strata are then recombined
## and SVA estimates whatever technical/biological structure remains.
##
## mean.only = TRUE for BOTH strata: male slide s6 carries a single sample, and
## ComBat auto-forces mean.only when any batch has n=1. Setting it explicitly
## keeps the two strata under an identical correction instead of a silent
## asymmetry. Slide-to-slide differences here are treated as location shifts.
##
## The smoking proxy is in both SVA models (mod and mod0), so smoking is an
## explicitly modeled covariate rather than something the SVs are hoped to
## absorb; the panel CpGs are excluded from the tested matrix.
suppressPackageStartupMessages({
  library(minfi); library(sva); library(data.table)
})
setwd("repo"); source("_setup.R")
set.seed(42)

g     <- readRDS(data_path("03_grs_filtered.rds"))
sva_o <- readRDS(data_path("05_sva_prev.rds"))
keep  <- sva_o$keep; mdk <- sva_o$mdk; propk <- sva_o$propk
beta  <- minfi::getBeta(g)[, keep, drop = FALSE]
cat("beta:", paste(dim(beta), collapse = " x "), "\n")

## ---- 1. smoking proxy -------------------------------------------------------
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
pv  <- unlist(panel, use.names = FALSE)
pv2 <- pv[pv %in% rownames(beta)]
cat("panel probes:", length(pv), "| present:", length(pv2), "\n")

Bp <- beta[pv2, , drop = FALSE]
Z  <- t(scale(t(Bp))); Z <- Z[rowSums(is.na(Z)) == 0, , drop = FALSE]
pc <- prcomp(t(Z), center = TRUE, scale. = FALSE)
smoke <- pc$x[, 1]
pve_panel <- pc$sdev^2 / sum(pc$sdev^2)
if (cor(smoke, Bp["cg05575921", ], use = "complete.obs") > 0) smoke <- -smoke
smoke <- as.numeric(scale(smoke))
cat("panel PC1 pve:", round(100 * pve_panel[1], 1), "% | anchor cor:",
    round(cor(smoke, Bp["cg05575921", ], use = "complete.obs"), 3), "\n")

## ---- 2. tested matrix, panel probes excluded --------------------------------
M_all <- log2(beta / (1 - beta))
excl  <- rownames(M_all) %in% pv2
M     <- M_all[!excl, , drop = FALSE]
M     <- M[rowSums(!is.finite(M)) == 0, , drop = FALSE]
rm(M_all, beta); invisible(gc())
cat("excluded panel probes:", sum(excl), "| tested probes:", nrow(M), "\n")

ptsd <- factor(mdk$ptsd); sex <- factor(mdk$sex); age <- as.numeric(mdk$age)
slide <- factor(as.character(mdk$slide))
levels(slide) <- paste0("s", seq_len(nlevels(slide)))
pos <- factor(substr(mdk$array_pos, 1, 3))
P <- as.data.frame(propk)
cell <- c("Neu","NK","CD4T","CD8T","Bcell","Mono")

r2 <- function(y, X) summary(lm(y ~ X))$r.squared
sex_r2_mean <- function(X, n = 20000) {
  idx <- seq_len(min(n, nrow(X)))
  mean(apply(X[idx, , drop = FALSE], 1, function(y) summary(lm(y ~ sex))$r.squared))
}

cat("\n--- batch structure ---\n")
cat("R2(sex ~ slide) cohort-wide:", round(r2(as.numeric(sex), slide), 3),
    "| R2(PTSD ~ slide):", round(r2(as.numeric(ptsd), slide), 3), "\n")

## ---- 3. ComBat within each sex stratum -------------------------------------
cat("\n--- stratified ComBat on slide (mean.only = TRUE) ---\n")
Mcb <- M
for (s in levels(sex)) {
  i  <- which(sex == s)
  bs <- droplevels(slide[i])
  ## sex is constant inside a stratum, so it is not (and cannot be) in cb_mod
  d  <- cbind(P[i, cell, drop = FALSE],
              ptsd = ptsd[i], age = age[i], smoke = smoke[i])
  cb_mod <- model.matrix(~ ptsd + age + smoke + Neu + NK + CD4T + CD8T + Bcell + Mono,
                         data = d)
  tb <- table(bs)
  cat(sprintf("  %s: n = %d | slides = %d | min batch = %d | cb_mod cols = %d\n",
              s, length(i), nlevels(bs), min(tb), ncol(cb_mod)))
  Mcb[, i] <- ComBat(dat = M[, i, drop = FALSE], batch = bs, mod = cb_mod,
                     par.prior = TRUE, mean.only = TRUE, prior.plots = FALSE)
}

## ---- 4. did it work, and did sex survive? ----------------------------------
cat("\n--- effect of the correction ---\n")
sub <- seq_len(min(20000, nrow(M)))
slide_r2_before <- mean(apply(M[sub, ],   1, function(y) r2(y, slide)))
slide_r2_after  <- mean(apply(Mcb[sub, ], 1, function(y) r2(y, slide)))
pos_r2_before   <- mean(apply(M[sub, ],   1, function(y) r2(y, pos)))
pos_r2_after    <- mean(apply(Mcb[sub, ], 1, function(y) r2(y, pos)))
sxb <- sex_r2_mean(M); sxa <- sex_r2_mean(Mcb)
ptb <- mean(apply(M[sub, ], 1, function(y) r2(y, model.matrix(~ptsd)[, -1])))
pta <- mean(apply(Mcb[sub, ], 1, function(y) r2(y, model.matrix(~ptsd)[, -1])))
cat(sprintf("  mean R2 vs slide:    %.4f -> %.4f\n", slide_r2_before, slide_r2_after))
cat(sprintf("  mean R2 vs position: %.4f -> %.4f  (not corrected; left to SVA)\n",
            pos_r2_before, pos_r2_after))
cat(sprintf("  mean R2 vs sex:      %.4f -> %.4f  (must be preserved)\n", sxb, sxa))
cat(sprintf("  mean R2 vs PTSD:     %.4f -> %.4f\n", ptb, pta))

## ---- 5. SVA on the recombined, slide-corrected matrix ----------------------
d_all <- cbind(P[, cell, drop = FALSE], ptsd = ptsd, sex = sex, age = age, smoke = smoke)
mod  <- model.matrix(~ ptsd + sex + age + smoke + Neu + NK + CD4T + CD8T + Bcell + Mono, data = d_all)
mod0 <- model.matrix(~        sex + age + smoke + Neu + NK + CD4T + CD8T + Bcell + Mono, data = d_all)
cat("\n--- SVA ---\n")
cat("mod cols:", ncol(mod), "| rank:", qr(mod)$rank, "| mod0 cols:", ncol(mod0), "\n")

nsv_pre  <- num.sv(M,   mod, method = "be")
nsv_post <- num.sv(Mcb, mod, method = "be")
cat("num.sv BEFORE ComBat (proxy in model):", nsv_pre, "\n")
cat("num.sv AFTER  ComBat (proxy in model):", nsv_post, "\n")
cat("num.sv originally (no proxy, no ComBat):", ncol(sva_o$SV), "\n")

svo <- sva(Mcb, mod, mod0, n.sv = nsv_post)
SV  <- svo$sv; colnames(SV) <- paste0("SV", seq_len(ncol(SV)))
cat("SV dims:", paste(dim(SV), collapse = " x "), "\n")

chip_r2 <- apply(SV, 2, function(s) r2(s, slide))
pos_r2  <- apply(SV, 2, function(s) r2(s, pos))
cat("chip R2 per SV:", paste(sprintf("SV%d=%.2f", seq_along(chip_r2), chip_r2), collapse = " "), "\n")
cat("position R2 per SV:", paste(sprintf("SV%d=%.2f", seq_along(pos_r2), pos_r2), collapse = " "), "\n")
cat("R2(smoke ~ SVs):", round(r2(smoke, SV), 3), "\n")
cat("max R2(SV ~ PTSD):", round(max(apply(SV, 2, function(s) r2(s, model.matrix(~ptsd)[, -1]))), 4), "\n")

np <- 1 + 1 + 1 + 1 + 1 + length(cell) + ncol(SV)   # int + ptsd + sex + age + smoke + cells + SVs
cat("combined design params:", np, "| residual df:", 87 - np, "\n")
cat("female-arm params:", np - 1, "| residual df:", 45 - (np - 1), "\n")
cat("male-arm   params:", np - 1, "| residual df:", 42 - (np - 1), "\n")

## ---- 6. save ---------------------------------------------------------------
saveRDS(list(SV = SV, n.sv = nsv_post, keep = keep, mdk = mdk, propk = propk,
             smoke = smoke, panel_probes = pv2, panel_pve = pve_panel,
             excluded_n = sum(excl), tested_probes = nrow(Mcb),
             prev_n_sv = ncol(sva_o$SV), nsv_pre_combat = nsv_pre,
             r2_smoke_sv = r2(smoke, SV), chip_r2 = chip_r2, pos_r2 = pos_r2,
             combat = list(mean_only = TRUE, batch = "slide", within = "sex",
                           slide_r2 = c(before = slide_r2_before, after = slide_r2_after),
                           pos_r2   = c(before = pos_r2_before,   after = pos_r2_after),
                           sex_r2   = c(before = sxb, after = sxa),
                           ptsd_r2  = c(before = ptb, after = pta))),
        data_path("05_sva.rds"))
saveRDS(list(smoke = smoke, panel = panel, present = pv2, pve = pve_panel,
             loadings = pc$rotation[, 1]),
        data_path("05_smoking_proxy.rds"))
saveRDS(Mcb, data_path("05_mvals_combat.rds"))
cat("\nsaved 05_sva.rds, 05_smoking_proxy.rds, 05_mvals_combat.rds\n")
