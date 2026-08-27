## Per-SV variance explained, for the ADOPTED design:
##   ComBat(slide) within sex strata, position a fixed covariate, smoke in mod & mod0
suppressPackageStartupMessages({library(sva); library(limma)})
setwd("repo"); source("_setup.R"); setwd("..")
set.seed(42)
sv  <- readRDS("repo/data/05_sva.rds")
Mcb <- readRDS("repo/data/05_mvals_combat.rds"); if (is.list(Mcb)) Mcb <- Mcb[[1]]
md <- sv$mdk; smoke <- sv$smoke; P <- as.data.frame(sv$propk)
ptsd <- relevel(factor(md$ptsd), ref="Control"); sex <- factor(md$sex); age <- as.numeric(md$age)
slide <- factor(as.character(md$slide)); levels(slide) <- paste0("s", seq_len(nlevels(slide)))
pos <- factor(substr(as.character(md$array_pos), 1, 3))
cells <- c("Neu","NK","CD4T","CD8T","Bcell","Mono")
base <- data.frame(ptsd, sex, age, smoke, P[, cells])
modC <- model.matrix(~ . + pos, base); mod0C <- model.matrix(~ . - ptsd + pos, base)

s15 <- sva(Mcb, modC, mod0C, n.sv = 15)
SV <- s15$sv
sub <- order(apply(Mcb, 1, var), decreasing = TRUE)[1:20000]
Ms  <- Mcb[sub, ]

## residualise on the FULL measured model, so "variance explained" means
## variance the SVs explain over and above ptsd/sex/age/smoke/cells/position
R <- residuals(lmFit(Ms, modC), Ms)
tot <- sum(R^2)

## marginal: each SV alone against the residual
marg <- sapply(seq_len(ncol(SV)), function(j) {
  H <- SV[, j, drop = FALSE]
  fitted <- R %*% H %*% solve(crossprod(H)) %*% t(H)
  sum(fitted^2) / tot
})
## sequential (orthogonalised, in SVA's native order): incremental contribution
seq_pve <- numeric(ncol(SV))
for (j in seq_len(ncol(SV))) {
  H <- SV[, 1:j, drop = FALSE]
  fitted <- R %*% H %*% solve(crossprod(H)) %*% t(H)
  seq_pve[j] <- sum(fitted^2) / tot - (if (j > 1) sum(seq_pve[1:(j-1)]) else 0)
}
r2 <- function(y, X) summary(lm(y ~ X))$r.squared
chip <- sapply(seq_len(ncol(SV)), function(j) r2(SV[, j], slide))
posr <- sapply(seq_len(ncol(SV)), function(j) r2(SV[, j], pos))
smk  <- sapply(seq_len(ncol(SV)), function(j) summary(lm(smoke ~ SV[, j]))$r.squared)
ptr  <- sapply(seq_len(ncol(SV)), function(j) r2(SV[, j], ptsd))

out <- data.frame(SV = seq_len(ncol(SV)),
                  pve_marginal = round(100*marg, 3),
                  pve_sequential = round(100*seq_pve, 3),
                  r2_slide = round(chip, 3), r2_position = round(posr, 3),
                  r2_smoke = round(smk, 3), r2_ptsd = round(ptr, 3))
cat("=== per-SV variance explained (% of residual variance after the measured model)\n")
cat("=== SVA's native order (NOT variance-ordered):\n")
print(out, row.names = FALSE)
cat(sprintf("\nSV block total (all 15 jointly): %.2f%% of residual variance\n",
    100 * { H <- SV; sum((R %*% H %*% solve(crossprod(H)) %*% t(H))^2)/tot }))
o <- order(out$pve_marginal, decreasing = TRUE)
cat("\n=== ranked by marginal pve, with cumulative:\n")
cum <- cumsum(out$pve_marginal[o])
for (i in seq_along(o)) cat(sprintf("  rank %2d: SV%-2d  pve %.3f%%  cum %.2f%%  (slide %.2f pos %.2f smoke %.2f)\n",
    i, out$SV[o][i], out$pve_marginal[o][i], cum[i],
    out$r2_slide[o][i], out$r2_position[o][i], out$r2_smoke[o][i]))
write.csv(out, "repo/data/05_sv_pve.csv", row.names = FALSE)
saveRDS(list(table = out, order = o, SV = SV), "repo/data/05_sv_pve.rds")
cat("\nsaved 05_sv_pve.csv\n")
