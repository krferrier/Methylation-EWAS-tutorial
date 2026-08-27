# Final chapter 05 SVA under the adopted design:
#   ComBat(slide) within sex strata  ->  position as fixed covariate  ->  SVA, k = 6
suppressPackageStartupMessages({library(sva); library(limma)})
set.seed(42)
t0 <- Sys.time()
setwd("repo"); source("_setup.R")

K <- 6
Mcb <- readRDS("data/05_mvals_combat.rds")     # slide-corrected within strata
prev <- readRDS("data/05_sva.rds")
mdk <- prev$mdk; propk <- prev$propk; smoke <- prev$smoke
stopifnot(identical(colnames(Mcb), as.character(mdk$sample_id)))

pos  <- factor(as.character(mdk$array_pos))
ptsd <- relevel(factor(mdk$ptsd), ref = "Control")
sex  <- factor(mdk$sex)

mod  <- model.matrix(~ ptsd + sex + mdk$age + smoke + pos + propk)
mod0 <- model.matrix(~        sex + mdk$age + smoke + pos + propk)
colnames(mod) <- make.names(colnames(mod))
cat("mod cols:", ncol(mod), " rank:", qr(mod)$rank, "\n")

sv <- sva(Mcb, mod, mod0, n.sv = K)
SV <- sv$sv; colnames(SV) <- paste0("SV", seq_len(K))
cat("n.sv used:", ncol(SV), "\n")

# what each retained SV tracks
slide_f <- factor(as.character(mdk$slide))
r2 <- function(y, x) summary(lm(y ~ x))$r.squared
tab <- data.frame(
  SV       = seq_len(K),
  r2_slide = sapply(seq_len(K), function(i) r2(SV[, i], slide_f)),
  r2_pos   = sapply(seq_len(K), function(i) r2(SV[, i], pos)),
  r2_smoke = sapply(seq_len(K), function(i) r2(SV[, i], smoke)),
  r2_ptsd  = sapply(seq_len(K), function(i) r2(SV[, i], ptsd)),
  r2_age   = sapply(seq_len(K), function(i) r2(SV[, i], mdk$age)),
  r2_sex   = sapply(seq_len(K), function(i) r2(SV[, i], sex)))
print(round(tab, 4))
cat("R2(smoke ~ SVs):", round(summary(lm(smoke ~ SV))$r.squared, 4), "\n")
cat("max R2(SV ~ PTSD):", round(max(tab$r2_ptsd), 4), "\n")

# df accounting under the final design
npar <- 1 + 1 + 1 + 1 + 1 + 7 + 6 + K   # int, ptsd, sex, age, smoke, 7 pos, 6 cells, K SVs
cat(sprintf("params combined: %d  resid df: %d (n=87)\n", npar, 87 - npar))
cat(sprintf("params stratum : %d  resid df F: %d (n=45)  M: %d (n=42)\n",
            npar - 1, 45 - (npar - 1), 42 - (npar - 1)))

out <- prev
out$SV <- SV; out$n.sv <- K; out$sv_table <- tab
out$k_selected <- K; out$k_candidates <- c(6, 8, 10, 15)
out$r2_smoke_sv <- summary(lm(smoke ~ SV))$r.squared
out$position_in_model <- TRUE
out$design_note <- "ComBat(slide) within sex strata; array position as fixed covariate; SVA k=6"
saveRDS(out, "data/05_sva.rds")
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")
