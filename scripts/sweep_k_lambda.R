# Select k on calibration + precision, not residual batch R2.
# For each k: sva(n.sv = k) -> SVs in SVA's NATIVE order -> full-probe limma.
suppressPackageStartupMessages({library(sva); library(limma)})
set.seed(42)
t0 <- Sys.time()
setwd("repo"); source("_setup.R")

Mcb <- readRDS("data/05_mvals_combat.rds")
s   <- readRDS("data/05_sva.rds")
mdk <- s$mdk; propk <- s$propk; smoke <- s$smoke
stopifnot(identical(colnames(Mcb), as.character(mdk$sample_id)))

ptsd <- relevel(factor(mdk$ptsd), ref = "Control")
sex  <- factor(mdk$sex); pos <- factor(as.character(mdk$array_pos))
slide_f <- factor(as.character(mdk$slide))

mod  <- model.matrix(~ ptsd + sex + mdk$age + smoke + pos + propk)
mod0 <- model.matrix(~        sex + mdk$age + smoke + pos + propk)
colnames(mod) <- make.names(colnames(mod))
coefn <- "ptsdCase"; stopifnot(coefn %in% colnames(mod))

lam <- function(p) median(qchisq(1 - p, 1), na.rm = TRUE) / qchisq(0.5, 1)
v <- apply(Mcb, 1, var); idx <- order(v, decreasing = TRUE)[1:20000]

# baseline: no SVs
base <- local({
  fit <- eBayes(lmFit(Mcb, mod))
  tt  <- limma::topTable(fit, coef = coefn, number = Inf, sort.by = "none")
  list(tt = tt, se = sqrt(fit$s2.post) * fit$stdev.unscaled[, coefn], npar = ncol(mod))
})

res <- list(); svlist <- list()
for (K in c(6, 8, 10, 15)) {
  set.seed(42)
  sv <- sva(Mcb, mod, mod0, n.sv = K)
  SV <- sv$sv; colnames(SV) <- paste0("SV", seq_len(K))   # native order, first k
  dm <- cbind(mod, SV)
  fit <- eBayes(lmFit(Mcb, dm))
  tt  <- limma::topTable(fit, coef = coefn, number = Inf, sort.by = "none")
  se  <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coefn]

  # residual slide structure after this model
  R <- t(resid(lm.fit(dm, t(Mcb[idx, ])))); tot <- rowSums(sweep(R, 1, rowMeans(R))^2)
  X <- model.matrix(~ slide_f)
  r2slide <- mean(1 - rowSums(t(resid(lm.fit(X, t(R))))^2) / tot)

  # strata rank check (drop sex inside a stratum)
  ok <- TRUE
  for (lv in levels(sex)) {
    i <- which(sex == lv)
    ds <- cbind(model.matrix(~ ptsd + mdk$age + smoke + pos + propk)[i, , drop = FALSE],
                SV[i, , drop = FALSE])
    ds <- ds[, apply(ds, 2, function(z) length(unique(z)) > 1 | all(z == 1)), drop = FALSE]
    if (qr(ds)$rank < ncol(ds)) ok <- FALSE
  }

  res[[as.character(K)]] <- data.frame(
    k = K, npar = ncol(dm), resid_df = nrow(mdk) - ncol(dm),
    lambda = lam(tt$P.Value), med_se = median(se),
    se_ratio = median(se / base$se), frac_gain = mean(se < base$se),
    n_p1e5 = sum(tt$P.Value < 1e-5), n_fdr = sum(tt$adj.P.Val < 0.05),
    n_bonf = sum(tt$P.Value < 0.05 / nrow(tt)),
    r2_slide = r2slide, r2_smoke_sv = summary(lm(smoke ~ SV))$r.squared,
    max_r2_sv_ptsd = max(sapply(seq_len(K), function(i) summary(lm(SV[, i] ~ ptsd))$r.squared)),
    strata_full_rank = ok,
    df_F = 45 - (ncol(dm) - 1), df_M = 42 - (ncol(dm) - 1))
  svlist[[as.character(K)]] <- SV
  cat("done k =", K, "\n"); flush.console()
}

tabl <- do.call(rbind, c(list(data.frame(
  k = 0, npar = base$npar, resid_df = nrow(mdk) - base$npar,
  lambda = lam(base$tt$P.Value), med_se = median(base$se), se_ratio = 1, frac_gain = NA,
  n_p1e5 = sum(base$tt$P.Value < 1e-5), n_fdr = sum(base$tt$adj.P.Val < 0.05),
  n_bonf = sum(base$tt$P.Value < 0.05 / nrow(base$tt)),
  r2_slide = NA, r2_smoke_sv = NA, max_r2_sv_ptsd = NA, strata_full_rank = TRUE,
  df_F = 45 - (base$npar - 1), df_M = 42 - (base$npar - 1))), res))
print(tabl, row.names = FALSE, digits = 4)
write.csv(tabl, "data/05_k_lambda_sweep.csv", row.names = FALSE)
saveRDS(list(tabl = tabl, svlist = svlist, base_se = base$se,
             base_tt = base$tt), "data/05_k_lambda_sweep.rds")
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")
