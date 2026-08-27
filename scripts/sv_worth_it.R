# Is the SV block earning its degrees of freedom?
# Compare the adopted design with and without the 6 SVs on the FULL probe set.
suppressPackageStartupMessages({library(limma)})
set.seed(42)
t0 <- Sys.time()
setwd("repo"); source("_setup.R")

Mcb <- readRDS("data/05_mvals_combat.rds")
s   <- readRDS("data/05_sva.rds")
mdk <- s$mdk; propk <- s$propk; smoke <- s$smoke; SV <- s$SV
stopifnot(ncol(SV) == 6, identical(colnames(Mcb), as.character(mdk$sample_id)))

ptsd <- relevel(factor(mdk$ptsd), ref = "Control")
sex  <- factor(mdk$sex); pos <- factor(as.character(mdk$array_pos))
slide_f <- factor(as.character(mdk$slide))

d0 <- model.matrix(~ ptsd + sex + mdk$age + smoke + pos + propk)
d6 <- cbind(d0, SV)
colnames(d0) <- make.names(colnames(d0)); colnames(d6) <- make.names(colnames(d6))
coefn <- grep("^ptsdCase$", colnames(d0), value = TRUE); stopifnot(length(coefn) == 1)

lam <- function(p) median(qchisq(1 - p, 1), na.rm = TRUE) / qchisq(0.5, 1)

run <- function(dm, tag) {
  fit <- eBayes(lmFit(Mcb, dm))
  tt  <- limma::topTable(fit, coef = coefn, number = Inf, sort.by = "none")
  se  <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coefn]
  cat(sprintf("%s: par %d  resid df %d  lambda %.4f  p<1e-5 %d  FDR<.05 %d  Bonf %d  median SE %.5f\n",
      tag, ncol(dm), nrow(mdk) - ncol(dm), lam(tt$P.Value),
      sum(tt$P.Value < 1e-5), sum(tt$adj.P.Val < 0.05),
      sum(tt$P.Value < 0.05 / nrow(tt)), median(se)))
  list(tt = tt, se = se, fit = fit, dm = dm)
}

r0 <- run(d0, "k=0 ")
r6 <- run(d6, "k=6 ")

# --- how much does the SV block move the inference?
cat("\ncor(t-stat)        :", round(cor(r0$tt$t, r6$tt$t), 4), "\n")
cat("cor(-log10 p)      :", round(cor(-log10(r0$tt$P.Value), -log10(r6$tt$P.Value)), 4), "\n")
cat("median |dlogFC|    :", signif(median(abs(r0$tt$logFC - r6$tt$logFC)), 3),
    " (median |logFC| =", signif(median(abs(r0$tt$logFC)), 3), ")\n")
cat("median SE ratio 6/0:", round(median(r6$se / r0$se), 4), "\n")
top0 <- rownames(r0$tt)[order(r0$tt$P.Value)][1:20]
top6 <- rownames(r6$tt)[order(r6$tt$P.Value)][1:20]
cat("top-20 overlap     :", length(intersect(top0, top6)), "/ 20\n")

# --- residual batch structure left behind by each model (20k most variable probes)
v  <- apply(Mcb, 1, var); idx <- order(v, decreasing = TRUE)[1:20000]
r2resid <- function(dm) {
  R <- t(resid(lm.fit(dm, t(Mcb[idx, ]))))
  ss <- function(f) { X <- model.matrix(~ f)
    1 - colSums(t(resid(lm.fit(X, t(R))))^2) / rowSums(sweep(R, 1, rowMeans(R))^2) }
  c(slide = mean(ss(slide_f)), pos = mean(ss(pos)), sex = mean(ss(sex)))
}
cat("\nresidual R2 k=0:", paste(names(r2resid(d0)), round(r2resid(d0), 4), collapse = "  "), "\n")
cat("residual R2 k=6:", paste(names(r2resid(d6)), round(r2resid(d6), 4), collapse = "  "), "\n")

saveRDS(list(tt0 = r0$tt, tt6 = r6$tt, se0 = r0$se, se6 = r6$se,
             lam0 = lam(r0$tt$P.Value), lam6 = lam(r6$tt$P.Value),
             r2_0 = r2resid(d0), r2_6 = r2resid(d6),
             npar0 = ncol(d0), npar6 = ncol(d6)),
        "data/05_sv_worth_it.rds")
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")
