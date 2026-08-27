suppressPackageStartupMessages({library(limma)})
setwd("repo"); source("_setup.R")
w   <- readRDS("data/05_sv_worth_it.rds")
Mcb <- readRDS("data/05_mvals_combat.rds")
s   <- readRDS("data/05_sva.rds")
mdk <- s$mdk; propk <- s$propk; smoke <- s$smoke; SV <- s$SV

ptsd <- relevel(factor(mdk$ptsd), ref="Control"); sex <- factor(mdk$sex)
pos  <- factor(as.character(mdk$array_pos)); slide_f <- factor(as.character(mdk$slide))
d0 <- model.matrix(~ ptsd + sex + mdk$age + smoke + pos + propk); d6 <- cbind(d0, SV)

v <- apply(Mcb, 1, var); idx <- order(v, decreasing=TRUE)[1:20000]
r2resid <- function(dm) {
  R  <- t(resid(lm.fit(dm, t(Mcb[idx, ]))))          # probes x samples
  tot <- rowSums(sweep(R, 1, rowMeans(R))^2)
  ss <- function(f) { X <- model.matrix(~ f)
    1 - rowSums(t(resid(lm.fit(X, t(R))))^2) / tot }
  c(slide=mean(ss(slide_f)), pos=mean(ss(pos)), sex=mean(ss(sex)))
}
a <- r2resid(d0); b <- r2resid(d6)
cat(sprintf("residual R2  k=0 : slide %.4f  pos %.4f  sex %.4f\n", a[1], a[2], a[3]))
cat(sprintf("residual R2  k=6 : slide %.4f  pos %.4f  sex %.4f\n", b[1], b[2], b[3]))

lam <- function(p) median(qchisq(1-p,1), na.rm=TRUE)/qchisq(0.5,1)
t0 <- w$tt0; t6 <- w$tt6
cat(sprintf("\nk=0 : par %d  resid df %d  lambda %.4f  p<1e-5 %d  FDR<.05 %d  Bonf %d  medSE %.5f\n",
  w$npar0, 87-w$npar0, lam(t0$P.Value), sum(t0$P.Value<1e-5), sum(t0$adj.P.Val<0.05),
  sum(t0$P.Value < 0.05/nrow(t0)), median(w$se0)))
cat(sprintf("k=6 : par %d  resid df %d  lambda %.4f  p<1e-5 %d  FDR<.05 %d  Bonf %d  medSE %.5f\n",
  w$npar6, 87-w$npar6, lam(t6$P.Value), sum(t6$P.Value<1e-5), sum(t6$adj.P.Val<0.05),
  sum(t6$P.Value < 0.05/nrow(t6)), median(w$se6)))
cat("\ncor(t)              :", round(cor(t0$t, t6$t),4), "\n")
cat("cor(-log10 p)       :", round(cor(-log10(t0$P.Value), -log10(t6$P.Value)),4), "\n")
cat("median |dlogFC|     :", signif(median(abs(t0$logFC-t6$logFC)),3),
    " vs median |logFC| ", signif(median(abs(t0$logFC)),3), "\n")
cat("median SE ratio 6/0 :", round(median(w$se6/w$se0),4), "\n")
o0 <- rownames(t0)[order(t0$P.Value)][1:20]; o6 <- rownames(t6)[order(t6$P.Value)][1:20]
cat("top-20 overlap      :", length(intersect(o0,o6)), "/20\n")
q <- c(0.5,0.9,0.99,0.999)
cat("quantiles |dt|      :", paste(round(quantile(abs(t0$t-t6$t), q),3), collapse="  "), "(50/90/99/99.9%)\n")
w$r2_0 <- a; w$r2_6 <- b; saveRDS(w, "data/05_sv_worth_it.rds")
