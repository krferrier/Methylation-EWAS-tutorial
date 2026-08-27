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
sub <- order(apply(Mcb, 1, var), decreasing = TRUE)[1:20000]
r2 <- function(y, X) summary(lm(y ~ X))$r.squared
res <- list()
for (k in c(6, 8, 10, 15)) {
  s <- sva(Mcb, modC, mod0C, n.sv = k)
  X <- cbind(modC, SV = s$sv); f <- lmFit(Mcb[sub, ], X)
  tt <- limma::topTable(eBayes(f), coef = "ptsdCase", number = Inf, sort.by = "none")
  lam <- median(qchisq(1 - tt$P.Value, 1)) / qchisq(0.5, 1)
  rr <- t(residuals(f, Mcb[sub, ]))
  lo <- mean(apply(rr, 2, function(y) r2(y, slide)))
  np <- ncol(X); nps <- np - 1   # sex drops inside a stratum
  ok <- all(sapply(c("F","M"), function(ss) { i <- sex==ss
        Xs <- cbind(1, as.numeric(ptsd[i])-1, age[i], smoke[i], as.matrix(P[i,cells]),
                    model.matrix(~ droplevels(pos[i]))[,-1], s$sv[i,,drop=FALSE])
        qr(Xs)$rank == ncol(Xs) }))
  cat(sprintf("k=%2d | par=%2d | lambda=%.3f | resid slide R2 %.4f | df comb %d F %d M %d | strata full rank %s | R2(smoke~SV) %.3f\n",
      k, np, lam, lo, 87-np, 45-nps, 42-nps, ok, r2(smoke, s$sv)))
  res[[as.character(k)]] <- list(sv=s$sv, lambda=lam, npar=np, slide_r2=lo, full_rank=ok)
}
saveRDS(res, "repo/data/05_k_sweep_posmod.rds"); cat("saved\n")
