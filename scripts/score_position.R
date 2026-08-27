suppressPackageStartupMessages({library(sva); library(limma)})
setwd("repo"); source("_setup.R"); setwd("..")
set.seed(42)
sv  <- readRDS("repo/data/05_sva.rds")
Mcb <- readRDS("repo/data/05_mvals_combat.rds"); if (is.list(Mcb)) Mcb <- Mcb[[1]]
Mb  <- readRDS("repo/data/05_mvals_combat_pos.rds")
md <- sv$mdk; smoke <- sv$smoke; P <- as.data.frame(sv$propk)
ptsd <- relevel(factor(md$ptsd), ref="Control")   # match ch06: Control is reference
sex <- factor(md$sex); age <- as.numeric(md$age)
slide <- factor(as.character(md$slide)); levels(slide) <- paste0("s", seq_len(nlevels(slide)))
pos <- factor(substr(as.character(md$array_pos), 1, 3))
cells <- c("Neu","NK","CD4T","CD8T","Bcell","Mono")
base <- data.frame(ptsd, sex, age, smoke, P[, cells])
modA <- model.matrix(~ ., base);          mod0A <- model.matrix(~ . - ptsd, base)
modC <- model.matrix(~ . + pos, base);    mod0C <- model.matrix(~ . - ptsd + pos, base)
stopifnot("ptsdCase" %in% colnames(modA), "ptsdCase" %in% colnames(modC))

sub <- order(apply(Mcb, 1, var), decreasing = TRUE)[1:20000]
r2 <- function(y, X) summary(lm(y ~ X))$r.squared
score <- function(Mat, mod, mod0, k, label) {
  s <- sva(Mat, mod, mod0, n.sv = k)
  X <- cbind(mod, SV = s$sv)
  f <- lmFit(Mat[sub, ], X)
  tt <- limma::topTable(eBayes(f), coef = "ptsdCase", number = Inf, sort.by = "none")
  lam <- median(qchisq(1 - tt$P.Value, 1)) / qchisq(0.5, 1)
  res <- t(residuals(f, Mat[sub, ]))
  lo <- mean(apply(res, 2, function(y) r2(y, slide)))
  po <- mean(apply(res, 2, function(y) r2(y, pos)))
  se <- mean(apply(res, 2, function(y) r2(y, sex)))
  np <- ncol(X)
  cat(sprintf("%s | k=%2d par=%2d | lambda=%.3f | resid R2: slide %.4f pos %.4f sex %.4f | df comb %d F %d M %d | n_p<1e-5 %d\n",
      label, k, np, lam, lo, po, se, 87-np, 45-(np-1), 42-(np-1), sum(tt$P.Value < 1e-5)))
  list(sv = s$sv, lambda = lam, npar = np, slide_r2 = lo, pos_r2 = po, sex_r2 = se)
}
cat("=== residual R2 of the INPUT matrices (no model), for reference:\n")
for (nm in c("Mcb","Mb")) {
  Mx <- get(nm); rr <- t(Mx[sub, ])
  cat(sprintf("  %s: slide %.4f pos %.4f sex %.4f\n", nm,
    mean(apply(rr,2,function(y) r2(y,slide))), mean(apply(rr,2,function(y) r2(y,pos))),
    mean(apply(rr,2,function(y) r2(y,sex)))))
}
cat("\n=== the three designs:\n")
A <- score(Mcb, modA, mod0A, 15, "A ComBat(slide), pos->SVA   ")
B <- score(Mb,  modA, mod0A, 25, "B ComBat(slide)+ComBat(pos) ")
C <- score(Mcb, modC, mod0C, 15, "C ComBat(slide), pos in mod ")
saveRDS(list(num_sv=c(A=15,B=25,C=15), A=A[-1], B=B[-1], C=C[-1], svC=C$sv, svA=A$sv),
        "repo/data/05_position_comparison.rds")
cat("\nsaved 05_position_comparison.rds\n")
