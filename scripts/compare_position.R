## Three ways to handle measured array position, compared on the same data.
##  A: ComBat(slide, within sex) only          -> position left to SVA
##  B: ComBat(slide) then ComBat(position)     -> sequential correction
##  C: ComBat(slide) + position as fixed covariate in the model
suppressPackageStartupMessages({library(sva); library(limma)})
setwd("repo"); source("_setup.R"); setwd("..")
set.seed(42)

sv  <- readRDS("repo/data/05_sva.rds")
Mcb <- readRDS("repo/data/05_mvals_combat.rds")
if (is.list(Mcb)) Mcb <- Mcb[[1]]
md  <- sv$mdk; smoke <- sv$smoke; P <- as.data.frame(sv$propk)
stopifnot(ncol(Mcb) == nrow(md))

ptsd <- factor(md$ptsd); sex <- factor(md$sex); age <- as.numeric(md$age)
slide <- factor(as.character(md$slide)); levels(slide) <- paste0("s", seq_len(nlevels(slide)))
pos <- factor(substr(as.character(md$array_pos), 1, 3))
cells <- c("Neu","NK","CD4T","CD8T","Bcell","Mono")

## ---- arm B: second ComBat pass on position, within sex stratum -------------
Mb <- Mcb
for (s in levels(sex)) {
  i <- which(sex == s)
  cbm <- model.matrix(~ ptsd[i] + age[i] + smoke[i] + as.matrix(P[i, cells]))
  Mb[, i] <- ComBat(dat = Mcb[, i], batch = droplevels(pos[i]), mod = cbm,
                    par.prior = TRUE, mean.only = FALSE)
}
saveRDS(Mb, "repo/data/05_mvals_combat_pos.rds")

## ---- num.sv for each arm ---------------------------------------------------
base  <- data.frame(ptsd, sex, age, smoke, P[, cells])
modA  <- model.matrix(~ ., base);                mod0A <- model.matrix(~ . - ptsd, base)
modC  <- model.matrix(~ . + pos, base);          mod0C <- model.matrix(~ . - ptsd + pos, base)
nsv <- list()
nsv$A <- sv$n.sv                                   # already computed = 15
nsv$B <- num.sv(Mb,  modA, method = "be")
nsv$C <- num.sv(Mcb, modC, method = "be")
cat(sprintf("num.sv  A=%d  B=%d  C=%d\n", nsv$A, nsv$B, nsv$C))

## ---- run SVA and score each arm on a common 20k subset ---------------------
sub <- order(apply(Mcb, 1, var), decreasing = TRUE)[1:20000]
r2 <- function(y, X) summary(lm(y ~ X))$r.squared
score <- function(Mat, mod, mod0, k, label, extra_df) {
  s  <- sva(Mat, mod, mod0, n.sv = k)
  X  <- cbind(mod, s$sv)
  fit <- eBayes(lmFit(Mat[sub, ], X))
  tt <- limma::topTable(fit, coef = "ptsdCase", number = Inf, sort.by = "none")
  lam <- median(qchisq(1 - tt$P.Value, 1)) / qchisq(0.5, 1)
  res <- t(residuals(lmFit(Mat[sub, ], X), Mat[sub, ]))
  lo <- mean(apply(res, 2, function(y) r2(y, slide)))
  po <- mean(apply(res, 2, function(y) r2(y, pos)))
  se <- mean(apply(res, 2, function(y) r2(y, sex)))
  npar <- ncol(X)
  cat(sprintf("%s | k=%2d par=%2d | lambda=%.3f | resid R2: slide %.4f pos %.4f sex %.4f | df comb %d F %d M %d\n",
      label, k, npar, lam, lo, po, se, 87-npar, 45-(npar-1+extra_df), 42-(npar-1+extra_df)))
  invisible(list(sv=s, lambda=lam, npar=npar))
}
## in-stratum parameter count: drop sex (constant), keep the rest
outA <- score(Mcb, modA, mod0A, nsv$A, "A ComBat(slide), pos->SVA  ", 0)
outB <- score(Mb,  modA, mod0A, nsv$B, "B ComBat(slide)+ComBat(pos)", 0)
outC <- score(Mcb, modC, mod0C, nsv$C, "C ComBat(slide), pos in mod", 0)

saveRDS(list(nsv=nsv, A=outA[c("lambda","npar")], B=outB[c("lambda","npar")],
             C=outC[c("lambda","npar")], svC=outC$sv$sv, svB=outB$sv$sv),
        "repo/data/05_position_comparison.rds")
