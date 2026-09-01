## Regenerate ONLY repo/data/06_ewas.rds under the published v8.1 design.
##
## This is the EWAS + BACON half of repo/scripts/recompute_06_final.R, with the
## side effects that would disturb already-published outputs removed:
##   - 06_bacon_summary.rds is NOT rewritten (the on-disk copy carries extra
##     fields a later pass added, and chapter 06 reads four of them)
##   - the Manhattan / volcano / BACON diagnostic images are NOT redrawn
##     (the on-disk copies are the restyled, published ones)
## Every published anchor is asserted against 06_bacon_summary.rds at the end.
suppressPackageStartupMessages({
  library(limma); library(bacon)
})
setwd("repo"); source("_setup.R")
set.seed(42)
t0 <- Sys.time()

Mcb <- readRDS(data_path("05_mvals_combat.rds"))
sva <- readRDS(data_path("05_sva.rds"))
mdk <- sva$mdk; propk <- sva$propk; smoke <- sva$smoke
K <- sva$k_selected
SV <- sva$SV[, seq_len(K), drop = FALSE]
stopifnot(ncol(Mcb) == nrow(mdk), K == 6)
cat("EWAS matrix:", nrow(Mcb), "probes x", ncol(Mcb), "samples | k =", K, "\n")
cat("panel probes excluded:", sva$excluded_n, "\n")

ptsd <- relevel(factor(mdk$ptsd), ref = "Control")
pos <- factor(as.character(mdk$array_pos))
cat("ptsd:", paste(levels(ptsd), collapse = "/"), "=",
    paste(table(ptsd), collapse = "/"), "| positions:", nlevels(pos), "\n")

dd <- data.frame(ptsd = ptsd, sex = factor(mdk$sex), age = mdk$age,
                 smoke = smoke, pos = pos, propk, SV = SV, check.names = FALSE)
design <- model.matrix(~ ptsd + sex + age + smoke + pos + CD8T + CD4T + NK +
                         Bcell + Mono + Neu + SV, data = dd)
cat("design:", ncol(design), "columns | rank:", qr(design)$rank,
    "| residual df:", ncol(Mcb) - ncol(design), "\n")
coef_name <- "ptsdCase"
stopifnot(coef_name %in% colnames(design))

fit <- eBayes(lmFit(Mcb, design))
tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
tt$probe <- rownames(tt)

## delta-beta on the ComBat-adjusted scale actually modeled
bcb <- 2^Mcb / (2^Mcb + 1)
is_case <- ptsd == "Case"
db <- rowMeans(bcb[, is_case, drop = FALSE], na.rm = TRUE) -
      rowMeans(bcb[, !is_case, drop = FALSE], na.rm = TRUE)
tt$delta_beta <- db[tt$probe]
rm(bcb); gc(verbose = FALSE)

lam <- function(p) median(qchisq(1 - p, 1), na.rm = TRUE) / qchisq(0.5, 1)
lambda <- lam(tt$P.Value)
cat("lambda:", round(lambda, 7), "\n")

tt <- tt[, c("probe", "logFC", "delta_beta", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

## ---- BACON ----
es_v <- tt$logFC
se_v <- tt$logFC / tt$t
bc <- bacon(NULL, effectsizes = es_v, standarderrors = se_v)
infl <- inflation(bc); bi <- bias(bc)
tt$bacon.p <- pval(bc)
tt$bacon.es <- es(bc)
tt$bacon.se <- se(bc)
tt$bacon.adj.P <- p.adjust(tt$bacon.p, method = "BH")
lambda_bacon <- lam(tt$bacon.p)
cat("inflation:", round(infl, 4), "| bias:", round(bi, 4),
    "| lambda_raw:", round(lambda, 4), "| lambda_bacon:", round(lambda_bacon, 4), "\n")

saveRDS(list(tt = tt, lambda = lambda, n = ncol(Mcb), design_ncol = ncol(design),
             resid_df = ncol(Mcb) - ncol(design), k_sv = K,
             ncase = sum(is_case), nctrl = sum(!is_case),
             n_tested = nrow(tt), panel_excluded = sva$excluded_n),
        data_path("06_ewas.rds"))

## ---- verify against the published BACON summary ----
bs <- readRDS(data_path("06_bacon_summary.rds"))
bonf <- 0.05 / nrow(tt)
got <- list(n = ncol(Mcb), resid_df = ncol(Mcb) - ncol(design),
            n_tested = nrow(tt), inflation = infl, bias = bi,
            lambda_raw = lambda, lambda_bacon = lambda_bacon,
            bonf_thresh = bonf,
            n_bonf_bacon = sum(tt$bacon.p < bonf),
            n_fdr_bacon = sum(tt$bacon.adj.P < 0.05),
            n_p1e5_bacon = sum(tt$bacon.p < 1e-5),
            n_bonf_raw = sum(tt$P.Value < bonf),
            n_fdr_raw = sum(tt$adj.P.Val < 0.05),
            n_p1e5_raw = sum(tt$P.Value < 1e-5),
            top_probe = tt$probe[which.min(tt$bacon.p)],
            top_p_limma = min(tt$P.Value),
            top_p_bacon = min(tt$bacon.p))
cat("\n=== published vs regenerated ===\n")
ok <- TRUE
for (nm in names(got)) {
  if (is.null(bs[[nm]])) next
  a <- bs[[nm]]; b <- got[[nm]]
  same <- if (is.character(a) || is.character(b)) identical(as.character(a), as.character(b))
          else isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-6))
  if (!same) ok <- FALSE
  cat(sprintf("%-14s published=%-24s regenerated=%-24s %s\n", nm,
              format(a), format(b), if (same) "OK" else "*** DIFFERS ***"))
}
cat("\nALL ANCHORS MATCH:", ok, "\n")
cat("top 6 by BACON p:\n")
print(head(tt[order(tt$bacon.p), c("probe","logFC","delta_beta","P.Value","bacon.p","bacon.adj.P")], 6))
cat("elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
