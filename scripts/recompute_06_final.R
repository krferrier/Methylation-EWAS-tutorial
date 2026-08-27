## Recompute chapter 06: limma EWAS on the stratified-ComBat matrix with the
## adopted v8.1 design:  ~ ptsd + sex + age + smoke + pos + 6 cell props + 6 SVs
suppressPackageStartupMessages({
  library(limma); library(bacon); library(ggplot2); library(minfi)
})
setwd("repo"); source("_setup.R")
set.seed(42)
t0 <- Sys.time()

Mcb <- readRDS(data_path("05_mvals_combat.rds"))
sva <- readRDS(data_path("05_sva.rds"))
mdk <- sva$mdk; propk <- sva$propk; smoke <- sva$smoke
K   <- sva$k_selected
SV  <- sva$SV[, seq_len(K), drop = FALSE]
stopifnot(ncol(Mcb) == nrow(mdk), K == 6)
cat("EWAS matrix:", nrow(Mcb), "probes x", ncol(Mcb), "samples | k =", K, "\n")
cat("panel probes excluded:", sva$excluded_n, "\n")

ptsd <- relevel(factor(mdk$ptsd), ref = "Control")
pos  <- factor(as.character(mdk$array_pos))
cat("ptsd:", paste(levels(ptsd), collapse="/"), "=", paste(table(ptsd), collapse="/"),
    "| positions:", nlevels(pos), "\n")

dd <- data.frame(ptsd = ptsd, sex = factor(mdk$sex), age = mdk$age,
                 smoke = smoke, pos = pos, propk, SV = SV, check.names = FALSE)
design <- model.matrix(~ ptsd + sex + age + smoke + pos + CD8T + CD4T + NK +
                         Bcell + Mono + Neu + SV, data = dd)
cat("design:", ncol(design), "columns | rank:", qr(design)$rank,
    "| residual df:", ncol(Mcb) - ncol(design), "\n")
coef_name <- "ptsdCase"
stopifnot(coef_name %in% colnames(design))

fit <- eBayes(lmFit(Mcb, design))
tt  <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
tt$probe <- rownames(tt)

## delta-beta on the ComBat-adjusted scale actually modeled
bcb <- 2^Mcb / (2^Mcb + 1)
is_case <- ptsd == "Case"
db <- rowMeans(bcb[, is_case, drop = FALSE], na.rm = TRUE) -
      rowMeans(bcb[, !is_case, drop = FALSE], na.rm = TRUE)
tt$delta_beta <- db[tt$probe]
rm(bcb); gc(verbose = FALSE)

## Genomic inflation factor, computed from p-values (the conventional definition,
## and the same one used for the k-selection sweep in chapter 05).
lam <- function(p) median(qchisq(1 - p, 1), na.rm = TRUE) / qchisq(0.5, 1)
lambda <- lam(tt$P.Value)
cat("lambda:", round(lambda, 7), "\n")

tt <- tt[, c("probe", "logFC", "delta_beta", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

## ---- BACON ----
es_v <- tt$logFC
se_v <- tt$logFC / tt$t
bc <- bacon(NULL, effectsizes = es_v, standarderrors = se_v)
infl <- inflation(bc); bi <- bias(bc)
tt$bacon.p  <- pval(bc)
tt$bacon.es <- es(bc)
tt$bacon.se <- se(bc)
tt$bacon.adj.P <- p.adjust(tt$bacon.p, method = "BH")
lambda_bacon <- lam(tt$bacon.p)
cat("inflation:", round(infl,4), "| bias:", round(bi,4),
    "| lambda_raw:", round(lambda,4), "| lambda_bacon:", round(lambda_bacon,4), "\n")

saveRDS(list(inflation = infl, bias = bi, lambda_raw = lambda,
             lambda_bacon = lambda_bacon, n = ncol(Mcb)),
        data_path("06_bacon_summary.rds"))
write.csv(tt, gzfile(data_path("06_ewas_bacon_toptable.csv.gz")), row.names = FALSE)

## Save 06_ewas.rds AFTER BACON so tt carries bacon.p / bacon.es / bacon.se /
## bacon.adj.P. Chapter 06 shows the BACON call with `eval: false`, so the
## rendered hits chunk must find those columns already present in ew$tt.
saveRDS(list(tt = tt, lambda = lambda, n = ncol(Mcb), design_ncol = ncol(design),
             resid_df = ncol(Mcb) - ncol(design), k_sv = K,
             ncase = sum(is_case), nctrl = sum(!is_case),
             n_tested = nrow(tt), panel_excluded = sva$excluded_n),
        data_path("06_ewas.rds"))

## ---- BACON diagnostics ----
dir.create(data_path("06_bacon"), showWarnings = FALSE)
jpeg(data_path("06_bacon/fit.jpg"), width = 1600, height = 1200, res = 200)
print(fit(bc, n = 100)); dev.off()
jpeg(data_path("06_bacon/qq.jpg"), width = 1600, height = 1200, res = 200)
print(plot(bc, type = "qq")); dev.off()
jpeg(data_path("06_bacon/traces.jpg"), width = 1600, height = 1400, res = 180)
traces(bc, burnin = FALSE); dev.off()
jpeg(data_path("06_bacon/posteriors.jpg"), width = 1600, height = 1200, res = 200)
posteriors(bc); dev.off()

## ---- Manhattan + volcano ----
grs <- readRDS(data_path("03_grs_filtered.rds"))
ann <- getAnnotation(grs)[tt$probe, c("chr", "pos")]
rm(grs); gc(verbose = FALSE)
pl <- data.frame(probe = tt$probe, chr = ann$chr, pos = ann$pos,
                 p = tt$P.Value, logFC = tt$logFC, db = tt$delta_beta)
pl$chrn <- as.integer(sub("chr", "", pl$chr))
pl <- pl[!is.na(pl$chrn), ]
pl <- pl[order(pl$chrn, pl$pos), ]
off <- cumsum(c(0, tapply(pl$pos, pl$chrn, max, na.rm = TRUE)))
names(off) <- c(names(tapply(pl$pos, pl$chrn, max)), "end")
pl$x <- pl$pos + off[as.character(pl$chrn)]
ctr <- tapply(pl$x, pl$chrn, function(z) (min(z)+max(z))/2)
bonf <- 0.05 / nrow(tt)

pm <- ggplot(pl, aes(x, -log10(p), color = factor(chrn %% 2))) +
  geom_point(size = 0.35, alpha = 0.75) +
  geom_hline(yintercept = -log10(bonf), color = "#C44E52", linewidth = 0.5) +
  geom_hline(yintercept = -log10(1e-5), color = "gray45", linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(values = c("0" = "#4C72B0", "1" = "#8FA8CC"), guide = "none") +
  scale_x_continuous(breaks = ctr, labels = names(ctr), expand = c(0.01, 0)) +
  labs(x = "Chromosome", y = expression(-log[10](p)),
       title = sprintf("PTSD EWAS, %s CpGs tested", format(nrow(tt), big.mark = ","))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
ggsave(data_path("06_manhattan.png"), plot = pm, width = 10, height = 4.5, dpi = 200)

pv <- ggplot(tt, aes(delta_beta * 100, -log10(P.Value))) +
  geom_point(size = 0.35, alpha = 0.55, color = "#4C72B0") +
  geom_hline(yintercept = -log10(1e-5), color = "gray45", linetype = "dashed", linewidth = 0.4) +
  labs(x = "Difference in methylation, case - control (percentage points)",
       y = expression(-log[10](p)), title = "Effect size versus significance") +
  theme_minimal(base_size = 11)
ggsave(data_path("06_volcano.png"), plot = pv, width = 7.5, height = 5, dpi = 200)

cat("\nprobes tested:", nrow(tt), "\n")
cat("Bonferroni:", sum(tt$P.Value < bonf), "| FDR<0.05:", sum(tt$adj.P.Val < 0.05),
    "| p<1e-5:", sum(tt$P.Value < 1e-5), "\n")
cat("BACON FDR<0.05:", sum(tt$bacon.adj.P < 0.05),
    "| BACON Bonferroni:", sum(tt$bacon.p < bonf), "\n")
cat("top 10:\n")
print(head(tt[, c("probe","logFC","delta_beta","P.Value","adj.P.Val","bacon.p")], 10))
cat("elapsed:", round(difftime(Sys.time(), t0, units="mins"),2), "min\n")
