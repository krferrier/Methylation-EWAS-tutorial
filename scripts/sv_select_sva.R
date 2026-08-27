## Choosing how many SVs to carry, using criteria appropriate to SVA.
##
## A scree/percent-variance rule does not transfer: SVA returns an unordered
## basis, so SV index carries no variance ranking and there is no elbow.
## The SVA-native dimensionality criterion is num.sv()'s permutation test,
## which already ran (method "be") and returned n.sv. Everything below is about
## whether a TRUNCATED subset is defensible for the low-df stratified arms,
## and if so which SVs to keep.
##
## SVs are meant to absorb technical batch AND unmeasured biology - here
## smoking above all, since GSE132203 ships no smoking variable. So rank SVs
## against both: measured technical factors, and a methylation-derived smoking
## proxy built from established smoking CpGs.
suppressPackageStartupMessages({
  library(minfi); library(data.table); library(ggplot2); library(sva)
})
setwd("repo"); source("_setup.R")
set.seed(42)

g   <- readRDS(data_path("03_grs_filtered.rds"))
sva_o <- readRDS(data_path("05_sva.rds"))
SV  <- sva_o$SV; mdk <- sva_o$mdk; propk <- sva_o$propk
beta <- minfi::getBeta(g)[, sva_o$keep, drop = FALSE]
cat("beta:", paste(dim(beta), collapse = " x "), "| n.sv from num.sv(be):", ncol(SV), "\n")

## ---- 1. smoking proxy from established smoking-associated CpGs -------------
## Well-replicated loci (Joehanes 2016 and predecessors). Sign of effect is NOT
## taken from memory: the panel's shared axis is extracted as PC1 and oriented
## using the single most robust direction in the literature - AHRR cg05575921
## is hypomethylated in smokers.
panel <- list(
  AHRR      = c("cg05575921","cg21161138","cg23576855","cg25648203","cg26703534","cg11902777"),
  F2RL3     = c("cg03636183","cg21911711"),
  GPR15     = c("cg19859270"),
  `2q37.1`  = c("cg05951221","cg21566642","cg01940273","cg03329539","cg06126421"),
  GFI1      = c("cg09935388","cg12876356","cg18146737","cg06338710"),
  MYO1G     = c("cg12803068","cg22132788","cg04180046","cg19089201"),
  PRSS23    = c("cg14391737"),
  RARA      = c("cg17739917")
)
pv <- unlist(panel, use.names = FALSE)
gene_of <- rep(names(panel), lengths(panel))
present <- pv %in% rownames(beta)
cat("smoking panel probes:", length(pv), "| present after v8.1 filtering:", sum(present), "\n")
if (any(!present)) cat("  absent:", paste(pv[!present], collapse = ", "), "\n")
pv2 <- pv[present]; gene2 <- gene_of[present]

Bp <- beta[pv2, , drop = FALSE]
Z  <- t(scale(t(Bp)))                       # z-score each CpG across samples
Z  <- Z[rowSums(is.na(Z)) == 0, , drop = FALSE]
pc <- prcomp(t(Z), center = TRUE, scale. = FALSE)
pve_panel <- pc$sdev^2 / sum(pc$sdev^2)
score <- pc$x[, 1]
## orient: smokers are hypomethylated at AHRR cg05575921, so the score must be
## NEGATIVELY correlated with that probe's beta
anchor <- "cg05575921"
if (anchor %in% rownames(Bp)) {
  if (cor(score, Bp[anchor, ], use = "complete.obs") > 0) score <- -score
  cat("anchor cg05575921 cor with oriented score:",
      round(cor(score, Bp[anchor, ], use = "complete.obs"), 3), "\n")
}
cat("panel PC1 explains", round(100 * pve_panel[1], 1), "% of panel variance (PC2:",
    round(100 * pve_panel[2], 1), "%)\n")

## does the panel cohere? loadings per gene on PC1
ld <- data.table(probe = rownames(Z), gene = gene2[match(rownames(Z), pv2)],
                 loading = pc$rotation[, 1])
ld[, oriented := loading * sign(cor(pc$x[,1], score))]
cat("\nPC1 loadings by locus (oriented so + = more smoking-like):\n")
print(ld[order(-oriented), .(probe, gene, loading = round(oriented, 3))])

## sanity: the proxy should track smoking-like biology, not the chip
r2 <- function(y, x) {
  d <- data.frame(y = y, x = x); d <- d[complete.cases(d), ]
  if (length(unique(d$x)) < 2) return(NA_real_)
  summary(lm(y ~ x, data = d))$r.squared
}
chip <- factor(as.character(mdk$slide))
apos <- factor(substr(mdk$array_pos, 1, 3))
cat("\nsmoking proxy vs chip R2:", round(r2(score, chip), 3),
    "| vs array position:", round(r2(score, apos), 3),
    "| vs PTSD:", round(r2(score, factor(mdk$ptsd)), 3),
    "| vs age:", round(r2(score, mdk$age), 3),
    "| vs sex:", round(r2(score, factor(mdk$sex)), 3), "\n")

## ---- 2. rank each SV against technical factors AND the smoking proxy ------
svv <- readRDS(data_path("05_sv_variance.rds"))
d <- data.table(sv = paste0("SV", seq_len(ncol(SV))), k = seq_len(ncol(SV)))
d[, pve      := 100 * svv$pve]
d[, chip     := vapply(seq_len(ncol(SV)), function(i) r2(SV[,i], chip), 0)]
d[, arraypos := vapply(seq_len(ncol(SV)), function(i) r2(SV[,i], apos), 0)]
d[, smoking  := vapply(seq_len(ncol(SV)), function(i) r2(SV[,i], score), 0)]
d[, PTSD     := vapply(seq_len(ncol(SV)), function(i) r2(SV[,i], factor(mdk$ptsd)), 0)]
d[, tech     := pmax(chip, arraypos)]
## a single "nuisance capture" score: the strongest thing an SV is absorbing
d[, capture  := pmax(tech, smoking)]
setorder(d, -capture)
d[, rank := .I]
cat("\nSVs ranked by strongest nuisance captured (technical or smoking):\n")
print(d[, .(rank, sv, pve = round(pve,2), chip = round(chip,3),
            arraypos = round(arraypos,3), smoking = round(smoking,3),
            capture = round(capture,3), PTSD = round(PTSD,3))])

cat("\nsmoking-leading SVs (smoking > tech):",
    paste(d[smoking > tech]$sv, collapse = ", "), "\n")
cat("total smoking R2 across all SVs:", round(sum(d$smoking), 3), "\n")
saveRDS(list(ranking = d, score = score, panel_loadings = ld,
             pve_panel = pve_panel),
        data_path("05_sv_selection_inputs.rds"))
fwrite(d, data_path("05_sv_ranking.csv"))
