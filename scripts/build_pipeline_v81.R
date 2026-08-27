## Rebuild the Snakemake pipeline inputs so chapter 07 fits EXACTLY the model
## chapter 06 fits:  ~ PTSD + sex + age + smoke + pos + 6 cell props + 6 SVs
## Inputs come from the stratified-ComBat matrix (05_mvals_combat.rds), the same
## matrix limma sees in chapter 06 -- not the raw filtered betas.
suppressPackageStartupMessages({ library(data.table) })
setwd("repo"); source("_setup.R")
t0 <- Sys.time()

Mcb <- readRDS(data_path("05_mvals_combat.rds"))
sva <- readRDS(data_path("05_sva.rds"))
mdk <- sva$mdk; propk <- sva$propk; smoke <- sva$smoke
K   <- sva$k_selected
stopifnot(K == 6, ncol(Mcb) == nrow(mdk))
SV  <- sva$SV[, seq_len(K), drop = FALSE]
colnames(SV) <- paste0("SV", seq_len(K))

out <- "../ewas_pipeline/data"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cat("ComBat matrix:", nrow(Mcb), "CpGs x", ncol(Mcb), "samples | k =", K, "\n")

## ---- pheno: every column after the first IS a model term in ewas.R ----------
## Column order here defines the glm formula, so it must mirror chapter 06.
ph <- data.table(SampleID = mdk$sample_id,
                 PTSD  = ifelse(mdk$ptsd == "Case", 1L, 0L),
                 sex   = mdk$sex,
                 age   = mdk$age,
                 smoke = as.numeric(smoke),
                 pos   = as.character(mdk$array_pos))
ph <- cbind(ph, as.data.table(propk)[, .(CD8T, CD4T, NK, Bcell, Mono, Neu)])
ph <- cbind(ph, as.data.table(SV))
stopifnot(identical(ph$SampleID, colnames(Mcb)))

fwrite(ph, file.path(out, "pheno.csv"))
cat("pheno.csv:", nrow(ph), "samples x", ncol(ph), "cols\n")
cat("cols:", paste(names(ph), collapse = ", "), "\n")
cat("PTSD table:", paste(table(ph$PTSD), collapse = "/"), "\n")
cat("sex table:", paste(paste(names(table(ph$sex)), table(ph$sex), sep = "="),
                       collapse = " "), "\n")
cat("positions:", length(unique(ph$pos)), "->",
    paste(sort(unique(ph$pos)), collapse = ","), "\n")

## combined-arm parameter count: intercept + PTSD + sex + age + smoke
## + (npos-1) dummies + 6 props + K SVs
npar <- 1 + 1 + 1 + 1 + 1 + (length(unique(ph$pos)) - 1) + 6 + K
cat("combined npar =", npar, "-> resid df =", nrow(ph) - npar,
    "  (chapter 06: 24 ->", nrow(ph) - 24, ")\n")

## per-stratum: sex is dropped by stratify.R; positions may be fewer
for (s in unique(ph$sex)) {
  sub <- ph[sex == s]
  np <- 1 + 1 + 1 + 1 + (length(unique(sub$pos)) - 1) + 6 + K
  cat(sprintf("  stratum %s: n=%d positions=%d npar=%d resid_df=%d\n",
              s, nrow(sub), length(unique(sub$pos)), np, nrow(sub) - np))
}

## ---- mvals: ewas.R expects SAMPLES IN ROWS, CpGs in columns ----------------
Mt <- t(Mcb)
fwrite(data.table(SampleID = rownames(Mt), as.data.table(Mt)),
       file.path(out, "mvals.csv.gz"))
cat("mvals.csv.gz:", nrow(Mt), "samples x", ncol(Mt), "CpGs |",
    round(file.size(file.path(out, "mvals.csv.gz")) / 1e6, 1), "MB\n")

cat("elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
