## Build the Snakemake pipeline inputs from the v8.1-refiltered matrix.
## Full probe set now -- no 20k subsample (user ruling).
suppressPackageStartupMessages({ library(minfi); library(data.table) })
setwd("repo"); source("_setup.R")
t0 <- Sys.time()

grs <- readRDS(data_path("03_grs_filtered.rds"))
sva <- readRDS(data_path("05_sva.rds"))
keep <- sva$keep; mdk <- sva$mdk; propk <- sva$propk

betas <- getBeta(grs)[, keep]
M <- log2((betas + 1e-3) / (1 - betas + 1e-3))
M <- M[complete.cases(M), ]
cat("mvals:", nrow(M), "CpGs x", ncol(M), "samples\n")

out <- "../ewas_pipeline/data"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
## ewas.R expects SAMPLES IN ROWS, CpGs in columns, first column = sample IDs.
Mt <- t(M)
fwrite(data.table(SampleID = rownames(Mt), as.data.table(Mt)),
       file.path(out, "mvals.csv.gz"))
cat("mvals written as", nrow(Mt), "samples x", ncol(Mt), "CpGs\n")
cat("wrote mvals.csv.gz:",
    round(file.size(file.path(out, "mvals.csv.gz"))/1e6, 1), "MB\n")

## pheno: every column after the first becomes a model covariate in ewas.R,
## so this file IS the model. It must match chapter 06's limma model exactly,
## otherwise the two chapters fit different models to the same data:
##   methylation ~ PTSD + sex + age + 6 cell proportions + all 15 SVs
## The SVs are NOT regressed out of mvals.csv.gz -- they are carried as
## covariates here, so the pipeline adjusts for technical batch and for
## unmeasured biology (smoking above all; GSE132203 ships no smoking variable).
## All 15 are kept: num.sv()'s permutation test chose 15, and the technical and
## smoking-linked SVs are nearly disjoint, so truncation drops one or the other.
old <- fread("../ewas_pipeline/data/pheno.csv.orig", header = TRUE)
cat("original pheno cols:", paste(names(old), collapse = ", "), "\n")

ph <- data.table(
  SampleID = mdk$sample_id,
  PTSD = ifelse(mdk$ptsd == "Case", 1L, 0L),
  sex  = mdk$sex,
  age  = mdk$age
)
ph <- cbind(ph, as.data.table(propk)[, .(Neu, CD8T, CD4T, NK, Bcell, Mono)])

SV <- sva$SV
colnames(SV) <- paste0("SV", seq_len(ncol(SV)))
stopifnot(nrow(SV) == nrow(ph))
ph <- cbind(ph, as.data.table(SV))
cat("SV columns added:", ncol(SV), "\n")
stopifnot(identical(ph$SampleID, colnames(M)))
fwrite(ph, file.path(out, "pheno.csv"))
cat("wrote pheno.csv:", nrow(ph), "samples x", ncol(ph), "cols\n")
cat("cols:", paste(names(ph), collapse = ", "), "\n")
cat("PTSD table:", paste(table(ph$PTSD), collapse = "/"), "\n")
cat("elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
