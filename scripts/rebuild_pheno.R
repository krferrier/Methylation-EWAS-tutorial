## Rebuild ONLY pheno.csv with the SV covariates (mvals.csv.gz is unchanged).
suppressPackageStartupMessages({ library(data.table) })
setwd("repo"); source("_setup.R")

sva <- readRDS(data_path("05_sva.rds"))
mdk <- sva$mdk; propk <- sva$propk; SV <- sva$SV
colnames(SV) <- paste0("SV", seq_len(ncol(SV)))
out <- "../ewas_pipeline/data"

ph <- data.table(SampleID = mdk$sample_id,
                 PTSD = ifelse(mdk$ptsd == "Case", 1L, 0L),
                 sex = mdk$sex, age = mdk$age)
ph <- cbind(ph, as.data.table(propk)[, .(Neu, CD8T, CD4T, NK, Bcell, Mono)])
ph <- cbind(ph, as.data.table(SV))

## sample order must still match the mvals header
hdr <- fread(file.path(out, "mvals.csv.gz"), select = 1L)
stopifnot(identical(as.character(hdr[[1]]), ph$SampleID))
cat("sample order matches mvals:", nrow(hdr), "samples\n")

fwrite(ph, file.path(out, "pheno.csv"))
cat("pheno.csv:", nrow(ph), "samples x", ncol(ph), "cols\n")
cat("cols:", paste(names(ph), collapse = ", "), "\n")
cat("PTSD table:", paste(table(ph$PTSD), collapse = "/"), "\n")
cat("sex table:", paste(names(table(ph$sex)), table(ph$sex), collapse = " "), "\n")
cat("model params (combined) =", 1 + 1 + 1 + 1 + 6 + ncol(SV),
    "-> resid df =", nrow(ph) - (1 + 1 + 1 + 1 + 6 + ncol(SV)), "\n")
