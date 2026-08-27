## Recompute chapter 08 Part 3: gometh, methylGSA, KYCG on the regenerated toptable.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setwd("repo"); source("_setup.R")
t0 <- Sys.time()
outd <- data_path("08_annotation")

tt <- fread(data_path("06_ewas_bacon_toptable.csv.gz"))
setorder(tt, bacon.p)
all.cpg <- tt$probe
sig.cpg <- tt$probe[1:1000]
cat("universe:", length(all.cpg), "| query:", length(sig.cpg), "\n")
suppressPackageStartupMessages(library(ggplot2))
## ---------------- 3. KYCG ----------------
## KYCG lives in `knowYourCG` in recent Bioconductor releases, but in sesame
## <= 1.16 `testEnrichment` is exported by sesame itself. Accept either route,
## otherwise the check reports "unavailable" on a machine where it works fine.
kycg_ok <- requireNamespace("knowYourCG", quietly = TRUE) ||
           "testEnrichment" %in% getNamespaceExports("sesame")
cat("\nKYCG available:", kycg_ok,
    "| via:", if (requireNamespace("knowYourCG", quietly = TRUE)) "knowYourCG" else "sesame", "\n")
if (kycg_ok) {
  if (requireNamespace("knowYourCG", quietly = TRUE)) {
    suppressPackageStartupMessages(library(knowYourCG))
  } else {
    suppressPackageStartupMessages(library(sesame))
  }
  suppressPackageStartupMessages({ library(sesameData); library(dbplyr); library(dplyr) })
  ## AnnotationHub 3.6.0 x dbplyr 2.5.0: collect(Inf) passed positionally.
  source("../patch_annotationhub_collect.R"); invisible(patch_ah_collect())
  suppressMessages(sesameDataCache())
  dbs <- c("KYCG.EPIC.CGI.20210713", "KYCG.EPIC.chromHMM.20211020",
           "KYCG.EPIC.TFBSconsensus.20211013", "KYCG.EPIC.HMconsensus.20211013")
  res <- try(testEnrichment(query = sig.cpg, databases = dbs,
                            universe = all.cpg, platform = "EPIC"), silent = TRUE)
  if (inherits(res, "try-error")) {
    cat("KYCG ERROR:", attr(res, "condition")$message, "\n")
  } else {
    res <- as.data.table(res)[order(FDR)]
    saveRDS(res, file.path(outd, "08_kycg.rds"))
    fwrite(res, file.path(outd, "08_kycg_results.csv"))
    cat("KYCG rows:", nrow(res), "| FDR<0.05:", sum(res$FDR < 0.05, na.rm = TRUE), "\n")
    ## The label column has been named `dbname` in some sesame versions and
    ## `group`/`db` in others; resolve it rather than assuming.
    labcol <- intersect(c("dbname", "group", "db", "gene_name"), names(res))[1]
    if (is.na(labcol)) labcol <- names(res)[1]
    cat("label column:", labcol, "| columns:", paste(names(res), collapse = ", "), "\n")
    print(head(res, 10))
    kp <- head(res[is.finite(estimate)], 12)
    kp[, lab := factor(get(labcol), levels = rev(get(labcol)))]
    ## Color by database group: the top hits are all FDR < 0.05, so a
    ## significance legend would carry no information here.
    kp[, grp := sub("^KYCG\\.EPIC\\.", "", sub("\\.[0-9]+$", "", group))]
    pk <- ggplot(kp, aes(estimate, lab, fill = grp)) +
      geom_col() +
      scale_fill_manual(values = c(chromHMM = "#4C72B0", TFBSconsensus = "#C44E52",
                                   HMconsensus = "#DD8452", CGI = "#55A868"),
                        name = NULL) +
      labs(x = expression(log[2]~"odds ratio"), y = NULL,
           title = "KYCG enrichment, top 1,000 CpGs by BACON p",
           subtitle = "all bars shown are FDR < 0.05") +
      theme_minimal(base_size = 10) + theme(legend.position = "top")
    ggsave(file.path(outd, "08_kycg_enrichment.png"), plot = pk, width = 8.5, height = 5, dpi = 200)
    cat("\n--- per-group summary ---\n")
    print(res[, .(n_tested = .N, n_sig = sum(FDR < 0.05, na.rm = TRUE),
                  best = get(labcol)[which.min(FDR)],
                  best_FDR = min(FDR, na.rm = TRUE),
                  best_log2OR = estimate[which.min(FDR)]), by = group])
    cat("\n--- CGI group, all rows ---\n")
    print(res[grepl("CGI", group), .(dbname, estimate, p.value, FDR, nQ, nD, overlap)])
    cat("\n--- chromHMM group, all rows ---\n")
    print(res[grepl("chromHMM", group), .(dbname, estimate, p.value, FDR, overlap)])
  }
}
cat("elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
