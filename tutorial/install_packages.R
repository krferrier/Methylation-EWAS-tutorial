## ---------------------------------------------------------------------------
## Install every R package this tutorial needs, without conda.
##
## Run once from the `tutorial/` directory:
##     Rscript install_packages.R
## or paste it into an R console.
##
## Requires R 4.2.x. BiocManager maps that to Bioconductor 3.16 automatically,
## which is the release that produced every checkpoint in the Zenodo record.
##
## Expect 20-40 minutes on a first run: several of these compile from source.
## ---------------------------------------------------------------------------

options(timeout = 1200)          # some Bioconductor tarballs are large
options(Ncpus = max(1L, parallel::detectCores() - 1L))   # parallel compiles

## --- 0. sanity check the R version -----------------------------------------
rv <- getRversion()
if (rv < "4.2" || rv >= "4.3") {
  warning(
    "This tutorial was built on R 4.2.x (Bioconductor 3.16). You are on ", rv,
    ".\nIt may still work, but package versions will differ from the ones that\n",
    "produced the published numbers.", call. = FALSE, immediate. = TRUE)
}

## --- 1. bootstrap the installers -------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
cat("Bioconductor release:", as.character(BiocManager::version()), "\n\n")

## --- 2. the matrixStats pin, FIRST -----------------------------------------
## matrixStats made `useNames = NA` defunct in 1.2.0, but Bioconductor 3.16's
## MatrixGenerics still passes it -- so a current matrixStats breaks
## detectionP() and preprocessFunnorm(). Install the pinned version before
## anything can pull in a newer one as a dependency.
need_pin <- !requireNamespace("matrixStats", quietly = TRUE) ||
  packageVersion("matrixStats") >= "1.2.0"
if (need_pin) {
  cat("Installing matrixStats 1.0.0 (pinned -- see SESSIONINFO.md)\n")
  ok <- tryCatch({
    remotes::install_version("matrixStats", version = "1.0.0",
                             repos = "https://cloud.r-project.org", upgrade = "never")
    TRUE
  }, error = function(e) { cat("  install failed:", conditionMessage(e), "\n"); FALSE })

  ## 1.0.0 is source-only, so this step is the one that needs a compiler.
  ## If it failed, the 1.3.0 binary plus a hidden option is a working fallback:
  ##   install.packages("matrixStats")
  ##   options(matrixStats.useNames.NA = "deprecated")   # 1.3.0 only; gone in 1.5.0
  ## Put that option in your .Rprofile so it applies to every session.
  if (!ok) {
    cat("\n  Could not build matrixStats 1.0.0 (it is source-only and needs a\n",
        "  compiler: Rtools on Windows, xcode-select --install on macOS).\n",
        "  Alternative -- use the binary and relax the check:\n",
        "      install.packages(\"matrixStats\")\n",
        "      options(matrixStats.useNames.NA = \"deprecated\")\n",
        "  See the Setup chapter for details.\n\n", sep = "")
  }
}

## --- 3. everything else -----------------------------------------------------
cran <- c("data.table", "ggplot2", "knitr", "DT")

bioc <- c(
  # array I/O and preprocessing
  "minfi", "illuminaio", "wateRmelon", "GEOquery", "Biobase",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  # cell composition
  "FlowSorted.Blood.EPIC", "genefilter",
  # batch effects, association testing, inflation
  "sva", "limma", "bacon",
  # annotation and enrichment
  "sesame", "sesameData", "GenomicRanges", "rtracklayer",
  "missMethyl", "methylGSA"
)

## `upgrade = "never"` is what protects the matrixStats pin from being
## silently bumped while resolving these.
BiocManager::install(c(cran, bioc), ask = FALSE, update = FALSE,
                     upgrade = "never")

## --- 4. report --------------------------------------------------------------
cat("\n--- installed versions ---\n")
for (p in c("matrixStats", cran, bioc)) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")
  cat(sprintf("  %-48s %s\n", p, v))
}

missing <- Filter(function(p) !requireNamespace(p, quietly = TRUE),
                  c("matrixStats", cran, bioc))
if (length(missing)) {
  cat("\nFAILED to install:", paste(missing, collapse = ", "), "\n")
  cat("Most first-time failures are missing system libraries; the error text\n",
      "above usually names the one to install.\n")
} else {
  cat("\nAll packages installed.\n")
}

ms <- as.character(packageVersion("matrixStats"))
if (utils::compareVersion(ms, "1.2.0") >= 0) {
  cat("\nWARNING: matrixStats is", ms, "-- chapters 01-02 will fail with\n",
      "'useNames = NA is defunct'. Re-run step 2 above.\n")
}
