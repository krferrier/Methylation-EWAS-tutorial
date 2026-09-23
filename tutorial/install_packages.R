## ---------------------------------------------------------------------------
## Install every R package this tutorial needs, without conda.
##
## Run once from the `tutorial/` directory:
##     Rscript install_packages.R
## or paste it into an R console.
##
## Requires R 4.4 or newer -- install the current release. BiocManager picks the
## Bioconductor release that matches your R (R 4.6 -> 3.23, R 4.5 -> 3.22), and
## all 28 package names resolve in every release from 3.19 on (knowYourCG first
## appears in 3.19, which is why R 4.3 and Bioconductor 3.18 are not enough). Nothing here
## needs a version pin: install the current matrixStats like everything else.
##
## The published numbers and the Zenodo checkpoints came from R 4.2.3 /
## Bioconductor 3.16, so a value may differ in its last digits on a newer
## stack. The code is unchanged. See SESSIONINFO.md.
##
## Budget 20-40 minutes on a first run. On Windows and macOS almost everything
## arrives as a pre-built binary -- there are just a lot of packages. On Linux
## they build from source, which is slower.
## ---------------------------------------------------------------------------

options(timeout = 1200)          # some Bioconductor tarballs are large
options(Ncpus = max(1L, parallel::detectCores() - 1L))   # parallel compiles

## --- 0. check the R version -------------------------------------------------
rv <- getRversion()
if (rv < "4.4") {
  stop("This tutorial needs R 4.4 or newer; you are on ", rv, ".\n",
       "knowYourCG, used in chapter 08, first appears in Bioconductor 3.19,\n",
       "which needs R 4.4. (Older still, R 4.2 selects Bioconductor 3.16, where\n",
       "MatrixGenerics passes `useNames = NA` -- an argument matrixStats made\n",
       "defunct in 1.2.0 -- so chapters 01-02 fail.) Install the current\n",
       "R from https://cran.r-project.org/ and re-run this script.",
       call. = FALSE)
}

## --- 1. bootstrap BiocManager -----------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_ver <- BiocManager::version()
cat("R", as.character(rv), "-> Bioconductor", as.character(bioc_ver), "\n\n")

## MatrixGenerics 1.13.1 switched every `useNames` default from NA to TRUE,
## which shipped in Bioconductor 3.18. Below that release, a current
## matrixStats breaks detectionP() and preprocessFunnorm(). knowYourCG needs
## 3.19 or later.
if (bioc_ver < "3.19") {
  stop("BiocManager selected Bioconductor ", bioc_ver, ", which predates\n",
       "knowYourCG (first released in 3.19). Update BiocManager with\n",
       "install.packages(\"BiocManager\"), or move to a newer R.", call. = FALSE)
}

## --- 2. the packages --------------------------------------------------------
cran <- c("data.table", "ggplot2", "knitr", "DT",
          "ggrepel", "ggtext",   # figure labels (chapters 06 and 08)
          "fst", "jsonlite")     # pipeline input file (07); UCSC API query (08)

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
  "missMethyl", "methylGSA", "knowYourCG", "ENmix"
)

## `upgrade = "never"` keeps this from rebuilding packages you already have
## while it resolves dependencies -- it makes a re-run cheap, not just a
## first run.
BiocManager::install(c(cran, bioc), ask = FALSE, update = FALSE,
                     upgrade = "never")

## --- 3. report --------------------------------------------------------------
cat("\n--- installed versions ---\n")
for (p in c(cran, bioc, "matrixStats", "MatrixGenerics")) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")
  cat(sprintf("  %-48s %s\n", p, v))
}

missing <- Filter(function(p) !requireNamespace(p, quietly = TRUE),
                  c(cran, bioc))
if (length(missing)) {
  cat("\nFAILED to install:", paste(missing, collapse = ", "), "\n")
  cat("Most first-time failures are missing system libraries; the error text\n",
      "above usually names the one to install.\n")
} else {
  cat("\nAll packages installed.\n")
}
