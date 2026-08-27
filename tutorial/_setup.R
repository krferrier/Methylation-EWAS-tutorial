# Shared setup sourced by every notebook.
# Ensures the writable library holding the EPIC data packages is on the path.
# Notebooks render from tutorial/, so also look one level up at the workspace libs.
.cand_libs <- c(".r-libs/methyl", "../.r-libs/methyl",
                "/home/krferrier/.claude-science/r-libs/28e92aec-1a42-49fe-9453-8f7c0cb91a8e/methyl",
                Sys.getenv("EWAS_EPIC_LIB", unset = NA))
.cand_libs <- .cand_libs[!is.na(.cand_libs) & nzchar(.cand_libs)]
.cand_libs <- .cand_libs[dir.exists(.cand_libs)]
if (length(.cand_libs)) .libPaths(c(.cand_libs, .libPaths()))

suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(ggplot2)
  library(data.table)
})

# Consistent ggplot theme for the tutorial
theme_ewas <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "right")
}
theme_set(theme_ewas())

# Palette for phenotype groups
pheno_pal <- c(Control = "#4C72B0", Case = "#C44E52",
               Female = "#DD8452", Male = "#4C72B0",
               neg = "#4C72B0", pos = "#C44E52")

data_path <- function(f) file.path("data", f)
