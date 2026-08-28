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

# ---------------------------------------------------------------------------
# Figure palette. These are the same hexes the site theme (theme.scss) uses, so
# plots and prose share one visual language. Teal/plum reads as a difference in
# both hue and lightness, which keeps the contrast legible in grayscale and
# under the common forms of color-vision deficiency.
# ---------------------------------------------------------------------------
ewas_col <- c(
  teal_dark  = "#0F3D43",
  teal       = "#1A6B75",
  teal_light = "#2A8F9B",
  sand       = "#B8873F",
  plum       = "#8C3A4A",
  gray_dark  = "#3A362F",
  gray       = "#6E675B",
  gray_light = "#B5AFA4"
)

# Palette for phenotype groups. Case/Control and Female/Male each get one
# teal and one warm contrast so the two groupings never look interchangeable.
pheno_pal <- c(Control = "#1A6B75", Case = "#8C3A4A",
               Female  = "#B8873F", Male = "#1A6B75",
               neg     = "#1A6B75", pos  = "#8C3A4A")

# Body typeface for figures, matching the site. The four Source Sans 3 faces
# ship in assets/fonts/otf/ so a fresh clone renders the same figures without
# installing anything system-wide; systemfonts::register_font() makes them
# visible to R's graphics devices. If the files are absent (or systemfonts is
# too old to register), .ewas_family stays "" and every plot silently falls
# back to the device default — nothing breaks, the type just changes.
.ewas_family <- local({
  fam  <- "Source Sans 3"
  dirs <- c("assets/fonts/otf", "../assets/fonts/otf")
  dir  <- dirs[dir.exists(dirs)][1]

  # Already visible to R (installed system-wide)? Then use it as-is.
  seen <- tryCatch(fam %in% systemfonts::system_fonts()$family,
                   error = function(e) FALSE)
  if (isTRUE(seen)) return(fam)

  if (is.na(dir)) return("")
  f <- file.path(dir, paste0("SourceSans3-",
                             c("Regular", "Semibold", "It", "SemiboldIt"),
                             ".otf"))
  if (!all(file.exists(f))) return("")

  ok <- tryCatch({
    systemfonts::register_font(fam, plain = f[1], bold = f[2],
                               italic = f[3], bolditalic = f[4])
    TRUE
  }, error = function(e) FALSE)

  if (isTRUE(ok)) fam else ""
})

# Consistent ggplot theme for the tutorial
theme_ewas <- function(base_size = 12, base_family = .ewas_family) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E4E1DA", linewidth = 0.35),
      panel.border     = element_rect(colour = "#B5AFA4", fill = NA, linewidth = 0.5),
      plot.title       = element_text(face = "bold", colour = "#0F3D43",
                                     size = rel(1.05)),
      plot.subtitle    = element_text(colour = "#6E675B", size = rel(0.92)),
      plot.caption     = element_text(colour = "#6E675B", size = rel(0.8),
                                      hjust = 0),
      axis.title       = element_text(colour = "#3A362F"),
      axis.text        = element_text(colour = "#6E675B"),
      strip.background = element_rect(fill = "#F2F0EC", colour = "#B5AFA4",
                                      linewidth = 0.5),
      strip.text       = element_text(colour = "#0F3D43", face = "bold",
                                      size = rel(0.9)),
      legend.position  = "right",
      legend.title     = element_text(colour = "#3A362F"),
      legend.key       = element_blank()
    )
}
theme_set(theme_ewas())

data_path <- function(f) file.path("data", f)
