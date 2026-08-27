kg <- read.delim("kegg_hsa_pathways.tsv", header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("id", "name"))
kg$id5  <- sub("^hsa", "", kg$id)
kg$name <- sub(" - Homo sapiens .human.$", "", kg$name)
lut <- setNames(kg$name, kg$id5)

d <- "repo/data/08_annotation"

g <- read.csv(file.path(d, "08_methylgsa_glm_top.csv"), stringsAsFactors = FALSE,
              colClasses = c(ID = "character"))
g$Description <- lut[g$ID]
cat("glm unmatched:", sum(is.na(g$Description)), "\n")
write.csv(g, file.path(d, "08_methylgsa_glm_top.csv"), row.names = FALSE)
print(g[, c("ID", "Description", "Size", "pvalue", "padj")], row.names = FALSE)

r <- read.csv(file.path(d, "08_methylgsa_rra_top.csv"), stringsAsFactors = FALSE,
              colClasses = c(KEGG = "character"))
r$Pathway <- lut[r$KEGG]
cat("\nrra unmatched:", sum(is.na(r$Pathway)), "\n")
write.csv(r, file.path(d, "08_methylgsa_rra_top.csv"), row.names = FALSE)
print(r[, c("KEGG", "Pathway", "Size", "NES", "pvalue", "padj")], row.names = FALSE)
cat("\nrra n padj<0.05:", sum(r$padj < 0.05, na.rm = TRUE),
    "| min padj:", min(r$padj, na.rm = TRUE), "\n")
cat("glm n padj<0.05:", sum(g$padj < 0.05, na.rm = TRUE),
    "| min padj:", min(g$padj, na.rm = TRUE), "\n")
