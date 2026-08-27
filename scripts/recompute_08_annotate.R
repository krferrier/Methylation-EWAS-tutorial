## Recompute chapter 08 Part 1: annotate the regenerated BACON toptable.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setwd("repo"); source("_setup.R")
t0 <- Sys.time()
outd <- data_path("08_annotation")
dir.create(outd, showWarnings = FALSE)

tt <- fread(data_path("06_ewas_bacon_toptable.csv.gz"))
cat("toptable:", nrow(tt), "rows |", paste(names(tt), collapse=", "), "\n")

## ---- gencode v41: gene names, biotypes, signed distance to TSS, coordinates ----
v41 <- fread(cmd = "zcat ../masks/EPIC.hg38.manifest.gencode.v41.tsv.gz")
nearest_tss <- function(x) {
  v <- suppressWarnings(as.numeric(strsplit(x, ";")[[1]]))
  if (all(is.na(v))) return(NA_real_)
  v[which.min(abs(v))]
}
v41[, distTSS := vapply(distToTSS, nearest_tss, numeric(1))]

## ---- v36 carries the CGI relation ----
v36 <- fread(cmd = "zcat ../masks/EPIC.hg38.manifest.gencode.v36.tsv.gz",
             select = c("probeID","CGIposition"))
v36[, island := fifelse(is.na(CGIposition) | CGIposition == "", "OpenSea", CGIposition)]

## v41 already carries hg38 coordinates keyed on probeID, so it is both the
## coordinate and the gene-model source (the standalone coord file is unkeyed).
annotated <- tt |>
  merge(v41[, .(probe = probeID, CpG_chrm, CpG_beg, CpG_end,
                genesUniq, transcriptTypes, distTSS)], by = "probe", all.x = TRUE) |>
  merge(v36[, .(probe = probeID, island)], by = "probe", all.x = TRUE)

annotated[, feature := fifelse(is.na(genesUniq) | genesUniq %in% c("", "NA"), "Intergenic",
                        fifelse(!is.na(distTSS) & abs(distTSS) <= 1500, "Promoter", "Gene body"))]
setorder(annotated, bacon.p)

n_tot   <- nrow(annotated)
n_coord <- sum(!is.na(annotated$CpG_chrm))
n_gene  <- sum(!(is.na(annotated$genesUniq) | annotated$genesUniq %in% c("","NA")))
cat("annotated:", n_tot, "| with coords:", n_coord,
    "| mapped to a gene:", n_gene, sprintf("(%.1f%%)", 100*n_gene/n_tot), "\n")

## ---- BIOS eQTM ----
eqtm <- fread(file.path(outd, "bios_eqtm_hgnc_annotated.tsv"))
eq <- eqtm[HGNCName_GRCh38 != "" & !is.na(HGNCName_GRCh38),
           .(BIOS_eQTM_genes  = paste(sort(unique(HGNCName_GRCh38)), collapse = ";"),
             BIOS_eQTM_minFDR = min(as.numeric(FDR), na.rm = TRUE)),
           by = .(probe = SNPName)]
annotated <- merge(annotated, eq, by = "probe", all.x = TRUE)
setorder(annotated, bacon.p)
n_eqtm <- sum(!is.na(annotated$BIOS_eQTM_genes))
cat("eQTM-annotated:", n_eqtm, sprintf("(%.1f%%)", 100*n_eqtm/n_tot), "\n")

fwrite(annotated, file.path(outd, "PTSD_ewas_annotated_zhou.csv.gz"))
saveRDS(annotated, file.path(outd, "annotated.rds"))

## ---- BED for comb-p ----
bed <- annotated[!is.na(CpG_chrm),
                 .(`#chrom` = CpG_chrm, start = CpG_beg, end = CpG_end,
                   pvals = bacon.p, cpgid = probe)]
setorder(bed, `#chrom`, start)
fwrite(bed, file.path(outd, "PTSD_ewas_annotated_results.bed"), sep = "\t")
cat("bed rows:", nrow(bed), "\n")

## ---- top-10 display table ----
t10 <- head(annotated, 10)[, .(
  probe, chr = CpG_chrm, pos = CpG_beg,
  gene = fifelse(is.na(genesUniq) | genesUniq == "", "-", genesUniq),
  distTSS = round(distTSS), feature, island,
  delta_beta = round(100 * delta_beta, 2),
  bacon_p = signif(bacon.p, 3), bacon_FDR = signif(bacon.adj.P, 3))]
setnames(t10, c("probe","chr","pos","gene","distToTSS","feature","island",
                "Δβ (pp)","BACON p","BACON FDR"))
fwrite(t10, file.path(outd, "top10_display.csv"))
print(t10)

## ---- eQTM display table ----
eqd <- annotated[!is.na(BIOS_eQTM_genes)][1:10]
eqd[, nearest := fifelse(mapply(function(g, e) {
      if (is.na(g) || g == "") return(FALSE)
      any(strsplit(e, ";")[[1]] %in% strsplit(g, ";")[[1]])
    }, genesUniq, BIOS_eQTM_genes), "yes", "no")]
eqdisp <- eqd[, .(probe, nearest_gene = fifelse(is.na(genesUniq)|genesUniq=="", "-", genesUniq),
                  eQTM_gene = BIOS_eQTM_genes, `nearest?` = nearest,
                  eQTM_FDR = signif(BIOS_eQTM_minFDR, 3),
                  bacon_p = signif(bacon.p, 3))]
fwrite(eqdisp, file.path(outd, "eqtm_display.csv"))
fwrite(annotated[!is.na(BIOS_eQTM_genes)][1:200], file.path(outd, "eqtm_top_hits.csv"))
print(head(eqdisp, 4))

## ---- feature / island distribution figures ----
TOPN <- 1000
top <- head(annotated, TOPN)
mk <- function(col, lev) {
  a <- data.table(set = "Top 1,000", grp = factor(top[[col]], levels = lev))
  b <- data.table(set = "All tested", grp = factor(annotated[[col]], levels = lev))
  d <- rbind(a, b)[!is.na(grp), .N, by = .(set, grp)]
  d[, pct := 100 * N / sum(N), by = set]
  d[, set := factor(set, levels = c("Top 1,000", "All tested"))]
  d
}
fd <- mk("feature", c("Promoter","Gene body","Intergenic"))
isl_lev <- c("Island","N_Shore","S_Shore","N_Shelf","S_Shelf","OpenSea")
id <- mk("island", isl_lev)

pf <- ggplot(fd, aes(grp, pct, fill = set)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_dodge(width = 0.9),
            vjust = -0.35, size = 3) +
  scale_fill_manual(values = c("Top 1,000" = "#C44E52", "All tested" = "#4C72B0"), name = NULL) +
  labs(x = NULL, y = "% of CpGs", title = "Gene-feature context") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(outd, "08_feature_distribution.png"), plot = pf, width = 8, height = 4.2, dpi = 200)

pi <- ggplot(id, aes(grp, pct, fill = set)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_dodge(width = 0.9),
            vjust = -0.35, size = 2.9) +
  scale_fill_manual(values = c("Top 1,000" = "#C44E52", "All tested" = "#4C72B0"), name = NULL) +
  labs(x = NULL, y = "% of CpGs", title = "CpG-island context") +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(outd, "08_island_distribution.png"), plot = pi, width = 8.5, height = 4.2, dpi = 200)

cat("\n--- feature distribution ---\n"); print(dcast(fd, grp ~ set, value.var = "pct"))
cat("\n--- island distribution ---\n"); print(dcast(id, grp ~ set, value.var = "pct"))
bonf <- 0.05 / n_tot
cat("\nBonferroni threshold:", signif(bonf, 3),
    "| smallest bacon p:", signif(min(annotated$bacon.p), 3), "\n")
cat("n Bonferroni:", sum(annotated$bacon.p < bonf),
    "| n FDR<0.05:", sum(annotated$bacon.adj.P < 0.05),
    "| n p<1e-5:", sum(annotated$bacon.p < 1e-5), "\n")
cat("n within 1.5kb of TSS in top 10:", sum(head(annotated,10)$feature == "Promoter"), "\n")
cat("top 10 probes:", paste(head(annotated$probe, 10), collapse = ", "), "\n")
cat("rank of cg16340178:", which(annotated$probe == "cg16340178"), "\n")
cat("elapsed:", round(difftime(Sys.time(), t0, units="mins"),2), "min\n")
