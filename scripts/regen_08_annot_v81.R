suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setwd("repo"); source("_setup.R")
t0 <- Sys.time()
outd <- data_path("08_annotation")

## authoritative v8.1 toptable
e  <- readRDS(data_path("06_ewas.rds"))
tt <- as.data.table(e$tt)
cat("toptable:", nrow(tt), "rows | n_tested field:", e$n_tested, "\n")
stopifnot(nrow(tt) == 756251L)
## re-emit the documented ch08 input, then read it back so the annotation runs
## through exactly the code path chapter 08 documents (fread of the csv.gz).
write.csv(tt, gzfile(data_path("06_ewas_bacon_toptable.csv.gz")), row.names = FALSE)
tt <- fread(data_path("06_ewas_bacon_toptable.csv.gz"))
cat("re-read from csv.gz:", nrow(tt), "rows |", paste(names(tt), collapse=", "), "\n")

v41 <- fread(cmd = "zcat ../masks/EPIC.hg38.manifest.gencode.v41.tsv.gz")
nearest_tss <- function(x) {
  v <- suppressWarnings(as.numeric(strsplit(x, ";")[[1]]))
  if (all(is.na(v))) return(NA_real_)
  v[which.min(abs(v))]
}
v41[, distTSS := vapply(distToTSS, nearest_tss, numeric(1))]

v36 <- fread(cmd = "zcat ../masks/EPIC.hg38.manifest.gencode.v36.tsv.gz",
             select = c("probeID","CGIposition"))
v36[, island := fifelse(is.na(CGIposition) | CGIposition == "", "OpenSea", CGIposition)]

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
cat("annotated:", n_tot, "| coords:", n_coord, "| gene-mapped:", n_gene,
    sprintf("(%.1f%%)", 100*n_gene/n_tot), "\n")

eqtm <- fread(file.path(outd, "bios_eqtm_hgnc_annotated.tsv"))
eq <- eqtm[HGNCName_GRCh38 != "" & !is.na(HGNCName_GRCh38),
           .(BIOS_eQTM_genes  = paste(sort(unique(HGNCName_GRCh38)), collapse = ";"),
             BIOS_eQTM_minFDR = min(as.numeric(FDR), na.rm = TRUE)),
           by = .(probe = SNPName)]
annotated <- merge(annotated, eq, by = "probe", all.x = TRUE)
setorder(annotated, bacon.p)
n_eqtm <- sum(!is.na(annotated$BIOS_eQTM_genes))
cat("eQTM-annotated:", n_eqtm, sprintf("(%.1f%%)", 100*n_eqtm/n_tot), "\n")
cat("ncol:", ncol(annotated), "\n")

fwrite(annotated, file.path(outd, "PTSD_ewas_annotated_zhou.csv.gz"))
saveRDS(annotated, file.path(outd, "annotated.rds"))

## BED written to a COMPARISON name; do not clobber the comb-p input yet
bed <- annotated[!is.na(CpG_chrm),
                 .(`#chrom` = CpG_chrm, start = CpG_beg, end = CpG_end,
                   pvals = bacon.p, cpgid = probe)]
setorder(bed, `#chrom`, start)
fwrite(bed, file.path(outd, "PTSD_ewas_bed_from_v41.bed"), sep = "\t")
cat("bed rows (v41 route):", nrow(bed), "\n")

## display tables -> comparison names
t10 <- head(annotated, 10)[, .(
  probe, chr = CpG_chrm, pos = CpG_beg,
  gene = fifelse(is.na(genesUniq) | genesUniq == "", "-", genesUniq),
  distTSS = round(distTSS), feature, island,
  delta_beta = round(100 * delta_beta, 2),
  bacon_p = signif(bacon.p, 3), bacon_FDR = signif(bacon.adj.P, 3))]
setnames(t10, c("probe","chr","pos","gene","distToTSS","feature","island",
                "\u0394\u03b2 (pp)","BACON p","BACON FDR"))
fwrite(t10, file.path(outd, "top10_display.csv"))
print(t10)

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
cat("\n--- feature distribution ---\n"); print(dcast(fd, grp ~ set, value.var = "pct"))
cat("\n--- island distribution ---\n"); print(dcast(id, grp ~ set, value.var = "pct"))

bonf <- 0.05 / n_tot
cat("\nBonferroni:", signif(bonf, 4), "| min bacon p:", signif(min(annotated$bacon.p), 4), "\n")
cat("n Bonf:", sum(annotated$bacon.p < bonf),
    "| n FDR<0.05:", sum(annotated$bacon.adj.P < 0.05),
    "| n p<1e-5:", sum(annotated$bacon.p < 1e-5), "\n")
cat("promoters in top10:", sum(head(annotated,10)$feature == "Promoter"), "\n")
cat("rank cg16340178:", which(annotated$probe == "cg16340178"), "\n")
cat("elapsed:", round(difftime(Sys.time(), t0, units="mins"),2), "min\n")
