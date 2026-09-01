suppressMessages({library(data.table)})
setwd("repo")
ord <- fread("../yame_data/EPIC/EPIC.ordering.tsv.gz", select="Probe_ID")
crd <- fread("../yame_data/EPIC/EPIC.hg38.coord.tsv.gz")
map <- data.table(probe=ord$Probe_ID, chrm=crd$CpG_chrm, beg=crd$CpG_beg)

e  <- readRDS("data/06_ewas.rds")
tt <- as.data.table(e$tt)
reg <- merge(tt, map, by="probe")[chrm=="chr6" & beg>=28633493 & beg<=28633701]
setorder(reg, beg)
stopifnot(nrow(reg)==11)

# GENCODE feature / CGI relation / gene / distTSS from the chapter 08 annotation layer
ann <- readRDS("data/08_annotation/annotated.rds")
ann <- as.data.table(ann)[, .(probe, feature, CGIposition = island, genesUniq,
                              transcriptTypes, distTSS, CpG_end,
                              BIOS_eQTM_genes)]
reg <- merge(reg, ann, by="probe", all.x=TRUE); setorder(reg, beg)
stopifnot(!any(is.na(reg$feature)), !any(is.na(reg$CGIposition)))
cat("genes in region:", paste(unique(reg$genesUniq), collapse=" | "), "\n")
cat("distTSS range:", paste(range(reg$distTSS, na.rm=TRUE), collapse=" .. "), "\n")

disp <- reg[, .(CpG = probe,
                `hg38 position (chr6)` = formatC(beg, big.mark=",", format="d"),
                `GENCODE feature` = sub(" \\(.*", "", feature),
                `CGI context` = CGIposition,
                `Δβ (case − control)` = round(delta_beta, 3),
                `BACON p` = signif(bacon.p, 3))]
fwrite(disp, "data/08_annotation/dmr_cpg_display.csv")
print(disp, row.names=FALSE)

full <- reg[, .(probe, CpG_beg = beg, CpG_end, genesUniq, transcriptTypes, distTSS,
                feature, island = CGIposition, BIOS_eQTM_genes,
                delta_beta = round(delta_beta, 4),
                bacon.es   = round(bacon.es, 3),
                bacon.p    = signif(bacon.p, 3))]
fwrite(full, "data/08_annotation/dmr_cpg_annotation.csv")

cat("\n=== anchors ===\n")
cat("span_bp:", max(reg$beg)-min(reg$beg)+2, "\n")
cat("n_cpgs:", nrow(reg), "\n")
cat("best_raw_bacon_p:", signif(min(reg$bacon.p),3), "\n")
cat("dbeta_pp_min:", round(100*min(reg$delta_beta),2), " dbeta_pp_max:", round(100*max(reg$delta_beta),2), "\n")
cat("neglog10_min:", round(min(-log10(reg$bacon.p)),2), " max:", round(max(-log10(reg$bacon.p)),2), "\n")
