suppressPackageStartupMessages(library(data.table))
d <- "repo/data/08_annotation"

## ---- 1. three-method GSA comparison, with the winning term and its FDR ----
g  <- readRDS(file.path(d, "08_gometh.rds"))
kg <- as.data.table(g$kegg, keep.rownames = "ID")[order(P.DE)]
gl <- fread(file.path(d, "08_methylgsa_glm_top.csv"), colClasses = c(ID = "character"))
rr <- fread(file.path(d, "08_methylgsa_rra_top.csv"), colClasses = c(KEGG = "character"))

cmp <- data.table(
  Method = c("gometh (Wallenius hypergeometric)",
             "methylglm (logistic regression)",
             "methylRRA-GSEA (rank aggregation)"),
  Input  = c("top 1,000 CpGs vs universe",
             "full named p-value vector",
             "full named p-value vector"),
  `Bias model` = c("probes per gene (Wallenius)",
                   "probe count as covariate",
                   "gene-level RRA score"),
  `KEGG FDR<0.05` = c(sum(kg$FDR < 0.05), sum(gl$padj < 0.05), sum(rr$padj < 0.05)),
  `Top KEGG pathway` = c(kg$Description[1], gl$Description[1], rr$Pathway[1]),
  `Its nominal p` = signif(c(kg$P.DE[1], gl$pvalue[1], rr$pvalue[1]), 3),
  `Its FDR`       = signif(c(kg$FDR[1],  gl$padj[1],   rr$padj[1]),   3)
)
fwrite(cmp, file.path(d, "08_gsa_method_comparison.csv"))
cat("--- GSA comparison ---\n"); print(cmp)

## ---- 2. KYCG per-database summary ----
k <- fread(file.path(d, "08_kycg_results.csv"))
setnames(k, c("group", "dbname"), c("db", "feature"))
k[, kb := sub("^KYCG\\.EPIC\\.", "", sub("\\.[0-9]+$", "", db))]

ksum <- k[, .(
  `Features tested` = .N,
  `Nominal p<0.05`  = sum(p.value < 0.05, na.rm = TRUE),
  `FDR<0.05`        = sum(FDR < 0.05, na.rm = TRUE),
  `Top feature`     = feature[which.min(FDR)],
  `Top log2(OR)`    = round(estimate[which.min(FDR)], 2),
  `Top FDR`         = signif(min(FDR, na.rm = TRUE), 3)
), by = .(Knowledgebase = kb)][order(-`FDR<0.05`)]
fwrite(ksum, file.path(d, "08_kycg_db_summary.csv"))
cat("\n--- KYCG per-knowledgebase ---\n"); print(ksum)

## ---- 3. KYCG top hits (all FDR<0.05, capped for display) ----
ktop <- k[FDR < 0.05][order(FDR, p.value),
  .(Knowledgebase = kb, Feature = feature,
    `log2(OR)` = round(estimate, 2), `Query CpGs` = overlap,
    `DB size` = nD, `p` = signif(p.value, 3), FDR = signif(FDR, 3))]
fwrite(ktop, file.path(d, "08_kycg_top.csv"))
cat("\n--- KYCG top 15 of", nrow(ktop), "---\n"); print(head(ktop, 15))

## ---- 4. CGI relation, the full five rows (now non-significant) ----
kcgi <- k[kb == "CGI"][order(FDR),
  .(Feature = feature, `log2(OR)` = round(estimate, 2),
    `Query CpGs` = overlap, `DB size` = nD,
    p = signif(p.value, 3), FDR = signif(FDR, 3))]
fwrite(kcgi, file.path(d, "08_kycg_cgi.csv"))
cat("\n--- KYCG CGI relation ---\n"); print(kcgi)

cat("\ntotals: features=", nrow(k), " nominal=", sum(k$p.value < 0.05),
    " fdr=", sum(k$FDR < 0.05), " expected_nominal=", round(0.05 * nrow(k)), "\n", sep = "")
