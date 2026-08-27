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

## ---------------- 1. missMethyl::gometh ----------------
suppressPackageStartupMessages(library(missMethyl))

## UPSTREAM BUG PATCH (missMethyl 1.32.0): .getKEGG() has TWO independent bugs,
## and both must be fixed or gometh(collection = "KEGG") fails.
##
## (1) Species code case. It calls getGeneKEGGLinks(species.KEGG = "hsa") but
##     getKEGGPathwayNames(species.KEGG = "Hsa"). The KEGG REST API is
##     case-sensitive and answers HTTP 400 for "Hsa".
##
## (2) Pathway-ID prefix mismatch. getGeneKEGGLinks() returns PathwayID WITH the
##     KEGG qualifier ("path:hsa00010"), while getKEGGPathwayNames(remove.
##     qualifier = TRUE) returns it WITHOUT ("hsa00010"). The subsequent
##     merge(by = "PathwayID") therefore matches ZERO rows, idList comes back
##     empty, and gsameth() dies with "subscript out of bounds" inside `%in%`.
##     Verified here: 0 rows as shipped vs 39,576 rows / 372 pathways once the
##     "path:" prefix is stripped. Fixing (1) alone only converts the HTTP 400
##     into this silent empty join.
local({
  ns <- asNamespace("missMethyl")
  f  <- get(".getKEGG", envir = ns)
  txt <- paste(deparse(body(f)), collapse = "\n")
  ## (1) species code
  txt <- gsub('species.KEGG = "Hsa"', 'species.KEGG = "hsa"', txt, fixed = TRUE)
  ## (2) strip the "path:" qualifier before the merge
  txt <- sub("GeneID.PathID <- merge(GeneID.PathID, PathID.PathName, by = \"PathwayID\")",
             paste0("GeneID.PathID$PathwayID <- sub(\"^path:\", \"\", GeneID.PathID$PathwayID)\n",
                    "    GeneID.PathID <- merge(GeneID.PathID, PathID.PathName, by = \"PathwayID\")"),
             txt, fixed = TRUE)
  ## parse(text=)[[1]] returns the `{...}` call itself. eval() would RUN it and
  ## assign the return value as the body, corrupting the namespace.
  body(f) <- parse(text = txt)[[1]]
  unlockBinding(".getKEGG", ns)
  assign(".getKEGG", f, envir = ns)
  lockBinding(".getKEGG", ns)
  ## Fail loudly rather than silently producing an empty KEGG result.
  chk <- f()
  stopifnot(length(chk$idList) > 100)
  cat("patched missMethyl:::.getKEGG: Hsa -> hsa, stripped 'path:' prefix |",
      length(chk$idList), "pathways\n")
})

go <- gometh(sig.cpg = sig.cpg, all.cpg = all.cpg, collection = "GO",
             array.type = "EPIC", plot.bias = FALSE)
go <- go[order(go$P.DE), ]
kegg <- gometh(sig.cpg = sig.cpg, all.cpg = all.cpg, collection = "KEGG",
               array.type = "EPIC")
kegg <- kegg[order(kegg$P.DE), ]
saveRDS(list(go = go, kegg = kegg), file.path(outd, "08_gometh.rds"))
cat("gometh GO   FDR<0.05:", sum(go$FDR   < 0.05), "| top:", rownames(go)[1],
    "-", go$TERM[1], "P =", signif(go$P.DE[1], 3), "FDR =", signif(go$FDR[1], 3), "\n")
cat("gometh KEGG FDR<0.05:", sum(kegg$FDR < 0.05), "| top:", kegg$Description[1],
    "P =", signif(kegg$P.DE[1], 3), "FDR =", signif(kegg$FDR[1], 3), "\n")

pd <- rbind(
  data.table(collection = "GO",   term = go$TERM[1:10],        P = go$P.DE[1:10],   FDR = go$FDR[1:10]),
  data.table(collection = "KEGG", term = kegg$Description[1:10], P = kegg$P.DE[1:10], FDR = kegg$FDR[1:10]))
pd[, term := ifelse(nchar(term) > 52, paste0(substr(term, 1, 49), "..."), term)]
pd[, lab := factor(paste0(term, "  (", collection, ")"),
                   levels = rev(paste0(term, "  (", collection, ")")))]
pg <- ggplot(pd, aes(-log10(P), lab, fill = collection)) +
  geom_col() +
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "gray40") +
  scale_fill_manual(values = c(GO = "#4C72B0", KEGG = "#C44E52"), name = NULL) +
  labs(x = expression(-log[10]~italic(P)~"(nominal)"), y = NULL,
       title = "gometh enrichment, top 1,000 CpGs by BACON p",
       subtitle = sprintf("No term survives FDR < 0.05 (smallest GO FDR = %.2f, KEGG FDR = %.2f)",
                          min(go$FDR), min(kegg$FDR))) +
  theme_minimal(base_size = 10) + theme(legend.position = "top")
ggsave(file.path(outd, "08_gometh_enrichment.png"), plot = pg, width = 9, height = 5.6, dpi = 200)

## ---------------- 2. methylGSA ----------------
suppressPackageStartupMessages(library(methylGSA))
pv <- setNames(tt$bacon.p, tt$probe)
mglm <- methylglm(cpg.pval = pv, array.type = "EPIC", GS.type = "KEGG",
                  minsize = 10, maxsize = 500)
mrra <- methylRRA(cpg.pval = pv, array.type = "EPIC", method = "GSEA",
                  GS.type = "KEGG", minsize = 10, maxsize = 500)
saveRDS(list(glm = mglm, rra = mrra), file.path(outd, "08_methylgsa.rds"))
n_glm <- sum(mglm$padj < 0.05); n_rra <- sum(mrra$padj < 0.05)
cat("methylglm KEGG FDR<0.05:", n_glm, "| methylRRA-GSEA KEGG FDR<0.05:", n_rra, "\n")

fwrite(as.data.table(mglm)[1:10], file.path(outd, "08_methylgsa_glm_top.csv"))
rra <- as.data.table(mrra)
setnames(rra, "ID", "KEGG", skip_absent = TRUE)
setnames(rra, "Description", "Pathway", skip_absent = TRUE)
rra_sig <- rra[padj < 0.05]
if (nrow(rra_sig) == 0) rra_sig <- head(rra, 5)
rra_sig[, `:=`(pvalue = signif(pvalue, 3), padj = signif(padj, 3))]
fwrite(rra_sig, file.path(outd, "08_methylgsa_rra_top.csv"))
cat("\n--- methylRRA top ---\n")
print(head(rra_sig[, .(KEGG, Pathway, Size, pvalue, padj)], 10))

cmp <- data.table(
  Method = c("gometh (hypergeometric, top-1,000 cutoff)",
             "methylglm (logistic regression, all p-values)",
             "methylRRA-GSEA (rank aggregation, all p-values)"),
  Input = c("top 1,000 CpGs vs universe", "full named p-value vector", "full named p-value vector"),
  `Bias model` = c("probes per gene (Wallenius)", "probe count as covariate", "gene-level RRA score"),
  `KEGG pathways FDR<0.05` = c(sum(kegg$FDR < 0.05), n_glm, n_rra))
fwrite(cmp, file.path(outd, "08_gsa_method_comparison.csv"))
cat("\n--- method comparison ---\n"); print(cmp)

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
  suppressPackageStartupMessages(library(sesameData))
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
    print(head(res[, .(dbname = get(grep("^db", names(res), value = TRUE)[1])), ], 0))
    print(head(res, 8))
    kp <- head(res[is.finite(estimate)], 12)
    kp[, lab := factor(dbname, levels = rev(dbname))]
    pk <- ggplot(kp, aes(estimate, lab, fill = FDR < 0.05)) +
      geom_col() +
      scale_fill_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "gray70"),
                        labels = c("FDR ≥ 0.05", "FDR < 0.05"), name = NULL) +
      labs(x = expression(log[2]~"odds ratio"), y = NULL,
           title = "KYCG enrichment, top 1,000 CpGs") +
      theme_minimal(base_size = 10) + theme(legend.position = "top")
    ggsave(file.path(outd, "08_kycg_enrichment.png"), plot = pk, width = 8.5, height = 5, dpi = 200)
  }
}
cat("elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
