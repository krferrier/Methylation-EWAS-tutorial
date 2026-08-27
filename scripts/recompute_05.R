## Recompute chapter 05 (batch PCA + SVA) on the v8.1-refiltered matrix.
suppressPackageStartupMessages({
  library(minfi); library(sva); library(ggplot2); library(data.table)
})
setwd("repo"); source("_setup.R")
set.seed(42)
t0 <- Sys.time()

grs <- readRDS(data_path("03_grs_filtered.rds"))
cc  <- readRDS(data_path("04_cc_full.rds"))
cat("filtered grs:", nrow(grs), "x", ncol(grs), "\n")

beta <- getBeta(grs)
M    <- log2(beta / (1 - beta))
M[!is.finite(M)] <- NA

## metadata, aligned to the array columns (reuse the prior alignment)
prev <- readRDS(data_path("05_batch_pca.rds"))
md   <- prev$md
stopifnot(identical(md$sample_id, colnames(beta)) || nrow(md) == ncol(beta))
props <- cc$prop[md$sample_id, , drop = FALSE]
stopifnot(identical(rownames(props), md$sample_id))

## ---- PCA on the 20k most variable probes (diagnostic only) ----
v    <- apply(M, 1, var, na.rm = TRUE)
top  <- names(sort(v, decreasing = TRUE))[1:20000]
Mtop <- M[top, ]
Mtop <- Mtop[complete.cases(Mtop), ]
cat("PCA input:", nrow(Mtop), "probes\n")
pr   <- prcomp(t(Mtop), center = TRUE, scale. = FALSE)
pcs  <- pr$x[, 1:10]
pve  <- (pr$sdev^2 / sum(pr$sdev^2))[1:10]

slide  <- factor(md$slide)
posrow <- substr(md$array_pos, 1, 3)
vars <- list(chip_slide = slide, array_position = factor(posrow),
             PTSD = factor(md$ptsd), sex = factor(md$sex), age = md$age,
             Neu = props[, "Neu"], CD4T = props[, "CD4T"], CD8T = props[, "CD8T"])
r2 <- t(sapply(vars, function(x)
  sapply(1:10, function(j) summary(lm(pcs[, j] ~ x))$r.squared)))
colnames(r2) <- paste0("PC", 1:10)

saveRDS(list(pcs = pcs, pve = pve, r2 = r2, md = md, props = props,
             slide = slide, posrow = posrow),
        data_path("05_batch_pca.rds"))
cat("PC1 pve:", round(pve[1], 4), "| Neu-PC1 r2:", round(r2["Neu","PC1"], 3), "\n")

## ---- SVA on the FULL refiltered matrix, 6 cell proportions ----
keep <- complete.cases(md$ptsd, md$sex, md$age) & complete.cases(props)
mdk  <- md[keep, ]; propk <- props[keep, , drop = FALSE]
Mk   <- M[, keep]
Mk   <- Mk[complete.cases(Mk), ]
cat("SVA input:", nrow(Mk), "probes x", ncol(Mk), "samples\n")

dfk  <- data.frame(ptsd = factor(mdk$ptsd), sex = factor(mdk$sex), age = mdk$age,
                   propk, check.names = FALSE)
mod  <- model.matrix(~ ptsd + sex + age + Neu + NK + CD4T + CD8T + Bcell + Mono, data = dfk)
mod0 <- model.matrix(~        sex + age + Neu + NK + CD4T + CD8T + Bcell + Mono, data = dfk)
n.sv <- num.sv(Mk, mod, method = "be")
cat("num.sv (be):", n.sv, "\n")
svobj <- sva(Mk, mod, mod0, n.sv = n.sv)
SV <- svobj$sv
colnames(SV) <- paste0("SV", seq_len(ncol(SV)))

saveRDS(list(SV = SV, n.sv = n.sv, keep = keep, mdk = mdk, propk = propk),
        data_path("05_sva.rds"))

chip_r2 <- sapply(1:ncol(SV), function(j) summary(lm(SV[,j] ~ factor(mdk$slide)))$r.squared)
ptsd_r2 <- sapply(1:ncol(SV), function(j) summary(lm(SV[,j] ~ factor(mdk$ptsd)))$r.squared)
cat("strongest chip SV: SV", which.max(chip_r2), " r2=", round(max(chip_r2),3), "\n", sep="")
cat("max SV-PTSD r2:", round(max(ptsd_r2), 4), "\n")

## ---- figures ----
hm <- as.data.frame(as.table(r2)); names(hm) <- c("variable","PC","r2")
p1 <- ggplot(hm, aes(PC, variable, fill = r2)) + geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r2)), size = 2.7,
            color = ifelse(hm$r2 > 0.5, "white", "gray20")) +
  scale_fill_gradient(low = "#F7F7F7", high = "#C44E52", limits = c(0, 1), name = expression(R^2)) +
  labs(x = NULL, y = NULL, title = "PC associations with technical and biological variables") +
  theme_minimal(base_size = 11) + theme(panel.grid = element_blank())
ggsave(data_path("05_pc_heatmap.png"), plot = p1, width = 8, height = 4, dpi = 200)

sc <- data.frame(PC1 = pcs[,1], PC2 = pcs[,2], Neu = props[,"Neu"])
p2 <- ggplot(sc, aes(PC1, PC2, color = Neu)) + geom_point(size = 2.6, alpha = 0.9) +
  scale_color_viridis_c(name = "Neutrophil\nfraction") +
  labs(title = "PC1-PC2, colored by estimated neutrophil fraction",
       x = sprintf("PC1 (%.1f%%)", 100*pve[1]), y = sprintf("PC2 (%.1f%%)", 100*pve[2])) +
  theme_minimal(base_size = 11)
ggsave(data_path("05_pc_scatter_neu.png"), plot = p2, width = 7, height = 5, dpi = 200)

bars <- rbind(
  data.frame(SV = factor(paste0("SV", 1:n.sv), levels = paste0("SV", 1:n.sv)),
             r2 = chip_r2, variable = "chip (technical)"),
  data.frame(SV = factor(paste0("SV", 1:n.sv), levels = paste0("SV", 1:n.sv)),
             r2 = ptsd_r2, variable = "PTSD (exposure)"))
p3 <- ggplot(bars, aes(SV, r2, fill = variable)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("chip (technical)" = "#C44E52", "PTSD (exposure)" = "#4C72B0"), name = NULL) +
  labs(title = "What each surrogate variable captures", x = NULL, y = expression(R^2)) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(data_path("05_sva_bars.png"), plot = p3, width = 8, height = 4, dpi = 200)

cat("elapsed:", round(difftime(Sys.time(), t0, units="mins"),2), "min\n")
