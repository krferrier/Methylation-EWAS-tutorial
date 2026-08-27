## How much residual variance does each surrogate variable actually carry?
## SVs are estimated from the matrix after the full model is projected out, so the
## honest denominator is the residual variance, not total variance.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setwd("repo"); source("_setup.R")
t0 <- Sys.time()

g   <- readRDS(data_path("03_grs_filtered.rds"))
sva <- readRDS(data_path("05_sva.rds"))

beta <- minfi::getBeta(g)
M <- log2(beta / (1 - beta)); M[!is.finite(M)] <- NA
M <- M[, sva$keep, drop = FALSE]
ok <- rowSums(is.na(M)) == 0
M <- M[ok, , drop = FALSE]
cat("M:", paste(dim(M), collapse = " x "), "\n")

## rebuild the same full model used for the EWAS (sva$mdk / sva$propk are aligned)
mdk <- sva$mdk; propk <- sva$propk
ptsd <- factor(mdk$ptsd, levels = c("Control", "Case"))
if (all(is.na(ptsd))) ptsd <- factor(mdk$ptsd)
mod <- model.matrix(~ ptsd + sex + age + Neu + NK + CD4T + CD8T + Bcell + Mono,
                    data = data.frame(ptsd = ptsd, sex = factor(mdk$sex), age = mdk$age,
                                      propk, check.names = FALSE))
cat("full model columns:", ncol(mod), "\n")

## residualize: R = M (I - H), H = hat matrix of the full model
H <- mod %*% solve(crossprod(mod), t(mod))
R <- M - (M %*% H)
tot <- sum(R^2)

SV <- sva$SV
pve <- apply(SV, 2, function(v) { v <- v / sqrt(sum(v^2)); sum((R %*% v)^2) / tot })
names(pve) <- paste0("SV", seq_along(pve))

## also the unconstrained scree of the residual space for context
ev <- eigen(crossprod(R) / tot, symmetric = TRUE, only.values = TRUE)$values
ev <- ev[ev > 0]

d <- data.table(k = seq_along(pve), sv = names(pve), pve = 100 * pve)
d[, cum := cumsum(pve)]
d[, share := 100 * pve / sum(pve)]
d[, cum_share := cumsum(share)]
print(d[, .(sv, pve = round(pve, 2), cum_resid = round(cum, 2),
            share_of_SV_block = round(share, 1), cum_share = round(cum_share, 1))])

cat("\nSV block captures", round(sum(d$pve), 1), "% of residual variance\n")
k80 <- d[cum_share >= 80][1]$k
cat("k for 80% of the SV-captured variance:", k80, "\n")

## elbow: largest drop in successive pve, and the point of diminishing returns
drops <- -diff(d$pve)
cat("successive drops (pp):", paste(round(drops, 2), collapse = " "), "\n")
cat("largest drop after SV:", which.max(drops), "\n")

## scree plot with both criteria marked
p <- ggplot(d, aes(k, pve)) +
  geom_line(color = "gray55") +
  geom_point(aes(fill = k <= k80), shape = 21, size = 3, stroke = 0.4) +
  geom_vline(xintercept = k80 + 0.5, linetype = 2, color = "#C44E52") +
  annotate("text", x = k80 + 0.7, y = max(d$pve) * 0.92, hjust = 0,
           label = sprintf("k = %d\n(%.0f%% of SV-captured\nresidual variance)", k80, d[k == k80]$cum_share),
           size = 3.1, color = "#C44E52") +
  scale_fill_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "gray75"), guide = "none") +
  scale_x_continuous(breaks = d$k) +
  labs(x = "surrogate variable", y = "% of residual variance",
       title = "Scree plot of the surrogate variables",
       subtitle = sprintf("15 SVs from num.sv(be); together they carry %.0f%% of the residual variance",
                          sum(d$pve))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(data_path("05_sv_scree.png"), plot = p, width = 7.5, height = 4.4, dpi = 200)

fwrite(d, data_path("05_sv_variance.csv"))
saveRDS(list(pve = pve, k80 = k80, resid_eigen = head(ev, 30)), data_path("05_sv_variance.rds"))

for (k in c(k80, 15)) {
  cat(sprintf("\nk=%2d -> Overall resid.df=%2d | Female=%2d | Male=%2d",
              k, 87 - (10 + k), 45 - (9 + k), 42 - (9 + k)))
}
cat("\nelapsed:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
