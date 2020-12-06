# Purpose: Plot GTEx results.
# Updated: 2020-12-05
library(cowplot)
library(data.table)
library(ggplot2)

# -----------------------------------------------------------------------------
# Compare effect sizes.
# -----------------------------------------------------------------------------

# Import data. 
data <- readRDS(file = "Data/SSN.rds")
data <- data[complete.cases(data), ]

# Restrict to genome-wide significant.
sig <- subset(
  x = data,
  (Uni.p <= 5e-8) | (Biv.p <= 5e-8)
)

# Effect sizes correlation.
r2 <- round(cor(sig$Uni.beta, sig$Biv.beta)^2, digits = 3)
r2_anno <- bquote(R^2==.(r2))

q1 <- ggplot(data = sig) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  geom_abline(
    intercept = 0, 
    slope = 1, 
    linetype = "dashed", 
    color = "#0073C2FF", 
    size = 1.2
  ) + 
  geom_point(
    aes(x = Uni.beta, y = Biv.beta), 
    color = "gray"
  ) + 
  labs(
    x = expression(Marginal ~ hat(beta)),
    y = expression(Spray ~ hat(beta)),
    title = "Estimated Effect Size"
  ) +
  lims(
    x = c(-2.25, 2),
    y = c(-2.25, 2)
  ) +
  annotate(
    geom = "text",
    x = -1.6,
    y = 1,
    label = r2_anno
  )

# -----------------------------------------------------------------------------
# Compare effect sizes.
# -----------------------------------------------------------------------------

q2 <- ggplot(data = sig) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  geom_abline(
    intercept = 0, 
    slope = 1, 
    linetype = "dashed", 
    color = "#0073C2FF", 
    size = 1.2
  ) + geom_point(
    aes(x = -log10(Uni.p), y = -log10(Biv.p)), 
    color = "gray"
  ) + labs(
    x = expression(-log[10] ~ "(" * Marginal ~ p * ")"),
    y = expression(-log[10] ~ "(" * Spray ~ p * ")"),
    title = "Comparison of P-values"
  ) + lims(
    x = c(0, 25),
    y = c(0, 50)
  )

q_top <- cowplot::plot_grid(q1, q2, nrow = 1, labels = c("A", "B"))

# -----------------------------------------------------------------------------
# Empirical relative efficiency by r2.
# -----------------------------------------------------------------------------

# Correlation between expression in brain and expression in blood.
rho <- readRDS(file = "Data/Corr.rds")
rho$r2 <- (rho$Pearson)^2
rho <- subset(
  x = rho,
  select = c("Ensembl", "r2")
)

# Add target-surrogate correlation.
sig <- merge(
  x = sig,
  y = rho,
  by = "Ensembl"
)

# Add empirical relative efficiency.
sig$re <- (sig$Uni.SE / sig$Biv.SE)^2

# Plot.
q_bottom <- ggplot(data = sig) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  geom_point(
    aes(x = r2, y = re), 
    alpha = 0.8, 
    color = "gray"
  ) + 
  stat_smooth(
    aes(x = r2, y = re),
    method = "loess", 
    size = 1.2,
    linetype = "dashed", 
    se = FALSE, 
    color = "#0073C2FF"
  ) + 
  labs(
    x = expression(r^2),
    y = "Empirical RE",
    title = expression(Empirical ~ Relative ~ Efficiency ~ by ~ Target - Surrogate ~ r^2)
  ) + 
  scale_y_continuous(
    limits = c(0.95, 2.0),
    breaks = c(1.00, 1.25, 1.50, 1.75, 2.00)
  ) + geom_hline(
    yintercept = 1.0,
    linetype = "dotted",
    color = "#0073C2FF"
  )

# Overdata plot
q_final <- cowplot::plot_grid(q_top, q_bottom, ncol = 1)

ggsave(
  plot = q_final, 
  filename = "Figures/gtex.png", 
  device = "png", 
  scale = 1, 
  width = 6.5, 
  height = 6.0, 
  units = "in", 
  dpi = 360
)

# -----------------------------------------------------------------------------
# Manhattan. 
# -----------------------------------------------------------------------------

# Extract chromosome.
data$chr <- as.numeric(gsub(pattern = "^([0-9]+)_.*", replacement = "\\1", x = data$SNP))

# Extract position.
data$pos <- as.numeric(gsub(pattern = "^[0-9]+_([0-9]+)_.*", replacement = "\\1", x = data$SNP))
data <- data[order(data$chr, data$pos), ]
data$coord <- seq_len(nrow(data))
data$alt <- factor((data$chr) %% 2)

# Marker positions.
pos <- tapply(data$coord, data$chr, median)
pos.even <- pos[seq(from = 2, to = 22, by = 2)]

# Bivariate Manhattan plot.
q1 <- ggplot(data = data) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) + 
  geom_point(
    aes(x = coord, y = -log10(Biv.p), color = alt), 
    show.legend = FALSE
  ) + 
  scale_color_manual(
    values = c("#0073C2FF", "gray")
  ) + scale_x_continuous(
    breaks = pos.even, 
    labels = seq(from = 2, to = 22, by = 2)
  ) + 
  labs(
    x = "Genomic Position",
    y = -log[10] ~ "(" * Spray ~ p * ")"
  ) + 
  ylim(c(0, 60))


# Univariate plot
q2 <- ggplot(data = data) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) + 
  geom_point(
    aes(x = coord, y = -log10(Uni.p), color = alt), 
    show.legend = FALSE
  ) + scale_color_manual(
    values = c("#0073C2FF", "gray")
  ) + 
  scale_x_continuous(
    breaks = pos.even, 
    labels = seq(from = 2, to = 22, by = 2)
  ) + 
  labs(
    x = NULL,
    y = -log[10] ~ "(" * Marginal ~ p * ")"
  ) + 
  scale_y_continuous(
    limits = c(60, 0), 
    trans = "reverse"
  )

# Mirror manhattan plots
q_final <- plot_grid(q1, q2, ncol = 1)
ggsave(
  plot = q_final, 
  filename = "Figures/Manhattan.png", 
  device = "png", 
  scale = 2,
  width = 6, 
  height = 3, 
  units = "in", 
  dpi = 360
)



