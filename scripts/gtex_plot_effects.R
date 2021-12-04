#' Purpose: Plot GTEx effect sizes and p-values.
#' Updated: 2021-11-25
library(cowplot)
library(dplyr)
library(data.table)
library(ggplot2)

# -----------------------------------------------------------------------------


#' Manhattan Plot
#' 
#' @param data Data.frame of results.
#' @param color Color.
#' @param min_p Minimum p-value for inclusion in plot.
#' @param y_lim Y-axis limits.
#' @return ggplot.

ManhattanPlot <- function(
  data, 
  color = "#4285F4", 
  min_p = 1e-3
) {
  data <- data[complete.cases(data), ]
  
  # Bonferroni threshold.
  threshold <- 0.05 / nrow(data)
  threshold <- -log10(threshold)
  
  # Threshold.
  data <- data %>%
    dplyr::filter(biv_p <= min_p | uni_p <= min_p)
  
  # Extract position.
  data$pos <- as.numeric(gsub(pattern = "^[0-9]+_([0-9]+)_.*", replacement = "\\1", x = data$snp))
  data <- data[order(data$chr, data$pos), ]
  data$coord <- seq_len(nrow(data))
  
  # Color alternator.
  data$alt <- factor((data$chr) %% 2)
  
  # Marker positions.
  pos <- tapply(data$coord, data$chr, median)
  pos_even <- pos[seq(from = 2, to = 22, by = 2)]
  pos_odd <- pos[seq(from = 1, to = 21, by = 2)]
  
  # Y-limits.
  y_lim <- 1.05 * max(-log10(data$biv_p), -log10(data$uni_p))
  
  # Bivariate Manhattan plot.
  q1 <- ggplot2::ggplot(data = data) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) + 
    ggplot2::geom_point(
      aes(x = coord, y = -log10(biv_p), color = alt), 
      show.legend = FALSE
    ) + 
    ggplot2::scale_color_manual(
      values = c(color, "gray")
    ) + 
    ggplot2::scale_x_continuous(
      breaks = pos_even, 
      labels = seq(from = 2, to = 22, by = 2),
      name = "Chromosome",
    ) + 
    ggplot2::scale_y_continuous(
      limits = c(0, y_lim),
      name = expression(-log[10] ~ "(" * Spray ~ p * ")")
    ) +
    ggplot2::geom_hline(
      yintercept = threshold,
      linetype = "dotted"
    ) +
    ggtitle("Manhattan Plots")
  
  # Univariate plot
  q2 <- ggplot2::ggplot(data = data) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) + 
    ggplot2::geom_point(
      aes(x = coord, y = -log10(uni_p), color = alt), 
      show.legend = FALSE
    ) + 
    ggplot2::scale_color_manual(
      values = c(color, "gray")
    ) + 
    ggplot2::scale_x_continuous(
      breaks = pos_odd, 
      labels = seq(from = 1, to = 21, by = 2),
      name = NULL,
      position = "top"
    ) + 
    ggplot2::scale_y_continuous(
      limits = c(y_lim, 0), 
      name = expression(-log[10] ~ "(" * Marginal ~ p * ")"),
      trans = "reverse"
    ) +
    ggplot2::geom_hline(
      yintercept = threshold,
      linetype = "dotted"
    )
  
  # Mirror Manhattan plots.
  q_final <- cowplot::plot_grid(q1, q2, ncol = 1, rel_heights = c(3, 2))
  return(q_final)
}


#' Create Effect Size Panel
#' 
#' @param data Data.frame.
#' @param color Point color.

PlotPanel <- function(data, color = "#4285F4") {
  data <- data[complete.cases(data), ]
  
  # Restrict to Bonferroni significant associations.
  threshold <- 0.05 / nrow(data)
  sig <- data %>%
    dplyr::filter(biv_p <= threshold | uni_p <= threshold)
  cat(nrow(sig), "significant associations.\n")
  
  # Effect sizes correlation.
  r2 <- round(cor(sig$uni_beta, sig$biv_beta)^2, digits = 3)
  r2_anno <- bquote(R^2==.(r2))
  
  # Plot bivarate vs. univariate effect sizes.
  q1 <- ggplot2::ggplot(data = sig) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(), 
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::geom_abline(
      intercept = 0, 
      slope = 1, 
      linetype = "dashed", 
      color = "gray", 
      size = 1.2
    ) + 
    ggplot2::geom_point(
      aes(x = uni_beta, y = biv_beta), 
      color = color,
      alpha = 0.8
    ) + 
    ggplot2::labs(
      x = expression(Marginal~Estimator),
      y = expression(Spray~Estimator),
      title = "Estimated Effect Size"
    ) +
    ggplot2::lims(
      x = c(-2.0, 2.0),
      y = c(-2.0, 2.0)
    ) +
    ggplot2::annotate(
      geom = "text",
      x = -1.6,
      y = 1,
      label = r2_anno,
      hjust = 0
    )
  
  # Compare p-values.
  q2 <- ggplot2::ggplot(data = sig) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(), 
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::geom_abline(
      intercept = 0, 
      slope = 1, 
      linetype = "dashed", 
      color = "gray", 
      size = 1.2
    ) + 
    ggplot2::geom_point(
      aes(x = -log10(uni_p), y = -log10(biv_p)), 
      color = color,
      alpha = 0.80
    ) + 
    ggplot2::labs(
      x = expression(-log[10] ~ "(" * Marginal ~ p * ")"),
      y = expression(-log[10] ~ "(" * Spray ~ p * ")"),
      title = "Comparison of P-values"
    ) + 
    ggplot2::lims(
      x = c(0, 25),
      y = c(0, 50)
    )
  
  q_top <- cowplot::plot_grid(q1, q2, nrow = 1, labels = c("A", "B"))
  q_bottom <- ManhattanPlot(data, color = color)
  
  # Overall plot
  q_final <- cowplot::plot_grid(
    q_top, q_bottom, 
    ncol = 1, 
    labels = c(NA, "C"),
    rel_heights = c(2, 3)
  )
  return(q_final)
}

# -----------------------------------------------------------------------------
# Compare effect sizes.
# -----------------------------------------------------------------------------

# Blood.
data <- readRDS(file = "results/ssn_blood_master.rds")
panel <- PlotPanel(data, color = "#DB4437")
ggplot2::ggsave(
  plot = panel,
  file = "figures/gtex_ssn_blood_panel.png",
  device = "png",
  width = 9.0,
  height = 9.0,
  units = "in",
  dpi = 480
)

# Cerebellum.
data <- readRDS(file = "results/ssn_cere_master.rds")
panel <- PlotPanel(data, color = "#4285F4")
ggplot2::ggsave(
  plot = panel,
  file = "figures/gtex_ssn_cere_panel.png",
  device = "png",
  width = 9.0,
  height = 9.0,
  units = "in",
  dpi = 480
)

# Muscle.
data <- readRDS(file = "results/ssn_muscle_master.rds")
panel <- PlotPanel(data, color = "#F4B400")
ggplot2::ggsave(
  plot = panel,
  file = "figures/gtex_ssn_mus_panel.png",
  device = "png",
  width = 9.0,
  height = 9.0,
  units = "in",
  dpi = 480
)

# -----------------------------------------------------------------------------
# Compare surrogates.
# -----------------------------------------------------------------------------

# Data.
blood <- readRDS(file = "results/ssn_blood_master.rds")
blood <- blood %>%
  dplyr::select(chr, gene, snp, biv_beta, biv_p) %>%
  dplyr::rename(blood_beta = biv_beta, blood_p = biv_p)
blood_threshold <- 0.05 / nrow(blood)

cere <- readRDS(file = "results/ssn_cere_master.rds")
cere <- cere %>%
  dplyr::select(chr, gene, snp, biv_beta, biv_p) %>%
  dplyr::rename(cere_beta = biv_beta, cere_p = biv_p)
cere_threshold <- 0.05 / nrow(cere)

data <- merge(x = blood, y = cere, on = c("chr", "gene", "snp"))

mus <- readRDS(file = "results/ssn_muscle_master.rds")
mus <- mus %>%
  dplyr::select(chr, gene, snp, biv_beta, biv_p) %>%
  dplyr::rename(mus_beta = biv_beta, mus_p = biv_p)
mus_threshold <- 0.05 / nrow(mus)

data <- merge(x = data, y = mus, on = c("chr", "gene", "snp"))
data <- data[complete.cases(data), ]

detections <- data %>%
  dplyr::mutate(
    sig_blood = (blood_p <= blood_threshold),
    sig_cere = (cere_p <= cere_threshold),
    sig_mus = (mus_p <= mus_threshold)
  ) %>%
  dplyr::select(chr, gene, snp, sig_blood, sig_cere, sig_mus) %>%
  dplyr::filter(sig_blood | sig_cere | sig_mus)

# Filter to eQTL significant in at least 1 analysis.
data <- data %>%
  dplyr::filter(blood_p <= blood_threshold | cere_p <= cere_threshold | mus_p <= mus_threshold)


#' Compare Effect Sizes
#' 
#' @param x Numeric.
#' @param x_lab X-axis label.
#' @param y Numeric.
#' @param y_lab Y-axis label.
#' @param color Point color.
#' @param title Plot title.
#' @return ggplot.

CompEffectSizes <- function(
  x, 
  x_lab,
  y, 
  y_lab,
  color = "#4285F4",
  plot_title = NULL
) {
  df <- data.frame(x = x, y = y)

  # Effect sizes correlation.
  r2 <- round(cor(df$x, df$y)^2, digits = 3)
  r2_anno <- bquote(R^2==.(r2))
  
  q <- ggplot2::ggplot(data = df) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) + 
    ggplot2::geom_abline(
      intercept = 0, 
      slope = 1, 
      linetype = "dashed", 
      color = "gray", 
      size = 1.2
    ) + 
    ggplot2::geom_point(
      aes(x = x, y = y), 
      color = color,
      alpha = 0.8
    ) + 
    ggplot2::labs(
      x = x_lab,
      y = y_lab,
      title = plot_title
    ) +
    ggplot2::lims(
      x = c(-2.0, 2.0),
      y = c(-2.0, 2.0)
    ) +
    ggplot2::annotate(
      geom = "text",
      x = -1.6,
      y = 1,
      label = r2_anno,
      hjust = 0
    ) 
  return(q)
}

#' Compare P-values
#' 
#' @param x Numeric.
#' @param x_lab X-axis label.
#' @param y Numeric.
#' @param y_lab Y-axis label.
#' @param color Point color.
#' @param title Plot title.
#' @return ggplot.

CompPvals <- function(
  x, 
  x_lab,
  y, 
  y_lab,
  color = "#DB4437",
  plot_title = NULL
) {
  df <- data.frame(x = x, y = y)
  df$x <- -log10(df$x)
  df$y <- -log10(df$y)
  tau <- max(c(df$x, df$y))
  
  q <- ggplot2::ggplot(data = df) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) + 
    ggplot2::geom_abline(
      intercept = 0, 
      slope = 1, 
      linetype = "dashed", 
      color = "gray", 
      size = 1.2
    ) + 
    ggplot2::geom_point(
      aes(x = x, y = y), 
      color = color,
      alpha = 0.8
    ) + 
    ggplot2::labs(
      x = x_lab,
      y = y_lab,
      title = plot_title
    ) +
    ggplot2::lims(
      x = c(0, tau),
      y = c(0, tau)
    )
  return(q)
}

# -----------------------------------------------------------------------------

# Effect sizes.
q_cere_blood <- CompEffectSizes(
  x = data$blood_beta,
  x_lab = expression(Blood~beta),
  y = data$cere_beta,
  y_lab = expression(Cerebellum~beta)
)

q_mus_blood <- CompEffectSizes(
  x = data$blood_beta,
  x_lab = expression(Blood~beta),
  y = data$mus_beta,
  y_lab = expression(Muscle~beta)
)

q_mus_cere <- CompEffectSizes(
  x = data$cere_beta,
  x_lab = expression(Cerebellum~beta),
  y = data$mus_beta,
  y_lab = expression(Muscle~beta)
)

q_top <- cowplot::plot_grid(
  q_cere_blood, q_mus_blood, q_mus_cere,
  nrow = 1
)

# -----------------------------------------------------------------------------

# Effect sizes.
q_cere_blood <- CompPvals(
  x = data$blood_p,
  x_lab = expression(-log[10]~(Blood~p)),
  y = data$cere_p,
  y_lab = expression(-log[10]~(Cerebellum~p))
)

q_mus_blood <- CompPvals(
  x = data$blood_p,
  x_lab = expression(-log[10]~(Blood~p)),
  y = data$mus_p,
  y_lab = expression(-log[10]~(Muscle~p))
)

q_mus_cere <- CompPvals(
  x = data$cere_p,
  x_lab = expression(-log[10]~(Cerebellum~p)),
  y = data$mus_p,
  y_lab = expression(-log[10]~(Muscle~p))
)

q_bottom <- cowplot::plot_grid(
  q_cere_blood, q_mus_blood, q_mus_cere,
  nrow = 1
)

# Overall figure.
q <- cowplot::plot_grid(
  q_top, q_bottom,
  ncol = 1,
  labels = c("A", "B")
)

ggplot2::ggsave(
  plot = q,
  file = "figures/gtex_surrogate_comparison.png",
  device = "png",
  width = 9.0,
  height = 6.0,
  units = "in",
  dpi = 480
)

