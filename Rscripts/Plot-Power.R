# Purpose: Plot results from power simulations
# Updated: 2020-12-04

# Library.
library(cowplot)
library(ggplot2)
library(optparse)

# Base directory.
base_dir <- getwd()

# -----------------------------------------------------------------------------

# Options.
opt_list <- list()

opt <- make_option(c("--dir"), type = "character", help = "directory", default = NULL)
opt_list <- c(opt_list, opt)

# Option parsing
parsed_opts <- OptionParser(option_list = opt_list)
params <- parse_args(object = parsed_opts)
setwd(params$dir)

# -----------------------------------------------------------------------------

# Type I error level.
k_alpha <- 0.05

# Color palettes.
rho <- c(0, 25, 50, 75)
nr <- length(rho)
color_ramp <- colorRampPalette(colors = c("white", "#0073C2FF"))
biv_palette <- color_ramp(n = 9)[c(3, 5, 7, 9)]

# -----------------------------------------------------------------------------

#' Data Import Loop.
#' 
#' @param file Current file.
#' @return Data.frame with all power results.

Loop <- function(file) {
  
  data <- readRDS(file = file)
  pwr <- mean(data$Bivariate <= k_alpha)
  rho <- as.numeric(gsub(pattern = ".*_R([0-9]+)_.*", replacement = "\\1", x = file))
  mt <- as.numeric(gsub(pattern = ".*_MT([0-9]+)_.*", replacement = "\\1", x = file))
  ms <- as.numeric(gsub(pattern = ".*_MS([0-9]+)_.*", replacement = "\\1", x = file))
  pve <- as.numeric(gsub(pattern = ".*_PVE([0-9]+.*)_.*", replacement = "\\1", x = file))
  out <- data.frame(
    "rho" = rho,
    "mt" = mt,
    "ms" = ms,
    "pve" = pve,
    "power" = pwr
  )
  return(out)
  
}

# -----------------------------------------------------------------------------

#' Missingness Coordinate.
#' 
#' @param data As returned by [Loop].
#' @return Data with missingness coordinate added.
MissCoord <- function(data) {
  px <- unique(data[, c("mt", "ms")])
  px$coord <- seq_len(nrow(px))
  data <- merge(
    x = data,
    y = px,
    by = c("mt", "ms"),
    all.x = TRUE
  )
  return(data)
}

#' Power Plot
#' 
#' @param data For a single missingness pattern.
#' @return ggplot

PlotPower <- function(data) {
  
  # Create factors
  mt <- as.numeric(unique(data$mt))
  ms <- as.numeric(unique(data$ms))
  data$rho <- factor(data$rho, levels = sort(unique(data$rho)), ordered = TRUE)
  data$mt <- factor(data$mt, levels = sort(unique(data$mt)), ordered = TRUE)
  data$ms <- factor(data$ms, levels = sort(unique(data$ms)), ordered = TRUE)
  
  # Plot title.
  title <- bquote(pi[T]:.(mt) * "%," ~ pi[S]:.(ms) * "%")

  # Bivariate plot
  q <- ggplot(data = data) +
    theme_bw() + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    ) +
    geom_point(aes(x = pve, y = power, color = rho)) +
    geom_line(aes(x = pve, y = power, color = rho)) +
    scale_color_manual(name = expression(rho ~ "(%)"), values = biv_palette) +
    scale_x_continuous(
      name = "Per-SNP Heritability (%)", 
      breaks = sort(unique(data$pve))
    ) +
    scale_y_continuous(
      name = "Power", 
      breaks = seq(from = 0, to = 1.0, by = 0.2), 
      limits = c(0, 1.0)
    ) +
    ggtitle(label = title)
    
  return(q)
}

# -----------------------------------------------------------------------------
# 1. Target-Missing
# -----------------------------------------------------------------------------

files <- dir()
files <- files[grepl(pattern = ".*_Master\\.rds", x = files)]

# Import data.
data <- lapply(files, Loop)
data <- do.call(rbind, data)

# Split by surrogate missingness.
data_sobs <- data[data$ms == 0, ]
data_smis <- data[data$ms  > 0, ]

data_sobs <- MissCoord(data_sobs)
data_smis <- MissCoord(data_smis)

# Plot surrogate observed.
data_sobs_split <- split(data_sobs, data_sobs$coord)
data_sobs_plots <- lapply(data_sobs_split, PlotPower)
sobs_plots <- cowplot::plot_grid(
  plotlist = data_sobs_plots,
  ncol = 2
)

# Plot surrogate missing.
data_smis_split <- split(data_smis, data_smis$coord)
data_smis_plots <- lapply(data_smis_split, PlotPower)
smis_plots <- cowplot::plot_grid(
  plotlist = data_smis_plots,
  ncol = 1
)

# Output.
setwd(base_dir)
ggsave(
  plot = smis_plots,
  filename = "Figures/power_smis.png",
  device = "png",
  width = 5.5,
  height = 7.5,
  units = "in",
  dpi = 360
)

ggsave(
  plot = sobs_plots,
  filename = "Figures/power_sobs.png",
  device = "png",
  width = 7.5,
  height = 5.5,
  units = "in",
  dpi = 360
)

