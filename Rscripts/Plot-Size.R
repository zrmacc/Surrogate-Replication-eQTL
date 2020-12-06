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

# Color palettes.
rho_levels <- c(0, 25, 50, 75)
nr <- length(rho_levels)
color_ramp <- colorRampPalette(colors = c("white", "#0073C2FF"))
biv_palette <- color_ramp(n = 9)[c(3, 5, 7, 9)]

# -----------------------------------------------------------------------------

#' QQ Plot
#' 
#' @param file Filename.
#' @return ggplot

PlotQQ <- function(file) {
  
  # Import data.
  data <- readRDS(file = file)
  data <- data$Bivariate
  
  obs <- -log10(sort(data))
  n <- length(obs)
  exp <- -log10(seq_len(n) / (n + 1))
  df <- data.frame(obs, exp)
  
  # Pointwise CIs.
  df$coord <- seq_len(n)
  df$lower <- -log10(qbeta(p = (1 - 0.05 / 2), df$coord, n - df$coord + 1))
  df$upper <- -log10(qbeta(p = (0.05 / 2), df$coord, n - df$coord + 1))
  
  # Simulation configuration.
  rho <- as.numeric(gsub(pattern = ".*_R([0-9]+)_.*", replacement = "\\1", x = file))
  mt <- as.numeric(gsub(pattern = ".*_MT([0-9]+)_.*", replacement = "\\1", x = file))
  ms <- as.numeric(gsub(pattern = ".*_MS([0-9]+)_.*", replacement = "\\1", x = file))
  
  # Plot title.
  title <- bquote(m[T]:.(mt) * "%," ~ m[S]:.(ms) * "%" ~ rho:.(rho))

  # Bivariate plot
  q <- ggplot(data = df) +
    theme_bw() + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "gray",
      linetype = "dashed"
    ) + 
    geom_ribbon(
      aes(x = exp, ymin = lower, ymax = upper),
      fill = biv_palette[which(rho == rho_levels)],
      alpha = 0.5
    ) + 
    geom_point(
      aes(x = exp, y = obs),
      color = biv_palette[which(rho == rho_levels)]
    ) +
    labs(
      x = "Expected Quantile",
      y = "Observed Quantile"
    ) + 
    ggtitle(label = title)
    
  return(q)
}


# -----------------------------------------------------------------------------
# 1. Target-Missing
# -----------------------------------------------------------------------------

files <- dir()
files <- files[grepl(pattern = ".*_Master\\.rds", x = files)]
ms <- as.numeric(gsub(pattern = ".*_MS([0-9]+)_.*", replacement = "\\1", x = files))

# Split into surrogate observed and surrogate missing.
files_sobs <- files[ms == 0]
files_smis <- files[ms >  0]

# Plots for surrogate observed.
plots_sobs <- lapply(files_sobs, PlotQQ)
plots_sobs <- cowplot::plot_grid(
  plotlist = plots_sobs,
  ncol = 4
)

# Plots for surrogate missing.
plots_smis <- lapply(files_smis, PlotQQ)
plots_smis <- cowplot::plot_grid(
  plotlist = plots_smis,
  ncol = 4
)

# Output.
setwd(base_dir)

ggsave(
  filename = "Figures/qq_sobs.png",
  device = "png",
  plot = plots_sobs,
  width = 10,
  height = 8,
  units = "in",
  dpi = 360
)

ggsave(
  filename = "Figures/qq_smis.png",
  device = "png",
  plot = plots_smis,
  width = 10,
  height = 6,
  units = "in",
  dpi = 360
)

