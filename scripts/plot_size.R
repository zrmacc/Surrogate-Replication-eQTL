# Purpose: Construct uniform quantile-quantile plots.
# Updated: 2021-10-14

# Library.
library(cowplot)
library(dplyr)
library(ggplot2)
library(optparse)

# Base directory.
base_dir <- getwd()

# -----------------------------------------------------------------------------

# Options.
opt_list <- list()

opt <- make_option(
  c("--dir"),
  type = "character", 
  help = "Directory", 
  default = NULL
)
opt_list <- c(opt_list, opt)

opt <- make_option(
  c("--out"), 
  type = "character", 
  help = "Output filename", 
  default = "Figures/qq"
)
opt_list <- c(opt_list, opt)

# Option parsing
parsed_opts <- OptionParser(option_list = opt_list)
params <- parse_args(object = parsed_opts)

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
  data <- readRDS(file = file) %>%
    dplyr::select(Bivariate) %>%
    dplyr::pull()
  
  # Number of data points before truncation.
  n <- length(data)
  
  # Truncate to facilitate plotting.
  pval <- 0.05
  n_kept <- sum(data <= pval)
  data <- data[data <= pval]
  
  obs <- -log10(sort(data))
  exp <- -log10(seq_len(n_kept) / (n + 1))
  df <- data.frame(obs, exp)
  
  # Pointwise CIs.
  df$coord <- seq_len(n_kept)
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
# Plotting. 
# -----------------------------------------------------------------------------

# Check for output.
sobs_file <- paste0(params$out, "_sobs.png")
smis_file <- paste0(params$out, "_mis.png")

# Proceed if output missing. 
if (any(!file.exists(sobs_file), !file.exists(smis_file))) {
  
  setwd(params$dir)
  files <- dir()
  files <- files[grepl(pattern = ".*_Master\\.rds", x = files)]
  ms <- as.numeric(gsub(pattern = ".*_MS([0-9]+)_.*", replacement = "\\1", x = files))
  
  # Split into surrogate observed and surrogate missing.
  files_sobs <- files[ms == 0]
  files_smis <- files[ms >  0]
  setwd(base_dir)
  
  # Generate plot for observed surrogate.
  if (!file.exists(sobs_file) & length(files_sobs) > 0) {
    
    setwd(params$dir)
    plots_sobs <- lapply(files_sobs, PlotQQ)
    plots_sobs <- cowplot::plot_grid(
      plotlist = plots_sobs,
      ncol = 4
    )
    
    setwd(base_dir)
    ggsave(
      filename = paste0(params$out, "_sobs.png"),
      device = "png",
      plot = plots_sobs,
      width = 10,
      height = 8,
      units = "in",
      dpi = 360
    )
    
  }
  
  
  # Generate plot for missing surrogate.
  if (!file.exists(smis_file) & length(files_smis) > 0) {
    
    setwd(params$dir)
    plots_smis <- lapply(files_smis, PlotQQ)
    plots_smis <- cowplot::plot_grid(
      plotlist = plots_smis,
      ncol = 4
    )
    
    setwd(base_dir)
    ggsave(
      filename = paste0(params$out, "_smis.png"),
      device = "png",
      plot = plots_smis,
      width = 10,
      height = 6,
      units = "in",
      dpi = 360
    )
    
  }
}