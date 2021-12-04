#' Purpose: Plot power curves.
#' Updated: 2021-10-17

# Library.
library(cowplot)
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
  default = "Figures/power"
)
opt_list <- c(opt_list, opt)

# Option parsing
parsed_opts <- OptionParser(option_list = opt_list)
params <- parse_args(object = parsed_opts)

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
  
  # Import data.
  data <- lapply(files, Loop)
  data <- do.call(rbind, data)
  
  # Split by surrogate missingness.
  data_sobs <- data[data$ms == 0, ]
  data_smis <- data[data$ms  > 0, ]
  
  setwd(base_dir)
  
  # Generate plot for observed surrogate.
  if (!file.exists(sobs_file) & nrow(data_sobs) > 0) {
    
    setwd(params$dir)
    data_sobs <- MissCoord(data_sobs)
    
    data_sobs_split <- split(data_sobs, data_sobs$coord)
    data_sobs_plots <- lapply(data_sobs_split, PlotPower)
    sobs_plots <- cowplot::plot_grid(
      plotlist = data_sobs_plots,
      ncol = 2
    )
    
    setwd(base_dir)
    ggsave(
      plot = sobs_plots,
      filename = paste0(params$out, "_sobs.png"),
      device = "png",
      height = 6,
      width = 10,
      units = "in",
      dpi = 360
    )
    
  }
  
  # Generate plot for missing surrogate.
  if (!file.exists(smis_file) & nrow(data_smis) > 0) {
    
    setwd(params$dir)
    data_smis <- MissCoord(data_smis)
    
    data_smis_split <- split(data_smis, data_smis$coord)
    data_smis_plots <- lapply(data_smis_split, PlotPower)
    smis_plots <- cowplot::plot_grid(
      plotlist = data_smis_plots,
      ncol = 1
    )
    
    setwd(base_dir)
    ggsave(
      plot = smis_plots,
      filename = paste0(params$out, "_smis.png"),
      device = "png",
      height = 6,
      width = 7.5,
      units = "in",
      dpi = 360
    )
    
  }
  
}
