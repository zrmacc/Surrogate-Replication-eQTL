# Purpose: Plot results from power simulations
# Updated: 2020-12-04

# Library.
library(data.table)
library(optparse)

# Base directory.
base_dir <- getwd()

# -----------------------------------------------------------------------------

# Options.
opt_list <- list()

opt <- make_option(c("--dir"), type = "character", help = "Directory", default = NULL)
opt_list <- c(opt_list, opt)

opt <- make_option(
  c("--out"), 
  type = "character", 
  help = "Output filename", 
  default = "Tables/output"
)
opt_list <- c(opt_list, opt)

# Option parsing
parsed_opts <- OptionParser(option_list = opt_list)
params <- parse_args(object = parsed_opts)
setwd(params$dir)

# -----------------------------------------------------------------------------

# Type I error levels.
k_alpha <- c(0.05, 1e-3, 1e-5)

# -----------------------------------------------------------------------------

#' Data Import Loop.
#' 
#' @param file Current file.
#' @return Data.frame with all power results.

Loop <- function(file) {
  
  data <- readRDS(file = file)
  pwr <- sapply(k_alpha, function(x){ mean(data$Bivariate <= x)})
  rho <- as.numeric(gsub(pattern = ".*_R([0-9]+)_.*", replacement = "\\1", x = file))
  mt <- as.numeric(gsub(pattern = ".*_MT([0-9]+)_.*", replacement = "\\1", x = file))
  ms <- as.numeric(gsub(pattern = ".*_MS([0-9]+)_.*", replacement = "\\1", x = file))
  if (grepl(pattern = "PVE", x = file)) {
    pve <- as.numeric(gsub(pattern = ".*_PVE([0-9]+.*)_.*", replacement = "\\1", x = file))
  } else {
    pve <- 0
  }
  out <- data.frame(
    "rho" = rho,
    "mt" = mt,
    "ms" = ms,
    "pve" = pve,
    "alpha" = k_alpha,
    "power" = pwr
  )
  return(out)
}

# -----------------------------------------------------------------------------
# Tabulate rejection probability.
# -----------------------------------------------------------------------------

files <- dir()
files <- files[grepl(pattern = ".*_Master\\.rds", x = files)]

# Import data.
data <- lapply(files, Loop)
data <- do.call(rbind, data)

# Export data.
setwd(base_dir)
data.table::fwrite(
  x = data,
  file = params$out,
  sep = "\t"
)
