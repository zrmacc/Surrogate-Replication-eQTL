# Purpose: Concatenate results for tabulation.
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

#' Data Import Loop.
#' 
#' @param file Current file.
#' @return Data.frame with estimation results.

Loop <- function(file) {
  
  data <- readRDS(file = file)
  n0 <- as.numeric(gsub(pattern = "^N([0-9]+)_.*", replacement = "\\1", x = file))
  rho <- as.numeric(gsub(pattern = ".*_R([0-9]+)_.*", replacement = "\\1", x = file))
  mt <- as.numeric(gsub(pattern = ".*_MT([0-9]+)_.*", replacement = "\\1", x = file))
  ms <- as.numeric(gsub(pattern = ".*_MS([0-9]+)_.*", replacement = "\\1", x = file))
  if (grepl(pattern = "PVE", x = file)) {
    pve <- as.numeric(gsub(pattern = ".*_PVE([0-9]+.*)_.*", replacement = "\\1", x = file))
  } else {
    pve <- 0
  }
  out <- cbind(
    "n0" = n0,
    "rho" = rho,
    "mt" = mt,
    "ms" = ms,
    "pve" = pve,
    data
  )
  return(out)
}

# -----------------------------------------------------------------------------
# Tabulate mean estimates.
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
