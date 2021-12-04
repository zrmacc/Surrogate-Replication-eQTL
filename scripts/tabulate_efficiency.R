# Purpose: Tabulate relative efficiency.
# Updated: 2020-12-04

# Library.
library(data.table)
library(optparse)

# Base directory.
base_dir <- getwd()

# -----------------------------------------------------------------------------

# Options.
opt_list <- list()

opt <- make_option(c("--tab"), type = "character", help = "Directory", default = NULL)
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

# -----------------------------------------------------------------------------

#' Configuration Coordinate.
#' 
#' @param data Tabulated results, containing "b.biv" and "b.uni" for each
#'   simulation setting.
#' @return Data with missingness coordinate added.
ConfigCoord <- function(data) {
  px <- unique(data[, c("rho", "mt", "ms", "pve")])
  px$coord <- seq_len(nrow(px))
  data <- merge(
    x = data,
    y = px,
    by = c("rho", "mt", "ms", "pve"),
    all.x = TRUE
  )
  return(data)
}


#' Asymptotic Relative Efficiency
#'
#' Calculates the ARE of the bivariate to the univariate
#' association test from the number of complete cases,
#' the target and surrogate missingness, and the
#' target-surrogate correlation. The variance of the
#' target and surrogates outcomes is assumed to be 1.
#'
#' @param n0 Number of complete cases
#' @param mt Target missingness
#' @param ms Surrogate missingness
#' @param r Target-surrogate correlation

ARE <- function(n0, mt, ms, rho) {
  
  # Total number of subjects.
  n <- n0 / (1 - mt - ms)
  
  # Number of subjects with target missingness.
  n1 <- mt * n
  
  # Number of subjects with surrogate missingness.
  n2 <- ms * n
  
  # Sigma.
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  sigma_inv <- solve(sigma)
  
  # Information matrices.
  ibb <- n0 * sigma_inv[1, 1] + n2 / sigma[1, 1]
  iba <- n0 * sigma_inv[1, 2]
  iaa <- n0 * sigma_inv[2, 2] + n1 / sigma[2, 2]
  
  # Efficient information.
  eff_info <- ibb - (iba)^2 / iaa
  
  # Univariate information.
  uni_info <- (n0 + n2) / sigma[1, 1]
  
  # ARE
  out <- eff_info / uni_info
  return(out)
}


#' Calculate Relative Efficiency
#' 
#' @param config Single configuration.
#' @return Data.frame with the observed and expected relative efficiency.

CalcRE <- function(config) {
  
  # Observed relative efficiency.
  obs_re <- config$EmpiricalVar[config$Parameter == "b.uni"] /
    config$EmpiricalVar[config$Parameter == "b.biv"]
  
  # Expected relative efficiency.
  n0 <- unique(config$n0) 
  rho <- unique(config$rho) / 100
  mt <- unique(config$mt) / 100
  ms <- unique(config$ms) / 100
  pve <- unique(config$pve) / 100

  # Output.
  out <- data.frame(
    "rho" = rho,
    "mt" = mt,
    "ms" = ms,
    "EmpirialRE" = obs_re,
    "ExpectedRE" = ARE(n0, mt, ms, rho)
  )
  return(out)
}

# -----------------------------------------------------------------------------
# Tabulate relative efficiency.
# -----------------------------------------------------------------------------

# Import data.
data <- data.table::fread(file = params$tab)
data$Bias <- NULL

# Calculate relative efficiency.
data <- subset(
  x = data,
  (Parameter == "b.biv") | (Parameter == "b.uni")
)
data <- ConfigCoord(data)

# Split by configuration coordinate.
split_data <- split(x = data, f = data$coord)
re_results <- lapply(split_data, CalcRE)
re_results <- do.call(rbind, re_results)

# Export data.
setwd(base_dir)
data.table::fwrite(
  x = re_results,
  file = params$out,
  sep = "\t"
)
