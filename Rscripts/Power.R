# Purpose: Power simulations for surrogate model
# Updated: 2020-12-02

# Packages.
library(optparse)
library(SurrogateRegression)
base_dir <- getwd()
setwd("Rscripts/")
source("SimLoops.R")

# -----------------------------------------------------------------------------
# Unpack simulation settings.
# -----------------------------------------------------------------------------

# Command line options.
opt_list <- list()

# Sample size.
opt <- make_option(c("--n0"), type = "integer", help = "Observed subjects", default = 1000)
opt_list <- c(opt_list, opt)

# Simulation settings.
opt <- make_option(c("--snps"), type = "integer", help = "SNPs", default = 100)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--reps"), type = "integer", help = "MC replicates", default = 2)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--idx"), type = "integer", help = "Job index", default = 1)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--out"), type = "character", help = "Output stem", default = "Test/")
opt_list <- c(opt_list, opt)

# Data generating process.
opt <- make_option(c("--mT"), type = "numeric", help = "Target Missingness", default = 0.00)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--mS"), type = "numeric", help = "Surrogate Missingness", default = 0.00)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--rho"), type = "numeric", help = "Correlation", default = 0.50)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--pve"), type = "numeric", help = "Heritability", default = 0.05)
opt_list <- c(opt_list, opt)

# Option parsing.
t0 <- proc.time()
parsed_opts <- OptionParser(option_list = opt_list)
params <- parse_args(object = parsed_opts)

# Output stem.
out_suffix <- paste0(
  "N", params$n0, 
  "_MT", params$mT * 100, 
  "_MS", params$mS * 100, 
  "_R", params$rho * 100, 
  "_PVE", params$pve * 100,
  "_", params$idx,
  ".rds"
)

# -----------------------------------------------------------------------------
# Simulations
# -----------------------------------------------------------------------------

sim <- lapply(seq_len(params$reps), function(b) {
  out <- PowerSimLoop(params = params)
  return(out)
})
sim <- data.frame(do.call(rbind, sim))

# -----------------------------------------------------------------------------
# Save summarized point estimates.
# -----------------------------------------------------------------------------

# Point estimates.
points <- c("point.b.biv", "point.b.uni")
means <- apply(sim[, points], 2, mean)

# Model variance.
ses <- c("se.b.biv", "se.b.uni")
var_model <- apply(sim[, ses], 2, function(x) {return(mean(x^2))})

# Empirical variance.
var_emp <- apply(sim[, points], 2, var)

# Output.
est <- data.frame(
  "Parameter" = c("b.biv", "b.uni"),
  "N" = nrow(sim),
  "Bias" = means,
  "ModelVar" =  var_model,
  "EmpiricalVar" = var_emp
)

setwd(base_dir)

out_stem <- paste0(params$out, "EstH1")
if (!dir.exists(out_stem)) {dir.create(out_stem, recursive = TRUE)}
out_file <- paste0(out_stem, "/", out_suffix)
saveRDS(object = est, file = out_file)


# -----------------------------------------------------------------------------
# Save p-values for power calculation.
# -----------------------------------------------------------------------------

pvals <- subset(
  x = sim,
  select = c("p.b.biv", "p.b.uni")
)
colnames(pvals) <- c("Bivariate", "Univariate")
out_stem <- paste0(params$out, "Power")
if (!dir.exists(out_stem)) {dir.create(out_stem, recursive = TRUE)}
out_file <- paste0(out_stem, "/", out_suffix)
saveRDS(object = pvals, file = out_file)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------
t1 <- proc.time()
show(t1-t0)
