#' Purpose: Estimation and type I error simulations.
#' Updated: 2021-12-04

# Packages.
library(optparse)
base_dir <- getwd()
setwd("scripts/")
source("sim_loops.R")

# -----------------------------------------------------------------------------
# Unpack simulation settings.
# -----------------------------------------------------------------------------

# Command line options.
opt_list <- list()

opt <- make_option(c("--array_idx"), type = "integer", help = "Job array index", default = 1)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--config"), type = "character", help = "Config file", default = 1)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--dist"), type = "character", help = "Distribution", default = "norm")
opt_list <- c(opt_list, opt)

opt <- make_option(c("--n0"), type = "integer", help = "Observed subjects", default = 1000)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--out"), type = "character", help = "Output stem", default = "Test/")
opt_list <- c(opt_list, opt)

opt <- make_option(c("--reps"), type = "integer", help = "MC replicates", default = 500)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--snps"), type = "integer", help = "SNPs", default = 200)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--testing"), type = "logical", help = "Testing?", default = TRUE)
opt_list <- c(opt_list, opt)

# Option parsing.
t0 <- proc.time()
parsed_opts <- OptionParser(option_list = opt_list)
params <- parse_args(object = parsed_opts)

# Load config file.
size_config <- data.table::fread(file = file.path(base_dir, params$config))
array_idx <- params$array_idx
params$idx <- size_config$idx[array_idx]
params$mT <- size_config$mT[array_idx]
params$mS <- size_config$mS[array_idx]
params$rho <- size_config$rho[array_idx]

if (params$testing) {
  params$n0 <- 100
  params$reps <- 10
  params$snps <- 10
}

# Output stem.
out_suffix <- paste0(
  "N", params$n0,
  "_", params$dist,
  "_MT", params$mT * 100, 
  "_MS", params$mS * 100, 
  "_R", params$rho * 100, 
  "_", params$idx,
  ".rds"
)

# -----------------------------------------------------------------------------
# Simulation.
# -----------------------------------------------------------------------------

sim <- lapply(seq_len(params$reps), function(b) {
  out <- SizeSimLoop(params = params)
})
sim <- data.frame(do.call(rbind, sim))

# -----------------------------------------------------------------------------
# Save summarized point estimates.
# -----------------------------------------------------------------------------

# Point estimates.
points <- c("point.b.biv", "point.b.uni", "point.a", "point.S11", "point.S12", "point.S22")
means <- apply(sim[, points], 2, mean)

# Model variance.
ses <- c("se.b.biv", "se.b.uni", "se.a", "se.S11", "se.S12", "se.S22")
var_model <- apply(sim[, ses], 2, function(x) {return(mean(x^2))})

# Empirical variance.
var_emp <- apply(sim[, points], 2, var)

# Output.
est <- data.frame(
  Parameter = c("b.biv", "b.uni", "a", "S11", "S12", "S22"),
  N = nrow(sim),
  Bias = means,
  ModelVar =  var_model,
  EmpiricalVar = var_emp
)

setwd(base_dir)

out_stem <- PathJoin(params$out, "est_h0")
if (!dir.exists(out_stem)) {dir.create(out_stem, recursive = TRUE)}
out_file <- PathJoin(out_stem, out_suffix)
saveRDS(object = est, file = out_file)

# -----------------------------------------------------------------------------
# Save p-values for type I error calculation.
# -----------------------------------------------------------------------------

pvals <- subset(
  x = sim,
  select = c("p.b.biv", "p.b.uni")
)
colnames(pvals) <- c("Bivariate", "Univariate")

out_stem <- PathJoin(params$out, "t1e")
if (!dir.exists(out_stem)) {dir.create(out_stem, recursive = TRUE)}
out_file <- PathJoin(out_stem, out_suffix)
saveRDS(object = pvals, file = out_file)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------
t1 <- proc.time()
show(t1-t0)
