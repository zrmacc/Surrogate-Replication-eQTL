# Purpose: Power simulations for surrogate model
# Updated: 2020-12-02

# Packages.
library(data.table)
library(optparse)
library(RNOmni)
library(SurrogateRegression)

# -----------------------------------------------------------------------------
# Unpack simulation settings.
# -----------------------------------------------------------------------------

# Command line options.
opt_list <- list()

opt <- make_option(c("--geno"), type = "character", help = "Genotype file", default = NULL)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--covar"), type = "character", help = "Covariate file", default = NULL)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--pheno"), type = "character", help = "Phenotype file", default = NULL)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--out"), type = "character", help = "Output file", default = NULL)
opt_list <- c(opt_list, opt)

# Option parsing.
t0 <- proc.time()
parsed_opts <- OptionParser(option_list = opt_list)
params <- parse_args(object = parsed_opts)

# -----------------------------------------------------------------------------
# Data import.
# -----------------------------------------------------------------------------

# Genotypes: all SNPs in *cis* to the gene. 
if (is.null(params$geno)) {
  geno <- replicate(n = 100, rbinom(n = 1000, size = 2, prob = 0.25))
} else {
  geno <- data.table::fread(file = params$geno)
}
geno <- data.matrix(geno)

# Covariates: Typically including age, sex, and genetic principal components. 
if (is.null(params$covar)) {
  covar <- replicate(n = 2, rnorm(n = 1000))
} else {
  covar <- data.table::fread(file = params$covar)
}
## Add intercept to the covariate matrix.
covar <- data.matrix(cbind(1, covar))

# Phenotypes: Expression in the target and surrogate tissues, respectively.
if (is.null(params$pheno)) {
  pheno <- replicate(n = 2, rnorm(n = 1000))
} else {
  pheno <- data.table::fread(file = params$pheno)
}
pheno <- data.matrix(pheno)

# Rank normal transform phenotypes.
pheno <- apply(pheno, 2, RNOmni::RankNorm)
t <- pheno[, 1]
s <- pheno[, 2]

# -----------------------------------------------------------------------------
# Association testing function.
# -----------------------------------------------------------------------------

#' Wrapper for Association Testing.
#' 
#' Subjects with missing genotypes are removed. The target phenotype `t`, 
#' surrogate phenotype `s`, and covariate matrix `covar`, with intercept 
#' included, are expected in the namespace. 
#' 
#' @param g Genotype at a single SNP.
#' @return Numeric vector containing the results from univariate and bivariate
#'   association testing. 

MapEQTL <- function(g) {
  
  # Remove subjects with missing genotypes at current SNP.
  is_obs <- !is.na(g)
  t.obs <- t[is_obs]
  s.obs <- s[is_obs]
  X.obs <- cbind(g[is_obs], covar[is_obs, ])
  
  # Test specification.
  test_spec <- c(TRUE, rep(FALSE, ncol(covar)))
  
  biv <- Fit.BNR(
    t = t.obs,
    s = s.obs,
    X = X.obs,
    Z = X.obs,
    eps = 1e-8
  )
  biv_out <- as.numeric(biv@Regression.tab[1, c("Point", "SE", "p")])
  
  uni <- coef(summary(lm(t.obs ~ 0 + X.obs)))
  uni_out <- as.numeric(uni[1, c("Estimate", "Std. Error", "Pr(>|t|)")])
  
  # Output
  out <- c(uni_out, biv_out)
  names(out) <- c("Uni.beta", "Uni.SE", "Uni.p", "Biv.beta", "Biv.SE", "Biv.p")
  return(out)
}


# -----------------------------------------------------------------------------
# Perform eQTL analysis.
# -----------------------------------------------------------------------------

# Analysis.
n_geno <- ncol(geno)
results <- lapply(seq_len(n_geno), function(i) {MapEQTL(geno[, i])})
results <- do.call(rbind, results)

# Output.
if (is.null(params$out)) {
  fout <- "Test/eQTL.txt"
} else {
  fout <- params$out
}

data.table::fwrite(
  x = data.frame(results),
  file = fout,
  sep = "\t"
)

# -----------------------------------------------------------------------------
# End.
# -----------------------------------------------------------------------------
t1 <- proc.time()
show(t1-t0)
