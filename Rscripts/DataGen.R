# Purpose: Generate data for surrogate model simulations.
# Updated: 2020-12-02
library(mvnfast)

# -----------------------------------------------------------------------------

#' Sample Sizes
#'
#' Given the number of observations with complete data, and the proportions of
#' observations with target and surrogate information only, calculates the
#' number of observations in each missingness category.
#'
#' @param n0 Observations with complete data.
#' @param mt Proportion of observations with target outcome missing.
#' @param ms Proportion of observations with surrogate outcome missing.
#'
#' @return A list containing:
#' \enumerate{
#'   \item 'n0', observations with complete data.
#'   \item 'n1', observations with surrogate data only.
#'   \item 'n2', observations with target data only.
#'   \item 'n', total number of observations.
#' }

SampleSizes <- function(n0, mt = 0, ms = 0) {
  
  # Input check.
  m <- mt + ms
  if (m < 0 || m >= 1) {
    stop("mt+ms must belong to the interval [0,1).")
  }

  # Sample size.
  out <- list()
  
  # Surrogate only.
  out$n1 <- ceiling(n0 * mt / (1 - mt - ms))
  
  # Target only.
  out$n2 <- ceiling(n0 * ms / (1 - mt - ms))
  
  # Overall sample size.
  out$n <- (n0 + out$n1 + out$n2)
  
  # output
  return(out)
}


# -----------------------------------------------------------------------------

#' Generate Covariates
#' 
#' @param n Total number of subjects.
#' @return nx2 covariate matrix, with columns represeting age and sex.

GenCovar <- function(n) {
  
  # Generative design matrix
  covars <- cbind(
    "Age" = rgamma(n = n, shape = 250, rate = 5),
    "Sex" = rbinom(n = n, size = 1, prob = 0.5)
  )
  covars <- scale(x = covars)
  return(covars)
}


# -----------------------------------------------------------------------------

#' Generate Genotypes
#' 
#' @param snps Number of SNPs.
#' @param n Total number of subjects.
#' @param maf Minor allele frequency.
#' @return Subject by SNP genotype matrix.

GenGeno <- function(snps, n, maf = 0.25) {
  geno <- replicate(
    n = snps, 
    expr = rbinom(n = n, size = 2, prob = maf)
  )
  storage.mode(geno) <- "numeric"
  geno <- scale(x = geno, center = TRUE, scale = TRUE)
  return(geno)
}


#' Get Principal Components
#' 
#' @param geno Subject by SNP genotype matrix.
#' @param n_pcs Number of principal components to retain.
#' @return Subject by principal component matrix.

GetPCs <- function(geno, n_pcs = 3) {
  decomp <- svd(x = geno, nu = n_pcs, nv = 0)
  pcs <- decomp$u
  pcs <- scale(x = pcs, center = FALSE, scale = TRUE)
  return(pcs)
}


# -----------------------------------------------------------------------------

#' Calculate Effect Sizes
#' 
#' @param n_covar Number of covariates (excluding PCs).
#' @param pve_covar Proportion of variation explained by covariates.
#' @param pve_geno Proportion of variation explained by genotype.
#' @param n_pcs Number of principal components.
#' @param pve_pcs Proportion f variation explained by principal components.
#' @param resid_var Residual variance. 
#' @return 

CalcEffectSizes <- function(
  n_covar, 
  pve_covar = 0.2, 
  pve_geno, 
  n_pcs, 
  pve_pcs = 0.05, 
  resid_var = 1.00
) {
  
  # Proportion variation explained by residual.
  pve_resid <- 1 - pve_covar - pve_geno - pve_pcs
  
  if (pve_resid < 0) {
    stop("Total variation explained by covars, geno, and PCs cannot exceed 1.")
  }
  
  # Coefficients. 
  bg <- sqrt(resid_var * pve_geno / pve_resid)
  bx <- sqrt(resid_var * pve_covar / (n_covar * pve_resid))
  bs <- sqrt(resid_var * pve_pcs / (n_pcs * pve_resid))
  
  # Output.
  out <- list(
    "bg" = bg,
    "bx" = rep(bx, times = n_covar),
    "bs" = rep(bs, times = n_pcs)
  )
  return(out)
}


#' Generate Linear Predictor.
#' 
#' Does *not* include the contribution of genotype.
#' 
#' @param covars Covariate matrix.
#' @param pcs Principal components. 
#' @param pve_geno Proportion of variation explained by genotype.
#' @param pve_pcs Proportion of variation explained by principal components.
#' @return List containing:
#' \itemize{
#'   \item 'betas', the regression coefficients.
#'   \item 'eta', the nx2 linear predictor matrix.
#' }

GetLinPred <- function(covars, pcs, pve_geno, pve_pcs) {
  n_covar <- ncol(covars)
  n_pcs <- ncol(pcs)
  
  # Calculate effect sizes.
  betas <- CalcEffectSizes(
    n_covar = n_covar, 
    pve_geno = pve_geno, 
    n_pcs = n_pcs,
    pve_pcs = pve_pcs
  )
  
  # Calculate linear predictor.
  eta <- (covars %*% betas$bx) + (pcs %*% betas$bs)
  eta <- cbind(eta, eta)
  
  # Output.
  out <- list(
    "betas" = betas,
    "eta" = eta
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Generate Phenotype
#' 
#' Does *not* include the contribution of genotype.
#'
#' @param eta nx2 linear predictor matrix.
#' @param n0 Number of complete obs.
#' @param n1 Number of obs with surrogate only.
#' @param n2 Number of obs with target only.
#' @param rho Target-surrogate correlation.
#'
#' @importFrom stats rnorm
#' @importFrom mvnfast rmvn

GenPheno <- function(eta, n0, n1 = 0, n2 = 0, rho) {
  
  # Residuals.
  sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  n <- n0 + n1 + n2
  resid <- mvnfast::rmvn(n = n, mu = c(0, 0), sigma = sigma)
  
  # Response.
  y <- eta + resid
  colnames(y) <- paste0("y", seq(1:2))

  # Missingness.
  if (n1 > 0) {
    y[1:n1, 1] <- NA
  }
  if (n2 > 0) {
    y[(n1 + 1):(n1 + n2), 2] <- NA
  }
  
  # Output.
  return(y)
}


# -----------------------------------------------------------------------------

#' Data Generating Process
#' 
#' @param n0 Number of complete obs.
#' @param n1 Number of obs with surrogate only.
#' @param n2 Number of obs with target only.
#' @param pve_geno Proportion of variation explained by genotype.
#' @param pve_pcs Proportion of variation explained by principal components.
#' @param rho Target-surrogate correlation. 
#' @param snps Number of SNPs.
#' @param 

GenData <- function(n0, n1, n2, pve_geno, pve_pcs, rho, snps) {
  
  # Total observations.
  n <- n0 + n1 + n2
  
  # Covariates.
  covars <- GenCovar(n = n)
  
  # Genotypes.
  geno <- GenGeno(snps = snps, n = n)
  
  # Principal components.
  pcs <- GetPCs(geno = geno)
  
  # Linear predictor, excluding contribution of genotype.
  lin_pred <- GetLinPred(
    covars = covars, 
    pcs = pcs, 
    pve_geno = pve_geno,
    pve_pcs = pve_pcs
  )
  eta <- lin_pred$eta
  
  # Generate phenotype.
  pheno <- GenPheno(eta = eta, n0 = n0, n1 = n1, n2 = n2, rho = rho)
  
  # Output.
  out <- list(
    "betas" = lin_pred$betas,
    "covars" = covars,
    "geno" = geno,
    "pcs" = pcs,
    "pheno" = pheno
  )
  return(out)
}