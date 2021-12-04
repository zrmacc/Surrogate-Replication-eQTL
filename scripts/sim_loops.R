#' Purpose: Simulation loops.
#' Updated: 2021-11-30
source("data_gen.R")

# -----------------------------------------------------------------------------

#' Join File Paths
#' 
#' @param a String.
#' @param b String.
#' @return File path.
PathJoin <- function(a, b) {
  a_last <- rev(unlist(strsplit(a, split = "")))[1]
  if (a_last == "/") {
    out <- paste0(a, b)
  } else {
    out <- paste0(a, "/", b)
  }
  return(out)
}

# -----------------------------------------------------------------------------

#' Rank Normal Transformation
#' 
#' @param u Numeric vector.
#' @return Numeric vector.

RankNorm <- function(u) {
  probs <- (order(u) - 0.5) / (length(u))
  z <- stats::qnorm(probs)
  return(z)
}


# -----------------------------------------------------------------------------
# Size Simulations.
# -----------------------------------------------------------------------------

#' Bivariate Test for Size Simulations.
#' 
#' @param g Genotype vector.
#' @param t Target outcome.
#' @param s Surrogate outcome.
#' @param design Design matrix.
#' @return Vector of point estimates and standard errors for parameters of 
#'   interest from the fitted bivariate regression model.

BivTestSize <- function(g, t, s, design) {
  
  geno_design <- cbind(g, design)
  
  fit <- SurrogateRegression::Fit.BNR(
    t = t, 
    s = s, 
    X = geno_design, 
    report = FALSE
  )
  
  # Extract parameters of interest.
  reg_param <- fit@Regression.tab
  reg_param <- subset(
    x = reg_param,
    Coefficient == "g",
    select = c("Point", "SE", "p")
  )
  cov_param <- fit@Covariance.tab
  
  # Output
  out <- c(
    point.b.biv = reg_param$Point[1],
    se.b.biv = reg_param$SE[1],
    p.b.biv = reg_param$p[1],
    point.a = reg_param$Point[2],
    se.a = reg_param$SE[2],
    p.a = reg_param$p[2],
    point.S11 = cov_param$Point[1],
    se.S11 = cov_param$SE[1],
    point.S12 = cov_param$Point[2],
    se.S12 = cov_param$SE[2],
    point.S22 = cov_param$Point[3],
    se.S22 = cov_param$SE[3]
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Univariate Test for Size Simulations.
#' 
#' @param g Genotype vector.
#' @param t Target outcome.
#' @param s Surrogate outcome.
#' @param design Design matrix.
#' @return Vector of point estimates and standard errors for parameters of 
#'   interest from the fitted bivariate regression model.

UniTestSize <- function(g, t, s, design) {
  
  # Subset to subjects with observed target outcome.
  keep <- !is.na(t)
  geno_design <- cbind(g, design)
  
  # Ordinary least squares inference.
  fit <- MatrixOps::fitOLS(y = t[keep], X = geno_design[keep, ])
  b <- fit$Beta[1]
  se <- sqrt(diag(MatrixOps::matInv(fit$Ibb)))[1]
  df <- sum(keep) - ncol(geno_design)
  p <- 2 * stats::pt(q = abs(b / se), df = df, lower.tail = FALSE)
  
  # Output
  out <- c(
    point.b.uni = b,
    se.b.uni = se,
    p.b.uni = p
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Size Simulation Loop
#' 
#' @param params List of simulation parameters.
#' @return Matrix of point estimates and standard errors for parameters of 
#'   interest from the fitted bivariate regression model.

SizeSimLoop <- function(params) {
  
  # Group sizes.
  sizes <- SampleSizes(n0 = params$n0, mt = params$mT, ms = params$mS)
  
  # Generate data. 
  data <- GenData(
    dist = params$dist,
    n0 = params$n0,
    n1 = sizes$n1,
    n2 = sizes$n2,
    pve_geno = 0,
    pve_pcs = 0.05,
    rho = params$rho,
    snps = params$snps
  )
  
  # If data generation process is non-normal, then rank normalize.
  if (params$dist != "norm") {
    pheno <- data$pheno
    y1 <- pheno[, 1]
    y2 <- pheno[, 2]
    y1[!is.na(y1)] <- RankNorm(y1[!is.na(y1)])
    y2[!is.na(y2)] <- RankNorm(y2[!is.na(y2)])
    pheno <- cbind(y1, y2)
    data$pheno <- pheno
    rm(pheno, y1, y2)
  }
  
  design <- cbind(1, data$covars, data$pcs)
  geno <- data.frame(data$geno)

  # Bivariate results.
  biv_results <- lapply(geno, function(g) {
    out <- BivTestSize(g = g, t = data$pheno[, 1], s = data$pheno[, 2], design = design)
    return(out)
  })
  biv_results <- do.call(rbind, biv_results)
  
  # Univariate results.
  uni_results <- lapply(geno, function(g) {
    out <- UniTestSize(g = g, t = data$pheno[, 1], design = design)
    return(out)
  })
  uni_results <- do.call(rbind, uni_results)
  
  # Output.
  results <- cbind(biv_results, uni_results)
  return(results)
}


# -----------------------------------------------------------------------------
# Power Simulations.
# -----------------------------------------------------------------------------

#' Bivariate Test Power.
#' 
#' @param g Genotype vector.
#' @param bg Genotypic effect.
#' @param t Target outcome.
#' @param s Surrogate outcome.
#' @param design Design matrix.
#' @return Vector of point estimates and standard errors for parameters of 
#'   interest from the fitted bivariate regression model.

BivTestPower <- function(g, bg, t, s, design) {
  
  # Add genotypic effect.
  tg <- t + bg * g
  
  # Fit association model.
  geno_design <- cbind(g, design)
  fit <- SurrogateRegression::Fit.BNR(
    t = tg, 
    s = s, 
    X = geno_design, 
    report = FALSE
  )
  
  # Output
  out <- c(
    point.b.biv = fit@Regression.tab$Point[1] - bg,
    se.b.biv = fit@Regression.tab$SE[1],
    p.b.biv = fit@Regression.tab$p[1]
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Power Calculation for Linear Regression
#' 
#' @param alpha Type I error level. 
#' @param beta Regression coefficient.
#' @param covar_var Marginal variance of the covariate.
#' @param n Sample size.
#' @param resid_var Residual variance of the response given the covariate.
#' @return Numeric estimated power of the two-sided chi-square test.
#' @export 

Power.LinReg <- function(alpha = 0.05, beta, covar_var, n, resid_var) {
  
  # Non-centrality parameter.
  ncp <- n * beta^2 * covar_var / resid_var
  
  # Power.
  power <- stats::pchisq(
    q = stats::qchisq(p = alpha, df = 1, lower.tail = FALSE),
    df = 1,
    ncp = ncp,
    lower.tail = FALSE
  ) 
  
  return(power)
}


#' University Test Power.
#' 
#' @param g Genotype vector.
#' @param bg Genotypic effect.
#' @param t Target outcome.
#' @param design Design matrix.
#' @return p-value for the effect of genotype on the target outcome from
#'   the univariate test.

UniTestPower <- function(g, bg, t, design) {
  
  # Add genotypic effect.
  tg <- t + bg * g
  
  # Subset to subjects with observed target outcome.
  keep <- !is.na(t)
  geno_design <- cbind(g, design)
  
  # Ordinary least squares inference.
  fit <- MatrixOps::fitOLS(y = tg[keep], X = geno_design[keep, ])
  b <- fit$Beta[1]
  se <- sqrt(diag(MatrixOps::matInv(fit$Ibb)))[1]
  df <- sum(keep) - ncol(geno_design)
  p <- 2 * stats::pt(q = abs(b / se), df = df, lower.tail = FALSE)
  
  # Output
  out <- c(
    point.b.uni = b - bg,
    se.b.uni = se,
    p.b.uni = p
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Power Simulation Loop
#' 
#' @param params List of simulation parameters.
#' @return Matrix of p-values for the bivariate and univariate association tests.

PowerSimLoop <- function(params) {
  
  # Group sizes.
  sizes <- SampleSizes(n0 = params$n0, mt = params$mT, ms = params$mS)
  
  # Generate data. 
  data <- GenData(
    dist = params$dist,
    n0 = params$n0,
    n1 = sizes$n1,
    n2 = sizes$n2,
    pve_geno = params$pve,
    pve_pcs = 0.05,
    rho = params$rho,
    snps = params$snps
  )
  
  bg <- data$betas$bg
  design <- cbind(1, data$covars, data$pcs)
  geno <- data.frame(data$geno)

  # Bivariate results.
  biv_results <- lapply(geno, function(g) {
    out <- BivTestPower(
      g = g, 
      bg = bg, 
      t = data$pheno[, 1], 
      s = data$pheno[, 2], 
      design = design
    )
    return(out)
  })
  biv_results <- do.call(rbind, biv_results)
  
  # Univariate results.
  uni_results <- lapply(geno, function(g) {
    out <- UniTestPower(g = g, bg = bg, t = data$pheno[, 1], design = design)
    return(out)
  })
  uni_results <- do.call(rbind, uni_results)
  
  # Output.
  results <- cbind(biv_results, uni_results)
  return(results)
}
