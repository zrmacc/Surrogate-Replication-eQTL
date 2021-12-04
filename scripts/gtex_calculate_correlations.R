#' Purpose: Calculate target-surrogate partial correlation, per transcript.
#' The correlation is calculated after residualizing with respect to
#' covariates.
#' Updated: 2021-08-05

# Packages.
suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})


# -----------------------------------------------------------------------------
# Unpack analysis settings.
# -----------------------------------------------------------------------------

# Command line options.
opt_list <- list()

opt <- make_option(c("--chr"), type = "integer", help = "Chromosome", default = 21)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--covar"), type = "character", help = "Covariate file", default = NULL)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--target_exp"), type = "character", help = "Target expression file", default = NULL)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--surrogate_exp"), type = "character", help = "Surrogate expression file", default = NULL)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--out"), type = "character", help = "Output file", default = NULL)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--testing"), type = "logical", help = "Testing?", default = TRUE)
opt_list <- c(opt_list, opt)

# Option parsing.
t0 <- proc.time()
parsed_opts <- optparse::OptionParser(option_list = opt_list)
params <- parse_args(object = parsed_opts)

if (FALSE) {
  params <- list(
    chr = 21,
    covar = "covariates/covar.rds",
    target_exp = "expression/ssn_exp.rds",
    surrogate_exp = "expression/mus_exp.rds",
    out = "results/",
    testing = TRUE
  )
}

# -----------------------------------------------------------------------------

#' Rank Normalize
#' 
#' @param u Numeric vector, possibly including NA.
#' @return Numeric vector. 

RankNorm <- function(u) {
  u[!is.na(u)] <- RNOmni::RankNorm(u[!is.na(u)])
  return(u)
}


#' Extract Tissue Name
#' 
#' @param x String.
#' @return String.

ExtractTissue <- function(x) {
  out <- basename(x)
  out <- gsub(pattern = "(.*)_exp\\.rds", replacement = "\\1", x = out)
  return(out)
}

# -----------------------------------------------------------------------------
# Data import.
# -----------------------------------------------------------------------------

# Covariate data.
all_covars <- readRDS(file = params$covar)

# Target data.
target_exp <- readRDS(file = params$target_exp) %>%
  dplyr::filter(chr == params$chr)

# Surrogate data. 
surrogate_exp <- readRDS(file = params$surrogate_exp) %>%
  dplyr::filter(chr == params$chr)

# Genes expressed in both tissues.
genes <- intersect(target_exp$gene, surrogate_exp$gene)
if (params$testing) {genes <- genes[1:2]}

# -----------------------------------------------------------------------------

# Output name.
out_file <- paste0(params$out, "chr", params$chr, ".txt")

if (!file.exists(out_file)) {
  
  # Loop over genes.
  rho <- lapply(genes, function(tx) {
    
    # ------------------------------------------------------
    # Prepare phenotype matrix.
    # ------------------------------------------------------
    
    target <- target_exp %>% dplyr::filter(gene == tx)
    target <- target %>%
      dplyr::select(-gene, -chr, -tx_start, -tx_stop) %>% 
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "id",
        values_to = "target"
      ) 
    
    surrogate <- surrogate_exp %>%
      dplyr::filter(gene == tx) %>%
      dplyr::select(-gene, -chr, -tx_start, -tx_stop) %>% 
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "id",
        values_to = "surrogate"
      ) 
    
    phenotype <- merge(target, surrogate, by = "id", all = TRUE) %>%
      dplyr::arrange(id) %>%
      dplyr::mutate(
        target = RankNorm(target),
        surrogate = RankNorm(surrogate)
      ) %>%
      dplyr::inner_join(all_covars, by = "id") %>%
      dplyr::filter(complete.cases(.)) %>%
      dplyr::select(-id)
    
    # ------------------------------------------------------
    # Calculate correlations.
    # ------------------------------------------------------
    
    # Residualize.
    fit_target <- lm(target ~ . - surrogate, data = phenotype)
    fit_surrogate <- lm(surrogate ~ . - target, data = phenotype)
    pheno_resid <- data.frame(
      target = resid(fit_target),
      surrogate = resid(fit_surrogate)
    )
    pearson_rho <- cor(pheno_resid$target, pheno_resid$surrogate, method = "pearson")
    spearman_rho <- cor(pheno_resid$target, pheno_resid$surrogate, method = "spearman") 
    
    # Output.
    out <- data.frame(
      chr = params$chr,
      target = ExtractTissue(params$target_exp),
      surrogate = ExtractTissue(params$surrogate_exp),
      gene = tx,
      pearson = pearson_rho,
      spearman = spearman_rho
    )
    return(out)
  })
  
  out <- do.call(rbind, rho)
  data.table::fwrite(
    x = out,
    file = out_file,
    sep = "\t"
  )
}

# -----------------------------------------------------------------------------
# End.
# -----------------------------------------------------------------------------
t1 <- proc.time()
cat(t1-t0)
