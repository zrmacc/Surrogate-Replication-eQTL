# Purpose: GTEx analysis script.
# Updated: 2021-08-04

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

opt <- make_option(c("--geno"), type = "character", help = "Genotype prefix", default = NULL)
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
  # Example of the structure of params.
  params <- list(
    chr = 21,
    geno = "genotypes/gtex_chr21",
    covar = "covariates/covar.rds",
    target_exp = "expression/ssn_exp.rds",
    surrogate_exp = "expression/cere_exp.rds",
    out = "testing/",
    testing = TRUE
  )
}

# Constants.
cis_range_bp <- 1e6

# -----------------------------------------------------------------------------

#' Rank Normalize
#' 
#' @param u Numeric vector, possibly including NA.
#' @return Numeric vector. 

RankNorm <- function(u) {
  u[!is.na(u)] <- RNOmni::RankNorm(u[!is.na(u)])
  return(u)
}

# -----------------------------------------------------------------------------
# Data import.
# -----------------------------------------------------------------------------

# Genotype meta-data.
bim <- data.table::fread(
  file = paste0(params$geno, ".bim"),
  select = c(1, 2, 4),
  col.names = c("chr", "snp", "bp")
) 

fam <- data.table::fread(
  file = paste0(params$geno, ".fam"), 
  select = c(1),
  col.names = c("id")
) %>%
  dplyr::mutate(
    id = gsub(pattern = "-", replacement = ".", x = id, fixed = TRUE)
  )

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

# Loop over genes.
for (tx in genes) {
  
  # Output name.
  out_file <- paste0(params$out, "chr", params$chr, "_", tx, ".txt")
  if (file.exists(out_file)) {next}
  
  # ------------------------------------------------------
  # Prepare phenotype matrix.
  # ------------------------------------------------------
  
  target <- target_exp %>% dplyr::filter(gene == tx)
  tx_start <- target$tx_start
  tx_stop <- target$tx_stop
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
    dplyr::arrange(id)
  
  # ------------------------------------------------------
  # Prepare covariate matrix.
  # ------------------------------------------------------
  
  covar <- all_covars[(all_covars$id %in% phenotype$id), ]
  covar <- covar[match(covar$id, phenotype$id), ]
  aligned <- all.equal(covar$id, phenotype$id)
  if (!aligned) {
    cat("Phenotype and covariates not aligned.\n")
    next
  }
  covar <- covar %>% dplyr::select(-id)
  
  # ------------------------------------------------------
  # Prepare genotype matrix.
  # ------------------------------------------------------ 
  
  bim_rows <- which(bim$bp >= tx_start - cis_range_bp & bim$bp <= tx_stop + cis_range_bp)
  fam_rows <- which(fam$id %in% phenotype$id)
  if (length(bim_rows) == 0) {
    cat(paste0(tx, " contains no cis-SNPs.\n"))
    next
  }
  
  geno <- ReadPlink::ReadGeno(
    stem = params$geno,
    bim_rows = bim_rows,
    fam_rows = fam_rows
  )
  rownames(geno) <- gsub(pattern = "-", replacement = ".", x = rownames(geno))
  aligned <- all.equal(rownames(geno), phenotype$id)
  if (!aligned) {
    cat("Phenotype and genotype not aligned.\n")
    next
  }
  
  # ------------------------------------------------------
  # Association testing.
  # ------------------------------------------------------ 
  
  results <- lapply(seq_len(ncol(geno)), function(j) {
    g <- geno[, j] 
    
    # Remove subjects with missing genotypes at current SNP.
    is_obs <- !is.na(g)
    if (var(g[is_obs]) == 0) {return(NULL)}
    
    t_obs <- RankNorm(phenotype$target[is_obs])
    s_obs <- RankNorm(phenotype$surrogate[is_obs])
    X_obs <- data.matrix(cbind(g[is_obs], 1, covar[is_obs, ]))
    
    # ------------------------------------------------------ 
    # Bivariate analysis.
    # ------------------------------------------------------ 
    
    biv <- try({
      SurrogateRegression::Fit.BNR(
        t = t_obs,
        s = s_obs,
        X = X_obs,
        Z = X_obs,
        report = FALSE
      )
    }, silent = TRUE)
    
    if (is(biv, "try-error")) {
      biv_out <- rep(NA, 3)
    } else {
      biv_out <- as.numeric(biv@Regression.tab[1, c("Point", "SE", "p")])
    }
    
    biv_out <- data.frame(
      biv_beta = biv_out[1],
      biv_se = biv_out[2],
      biv_p = biv_out[3]
    )
    
    # ------------------------------------------------------ 
    # Univariate analysis.
    # ------------------------------------------------------ 
    
    uni <- try({
      coef(summary(lm(t_obs ~ 0 + X_obs)))
    }, silent = TRUE)
    
    if (is(uni, "try-error")) {
      uni_out <- rep(NA, 3)
    } else {
      uni_out <- as.numeric(uni[1, c("Estimate", "Std. Error", "Pr(>|t|)")])
    }
    
    uni_out <- data.frame(
      uni_beta = uni_out[1],
      uni_se = uni_out[2],
      uni_p = uni_out[3]
    )
    
    # Output
    out <- data.frame(
      chr = params$chr,
      gene = tx,
      snp = colnames(geno)[j],
      biv_out,
      uni_out
    )
    return(out)
  })
  results <- do.call(rbind, results)
  
  # Output.
  if (is(results, "data.frame")) {
    data.table::fwrite(
      x = results,
      file = out_file,
      sep = "\t"
    )
  }
}

# -----------------------------------------------------------------------------
# End.
# -----------------------------------------------------------------------------
t1 <- proc.time()
cat(t1-t0)
