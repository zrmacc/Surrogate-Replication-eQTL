#' Purpose: Count number of subjects with genotype and covariate data in each tissue.
#' Updated: 2021-11-22
library(dplyr)

setwd("~/Documents/Lab/Projects/Surrogates/Cluster/")

# Covariates.
covar <- readRDS(file = "covariates/covar.rds")
covar_ids <- covar$id

# Genotpyes.
geno <- data.table::fread(file = "genotypes/gtex_chr21.fam", select = 1)
geno_ids <- gsub("-", ".", geno$V1)

# Geno and covar.
geno_covar_ids <- intersect(geno_ids, covar_ids)
  
# Blood.
blood <- readRDS(file = "expression/blood_exp.rds")
blood_ids <- colnames(blood)[5:ncol(blood)]
blood_geno_covar_ids <- intersect(blood_ids, geno_covar_ids)

# Muscle.
muscle <- readRDS(file = "expression/mus_exp.rds")
muscle_ids <- colnames(muscle)[5:ncol(muscle)]
muscle_geno_covar_ids <- intersect(muscle_ids, geno_covar_ids)

# Cerebellum.
cere <- readRDS(file = "expression/cere_exp.rds")
cere_ids <- colnames(cere)[5:ncol(cere)]
cere_geno_covar_ids <- intersect(cere_ids, geno_covar_ids)

# SSN.
ssn <- readRDS(file = "expression/ssn_exp.rds")
ssn_ids <- colnames(ssn)[5:ncol(ssn)]
ssn_geno_covar_ids <- intersect(ssn_ids, geno_covar_ids)

# -----------------------------------------------------------------------------

# Counts: SSN + blood.
n0 <- length(intersect(ssn_geno_covar_ids, blood_geno_covar_ids))
n1 <- length(setdiff(blood_geno_covar_ids, ssn_geno_covar_ids))
n2 <- length(setdiff(ssn_geno_covar_ids, blood_geno_covar_ids))
n <- n0 + n1 + n2

# Counts: SSN + muscle.
n0 <- length(intersect(ssn_geno_covar_ids, muscle_geno_covar_ids))
n1 <- length(setdiff(muscle_geno_covar_ids, ssn_geno_covar_ids))
n2 <- length(setdiff(ssn_geno_covar_ids, muscle_geno_covar_ids))
n <- n0 + n1 + n2

# Counts: Cerebellum + muscle.
n0 <- length(intersect(ssn_geno_covar_ids, cere_geno_covar_ids))
n1 <- length(setdiff(cere_geno_covar_ids, ssn_geno_covar_ids))
n2 <- length(setdiff(ssn_geno_covar_ids, cere_geno_covar_ids))
n <- n0 + n1 + n2
