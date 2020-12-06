# Purpose: Plot GTEx results.
# Updated: 2020-12-05

# -----------------------------------------------------------------------------
# Compare effect sizes.
# -----------------------------------------------------------------------------

# Import data. 
data <- readRDS(file = "Data/SSN.rds")
data <- data[complete.cases(data), ]

# Restrict to genome-wide significant.
sig <- subset(
  x = data,
  (Uni.p <= 5e-8) | (Biv.p <= 5e-8)
)

# Counts for the univariate analysis:
cat("Significant in univariate:\n")
sum(sig$Uni.p <= 5e-8)

cat("Genes: ")
length(unique(sig$Ensembl[sig$Uni.p <= 5e-8]))

cat("SNPs: ")
length(unique(sig$SNP[sig$Uni.p <= 5e-8]))
cat("\n")

# Counts for the bivariate analysis:
cat("Significant in bivariate:\n")
sum(sig$Biv.p <= 5e-8)

cat("Genes: ")
length(unique(sig$Ensembl[sig$Biv.p <= 5e-8]))

cat("SNPs: ")
length(unique(sig$SNP[sig$Biv.p <= 5e-8]))
cat("\n")

# Discovery gain.
biv_only <- sum((sig$Biv.p <= 5e-8) & (sig$Uni.p > 5e-8))
uni_only <- sum((sig$Biv.p > 5e-8) & (sig$Uni.p <= 5e-8))
either <- sum((sig$Biv.p <= 5e-8) | (sig$Uni.p <= 5e-8))
cat("Discovery Gain:\n")
round((biv_only - uni_only) / (either) * 100, digits = 1)
cat("\n")

# Average Efficiency gain.
biv_var <- mean(sig$Biv.SE^2)
uni_var <- mean(sig$Uni.SE^2)
cat("Average Efficiency Gain:\n")
round((uni_var / biv_var - 1) * 100, digits = 1)
cat("\n")
