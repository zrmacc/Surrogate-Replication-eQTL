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

# Discovery gain.
biv_only <- sum((sig$Biv.p <= 5e-8) & (sig$Uni.p > 5e-8))
uni_only <- sum((sig$Biv.p > 5e-8) & (sig$Uni.p <= 5e-8))
either <- sum((sig$Biv.p <= 5e-8) | (sig$Uni.p <= 5e-8))
cat("Discovery Gain:\n")
round((biv_only - uni_only) / (either) * 100, digits = 1)

# Average Efficiency gain.
biv_var <- mean(sig$Biv.SE^2)
uni_var <- mean(sig$Uni.SE^2)
cat("Average Efficiency Gain:\n")
round((uni_var / biv_var - 1) * 100, digits = 1)
