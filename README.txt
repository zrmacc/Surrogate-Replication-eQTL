# Ensure dependencies in Rscripts/Dependencies.R are installed.

# Type I error and H0 estimation simulations:
sh SizeSim.sh

# Power, relative efficiency, and H1 estimation simulations:
sh PowerSim.sh

# Prepare GTEx figures and tables:
sh GTEx.sh

# SSN eQTL summary statistics are located under Data/SSN.rds, and descriptions of the fields are provided in Data/Description.tsv.