# Replication Materials for *Leveraging a surrogate outcome to improve inference on a partially missing target outcome*.
By: Zachary R. McCaw <br>
Updated: 2021-12-04


# Dependencies

Ensure the dependencies under `scripts/dependencies.R` are installed.


# Type I error simulations

`size_sim.sh` provides a working example of simulations under the null hypothesis that genotype has no effect on the target outcome. Running `sh size_sim.sh` as is will perform 1000 simulations, at each of 16 different configurations, with a sample size of 100, tabulate the results, then construct a uniform qq plot for the p-values. The code may be scaled up by setting `--testing FALSE`. 


# Power simulations
`power_sim.sh` provides a working example of simulations under the alternative hypothesis, wherein genotype has an effect on the target outcome. Running `sh power_sim.sh` as is will perform 100 simulations at each of 160 different configurations, with a sample size of 100, tabulate the results, then construct power curves. The code may be scaled up by setting `--testing FALSE`. 


# eQTL analyses

eQTL summary statistics for SNP-transcript pairs with p-values <= 1e-3 are located under `data`. In each case the target tissue is the substantia nigra. The surrogate tissues include whole blood, skeletal muscle, and cerebellum. `data/data_description.tsv` provides a description of the columns of each `.rds` file.
