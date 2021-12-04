#!/bin/bash
config_file=configs/power_unilateral_10x.tsv
test_dir=test/norm

idx=0
sed 1d ${config_file} | while read line
do
	idx=$(($idx+1))
	echo "Config: ${idx}"
	Rscript scripts/sim_power.R \
		--array_idx ${idx} \
		--config ${config_file} \
		--dist "norm" \
		--testing TRUE \
		--out "${test_dir}";
done

# Merge results.
Rscript scripts/merge_concat.R --dir "${test_dir}/power";
Rscript scripts/merge_mean.R --dir "${test_dir}/est_h1";

# Tabulate.
Rscript scripts/tabulate_rejection.R --dir "${test_dir}/power" --out "${test_dir}/power.txt";
Rscript scripts/tabulate_concat.R --dir "${test_dir}/est_h1" --out "${test_dir}/est_h1.txt";

# Plot.
Rscript scripts/plot_power.R --dir "${test_dir}/power" --out "${test_dir}/power";