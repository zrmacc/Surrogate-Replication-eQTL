#!/bin/bash
config_file=configs/size_unilateral_10x.tsv
test_dir=test/norm

idx=0
sed 1d ${config_file} | while read line
do
	idx=$(($idx+1))
	echo "Config: ${idx}"
	Rscript scripts/sim_size.R \
		--array_idx ${idx} \
		--config ${config_file} \
		--dist "norm" \
		--testing TRUE \
		--out "${test_dir}";
done

# Merge results. 
Rscript scripts/merge_concat.R --dir "${test_dir}/t1e";
Rscript scripts/merge_mean.R --dir "${test_dir}/est_h0";

# Tabulate.
Rscript scripts/tabulate_rejection.R --dir "${test_dir}/t1e" --out "${test_dir}/t1e.txt";
Rscript scripts/tabulate_concat.R --dir "${test_dir}/est_h0" --out "${test_dir}/est_h0.txt";

# Plot.
Rscript scripts/plot_size.R --dir "${test_dir}/t1e" --out "${test_dir}/qq";