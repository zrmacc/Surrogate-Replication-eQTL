#!/bin/bash
fin=Configs/PowerConfig.txt

sed 1d ${fin} | while read line
do
	# Read parameter configuration.
	idx=$(echo ${line} | awk "{print \$1}")
	mT=$(echo ${line} | awk "{print \$2}")
	mS=$(echo ${line} | awk "{print \$3}")
	rho=$(echo ${line} | awk "{print \$4}")
	pve=$(echo ${line} | awk "{print \$5}")
	
	# Run simulation.
	Rscript Rscripts/Power.R --n0 1000 --snps 10 --reps 10 --idx ${idx} --out "Test/" --mT ${mT} --mS ${mS} --rho ${rho} --pve ${pve};
done

# Merge results.
Rscript Rscripts/Merge-Concat.R --dir "Test/Power";
Rscript Rscripts/Merge-Mean.R --dir "Test/EstH1";

# Tabulate results.
Rscript Rscripts/Tabulate-Concat.R --dir "Test/EstH1" --out "Tables/EstH1.txt";
Rscript Rscripts/Tabulate-Rejection.R --dir "Test/Power" --out "Tables/Power.txt";

# Calculate relative efficiency.
Rscript Rscripts/Tabulate-Efficiency.R --tab "Tables/EstH1.txt" --out "Tables/RE_H1.txt";

# Plot power curves.
Rscript Rscripts/Plot-Power.R --dir "Test/Power";