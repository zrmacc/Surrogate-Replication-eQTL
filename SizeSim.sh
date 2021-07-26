#!/bin/bash
fin=Configs/SizeConfig.txt

sed 1d ${fin} | while read line
do
	# Read parameter configuration.
	idx=$(echo ${line} | awk "{print \$1}")
	mT=$(echo ${line} | awk "{print \$2}")
	mS=$(echo ${line} | awk "{print \$3}")
	rho=$(echo ${line} | awk "{print \$4}")
	

	# Run simulation.
	Rscript Rscripts/Size.R --dist "lognorm" --n0 1000 --snps 10 --reps 10 --idx ${idx} --out "Test/LogNorm/" --mT ${mT} --mS ${mS} --rho ${rho};
done

# Merge results. 
Rscript Rscripts/Merge-Concat.R --dir "Test/TypeIerror";
Rscript Rscripts/Merge-Mean.R --dir "Test/EstH0";

# Tabulate results.
Rscript Rscripts/Tabulate-Concat.R --dir "Test/EstH0" --out "Tables/EstH0.txt";
Rscript Rscripts/Tabulate-Rejection.R --dir "Test/TypeIerror" --out "Tables/Size.txt";

# Calculate relative efficiency.
#Rscript Rscripts/Tabulate-Efficiency.R --tab "Tables/EstH0.txt" --out "Tables/RE_H0.txt";

# Quantile quantile plots.
#Rscript Rscripts/Plot-Size.R --dir "Test/TypeIerror";