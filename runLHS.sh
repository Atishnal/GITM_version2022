#!/bin/bash
# Bash script to run all the cases in reaction_rate.csv in GITM

# Count number of columns
NumCol=$(head -n1 reaction_rates_master.csv |tr '\,' '\n' |wc -l)

for i in $(seq 1 $[$NumCol-1]);do
	NumPad=$(printf "%03d" $i)
	RUNDIR=run$NumPad
make rundir
Col=$[$i+1]
cut -d, -f1,$Col reaction_rates_master.csv > run$NumPad/reaction_rates.csv
cd run$NumPad
qsub rungitm.pbs
cd ..
done



