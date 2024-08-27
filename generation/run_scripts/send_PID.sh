#!/bin/bash

#particle="e-"
#particle="mu-"
for particle in "e-"
do
    for energy in 1500 3500 5500 7500 9500 11500 13500 15000
    do
	#source generic_condor.sh $particle $energy $conf "grid_-40-40_"$particle$energy"GeV.mac"
	# sbatch -N1 -t 1-0 -o slurm_out/slurm-%j.out --partition=htc ./generic_condor.sh $particle $energy $conf 
	./generic_condor_single.sh $particle $energy 
	#break
    done
done
