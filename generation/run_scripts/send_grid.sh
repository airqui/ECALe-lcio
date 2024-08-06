#!/bin/bash

#particle="e-"
#particle="mu-"
for particle in "pi-" "mu-" "gamma"
do
    for energy in 2 50 100
    do
	#source generic_condor.sh $particle $energy $conf "grid_-40-40_"$particle$energy"GeV.mac"
	# sbatch -N1 -t 1-0 -o slurm_out/slurm-%j.out --partition=htc ./generic_condor.sh $particle $energy $conf 
	./generic_condor.sh $particle $energy 
	#break
    done
done
