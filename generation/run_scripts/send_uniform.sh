#!/bin/bash

#particle="e-"
#particle="mu-"
for particle in "neutron" "gamma" "pi-"
do
    source uniform_condor.sh $particle 0.5 10
done
