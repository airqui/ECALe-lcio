#!/bin/bash

echo -n "Also remove all jobs of ${USER} in the condor queue? (y/n)"
read flag_condor_rm

if [ ${flag_condor_rm} = "y" ] || [ ${flag_condor_rm} = "Y" ]; then
    echo "Proceeding with condor jobs removed."
    condor_rm ${USER}
else
    echo "Proceeding with condor jobs intact.."
fi
rm -r steer
rm -r data
rm -r log
rm -r macros

