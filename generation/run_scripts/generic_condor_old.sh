#!/bin/bash

if [ ! -e steer ]; then
    mkdir steer
fi
if [ ! -e data ]; then
    mkdir data
fi
if [ ! -e log ]; then
    mkdir log
fi
if [ ! -e macros ]; then
    mkdir macros 
fi

particle=$1
energy=$2
#conf=$3

ilcsoft_path="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/"
local=$PWD
geometry_folder="${PWD}/../geometry/ECALe_luxe/"
data_path="${local}/data"
steer_path="${local}/steer"
log_path="${local}/log"

#macfile=$4
macfile=grid_-0.0-0.0_${particle}_${energy}GeV.mac
#note: -0-0 are the beam position -x-y in mm

nevt=10000
# nevt=10

cat > ${local}/macros/$macfile <<EOF

/gps/verbose 1
/gps/particle ${particle}
/gps/direction 0 0 1
/gps/pos/type Beam
/gps/pos/shape Circle
/gps/pos/centre 0. 0. -10000. mm
/gps/pos/sigma_x 5 mm
/gps/pos/sigma_y 5 mm

/gps/ene/type Mono
/gps/ene/mono ${energy} GeV
/run/beamOn ${nevt}

EOF


#physl=("QGSP_BERT" "FTFP_BERT")
physl=("QGSP_BERT")

#for energy in ${ens[@]}; do
for physlist in ${physl[@]}; do
  for it in {1..5}; do
    echo $energy $particle $it
    
    label=${physlist}_${particle}_${energy}GeV_${it}
    echo $label
    
    scriptname=runddsim_${label}.py
    condorsh=runddsim_${label}.sh
    condorsub=runddsim_${label}.sub 
    condorfile=runddsim_${label}
    
    #Write the steering file (careful with the compact file name)
    cat > ${local}/steer/$scriptname <<EOF
from DDSim.DD4hepSimulation import DD4hepSimulation
#from SystemOfUnits import mm, GeV, MeV
from g4units import GeV, mm, MeV

SIM = DD4hepSimulation()

SIM.runType = "run"
# Number of events defined in macro file
#SIM.numberOfEvents = $nevt

SIM.skipNEvents = 0
SIM.outputFile = "${data_path}/ECALe_luxe_v0_${label}.slcio"
# SIM.outputFile = "${data_path}/ECALe_luxe_v0_${label}.slcio"

SIM.compactFile = "${geometry_folder}/ECALe_luxe_v0.xml"
SIM.dumpSteeringFile = "${local}/steer/dumpSteering.xml"

SIM.field.eps_min = 0.0001*mm
SIM.part.minimalKineticEnergy = 0.2*MeV
SIM.physicsList = "${physlist}"
SIM.enableDetailedShowerMode=True
EOF
    
    # Here should add mac file creation
    
    #echo "Condor sh ${local}/steer/$condorsh"
    
    cat > ${local}/steer/$condorsh <<EOF
#!/bin/bash
cp -r ${local}/steer/runddsim_${label}.* .
source ${ilcsoft_path}/init_ilcsoft.sh
ddsim --enableG4GPS --macroFile ${local}/macros/${macfile} --steeringFile ${local}/steer/$scriptname
&> ${local}/log/${label}.log
#tar czvf ${local}/ECALe_luxe_v0_${label}.slcio.tar.gz ECALe_luxe_v0_${label}.slcio 
#rm ${local}/log/errors_${condorfile}* ${local}/log/outfile_${condorfile}*
mv ${steer_path}/*$condorfile.txt ${log_path}/.
EOF
    
    cat > ${local}/steer/$condorsub <<EOF
# Unix submit description file
# simple Marlin job
universe = vanilla 
executable              = $condorsh
log                     = $condorfile.log
output                  = outfile_$condorfile.txt
error                   = errors_$condorfile.txt

#should_transfer_files   = Yes
#when_to_transfer_output = ON_EXIT
+JobFlavour = "largo"
#For longer jobs, use other flavours
queue
EOF

    cd ${local}/steer/

    # in2p3 server specific
    chmod +x $condorsh
    # srun -lN1 -t 1-0 --partition=htc ./$condorsh
    # sbatch -N1 -t 1-0 -o ../slurm_out/slurm-%j.out --partition=htc ./$condorsh
    condor_submit $condorsub
    
    cd -
  done
done

