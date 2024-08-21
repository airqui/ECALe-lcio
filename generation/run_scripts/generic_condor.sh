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

# Default options
particle="e-"
energy="1."
nevt=10
nrun0=0
nrun=0
tan_px=0
tan_py=0
gun_x=0
gun_y=0

# Bash options
POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
      echo "USAGE: $0 {particle} {energy in GeV} [-N {number of events} -r {min run id} -R {max run id}]"
      exit 0
      ;;
    -N|--nevt)
      nevt="$2"
      shift # past argument
      shift # past value
      ;;
    -r|--runid)
      nrun0="$2"
      shift # past argument
      shift # past value
      ;;
    -R|--RUNID)
      nrun="$2"
      shift # past argument
      shift # past value
      ;;
    -gx|--gunX)
      gun_x="$2"
      shift # past argument
      shift # past value
      ;;
    -gy|--gunY)
      gun_y="$2"
      shift # past argument
      shift # past value
      ;;
    -ax|--angleX)
      tan_px="$2"
      shift # past argument
      shift # past value
      ;;    
    -ay|--angleY)
      tan_py="$2"
      shift # past argument
      shift # past value
      ;;
    -*|--*)
      echo "Unknown option $1"
    #   echo $( $0 --help )
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done
if [ ${#POSITIONAL_ARGS[@]} -eq 2 ]; then
    particle=${POSITIONAL_ARGS[0]}
    energy=${POSITIONAL_ARGS[1]}
elif [ ${#POSITIONAL_ARGS[@]} -eq 0 ]; then
    echo "No inputs of particle and its energy! Run default options as a test."
else
    echo "Require particle and its energy!"
    # echo $( $0 --help )
    exit -1
fi


# ilcsoft_path="/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/"
run_local=$PWD

source_file="${run_local}/../../init_key4hep.sh"
geometry_folder="${run_local}/../geometry/ECALe_luxe/"
data_path="${run_local}/data"
steer_path="${run_local}/steer"
log_path="${run_local}/log"

energy_in_mev=$( echo "${energy}*1000" | bc )
energy_in_mev=$( printf "%.0f" ${energy_in_mev} )
macfile=ECALe_luxe_${particle}_${energy_in_mev}MeV_x${gun_x}y${gun_y}_ax${tan_px}ay${tan_py}.mac
#note: -0-0 are the beam position -x-y in mm


#Write the G4 mac file
cat > ${run_local}/macros/$macfile <<EOF

/gps/verbose 1
/gps/particle ${particle}
/gps/direction ${tan_px} ${tan_py} 1
/gps/pos/type Beam
/gps/pos/shape Circle
/gps/pos/centre ${gun_x} ${gun_y} -1000. mm
/gps/pos/sigma_x 10 mm
/gps/pos/sigma_y 10 mm

/gps/ene/type Mono
/gps/ene/mono ${energy} GeV
/run/beamOn ${nevt}

EOF


#physl=("QGSP_BERT" "FTFP_BERT")
physl=("QGSP_BERT")

#for energy in ${ens[@]}; do
for physlist in ${physl[@]}; do
  for it in $( eval echo {${nrun0}..${nrun}} ); do
    echo $energy $particle $it
    
    label=${physlist}_${particle}_${energy_in_mev}MeV_${it}
    echo $label
    
    scriptname=runddsim_${label}.py
    condorsh=runddsim_${label}.sh
    condorsub=runddsim_${label}.sub 
    condorfile=runddsim_${label}
    
    #Write the steering file (careful with the compact file name)
    cat > ${run_local}/steer/${scriptname} <<EOF
from DDSim.DD4hepSimulation import DD4hepSimulation
#from SystemOfUnits import mm, GeV, MeV
from g4units import GeV, mm, MeV

SIM = DD4hepSimulation()

SIM.runType = "run"
# Number of events defined in macro file
#SIM.numberOfEvents = ${nevt}

SIM.skipNEvents = 0
SIM.outputFile = "${data_path}/ECALe_luxe_v0_${label}.slcio"

SIM.compactFile = "${geometry_folder}/ECALe_luxe_v0.xml"
SIM.dumpSteeringFile = "${run_local}/steer/dumpSteering.xml"

SIM.field.eps_min = 0.0001*mm
SIM.part.minimalKineticEnergy = 0.2*MeV
SIM.physicsList = "${physlist}"
SIM.enableDetailedShowerMode=True
EOF
    
    # Here should add mac file creation
    
    #echo "Condor sh ${run_local}/steer/$condorsh"
    
    cat > ${run_local}/steer/${condorsh} <<EOF
#!/bin/bash

# cp -r ${run_local}/steer/runddsim_${label}.* .
source ${source_file}
ddsim --enableG4GPS --macroFile ${run_local}/macros/${macfile} --steeringFile ${run_local}/steer/${scriptname}
#&> ${run_local}/log/${label}.log
#tar czvf ${run_local}/ECALe_luxe_v0_${label}.slcio.tar.gz ECALe_luxe_v0_${label}.slcio 
# mv ${steer_path}/*${condorfile}.log ${log_path}/
#rm ${run_local}/log/errors_${condorfile}* ${run_local}/log/outfile_${condorfile}*
EOF
    
    cat > ${run_local}/steer/$condorsub <<EOF
# Unix submit description file
# simple Marlin job
universe = vanilla 
executable              = ${run_local}/steer/$condorsh
log                     = ${run_local}/log/condor_$condorfile.log
output                  = ${run_local}/log/outfile_$condorfile.log
error                   = ${run_local}/log/errors_$condorfile.log

#should_transfer_files   = Yes
#when_to_transfer_output = ON_EXIT
+JobFlavour = "largo"
#For longer jobs, use other flavours
queue
EOF

    cd ${run_local}/steer/

    # in2p3 server specific
    chmod +x $condorsh
    # srun -lN1 -t 1-0 --partition=htc ./$condorsh
    # sbatch -N1 -t 1-0 -o ../slurm_out/slurm-%j.out --partition=htc ./$condorsh
    condor_submit $condorsub
    
    cd -
  done
done

