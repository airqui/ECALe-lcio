
#DO THIS ONLY ONCE!!
source ../../init_key4hep.sh
export MARLIN_DLL="$MARLIN_DLL:$PWD/../../PixelizationProcessor/build/lib/libPixelizationProcessor.so"
export MARLIN_DLL="$MARLIN_DLL:$PWD/../../SimCalo2CaloHit/build/lib/libSimCalo2CaloHit.so"
export MARLIN_DLL="$MARLIN_DLL:$PWD/../../NNClusteringProcessor/build/lib/libNNClusteringProcessor.so"
export MARLIN_DLL="$MARLIN_DLL:/lhome/ific/a/airqui/ECALe-Arbor/lib/libRanger.so"


## run the Marlin 
Marlin test.xml
