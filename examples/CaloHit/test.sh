
#DO THIS ONLY ONCE!!
source ../../init_ilcsoft.sh
export MARLIN_DLL="$MARLIN_DLL:$PWD/../../PixelizationProcessor/build/lib/libPixelizationProcessor.so"
export MARLIN_DLL="$MARLIN_DLL:$PWD/../../SimCalo2CaloHit/build/lib/libSimCalo2CaloHit.so"

## run the Marlin 
Marlin test.xml