
#DO THIS ONLY ONCE!!
source ../../init_ilcsoft.sh
export MARLIN_DLL="$MARLIN_DLL:$PWD/../build/lib/libPixelizationProcessor.so"


## run the Marlin 
Marlin test.xml