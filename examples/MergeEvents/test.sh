#DO THIS ONLY ONCE!!
source ../../init_key4hep.sh
export MARLIN_DLL="$MARLIN_DLL:$PWD/../../MergeEvents/build/lib/libMergeEvents.so"

## run the Marlin 
Marlin test.xml