# A. Irles
# IFIC - CSIC/UV
# https://aitanatop.ific.uv.es/aitanatop/luxe/
# https://aitanatop.ific.uv.es/aitanatop/detector-rd/
adrian.irles_at_ific.uv.es
2024/07/31


# Example Processor
Simple processor that reads LCIO files with SimCalorimeterHit and writes in screen the energy deposition and coordinates of each hit

# Compile:

> cd XXXXX/ExampleProcessor
> 
> source ../init_key4hep.sh

> mkdir build
> 
> cd build
> 
> cmake ..
> 
> make -j3
> 
> make install

# Run: 
scripts_condor are copied from other example, not usable for us unless few modifications are done.

> cd scripts/
> 
> export MARLIN_DLL="$MARLIN_DLL:$PWD/../lib/libExampleProcessor.so"
> 
(the export MARLIN_DLL is done only once per session, same as the source init_ilcsoft.sh)
> Marlin test.xml
> 
