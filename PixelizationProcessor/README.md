# A. Irles
# IFIC - CSIC/UV
# https://aitanatop.ific.uv.es/aitanatop/luxe/
# https://aitanatop.ific.uv.es/aitanatop/detector-rd/
# August 2024

# PixelizationProcessor

The ECALe layers have an internal structure due to the wafer layout.
They have 305um of dead space at the edges of each sensor and there is a 10um gap between pixels.
In addition, we should assume some tolerance in the separation between sensors. 100um seems reasonable.
Ideally we would tune the simulation to not have silicon in these 100um, for the moment we just ignore the hits in that region but not change the material content.

The pixelization assumes tha tthe simulations gives a long layer of pure silicon (one single readout-channel of the full surface of the layer).
For LUXE, we are assuming the ECALe to have a size of 4 wafers (in X) and 2 wafers (in Y), i.e. 2 ASUs

# Compile:

> cd XXXXX/PixelizationProcessor
>
> source init_ilcsoft.sh

(or use an updated one from /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03 or similar folders!)

> mkdir build
> 
> cd build
> 
> cmake -C $ILCSOFT/ILCSoft.cmake ..
> 
> make -j3
> 
> make install

# Run: 
scripts_condor are copied from other example, not usable for us unless few modifications are done.

> cd scripts/
> 
> export MARLIN_DLL="$MARLIN_DLL:$PWD/../lib/libPixelizationProcessor.so"
> 
(the export MARLIN_DLL is done only once per session, same as the source init_ilcsoft.sh)
> Marlin test.xml
> 
