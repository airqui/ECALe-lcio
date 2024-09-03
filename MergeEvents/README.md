# A. Irles
# IFIC - CSIC/UV
# https://aitanatop.ific.uv.es/aitanatop/luxe/
# https://aitanatop.ific.uv.es/aitanatop/detector-rd/
# August 2024

# MergeEvents

Processor that will read collections from two different files and merge all their elements in a single collection (so we can add two particles to the same event) 

Starting point: the Overlay.cc of /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Overlay/v00-22-03/

Work-in-progress

# Compile

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
