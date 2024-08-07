# A. Irles
# IFIC - CSIC/UV
# https://aitanatop.ific.uv.es/aitanatop/luxe/
# https://aitanatop.ific.uv.es/aitanatop/detector-rd/
adrian.irles_at_ific.uv.es
2024/07/31

## Instruction to run the CED event display.

It assumes that you have run the Clustering example.

# DISCLAIMER

For some reason, the libglut3 libraries are not installed in glui01, but you cand download them from (last version) https://freeglut.sourceforge.net/docs/install.php and copy it to your local folder in glui01.ific.uv.es

--> BUT IT WORKS IF YOU usin init_key4hep.sh instead of init_ilcsoft.sh


# INSTALLATION - ilcsoft (please, consider using key4hep instead)

In your /lhome/ or wherever you decide:

> source /cvmfs/ilc.desy.de/sw/x86_64_gcc103_centos7/v02-03-03/init_ilcsoft.sh
>
> tar xzvf freeglut-X.Y.Z.tar.gz (with X.Y.Z being the version of your download)
>
> cd freeglut-X.Y.Z
>
> mkdir build
>
> cd build
>
> cmake ..
>
> make -j3
>
> export LD_LIBRARY_PATH=${PWD}/freeglut-3.6.0/build/lib:$LD_LIBRARY_PATH

--> this should be done everytime, at the same time than source /cvmfs/ilc.desy.de/sw/x86_64_gcc103_
centos7/v02-03-03/init_ilcsoft.sh -> my  advice is that you make a local copy of        the init_ilcsoft.sh in your  own folder and include that line in the init_ilcsoft.sh   (you can run it every time you open a terminal by including it in your ~/.bashrc or just do it manually).


# INSTALLATION - keyforhep 

No need to do anything except

source init_key4hep.sh

# to make it run:

> glced &
>
> Marlin display.xml

