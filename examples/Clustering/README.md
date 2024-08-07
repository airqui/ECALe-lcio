# A. Irles
# IFIC - CSIC/UV
# https://aitanatop.ific.uv.es/aitanatop/luxe/
# https://aitanatop.ific.uv.es/aitanatop/detector-rd/
adrian.irles_at_ific.uv.es
2024/08/07

## Instruction to run the Full Chain until the reconstruction of clusters

It assumes that you have run the Clustering example.

Assumes that you have generated a file with the geometry described in ../../generation

Processors:

    <processor name="MyPixelizationProcessor"/> 
    <processor name="MySimCalo2CaloHit"/>
    <processor name="MyNNClustering"/>
    <processor name="DSTOutput"/>

# DISCLAIMER

Follow the instructions in
../generation
../PixelizationProcessor
../SimCalo2CaloHit
../NNClustering

before starting with this example.