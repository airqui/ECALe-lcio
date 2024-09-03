# A. Irles
# IFIC - CSIC/UV
# https://aitanatop.ific.uv.es/aitanatop/luxe/
# https://aitanatop.ific.uv.es/aitanatop/detector-rd/
adrian.irles_at_ific.uv.es
2024/09/03

## Instruction to run the a processor to merge two events from different files: we merge all elements for collections with the same name
## fromt o different files: the default input file (in global) and the  InputFileNames in the processor input

Assumes that you have generated a file with the geometry described in ../../generation

Processors:
    
    <processor name="MyMergeEvents"/>
    <processor name="DSTOutput"/>

# DISCLAIMER

processor created from the original ILCSOft/Overlay
