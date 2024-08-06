# ECALe-lcio
## AITANA - IFIC
https://aitanatop.ific.uv.es/aitanatop/detector-rd/ 
https://aitanatop.ific.uv.es/aitanatop/luxe/

Started in August  2024

## Intro
Code and processors for ECALe for LUXE simulations.
It will be a variation of the SiW-ECAL prototype adapted to LUXE: 15 layers of 2ASUs.
This concept can be used for the ECALe of LUXE and/or the LUXE-NPOD calorimeter

* ASU= Active Signal Unit = 2x2 wafer readout system
  
### generation
Tools for generation of events using ddsim
We generate one single active surface per layer (no pixelization).

### ExampleProcessor
Example Processor to read the SimCalorimeterHits.

### PixelizationProcessor
Processor to read the sim-files and perform proper pixelization and wafer boundary definition for each active  layer.
It assumes the 2xASU geometry.
It takes SimCalorimeterHit and returns SimCalorimeterHit

### MIPCalibration
Work in progress.
It should be run after PixelizationProcessor and use as input the collection created there (not the one from the generation!!).
Takes SimCalorimeterHit and returns a vector of numbers. Should be used only once.
These numbers are the energy2MIP conversion for each layer.

### SimCalo2CaloHit
Work in progress
It requires the output of PixelizationProcessor and knowing the energy2MIP value of each layer. 
It will convert SimCalorimeterHits to CalorimeterHits and will include some level of digitization.
