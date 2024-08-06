# Orignally created for the CALICE COllaboration - SiWECAL beam tests
# Generation of beam events and the SiW-ECAL

The generation is based on the [DD4hep](https://github.com/iLCSoft/lcgeo) models for detectors.

*IMPORTANT* There is no pixelization for this simulation: each layer corresponds to one single active area. The pixelization is done a posteriori (see PixelizationProcessor)

## Setup

The software here has been tested with key4hep : `source init_ky4hep.sh` that can be found in the main folder of this code.

## Structure

- `geometry`: Contains the definition of the detector model geometries and configurations.
- `run_scripts`: A set of tools for running the simulations with different beam conditions. Includes scripts for sending batch jobs.


### Visualization

Run 
```bash
geoDisplay file.xml # or any other geometry
```
`To check materials, distances and possible overlappings`:
"materialScan compactgeometryfile.xml x0 y0 z0 x1 y1 z1"
it will display a list of materials moving in a straight line from (x0,y0,z0) to (x1,y1,z1)


## Previous work

This part of the repository is based on [Daniel Jeans' work](https://gitlab.cern.ch/calice/calice_dd4heptestbeamsim/-/tree/master/2017_SiECAL_DESY/).
