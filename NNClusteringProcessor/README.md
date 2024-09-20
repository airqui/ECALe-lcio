# A. Irles
# IFIC - CSIC/UV
# https://aitanatop.ific.uv.es/aitanatop/luxe/
# https://aitanatop.ific.uv.es/aitanatop/detector-rd/
# August 2024

# NNClustering

 File modified from the original MarlinReco NN Cluster:
https://github.com/iLCSoft/MarlinReco/blob/master/Clustering/NNClustering
and
https://github.com/iLCSoft/MarlinReco/blob/master/Analysis/AddClusterProperties/src/AddClusterProperties.cc
 
Example processor that does a simple nearest neighbour (NN) clustering on one or more CalorimeterHit 

Work-in-progress

# Notes
`NNClustering` gives:
- Position of the cluster based on its centre of gravity with the energy deposit as the weight for each hit;
- Direction of the cluster based on its covariance matrix.
These information are obtained using methods in MarlinUtil::WeightedPoints3D.

# Compile:

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
