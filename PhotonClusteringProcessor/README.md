PhotonClusteringProcessor
=========================
August 2024

by **A. Irles** (_IFIC - CSIC/UV_)
and
AITANA group (projects
[LUXE](https://aitanatop.ific.uv.es/aitanatop/luxe/)
[DRD](https://aitanatop.ific.uv.es/aitanatop/detector-rd/))

# NNClustering

 File modified from the original MarlinReco NN Cluster:
https://github.com/iLCSoft/MarlinReco/blob/master/Clustering/NNClustering
and
https://github.com/iLCSoft/MarlinReco/blob/master/Analysis/AddClusterProperties/src/AddClusterProperties.cc
 
Example processor that does a simple nearest neighbour (NN) clustering on one or more CalorimeterHit 

Work-in-progress

## Routines
<!-- 1. Read the MC particle information
    - storing only the  -->
1. Read from the original cluster collection
    - listing the clusters by energy
    - taking the top two as the seed clusters (*NOTE*: this step requires for at least two clusters)
    - connecting the soft clusters to the seeds if demanded
2. Dump out the seed clusters (**strategy 0**)
3. Recluster the hits from the seeds with four kinds of strategies:
    - **Simple plane division** `| c1 | c2 |` A plane at `x==xm` divides all hits into two clusters, where `xm` is the midpoint of the two major clusters
    - **Simple cylindrical collection** `|(c1)|(c2)|` Based on strategy 1, only hits within a cylinder along the centre of cluster are collected (instead of everything in the half ECAL)
    - **Moliere cylindrical collection**: Based on stategy 2, only hits within the radius of 90% of the cluster energy are collected. The radii are upto the point with the two clusters begin to overlap.
    - **Limited cylindrical collection** Based on stategy 3, instead of collecting the whole cylinder, rear layers are discarded if the number of hits in those layers is less than 1/2 of the maximum point

## Tunable parameters


## Connecting


# Usage
```shell
source ../init_key4hep.sh
mkdir build && cd build
cmake ..
make
make install
```