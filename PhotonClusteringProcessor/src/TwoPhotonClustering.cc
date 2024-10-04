#include "TwoPhotonClustering.hh"
#include <iostream>
#include "WeightedPoints3D.h"

#include "time.h"

#include <algorithm>
// fitting stuff
#include <TMath.h>
#include <Math/Point2D.h>
#include <Math/Vector2D.h>
#include <TROOT.h>


using namespace lcio;
using namespace marlin;
using namespace std;

TwoPhotonClustering aTwoPhotonClustering;

TwoPhotonClustering::TwoPhotonClustering() : Processor("TwoPhotonClustering") {

    // modify processor description
    _description = "TwoPhotonClustering : simple nearest neighbour clustering";

    registerInputCollection(LCIO::CALORIMETERHIT,
                            "CaloHitInputCollection",
                            "Name of the input CaloHit collection",
                            _CaloHitCol,
                            std::string("EcalCollection"));

    registerInputCollection(LCIO::CLUSTER,
                            "NNClusterInputCollection",
                            "Name of the input NNCluster collection",
                            _NNClusterCol,
                            std::string("NNClusters"));

    registerOutputCollection(LCIO::CLUSTER,
                            "OutputCollection",
                            "Name of the output collections",
                            _outputColName,
                            std::string("TwoPhotonClusters"));

    // For the moment keep this as an external input.
    /*registerProcessorParameter("NPhotons",
                                 "This is the number of photons that are expected",
                                 _nPhotonstoReconstruct,
                                 (int)2);*/

    registerProcessorParameter("DoReclustering",
                               "dedice if you want to merge the soft clusters with the main one",
                               _doRecluster,
                               (bool)false);
    registerProcessorParameter("Strategy",
                               "strategy to follow to separate clusters and measure direction",
                               _strategytofollow,
                               (int)4);

    registerProcessorParameter("RadiusFactorClusteringCylinder",
                               "Factor between 0 and 1 dictates the size of clustering cylinder", // 1 = two cylinders touching at the midpoint
                               DISTANCE_RATIO,
                               (float) 1);
    registerProcessorParameter("MolierePercentage",
                               "Factor between 0 and 1 dictates the percentage of energy recovered for a cluster",
                               ENERGY_FACTOR,
                               (float) 0.7);

    registerProcessorParameter("LastLayer",
                               "Last layer to be considered... WIP!!",
                               _lastlayer,
                               (int) 10);

    registerProcessorParameter("AnalogueFit",
                                "simplistic implementation of fit with weighted by energy hits",
                                _boolAnalogue,
                                (bool) false);
}


void TwoPhotonClustering::init() {
    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;

    AIDAProcessor::tree(this);
    _x_res_origin = new TH1D("x_residual_vertex","Residual of X-position at initial Vtx; x_{MC}-x_{cl} [mm]; N", 2001, -1001, 1001);
    _y_res_origin = new TH1D("y_residual_vertex","Residual of Y-position at initial Vtx; y_{MC}-y_{cl} [mm]; N", 2001, -1001, 1001);
    _r_res_origin = new TH1D("r_residual_vertex","Residual of XY-position at initial Vtx; #sqrt{(x_{MC}-x_{cl})^{2}+(y_{MC}-y_{cl})^{2}}) [mm]; N", 2001, -1001, 1001);

    //_x_res_layer0 = new TH1D("x_residual_layer0","Residual of X-position at  the front of the ECAL; x_{MC}-x_{cl} [mm]; N", 2001, -1001, 1001);
    //_y_res_layer0 = new TH1D("y_residual_layer0","Residual of Y-position at  the front of the ECAL; y_{MC}-y_{cl} [mm]; N", 2001, -1001, 1001);
    //_r_res_layer0 = new TH1D("r_residual_layer0","Residual of XY-position at  the front of the ECAL; #sqrt{(x_{MC}-x_{cl})^{2}+(y_{MC}-y_{cl})^{2}}) [mm]; N", 20001, -1001, 1001);

}

void TwoPhotonClustering::processRunHeader(LCRunHeader * /*run*/) {
    _nRun++;
}

void TwoPhotonClustering::processEvent(LCEvent *evt) {
    clock_t start = clock(); // Implementation of real time consumption
    streamlog_out(MESSAGE0) << " ---- TwoPhotonClustering::processEvent() - evt: "
                            << evt->getRunNumber() << ", " << evt->getEventNumber()
                            << std::endl;

    IMPL::LCCollectionVec *output_cluster = new IMPL::LCCollectionVec(LCIO::CLUSTER);
    
    /* **READ FROM INPUT CLUSTER COLLECTION** */
    LCCollection *input_caloclusters = evt->getCollection(_NNClusterCol);
    int n_clusters = input_caloclusters->getNumberOfElements();
    LCCollection *input_calohits = evt->getCollection(_CaloHitCol);
    int n_hits = input_calohits->getNumberOfElements();
    CellIDDecoder<EVENT::CalorimeterHit> cd(input_calohits);


    if (n_clusters < 2) {
        streamlog_out(MESSAGE0) << "Only found "<< n_clusters << "cluster(s) - fewer than required 2" << std::endl
                                << "Will dump an empty collection in this event." <<std::endl;
    } else {
        /* **Begin mapping the clusters by energy** */
        std::map<float, Cluster *> input_clusters_map;

        for (int i = 0; i < n_clusters; i++) {
            Cluster *clus = dynamic_cast<Cluster *>(input_caloclusters->getElementAt(i));
            float energy_clus = clus->getEnergy();
            // fill the map
            auto it = input_clusters_map.find(energy_clus);
            if (it == input_clusters_map.end()) { // cannot find a match in the mapping
                input_clusters_map.insert(pair<float, Cluster *>(energy_clus, clus));
                // it->second.push_back(px);
            } else { // found a match in the mapping, suprise suprise!
                streamlog_out(ERROR) << " WHAT ??? there are two clusters with same energy ???? i:" << i << " e:" << energy_clus << std::endl;
                input_clusters_map.insert(pair<float, Cluster *>(energy_clus+0.001, clus)); // Add the cluster anyway but under a different key
                // TODO: Considering to put the previous line with a do-while-loop,
                // just in the case there would be three clusters with the same energy...
            }
        }
        // sorting the map of clusters according to energy, in descending mode
        std::vector<pair<float, Cluster *>> pairs;
        for (auto &it : input_clusters_map) {
            pairs.push_back(it);
        }
        sort(pairs.begin(), pairs.end(), [](auto &a, auto &b)
                 { return a.first > b.first; }); // Sorting from large to small
        /* **Finish mapping the clusters by energy** */

        // fill two vector of clusters:_ the seed ones and the soft ones.
        std::vector<Cluster *> seed_clusters; // First two clusters
        std::vector<Cluster *> soft_clusters; // All the rest ones
        int counter = 0;
        for (auto &pair : pairs) {
            if (counter < 2) {
                seed_clusters.push_back(pair.second);
            } else {
                soft_clusters.push_back(pair.second);
            }
            counter++;
        }
        streamlog_out(MESSAGE) << "N main clusters: "<<seed_clusters.size()<<", rest of clusters: "<<soft_clusters.size() <<endl;
        streamlog_out(MESSAGE) << "energies: "<< seed_clusters.at(0)->getEnergy()<<" "<<seed_clusters.at(1)->getEnergy();
        for(int ix=0; ix<soft_clusters.size(); ix++) streamlog_out(MESSAGE)<<" "<<soft_clusters.at(ix)->getEnergy();
        streamlog_out(MESSAGE) <<endl;

        float r0_90;
        float r1_90;
        //The FindZEMShower needs to be revised, it returns always a value of 10000 
        //float last_layer_r0_90=FindZEndShower(input_calohits,_x0,_y0,r0_90); // last layer with 1/2 of the maximum hits
        //float last_layer_r1_90=FindZEndShower(input_calohits,_x1,_y1,r1_90);
        float last_layer_r0_90;//in units of mm, assuming 15 layers... REVISE THIS
        float last_layer_r1_90;
            
            
        // START WITH THE FIT TESTS
        // FOR THE MOMENT THIS IS DONE ASSUMING 2 seed CLUSTERS
        const int n = n_hits; // Number of total hits
        int n0 = 0;
        int n1 = 0;
        /* Change the following lines to double because of WeightedPoints3D,
           hopefully will not create a fuss... */
        double x0[n], y0[n], z0[n], e0[n]; // all hits near cluster seed 0
        double x1[n], y1[n], z1[n], e1[n]; // all hits near cluster seed 1

        // std::vector<float> two_cluster_coordinates=GetCenterCoordinates(seed_clusters,soft_clusters);
        // float _x0=two_cluster_coordinates.at(0);
        // float _y0=two_cluster_coordinates.at(1);
        // float _e0=two_cluster_coordinates.at(2);
        // float _x1=two_cluster_coordinates.at(3);
        // float _y1=two_cluster_coordinates.at(4);
        // float _e1=two_cluster_coordinates.at(5);

        std::vector<ClusterInfo> twoPhotonClusters = GetTwoClustersArray(seed_clusters, soft_clusters);
        double cluster_e_ratio = twoPhotonClusters[0].energy / twoPhotonClusters[1].energy;
        if (cluster_e_ratio < 0.5 || cluster_e_ratio > 2.0) {
            std::cout<<"NNClustering Failed! Energy difference is larger than 2 (should be 1)..." <<endl;
            return;
        }
        /* **STRATEGY 0** GOES OUT DIRECTLY */
        if (_strategytofollow != 0) {
            float _x0 = twoPhotonClusters[0].position[0];
            float _y0 = twoPhotonClusters[0].position[1];
            float _e0 = twoPhotonClusters[0].energy;
            float _x1 = twoPhotonClusters[1].position[0];
            float _y1 = twoPhotonClusters[1].position[1];
            float _e1 = twoPhotonClusters[1].energy;
            
            float xmean = (_x0 + _x1) / 2.;
            float ymean = (_y0 + _y1) / 2.;
            float rmean = sqrt((_x0 - xmean) * (_x0 - xmean) + (_y0 - ymean) * (_y0 - ymean));
   
            // FIRST RECLUSTER OF HITS...  for a simple finding of the center of the shower layer per layer.. 
            // so we can play later with the moliere radius in a cleaner way -- strategy 3 and 4
            int nn0[40]={0};
            int nn1[40]={0};
            double xx0[40][n]={0},yy0[40][n]={0},zz0[40][n]={0},ee0[40][n]={0};
            double xx1[40][n]={0},yy1[40][n]={0},zz1[40][n]={0},ee1[40][n]={0};
            float _x0_layer[20]={0};
            float _y0_layer[20]={0};
            float _x1_layer[20]={0};
            float _y1_layer[20]={0};
            float last_layer_r0_90;
            float last_layer_r1_90;

            if (_strategytofollow == 3 || _strategytofollow ==4) {

                for (int i = 0; i < n_hits; i++) {
                    // get all simcalorimeter hits
                    CalorimeterHit *hit = dynamic_cast<CalorimeterHit *>(input_calohits->getElementAt(i));
                    float position[3];
                    position[0] = hit->getPosition()[0];
                    position[1] = hit->getPosition()[1];
                    position[2] = hit->getPosition()[2];

                    int k = cd(hit)["K"]-1;

                    if (fabs(position[0] - _x0) < fabs(position[0] - _x1))
                    {
                        xx0[k][nn0[k]] = position[0];
                        yy0[k][nn0[k]] = position[1];
                        zz0[k][nn0[k]] = position[2];
                        ee0[k][nn0[k]] = hit->getEnergy();
                        nn0[k]++;
                    }
                    else
                    {
                        xx1[k][nn1[k]] = position[0];
                        yy1[k][nn1[k]] = position[1];
                        zz1[k][nn1[k]] = position[2];
                        ee1[k][nn1[k]] = hit->getEnergy();
                        nn1[k]++;
                    }

                }
                for(int ilayer=0; ilayer<15; ilayer ++) {
                    int ntemp=nn0[ilayer];
                    double xtemp0[n],ytemp0[n], ztemp0[n], etemp0[n];
                    for(int in=0; in<ntemp; in++) {xtemp0[in]=xx0[ilayer][in]; ytemp0[in]=yy0[ilayer][in]; ztemp0[in]=zz0[ilayer][in]; etemp0[in]=ee0[ilayer][in];}
                    std::vector<float> _xy0_layer= WeightedCenter( ntemp, xtemp0, ytemp0, ztemp0, etemp0);
                    if( isnan( _xy0_layer.at(0)  ) || isnan( _xy0_layer.at(1)  )) {
                        _x0_layer[ilayer]=_x0;
                        _y0_layer[ilayer]=_y0;
                    } else {
                        _x0_layer[ilayer]=_xy0_layer.at(0);
                        _y0_layer[ilayer]=_xy0_layer.at(1);
                    }
                    
                    ntemp=nn1[ilayer];
                    double xtemp1[n],ytemp1[n], ztemp1[n], etemp1[n];
                    for(int in=0; in<ntemp; in++) {xtemp1[in]=xx1[ilayer][in]; ytemp1[in]=yy1[ilayer][in]; ztemp1[in]=zz1[ilayer][in]; etemp1[in]=ee1[ilayer][in];}
                    std::vector<float> _xy1_layer= WeightedCenter( ntemp, xtemp1, ytemp1, ztemp1, etemp1);                    
                    if( isnan( _xy1_layer.at(0)  ) || isnan( _xy1_layer.at(1)  )) {
                        _x1_layer[ilayer]=_x1;
                        _y1_layer[ilayer]=_y1;
                    }else {
                        _x1_layer[ilayer]=_xy1_layer.at(0);
                        _y1_layer[ilayer]=_xy1_layer.at(1);
                    }
                }

                r0_90=FindRadius90(input_calohits,_x0_layer,_y0_layer,rmean,_e0);
                r1_90=FindRadius90(input_calohits,_x1_layer,_y1_layer,rmean,_e1);

                //The FindZEMShower needs to be revised !!!
                last_layer_r0_90=FindZEndShower(input_calohits,_x0_layer,_y0_layer,r0_90); // last layer with 1/2 of the maximum hits
                last_layer_r1_90=FindZEndShower(input_calohits,_x1_layer,_y1_layer,r1_90);

                cout<<"Moliere: "<<rmean*DISTANCE_RATIO<<"\t"<<r0_90<<"\t"<<r1_90<<endl;
                cout<<"Clusters: "<<_e0*0.0141<<"\t"<<_e1*0.0141<<endl;
                cout<<"Calculated : LastLayers, in z-mm: "<<last_layer_r0_90<<"\t"<<last_layer_r1_90<<endl;
                last_layer_r0_90=15.*float(_lastlayer);//in units of mm, assuming 15 layers... REVISE THIS
                last_layer_r1_90=15.*float(_lastlayer);
                cout<<"Used (by steering file): LastLayers, in z-mm: "<<last_layer_r0_90<<"\t"<<last_layer_r1_90<<endl;
            }


            /* **RECLUSTER THE HITS** */
            for (int i = 0; i < n_hits; i++) {
                // get all simcalorimeter hits
                CalorimeterHit *hit = dynamic_cast<CalorimeterHit *>(input_calohits->getElementAt(i));
                float position[3];
                position[0] = hit->getPosition()[0];
                position[1] = hit->getPosition()[1];
                position[2] = hit->getPosition()[2];
                int k = cd(hit)["K"]-1;

                // ********************** strategy 1, I separate the two clusters with a plane in x=xmean between clus0 and clus1
                if (_strategytofollow == 1)
                {
                    if (fabs(position[0] - _x0) < fabs(position[0] - _x1))
                    {
                        x0[n0] = position[0];
                        y0[n0] = position[1];
                        z0[n0] = position[2];
                        e0[n0] = hit->getEnergy();
                        n0++;
                    }
                    else
                    {
                        x1[n1] = position[0];
                        y1[n1] = position[1];
                        z1[n1] = position[2];
                        e1[n1] = hit->getEnergy();
                        n1++;
                    }
                } 

                // ********************** strategy 2, I separate the two clusters with a cylinder of rmean = 2
                if (_strategytofollow == 2)
                {
                    if (sqrt( (position[0] - _x0)*(position[0] - _x0)+(position[1] - _y0)*(position[1] - _y0)) < rmean*DISTANCE_RATIO)
                    {
                        x0[n0] = position[0];
                        y0[n0] = position[1];
                        z0[n0] = position[2];
                        e0[n0] = hit->getEnergy();
                        n0++;
                    }
                    if( sqrt( (position[0] - _x1)*(position[0] - _x1)+(position[1] - _y1)*(position[1] - _y1)) < rmean*DISTANCE_RATIO)
                    {
                        x1[n1] = position[0];
                        y1[n1] = position[1];
                        z1[n1] = position[2];
                        e1[n1] = hit->getEnergy();
                        n1++;
                    }
                } 

                if (_strategytofollow==3 )
                {
                    // TO BE STUDIED: if I usw the layer per layer center position, it does not work... why ???
                    // if (sqrt( (position[0] - _x0_layer[k])*(position[0] - _x0_layer[k])+(position[1] - _y0_layer[k])*(position[1] - _y0_layer[k]))< r0_90 )  
                    if (sqrt( (position[0] - _x0)*(position[0] - _x0)+(position[1] - _y0)*(position[1] - _y0))< r0_90 )  
                    {
                        x0[n0] = position[0];
                        y0[n0] = position[1];
                        z0[n0] = position[2];
                        e0[n0] = hit->getEnergy();
                        n0++;
                    }
                    //if( sqrt( (position[0] - _x1_layer[k])*(position[0] - _x1_layer[k])+(position[1] - _y1_layer[k])*(position[1] - _y1_layer[k]) )<r1_90 )
                    if( sqrt( (position[0] - _x1)*(position[0] - _x1)+(position[1] - _y1)*(position[1] - _y1) )< r1_90)
                    {
                        x1[n1] = position[0];
                        y1[n1] = position[1];
                        z1[n1] = position[2];
                        e1[n1] = hit->getEnergy();
                        n1++;
                    }
                } 

                // ********************** strategy 4, we use the cylinder that contains the 90% of the energy and we remove the last layers with only few hits
                if (_strategytofollow == 4)
                {
                    // if (sqrt( (position[0] - _x0_layer[k])*(position[0] - _x0_layer[k])+(position[1] - _y0_layer[k])*(position[1] - _y0_layer[k]))< r0_90 && position[2]<last_layer_r0_90)  
                    if (sqrt( (position[0] - _x0)*(position[0] - _x0)+(position[1] - _y0)*(position[1] - _y0))< r0_90 && position[2]<last_layer_r0_90)  
                    {
                        x0[n0] = position[0];
                        y0[n0] = position[1];
                        z0[n0] = position[2];
                        e0[n0] = hit->getEnergy();
                        n0++;
                    }
                    //if( sqrt( (position[0] - _x1_layer[k])*(position[0] - _x1_layer[k])+(position[1] - _y1_layer[k])*(position[1] - _y1_layer[k]) )<r1_90 && position[2]<last_layer_r1_90)
                    if( sqrt( (position[0] - _x1)*(position[0] - _x1)+(position[1] - _y1)*(position[1] - _y1) )< r1_90 && position[2]<last_layer_r1_90)
                    {
                        x1[n1] = position[0];
                        y1[n1] = position[1];
                        z1[n1] = position[2];
                        e1[n1] = hit->getEnergy();
                        n1++;
                    }
                }
            }

 
            for (int iphotons = 0; iphotons < 2; iphotons++) {

                std::vector<float> position_direction_clustered_photon;
                if(iphotons==0) position_direction_clustered_photon=GetFittedCoordinates(n0,x0,y0,z0,e0);
                else position_direction_clustered_photon=GetFittedCoordinates(n1,x1,y1,z1,e1);
                // twoPhotonClusters[iphotons].theta = cl_direction_polar.Theta();
                for(int j=0; j<3; j++) twoPhotonClusters[iphotons].position[j]=position_direction_clustered_photon.at(j);
                twoPhotonClusters[iphotons].theta = position_direction_clustered_photon.at(3);
                twoPhotonClusters[iphotons].phi = position_direction_clustered_photon.at(4);

            }
        }
        ClusterImpl* newClusters[2];
        for (int iphotons = 0; iphotons < 2; iphotons++) {
            newClusters[iphotons] = new ClusterImpl(); // new object
        
            newClusters[iphotons]->setEnergy(twoPhotonClusters[iphotons].energy);
            // newClusters[iphotons]->setEnergyError(twoPhotonClusters[iphotons].error_energy);
            newClusters[iphotons]->setPosition(twoPhotonClusters[iphotons].position);
            // newClusters[iphotons]->setPositionError(twoPhotonClusters[iphotons].error_position);
            newClusters[iphotons]->setITheta(twoPhotonClusters[iphotons].theta);
            newClusters[iphotons]->setIPhi(twoPhotonClusters[iphotons].phi);
            // newClusters[iphotons]->setDirectionError(twoPhotonClusters[iphotons].error_direction);
            
            output_cluster->addElement(newClusters[iphotons]);
        }
    }
    
    evt->addCollection( output_cluster , _outputColName ) ;
    _nEvt++;




    clock_t end = clock();
    streamlog_out(MESSAGE0) << " ---- PhotonClusteringProcessor::processEvent() - time: "
                            << double(end - start) / double(CLOCKS_PER_SEC)
                            << std::endl;
}

void TwoPhotonClustering::check(LCEvent * /*evt*/)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void TwoPhotonClustering::end()
{

    //   streamlog_out(MESSAGE0) << "NNClusterProcessor::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}

std::vector<float> TwoPhotonClustering::GetCenterCoordinates(std::vector<Cluster*> seed_clusters, std::vector<Cluster*> soft_clusters) {
/*  Return the central coordinates (2D) of the two main clusters,
    with the option of reclustering, where the `soft_clusters` are required.
    Needs further clean-up. */

    float _x0 = seed_clusters.at(0)->getPosition()[0];
    float _y0 = seed_clusters.at(0)->getPosition()[1];
    // float _z0 = seed_clusters.at(0)->getPosition()[2];
    float _e0 = seed_clusters.at(0)->getEnergy();
    float _e0sum = seed_clusters.at(0)->getEnergy();
    // float _e0invsum=1./_e0;
    auto cluster_0 = ROOT::Math::XYPoint(_x0, _y0);

    float _x1 = seed_clusters.at(1)->getPosition()[0];
    float _y1 = seed_clusters.at(1)->getPosition()[1];
    // float _z1 = seed_clusters.at(1)->getPosition()[2];
    float _e1 = seed_clusters.at(1)->getEnergy();
    float _e1sum = seed_clusters.at(1)->getEnergy();
    // float _e1invsum=1./_e1;
    auto cluster_1 = ROOT::Math::XYPoint(_x1, _y1);
    
    // Geometric centre of the two clusters
    // float xmean = (_x0 + _x1) / 2.;
    // float ymean = (_y0 + _y1) / 2.;
    // float rmean = sqrt((_x0 - xmean) * (_x0 - xmean) + (_y0 - ymean) * (_y0 - ymean));
    float dmean = (cluster_0-cluster_1).R(); // dmean == 2*rmean, just for clarification

    //If _doRecluster == true and we have more than 0 soft-clusters, we recluster them if they are close to the main cluster.
    // ad-hoc criteria: close is defined as in a radius (x-y) smaller than half the distance between the two main clusters.
    // shall we tune this? shall we add more criteria like energy cuts ?
    float RECLUSTER_RATIO = 0.5; // Should be a factor <1/2, otherwise will cause dilemma of soft cluster assignment

    // Centre of a cluster is the weighted average of its hits
    if(soft_clusters.size()>0 && _doRecluster==true) {
        // Weighted positions of the clusters
        _x0 *= _e0;
        _y0 *= _e0;
        // _z0 *= _e0;
        _x1 *= _e1;
        _y1 *= _e1;
        // _z1 *= _e1;
        

        for(int isoft=0; isoft<soft_clusters.size(); isoft++) {
            float position[3];
            position[0] = soft_clusters.at(isoft)->getPosition()[0];
            position[1] = soft_clusters.at(isoft)->getPosition()[1];
            // position[2] = soft_clusters.at(isoft)->getPosition()[2];
            auto soft_position = ROOT::Math::XYPoint(position[0], position[1]);
            float _e = soft_clusters.at(isoft)->getEnergy();

            if ((cluster_0 - soft_position).R() < dmean*RECLUSTER_RATIO) {
                _x0 += position[0] * _e;
                _y0 += position[1] * _e;
                // _z0 += position[2] * _e;
                // _e0 += _e;
                // _e0invsum += 1./_e;
                _e0sum += _e;
            }
            if ((cluster_1 - soft_position).R() < dmean*RECLUSTER_RATIO) {
                _x1 += position[0] * _e;
                _y1 += position[1] * _e;
                // _z1 += position[2] * _e;
                // _e1 += _e;
                // _e1invsum+=1./_e;
                _e1sum += _e;
            }
            
        }


        _x0 /= _e0sum;
        _y0 /= _e0sum;
        // _z0 /= _e0sum;

        _x1 /= _e1sum;
        _y1 /= _e1sum;
        // _z1 /= _e1sum;
    }


    std::vector<float> result;
    result.push_back(_x0);
    result.push_back(_y0);
    result.push_back(_e0sum);
    result.push_back(_x1);
    result.push_back(_y1);
    result.push_back(_e1sum);

    return result;

}

float TwoPhotonClustering::FindZEndShower(LCCollection* input_calohits, float _x0[],float _y0[], float r){
/*  For a given cone, Return the last z coordinate associated to it, which has more than 1/2 of maximum hits
    Needs further optimise. */
    
    int n_hits = input_calohits->getNumberOfElements();
    CellIDDecoder<EVENT::CalorimeterHit> cd_RM(input_calohits);

    //  Fill map of calorimeter hits vector per every Z
    std::map<float, std::vector<CalorimeterHit *> > input_hits_map;
    for (int i = 0; i < n_hits; i++)
    {   
        CalorimeterHit *hit = dynamic_cast<CalorimeterHit *>(input_calohits->getElementAt(i));
        float hit_x = hit->getPosition()[0];
        float hit_y = hit->getPosition()[1];
        float hit_z = hit->getPosition()[2];
        int k = cd_RM(hit)["K"]-1;

        // fill the map

        streamlog_out(DEBUG) <<"radius= "<<r<<" "<<sqrt( (hit_x - _x0[k])*(hit_x - _x0[k])+(hit_y - _y0[k])*(hit_y - _y0[k]))<<endl;
        if (sqrt( (hit_x - _x0[k])*(hit_x - _x0[k])+(hit_y - _y0[k])*(hit_y - _y0[k])) > r)  continue;
        auto it = input_hits_map.find(hit_z);
        if (it == input_hits_map.end()) { // did not previously have record in the map
            std::vector<CalorimeterHit *> temp;
            temp.push_back(hit);
            input_hits_map.insert(pair<float, std::vector<CalorimeterHit *>>(hit_z, temp));
            // it->second.push_back(px);
        } else {
            it->second.push_back(hit);
        }
    }
    // sorting the map of hit vectors according to z,
    std::vector<pair<float, std::vector<CalorimeterHit *> > > pairs;
    for (auto &it : input_hits_map)
    {
        pairs.push_back(it);
    }
    sort(pairs.begin(), pairs.end(), [](auto &a, auto &b)
                { return a.first > b.first; }); // Sorting from large to small


    // Find last layer with Nhit > nhits_max/2.
    int max_nhit=0;
    for (auto &pair : pairs) {
        if(pair.second.size()>max_nhit) {
            max_nhit=pair.first;
        };
    }
    float last_layer_z=10000.;
    for (auto &pair : pairs) {
        if(pair.second.size()<max_nhit/4.) {
            last_layer_z=pair.first;
            break;
        };
    }
    // Find last layer with N>4 hits


    return last_layer_z;

}

float TwoPhotonClustering::FindRadius90(LCCollection* input_calohits, float _x0[],float _y0[], float rmean, float energy_cl){
/*  Return a radius in which contains 90% of the cluster's energy (from its centre);
    or the radius till the halfway towards the neighbouring cluster.
    Needs further optimise. */
    int n_hits = input_calohits->getNumberOfElements();

    float step = rmean/100.;
    float radius_90 = rmean; // Maximum if cannot reach 90%
    CellIDDecoder<EVENT::CalorimeterHit> cd_RM(input_calohits);

    //recalculate energy pre-cluster
    float energy_precluster=0;
    for (int i = 0; i < n_hits; i++) {
        // get all simcalorimeter hits
        CalorimeterHit *hit = dynamic_cast<CalorimeterHit *>(input_calohits->getElementAt(i));
        float position[3];
        position[0] = hit->getPosition()[0];
        position[1] = hit->getPosition()[1];
        // position[2] = hit->getPosition()[2];
        float energy = hit->getEnergy();
        int k = cd_RM(hit)["K"]-1;
        if (sqrt( (position[0] - _x0[k])*(position[0] - _x0[k])+(position[1] - _y0[k])*(position[1] - _y0[k]))< rmean)  {
            energy_precluster+=energy;
        }
    }
    for(int istep=1; istep<100; istep++) {

        float energy_sum=0;

        float radius=istep*step;

        // The following loop can be optimised using a sorting in advance
        for (int i = 0; i < n_hits; i++) {
            // get all simcalorimeter hits
            CalorimeterHit *hit = dynamic_cast<CalorimeterHit *>(input_calohits->getElementAt(i));
            float position[3];
            position[0] = hit->getPosition()[0];
            position[1] = hit->getPosition()[1];
            // position[2] = hit->getPosition()[2];
            float energy = hit->getEnergy();
            int k = cd_RM(hit)["K"]-1;

            if (sqrt( (position[0] - _x0[k])*(position[0] - _x0[k])+(position[1] - _y0[k])*(position[1] - _y0[k]))< radius)  {
                energy_sum+=energy;
            }
         }
         if(energy_sum/energy_precluster > ENERGY_FACTOR) {
            radius_90=radius;
            break;
         }
    }

    return radius_90;

}

// define the parametric line equation
void TwoPhotonClustering::line(double t, const double *p, double &x, double &y, double &z) {
     // a parametric line is define from 6 parameters but 4 are independent
     // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
     // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
     x = p[0] + p[1]*t;
     y = p[2] + p[3]*t;
     z = t;
}


std::vector<float> TwoPhotonClustering::WeightedCenter( int n, double x0[], double y0[], double z0[], double e0[]) {

    float x=0, y=0, z=0;
    float wsum=0;
    for(int i=0; i<n; i++) {
            wsum+=e0[i];
            float w=e0[i];
            x+=x0[i]*w;
            y+=y0[i]*w;
            z+=z0[i]*w;
    }
    x/=wsum;
    y/=wsum;
    z/=wsum;

    std::vector<float> result;
    result.push_back(x);
    result.push_back(y);
    result.push_back(z);
    return result;


}

std::vector<TwoPhotonClustering::ClusterInfo> TwoPhotonClustering::GetTwoClustersArray(std::vector<Cluster *> seed_clusters, std::vector<Cluster *> soft_clusters)
{
    /*  Return a struct array of two clusters containing their properties
    with the option of reclustering, where the `soft_clusters` are required.
    */

   /* Does it works for non perpendicular clusters ? */

    double _x[2], _y[2], _z[2], _esum[2];
    double _theta[2], _phi[2]; // TODO: NOT Implemented any change

    for (int i = 0; i < 2; i++)
    {
        // Centre of a cluster is the weighted average of its hits
        _x[i] = seed_clusters.at(i)->getPosition()[0];
        _y[i] = seed_clusters.at(i)->getPosition()[1];
        _z[i] = seed_clusters.at(i)->getPosition()[2];
        _theta[i] = seed_clusters.at(i)->getITheta();
        _phi[i] = seed_clusters.at(i)->getIPhi();
        _esum[i] = seed_clusters.at(i)->getEnergy();
    }
    auto cluster_0 = ROOT::Math::XYPoint(_x[0], _y[0]);
    auto cluster_1 = ROOT::Math::XYPoint(_x[1], _y[1]);

    double rmean = (cluster_0 - cluster_1).R();

    // If _doRecluster == true and we have more than 0 soft-clusters, we recluster them if they are close to the main cluster.
    //  ad-hoc criteria: close is defined as in a radius (x-y) smaller than half the distance between the two main clusters.
    //  shall we tune this? shall we add more criteria like energy cuts ?
    double RECLUSTER_RATIO = 1.0;

    if (soft_clusters.size() > 0 && _doRecluster == true)
    {
        // Weighted positions of the clusters
        _x[0] *= _esum[0];
        _y[0] *= _esum[0];
        _z[0] *= _esum[0];
        _x[1] *= _esum[1];
        _y[1] *= _esum[1];
        _z[1] *= _esum[1];

        for (int isoft = 0; isoft < soft_clusters.size(); isoft++)
        {
            double position[3];
            for (int i = 0; i < 3; i++)
            {
                position[i] = soft_clusters.at(isoft)->getPosition()[i];
            }
            auto soft_position = ROOT::Math::XYPoint(position[0], position[1]);
            double _e = soft_clusters.at(isoft)->getEnergy();

            // Merge the soft cluster into the closer one; else drop it
            if ((cluster_0 - soft_position).R() < rmean * RECLUSTER_RATIO)
            {
                _x[0] += position[0] * _e;
                _y[0] += position[1] * _e;
                _z[0] += position[2] * _e;
                _esum[0] += _e;
            }
            else if ((cluster_1 - soft_position).R() < rmean * RECLUSTER_RATIO)
            {
                _x[1] += position[0] * _e;
                _y[1] += position[1] * _e;
                _z[1] += position[2] * _e;
                _esum[1] += _e;
            }
        }
        _x[0] /= _esum[0];
        _y[0] /= _esum[0];
        _z[0] /= _esum[0];
        _x[1] /= _esum[1];
        _y[1] /= _esum[1];
        _z[1] /= _esum[1];
    }

    std::vector<ClusterInfo> results;
    for (int i = 0; i < 2; i++)
    {
        ClusterInfo result;
        result.position[0] = _x[i];
        result.position[1] = _y[i];
        result.position[2] = _z[i];
        result.energy = _esum[i];
        result.theta = _theta[i];
        result.phi = _phi[i];
        results.push_back(result);
    }
    return results;
}

std::vector<float> TwoPhotonClustering::GetFittedCoordinates(int n, double x[], double y[], double z[], double e[])
{

    std::vector<float> vector_result;
    float energy_cl = 0;
    float e_energy_cl = 0;

    // int n=0;
    // if(iphotons==0) n=n0;
    // else n=n1;

    // energy_cl=two_cluster_coordinates.at(2+3*iphotons);
    // twoPhotonClusters[iphotons].energy = energy_cl;
    // newCluster->setEnergyError(0);
    std::vector<float> center;
    center = WeightedCenter(n, x, y, z, e);

    for (int j = 0; j < 3; j++)
        vector_result.push_back(center.at(j)); // twoPhotonClusters[iphotons].position[j] = center.at(j);

    TGraph2D *g;

    if (_boolAnalogue == false)
    {
        g = new TGraph2D(n, x, y, z);
    }
    else
    {
        double xxx[10000], yyy[10000], zzz[10000];
        int nnn = 0;

        int nn = n;
        for (int ii = 0; ii < nn; ii++)
        {
            for (int jj = 0; jj < int(e[ii]) + 1; jj++)
            {
                xxx[nnn] = x[ii];
                yyy[nnn] = y[ii];
                zzz[nnn] = z[ii];
                nnn++;
            }
        }
        g = new TGraph2D(nnn, xxx, yyy, zzz);
    }
    // #######################################3
    ROOT::Fit::Fitter fitter;

    // make the functor objet
    SumDistance2 sdist(g);
    ROOT::Math::Functor fcn(sdist, 4);
    // set the function and the initial parameter values
    double pStart[4] = {1, 1, 1, 1};
    fitter.SetFCN(fcn, pStart);
    // set step sizes different than default ones (0.3 times parameter values)
    for (int i = 0; i < 4; ++i)
        fitter.Config().ParSettings(i).SetStepSize(0.01);

    bool ok = fitter.FitFCN();
    if (!ok)
    {
        Error("line3Dfit", "Line3D Fit failed");
    }
    const ROOT::Fit::FitResult &result = fitter.Result();
    // result.Print(std::cout);

    // get fit parameters
    const double *parFit = result.GetParams();
    ROOT::Math::Polar3DVector cl_direction_polar;
    float cl_direction_x = 0, cl_direction_y = 0, cl_direction_z = 0;

    // draw the fitted line
    double t0 = -2500;
    double dt = 1000;
    for (int iangle = 0; iangle < 2; ++iangle)
    {
        double t = t0 + dt * iangle;
        double xx, yy, zz;
        line(t, parFit, xx, yy, zz);

        if (iangle == 0 && ok == true)
        {
            cl_direction_x = -xx;
            cl_direction_y = -yy;
            cl_direction_z = -zz;
            _x_res_origin->Fill(-xx);
            _y_res_origin->Fill(30 - yy);
            _r_res_origin->Fill(sqrt(xx * xx + (30. - yy) * (30. - yy)));
        }
        if (iangle == 1)
        {
            cl_direction_x += xx;
            cl_direction_y += yy;
            cl_direction_z += zz;
        }
    }
    cl_direction_polar.SetXYZ(cl_direction_x, cl_direction_y, cl_direction_z);
    vector_result.push_back(cl_direction_polar.Theta()); // twoPhotonClusters[iphotons].theta = cl_direction_polar.Theta();
    vector_result.push_back(cl_direction_polar.Phi());   // twoPhotonClusters[iphotons].phi = cl_direction_polar.Phi();

    // output_cluster->addElement(newCluster);

    // Cleanup TGraph2D *g in TDirectory; TGraph2D object by default has name "Graph2D"
    gROOT->Delete("Graph2D");

    return vector_result;
}