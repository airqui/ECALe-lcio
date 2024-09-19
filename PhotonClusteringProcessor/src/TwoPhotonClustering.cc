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
                               (float) 0.9);
}


void TwoPhotonClustering::init() {
    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;
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


            
            
        // START WITH THE FIT TESTS
        // FOR THE MOMENT THIS IS DONE ASSUMING 2 seed CLUSTERS
        const int n = n_hits; // Number of total hits
        int n0 = 0;
        int n1 = 0;
        /* Change the following lines to double because of WeightedPoints3D,
           hopefully will not create a fuss... */
        double x0[n], y0[n], z0[n], e0[n]; // all hits near cluster seed 0
        double x1[n], y1[n], z1[n], e1[n]; // all hits near cluster seed 1

        std::vector<float> two_cluster_coordinates=GetCenterCoordinates(seed_clusters,soft_clusters);
        // float _x0=two_cluster_coordinates.at(0);
        // float _y0=two_cluster_coordinates.at(1);
        // float _e0=two_cluster_coordinates.at(2);
        // float _x1=two_cluster_coordinates.at(3);
        // float _y1=two_cluster_coordinates.at(4);
        // float _e1=two_cluster_coordinates.at(5);

        ClusterInfo* twoPhotonClusters = GetTwoClustersArray(seed_clusters, soft_clusters);
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
            float r0_90=FindRadius90(input_calohits,_x0,_y0,rmean,_e0);
            float r1_90=FindRadius90(input_calohits,_x1,_y1,rmean,_e1);
            float last_layer_r0_90=FindZEndShower(input_calohits,_x0,_y0,r0_90); // last layer with 1/2 of the maximum hits
            float last_layer_r1_90=FindZEndShower(input_calohits,_x1,_y1,r1_90);

            /* **RECLUSTER THE HITS** */
            for (int i = 0; i < n_hits; i++) {
                // get all simcalorimeter hits
                CalorimeterHit *hit = dynamic_cast<CalorimeterHit *>(input_calohits->getElementAt(i));
                float position[3];
                position[0] = hit->getPosition()[0];
                position[1] = hit->getPosition()[1];
                position[2] = hit->getPosition()[2];

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
                } else

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
                } else

                // ********************** strategy 3, we use the cylinder that contains the 90% of the energy
                if (_strategytofollow == 3)
                {
                    if (sqrt( (position[0] - _x0)*(position[0] - _x0)+(position[1] - _y0)*(position[1] - _y0))< r0_90)  
                    {
                        x0[n0] = position[0];
                        y0[n0] = position[1];
                        z0[n0] = position[2];
                        e0[n0] = hit->getEnergy();
                        n0++;
                    }
                    if( sqrt( (position[0] - _x1)*(position[0] - _x1)+(position[1] - _y1)*(position[1] - _y1) )<r1_90)
                    {
                        x1[n1] = position[0];
                        y1[n1] = position[1];
                        z1[n1] = position[2];
                        e1[n1] = hit->getEnergy();
                        n1++;
                    }
                } else

            // ********************** strategy 4, we use the cylinder that contains the 90% of the energy and we remove the last layers with only few hits
                if (_strategytofollow == 4)
                {
                    if (sqrt( (position[0] - _x0)*(position[0] - _x0)+(position[1] - _y0)*(position[1] - _y0))< r0_90 && position[2]<last_layer_r0_90)  
                    {
                        x0[n0] = position[0];
                        y0[n0] = position[1];
                        z0[n0] = position[2];
                        e0[n0] = hit->getEnergy();
                        n0++;
                    }
                    if( sqrt( (position[0] - _x1)*(position[0] - _x1)+(position[1] - _y1)*(position[1] - _y1) )<r1_90 && position[2]<last_layer_r1_90)
                    {
                        x1[n1] = position[0];
                        y1[n1] = position[1];
                        z1[n1] = position[2];
                        e1[n1] = hit->getEnergy();
                        n1++;
                    }
                }
            }

            bool flagUsingFit = false;
            if (flagUsingFit) {
                for (int iphotons = 0; iphotons < 2; iphotons++) {

                    ClusterImpl *newCluster = new ClusterImpl(); // new object
                    float energy_cl=0;
                    float e_energy_cl=0;
                

                    // int n=0;
                    // if(iphotons==0) n=n0;
                    // else n=n1;
                    
                    energy_cl=two_cluster_coordinates.at(2+3*iphotons);
                    newCluster->setEnergy(energy_cl);
                    newCluster->setEnergyError(0);
                    std::vector<float> center;
                    if(iphotons==0) center=WeightedCenter(n0,x0,y0,z0,e0);
                    else center=WeightedCenter(n1,x1,y1,z1,e1);

                    float position[3]={center.at(0),center.at(1),center.at(2)};
                    // float eposition[3]={center.at(3),center.at(4),center.at(5)}; // NOT DONE!
                    newCluster->setPosition(position);
                    //cluster->setPositionError(eposition);

                    TGraph2D *g;
                    if(iphotons==0) g = new TGraph2D(n0, x0, y0, z0);
                    else g = new TGraph2D(n1, x1, y1, z1);
                    //#######################################3
                    ROOT::Fit::Fitter  fitter;

                    // make the functor objet
                    SumDistance2 sdist(g);
                    ROOT::Math::Functor fcn(sdist,4);
                    // set the function and the initial parameter values
                    double pStart[4] = {1,1,1,1};
                    fitter.SetFCN(fcn,pStart);
                    // set step sizes different than default ones (0.3 times parameter values)
                    for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);

                    bool ok = fitter.FitFCN();
                    if (!ok) {
                        Error("line3Dfit","Line3D Fit failed");
                    }
                    const ROOT::Fit::FitResult & result = fitter.Result();
                    result.Print(std::cout);

                    // get fit parameters
                    const double * parFit = result.GetParams();
                    ROOT::Math::Polar3DVector cl_direction_polar;
                    float cl_direction_x=0,cl_direction_y=0, cl_direction_z=0;

                    // draw the fitted line
                    double t0 = 0;
                    double dt = 1000;
                    for (int iangle = 0; iangle <2;++iangle) {
                        double t = t0+ dt*iangle;
                        double xx,yy,zz;
                        line(t,parFit,xx,yy,zz);
                
                        
                        if(iangle==0) {cl_direction_x=-xx;cl_direction_y=-yy;cl_direction_z=-zz;}
                        if(iangle==1) {cl_direction_x+=xx;cl_direction_y+=yy;cl_direction_z+=zz;}
                    }
                    cl_direction_polar.SetXYZ(cl_direction_x, cl_direction_y, cl_direction_z);
                    float cl_theta = cl_direction_polar.Theta();
                    float cl_phi = cl_direction_polar.Phi();
                    newCluster->setITheta(cl_theta);
                    newCluster->setIPhi(cl_phi);
                
                    output_cluster->addElement(newCluster);

                    // Cleanup TGraph2D *g in TDirectory; TGraph2D object by default has name "Graph2D"
                    gROOT->Delete("Graph2D");
                }
            } else {
                WeightedPoints3D *cl_hits;

                for (int iphotons = 0; iphotons < 2; iphotons++) {
                    if (iphotons==0) {
                        cl_hits = new WeightedPoints3D(n0, &e0[0], &x0[0], &y0[0], &z0[0]);
                    } else if (iphotons==1) {
                        cl_hits = new WeightedPoints3D(n1, &e1[0], &x1[0], &y1[0], &z1[0]);
                    } else streamlog_out(ERROR) << "WHAT??? Chech index `iphotons`" <<endl;
                    
                    twoPhotonClusters[iphotons].energy = cl_hits->getTotalWeight();
                    
                    double* cog = cl_hits->getCentreOfGravity();
                    for (int j=0; j<3; j++) twoPhotonClusters[iphotons].position[j] = cog[j];

                    double* evpv = cl_hits->getEigenVecPolar();
                    double evp[2][3]; for( int iii=0 ; iii<2 ; iii++ ) { for ( int jjj=0 ; jjj<3 ; jjj++ ) { evp[iii][jjj]=evpv[jjj+iii*3] ; } } ;
                    
                    twoPhotonClusters[iphotons].theta = evp[0][2];
                    twoPhotonClusters[iphotons].phi = evp[1][2];
                }
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

float TwoPhotonClustering::FindZEndShower(LCCollection* input_calohits, float _x0,float _y0, float r){
/*  For a given cone, Return the last z coordinate associated to it, which has more than 1/2 of maximum hits
    Needs further optimise. */
    
    int n_hits = input_calohits->getNumberOfElements();

    //  Fill map of calorimeter hits vector per every Z
    std::map<float, std::vector<CalorimeterHit *> > input_hits_map;
    for (int i = 0; i < n_hits; i++)
    {   
        CalorimeterHit *hit = dynamic_cast<CalorimeterHit *>(input_calohits->getElementAt(i));
        float hit_x = hit->getPosition()[0];
        float hit_y = hit->getPosition()[1];
        float hit_z = hit->getPosition()[2];

        // fill the map

        streamlog_out(DEBUG) <<"radius= "<<r<<" "<<sqrt( (hit_x - _x0)*(hit_x - _x0)+(hit_y - _y0)*(hit_y - _y0))<<endl;
        if (sqrt( (hit_x - _x0)*(hit_x - _x0)+(hit_y - _y0)*(hit_y - _y0)) > r)  continue;
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
        if(pair.second.size()>max_nhit/2.) {
            last_layer_z=pair.first;
            break;
        };
    }
    // Find last layer with N>4 hits


    return last_layer_z;

}

float TwoPhotonClustering::FindRadius90(LCCollection* input_calohits, float _x0,float _y0, float rmean, float energy_cl){
/*  Return a radius in which contains 90% of the cluster's energy (from its centre);
    or the radius till the halfway towards the neighbouring cluster.
    Needs further optimise. */
    int n_hits = input_calohits->getNumberOfElements();

    float step = rmean/100.;
    float radius_90 = rmean; // Maximum if cannot reach 90%

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

            if (sqrt( (position[0] - _x0)*(position[0] - _x0)+(position[1] - _y0)*(position[1] - _y0))< radius)  {
                energy_sum+=energy;
            }
         }
         if(energy_sum/energy_cl > ENERGY_FACTOR) {
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

TwoPhotonClustering::ClusterInfo* TwoPhotonClustering::GetTwoClustersArray(std::vector<Cluster*> seed_clusters, std::vector<Cluster*> soft_clusters) {
    /*  Return a struct array of two clusters containing their properties
    with the option of reclustering, where the `soft_clusters` are required.
    */

    double _x[2],_y[2],_z[2],_esum[2];
    double _theta[2],_phi[2]; // TODO: NOT Implemented any change
    
    for (int i=0; i<2; i++) {
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
    
    double rmean = (cluster_0-cluster_1).R();

    //If _doRecluster == true and we have more than 0 soft-clusters, we recluster them if they are close to the main cluster.
    // ad-hoc criteria: close is defined as in a radius (x-y) smaller than half the distance between the two main clusters.
    // shall we tune this? shall we add more criteria like energy cuts ?
    double RECLUSTER_RATIO = 1.0;

    if(soft_clusters.size()>0 && _doRecluster==true) {
        // Weighted positions of the clusters
        _x[0] *= _esum[0]; _y[0] *= _esum[0]; _z[0] *= _esum[0];
        _x[1] *= _esum[1]; _y[1] *= _esum[1]; _z[1] *= _esum[1];
        
        for(int isoft=0; isoft<soft_clusters.size(); isoft++) {
            double position[3];
            for (int i=0; i<3; i++) {
                position[i] = soft_clusters.at(isoft)->getPosition()[i];
            }
            auto soft_position = ROOT::Math::XYPoint(position[0], position[1]);
            double _e = soft_clusters.at(isoft)->getEnergy();

            // Merge the soft cluster into the closer one; else drop it
            if ((cluster_0 - soft_position).R() < rmean*RECLUSTER_RATIO) {
                _x[0] += position[0] * _e;
                _y[0] += position[1] * _e;
                _z[0] += position[2] * _e;
                _esum[0] += _e;
            } else if ((cluster_1 - soft_position).R() < rmean*RECLUSTER_RATIO) {
                _x[1] += position[0] * _e;
                _y[1] += position[1] * _e;
                _z[1] += position[2] * _e;
                _esum[1] += _e;
            }
            
        }
        _x[0] /= _esum[0]; _y[0] /= _esum[0]; _z[0] /= _esum[0];
        _x[1] /= _esum[1]; _y[1] /= _esum[1]; _z[1] /= _esum[1];
    }
    
    ClusterInfo result[2];
    for (int i=0; i<2; i++) {
        result[i].position[0] = _x[i];
        result[i].position[1] = _y[i];
        result[i].position[2] = _z[i];
        result[i].energy = _esum[i];
        result[i].theta = _theta[i];
        result[i].phi = _phi[i];
    }
    return result;
}