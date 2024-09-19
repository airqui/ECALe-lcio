#ifndef TwoPhotonClustering_h
#define TwoPhotonClustering_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <EVENT/Cluster.h>
#include <EVENT/CalorimeterHit.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/ClusterImpl.h>

#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>

#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <Math/Cartesian3D.h>
#include <Math/RotationZYX.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TRotation.h>

using namespace lcio ;
using namespace marlin ;
//using VectorXYZ   = ROOT::Math::XYZVectorF ;

// Processor to perform a two photon reclustering, starting from single layer clusers (NNClusters)

class TwoPhotonClustering : public Processor {
public:
  
    virtual Processor*  newProcessor() { return new TwoPhotonClustering ; }
  
  
    TwoPhotonClustering() ;
  
    /** Called at the begin of the job before anything is read.
        * Use to initialize the processor, e.g. book histograms.
    */
    virtual void init() ;
  
    /** Called for every run.
    */
    virtual void processRunHeader( LCRunHeader* run ) ;
  
    /** Called for every event - the working horse.
    */
    virtual void processEvent( LCEvent * evt ) ; 
  
  
    virtual void check( LCEvent * evt ) ; 
  
  
    /** Called after data processing for clean up.
    */
    virtual void end() ;

    // define the parametric line equation
    
private:

    virtual std::vector<float> GetCenterCoordinates( std::vector<Cluster*> cl1, std::vector<Cluster*> cl2);
    virtual float FindRadius90(LCCollection* input_calohits, float x0,float y0, float rmean, float energy_cl);
    virtual float FindZEndShower(LCCollection* input_calohits, float _x0,float _y0, float r);
    virtual void line(double t, const double *p, double &x, double &y, double &z) ;
    virtual std::vector<float> WeightedCenter( int n, double x0[], double y0[], double z0[], double e0[]);


// function Object to be minimized
    struct SumDistance2 {
    public:
        // the TGraph is a data member of the object
        TGraph2D *fGraph;

        SumDistance2(TGraph2D *g) : fGraph(g) {};

        // calculate distance line-point
        double distance2(double x,double y,double z, const double *p) {
            // distance line point is D= | (xp-x0) cross  ux |
            // where ux is direction of line and x0 is a point in the line (like t = 0)
            ROOT::Math::XYZVector xp(x,y,z);
            ROOT::Math::XYZVector x0(p[0], p[2], 0. );
            ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
            ROOT::Math::XYZVector u = (x1-x0).Unit();
            double d2 = ((xp-x0).Cross(u)).Mag2();
            return d2;
        };

        // implementation of the function to be minimized
        double operator() (const double *par) {
            assert(fGraph != 0);
            double * x = fGraph->GetX();
            double * y = fGraph->GetY();
            double * z = fGraph->GetZ();
            int npoints = fGraph->GetN();
            double sum = 0;
            for (int i  = 0; i < npoints; ++i) {
                double d = distance2(x[i],y[i],z[i],par);
                sum += d;
            }

            return sum;
        }
    };



    struct ClusterInfo {
    public:
        float position[3];
        double energy;
        double phi, theta;
        double error_position[3];   // TODO: NO implementation yet
        double error_energy;        // TODO: NO implementation yet
        double *error_direction;    // TODO: NO implementation yet
    };
    virtual std::vector<ClusterInfo> GetTwoClustersArray(std::vector<Cluster *> seed, std::vector<Cluster *> soft);
 
 protected:

    /** Input collection name.
    */

    std::string _CaloHitCol;
    std::string _NNClusterCol;
    std::string _outputColName;

    //int _nPhotonstoReconstruct;
    int _strategytofollow;
    bool _doRecluster;
    float ENERGY_FACTOR;
    float DISTANCE_RATIO;

    //float _distCut{};
    //float _eCut{};

    //int _nThetaPhi{};

    int _nRun{};
    int _nEvt{};
  

    //   NNClusteringer* _clusterer ;

} ;

#endif