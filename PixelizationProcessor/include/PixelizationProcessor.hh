#ifndef PixelizationProcessor_h
#define PixelizationProcessor_h 1
#include <iomanip>
#include "marlin/Processor.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/SimCalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <cmath>
#include <string>
#include <vector>

#include <TF1.h>
#include <TRandom.h>
#include <TFile.h>
#include <TTree.h>
// #include <TLorentzVector.h>
using namespace lcio;
using namespace marlin;
using namespace std;

class PixelizationProcessor : public Processor
{

public:
  virtual Processor *newProcessor() { return new PixelizationProcessor; }

  PixelizationProcessor();
  virtual ~PixelizationProcessor();

  /** Called at the begin of the job before anything is read.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader *run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent *evt);

  virtual void check(LCEvent *evt);

  virtual void end();

  // structure to store hits from G4
  struct pixel_hit{
    public:
      bool hit; //is hit or not
      float e;//energy
      float t;
      float x;//pad-coordinates in world G4 units, mm.
      float y;
      float z;
      int I;//coordinates I, J , K, with I=column, J=row, K= layer
      int J;
      int K;
      int wafer;
      float length;
      int pdg;
      float step_x;//coordinates of G4 step contribution
      float step_y;
      float step_z;
      void init() {
        hit=0; //is hit or not
        e=-1;//energy
        t=-1;
        x=-10000;//coordinates in G4 units, mm.
        y=-10000;
        z=-10000;
        I=-10000;//coordinates I, J , K, with I=column, J=row, K= layer
        J=-10000;
        K=-10000;
        wafer=-1;
        pdg=0;
        length=-1;
        step_x=0; step_y=0; step_z=0;

      }
      bool goodhit() {
        if(hit==1 && e>0 && x>-10000 
          && y>-10000
          && z>-10000
          && I>-10000
          && J>-10000
          && K>-10000
          && wafer>-1) return true;
        else return false;
      }
      void write(){
        streamlog_out(DEBUG)<<"pixel_hit: "<< hit<<", E "<<e<<" T "<<t<<", xyz "<<x<<" "<<y<<" "<<z<<", IJK "<<I<<" "<<J<<" "<<K<<" wafer:"<<wafer<<endl;
      }
    
  };
  virtual int GetWafer(float , float);
  virtual array<float, 2> GetWaferReferencePadXY(float , float, int);
  void GetPixelHit(float, float, float, float, int);

private:

  string _inputColName;
  string _outputColName;

  // WAFER INFO --> this could be given as input to the processor.
  // safety margin between two wafers
  float space_between_wafers = 100;//um
  // wafer size. All capitalised variables from HMMS rating doc K30-B70125
  float CHIP_SIZE = 89700;//um
  float ACTIVE_AREA = 88480;//um  -- only active area!
  float dead_wafer_space = 1/2.*(CHIP_SIZE - ACTIVE_AREA);// (in each side)
  // pixel size
  float PIXEL_PITCH = 5530;//um
  float PIXEL_GAP = 10;//um
  float pixel_size = PIXEL_PITCH - PIXEL_GAP;
  int npix_col = 16;//number of pixels per column of a wafer
  int npix_row = 16;//number of pixels per row of a wafer

  //we are assuming an 8 sensor plane detector (two ASU's)
  //seen from the beam point of view, the particle will see:
  // W1 W2 W3 W4 
  // W5 W6 W7 W8
  //--> (0.0) is between W2 and W7.
  //--> (max,max) is in W4
  //--> (-max,-max) is in W5

  struct pixel_hit px;



};

#endif
