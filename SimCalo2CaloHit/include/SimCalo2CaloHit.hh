#ifndef SimCalo2CaloHit_h
#define SimCalo2CaloHit_h 1
#include <iomanip>
#include "marlin/Processor.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/SimCalorimeterHitImpl.h>
#include <EVENT/CalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
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

class SimCalo2CaloHit : public Processor
{

public:
  virtual Processor *newProcessor() { return new SimCalo2CaloHit; }

  SimCalo2CaloHit();
  virtual ~SimCalo2CaloHit();

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

private:

  string _inputColName;
  string _outputColName;

  vector<float> _gev2mip;
  float _adc_s_over_n;
  float _trig_s_over_n;
  float _mip_threshold;
  bool _digitization;
  bool _convert_to_digital;


};

#endif
