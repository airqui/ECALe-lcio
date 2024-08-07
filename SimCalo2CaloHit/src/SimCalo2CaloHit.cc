#include "SimCalo2CaloHit.hh"


#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

// ----- include for verbosity dependent logging ---------
// #include "marlin/VerbosityLevels.h"
// #include "marlin/StringParameters.h"
// #define SLM streamlog_out(MESSAGE)

// using namespace std;
// using namespace lcio ;
// using namespace marlin ;

SimCalo2CaloHit aSimCalo2CaloHit;

SimCalo2CaloHit::SimCalo2CaloHit() : Processor("SimCalo2CaloHit")
{

	// modify processor description
	_description = "";

	// input collections
	registerInputCollection(LCIO::SIMCALORIMETERHIT,
							"PixelECALCollectionName",
							"Name of the Sim ECAL Collection with realistic pixel and wafer sizes",
							_inputColName,
							string("PixelSiEcalCollection"));
	//output
	registerInputCollection(LCIO::CALORIMETERHIT,
							"ECALCollectionName",
							"Name of the ECAL Collection with realistic pixel and wafer sizes, MIP calibrated, preliminary digitized",
							_outputColName,
							string("ECAL"));

	vector<float> gev2mip = {0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014,0.014};//dummy values
	registerProcessorParameter( "GeV2MIPfactor",
			"correction factor from GeV2MIP",
			_gev2mip,
			gev2mip);
	registerProcessorParameter( "ADC_Signal_over_Noise",
			"ADC signal over noise value for a dummy digitization, in MIP units",
			_adc_s_over_n,
			(float)10.);
	registerProcessorParameter( "TRIG_Signal_over_Noise",
			"Trigger signal over noise value for a dummy digitization, in MIP units",
			_trig_s_over_n,
			(float)5.);
	registerProcessorParameter( "MIP_threshold",
			"threshold for trigger on MIP units",
			_mip_threshold,
			(float)0.5);
	registerProcessorParameter( "DoDigitization",
			"true= if we want to perform digitization",
			_digitization,
			true);
	
}

SimCalo2CaloHit::~SimCalo2CaloHit() {}

void SimCalo2CaloHit::init()
{

	printParameters();
}

void SimCalo2CaloHit::processRunHeader(LCRunHeader *run)
{
}

void SimCalo2CaloHit::processEvent(LCEvent *evt)
{

	try
	{

		LCCollection *input_ecal = evt->getCollection(_inputColName);
		int n_el = input_ecal->getNumberOfElements();
		CellIDDecoder<EVENT::SimCalorimeterHit> cd(input_ecal);

		IMPL::LCCollectionVec *output_ecal = new IMPL::LCCollectionVec(LCIO::CALORIMETERHIT);
	
		IMPL::LCFlagImpl flag ;
    	flag.setBit( EVENT::LCIO::CHBIT_LONG ) ;
    	flag.setBit( EVENT::LCIO::CHBIT_ID1 ) ;
		flag.setBit( EVENT::LCIO::CHBIT_STEP ) ;
		flag.setBit( EVENT::LCIO::RCHBIT_ENERGY_ERROR) ;
		output_ecal->setFlag(flag.getFlag());

      	CellIDEncoder<IMPL::CalorimeterHitImpl> cd_output("I:9,J:9,K:6,M:3,S-1:3,K-1:6",output_ecal);

		for (int i = 0; i < n_el; i++)
		{

			// get all simcalorimeter hits
			SimCalorimeterHit *ecalhit = dynamic_cast<SimCalorimeterHit *>(input_ecal->getElementAt(i));

			int cd_i = cd(ecalhit)["I"];
			int cd_j = cd(ecalhit)["J"];
			int cd_k = cd(ecalhit)["K"];

			float x = ecalhit->getPosition()[0];//x-position of the hit in mm
			float y = ecalhit->getPosition()[1];//y
			float z = ecalhit->getPosition()[2];//y
			float e=0;
		
			int error= sizeof(_gev2mip); 
			if(error < (cd_k-2) ) {
			streamlog_out(ERROR) << "MIP CALIBRATION FACTOR NOT EXISTING for this layer! K="
				<<cd_k<<" size mip value array="<<sizeof(_gev2mip)/sizeof(_gev2mip[0])<<endl;
			}
			else e=ecalhit->getEnergy()/_gev2mip[cd_k-1];

			float t= 0;//ecalhit->getTime();


			float e_error = fabs(CLHEP::RandGauss::shoot( 1.0, 1./_adc_s_over_n));//error in energy is the only the pedestal width, in mip units
			float e_digitized = CLHEP::RandGauss::shoot( e, 1./_adc_s_over_n);

			float e_digitized_trig = CLHEP::RandGauss::shoot( e, 1./_trig_s_over_n);
			if(e_digitized_trig>_mip_threshold) {
				CalorimeterHitImpl *new_ecalhit = new CalorimeterHitImpl(); // new object

				float position[3]={x,y,z};
				new_ecalhit->setPosition(position);
				new_ecalhit->setEnergy(e_digitized);
				new_ecalhit->setEnergyError(e_error);
				new_ecalhit->setTime(t);
				cd_output["I"]=cd_i;
				cd_output["J"]=cd_j;
				cd_output["K"]=cd_k;
				cd_output["K-1"]=cd_k-1;
				cd_output["S-1"]=0;
				cd_output["M"]=1;

				cd_output.setCellID(new_ecalhit);
				output_ecal->addElement(new_ecalhit);
			}
		}
		evt->addCollection(output_ecal, _outputColName);

	}catch (DataNotAvailableException &e){
		streamlog_out(DEBUG) << "Whoops!....\n";
		streamlog_out(DEBUG) << e.what();
	}
}

void SimCalo2CaloHit::check(LCEvent *evt)
{
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}

void SimCalo2CaloHit::end()
{
}
