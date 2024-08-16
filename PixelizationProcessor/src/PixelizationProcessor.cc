#include "PixelizationProcessor.hh"

// ----- include for verbosity dependent logging ---------
// #include "marlin/VerbosityLevels.h"
// #include "marlin/StringParameters.h"
// #define SLM streamlog_out(MESSAGE)

// using namespace std;
// using namespace lcio ;
// using namespace marlin ;

PixelizationProcessor aPixelizationProcessor;

PixelizationProcessor::PixelizationProcessor() : Processor("PixelizationProcessor")
{

	// modify processor description
	_description = "";

	// input collections
	registerInputCollection(LCIO::SIMCALORIMETERHIT,
							"SimECALCollectionName",
							"Name of the Sim ECAL Collection",
							_inputColName,
							string("SiEcalCollection"));

	registerInputCollection(LCIO::SIMCALORIMETERHIT,
							"PixelECALCollectionName",
							"Name of the Sim ECAL Collection with realistic pixel and wafer sizes",
							_outputColName,
							string("PixelSiEcalCollection"));
}

PixelizationProcessor::~PixelizationProcessor() {}

void PixelizationProcessor::init()
{

	printParameters();
}

int PixelizationProcessor::GetWafer(float x, float y)
{

	// INPUTS ARE EXPECTED IN um

	// assuming SiWECAL ASUs made of 4 wafers
	// wafer_size = 89700um,
	// space_between_wafers = 100um
	// dead_wafer_space = 305um (in each side)

	// W1 W2 W3 W4 --> (0.0) is between W2 and W7.
	// W5 W6 W7 W8
	//--> (0.0) is between W2 and W7.
	//--> (max,max) is in W4
	//--> (-max,-max) is in W5

	float distance_between_wafers = (space_between_wafers + dead_wafer_space + wafer_size + dead_wafer_space );

	streamlog_out(DEBUG)<<"dist_bet_wafers: "<<distance_between_wafers<<endl;
	streamlog_out(DEBUG)<<"x-G4hit: "<<x<<" y-G4hit:"<<y<<endl;
	// W1,2,3,4
	if (y > 0)
	{
		// W1, W2
		if (x < 0)
		{
		if(fabs(x)>distance_between_wafers) return 1;
		else return 2;
		}
		else
		{ // W3 W4
		if(x>distance_between_wafers) return 4;
		else return 3;
		}
	}
	else
	{ // W5, W6, W7, W8
		// W5, W6
		if (x < 0)
		{
		if(fabs(x)>distance_between_wafers) return 5;
		else return 6;
		}
		else
		{ // W7 W8
		if(x>distance_between_wafers) return 8;
		else return 7;
		}
	}

	return -1;
}

array<float, 2> PixelizationProcessor::GetWaferReferencePadXY(float x, float y, int wafer)
{

	// function to calculate the center position of the most letp-up pixel in the wafer.

	float shift_x = -100;
	float shift_y = -100;
	if (wafer == 1 || wafer == 5)
		shift_x = -2;
	if (wafer == 2 || wafer == 6)
		shift_x = -1;
	if (wafer == 3 || wafer == 7)
		shift_x = 0;
	if (wafer == 4 || wafer == 8)
		shift_x = 1;

	if (wafer < 5)
		shift_y = 0;
	if (wafer > 4)
		shift_y = -1;

	//  WE GET THE center of the REFERENCE PAD for EACH Wafer, which is the one in the left, upper corner
	float first_pad_x = shift_x * (space_between_wafers / 2 + dead_wafer_space + wafer_size + dead_wafer_space + space_between_wafers / 2) + dead_wafer_space + space_between_wafers + pixel_gap / 2 + pixel_size / 2.;
	float first_pad_y = shift_y * (space_between_wafers / 2 + dead_wafer_space + wafer_size + dead_wafer_space + space_between_wafers / 2) + dead_wafer_space + space_between_wafers + pixel_gap / 2 + pixel_size / 2.;

	streamlog_out(DEBUG)<<"GetWaferReferencePadXY: x:"<<x<<" y:"<<y<<" w:"<<wafer<< " fp_x:"<<first_pad_x<<" fp_y:"<<first_pad_y<<endl;

	array<float, 2> result;
	result[0] = first_pad_x;
	result[1] = first_pad_y;

	return result;
}

void PixelizationProcessor::GetPixelHit(float x, float y, float ref_pad_x, float ref_pad_y, int wafer)
{

	// knowing which one is the left-upper pad of the wafer and knowing the pixel-to-pixel distance and numer of pixels per row(column)
	// we check if the hit is inside a pixel or in the gaps or in the dead areas of the wafer

	streamlog_out(DEBUG)<<"GetPixelHit: x:"<<x<<" y:"<<y<<" ref_pad_x:"<<ref_pad_x<<" ref_pad_y:"<<ref_pad_y<<" w:"<<wafer<<endl;

	int waferI = 0;
	if (wafer < 5)
		waferI = (wafer - 1) * npix_row+1;
	else
		waferI = (wafer - 5) * npix_row+1;

	int waferJ = 0;
	if (wafer < 5)
		waferJ = (npix_row + 1);
	else waferJ=1;

	float centers_x[100];
	float centers_y[100];

	for (int i = 0; i < npix_row; i++)
	{
		centers_x[i] = ref_pad_x + i * (pixel_size + pixel_gap);
		centers_y[i] = ref_pad_y + i * (pixel_size + pixel_gap);

	}

	//struct pixel_hit px;
	px.init();

	bool b_hitx = 0;//check if the hit is in an active area of the wafer, x-coordinate
	bool b_hity = 0;
	float hitx = -10000;//world coordinates
	float hity = -10000;
	float hitI = -10000;//IJK coordinates
	float hitJ = -10000;
	for (int i = 0; i < npix_row; i++)
	{

		if (fabs(x - centers_x[i]) < pixel_size / 2.)
		{
			b_hitx = 1;
			hitx = centers_x[i];
			hitI = i;
		}
		if (fabs(y - centers_y[i]) < pixel_size / 2.)
		{
			b_hity = 1;
			hity = centers_y[i];
			hitJ = i;
		}
	}

	if (b_hitx == 1 && b_hity == 1)
	{
		px.hit = 1;
		px.x = hitx/1000.;//in mm!!!
		px.y = hity/1000.;
		px.I = hitI + waferI;
		px.J = hitJ + waferJ;
		//notice that I,J start from 1, not from 0
	}

	//return result;
}

void PixelizationProcessor::processRunHeader(LCRunHeader *run)
{
}

void PixelizationProcessor::processEvent(LCEvent *evt)
{

	try
	{

		LCCollection *input_ecal = evt->getCollection(_inputColName);
		int n_el = input_ecal->getNumberOfElements();
		CellIDDecoder<EVENT::SimCalorimeterHit> cd(input_ecal);

		IMPL::LCCollectionVec *output_ecal = new IMPL::LCCollectionVec(LCIO::SIMCALORIMETERHIT);
	
		IMPL::LCFlagImpl flag ;
    	flag.setBit( EVENT::LCIO::CHBIT_LONG ) ;
    	flag.setBit( EVENT::LCIO::CHBIT_ID1 ) ;
		flag.setBit( EVENT::LCIO::CHBIT_STEP ) ;
		output_ecal->setFlag(flag.getFlag());

      	CellIDEncoder<IMPL::SimCalorimeterHitImpl> cd_output("I:8,J:8,K:8,wafer:4",output_ecal);

		for (int i = 0; i < n_el; i++)
		{

			// since we simulate one large silicon structure simulating the full layer, there is only one hit per layer.
			// we want to look at the sub-hits (MCParticles) to study in detail the position of each hit
			// all hits hitting in the dead areas are neglected for the calculation of the energy per pixel

			// For the moment we are ignoring the time information !!

			SimCalorimeterHit *ecalhit = dynamic_cast<SimCalorimeterHit *>(input_ecal->getElementAt(i));

			int k = cd(ecalhit)["layer"]+1;

			int n_cont = ecalhit->getNMCParticles();

			map<int, vector<pixel_hit>> map_hits_per_layer;

			bool save = 0;
			for (int j = 0; j < n_cont; j++)
			{

				float x = ecalhit->getStepPosition(j)[0]*1000.;//x-position of the hit in um
				float y = ecalhit->getStepPosition(j)[1]*1000.;//y

				int wafer = GetWafer(x, y);	
				streamlog_out(DEBUG)<<" Wafer: "<<wafer<<endl;	  // I calculate which wafer is hit

				array<float, 2> reference_pad = GetWaferReferencePadXY(x, y, wafer); // I calculate the reference pad of that wafer

				GetPixelHit(x, y, reference_pad[0], reference_pad[1],wafer);//here we find the pixel in which the hit falls
				px.write();

				//extra information about the hit in the pixel
				if (px.hit == 1)
				{
					px.e = ecalhit->getEnergyCont(j);
					px.K = k;
					px.z = ecalhit->getPosition()[2];
					px.t = (float)ecalhit->getTimeCont(j);// not yet used for anything...
					px.wafer=wafer;
					px.step_x=ecalhit->getStepPosition(j)[0];
					px.step_y=ecalhit->getStepPosition(j)[1];
					px.step_z=ecalhit->getStepPosition(j)[2];
					px.pdg=ecalhit->getPDGCont(j);
					px.length=ecalhit->getLengthCont(j);

					bool store=px.goodhit();

					// I build the map of hits in this layer, I want to separate them in different I,J,K.
					if(store==true) {
						int counter=1000000*px.K + 10000*px.J + px.I;
						auto it = map_hits_per_layer.find(counter);
						if( it== map_hits_per_layer.end()) {
							vector<pixel_hit> temp;
							temp.push_back(px);
							map_hits_per_layer.insert(pair<int,vector<pixel_hit>>(counter,temp));
							//it->second.push_back(px);
						} else {
							it->second.push_back(px);
						}
					}
				}
			}

			//loop over the map of hits per readout pchannel (pixel). We will dump these hits into SimCalorimeterHits. 
			//We store only one hit per channel/pixel and we keep all the MC particles (G4) associated to that hit.
			//the ones not falling inside a pixel ae discarded.
			
			for (auto it = map_hits_per_layer.begin(); it!= map_hits_per_layer.end(); it++) {
				bool addelement=0;
				SimCalorimeterHitImpl *new_ecalhit = new SimCalorimeterHitImpl(); // new object

				float position[3]={it->second.at(0).x,it->second.at(0).y,it->second.at(0).z};

				addelement=1;
				for(int i=0; i<it->second.size(); i++) {
					float position_step[3]={it->second.at(i).step_x,it->second.at(i).step_y,it->second.at(i).step_z};

					MCParticleImpl *pixelhit = new MCParticleImpl(); 
	
					new_ecalhit->addMCParticleContribution(pixelhit, (float)it->second.at(i).e, (float)it->second.at(i).t,(float)it->second.at(i).length,(float)it->second.at(i).pdg,position_step);
					cd_output["I"]=it->second.at(i).I;
					cd_output["J"]=it->second.at(i).J;
					cd_output["K"]=it->second.at(i).K;
					cd_output["wafer"]=it->second.at(i).wafer;
					cd_output.setCellID(new_ecalhit);
				}
				new_ecalhit->setPosition(position);//we store the X,Y,Z position of only one hit, since the granularity is exactly the pad-size
				
				if (addelement == 1)
				{
					output_ecal->addElement(new_ecalhit);
				}
			}
		}
		evt->addCollection(output_ecal, _outputColName);
	}catch (DataNotAvailableException &e){
		streamlog_out(DEBUG) << "Whoops!....\n";
		streamlog_out(DEBUG) << e.what();
	}
}

void PixelizationProcessor::check(LCEvent *evt)
{
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}

void PixelizationProcessor::end()
{
}
