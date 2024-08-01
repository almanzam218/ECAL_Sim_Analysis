#include "ExampleProcessor.hh"

// ROOT

// ----- include for verbosity dependent logging ---------
// #include "marlin/VerbosityLevels.h"
// #include "marlin/StringParameters.h"
// #define SLM streamlog_out(MESSAGE)

// using namespace std;
using namespace lcio ;
using namespace marlin ;
using EVENT::LCCollection;
using EVENT::MCParticle;
using EVENT::ReconstructedParticle;
using EVENT::Track;
using EVENT::Vertex;
using IMPL::LCRelationImpl;
using IMPL::ReconstructedParticleImpl;
using IMPL::TrackImpl;
using IMPL::TrackStateImpl;
using std::string;
using std::vector;
using UTIL::LCRelationNavigator;

ExampleProcessor aExampleProcessor;

ExampleProcessor::ExampleProcessor() : Processor("ExampleProcessor")
{

	// modify processor description
	_description = "";

	// input collections
	registerInputCollection(LCIO::MCPARTICLE,"MCCollectionName",
				"Name of the MC collection",
							_MCColName,
							std::string("MCParticles"));
	registerInputCollection(LCIO::SIMCALORIMETERHIT,
							"ECALCollection",
							"Name of the Sim ECAL Collection",
							_ECALColName,
							std::string("SiEcalCollection"));
}

ExampleProcessor::~ExampleProcessor() {}

void ExampleProcessor::init()
{
	int testVariable = 1;
	int test2 = 10;
	printParameters();
	_xHist = new TH1D("_xHist","X Distribution",30, -0.5, 29.5);
	_yHist = new TH1D("_yHist","Y Distribution",30, -0.5, 29.5);
	_zHist = new TH1D("_zHist","Z Distribution",30, -0.5, 29.5);
	_cellEnergyHist = new TH1F("_cellEnergyHist","Energy deposited in cells Distribution",100, 0, 0.01);
	_evEnergyHist = new TH1F("_evEnergyHist","Energy of shower Distribution",100, 0, 10);
	
	_xyHist = new TH2D("_xyHist","XY view all events",30,-0.5,29.5,30, -0.5, 29.5);

}


void ExampleProcessor::ShowMCInfo(EVENT::LCCollection *myCollection)
{
  int number = myCollection->getNumberOfElements();
  
    for (int i = 0; i < number; i++)
    {

      MCParticle *particle = dynamic_cast<MCParticle *>(myCollection->getElementAt(i));
      vector<MCParticle *> daughters = particle->getDaughters();

      streamlog_out(MESSAGE) << "\n MCCollection, particle:"  ;
      streamlog_out(MESSAGE) << " pdg=" << particle->getPDG();
      streamlog_out(MESSAGE) << " satus=" << particle->getGeneratorStatus();
      streamlog_out(MESSAGE) << " Ndaughters=" << daughters.size();
      streamlog_out(MESSAGE) << " E=" << particle->getEnergy();
      streamlog_out(MESSAGE) << " px=" << particle->getMomentum()[0];
      streamlog_out(MESSAGE) << " py=" << particle->getMomentum()[1];
      streamlog_out(MESSAGE) << " pz=" << particle->getMomentum()[2];
      streamlog_out(MESSAGE) << " m=" << particle->getMass();
      streamlog_out(MESSAGE) << " charge=" << particle->getCharge();

    }
  	//streamlog_out(MESSAGE) << "this treeed" << std::endl;

}

 void ExampleProcessor::ShowECALInfo(EVENT::LCCollection *myCollection)
{
  int number = myCollection->getNumberOfElements();
  CellIDDecoder<EVENT::SimCalorimeterHit> cd(myCollection);
	double totalEnergy = 0;

  for (int i = 0; i < number; i++)
    {

      SimCalorimeterHit *ecalhit = dynamic_cast<SimCalorimeterHit *>(myCollection->getElementAt(i));
      int x_in_IJK_coordinates = cd(ecalhit)["x"];
      int y_in_IJK_coordinates = cd(ecalhit)["y"];
      int z_in_IJK_coordinates = cd(ecalhit)["layer"];
		//_coordinateZ=z_in_IJK_coordinates;
      streamlog_out(MESSAGE) << "\n SimCalorimeterHit, :" << i;
      streamlog_out(MESSAGE) << " cellID-encoded=" << ecalhit->getCellID0();
      streamlog_out(MESSAGE) << " x_in_IJK_coordinates=" << x_in_IJK_coordinates;
      streamlog_out(MESSAGE) << " y_in_IJK_coordinates=" << y_in_IJK_coordinates;
      streamlog_out(MESSAGE) << " z_in_IJK_coordinates=" << z_in_IJK_coordinates;
      streamlog_out(MESSAGE) << " energy=" << ecalhit->getEnergy();
		if (z_in_IJK_coordinates<9)
		{
		totalEnergy=totalEnergy+(ecalhit->getEnergy()*(1+47.894));		
      	streamlog_out(MESSAGE) << " LESSSSSSSSSSSSSSSSSSSSSSSSSSSSS";
		}
		else
		{
		totalEnergy=totalEnergy+(ecalhit->getEnergy()*(1+(47.894*4/3)));		
      	streamlog_out(MESSAGE) << " MOREEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";
		}
      streamlog_out(MESSAGE) << " energy TUNGSTEN=" << totalEnergy;
		
		
		_xHist->Fill(x_in_IJK_coordinates);
		_yHist->Fill(y_in_IJK_coordinates);
		_zHist->Fill(z_in_IJK_coordinates);
		_cellEnergyHist->Fill(ecalhit->getEnergy());
		_xyHist->Fill(x_in_IJK_coordinates,y_in_IJK_coordinates);

    }
	_evEnergyHist->Fill(totalEnergy);
	AIDAProcessor::tree(this);

}


void ExampleProcessor::processRunHeader(LCRunHeader *run)
{
}

void ExampleProcessor::processEvent(LCEvent *evt)
{

	try
	  {
	    streamlog_out(MESSAGE) << "\n ----------------------------------------- ";
	    LCCollection *mccol = evt->getCollection(_MCColName);
	    ShowMCInfo(mccol);
	    
	    LCCollection *ecal = evt->getCollection(_ECALColName);
	    ShowECALInfo(ecal);

	
	}catch (DataNotAvailableException &e)
	{
		streamlog_out(DEBUG) << "Whoops!....\n";
		streamlog_out(DEBUG) << e.what();
	}
	//AIDAProcessor::tree(this);

}

	void ExampleProcessor::check(LCEvent * evt)
	{
		// nothing to check here - could be used to fill checkplots in reconstruction processor
	}

	void ExampleProcessor::end()
	{
	}
