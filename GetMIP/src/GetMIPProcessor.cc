#include "GetMIPProcessor.hh"
#include "langaus.C"

// ROOT
#include "TStyle.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TMath.h"

#include <math.h>
#include <iostream>
#include <fstream>

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

GetMIPProcessor aGetMIPProcessor;

GetMIPProcessor::GetMIPProcessor() : Processor("GetMIPProcessor")
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

GetMIPProcessor::~GetMIPProcessor() {}

void GetMIPProcessor::init()
{
	int testVariable = 1;
	int test2 = 10;
	printParameters();
	_xHist = new TH1D("_xHist","X Distribution",32, -0.5, 31.5);
	_yHist = new TH1D("_yHist","Y Distribution",32, -0.5, 31.5);
	_zHist = new TH1D("_zHist","Z Distribution",15, -0.5, 14.5);
	_cellEnergyHist = new TH1F("_cellEnergyHist","Energy deposited in cells Distribution",200, 0, 0.002);
	_evEnergyHist = new TH1F("_evEnergyHist","Energy of shower Distribution",100, 0, 25);
	
	_xyHist = new TH2D("_xyHist","XY view all events",30,-0.5,29.5,30, -0.5, 29.5);
	

	for (int i = 0; i < 15; i++)
	{
		_energyInLayerSi[i] = new TH1F(Form("_energyInLayerSi_%d",i+1),"Energy deposited in layer ",200, 0, 0.002);
		_energyInLayerSi[i]->SetTitle(Form("Total energy in layer %d",i+1));
	}
	AIDAProcessor::tree(this);
	
}


void GetMIPProcessor::ShowMCInfo(EVENT::LCCollection *myCollection)
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

 void GetMIPProcessor::ShowECALInfo(EVENT::LCCollection *myCollection)
{
  int number = myCollection->getNumberOfElements();
  CellIDDecoder<EVENT::SimCalorimeterHit> cd(myCollection);
	double totalEnergy = 0;
	double totalEnergyLayerSi[15] = {0};
	int hitsInLayer[15] = {0};
  for (int i = 0; i < number; i++)
    {

      SimCalorimeterHit *ecalhit = dynamic_cast<SimCalorimeterHit *>(myCollection->getElementAt(i));
      int x_in_IJK_coordinates = cd(ecalhit)["x"];
      int y_in_IJK_coordinates = cd(ecalhit)["y"];
      int z_in_IJK_coordinates = cd(ecalhit)["layer"];
		//_coordinateZ=z_in_IJK_coordinates;
 /*     streamlog_out(MESSAGE) << "\n SimCalorimeterHit, :" << i;
      streamlog_out(MESSAGE) << " cellID-encoded=" << ecalhit->getCellID0();
      streamlog_out(MESSAGE) << " x_in_IJK_coordinates=" << x_in_IJK_coordinates;
      streamlog_out(MESSAGE) << " y_in_IJK_coordinates=" << y_in_IJK_coordinates;
      streamlog_out(MESSAGE) << " z_in_IJK_coordinates=" << z_in_IJK_coordinates;
      streamlog_out(MESSAGE) << " energy=" << ecalhit->getEnergy();
      streamlog_out(MESSAGE) << " NUMBER=" << number;
*/		
		if (z_in_IJK_coordinates<9)
		{
			if (z_in_IJK_coordinates<5)
			{
				totalEnergy=totalEnergy+(ecalhit->getEnergy()*(1+47.894*4.2/0.650));		
			}
			else
			{
				totalEnergy=totalEnergy+(ecalhit->getEnergy()*(1+47.894*4.2/0.500));		
			}
			
		}
		else
		{
			if (z_in_IJK_coordinates<11)
			{
				totalEnergy=totalEnergy+(ecalhit->getEnergy()*(1+(47.894*5.6/500)));		
			}
			else
			{
				totalEnergy=totalEnergy+(ecalhit->getEnergy()*(1+47.894*5.6/0.320));		
			}
		}
		
		totalEnergyLayerSi[z_in_IJK_coordinates]=totalEnergyLayerSi[z_in_IJK_coordinates]+ecalhit->getEnergy();
		hitsInLayer[z_in_IJK_coordinates] = hitsInLayer[z_in_IJK_coordinates] + 1;
		_xHist->Fill(x_in_IJK_coordinates);
		_yHist->Fill(y_in_IJK_coordinates);
		_zHist->Fill(z_in_IJK_coordinates);
		_cellEnergyHist->Fill(ecalhit->getEnergy());
		_xyHist->Fill(x_in_IJK_coordinates,y_in_IJK_coordinates);

    }
    streamlog_out(MESSAGE) << " energy TUNGSTEN=" << totalEnergy;
	_evEnergyHist->Fill(totalEnergy);
	for (int i = 0; i < 15; i++)
	{
		if (hitsInLayer[i]==1)
		{
			_energyInLayerSi[i]->Fill(totalEnergyLayerSi[i]);
			totalEnergyLayerSi[i]=0;
		}
		
	}

	totalEnergy=0;

}


void GetMIPProcessor::processRunHeader(LCRunHeader *run)
{
}

void GetMIPProcessor::processEvent(LCEvent *evt)
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

	void GetMIPProcessor::check(LCEvent * evt)
	{
		// nothing to check here - could be used to fill checkplots in reconstruction processor
	}

	void GetMIPProcessor::end()
	{
		for (int i = 0; i < 15; i++)
		{
				// Fitting SNR histo
			printf("Fitting...\n");
			
			_energyInLayerSi[i]->Fit("landau");	
			TF1 *fit = (TF1*)_energyInLayerSi[i]->GetListOfFunctions()->FindObject("landau");
			gStyle->SetOptFit(1111);
			printf("Fitting done\nPlotting results...\n");
			for (int j = 0; j < 3; j++)
			{
				_layerFitParams[i][j] = fit->GetParameter(j);
			}
			
		}
		for (int j = 0; j < 15; j++)
		{
			streamlog_out(MESSAGE) << "\n Fit PARAMS for layer " << j;
			for (int h = 0; h < 3; h++)
			{
				streamlog_out(MESSAGE) << "\n par["<< h <<"] = "<< _layerFitParams[j][h];
			
			}
		}	
	}
