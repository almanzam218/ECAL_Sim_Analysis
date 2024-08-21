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
							std::string("PixelSiEcalCollection"));//Name of collection after using the Pixelization Processor, giving coordinates of hit in I,J,K starting from 1
}

GetMIPProcessor::~GetMIPProcessor() {}

void GetMIPProcessor::init()
{
	printParameters();
	AIDAProcessor::tree(this);//Using the AIDAProcessor to save the histograms created in init() in a root file
	_xHist = new TH1D("_xHist","X Distribution",64, 0.5, 64.5);//Histogram of X distribution of hits in ECAL pixel coordinates
	_yHist = new TH1D("_yHist","Y Distribution",32, 0.5, 32.5);//Histogram of Y distribution of hits
	_zHist = new TH1D("_zHist","Z Distribution",15, 0.5, 15.5);//Histogram of Z (Layer) distribution of hits
	_cellEnergyHist = new TH1F("_cellEnergyHist","Energy deposited in cells Distribution",200, 0, 0.002);//Histogram of the energy deposition in all cell for all events
	
	_xyHist = new TH2D("_xyHist","XY view all events",64,0.5,64.5,32, 0.5, 32.5);//Front view of the ECAL, XY distribution of hits
	

	for (int i = 0; i < 15; i++)
	{
		_energyInLayerSi[i] = new TH1F(Form("_energyInLayerSi_%d",i+1),"Energy deposited in layer ",200, 0, 0.002);//Histogram for the energy deposited in a layer for each event
		_energyInLayerSi[i]->SetTitle(Form("Total energy in layer %d",i+1));
	}

}


void GetMIPProcessor::ShowMCInfo(EVENT::LCCollection *myCollection)
{
  int number = myCollection->getNumberOfElements();
  
    for (int i = 0; i < number; i++)//Loop through the MC Particle collection for one event
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

	double totalEnergyLayerSi[15] = {0};
	int hitsInLayer[15] = {0};
  for (int i = 0; i < number; i++)//Loop through the ECAL Hits in one event (after pixelization)
    {

      SimCalorimeterHit *ecalhit = dynamic_cast<SimCalorimeterHit *>(myCollection->getElementAt(i));
		int x_in_IJK_coordinates = cd(ecalhit)["I"];//Reading the x, y and layer coordinate of one hit in the collection
    	int y_in_IJK_coordinates = cd(ecalhit)["J"];
        int z_in_IJK_coordinates = cd(ecalhit)["K"];
		//_coordinateZ=z_in_IJK_coordinates;
 /*     streamlog_out(MESSAGE) << "\n SimCalorimeterHit, :" << i;
      streamlog_out(MESSAGE) << " cellID-encoded=" << ecalhit->getCellID0();
      streamlog_out(MESSAGE) << " x_in_IJK_coordinates=" << x_in_IJK_coordinates;
      streamlog_out(MESSAGE) << " y_in_IJK_coordinates=" << y_in_IJK_coordinates;
      streamlog_out(MESSAGE) << " z_in_IJK_coordinates=" << z_in_IJK_coordinates;
      streamlog_out(MESSAGE) << " energy=" << ecalhit->getEnergy();
      streamlog_out(MESSAGE) << " NUMBER=" << number;
*/		
		
		totalEnergyLayerSi[z_in_IJK_coordinates-1]=totalEnergyLayerSi[z_in_IJK_coordinates-1]+ecalhit->getEnergy();//Adding the energy of the event separated by layers
		hitsInLayer[z_in_IJK_coordinates-1]++;//Counting the hits per layer in one event
		_xHist->Fill(x_in_IJK_coordinates);//Fill the x, y, z and energy distribution histograms
		_yHist->Fill(y_in_IJK_coordinates);
		_zHist->Fill(z_in_IJK_coordinates);
		_cellEnergyHist->Fill(ecalhit->getEnergy());
		_xyHist->Fill(x_in_IJK_coordinates,y_in_IJK_coordinates);

    }


	for (int i = 0; i < 15; i++)
	{
		if (hitsInLayer[i]==1)//Filling energy in layer histograms only with events with one hit per layer (muon) to calculate the MIP
		{
			_energyInLayerSi[i]->Fill(totalEnergyLayerSi[i]);
			totalEnergyLayerSi[i]=0;
		}
		
	}

	

}//By this point all histograms are filled for one event, this is repeated for all the events in the collection


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

}

	void GetMIPProcessor::check(LCEvent * evt)
	{
		// nothing to check here - could be used to fill checkplots in reconstruction processor
	}

	void GetMIPProcessor::end()//Using this function to fit the energy in layer histograms after all the events have been processed 
	{
		std::vector<float> MIP;
		for (int i = 0; i < 15; i++)//For each layer
		{
				// Fitting SNR histo
			printf("Fitting...\n");
			
			_energyInLayerSi[i]->Fit("landau");	//Fit a landau to the distribution
			TF1 *fit = (TF1*)_energyInLayerSi[i]->GetListOfFunctions()->FindObject("landau");
			gStyle->SetOptFit(1111);//Set to 1 to show and save the fit with the histogram in the root file generated by the AIDAProcessor
			printf("Fitting done\nPlotting results...\n");
			for (int j = 0; j < 3; j++)//Save all fit parameters in the array
			{
				_layerFitParams[i][j] = fit->GetParameter(j);
			}
			
		}
		for (int j = 0; j < 15; j++)
		{
			streamlog_out(MESSAGE) << "\n Fit PARAMS for layer " << j;//Printing all fit parameters for each layer
			for (int h = 0; h < 3; h++)
			{
				streamlog_out(MESSAGE) << "\n par["<< h <<"] = "<< _layerFitParams[j][h];
			}
				MIP.push_back(_layerFitParams[j][1]);
		}	
    std::cout << '\n';
	
	// Print out the vector
    for (float n : MIP)
        std::cout << n << ' ';
    std::cout << '\n';

	}
