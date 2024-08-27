#include "GetEnergyResolutionProcessor.hh"
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

GetEnergyResolutionProcessor aGetEnergyResolutionProcessor;

GetEnergyResolutionProcessor::GetEnergyResolutionProcessor() : Processor("GetEnergyResolutionProcessor")
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
							std::string("PixelSiEcalCollection"));
}

GetEnergyResolutionProcessor::~GetEnergyResolutionProcessor() {}

void GetEnergyResolutionProcessor::init()
{
	AIDAProcessor::tree(this);
	printParameters();
	_xHist = new TH1F("_xHist","X Distribution",64, 0.5, 64.5);
	_yHist = new TH1F("_yHist","Y Distribution",32, 0.5, 32.5);
	_zHist = new TH1D("_zHist","Z Distribution",15, 0.5, 15.5);
	_cellEnergyHist = new TH1F("_cellEnergyHist","Energy deposited in cells Distribution",100, 0, 0.05);
	_evEnergyHist = new TH1F("_evEnergyHist","Energy of shower Distribution",100, 0, 0.2);
	_evHitsHist = new TH1D("_evHitsHist","Number of hits Distribution", 62, -5, 305);
	
	_xyHist = new TH2D("_xyHist","XY view all events",64,0.5,64.5,32, 0.5, 32.5);
	_zxHist = new TH2D("_zxHist","ZX view all events",15,0.5,15.5,64, 0.5, 64.5);
	_zyHist = new TH2D("_zyHist","ZY view all events",15,0.5,15.5,32, 0.5, 32.5);
	

	for (int i = 0; i < 15; i++)
	{
		_energyInLayerSi[i] = new TH1F(Form("_energyInLayerSi_%d",i+1),"Energy deposited in layer ",100, 0, 0.05);
		_energyInLayerSi[i]->SetTitle(Form("Total energy in layer %d",i+1));
		
		_hitsInLayer[i] = new TH1D(Form("_HitsInLayer_%d",i+1),"Hits in layer ",22, -5, 105);
		_hitsInLayer[i]->SetTitle(Form("Total hits in layer %d",i+1));
	}
	
}


void GetEnergyResolutionProcessor::ShowMCInfo(EVENT::LCCollection *myCollection)
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

 void GetEnergyResolutionProcessor::ShowECALInfo(EVENT::LCCollection *myCollection)
{
  int number = myCollection->getNumberOfElements();
  CellIDDecoder<EVENT::SimCalorimeterHit> cd(myCollection);
	double totalEnergy = 0;
	double totalEnergyLayerSi[15] = {0};
	int hitsInLayer[15] = {0};
  for (int i = 0; i < number; i++)
    {

      SimCalorimeterHit *ecalhit = dynamic_cast<SimCalorimeterHit *>(myCollection->getElementAt(i));
    int x_in_IJK_coordinates = cd(ecalhit)["I"];
    int y_in_IJK_coordinates = cd(ecalhit)["J"];
      int z_in_IJK_coordinates = cd(ecalhit)["K"];
		//_coordinateZ=z_in_IJK_coordinates;
      streamlog_out(MESSAGE) << "\n SimCalorimeterHit, :" << i;
      streamlog_out(MESSAGE) << " cellID-encoded=" << ecalhit->getCellID0();
      streamlog_out(MESSAGE) << " x_in_IJK_coordinates=" << x_in_IJK_coordinates;
      streamlog_out(MESSAGE) << " y_in_IJK_coordinates=" << y_in_IJK_coordinates;
      streamlog_out(MESSAGE) << " z_in_IJK_coordinates=" << z_in_IJK_coordinates;
      streamlog_out(MESSAGE) << " energy=" << ecalhit->getEnergy();
      streamlog_out(MESSAGE) << " NUMBER=" << number;
		
		totalEnergy=totalEnergy+ecalhit->getEnergy();
		noHits++;
		totalEnergyLayerSi[z_in_IJK_coordinates-1]=totalEnergyLayerSi[z_in_IJK_coordinates-1]+ecalhit->getEnergy();
		hitsInLayer[z_in_IJK_coordinates-1] = hitsInLayer[z_in_IJK_coordinates-1] + 1;
		_xHist->Fill(x_in_IJK_coordinates);
		_yHist->Fill(y_in_IJK_coordinates);
		_zHist->Fill(z_in_IJK_coordinates);
		_cellEnergyHist->Fill(ecalhit->getEnergy());
		_xyHist->Fill(x_in_IJK_coordinates,y_in_IJK_coordinates);
		_zxHist->Fill(z_in_IJK_coordinates,x_in_IJK_coordinates);
		_zyHist->Fill(z_in_IJK_coordinates,y_in_IJK_coordinates);

    }
    streamlog_out(MESSAGE) << " energy TUNGSTEN=" << totalEnergy;
	_evEnergyHist->Fill(totalEnergy);
	_evHitsHist->Fill(noHits);
	noHits=0;
	totalEnergy=0;
	for (int i = 0; i < 15; i++)
	{
		if (hitsInLayer[i]>0)
		{
			_energyInLayerSi[i]->Fill(totalEnergyLayerSi[i]);
		}
		_hitsInLayer[i]->Fill(hitsInLayer[i]);
		totalEnergyLayerSi[i]=0;
		
	}


}


void GetEnergyResolutionProcessor::processRunHeader(LCRunHeader *run)
{
}

void GetEnergyResolutionProcessor::processEvent(LCEvent *evt)
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

	void GetEnergyResolutionProcessor::check(LCEvent * evt)
	{
		// nothing to check here - could be used to fill checkplots in reconstruction processor
	}

	void GetEnergyResolutionProcessor::end()
	{
			_evHitsHist->Fit("gaus");	//Fit a gaussian to the distribution
			TF1 *fit1 = (TF1*)_evHitsHist->GetListOfFunctions()->FindObject("gaus");
			gStyle->SetOptFit(1111);//Set to 1 to show and save the fit with the histogram in the root file generated by the AIDAProcessor
			_evHitsHist->Fit("gaus","","",fit1->GetParameter(1)-fit1->GetParameter(2),fit1->GetParameter(1)+fit1->GetParameter(2));	//Fit a gaussian to the distribution
			TF1 *fit2 = (TF1*)_evHitsHist->GetListOfFunctions()->FindObject("gaus");
			gStyle->SetOptFit(1111);//Set to 1 to show and save the fit with the histogram in the root file generated by the AIDAProcessor

			_evEnergyHist->Fit("gaus");	//Fit a gaussian to the distribution
			TF1 *fit = (TF1*)_evEnergyHist->GetListOfFunctions()->FindObject("gaus");
			gStyle->SetOptFit(1111);//Set to 1 to show and save the fit with the histogram in the root file generated by the AIDAProcessor

			printf("Fitting done\nPlotting results...\n");
			streamlog_out(MESSAGE) << "\n Fit " << fit->GetParameter(2);

		for (int i = 0; i < 15; i++)
		{
			// Fitting SNR histo
			printf("Fitting...\n");

			_hitsInLayer[i]->Fit("gaus");	
			TF1 *fitLayer1 = (TF1*)_hitsInLayer[i]->GetListOfFunctions()->FindObject("gaus");
			gStyle->SetOptFit(1111);
			
			_energyInLayerSi[i]->Fit("gaus");	
			TF1 *fitLayer = (TF1*)_energyInLayerSi[i]->GetListOfFunctions()->FindObject("gaus");
			gStyle->SetOptFit(1111);
			printf("Fitting done\nPlotting results...\n");
		}
			
	}
