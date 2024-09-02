#include "GetMIPProcessor.hh"
// #include "langaus.C"

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
// using namespace lcio ;
// using namespace marlin ;
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
	registerInputCollection(LCIO::MCPARTICLE, "MCParticle",
							"Muon MC collection",
							_MCColName,
							std::string("MCParticle"));
    registerInputCollection(LCIO::SIMCALORIMETERHIT, "SiEcalCollection",
                            "Sim ECAL Monolithic Collection",
                            _ECALColName,
                            std::string("SiEcalCollection"));
	registerInputCollection(LCIO::SIMCALORIMETERHIT, "PixelSiEcalCollection",
                            "Sim ECAL Pixelised Collection",
                            _pECALColName,
							std::string("PixelSiEcalCollection"));//Name of collection after using the Pixelization Processor, giving coordinates of hit in I,J,K starting from 1
}

GetMIPProcessor::~GetMIPProcessor() {}

void GetMIPProcessor::init()
{
	printParameters();
	AIDAProcessor::tree(this);//Using the AIDAProcessor to save the histograms created in init() in a root file
    _xHist = new TH1D("xHist","X Distribution; x [cell]; Number of hit", NUMBER_OF_CELLX, 0.5, NUMBER_OF_CELLX +.5);//Histogram of X distribution of hits in ECAL pixel coordinates
    _yHist = new TH1D("yHist","Y Distribution; y [cell]; Number of hit", NUMBER_OF_CELLY, 0.5, NUMBER_OF_CELLY +.5);//Histogram of Y distribution of hits
    _zHist = new TH1D("zHist","Z Distribution; z [layer]; Number of hit", NUMBER_OF_LAYER, 0.5, NUMBER_OF_LAYER +.5);//Histogram of Z (Layer) distribution of hits
    _cellEnergyHist = new TH1F("cellEnergyHist","Energy deposited in cells Distribution", 200, FIT_RANGE_MIN, FIT_RANGE_MAX);//Histogram of the energy deposition in all cell for all events
	
    _xyHist = new TH2D("xyHist","XY view all events", NUMBER_OF_CELLX, 0.5, NUMBER_OF_CELLX +.5, NUMBER_OF_CELLY, 0.5, NUMBER_OF_CELLY +.5);//Front view of the ECAL, XY distribution of hits
	

    for (int i = 0; i < NUMBER_OF_LAYER; i++) {
        _energyInLayerSi[i] = new TH1F(Form("energyInLayerSi_%d",i+1),"Energy deposited in layer ", 200, FIT_RANGE_MIN, FIT_RANGE_MAX);//Histogram for the energy deposited in a layer for each event
        _energyInLayerSi[i]->SetTitle(Form("Total energy in layer %d",i+1));
        _energyInPixelLayerSi[i] = new TH1F(Form("energyInPixelLayerSi_%d",i+1),"Energy deposited in pixelised layer ", 200, FIT_RANGE_MIN, FIT_RANGE_MAX);//Histogram for the energy deposited in a layer for each event
        _energyInPixelLayerSi[i]->SetTitle(Form("Total energy in pixelised layer %d",i+1));
	}
    _fittedMIP = new TH1D("fitMIP", "Fitted MIP based on hits per layer", NUMBER_OF_LAYER, 0.5, NUMBER_OF_LAYER +.5);
    _fittedPixelMIP = new TH1D("fitPixelMIP", "Fitted MIP based on hits per pixel", NUMBER_OF_LAYER, 0.5, NUMBER_OF_LAYER +.5);


}


void GetMIPProcessor::ShowMCInfo(EVENT::LCCollection *myCollection) {
    int number = myCollection->getNumberOfElements();
  
    for (int i = 0; i < number; i++) {//Loop through the MC Particle collection for one event
        MCParticle *particle = dynamic_cast<MCParticle *>(myCollection->getElementAt(i));
        vector<MCParticle *> daughters = particle->getDaughters();
    
        streamlog_out(DEBUG) << "\n MCCollection, particle:" << i;
        streamlog_out(DEBUG) << " pdg = " << particle->getPDG() <<",";
        streamlog_out(DEBUG) << " status = " << particle->getGeneratorStatus() <<",";
        streamlog_out(DEBUG) << " N_daughters = " << daughters.size() <<",";
        streamlog_out(DEBUG) << " E = " << particle->getEnergy() <<" GeV,";
        streamlog_out(DEBUG) << " px = " << particle->getMomentum()[0] <<" GeV,";
        streamlog_out(DEBUG) << " py = " << particle->getMomentum()[1] <<" GeV,";
        streamlog_out(DEBUG) << " pz = " << particle->getMomentum()[2] <<" GeV,";
        streamlog_out(DEBUG) << " m = " << particle->getMass() <<" GeV,";
        streamlog_out(DEBUG) << " charge = " << particle->getCharge() <<".";
    }
  	streamlog_out(DEBUG) << std::endl;
}

void GetMIPProcessor::ShowECALInfo(EVENT::LCCollection *myCollection) {
    int number = myCollection->getNumberOfElements();
    CellIDDecoder<EVENT::SimCalorimeterHit> cd(myCollection);
    
    float totalEnergyLayerSi[NUMBER_OF_LAYER];
    int hitsInLayer[NUMBER_OF_LAYER];
    std::fill(std::begin(totalEnergyLayerSi), std::end(totalEnergyLayerSi), 0.0);
    std::fill(std::begin(hitsInLayer), std::end(hitsInLayer), 0.0);

    for (int i = 0; i < number; i++) {
        SimCalorimeterHit *ecalhit = dynamic_cast<SimCalorimeterHit *>(myCollection->getElementAt(i));

        int xyz_x = cd(ecalhit)["x"];
        int xyz_y = cd(ecalhit)["y"];
        int xyz_z = cd(ecalhit)["layer"];
        float hit_energy = ecalhit->getEnergy();

        streamlog_out(DEBUG) << "\n SimCalorimeterHit, :" << i;
        streamlog_out(DEBUG) << " cellID-encoded=" << ecalhit->getCellID0();
        streamlog_out(DEBUG) << " x = " << xyz_x <<" mm,";
        streamlog_out(DEBUG) << " y = " << xyz_y <<" mm,";
        streamlog_out(DEBUG) << " z = " << xyz_z <<" layer,";
        streamlog_out(DEBUG) << " energy = " << hit_energy <<"GeV.\n";
        totalEnergyLayerSi[xyz_z] += hit_energy;
        hitsInLayer[xyz_z]++;
    }
    // return totalEnergyLayerSi;
    for (int i = 0; i < NUMBER_OF_LAYER; i++) {
        if (hitsInLayer[i]==1) {//Filling energy in layer histograms only with events with one hit per layer (muon) to calculate the MIP
            _energyInLayerSi[i]->Fill(totalEnergyLayerSi[i]);
        }
	}
}//By this point all histograms are filled for one event, this is repeated for all the events in the collection

void GetMIPProcessor::ShowPixelECALInfo(EVENT::LCCollection *myCollection) {
    int number = myCollection->getNumberOfElements();
    CellIDDecoder<EVENT::SimCalorimeterHit> cd(myCollection);

    float totalEnergyLayerSi[NUMBER_OF_LAYER];
    int hitsInLayer[NUMBER_OF_LAYER];
    std::fill(std::begin(totalEnergyLayerSi), std::end(totalEnergyLayerSi), 0.0);
    std::fill(std::begin(hitsInLayer), std::end(hitsInLayer), 0);

    for (int i = 0; i < number; i++) {
        SimCalorimeterHit *ecalhit = dynamic_cast<SimCalorimeterHit *>(myCollection->getElementAt(i));

        int IJK_I = cd(ecalhit)["I"];
        int IJK_J = cd(ecalhit)["J"];
        int IJK_K = cd(ecalhit)["K"];
        float hit_energy = ecalhit->getEnergy();

        streamlog_out(DEBUG) << "\n SimCalorimeterHit, :" << i;
        streamlog_out(DEBUG) << " cellID-encoded=" << ecalhit->getCellID0();
        streamlog_out(DEBUG) << " I = " << IJK_I <<" mm,";
        streamlog_out(DEBUG) << " J = " << IJK_J <<" mm,";
        streamlog_out(DEBUG) << " K = " << IJK_K <<" layer,";
        streamlog_out(DEBUG) << " energy = " << hit_energy <<" GeV.\n";
        totalEnergyLayerSi[IJK_K-1] += hit_energy;
        hitsInLayer[IJK_K-1]++;
        _xHist->Fill(IJK_I);//Fill the x, y, z and energy distribution histograms
        _yHist->Fill(IJK_J);
        _zHist->Fill(IJK_K);
        _cellEnergyHist->Fill(hit_energy);
        _xyHist->Fill(IJK_I, IJK_J);
    }
    // return totalEnergyLayerSi;
    for (int i = 0; i < NUMBER_OF_LAYER; i++) {
        if (hitsInLayer[i]==1) {//Filling energy in layer histograms only with events with one hit per layer (muon) to calculate the MIP
            _energyInPixelLayerSi[i]->Fill(totalEnergyLayerSi[i]);
        }
	}
}//By this point all histograms are filled for one event, this is repeated for all the events in the collection


void GetMIPProcessor::processRunHeader(LCRunHeader *run) {
}

void GetMIPProcessor::processEvent(LCEvent *evt)
{

	try {
	    streamlog_out(DEBUG) << "\n ----------------------------------------- ";
	    LCCollection *mccol = evt->getCollection(_MCColName);
	    ShowMCInfo(mccol);
	    LCCollection *ecal = evt->getCollection(_ECALColName);
	    ShowECALInfo(ecal);
	    LCCollection *pecal = evt->getCollection(_pECALColName);
	    ShowPixelECALInfo(pecal);
    } catch (DataNotAvailableException &e) {
		streamlog_out(DEBUG) << "Whoops!....\n";
		streamlog_out(DEBUG) << e.what();
	}

}

void GetMIPProcessor::check(LCEvent * evt) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void GetMIPProcessor::end() {//Using this function to fit the energy in layer histograms after all the events have been processed 
    std::vector<float> MIP;
    std::vector<float> pMIP;

    landauFunc->SetParLimits(0, 0, FIT_INIT_AMP*10.);
    landauFunc->SetParLimits(1, 0, FIT_INIT_MPV*10.);
    // landauFunc->SetParLimits(2, FIT_INIT_SIG/2., FIT_INIT_SIG*2.);
    for (int i = 0; i < NUMBER_OF_LAYER; i++) {//For each layer
        // Fitting SNR histo
        streamlog_out(MESSAGE) << "Fitting 1..." <<std::endl;
        streamlog_out(MESSAGE) << "Items inside the histogram " << _energyInLayerSi[i]->GetEntries() <<std::endl;

        landauFunc->SetParameters(FIT_INIT_AMP, FIT_INIT_MPV, FIT_INIT_SIG);
        _energyInLayerSi[i]->Fit("fitLandauFunc");	//Fit a landau to the distribution
        TF1 *landaufit = (TF1*)_energyInLayerSi[i]->GetListOfFunctions()->FindObject("fitLandauFunc");
        gStyle->SetOptFit(1111);//Set to 1 to show and save the fit with the histogram in the root file generated by the AIDAProcessor
        streamlog_out(MESSAGE) << "Fitting done" <<std::endl << "Plotting results..." <<std::endl;
        for (int j = 0; j < 3; j++) {//Save all fit parameters in the array
            _layerFitParams[i][j] = landaufit->GetParameter(j);
        }
        
        streamlog_out(MESSAGE) << "Fitting 2..." <<std::endl;
        streamlog_out(MESSAGE) << "Items inside the histogram " << _energyInPixelLayerSi[i]->GetEntries() <<std::endl;

        landauFunc->SetParameters(FIT_INIT_AMP, FIT_INIT_MPV, FIT_INIT_SIG);
        _energyInPixelLayerSi[i]->Fit("fitLandauFunc");	//Fit a landau to the distribution
        TF1 *landaufit2 = (TF1*)_energyInPixelLayerSi[i]->GetListOfFunctions()->FindObject("fitLandauFunc");
        gStyle->SetOptFit(1111);//Set to 1 to show and save the fit with the histogram in the root file generated by the AIDAProcessor
        streamlog_out(MESSAGE) << "Fitting done" <<std::endl << "Plotting results..." <<std::endl;
        for (int j = 0; j < 3; j++) {//Save all fit parameters in the array
            _pixelLayerFitParams[i][j] = landaufit2->GetParameter(j);
        }
    }

    for (int i = 0; i < NUMBER_OF_LAYER; i++) {
        streamlog_out(MESSAGE) << "\n Fit PARAMS for layer " << i;//Printing all fit parameters for each layer
        for (int j = 0; j < 3; j++) {
            streamlog_out(MESSAGE) << "\n par["<< j <<"] = "<< _layerFitParams[i][j];
        }
        for (int j = 0; j < 3; j++) {
            streamlog_out(MESSAGE) << "\n par["<< j <<"] = "<< _pixelLayerFitParams[i][j];
        }
        MIP.push_back(_layerFitParams[i][1]);
        _fittedMIP->SetBinContent(i+1, _layerFitParams[i][1]); // bin number counts from 1
        pMIP.push_back(_pixelLayerFitParams[i][1]);
        _fittedPixelMIP->SetBinContent(i+1, _pixelLayerFitParams[i][1]);
    }
    std::cout << '\n';

    // Print out the vector
    for (float n : MIP) std::cout << n << ' ';
    std::cout << '\n';
    for (float n : pMIP) std::cout << n << ' ';
    std::cout << '\n';
}
