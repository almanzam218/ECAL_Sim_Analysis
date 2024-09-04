#include "GetEnergyResolutionProcessor.hh"
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

GetEnergyResolutionProcessor::GetEnergyResolutionProcessor() : Processor("GetEnergyResolutionProcessor") {

	// modify processor description
	_description = "";

	// input collections
	registerInputCollection(LCIO::MCPARTICLE,"MCCollectionName",
                            "Name of the MC collection",
							_MCColName,
							// std::string("MCParticle"));
							std::string("Not configured in xml file"));
	registerInputCollection(LCIO::SIMCALORIMETERHIT, "ECALCollectionName",
							"Monolithic ECAL Hits Collection",
							_ECALColName,
							// std::string("SiEcalCollection"));
							std::string("Not configured in xml file"));
	registerInputCollection(LCIO::SIMCALORIMETERHIT, "PixelisedECALCollectionName",
							"Pixelised ECAL Hits Collection",
							_pECALColName,
							// std::string("PixelSiEcalCollection"));
							std::string("Not configured in xml file"));
	registerInputCollection(LCIO::CALORIMETERHIT, "DigitisedECALCollectionName",
							"Digitised ECAL Hits Collection",
							_dECALColName,
							// std::string("EcalCollection"));
							std::string("Not configured in xml file"));
}

GetEnergyResolutionProcessor::~GetEnergyResolutionProcessor() {}

void GetEnergyResolutionProcessor::init() {
	AIDAProcessor::tree(this);
	printParameters();
    if (_MCColName!="Not configured in xml file") {_flagMcCol = true;}
    if (_ECALColName!="Not configured in xml file") {_flagEcalCol = true;}
    if (_pECALColName!="Not configured in xml file") {_flagPixelEcalCol = true;}
    if (_dECALColName!="Not configured in xml file") {_flagDigitEcalCol = true;}
    // MC particle
    if (_flagMcCol) {
        _runEnergy = new TH1F("_runInfo", "Run Information", 1, 0, 1); // Bin 1: Beam energy‘        
    }

    // Monolithic calorimeter hits
    if (_flagEcalCol) {
        _zMonoHist = new TH1D("mono_zHist","Z Distribution; z [layer]; Number of hit", NUMBER_OF_LAYER, 0.5, NUMBER_OF_LAYER +.5);//Histogram of Z (Layer) distribution of hits
        _evMonoEnergyHist = new TH1D("mono_evEnergyHist","Energy of shower Distribution; E_{dep} [GeV]", evHistBins, 0, evHistBins);
        for (int i = 0; i < NUMBER_OF_LAYER; i++) {
            _energyInMonoLayerSi[i] = new TH1F(Form("mono_energyInLayerSi_%d",i+1),"Energy deposited in monolithic layer; E_{dep} [GeV];",100, 0, 0.05);
        }
    }

    // Pixelised calorimeter hits
    if (_flagPixelEcalCol) {
        _xHist = new TH1D("pixel_xHist","X Distribution; x [cell]; Number of hit", NUMBER_OF_CELLX, 0.5, NUMBER_OF_CELLX +.5);//Histogram of X distribution of hits in ECAL pixel coordinates
        _yHist = new TH1D("pixel_yHist","Y Distribution; y [cell]; Number of hit", NUMBER_OF_CELLY, 0.5, NUMBER_OF_CELLY +.5);//Histogram of Y distribution of hits
        _zHist = new TH1D("pixel_zHist","Z Distribution; z [layer]; Number of hit", NUMBER_OF_LAYER, 0.5, NUMBER_OF_LAYER +.5);//Histogram of Z (Layer) distribution of hits
        _xyHist = new TH2F("pixel_xyHist","XY view all events", NUMBER_OF_CELLX, 0.5, NUMBER_OF_CELLX +.5, NUMBER_OF_CELLY, 0.5, NUMBER_OF_CELLY +.5);
        _zxHist = new TH2F("pixel_zxHist","ZX view all events", NUMBER_OF_LAYER, 0.5, NUMBER_OF_LAYER +.5, NUMBER_OF_CELLX, 0.5, NUMBER_OF_CELLX +.5);
        _zyHist = new TH2F("pixel_zyHist","ZY view all events", NUMBER_OF_LAYER, 0.5, NUMBER_OF_LAYER +.5, NUMBER_OF_CELLY, 0.5, NUMBER_OF_CELLY +.5);
        _cellEnergyHist = new TH1F("pixel_cellEnergyHist","Energy deposited in cells Distribution; E_{dep} [GeV]; Number of hit", 200, E_RANGE_MIN, E_RANGE_MAX);//Histogram of the energy deposition in all cell for all events
        // The histogram will instead be declared and filled at the ending stage
        // Bin ranges will be changed at the final stage
        _evEnergyHist = new TH1D("pixel_evEnergyHist","Energy of shower Distribution; E_{dep} [GeV]", evHistBins, 0, evHistBins);
        _evHitsHist = new TH1D("pixel_evHitsHist","Number of hits Distribution; Hit", evHistBins, 0, evHistBins);
        for (int i = 0; i < NUMBER_OF_LAYER; i++) {
            _energyInLayerSi[i] = new TH1F(Form("pixel_energyInLayerSi_%d",i+1),"Energy deposited in layer; E_{dep} [GeV];",100, 0, 0.05);
            _energyInLayerSi[i]->SetTitle(Form("Total energy in layer %d",i+1));
            _hitsInLayer[i] = new TH1F(Form("pixel_HitsInLayer_%d",i+1),"Hits in layer; Hit;",100, -0.5, 99.5);
            _hitsInLayer[i]->SetTitle(Form("Total hits in layer %d",i+1));
        }
    }

    // Digitised calorimeter hits
    if (_flagDigitEcalCol) {
        _xDigitHist = new TH1D("digit_xHist","X Distribution; x [cell]; Number of hit", NUMBER_OF_CELLX, 0.5, NUMBER_OF_CELLX +.5);//Histogram of X distribution of hits in ECAL pixel coordinates
        _yDigitHist = new TH1D("digit_yHist","Y Distribution; y [cell]; Number of hit", NUMBER_OF_CELLY, 0.5, NUMBER_OF_CELLY +.5);//Histogram of Y distribution of hits
        _zDigitHist = new TH1D("digit_zHist","Z Distribution; z [layer]; Number of hit", NUMBER_OF_LAYER, 0.5, NUMBER_OF_LAYER +.5);//Histogram of Z (Layer) distribution of hits
        _xyDigitHist = new TH2F("digit_xyHist","XY view all events", NUMBER_OF_CELLX, 0.5, NUMBER_OF_CELLX +.5, NUMBER_OF_CELLY, 0.5, NUMBER_OF_CELLY +.5);
        _zxDigitHist = new TH2F("digit_zxHist","ZX view all events", NUMBER_OF_LAYER, 0.5, NUMBER_OF_LAYER +.5, NUMBER_OF_CELLX, 0.5, NUMBER_OF_CELLX +.5);
        _zyDigitHist = new TH2F("digit_zyHist","ZY view all events", NUMBER_OF_LAYER, 0.5, NUMBER_OF_LAYER +.5, NUMBER_OF_CELLY, 0.5, NUMBER_OF_CELLY +.5);
        _cellDigitEnergyHist = new TH1F("digit_cellEnergyHist","Energy deposited in cells Distribution; E_{dep} [MIP]; Number of hit", 200, E_RANGE_MIN, E_RANGE_MAX);//Histogram of the energy deposition in all cell for all events
        // The histogram will instead be declared and filled at the ending stage
        // Bin ranges will be changed at the final stage
        _evDigitEnergyHist = new TH1D("digit_evEnergyHist","Energy of shower Distribution; E_{dep} [MIP]", evHistBins, 0, evHistBins);
        _evDigitHitsHist = new TH1D("digit_evHitsHist","Number of hits Distribution; Hit", evHistBins, 0, evHistBins);
        for (int i = 0; i < NUMBER_OF_LAYER; i++) {
            _energyInDigitLayerSi[i] = new TH1F(Form("digit_energyInLayerSi_%d",i+1),"Energy deposited in layer; E_{dep} [MIP];",100, 0, 0.05);
            _energyInDigitLayerSi[i]->SetTitle(Form("Total energy in layer %d",i+1));
            _hitsInDigitLayer[i] = new TH1F(Form("digit_HitsInLayer_%d",i+1),"Hits in layer; Hit;",100, -0.5, 99.5);
            _hitsInDigitLayer[i]->SetTitle(Form("Total hits in layer %d",i+1));
        }
    }
}


void GetEnergyResolutionProcessor::ShowMCInfo(EVENT::LCCollection *myCollection) {
    int number = myCollection->getNumberOfElements();
    
    if (runEnergy==-1) {
        for (int i = 0; i < number; i++) {//Loop through the MC Particle collection for one event
            MCParticle *particle = dynamic_cast<MCParticle *>(myCollection->getElementAt(i));
            vector<MCParticle *> daughters = particle->getDaughters();
            runEnergy = particle->getEnergy();
            
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
}

void GetEnergyResolutionProcessor::ShowECALInfo(EVENT::LCCollection *myCollection) {
    int number = myCollection->getNumberOfElements();
    streamlog_out(DEBUG) << "TOTAL NUMBER OF HITS: " << number <<endl;
    CellIDDecoder<EVENT::SimCalorimeterHit> cd(myCollection);
    
	double totalEnergy = 0;
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
        totalEnergy += hit_energy;
        totalEnergyLayerSi[xyz_z] += hit_energy;
        hitsInLayer[xyz_z]++;
        _zMonoHist->Fill(xyz_z+1);
    }
    // return totalEnergyLayerSi;
    _evMonoEnergyVec.push_back(totalEnergy);
	for (int i = 0; i < NUMBER_OF_LAYER; i++) {
		if (hitsInLayer[i]>0) {
			_energyInMonoLayerSi[i]->Fill(totalEnergyLayerSi[i]);
		}
	}
}//By this point all histograms are filled for one event, this is repeated for all the events in the collection

void GetEnergyResolutionProcessor::ShowPixelECALInfo(EVENT::LCCollection *myCollection) {
    int number = myCollection->getNumberOfElements();
    streamlog_out(DEBUG) << "TOTAL NUMBER OF HITS: " << number <<endl;
    CellIDDecoder<EVENT::SimCalorimeterHit> cd(myCollection);
    
	double totalEnergy = 0;
    int totalHits = 0;
	double totalEnergyLayerSi[NUMBER_OF_LAYER];
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
		totalEnergy += hit_energy;
        totalEnergyLayerSi[IJK_K-1] += hit_energy;
        totalHits++;
        hitsInLayer[IJK_K-1]++;
        _xHist->Fill(IJK_I);
        _yHist->Fill(IJK_J);
        _zHist->Fill(IJK_K);
		_cellEnergyHist->Fill(hit_energy);
		_xyHist->Fill(IJK_I,IJK_J);
		_zxHist->Fill(IJK_K,IJK_I);
		_zyHist->Fill(IJK_K,IJK_J);
    }
    streamlog_out(DEBUG) << "Total energy deposit: " << totalEnergy << " GeV" <<endl;
    // Instead of filling the histograms now, we store the numbers in vectors first, then decide the binsize later
	// _evEnergyHist->Fill(totalEnergy);
	// _evHitsHist->Fill(totalHits);
    _evEnergyVec.push_back(totalEnergy);
    _evHitsVec.push_back(totalHits);
	for (int i = 0; i < NUMBER_OF_LAYER; i++) {
		if (hitsInLayer[i]>0) {
			_energyInLayerSi[i]->Fill(totalEnergyLayerSi[i]);
		}
		_hitsInLayer[i]->Fill(hitsInLayer[i]);
	}
}

void GetEnergyResolutionProcessor::ShowDigitECALInfo(EVENT::LCCollection *myCollection) {
    int number = myCollection->getNumberOfElements();
    streamlog_out(DEBUG) << "TOTAL NUMBER OF HITS: " << number <<endl;
    CellIDDecoder<EVENT::CalorimeterHit> cd(myCollection);

	double totalEnergy = 0;
    int totalHits = 0;
	double totalEnergyLayerSi[NUMBER_OF_LAYER];
	int hitsInLayer[NUMBER_OF_LAYER];
    std::fill(std::begin(totalEnergyLayerSi), std::end(totalEnergyLayerSi), 0.0);
    std::fill(std::begin(hitsInLayer), std::end(hitsInLayer), 0);
    
    for (int i = 0; i < number; i++) {
        CalorimeterHit *ecalhit = dynamic_cast<CalorimeterHit *>(myCollection->getElementAt(i));
        
        int IJK_I = cd(ecalhit)["I"];
        int IJK_J = cd(ecalhit)["J"];
        int IJK_K = cd(ecalhit)["K"];
        float hit_energy = ecalhit->getEnergy();

        streamlog_out(DEBUG) << "\n CalorimeterHit, :" << i;
        streamlog_out(DEBUG) << " cellID-encoded=" << ecalhit->getCellID0();
        streamlog_out(DEBUG) << " I = " << IJK_I <<" mm,";
        streamlog_out(DEBUG) << " J = " << IJK_J <<" mm,";
        streamlog_out(DEBUG) << " K = " << IJK_K <<" layer,";
        streamlog_out(DEBUG) << " energy = " << hit_energy <<" GeV.\n";
		totalEnergy += hit_energy;
        totalEnergyLayerSi[IJK_K-1] += hit_energy;
        totalHits++;
        hitsInLayer[IJK_K-1]++;
        _xDigitHist->Fill(IJK_I);
        _yDigitHist->Fill(IJK_J);
        _zDigitHist->Fill(IJK_K);
		_cellDigitEnergyHist->Fill(hit_energy);
		_xyDigitHist->Fill(IJK_I,IJK_J);
		_zxDigitHist->Fill(IJK_K,IJK_I);
		_zyDigitHist->Fill(IJK_K,IJK_J);
    }
    streamlog_out(DEBUG) << "Total energy deposit: " << totalEnergy << " GeV" <<endl;
    // Instead of filling the histograms now, we store the numbers in vectors first, then decide the binsize later
	// _evEnergyHist->Fill(totalEnergy);
	// _evHitsHist->Fill(totalHits);
    _evDigitEnergyVec.push_back(totalEnergy);
    _evDigitHitsVec.push_back(totalHits);
	for (int i = 0; i < NUMBER_OF_LAYER; i++) {
		if (hitsInLayer[i]>0) {
			_energyInDigitLayerSi[i]->Fill(totalEnergyLayerSi[i]);
		}
		_hitsInDigitLayer[i]->Fill(hitsInLayer[i]);
		totalEnergyLayerSi[i]=0;
	}
}


void GetEnergyResolutionProcessor::processRunHeader(LCRunHeader *run)
{
}

void GetEnergyResolutionProcessor::processEvent(LCEvent *evt) {
    try {
        streamlog_out(DEBUG) << "\n ----------------------------------------- ";
        if (_flagMcCol) {
            LCCollection *mccol = evt->getCollection(_MCColName);
            ShowMCInfo(mccol);
        }
        if (_flagEcalCol) {
            LCCollection *ecal = evt->getCollection(_ECALColName);
            ShowECALInfo(ecal);
        }
        if (_flagPixelEcalCol) {
            LCCollection *pecal = evt->getCollection(_pECALColName);
            ShowPixelECALInfo(pecal);
        }
        if (_flagDigitEcalCol) {
            LCCollection *decal = evt->getCollection(_dECALColName);
            ShowDigitECALInfo(decal);
        }
    } catch (DataNotAvailableException &e) {
        streamlog_out(DEBUG) << "Whoops!....\n";
        streamlog_out(DEBUG) << e.what();
    }
}

void GetEnergyResolutionProcessor::check(LCEvent * evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void GetEnergyResolutionProcessor::end() {
    streamlog_out(MESSAGE) << "Event loop finished. Starting the fit..." <<endl;
    if (_flagMcCol) {
        _runEnergy->SetBinContent(1, runEnergy);
    }
    if (_flagEcalCol) {
        // Using the full statistical information to determine the binsize
        
        float meanMonoEnergy = TMath::Mean(_evMonoEnergyVec.begin(), _evMonoEnergyVec.end());
        float sigmaMonoEnergy = TMath::RMS(_evMonoEnergyVec.begin(), _evMonoEnergyVec.end());
        // Now declaring the histograms
        float minMonoEnergy = meanMonoEnergy - radiusOverSigma*sigmaMonoEnergy;
        float maxMonoEnergy = meanMonoEnergy + radiusOverSigma*sigmaMonoEnergy;
        _evMonoEnergyHist->SetBins(evHistBins, minMonoEnergy, maxMonoEnergy);
        streamlog_out(DEBUG) << _evMonoEnergyVec.size() <<"\t"<< minMonoEnergy <<","<< maxMonoEnergy <<endl;
        // Now filling the histograms
        for (double energy:_evMonoEnergyVec) {_evMonoEnergyHist->Fill(energy);}

        // Fitting
        _evMonoEnergyHist->Fit("gaus");
        TF1 *mono_fit = (TF1*) _evMonoEnergyHist->GetListOfFunctions()->FindObject("gaus");
        gStyle->SetOptFit(1111);
        streamlog_out(MESSAGE) << "\n Mono Fit Res. " << mono_fit->GetParameter(2) / mono_fit->GetParameter(1) <<endl;
    }
    
    if (_flagPixelEcalCol) {
        // Using the full statistical information to determine the binsize
        float meanEnergy = TMath::Mean(_evEnergyVec.begin(), _evEnergyVec.end());
        float sigmaEnergy = TMath::RMS(_evEnergyVec.begin(), _evEnergyVec.end());
        float meanHits = TMath::Mean(_evHitsVec.begin(), _evHitsVec.end());
        float sigmaHits = TMath::RMS(_evHitsVec.begin(), _evHitsVec.end());
        // Now declaring the histograms
        float minEnergy = meanEnergy - radiusOverSigma*sigmaEnergy;
        float maxEnergy = meanEnergy + radiusOverSigma*sigmaEnergy;
        float minHits = meanHits - radiusOverSigma*sigmaHits;
        float maxHits = meanHits + radiusOverSigma*sigmaHits;
        _evEnergyHist->SetBins(evHistBins, minEnergy, maxEnergy);
        _evHitsHist->SetBins(evHistBins, minHits, maxHits);
        streamlog_out(DEBUG) << _evEnergyVec.size() <<"\t"<< minEnergy <<","<< maxEnergy <<endl;
        // Now filling the histograms
        for (double energy:_evEnergyVec) {_evEnergyHist->Fill(energy);}
        for (int hits:_evHitsVec) {_evHitsHist->Fill(hits);}

        
        _evHitsHist->Fit("gaus");	//Fit a gaussian to the distribution
        TF1 *fit1 = (TF1*)_evHitsHist->GetListOfFunctions()->FindObject("gaus");
        gStyle->SetOptFit(1111);//Set to 1 to show and save the fit with the histogram in the root file generated by the AIDAProcessor
        // _evHitsHist->Fit("gaus","","",fit1->GetParameter(1)-fit1->GetParameter(2),fit1->GetParameter(1)+fit1->GetParameter(2));	//Fit a gaussian to the distribution
        // TF1 *fit2 = (TF1*)_evHitsHist->GetListOfFunctions()->FindObject("gaus");
        // gStyle->SetOptFit(1111);//Set to 1 to show and save the fit with the histogram in the root file generated by the AIDAProcessor
        _evEnergyHist->Fit("gaus");	//Fit a gaussian to the distribution
        TF1 *fit = (TF1*)_evEnergyHist->GetListOfFunctions()->FindObject("gaus");
        gStyle->SetOptFit(1111);//Set to 1 to show and save the fit with the histogram in the root file generated by the AIDAProcessor

        streamlog_out(MESSAGE) << "\n Pixel Fit Res. " << fit->GetParameter(2) / fit->GetParameter(1) <<endl;

        // Doing the fitting layer-by-layer seems unnecessary,
        // as the eventual Gaussian distribution comes from the large number of Laudau sampling.
        // for (int i = 0; i < NUMBER_OF_LAYER; i++) {
        //     // Fitting SNR histo
        //     printf("Fitting...\n");

        //     _hitsInLayer[i]->Fit("gaus");	
        //     TF1 *fitLayer1 = (TF1*)_hitsInLayer[i]->GetListOfFunctions()->FindObject("gaus");
        //     gStyle->SetOptFit(1111);
            
        //     _energyInLayerSi[i]->Fit("gaus");	
        //     TF1 *fitLayer = (TF1*)_energyInLayerSi[i]->GetListOfFunctions()->FindObject("gaus");
        //     gStyle->SetOptFit(1111);
        //     printf("Fitting done\nPlotting results...\n");
        // }
    }

    if (_flagDigitEcalCol) {
        // Using the full statistical information to determine the binsize
        float meanDigitEnergy = TMath::Mean(_evDigitEnergyVec.begin(), _evDigitEnergyVec.end());
        float sigmaDigitEnergy = TMath::RMS(_evDigitEnergyVec.begin(), _evDigitEnergyVec.end());
        float meanDigitHits = TMath::Mean(_evDigitHitsVec.begin(), _evDigitHitsVec.end());
        float sigmaDigitHits = TMath::RMS(_evDigitHitsVec.begin(), _evDigitHitsVec.end());
        // Now declaring the histograms
        float minDigitEnergy = meanDigitEnergy - radiusOverSigma*sigmaDigitEnergy;
        float maxDigitEnergy = meanDigitEnergy + radiusOverSigma*sigmaDigitEnergy;
        float minDigitHits = meanDigitHits - radiusOverSigma*sigmaDigitHits;
        float maxDigitHits = meanDigitHits + radiusOverSigma*sigmaDigitHits;
        _evDigitEnergyHist->SetBins(evHistBins, minDigitEnergy, maxDigitEnergy);
        _evDigitHitsHist->SetBins(evHistBins, minDigitHits, maxDigitHits);
        streamlog_out(MESSAGE) << _evDigitEnergyVec.size() <<"\t"<< minDigitEnergy <<","<< maxDigitEnergy <<endl;
        // Now filling the histograms
        for (double energy:_evDigitEnergyVec) {_evDigitEnergyHist->Fill(energy);}
        for (int hits:_evDigitHitsVec) {_evDigitHitsHist->Fill(hits);}

        
        _evDigitHitsHist->Fit("gaus");	//Fit a gaussian to the distribution
        TF1 *digit_fit1 = (TF1*)_evDigitHitsHist->GetListOfFunctions()->FindObject("gaus");
        gStyle->SetOptFit(1111);//Set to 1 to show and save the fit with the histogram in the root file generated by the AIDAProcessor
        _evDigitEnergyHist->Fit("gaus");	//Fit a gaussian to the distribution
        TF1 *digit_fit = (TF1*)_evDigitEnergyHist->GetListOfFunctions()->FindObject("gaus");
        gStyle->SetOptFit(1111);//Set to 1 to show and save the fit with the histogram in the root file generated by the AIDAProcessor

        streamlog_out(MESSAGE) << "\n Digit Fit Res. " << digit_fit->GetParameter(2) / digit_fit->GetParameter(1) <<endl;
    }
    streamlog_out(MESSAGE) << "Fitting done\nPlotting results..." <<endl;
}
