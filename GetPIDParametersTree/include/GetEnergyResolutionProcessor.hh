#ifndef GetEnergyResolutionProcessor_h
#define GetEnergyResolutionProcessor_h 1
#include <iomanip>
#include <EVENT/LCRelation.h>
#include "marlin/Processor.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/PIDHandler.h>
#include <UTIL/CellIDDecoder.h>
#include <cmath>
#include <string>
#include <vector>
#include <array>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include <TF1.h>
#include <TRandom.h>
#include <TFile.h>
#include <TTree.h>
// #include <TLorentzVector.h>
//MAS: Included this three lines
#ifdef MARLIN_AIDA //AIDA
#include <marlin/AIDAProcessor.h>
#endif

using namespace lcio;
using namespace marlin;
using namespace std;

class GetEnergyResolutionProcessor : public Processor {

private:
    const static int NUMBER_OF_LAYER = 15;
    const static int NUMBER_OF_CELLX = 64;
    const static int NUMBER_OF_CELLY = 32;
    float FIT_INIT_AMP = 5000;//events
    float FIT_INIT_MPV = 0.1e-3;//GeV
    float FIT_INIT_SIG = 0.01e-3;//GeV
    float E_RANGE_MIN = 0;
    float E_RANGE_MAX = 0.05;//GeV

public:
  virtual Processor *newProcessor() { return new GetEnergyResolutionProcessor; }

  GetEnergyResolutionProcessor();
  virtual ~GetEnergyResolutionProcessor();

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

  // Histogram definitions for ECALHit class

    // MC particle
    TH1* _runEnergy;
    // Monolithic calorimeter hits
    TH1* _zMonoHist;
    TH1* _evMonoEnergyHist;
    TH1* _energyInMonoLayerSi[NUMBER_OF_LAYER];
    // Pixelised calorimeter hits
    TH1* _xHist;
    TH1* _yHist;
    TH1* _zHist;
    TH2* _xyHist;
    TH2* _zxHist;
    TH2* _zyHist;
    TH1* _cellEnergyHist;
    TH1* _evEnergyHist;
    TH1* _evHitsHist;
    TH1* _energyInLayerSi[NUMBER_OF_LAYER];
    TH1* _hitsInLayer[NUMBER_OF_LAYER];
    // Digitised calorimeter hits
    TH1* _xDigitHist;
    TH1* _yDigitHist;
    TH1* _zDigitHist;
    TH2* _xyDigitHist;
    TH2* _zxDigitHist;
    TH2* _zyDigitHist;
    TH1* _cellDigitEnergyHist;
    TH1* _evDigitEnergyHist;
    TH1* _evDigitHitsHist;
    TH1* _energyInDigitLayerSi[NUMBER_OF_LAYER];
    TH1* _hitsInDigitLayer[NUMBER_OF_LAYER];
    // double _layerFitParams[NUMBER_OF_LAYER][4];
    // double energyRes;
    

private:
  virtual void ShowMCInfo(LCCollection *col);
  virtual void ShowECALInfo(LCCollection *col);
    virtual void ShowPixelECALInfo(LCCollection *col);
    virtual void ShowDigitECALInfo(LCCollection *col);

  std::string _MCColName;
  std::string _ECALColName;
    std::string _pECALColName;
    std::string _dECALColName;
    
    uint evHistBins = 31;
    float radiusOverSigma = 2.0;//Fitting range == mean +- radius
    vector<double> _evMonoEnergyVec;
    vector<double> _evEnergyVec;
    vector<int> _evHitsVec;
    vector<double> _evDigitEnergyVec;
    vector<int> _evDigitHitsVec;
    
    bool _flagMcCol = false;
    bool _flagEcalCol = false;
    bool _flagPixelEcalCol = false;
    bool _flagDigitEcalCol = false;
    float runEnergy = -1;


};

#endif
